/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2014, University of Twente
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors
 may be used to endorse or promote products derived from this software without
 specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <chrono>
#include <random>
#include "petsc.h"

#include "Base/CommandLineOptions.h"
#include "Base/H1ConformingTransformation.h"
#include "Base/MeshManipulator.h"
#include "Integration/ElementIntegral.h"
#include "Output/VTKSpecificTimeWriter.h"
#include "Utilities/GlobalIndexing.h"
#include "Utilities/GlobalMatrix.h"
#include "Utilities/GlobalVector.h"

#include "L2ProjectionErrorQualityMetric.h"
#include "SIPGErrorQualityMetric.h"

#include "Logger.h"

using namespace hpgem;

auto& name = Base::register_argument<std::string>(
    'm', "meshFile", "name of the mesh file", true);

template <std::size_t>
void runWithDimension();

int main(int argc, char** argv) {
    Base::parse_options(argc, argv);

    auto startTime = std::chrono::system_clock::now();

    runWithDimension<2>();

    auto endTime = std::chrono::system_clock::now();
    std::chrono::duration<double> duration = endTime - startTime;
    std::cout << "Elapsed time " << duration.count() << std::endl;

    return 0;
}

template <std::size_t dim>
std::vector<std::unique_ptr<QualityMetricComputation<dim>>> getMetrics(
    const Base::MeshManipulator<dim>& mesh);

template <std::size_t DIM>
void runTestsWithDGBasis(Base::MeshManipulator<DIM>& mesh);
template <std::size_t DIM>
void runTestsWithConformingBasis(Base::MeshManipulator<DIM>& mesh);

template <std::size_t DIM>
void runWithDimension() {
    std::size_t numberOfVectors = 1;
    Base::ConfigurationData configData(1, 0);

    Base::MeshManipulator<DIM> mesh(&configData, 1, 1, 1, 1);
    mesh.readMesh(name.getValue());

    std::vector<std::unique_ptr<QualityMetricComputation<DIM>>> metrics =
        getMetrics<DIM>(mesh);
    Output::VTKSpecificTimeWriter<DIM> plotter("meshQuality", &mesh, 0, 1);
    for (auto& metric : metrics) {
        metric->computeAndPlotMetric(mesh, plotter, "meshQuality-");
    }
}

template <std::size_t DIM>
double testingFunction(const Geometry::PointPhysical<DIM>& p) {
    double functionValue = 0;
    for (std::size_t i = 0; i < DIM; ++i) {
        functionValue += (i + 1) * p[i];
    }
    return functionValue;
}

template <std::size_t DIM>
void runTestsWithDGBasis(Base::MeshManipulator<DIM>& mesh) {
    const std::size_t ELEM_MAT_ID = 0;
    const std::size_t ELEM_VEC_ID = 0;

    // Linear patch test
    mesh.useDefaultDGBasisFunctions(1);
    Integration::ElementIntegral<DIM> integral;
    integral.setTransformation(
        std::shared_ptr<Base::CoordinateTransformation<DIM>>(
            new Base::H1ConformingTransformation<DIM>()));

    std::vector<double> phis;

    for (Base::Element* element : mesh.getElementsList()) {
        LinearAlgebra::MiddleSizeMatrix mat;
        // Compute mass matrix
        mat = integral.integrate(
            element, [&phis](Base::PhysicalElement<DIM>& pelemnt) {
                std::size_t n = pelemnt.getNumberOfBasisFunctions();
                LinearAlgebra::MiddleSizeMatrix mat(n, n);
                phis.resize(n);
                for (std::size_t i = 0; i < n; ++i) {
                    phis[i] = pelemnt.basisFunction(i);
                }
                for (std::size_t i = 0; i < n; ++i) {
                    for (std::size_t j = i; j < n; ++j) {
                        mat(i, j) = phis[i] * phis[j];
                        mat(j, i) = mat(i, j);
                    }
                }
                return mat;
            });
        element->setElementMatrix(mat, ELEM_MAT_ID);
        // Compute source vector
        LinearAlgebra::MiddleSizeVector vec;
        vec = integral.integrate(
            element, [](Base::PhysicalElement<DIM>& pelement) {
                std::size_t n = pelement.getNumberOfBasisFunctions();
                LinearAlgebra::MiddleSizeVector vec(n);
                // Position for the patch test
                // x + 2y + 3z etc.
                const Geometry::PointPhysical<DIM>& p =
                    pelement.getPointPhysical();
                double functionValue = testingFunction(p);
                for (std::size_t i = 0; i < n; ++i) {
                    vec[i] = functionValue * pelement.basisFunction(i);
                }
                return vec;
            });
        element->setElementVector(vec, ELEM_VEC_ID);
    }
    // Setup global system
    Utilities::GlobalIndexing indexing(&mesh);
    Utilities::GlobalPetscMatrix massMatrix(indexing, ELEM_MAT_ID, -1);
    Utilities::GlobalPetscVector loadVector(indexing, ELEM_VEC_ID, -1);
    Utilities::GlobalPetscVector resultVector(indexing, -1, -1);
    // Solve system
    KSP ksp;
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPPREONLY);

    KSPSetOperators(ksp, massMatrix, massMatrix);
    KSPSolve(ksp, loadVector, resultVector);

    // Reserve time integration vector
    for (Base::Element* element : mesh.getElementsList()) {
        element->setNumberOfTimeIntegrationVectors(1);
    }
    // Computing the error
    resultVector.writeTimeIntegrationVector(0);

    std::map<const Base::Element*, double> errors;

    for (Base::Element* element : mesh.getElementsList()) {
        LinearAlgebra::MiddleSizeVector coefficients =
            element->getTimeIntegrationVector(0);

        double localError = integral.integrate(
            element, [&coefficients](Base::PhysicalElement<DIM>& pelement) {
                std::size_t n = pelement.getNumberOfBasisFunctions();
                const Geometry::PointPhysical<DIM>& p =
                    pelement.getPointPhysical();
                LinearAlgebra::MiddleSizeVector::type error =
                    testingFunction(p);
                for (std::size_t i = 0; i < n; ++i) {
                    error -= coefficients[i] * pelement.basisFunction(i);
                }
                double result = std::real(error) * std::real(error);
                result += std::imag(error) * std::imag(error);
                return result;
            });
        errors[element] = localError;
    }

    // Plotting
    Output::VTKSpecificTimeWriter<DIM> writer("meshquality", &mesh);
    writer.write(
        [&errors](Base::Element* element, const Geometry::PointReference<DIM>&,
                  std::size_t) { return errors[element]; },
        "error");
}

template <std::size_t dim>
class SinCosFunc : public Func<dim> {

   public:
    SinCosFunc(LinearAlgebra::SmallVector<dim>& k, double phase, bool sinFunc)
        : k_(k), sinFunc_(sinFunc), phase_ (phase){};

    double value(const Geometry::PointPhysical<dim>& p) const final {
        if (sinFunc_) {
            return std::sin(p.getCoordinates() * k_ + phase_);
        } else {
            return std::cos(p.getCoordinates() * k_ + phase_);
        }
    }

    double laplacian(const Geometry::PointPhysical<dim>& p) const final {
        return -k_.l2NormSquared() * value(p);
    }

   private:
    bool sinFunc_;
    LinearAlgebra::SmallVector<dim> k_;
    double phase_;
};

template <std::size_t dim>
class TestingFunc : public Func<dim> {
   public:
    double value(const Geometry::PointPhysical<dim>& p) const final {
        return p[0] * (1 - p[0]) * p[1] * (1 - p[1]);
        //        return 1.0;
    }
    double laplacian(const Geometry::PointPhysical<dim>& p) const final {
        return -2 * (p[0] * (1 - p[0]) + p[1] * (1 - p[1]));
        //        return -0.0;
    }
};

template <std::size_t dim>
std::vector<std::unique_ptr<QualityMetricComputation<dim>>> getMetrics(
    const Base::MeshManipulator<dim>& mesh) {
    std::vector<std::unique_ptr<QualityMetricComputation<dim>>> metrics;

    // L2-projection for exp(ik.x), using both dg and cg basis functions.  Split
    // over two source functions, sin and cos.
    {
        using Function =
            typename L2ProjectionErrorQualityMetric<dim>::TestingFunction;
        std::vector<Function> functions;
        LinearAlgebra::SmallVector<dim> k;
        // Compute k such that we have exactly 1 period over the mesh bounding
        // box
        LinearAlgebra::SmallVector<dim> minVec, maxVec;
        minVec.set(std::numeric_limits<double>::max());
        maxVec.set(std::numeric_limits<double>::min());
        for (const Geometry::PointPhysical<dim>& p :
             mesh.getNodeCoordinates()) {
            for (std::size_t i = 0; i < dim; ++i) {
                minVec[i] = std::min(minVec[i], p[i]);
                maxVec[i] = std::max(maxVec[i], p[i]);
            }
        }
        k = maxVec - minVec;
        double maxK = 0.0;
        for (std::size_t i = 0; i < dim; ++i) {
            double value = 2 * M_PI / k[i];
            k[i] = value;
            maxK = std::max(value, maxK);
        }

        functions.push_back([k](const Geometry::PointPhysical<dim>& p) {
            return std::cos(p.getCoordinates() * k);
        });
        functions.push_back([k](const Geometry::PointPhysical<dim>& p) {
            return std::sin(p.getCoordinates() * k);
        });
        using BFType =
            typename L2ProjectionErrorQualityMetric<dim>::BasisFunctionType;

        metrics.emplace_back(new L2ProjectionErrorQualityMetric<dim>(
            BFType::DISCONTINUOUS, 1, "L2-projection-dg1", functions));
        metrics.emplace_back(new L2ProjectionErrorQualityMetric<dim>(
            BFType::CONFORMING, 1, "L2-projection-cg1", functions));
        metrics.emplace_back(new L2ProjectionErrorQualityMetric<dim>(
            BFType::DISCONTINUOUS, 2, "L2-projection-dg2", functions));
        metrics.emplace_back(new L2ProjectionErrorQualityMetric<dim>(
            BFType::CONFORMING, 2, "L2-projection-cg2", functions));

        // Several pairs of the form {sin(kx), cos(kx)} for different k.

        const std::size_t NUM_DIRECTIONS = 10;
        std::vector<std::shared_ptr<Func<dim>>> sipgFunctions;

        std::mt19937 generator;  // Use fixed default seed
        std::uniform_real_distribution<double> distribution(-1.0, 1.0);

        for (std::size_t i = 0; i < NUM_DIRECTIONS; ++i) {
            LinearAlgebra::SmallVector<dim> localK;
            for (std::size_t j = 0; j < dim; ++j) {
                localK[j] = 2 * maxK * distribution(generator);
            }
            double phase = M_PI*distribution(generator);
            sipgFunctions.push_back(
                std::make_shared<SinCosFunc<dim>>(localK, phase, true));
            sipgFunctions.push_back(
                std::make_shared<SinCosFunc<dim>>(localK, phase, false));
        }
        metrics.emplace_back(
            new SIPGErrorQualityMetric<dim>(1, "SIPG-linear", sipgFunctions));
    }

    return metrics;
}