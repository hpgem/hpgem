/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2021, University of Twente
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
#ifndef HPGEM_GLOBALSOLVEQUALITYMETRIC_H
#define HPGEM_GLOBALSOLVEQUALITYMETRIC_H

#include "QualityMetricComputation.h"
#include <Utilities/GlobalMatrix.h>
#include <Utilities/GlobalVector.h>

/**
 * Error metric that is based on solving a global system.
 *
 * This abstract class for an error metric assumes that the approach satisfies
 * several properties:
 *
 *  1. The metric is based on solving a system of the form Sx = b for one or
 *     more vectors b, a FEM matrix S (e.g. mass or stiffness) and then
 * computing
 *  2. The error is the L2 error between the computed solution and a known
 *     solution.
 *  3. When multiple right hand sides are in use the total error is the l2-norm
 *     of the individual errors.
 *
 *
 * @tparam dim
 */
template <std::size_t dim>
class GlobalSolveQualityMetric : public QualityMetricComputation<dim> {
   protected:
    // Time integration vector to use for the solution
    const std::size_t SOLUTION_VECTOR_ID = 0;
    const std::size_t ELEMENT_MATRIX_ID = 0;
    const std::size_t ELEMENT_VECTOR_ID = 0;
    const std::size_t FACE_MATRIX_ID = 0;
    const std::size_t FACE_VECTOR_ID = 0;
    // Type of the value from solutions, may be complex.
    using VALUE_TYPE = hpgem::LinearAlgebra::MiddleSizeVector::type;

   public:
    void computeAndPlotMetric(
        hpgem::Base::MeshManipulator<dim>& mesh,
        hpgem::Output::VTKSpecificTimeWriter<dim>& plotter,
        const std::string& filePrefix) final;

    /**
     * @return  The order of the basis functions
     */
    std::size_t getOrder() const { return order_; }

   protected:
    GlobalSolveQualityMetric(std::size_t order, std::string name)
        : name_(std::move(name)), order_(order){};

    /**
     * @return Whether face matrices and/or vectors are in use, or only element
     * matrices and vectors are used.
     */
    virtual bool usesFaceMatrixAndVector() const = 0;
    /**
     * @return  The number of right hand sides and thus solves.
     */
    virtual std::size_t numberOfSolves() const = 0;
    /**
     * Compute the local matrices, both element and if in use face.
     *
     * NOTE: Requires initializing the basis functions.
     * @param mesh The mesh to use
     */
    virtual void computeLocalMatrices(
        hpgem::Base::MeshManipulator<dim>& mesh) = 0;
    /**
     * Compute the local vectors
     * @param solve The id of the solve
     * @param mesh The mesh to use
     */
    virtual void computeLocalVectors(
        std::size_t solve, hpgem::Base::MeshManipulator<dim>& mesh) = 0;

    /**
     * Compute the square of the L2 error on an element
     * @param solve The index of the solve
     * @param element The element to compute on
     * @return The square of the L2 error
     */
    double computeElementSquareError(std::size_t solve,
                                     hpgem::Base::Element* element);

    /**
     * Compute the solution at a point in an element
     * @param element The element to compute the solution on
     * @param p The reference coordinate on that element
     * @return The solution at that point
     */
    virtual double computeSolutionValue(
        hpgem::Base::Element* element,
        const hpgem::Geometry::PointReference<dim>& p) const = 0;
    /**
     * Compute the theoretical solution
     * @param solve The index of the solve
     * @param p The point
     * @return The function value at the point
     */
    virtual double computeFunctionValue(
        std::size_t solve,
        const hpgem::Geometry::PointPhysical<dim>& p) const = 0;
    /**
     * @return The integrator to be used for the L2 error
     */
    virtual hpgem::Integration::ElementIntegral<dim>& getIntegrator() = 0;

   private:
    std::string name_;
    std::size_t order_;
};

template <std::size_t dim>
void GlobalSolveQualityMetric<dim>::computeAndPlotMetric(
    hpgem::Base::MeshManipulator<dim>& mesh,
    hpgem::Output::VTKSpecificTimeWriter<dim>& plotter,
    const std::string& filePrefix) {

    using namespace hpgem;

    if (numberOfSolves() == 0) {
        return;
    }

    logger(INFO, "Computing local matrices");

    computeLocalMatrices(mesh);

    // Initialize vectors with appropriately sized zero vectors so that global
    // vector assembly works.
    for (Base::Element* element : mesh.getElementsList()) {
        LinearAlgebra::MiddleSizeVector vec(
            element->getTotalNumberOfBasisFunctions());
        element->setElementVector(vec, ELEMENT_VECTOR_ID);
    }
    if (usesFaceMatrixAndVector()) {
        for (Base::Face* face : mesh.getFacesList()) {
            std::size_t n = face->getNumberOfBasisFunctions();
            LinearAlgebra::MiddleSizeVector vec(n);
            face->setFaceVector(vec, FACE_VECTOR_ID);
        }
    }

    // Build global matrices and vectors
    logger(INFO, "Assembling global matrix");
    Utilities::GlobalIndexing indexing(&mesh);
    Utilities::GlobalPetscMatrix globalMatrix(
        indexing, ELEMENT_MATRIX_ID,
        usesFaceMatrixAndVector() ? FACE_MATRIX_ID : -1);
    Utilities::GlobalPetscVector loadVector(
        indexing, ELEMENT_VECTOR_ID,
        usesFaceMatrixAndVector() ? FACE_VECTOR_ID : -1);
    Utilities::GlobalPetscVector resultVector(indexing, -1, -1);

    // Setup the solver
    KSP ksp;
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPGMRES);
    KSPSetOperators(ksp, globalMatrix, globalMatrix);
    KSPSetTolerances(ksp, 1e-12, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSetFromOptions(ksp);
    KSPSetUp(ksp);

    std::map<const Base::Element*, double> totalElementErrorSquared;
    std::map<const Base::Element*, double> elementErrorSquared;

    for (std::size_t sid = 0; sid < numberOfSolves(); ++sid) {
        logger(INFO, "Creating sources for %-%", name_, sid);
        computeLocalVectors(sid, mesh);

        loadVector.reinit();

        globalMatrix.writeMatlab("Stiff.m");
        loadVector.writeMatlab("Load.m");

        logger(INFO, "Solving for %-%", name_, sid);
        KSPSetOperators(ksp, globalMatrix, globalMatrix);
        KSPSolve(ksp, loadVector, resultVector);
        {
            // Basic solve diagnostics
            KSPConvergedReason reasonId;
            const char* reason;
            KSPGetConvergedReason(ksp, &reasonId);  // For convergence or not
            KSPGetConvergedReasonString(ksp, &reason);
            PetscInt iterations;
            KSPGetTotalIterations(ksp, &iterations);
            if (reasonId > 0) {
                logger(INFO, "Converged in % iterations with %", iterations,
                       reason);
            } else {
                logger(ERROR, "Not converged after % iterations with %",
                       iterations, reason);
            }
        }

        resultVector.writeTimeIntegrationVector(SOLUTION_VECTOR_ID);
        for (Base::Element* element : mesh.getElementsList()) {
            double error = computeElementSquareError(sid, element);
            totalElementErrorSquared[element] += error;
            elementErrorSquared[element] = error;
        }

        // Plot details

        std::stringstream detailPlotFile;
        detailPlotFile << filePrefix << name_ << "-" << sid;
        Output::VTKSpecificTimeWriter<dim> detailPlotter(detailPlotFile.str(),
                                                         &mesh, 0, order_);

        // The solution
        detailPlotter.write(
            [&sid, this](Base::Element* element,
                         const Geometry::PointReference<dim>& p, std::size_t) {
                return computeSolutionValue(element, p);
            },
            "solution");
        // The original function
        detailPlotter.write(
            [&sid, this](Base::Element* element,
                         const Geometry::PointReference<dim>& p, std::size_t) {
                return computeFunctionValue(sid,
                                            element->referenceToPhysical(p));
            },
            "function");
        // The error in the solution
        detailPlotter.write(
            [&sid, this](Base::Element* element,
                         const Geometry::PointReference<dim>& p, std::size_t) {
                double error = computeSolutionValue(element, p);
                error -=
                    computeFunctionValue(sid, element->referenceToPhysical(p));
                return error;
            },
            "error");

        // Local error
        detailPlotter.write(
            [&elementErrorSquared](Base::Element* element,
                                   const Geometry::PointReference<dim>&,
                                   std::size_t) {
                return std::sqrt(elementErrorSquared[element]);
            },
            "elementError");
    }
    // Plot the error in the result
    plotter.write(
        [&totalElementErrorSquared](Base::Element* element,
                                    const Geometry::PointReference<dim>&,
                                    std::size_t) {
            return std::sqrt(totalElementErrorSquared[element]);
        },
        name_);
}

template <std::size_t dim>
double GlobalSolveQualityMetric<dim>::computeElementSquareError(
    std::size_t solve, hpgem::Base::Element* element) {
    using namespace hpgem;
    return getIntegrator().integrate(
        element,
        [&](Base::PhysicalElement<dim>& pelement) -> double {
            double error =
                computeSolutionValue(element, pelement.getPointReference());
            error -= computeFunctionValue(solve, pelement.getPointPhysical());
            return error * error;
        },
        // Higher quadrature order than needed to improve accuracy
        QuadratureRules::AllGaussQuadratureRules::instance().getRule(
            element->getReferenceGeometry(), 2 * order_ + 4));
}

#endif  // HPGEM_GLOBALSOLVEQUALITYMETRIC_H
