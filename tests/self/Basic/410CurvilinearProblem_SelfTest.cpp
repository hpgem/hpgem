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

#include <chrono>

#include <Logger.h>
#include <Base/CommandLineOptions.h>
#include <Base/HpgemAPILinearSteadyState.h>

#include "../ConvergenceTest.h"
#include "../TestMeshes.h"

using namespace hpgem;

class Problem : public Base::HpgemAPILinearSteadyState<2> {
   public:
    Problem(std::size_t order, double penalty)
        : Base::HpgemAPILinearSteadyState<2>(1, order, true, true),
          penalty_(penalty) {
        writeTec_ = false;
    };

    void readMesh(const std::string meshName) override {
        HpgemAPILinearSteadyState::readMesh(meshName);
        std::size_t elementMatrices = 1;
        std::size_t elementVectors = 1;
        std::size_t faceMatrices = 1;
        std::size_t faceVectors = 1;
        // No faceVectors as nothing at the boundary
        this->addMesh(meshName, elementMatrices, elementVectors, faceMatrices,
                      faceVectors);
        this->meshes_[0]->useDefaultDGBasisFunctions(polynomialOrder_);
        // Only needed for the solution
        this->setNumberOfTimeIntegrationVectorsGlobally(1);
    }

    LinearAlgebra::MiddleSizeMatrix computeIntegrandStiffnessMatrixAtElement(
        Base::PhysicalElement<2>& element) override {
        auto nbasis = element.getNumberOfBasisFunctions();
        LinearAlgebra::MiddleSizeMatrix& integrandVal =
            element.getResultMatrix();

        for (std::size_t i = 0; i < nbasis; ++i) {
            for (std::size_t j = 0; j < nbasis; ++j) {
                integrandVal(i, j) = element.basisFunctionDeriv(i) *
                                     element.basisFunctionDeriv(j);
            }
        }
        return integrandVal;
    }

    Base::FaceMatrix computeIntegrandStiffnessMatrixAtFace(
        Base::PhysicalFace<2>& face) override {
        std::size_t numberOfBasisFunctions =
            face.getFace()->getNumberOfBasisFunctions();

        // Create the FaceMatrix integrandVal with the correct size.
        Base::FaceMatrix& integrandVal = face.getResultMatrix();

        LinearAlgebra::SmallVector<2> phiNormalI, phiNormalJ, phiDerivI,
            phiDerivJ;

        double penalty = penalty_ / face.getFace()->getDiameter();
        double factor = face.isInternal() ? 0.5 : 1.0;

        for (std::size_t i = 0; i < numberOfBasisFunctions; ++i) {
            phiNormalI = face.basisFunctionUnitNormal(i);
            phiDerivI = face.basisFunctionDeriv(i);

            for (std::size_t j = 0; j < numberOfBasisFunctions; ++j) {
                phiNormalJ = face.basisFunctionUnitNormal(j);
                phiDerivJ = face.basisFunctionDeriv(j);

                integrandVal(j, i) =
                    -(phiNormalI * phiDerivJ + phiNormalJ * phiDerivI) *
                        factor +
                    penalty * phiNormalI * phiNormalJ;
            }
        }

        return integrandVal;
    }

    LinearAlgebra::MiddleSizeVector computeIntegrandSourceTermAtFace(
        Base::PhysicalFace<2>& face) final {
        auto numdofs = face.getNumberOfBasisFunctions();
        LinearAlgebra::MiddleSizeVector res(numdofs);

        if (face.isInternal()) {
            return res;
        }

        double penalty = penalty_ / face.getFace()->getDiameter();
        const auto& normal = face.getUnitNormalVector();

        double value = 1.0;

        for (std::size_t i = 0; i < numdofs; ++i) {
            res(i) = normal *
                     (penalty * face.basisFunctionUnitNormal(i) -
                      face.basisFunctionDeriv(i)) *
                     value;
        }
        // By strange convention this there is no -1 in front of the Laplacian
        // requiring compensation on the right hand side
        res *= -1.0;
        return res;
    }

    LinearAlgebra::MiddleSizeVector getExactSolution(
        const PointPhysicalT& pPhys) override {
        // 1 - r^2
        return {2.0 - pPhys.l2NormSquared()};
    }

    LinearAlgebra::MiddleSizeVector getSourceTerm(
        const PointPhysicalT& pPhys) override {
        // Laplacian -r^2 = 4
        // By strange convention this there is no -1 in front of the Laplacian
        // requiring compensation on the right hand side
        return {-4.0};
    }

   private:
    double penalty_;
};

double solve(const std::string& filename, std::size_t level) {

    std::stringstream outputFile;
    outputFile << "solution-" << level << "-";

    Problem problem(3, 10.0);
    problem.setOutputFileName(outputFile.str());

    problem.readMesh(filename);
    problem.solveSteadyStateWithPetsc(false);
    return problem.computeTotalError(0, 0).real();
}

int main(int argc, char** argv) {
    using namespace std::string_literals;
    Base::parse_options(argc, argv);

    // Define clocks for measuring simulation time.
    std::chrono::time_point<std::chrono::system_clock> startClock, endClock;
    startClock = std::chrono::system_clock::now();

    // For recomputing the error tables
    bool ignoreFailures = true;

    // Convergence of 2-r^2 on a unit circle with second order basis functions.
    // Were it not for the curvilinear elements, the solution would be perfectly
    // representable. The error decreases with a factor of about 9-11,
    // quicker than the 2^{p+1} that is expected in general.
    ConvergenceTestSet set = {getUnitCircleQuadraticTriangleMeshes(),
                              {
                                  5.12845861e-03,  //------
                                  5.51975907e-04,  //  9.29
                                  5.52196123e-05,  // 10.00
                                  5.18265586e-06,  // 10.65
                              }};
    runConvergenceTest(set, ignoreFailures, solve);

    endClock = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = endClock - startClock;
    logger(INFO, "Elapsed time %", elapsed_seconds.count());

    return 0;
}