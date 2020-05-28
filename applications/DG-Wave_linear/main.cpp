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

#include "Base/HpgemAPIBase.h"
#include "Base/L2Norm.h"
#include "Output/TecplotSingleElementWriter.h"
#include <fstream>
#include "Utilities/BasisFunctions2DH1ConformingSquare.h"
#include "Integration/ElementIntegrandBase.h"
#include "Integration/FaceIntegrandBase.h"
#include "Integration/FaceIntegral.h"
#include "Integration/ElementIntegral.h"
#include "Base/GlobalData.h"
#include "Base/ConfigurationData.h"
#include "Base/Element.h"
#include "Base/Face.h"
#include "Output/TecplotDiscontinuousSolutionWriter.h"
#include "Output/VTKTimeDependentWriter.h"
#include <cmath>

#include "Utilities/GlobalMatrix.h"
#include "Utilities/GlobalVector.h"
#include "petscksp.h"

const std::size_t DIM = 2;

class DGWave : public Base::HpgemAPIBase<DIM>,
               public Output::TecplotSingleElementWriter<DIM> {
   public:
    DGWave(std::size_t p)
        : HpgemAPIBase<DIM>(new Base::GlobalData(),
                            new Base::ConfigurationData(2)),
          p_(p) {}

    /// set up the mesh
    bool initialise() {
        addMesh("mesh.hpgem", 1, 1, 1, 1);
        setNumberOfTimeIntegrationVectorsGlobally(1);
        meshes_[0]->setDefaultBasisFunctionSet(
            Utilities::createInteriorBasisFunctionSet2DH1Square(p_));
        std::vector<const Base::BasisFunctionSet*> bFsets;
        bFsets = Utilities::createVertexBasisFunctionSet2DH1Square(p_);
        meshes_[0]->addVertexBasisFunctionSet(bFsets);
        std::vector<const Base::OrientedBasisFunctionSet*> oBFsets;
        oBFsets = Utilities::createFaceBasisFunctionSet2DH1Square(p_);
        meshes_[0]->addFaceBasisFunctionSet(oBFsets);
        return true;
    }

    class : public Integration::ElementIntegrandBase<
                LinearAlgebra::MiddleSizeMatrix, DIM> {
       public:
        void elementIntegrand(
            Base::PhysicalElement<DIM>& element,
            LinearAlgebra::MiddleSizeMatrix& ret) override final {
            std::size_t numberOfBasisFunctions =
                element.getElement()->getNumberOfBasisFunctions();
            ret.resize(numberOfBasisFunctions *
                           element.getElement()->getNumberOfUnknowns(),
                       numberOfBasisFunctions *
                           element.getElement()->getNumberOfUnknowns());
            for (std::size_t i = 0; i < numberOfBasisFunctions; ++i) {
                for (std::size_t j = 0; j <= i; ++j) {
                    ret(element.getElement()->convertToSingleIndex(i, 1),
                        element.getElement()->convertToSingleIndex(j, 1)) =
                        ret(element.getElement()->convertToSingleIndex(j, 1),
                            element.getElement()->convertToSingleIndex(i, 1)) =
                            element.basisFunctionDeriv(i) *
                            element.basisFunctionDeriv(j);
                }
            }
        }
    } stifnessIntegrand;

    class : public Integration::FaceIntegrandBase<
                LinearAlgebra::MiddleSizeMatrix, DIM> {
       public:
        void faceIntegrand(
            Base::PhysicalFace<DIM>& face,
            LinearAlgebra::MiddleSizeMatrix& ret) override final {
            std::size_t numberOfBasisFunctions =
                face.getFace()->getNumberOfBasisFunctions();
            ret.resize(
                numberOfBasisFunctions *
                    face.getFace()->getPtrElementLeft()->getNumberOfUnknowns(),
                numberOfBasisFunctions *
                    face.getFace()->getPtrElementLeft()->getNumberOfUnknowns());
            for (std::size_t i = 0; i < numberOfBasisFunctions; ++i) {
                for (std::size_t j = 0; j < numberOfBasisFunctions;
                     ++j) {  // the basis functions belonging to internal
                             // parameters are 0 on the free surface anyway.
                    for (std::size_t k = 0; k < face.getFace()
                                                    ->getPtrElementLeft()
                                                    ->getNumberOfUnknowns();
                         ++k) {
                        ret(face.getFace()
                                ->getPtrElementLeft()
                                ->convertToSingleIndex(i, k),
                            face.getFace()
                                ->getPtrElementLeft()
                                ->convertToSingleIndex(j, k)) =
                            face.basisFunction(i) * face.basisFunction(j);
                    }
                }
            }
        }
    } massIntegrand;

    static void exactSolution(const double t, const PointPhysicalT& p,
                              LinearAlgebra::MiddleSizeVector& ret) {
        ret.resize(2);
        // ret[0]=cos(2*M_PI*p[0])*cos(sqrt(2*M_PI*tanh(2*M_PI))*t)*cosh(2*M_PI*(p[DIM-1]+1))/cosh(2*M_PI);//standing
        // wave
        // ret[1]=cos(2*M_PI*p[0])*sin(-sqrt(2*M_PI*tanh(2*M_PI))*t)/sqrt(2*M_PI*tanh(2*M_PI))*cosh(2*M_PI*(p[DIM-1]+1))/cosh(2*M_PI);
        double k = 2. * M_PI;
        ret[0] = cosh(k * (p[DIM - 1] + 1)) *
                 sin(sqrt(k * tanh(k)) * t - k * p[0]) * sqrt(k * tanh(k)) /
                 cosh(k) * 0.001;  // moving wave
        ret[1] = cosh(k * (p[DIM - 1] + 1)) *
                 cos(sqrt(k * tanh(k)) * t - k * p[0]) / cosh(k) * 0.001;
    }

    static void initialConditions(const PointPhysicalT& p,
                                  LinearAlgebra::MiddleSizeVector& ret) {
        exactSolution(0, p, ret);
        // ret.resize(2);
        // ret[0]=cos(2*M_PI*p[0])*cosh(2*M_PI*(p[DIM-1]+1))/cosh(2*M_PI);//standing
        // wave ret[1]=0; double k = 2. * M_PI; ret[0] = cosh(k * (p[DIM - 1] +
        // 1)) * sin(-k * p[0]) * sqrt(k * tanh(k)) / cosh(k) * 0.001; //moving
        // wave ret[1] = cosh(k * (p[DIM - 1] + 1)) * cos(-k * p[0]) / cosh(k) *
        // 0.001;
    }

    class : public Integration::FaceIntegrandBase<
                LinearAlgebra::MiddleSizeVector, DIM> {
       public:
        void faceIntegrand(
            Base::PhysicalFace<DIM>& face,
            LinearAlgebra::MiddleSizeVector& ret) override final {
            const PointPhysicalT& pPhys = face.getPointPhysical();
            ret.resize(2 * face.getFace()->getNumberOfBasisFunctions());
            LinearAlgebra::MiddleSizeVector data;
            initialConditions(pPhys, data);
            for (std::size_t i = 0;
                 i < face.getFace()->getNumberOfBasisFunctions(); ++i) {
                ret(face.getFace()->getPtrElementLeft()->convertToSingleIndex(
                    i, 0)) = face.basisFunction(i) * data[0] * 2;
                ret(face.getFace()->getPtrElementLeft()->convertToSingleIndex(
                    i, 1)) = face.basisFunction(i) * data[1] * 2;
            }
        }
    } interpolator;

    class : public Integration::FaceIntegrandBase<
                LinearAlgebra::MiddleSizeVector, DIM> {
       public:
        void faceIntegrand(
            Base::PhysicalFace<DIM>& face,
            LinearAlgebra::MiddleSizeVector& ret) override final {
            const PointPhysicalT& pPhys = face.getPointPhysical();
            ret.resize(2);
            static LinearAlgebra::MiddleSizeVector exact(2);
            if (std::abs(pPhys[DIM - 1]) < 1e-9) {
                ret = face.getSolution(Base::Side::LEFT);
                exactSolution(t, pPhys, exact);
                ret -= exact;
                ret[0] *= ret[0];
                ret[1] *= ret[1];
            }
        }
    } error;

    class : public Integration::FaceIntegrandBase<
                LinearAlgebra::MiddleSizeVector, DIM> {
       public:
        void faceIntegrand(
            Base::PhysicalFace<DIM>& face,
            LinearAlgebra::MiddleSizeVector& ret) override final {
            const PointPhysicalT& pPhys = face.getPointPhysical();
            ret.resize(1);
            static LinearAlgebra::MiddleSizeVector dummySolution(2);
            if (std::abs(pPhys[DIM - 1]) < 1e-9) {
                dummySolution = face.getSolution(Base::Side::LEFT);
                ret[0] = dummySolution[0] * dummySolution[0];
                ret[0] /= 2;  // assumes g=1
            }
        }
    } faceEnergy;

    class : public Integration::ElementIntegrandBase<
                LinearAlgebra::MiddleSizeVector, DIM> {
       public:
        void elementIntegrand(
            Base::PhysicalElement<DIM>& element,
            LinearAlgebra::MiddleSizeVector& ret) override final {
            ret.resize(1);
            LinearAlgebra::SmallVector<DIM> gradPhi =
                element.getSolutionDeriv()[1];
            ret[0] = (gradPhi * gradPhi) / 2;
        }
    } elementEnergy;

    void writeToTecplotFile(const Base::Element* element,
                            const PointReferenceT& p,
                            std::ostream& out) override final {
        LinearAlgebra::MiddleSizeVector value = element->getSolution(0, p);
        out << value[0] << " " << value[1];
        PointPhysicalT pPhys = element->referenceToPhysical(p);
        exactSolution(t, pPhys, value);
        out << " " << value[0] << " " << value[1];
    }

    // assumes that 'all' means 'all relevant'; computes the mass matrix at the
    // free surface
    void doAllFaceIntegration() {
        Integration::FaceIntegral<DIM> integral;
        PointPhysicalT pPhys;
        LinearAlgebra::MiddleSizeMatrix result;
        LinearAlgebra::MiddleSizeVector initialconditions;
        for (Base::Face* face : meshes_[0]->getFacesList()) {
            const PointReferenceOnFaceT& p =
                face->getReferenceGeometry()->getCenter();
            pPhys = face->referenceToPhysical(p);
            if (std::abs(pPhys[DIM - 1]) < 1e-9) {
                result = integral.integrate(face, &massIntegrand);
                initialconditions = integral.integrate(face, &interpolator);
                face->getPtrElementLeft()->setTimeIntegrationVector(
                    0, initialconditions);
            } else {
                result.resize(face->getNumberOfBasisFunctions() * 2,
                              face->getNumberOfBasisFunctions() * 2);
                for (std::size_t i = 0; i < result.size(); ++i) {
                    result[i] = 0;
                }
            }
            face->setFaceMatrix(result);
        }
    }

    void doAllElementIntegration() {
        Integration::ElementIntegral<DIM> integral;
        LinearAlgebra::MiddleSizeMatrix result;
        LinearAlgebra::MiddleSizeVector zero;
        for (Base::Element* element : meshes_[0]->getElementsList()) {
            zero.resize(element->getNumberOfUnknowns() *
                        element->getNumberOfBasisFunctions());
            result = integral.integrate(element, &stifnessIntegrand);
            element->setElementMatrix(result);
            element->setTimeIntegrationVector(0, zero);
        }
    }

    // messy routine that collects the row numbers of basisfunctions active on
    // the free surface for use with PETSc
    void getSurfaceIS(Utilities::GlobalPetscMatrix& S, IS* surface, IS* rest) {
        PointPhysicalT pPhys;
        std::size_t numberOfBasisFunctions(0);
        std::vector<int> facePositions;
        for (const Base::Face* face : meshes_[0]->getFacesList()) {
            const PointReferenceOnFaceT& p =
                face->getReferenceGeometry()->getCenter();
            pPhys = face->referenceToPhysical(p);
            if (std::abs(pPhys[DIM - 1]) < 1e-9) {
                S.getMatrixBCEntries(face, numberOfBasisFunctions,
                                     facePositions);
            }
        }
        ISCreateGeneral(PETSC_COMM_WORLD, numberOfBasisFunctions,
                        facePositions.data(), PETSC_COPY_VALUES, surface);
        int totalFuncs;
        MatGetSize(S, &totalFuncs, PETSC_NULL);
        ISSortRemoveDups(*surface);
        ISComplement(*surface, 0, totalFuncs, rest);
    }

    void printError() {
        Integration::FaceIntegral<DIM> integral;
        LinearAlgebra::MiddleSizeVector totalError(2), contribution(2);
        for (Base::Face* face : meshes_[0]->getFacesList()) {
            contribution = integral.integrate(face, &error);
            totalError += contribution;
            contribution[0] = 0;
            contribution[1] = 0;
        }
        std::cout << "t: " << t << " eta: " << sqrt(totalError[0])
                  << " phi: " << sqrt(totalError[1]) << std::endl;
    }

    void computeEnergy() {
        Integration::ElementIntegral<DIM> elIntegral;
        Integration::FaceIntegral<DIM> faIntegral;
        LinearAlgebra::MiddleSizeVector totalEnergy(1), contribution(1);
        for (Base::Face* face : meshes_[0]->getFacesList()) {
            contribution = faIntegral.integrate(face, &faceEnergy);
            totalEnergy += contribution;
            contribution[0] = 0;
        }
        for (Base::Element* element : meshes_[0]->getElementsList()) {
            contribution = elIntegral.integrate(element, &elementEnergy);
            totalEnergy += contribution;
            contribution[0] = 0;
        }
        std::cout << "Energy: " << totalEnergy[0] << std::endl;
    }

    void setUpSwap(Mat swap) {
        PetscInt n;
        MatGetSize(swap, &n, PETSC_NULL);
        if (n > 0) {
            for (std::size_t i = 0; i < static_cast<std::size_t>(n - 1);
                 i += 2) {
                MatSetValue(swap, i, i + 1, 1., INSERT_VALUES);
                MatSetValue(swap, i + 1, i, 1., INSERT_VALUES);
            }
        }
        MatAssemblyBegin(swap, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(swap, MAT_FINAL_ASSEMBLY);
    }

    bool solve() {
        std::cout.precision(14);
        doAllElementIntegration();
        doAllFaceIntegration();

        Utilities::GlobalIndexing indexing(meshes_[0]);
        Utilities::GlobalPetscVector eta(indexing), phi(indexing);
        Utilities::GlobalPetscMatrix M(indexing, -1, 0), S(indexing, 0, -1);
        Vec phiS, phiOther, etaActually, interiorRHS, surfaceRHS, surfaceExtra;
        Mat surfaceMass, interiorStifness, surfaceStifness, mixStifness,
            backStiffness, swapSurfaceVars;

        IS isSurface, isRest;
        getSurfaceIS(S, &isSurface, &isRest);

        std::ofstream outFile("output.dat");
        Output::TecplotDiscontinuousSolutionWriter<DIM> writeFunc(
            outFile, "test", "01", "eta_num,phi_num,eta_exact,phi_exact");
        Output::VTKTimeDependentWriter<DIM> paraWrite("output", meshes_[0]);

        // deal with initial conditions
        double g(1.), dt(M_PI / 16);
        eta.constructFromTimeIntegrationVector(0, 0);
        phi.constructFromTimeIntegrationVector(0, 1);

        std::size_t numberOfSnapshots(65);  // placeholder parameters
        std::size_t numberOfTimeSteps(1);

        MatGetSubMatrix(M, isSurface, isSurface, MAT_INITIAL_MATRIX,
                        &surfaceMass);
        MatDuplicate(surfaceMass, MAT_DO_NOT_COPY_VALUES, &swapSurfaceVars);
        setUpSwap(swapSurfaceVars);
        MatGetSubMatrix(S, isRest, isRest, MAT_INITIAL_MATRIX,
                        &interiorStifness);
        MatGetSubMatrix(S, isSurface, isSurface, MAT_INITIAL_MATRIX,
                        &surfaceStifness);
        MatGetSubMatrix(S, isRest, isSurface, MAT_INITIAL_MATRIX, &mixStifness);
        MatGetSubMatrix(S, isSurface, isRest, MAT_INITIAL_MATRIX,
                        &backStiffness);

        KSP interior, surface;
        KSPCreate(PETSC_COMM_WORLD, &interior);
        KSPSetOperators(interior, interiorStifness, interiorStifness);
        KSPSetInitialGuessNonzero(interior, PETSC_TRUE);
        KSPSetFromOptions(interior);
        KSPCreate(PETSC_COMM_WORLD, &surface);
        KSPSetOperators(surface, surfaceMass, surfaceMass);
        KSPSetFromOptions(surface);
        KSPConvergedReason conferge;
        int iterations;

        VecGetSubVector(phi, isSurface, &phiS);
        VecGetSubVector(phi, isRest, &phiOther);
        VecGetSubVector(eta, isSurface, &etaActually);
        VecDuplicate(phiS, &surfaceRHS);
        VecDuplicate(phiS, &surfaceExtra);
        VecDuplicate(phiOther, &interiorRHS);
        VecCopy(phiS, surfaceRHS);
        KSPSolve(surface, surfaceRHS, phiS);
        KSPGetConvergedReason(surface, &conferge);
        KSPGetIterationNumber(surface, &iterations);
        std::cout
            << "Finalizing interpolation (1): KSP solver ended because of "
            << KSPConvergedReasons[conferge] << " in " << iterations
            << " iterations." << std::endl;
        VecCopy(etaActually, surfaceRHS);
        KSPSolve(surface, surfaceRHS, etaActually);
        KSPGetConvergedReason(surface, &conferge);
        KSPGetIterationNumber(surface, &iterations);
        std::cout
            << "Finalizing interpolation (2): KSP solver ended because of "
            << KSPConvergedReasons[conferge] << " in " << iterations
            << " iterations." << std::endl;
        VecRestoreSubVector(eta, isSurface, &etaActually);
        VecRestoreSubVector(phi, isSurface, &phiS);
        VecRestoreSubVector(phi, isRest, &phiOther);

        eta.writeTimeIntegrationVector(0, 0);
        phi.writeTimeIntegrationVector(0, 1);

        writeFunc.write(meshes_[0], "solution at time 0", false, this);
        std::function<double(Base::Element*,
                             const Geometry::PointReference<DIM>&, std::size_t)>
            funcEta = [](Base::Element* element,
                         const Geometry::PointReference<DIM>& p,
                         std::size_t timeIntegrationVectorId) -> double {
            return element->getSolution(timeIntegrationVectorId, p)[0];
        };
        std::function<double(Base::Element*,
                             const Geometry::PointReference<DIM>&, std::size_t)>
            funcPhi = [](Base::Element* element,
                         const Geometry::PointReference<DIM>& p,
                         std::size_t timeIntegrationVectorId) -> double {
            return element->getSolution(timeIntegrationVectorId, p)[1];
        };
        paraWrite.write(funcEta, "eta", t, 0);
        paraWrite.write(funcPhi, "phi", t, 0);
        for (std::size_t i = 1; i < numberOfSnapshots; ++i) {
            VecGetSubVector(phi, isSurface, &phiS);
            VecGetSubVector(phi, isRest, &phiOther);
            VecGetSubVector(eta, isSurface, &etaActually);
            for (std::size_t j = 0; j < numberOfTimeSteps; ++j) {
                t += dt;
                MatMult(swapSurfaceVars, etaActually, surfaceRHS);
                VecAXPY(phiS, -g * dt / 2,
                        surfaceRHS);  //(can in principle be combined with final
                                      //step, but this makes output easier)
                MatMult(mixStifness, phiS, interiorRHS);
                VecScale(interiorRHS, -1);
                KSPSolve(interior, interiorRHS, phiOther);
                KSPGetConvergedReason(interior, &conferge);
                KSPGetIterationNumber(interior, &iterations);
                std::cout << "Laplace problem: KSP solver ended because of "
                          << KSPConvergedReasons[conferge] << " in "
                          << iterations << " iterations." << std::endl;

                MatMult(backStiffness, phiOther, surfaceExtra);
                MatMultAdd(surfaceStifness, phiS, surfaceExtra, surfaceExtra);
                MatMult(swapSurfaceVars, surfaceExtra, surfaceRHS);
                VecScale(surfaceRHS, dt);
                MatMultAdd(surfaceMass, etaActually, surfaceRHS, surfaceRHS);
                KSPSolve(surface, surfaceRHS, etaActually);
                KSPGetConvergedReason(surface, &conferge);
                KSPGetIterationNumber(surface, &iterations);
                std::cout << "Updating \\eta: KSP solver ended because of "
                          << KSPConvergedReasons[conferge] << " in "
                          << iterations << " iterations." << std::endl;

                MatMult(swapSurfaceVars, etaActually, surfaceRHS);
                VecAXPY(phiS, -g * dt / 2, surfaceRHS);
            }
            VecRestoreSubVector(phi, isSurface, &phiS);
            VecRestoreSubVector(phi, isRest, &phiOther);
            VecRestoreSubVector(eta, isSurface, &etaActually);
            eta.writeTimeIntegrationVector(0, 0);
            phi.writeTimeIntegrationVector(0, 1);
            writeFunc.write(meshes_[0], "solution at time t", false, this);
            paraWrite.write(funcEta, "eta", t, 0);
            paraWrite.write(funcPhi, "phi", t, 0);
            printError();
            computeEnergy();
        }
        return true;
    }

   private:
    // polynomial order of the approximation
    std::size_t p_;
    static double t;
};

double DGWave::t = 0;
auto& p =
    Base::register_argument<std::size_t>('p', "poly", "Polynomial order", true);

int main(int argc, char** argv) {
    Base::parse_options(argc, argv);
    DGWave test(p.getValue());
    test.initialise();
    test.solve();
    return 0;
}
