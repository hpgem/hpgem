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

#include "HpgemAPINonLinearSteadyState.h"
#include "Utilities/GlobalVector.h"

#if defined(HPGEM_USE_SUNDIALS)
#include "kinsol/kinsol.h"
#include "kinsol/kinsol_spgmr.h"
#include "nvector/nvector_serial.h"
#endif

#include "Logger.h"
#include <iostream>  //remove this after this API is official
#include <fstream>   //remove this after this API is official

namespace Base {
// todo: Possibly add restart from data file to create a p-multigrid solution
// strategy todo: Add function that can add setup functions to the solve function

// note: compute error is always on, compute both faces is always on
template <std::size_t DIM>
HpgemAPINonLinearSteadyState<DIM>::HpgemAPINonLinearSteadyState(
    const std::size_t numberOfVariables, const std::size_t polynomialOrder,
    const bool computeBothFaces)
    : HpgemAPISimplified<DIM>(numberOfVariables, polynomialOrder, 4, 0,
                              computeBothFaces),
      sourceElementVectorID_(0),
      sourceFaceVectorID_(0) {}

template <std::size_t DIM>
void HpgemAPINonLinearSteadyState<DIM>::computeJacobian() {
    // Compute new Jacobian element matrices for J(u)
    for (Base::Element *ptrElement : this->meshes_[0]->getElementsList()) {
        // Get the solution coefficients required for the calculation
        LinearAlgebra::MiddleSizeVector &solutionCoefficients(
            ptrElement->getTimeIntegrationVector(this->solutionVectorId_));
        // Calculate and store the matrix
        ptrElement->setElementMatrix(
            computeJacobianAtElement(ptrElement, solutionCoefficients, 0),
            jacobianElementMatrixID_);
    }

    // Compute the local Jacobian Matrix: Face integrals
    for (Base::Face *ptrFace : this->meshes_[0]->getFacesList()) {
        // grab the coefficients
        LinearAlgebra::MiddleSizeVector solutionCoefficientsLeft(
            ptrFace->getPtrElementLeft()->getTimeIntegrationVector(
                this->solutionVectorId_));
        LinearAlgebra::MiddleSizeVector solutionCoefficientsRight(
            ptrFace->getPtrElementRight()->getTimeIntegrationVector(
                this->solutionVectorId_));

        // Calculate face Matrix
        int numberOfDOFLeft =
            ptrFace->getPtrElement(Base::Side::LEFT)
                ->getNumberOfBasisFunctions() *
            ptrFace->getPtrElement(Base::Side::LEFT)->getNumberOfUnknowns();
        int numberOfDOFRight =
            ptrFace->getPtrElement(Base::Side::RIGHT)
                ->getNumberOfBasisFunctions() *
            ptrFace->getPtrElement(Base::Side::RIGHT)->getNumberOfUnknowns();
        Base::FaceMatrix faceMatrix(numberOfDOFLeft, numberOfDOFRight);
        // Left element. Side indexes are: elementSide, derivativeSide
        faceMatrix.setElementMatrix(
            computeJacobianAtFace(ptrFace, solutionCoefficientsLeft,
                                  solutionCoefficientsRight, 0,
                                  Base::Side::LEFT, Base::Side::LEFT),
            Base::Side::LEFT, Base::Side::LEFT);
        faceMatrix.setElementMatrix(
            computeJacobianAtFace(ptrFace, solutionCoefficientsLeft,
                                  solutionCoefficientsRight, 0,
                                  Base::Side::LEFT, Base::Side::RIGHT),
            Base::Side::LEFT, Base::Side::RIGHT);
        // Right element
        faceMatrix.setElementMatrix(
            computeJacobianAtFace(ptrFace, solutionCoefficientsLeft,
                                  solutionCoefficientsRight, 0,
                                  Base::Side::RIGHT, Base::Side::LEFT),
            Base::Side::RIGHT, Base::Side::LEFT);
        faceMatrix.setElementMatrix(
            computeJacobianAtFace(ptrFace, solutionCoefficientsLeft,
                                  solutionCoefficientsRight, 0,
                                  Base::Side::RIGHT, Base::Side::RIGHT),
            Base::Side::RIGHT, Base::Side::RIGHT);

        // Store face Matrix
        ptrFace->setFaceMatrix(faceMatrix, jacobianFaceMatrixID_);  // time=0
    }
}

#if defined(HPGEM_USE_SUNDIALS)
// static int KINSol function required function to compute the right hand side.
// The user_data is recast into the API class and a corresponding memeber
// function is then called to compute the RHS
template <std::size_t DIM>
int HpgemAPINonLinearSteadyState<DIM>::func(N_Vector u, N_Vector fval,
                                            void *user_data) {
    return static_cast<HpgemAPINonLinearSteadyState<DIM> *>(user_data)
        ->computeRHS(u, fval);
}

template <std::size_t DIM>
int HpgemAPINonLinearSteadyState<DIM>::jtimes(N_Vector v, N_Vector Jv,
                                              N_Vector u, booleantype *new_u,
                                              void *user_data) {
    return static_cast<HpgemAPINonLinearSteadyState<DIM> *>(user_data)
        ->computeJacTimesVector(v, Jv, u, new_u);
}

// Member function computing the rhs
template <std::size_t DIM>
int HpgemAPINonLinearSteadyState<DIM>::computeRHS(N_Vector u, N_Vector fval) {
    // Pass the solution vector u to the GlobalSundialsVector
    globalVector_->setVector(u);
    // The GlobalSundialsVector will then put the data in the hpGEM data
    // structure
    globalVector_->writeTimeIntegrationVector(this->solutionVectorId_);

    // Compute RHS with the new data
    // solutionCoefficients are stored at this->solutionVectorId_ (=0 by
    // default) and the rhs at 1, the time = 0
    this->computeRightHandSide(this->solutionVectorId_, this->solutionRHS_, 0);
    // Set the correct N_Vector pointer
    globalVector_->setVector(fval);
    // Pass the new solution back to KINsol
    globalVector_->constructFromTimeIntegrationVector(
        this->solutionRHS_);  // 1 corresponds to the resultId vector

    // Write intermediate solutions
    if (doOutputIntermediateSolutions_) {
        // Update the nstep number
        // NOTE: step is not a real step, just the number of times the function
        // is called
        step_++;
        logger(INFO, "Writing intermediate solution: %", step_);
        tecplotWriter_->write(this->meshes_[0], this->solutionTitle_, false,
                              this, step_);
        this->VTKWrite(*(this->VTKWriter_), step_, this->solutionVectorId_);
    }
    return 0;
}

// todo: If required, make this function work for conforming
template <std::size_t DIM>
int HpgemAPINonLinearSteadyState<DIM>::computeJacTimesVector(
    N_Vector v, N_Vector Jv, N_Vector u, booleantype *new_u) {
    // Pass the solution vector u to the GlobalSundialsVector
    globalVector_->setVector(u);
    // The GlobalSundialsVector will then put the data in the hpGEM data
    // structure
    globalVector_->writeTimeIntegrationVector(this->solutionVectorId_);

    // If a new solution exists, new Jacobian values need to be calculated
    if (*new_u == true) {
        logger(INFO, "Computing new Jacobian Matrix.");
        computeJacobian();
        // Set new_u to false since we have computed the new matrices
        *new_u = false;

        // DEBUG
        // todo: Move this to a function, that is there just for debug purposes
        // Create a very awesome big matrix
        LinearAlgebra::MiddleSizeMatrix matrix;
        std::size_t totalNumberOfDOF = NV_LENGTH_S(u);
        std::cout << "totalNumberOfDOF: " << totalNumberOfDOF << std::endl;
        LinearAlgebra::MiddleSizeMatrix jacobianMatrixFull(totalNumberOfDOF,
                                                           totalNumberOfDOF);
        std::size_t numberOfDOFLocal;
        std::size_t numberOfDOFNonLocal;
        int startLocal;     // startposition of the local element
        int startNonLocal;  // start position of the non local element
        std::vector<int> startPositionsOfElementsInTheVector =
            globalVector_->getStartPositionsOfElementsInTheVector();

        // Fill the very big awesome matrix
        for (Base::Element *ptrElement : this->meshes_[0]->getElementsList()) {
            // Add the element matrix
            numberOfDOFLocal = ptrElement->getNrOfBasisFunctions() *
                               ptrElement->getNrOfUnknowns();
            matrix = ptrElement->getElementMatrix(jacobianElementMatrixID_);
            startLocal =
                startPositionsOfElementsInTheVector[ptrElement->getID()];
            for (std::size_t i = 0; i < numberOfDOFLocal; i++) {
                for (std::size_t j = 0; j < numberOfDOFLocal; j++) {
                    jacobianMatrixFull(i + startLocal, j + startLocal) +=
                        matrix(i, j);
                }
            }

            // Face matrices
            for (const Base::Face *ptrFace : ptrElement->getFacesList()) {
                Base::Side elementSide;
                Base::Side neighbourElementSide;
                // Check if current element is the left or the right element
                if (ptrFace->getPtrElementLeft() == ptrElement) {
                    elementSide = Base::Side::LEFT;
                    neighbourElementSide = Base::Side::RIGHT;
                } else {
                    elementSide = Base::Side::RIGHT;
                    neighbourElementSide = Base::Side::LEFT;
                }
                numberOfDOFLocal = ptrFace->getPtrElement(elementSide)
                                       ->getNrOfBasisFunctions() *
                                   ptrElement->getNrOfUnknowns();
                numberOfDOFNonLocal =
                    ptrFace->getPtrElement(neighbourElementSide)
                        ->getNrOfBasisFunctions() *
                    ptrElement->getNrOfUnknowns();
                startLocal = startPositionsOfElementsInTheVector
                    [ptrFace->getPtrElement(elementSide)->getID()];
                startNonLocal = startPositionsOfElementsInTheVector
                    [ptrFace->getPtrElement(neighbourElementSide)->getID()];

                // Compute local face contribution
                matrix = ptrFace->getFaceMatrix(this->jacobianFaceMatrixID_)
                             .getElementMatrix(elementSide, elementSide);
                for (std::size_t i = 0; i < numberOfDOFLocal; i++) {
                    for (std::size_t j = 0; j < numberOfDOFLocal; j++) {
                        jacobianMatrixFull(i + startLocal, j + startLocal) +=
                            matrix(i, j);
                    }
                }

                // Compute non local face contribution
                matrix =
                    ptrFace->getFaceMatrix(jacobianFaceMatrixID_)
                        .getElementMatrix(elementSide, neighbourElementSide);
                for (std::size_t i = 0; i < numberOfDOFLocal; i++) {
                    for (std::size_t j = 0; j < numberOfDOFNonLocal; j++) {
                        jacobianMatrixFull(i + startLocal, j + startNonLocal) +=
                            matrix(i, j);
                    }
                }
            }
        }

        // Write matrix to file
        remove("matrix");
        std::ofstream myfile("matrix");
        if (myfile.is_open()) {
            for (std::size_t i = 0; i < totalNumberOfDOF; i++) {
                // print a line
                for (std::size_t j = 0; j < totalNumberOfDOF - 1; j++) {
                    myfile << jacobianMatrixFull(i, j) << ",";
                }
                myfile << jacobianMatrixFull(i, totalNumberOfDOF - 1);
                myfile << "\n";
            }
            myfile.close();
        } else
            cout << "Unable to open file";
        std::exit(-1);
        // DEBUG
    }

    // Compute J(u) times v
    // Put v in the hpGEM structure
    globalVector_->setVector(v);
    globalVector_->writeTimeIntegrationVector(this->VectorVID_);
    // logger(INFO,"Computing new Jacobian multiplication.");
    // Compute J(u) times v for every element and place it in the hpGEM
    // structure
    for (Base::Element *ptrElement : this->meshes_[0]->getElementsList()) {
        LinearAlgebra::MiddleSizeVector jTimesV(
            ptrElement->getNumberOfBasisFunctions() *
            ptrElement->getNumberOfUnknowns());
        LinearAlgebra::MiddleSizeMatrix matrix;
        LinearAlgebra::MiddleSizeVector vector;

        // Element matrix
        matrix = ptrElement->getElementMatrix(jacobianElementMatrixID_);
        vector = ptrElement->getTimeIntegrationVector(this->VectorVID_);
        jTimesV += matrix * vector;

        // Face matrices
        for (const Base::Face *ptrFace : ptrElement->getFacesList()) {
            Base::Side elementSide;
            Base::Side neighbourElementSide;
            // Check if current element is the left or the right element
            if (ptrFace->getPtrElementLeft() == ptrElement) {
                elementSide = Base::Side::LEFT;
                neighbourElementSide = Base::Side::RIGHT;
            } else {
                elementSide = Base::Side::RIGHT;
                neighbourElementSide = Base::Side::LEFT;
            }

            // DIT IS EEN BUG?
            // Compute local face contribution
            matrix = ptrFace->getFaceMatrix(this->jacobianFaceMatrixID_)
                         .getElementMatrix(elementSide, elementSide);
            vector = ptrFace->getPtrElement(elementSide)
                         ->getTimeIntegrationVector(this->VectorVID_);
            jTimesV += matrix * vector;

            // Compute non local face contribution
            matrix = ptrFace->getFaceMatrix(jacobianFaceMatrixID_)
                         .getElementMatrix(elementSide, neighbourElementSide);
            vector = ptrFace->getPtrElement(neighbourElementSide)
                         ->getTimeIntegrationVector(this->VectorVID_);
            jTimesV += matrix * vector;
        }

        // Set resulting vector
        ptrElement->setTimeIntegrationVector(jTimesvecID_, jTimesV);
    }

    // Write the data from hpGEM to KINSOL
    globalVector_->setVector(Jv);
    globalVector_->constructFromTimeIntegrationVector(jTimesvecID_);
    return 0;
}

#endif

template <std::size_t DIM>
bool HpgemAPINonLinearSteadyState<DIM>::solve(bool doComputeInitialCondition,
                                              bool doComputeError,
                                              bool doUseJacobian) {
#if defined(HPGEM_USE_SUNDIALS)
    int flag;
    int maxl;
    int globalStrategy = KIN_LINESEARCH;  // For now: nothing special, in future
                                          // maybe KIN_LINESEARCH

    // Create output files for Paraview.
    std::string outputFileNameVTK = this->outputFileName_;

    // todo: Maybe add scaling

    this->registerVTKWriteFunctions();
    Output::VTKTimeDependentWriter<DIM> VTKWriter(outputFileNameVTK,
                                                  this->meshes_[0]);
    VTKWriter_ = &VTKWriter;

    // Create output files for Tecplot.
#ifdef HPGEM_USE_MPI
    std::string outputFileName =
        this->outputFileName_ + "." +
        std::to_string(Base::MPIContainer::Instance().getProcessorID());
#else
    std::string outputFileName = this->outputFileName_;
#endif
    std::string outputFileNameTecplot = outputFileName + ".dat";
    std::string dimensionsToWrite = "";
    for (std::size_t i = 0; i < this->configData_->dimension_; i++) {
        dimensionsToWrite = dimensionsToWrite + std::to_string(i);
    }

    std::string variableString = this->variableNames_[0];
    for (std::size_t iV = 1; iV < this->variableNames_.size(); iV++) {
        variableString = variableString + "," + this->variableNames_[iV];
    }

    std::ofstream outputFile(outputFileNameTecplot);
    Output::TecplotDiscontinuousSolutionWriter<DIM> tecplotWriter(
        outputFile, this->internalFileTitle_, dimensionsToWrite,
        variableString);
    tecplotWriter_ = &tecplotWriter;

    // Create KINsol memory block
    void *kmem = KINCreate();
    logger.assert_always(kmem != nullptr, "Failed creating KINmem.");

    // Initialise Global Solution Vector and template N_Vector u
    Utilities::GlobalSundialsVector globalVector(HpgemAPIBase<DIM>::meshes_[0],
                                                 sourceElementVectorID_,
                                                 sourceFaceVectorID_);
    globalVector_ = &globalVector;
    std::size_t numberOfDOF = globalVector_->getTotalNumberOfDOF();
    N_Vector u = N_VNew_Serial(numberOfDOF);

    // Create the scale vector. For now this is just all set to 1
    N_Vector scale = N_VNew_Serial(numberOfDOF);
    double *ptrScale = NV_DATA_S(scale);
    for (std::size_t i = 0; i < numberOfDOF; i++) {
        ptrScale[i] = 1.0;
    }

    // If an initial condition is supplied, compute it and initialise u with
    // this vector, else the vector is initialised as 0
    if (doComputeInitialCondition == true) {
        // Let hpGEM compute the initial coefficients
        this->setInitialSolution(this->solutionVectorId_, 0, 0);

        // obtain globalvector for initial conditions
        globalVector_->setVector(u);
        globalVector_->constructFromTimeIntegrationVector(0);
    }

    // Give userdata to KINSol, as a void pointer
    void *userdata = static_cast<void *>(this);
    flag = KINSetUserData(kmem, userdata);
    logger.assert_always(flag >= 0,
                         "User data not correctly passed to KINSol.");

    // Initialise kmem
    flag = KINInit(kmem, func, u);
    logger.assert_always(flag >= 0, "Initialisation failed with flag %.", flag);

    // set number of iterations
    flag = KINSetNumMaxIters(kmem, 8000);

    // Set FTOL and STOL
    flag = KINSetFuncNormTol(kmem, 1e-7);
    flag = KINSetScaledStepTol(kmem, 1e-13);

    // Set type of solver
    // todo: Make a switch between different types of solver
    maxl = 15;
    flag = KINSpgmr(kmem, maxl);
    logger.assert_debug(flag >= 0, "Failed to attach GMRS solver.");

    // Set jacobian times vector function, if specified
    if (doUseJacobian == true) {
        flag = KINSpilsSetJacTimesVecFn(kmem, jtimes);
        logger.assert_debug(flag >= 0,
                            "Failed to attach jtimes function to KINSol.");
    }

    // Set additional parameters

    this->tasksBeforeSolving();

    // Start solving
    flag = KINSol(kmem, u, globalStrategy, scale, scale);
    if ((flag == 0) || (flag == 1)) {
        logger(INFO, "Solution obtained");
    } else {
        logger(INFO, "No solution found0");
        // logger.assert_debug(flag >= 0, "Failed to solve the problem with flag
        // %.",flag);
    }

    // Write final solution to VTK
    step_++;
    tecplotWriter.write(this->meshes_[0], this->solutionTitle_, false, this,
                        step_);
    this->VTKWrite(VTKWriter, step_, this->solutionVectorId_);

    // DEBUG
    N_Vector temp = N_VNew_Serial(numberOfDOF);
    globalVector_->setVector(temp);
    // Pass the new solution back to KINsol
    globalVector_->constructFromTimeIntegrationVector(1);  // 1 corresponds to
                                                           // the resultId vector
    double *tempData = NV_DATA_S(temp);
    double max = 0;
    for (std::size_t i = 0; i < numberOfDOF; i++) {
        if (std::abs(tempData[i]) > max) {
            max = tempData[i];
        }
    }
    std::cout << "Max RHS value: " << max << std::endl;
    // DEBUG

    // Compute the energy norm of the error
    if (doComputeError) {
        LinearAlgebra::MiddleSizeVector::type totalError =
            this->computeTotalError(this->solutionVectorId_, 0);
        logger(INFO, "Total error: %.", totalError);
        LinearAlgebra::MiddleSizeVector maxError =
            this->computeMaxError(this->solutionVectorId_, 0);
        logger.assert_debug(
            maxError.size() == this->configData_->numberOfUnknowns_,
            "Size of maxError (%) not equal to the number of variables (%)",
            maxError.size(), this->configData_->numberOfUnknowns_);
        for (std::size_t iV = 0; iV < this->configData_->numberOfUnknowns_;
             iV++) {
            logger(INFO, "Maximum error %: %", this->variableNames_[iV],
                   maxError(iV));
        }
    }

    // Free memory
    KINFree(&kmem);
    N_VDestroy_Serial(u);
    N_VDestroy_Serial(scale);

    return false;
#endif
    logger(ERROR,
           "Sundials is required to solve the non linear steady state problem "
           "using this function.");
    return false;
}
}  // namespace Base
