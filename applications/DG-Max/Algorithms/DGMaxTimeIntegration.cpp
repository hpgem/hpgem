/*
This file forms part of hpGEM. This package has been developed over a number of
years by various people at the University of Twente and a full list of
contributors can be found at http://hpgem.org/about-the-code/team

This code is distributed using BSD 3-Clause License. A copy of which can found
below.


Copyright (c) 2018, Univesity of Twenete
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "DGMaxTimeIntegration.h"

#include "petscksp.h"

#include "Utilities/GlobalMatrix.h"
#include "Utilities/GlobalVector.h"

using namespace hpgem;

template <std::size_t DIM>
DGMaxTimeIntegration<DIM>::DGMaxTimeIntegration(
    Base::MeshManipulator<DIM>& mesh, std::size_t order, double stab)
    : mesh_(mesh), discretization(order, stab), snapshotTime(nullptr) {
    discretization.initializeBasisFunctions(mesh_);
}

template <std::size_t DIM>
DGMaxTimeIntegration<DIM>::~DGMaxTimeIntegration() {
    if (snapshotTime != nullptr) {
        delete[] snapshotTime;
    }
}

template <std::size_t DIM>
void DGMaxTimeIntegration<DIM>::solve(
    const SeparableTimeIntegrationProblem<DIM>& input,
    TimeIntegrationParameters<DIM> parameters) {
    // TODO: It might be better to separate some parts of the CO2 and CO4
    // algorithm
    using namespace std::placeholders;

    if (snapshotTime != nullptr) {
        delete[] snapshotTime;
        snapshotTime = nullptr;
    }
    snapshotTime = new double[parameters.numberOfSnapshots()];
    for (Base::Element* element :
         mesh_.getElementsList(Base::IteratorType::GLOBAL)) {
        element->setNumberOfTimeIntegrationVectors(
            parameters.numberOfSnapshots());
    }

    PetscErrorCode error;

    std::cout << "doing a time dependent simulation" << std::endl;
    std::map<std::size_t, typename DGMaxDiscretization<DIM>::InputFunction>
        elementVectors;

    elementVectors[DGMaxDiscretizationBase::ELEMENT_VECTOR_ID] =
        std::bind(&SeparableTimeIntegrationProblem<DIM>::sourceTermRef,
                  std::ref(input), _1);
    elementVectors[DGMaxDiscretizationBase::INITIAL_CONDITION_VECTOR_ID] =
        std::bind(&TimeIntegrationProblem<DIM>::initialCondition,
                  std::ref(input), _1);
    elementVectors
        [DGMaxDiscretizationBase::INITIAL_CONDITION_DERIVATIVE_VECTOR_ID] =
            std::bind(&TimeIntegrationProblem<DIM>::initialConditionDerivative,
                      std::ref(input), _1);

    // Using this is not supported
    const double dispersionFrequency =
        std::numeric_limits<double>::signaling_NaN();
    discretization.setMatrixHandling(DGMaxDiscretizationBase::INVERT);
    discretization.computeElementIntegrals(mesh_, elementVectors,
                                           dispersionFrequency);

    std::map<std::size_t, typename DGMaxDiscretization<DIM>::FaceInputFunction>
        faceVectors;
    faceVectors[DGMaxDiscretizationBase::FACE_VECTOR_ID] =
        std::bind(&SeparableTimeIntegrationProblem<DIM>::boundaryConditionRef,
                  std::ref(input), _1);

    discretization.computeFaceIntegrals(mesh_, faceVectors,
                                        dispersionFrequency);
    //    MHasToBeInverted_ = true;
    //    assembler->fillMatrices(this);

    Utilities::GlobalIndexing indexing(
        &mesh_, Utilities::GlobalIndexing::BLOCKED_PROCESSOR);
    Utilities::GlobalPetscMatrix massMatrix(
        indexing, DGMaxDiscretization<DIM>::MASS_MATRIX_ID, -1),
        stiffnessMatrix(indexing, DGMaxDiscretization<DIM>::STIFFNESS_MATRIX_ID,
                        DGMaxDiscretization<DIM>::FACE_STIFFNESS_MATRIX_ID);
    std::cout << "GlobalPetscMatrix initialised" << std::endl;
    Utilities::GlobalPetscVector resultVector(
        indexing, DGMaxDiscretization<DIM>::INITIAL_CONDITION_VECTOR_ID, -1),
        derivative(
            indexing,
            DGMaxDiscretization<DIM>::INITIAL_CONDITION_DERIVATIVE_VECTOR_ID,
            -1),
        rhsBoundary(indexing, -1, DGMaxDiscretization<DIM>::FACE_VECTOR_ID),
        rhsSource(indexing, DGMaxDiscretization<DIM>::ELEMENT_VECTOR_ID, -1);
    std::cout << "GlobalPetscVector initialised" << std::endl;
    resultVector.assemble();
    std::cout << "resultVector assembled" << std::endl;
    rhsBoundary.assemble();
    rhsSource.assemble();
    std::cout << "RHS_ assembled" << std::endl;
    derivative.assemble();
    std::cout << "derivative_ assembled" << std::endl;

    Vec dummy, dummy2, dummy3;
    error = VecDuplicate(derivative, &dummy);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    if (parameters.method == CO4) {
        error = VecDuplicate(derivative, &dummy3);
        CHKERRABORT(PETSC_COMM_WORLD, error);
    }
    error = VecDuplicate(derivative, &dummy2);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = MatMult(massMatrix, derivative, dummy);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = VecCopy(dummy, derivative);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = MatMult(massMatrix, resultVector, dummy);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = VecCopy(dummy, resultVector);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    double tau = parameters.timeStepSize;

    double t = 0;
    snapshotTime[0] = t;

    double scale0, scale1;
    if (parameters.method == CO2) {
        scale0 = (1 - input.conductivity() * tau / 2);
        scale1 = 1 / (1 + input.conductivity() * tau / 2);
    }

    LinearAlgebra::SmallVector<6> alpha, beta, alpha_sum, beta_sum,
        scale0vector, scale1vector;
    double sourceTime;
    if (parameters.method == CO4) {
        getCoeffCO4(alpha, beta, alpha_sum, beta_sum, scale0vector,
                    scale1vector);
    }
    int snapshotIndex = 0;

    std::cout << tau << " " << parameters.numberOfSteps << std::endl;

    for (int i = 0; i < parameters.numberOfSteps; ++i) {
        if (i % parameters.snapshotStride == 0) {
            logger(DEBUG, "Writing snapshot at timestep %, (t=%)", i, t);
            resultVector.writeTimeIntegrationVector(snapshotIndex);
            snapshotTime[snapshotIndex] = t;
            snapshotIndex++;
        }
        if (parameters.method == CO2) {
            // leap-frog sceme (Yee)
            error = VecAXPY(resultVector, tau / 2, derivative);
            CHKERRABORT(PETSC_COMM_WORLD, error);
            error = MatMult(stiffnessMatrix, resultVector, dummy);
            CHKERRABORT(
                PETSC_COMM_WORLD,
                error);  // starting here:dummy contaings partial update of x

            // compute dummy2, a time scaled version of the RHS contains
            // time-scaled version of RHS
            error = VecCopy(rhsBoundary, dummy2);
            CHKERRABORT(PETSC_COMM_WORLD, error);
            error = VecScale(dummy2, input.timeScalingBoundary(t));
            CHKERRABORT(PETSC_COMM_WORLD, error);
            error = VecAXPY(dummy2, input.timeScalingSource(t), rhsSource);
            CHKERRABORT(PETSC_COMM_WORLD, error);

            CHKERRABORT(PETSC_COMM_WORLD, error);
            error = VecScale(dummy2, 0.5 * tau);
            CHKERRABORT(PETSC_COMM_WORLD, error);
            error = VecAYPX(dummy, -tau, dummy2);
            CHKERRABORT(PETSC_COMM_WORLD, error);
            // Stage 2
            error = VecCopy(rhsBoundary, dummy2);
            CHKERRABORT(PETSC_COMM_WORLD, error);
            error = VecScale(dummy2, input.timeScalingBoundary(t + tau));
            CHKERRABORT(PETSC_COMM_WORLD, error);
            error =
                VecAXPY(dummy2, input.timeScalingSource(t + tau), rhsSource);
            CHKERRABORT(PETSC_COMM_WORLD, error);

            error = VecAXPY(dummy, 0.5 * tau, dummy2);
            CHKERRABORT(PETSC_COMM_WORLD, error);
            error = VecScale(dummy, scale1);
            CHKERRABORT(PETSC_COMM_WORLD, error);
            error = MatMult(massMatrix, dummy, dummy2);
            CHKERRABORT(
                PETSC_COMM_WORLD,
                error);  // starting here: dummy2 contiains partial update of x
            error = VecAYPX(derivative, scale0 * scale1, dummy2);
            CHKERRABORT(PETSC_COMM_WORLD, error);
            error = VecAXPY(resultVector, tau / 2, derivative);
            CHKERRABORT(PETSC_COMM_WORLD, error);
            t = t + tau;
        }

        if (parameters.method == CO4) {
            for (unsigned int k = 1; k < 6; ++k) {
                error = VecAXPY(resultVector, tau * (alpha[k - 1] + beta[k]),
                                derivative);
                CHKERRABORT(PETSC_COMM_WORLD, error);

                error = MatMult(stiffnessMatrix, resultVector, dummy);
                CHKERRABORT(PETSC_COMM_WORLD, error);

                sourceTime = t + (alpha_sum[k - 1] + beta_sum[k - 1]) * tau;
                double eta = input.timeScalingBoundary(sourceTime);
                // error = VecCopy(RHS_, dummy2);
                CHKERRABORT(PETSC_COMM_WORLD, error);
                error = VecScale(dummy2, beta[k] * tau * eta);
                CHKERRABORT(PETSC_COMM_WORLD, error);
                VecAYPX(dummy, -(beta[k] + alpha[k]) * tau, dummy2);
                CHKERRABORT(PETSC_COMM_WORLD, error);

                sourceTime = t + (alpha_sum[k] + beta_sum[k]) * tau;
                eta = input.timeScalingBoundary(sourceTime);
                // error = VecCopy(RHS_, dummy2);
                CHKERRABORT(PETSC_COMM_WORLD, error);
                error = VecScale(dummy2, alpha[k] * tau * eta);
                CHKERRABORT(PETSC_COMM_WORLD, error);
                error = VecAXPY(dummy, 1.0, dummy2);
                CHKERRABORT(PETSC_COMM_WORLD, error);

                error = VecScale(dummy, scale1vector[k]);
                CHKERRABORT(PETSC_COMM_WORLD, error);
                error = MatMult(massMatrix, dummy, dummy3);
                CHKERRABORT(PETSC_COMM_WORLD, error);
                VecAYPX(derivative, scale0vector[k] * scale1vector[k], dummy3);
                CHKERRABORT(PETSC_COMM_WORLD, error);
            }
            error = VecAXPY(resultVector, tau * alpha[5], derivative);
            CHKERRABORT(PETSC_COMM_WORLD, error);
            t += tau;
        }
        // getArrayRead only gets local values, so make the entire std::vector
        // local
    }
    // Save the last snapshot if needed.
    if (parameters.numberOfSteps % parameters.snapshotStride == 0) {
        logger(DEBUG, "Writing snapshot at timestep %, (t=%)",
               parameters.numberOfSteps, t);
        resultVector.writeTimeIntegrationVector(snapshotIndex);
        snapshotTime[snapshotIndex] = t;
        snapshotIndex++;
    }
    // Save the number of snapshots for writing the output.
    numberOfSnapshots = snapshotIndex;
}

template <std::size_t DIM>
void DGMaxTimeIntegration<DIM>::getCoeffCO4(
    LinearAlgebra::SmallVector<6>& alpha, LinearAlgebra::SmallVector<6>& beta,
    LinearAlgebra::SmallVector<6>& alpha_sum,
    LinearAlgebra::SmallVector<6>& beta_sum,
    LinearAlgebra::SmallVector<6>& scale0,
    LinearAlgebra::SmallVector<6>& scale1) const {
    alpha[0] = 0.0;
    beta[0] = 0.0;
    scale0[0] = 0.0;
    scale1[0] = 0.0;

    alpha[1] = (146 + 5 * std::sqrt(19.0)) / 540;
    beta[5] = alpha[1];

    alpha[2] = (-2 + 10 * std::sqrt(19.0)) / 135;
    beta[4] = alpha[2];

    alpha[3] = 1.0 / 5.0;
    beta[3] = alpha[3];

    alpha[4] = (-23 - 20 * std::sqrt(19.0)) / 270;
    beta[2] = alpha[4];

    alpha[5] = (14 - std::sqrt(19.0)) / 108;
    beta[1] = alpha[5];

    for (unsigned int i = 0; i < 6; ++i) {
        alpha_sum[i] = 0.0;
        beta_sum[i] = 0.0;
        for (unsigned int j = 0; j <= i; ++j) {
            alpha_sum[i] += alpha[j];
            beta_sum[i] += beta[j];
        }
        scale0[i] = 1;
        scale1[i] = 1;
    }
    return;
}

template <std::size_t DIM>
void DGMaxTimeIntegration<DIM>::writeTimeSnapshots(std::string fileName) const {
    std::ofstream fileWriter;
    fileWriter.open(fileName);

    std::string dimensions;
    std::string variables;
    if (DIM == 2) {
        dimensions = "01";
        variables = "E0,E1,H0,H1";
    } else if (DIM == 3) {
        dimensions = "012";
        variables = "E0,E1,E2,H0,H1,H2";
    } else {
        logger.assert_debug(false, "Only supporting 2D and 3D Maxwell");
    }
    Output::TecplotDiscontinuousSolutionWriter<DIM> tecplotWriter(
        fileWriter, "The electric field", dimensions, variables);

    std::stringstream timeString;
    for (std::size_t timeLevel = 0; timeLevel < numberOfSnapshots;
         ++timeLevel) {
        writeTimeLevel(tecplotWriter, timeLevel, timeLevel == 0);
    }
    fileWriter.close();
}

template <std::size_t DIM>
void DGMaxTimeIntegration<DIM>::writeTimeLevel(
    Output::TecplotDiscontinuousSolutionWriter<DIM>& writer,
    std::size_t timeLevel, bool firstLevel) const {
    std::stringstream zoneName("t=");
    zoneName << snapshotTime[timeLevel];
    writer.write(
        &mesh_, zoneName.str(), !firstLevel,
        [&](const Base::Element* element,
            const Geometry::PointReference<DIM>& point, std::ostream& stream) {
            const LinearAlgebra::MiddleSizeVector coefficients =
                element->getTimeIntegrationVector(timeLevel);
            auto fields =
                discretization.computeFields(element, point, coefficients);
            LinearAlgebra::SmallVector<DIM> electricField =
                fields.electricField;
            LinearAlgebra::SmallVector<DIM> curlField =
                fields.electricFieldCurl;
            if (DIM == 2) {
                stream << electricField[0] << " " << electricField[1] << " "
                       << curlField[0] << " " << curlField[1] << std::endl;
            } else if (DIM == 3) {
                stream << electricField[0] << " " << electricField[1] << " "
                       << electricField[2] << " " << curlField[0] << " "
                       << curlField[1] << " " << curlField[2] << std::endl;
            }
        });
}

template <std::size_t DIM>
void DGMaxTimeIntegration<DIM>::printErrors(
    const std::vector<typename DGMaxDiscretization<DIM>::NormType>& norms,
    const typename DGMaxDiscretization<DIM>::TimeFunction& exactField,
    const typename DGMaxDiscretization<DIM>::TimeFunction& exactCurl) {
    using NormType = typename DGMaxDiscretization<DIM>::NormType;
    std::set<NormType> normSet;

    std::cout << "t";
    for (auto norm : norms) {
        std::cout << "\t" << DGMaxDiscretization<DIM>::normName(norm);
        normSet.emplace(norm);
    }
    std::cout << std::endl;

    for (std::size_t level = 0; level < numberOfSnapshots; ++level) {
        double time = snapshotTime[level];
        std::map<NormType, double> normValues = discretization.computeError(
            mesh_, level, std::bind(exactField, std::placeholders::_1, time),
            std::bind(exactCurl, std::placeholders::_1, time), normSet);
        std::cout << time;
        for (auto norm : norms) {
            std::cout << "\t" << normValues[norm];
        }
        std::cout << std::endl;
    }
}

template <std::size_t DIM>
void DGMaxTimeIntegration<DIM>::printErrors(
    const std::vector<typename DGMaxDiscretization<DIM>::NormType>& norms,
    const ExactTimeIntegrationProblem<DIM>& problem) {
    printErrors(norms,
                std::bind(&ExactTimeIntegrationProblem<DIM>::exactSolution,
                          std::ref(problem), std::placeholders::_1,
                          std::placeholders::_2),
                std::bind(&ExactTimeIntegrationProblem<DIM>::exactSolutionCurl,
                          std::ref(problem), std::placeholders::_1,
                          std::placeholders::_2));
}

template class DGMaxTimeIntegration<2>;
template class DGMaxTimeIntegration<3>;
