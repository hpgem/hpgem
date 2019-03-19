/*
This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
http://hpgem.org/about-the-code/team

This code is distributed using BSD 3-Clause License. A copy of which can found below.


Copyright (c) 2018, Univesity of Twenete
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef ALGORITHMS_DGMAXTIMEINTEGRATION_h
#define ALGORITHMS_DGMAXTIMEINTEGRATION_h

#include "Output/TecplotDiscontinuousSolutionWriter.h"

#include "../ProblemTypes/TimeIntegrationProblem.h"

#include "DGMaxDiscretization.h"
#include "../BaseExtended.h"

template<std::size_t DIM>
struct TimeIntegrationParameters;

/// \brief TimeIntegrationProblem solver using the tradition DGMax discretization
template <std::size_t DIM>
class DGMaxTimeIntegration
{

public:
    enum IntegrationMethod
    {
        CO2, CO4
    };

    DGMaxTimeIntegration(hpGemUIExtentions<DIM>& base);
    ~DGMaxTimeIntegration();
    void solve(const SeparableTimeIntegrationProblem<DIM>& input, TimeIntegrationParameters<DIM> parameters);
    void writeTimeSnapshots(std::string fileName) const;
    void printErrors(const std::vector<typename DGMaxDiscretization<DIM>::NormType>& norms,
            const typename DGMaxDiscretization<DIM>::TimeFunction& exactField,
            const typename DGMaxDiscretization<DIM>::TimeFunction& exactCurl) const;
    void printErrors(const std::vector<typename DGMaxDiscretization<DIM>::NormType>& norms,
            const ExactTimeIntegrationProblem<DIM>& problem) const;
private:

    /// \brief Coefficients for the CO4 algorithm.
    void getCoeffCO4(LinearAlgebra::SmallVector<6>& alpha, LinearAlgebra::SmallVector<6>& beta,
                     LinearAlgebra::SmallVector<6>& alpha_sum, LinearAlgebra::SmallVector<6>& beta_sum,
                     LinearAlgebra::SmallVector<6>& scale0, LinearAlgebra::SmallVector<6>& scale1) const;
    void writeTimeLevel(Output::TecplotDiscontinuousSolutionWriter<DIM>& writer,
                        std::size_t timeLevel, bool firstLevel) const;
    hpGemUIExtentions<DIM>& base_;
    DGMaxDiscretization<DIM> discretization;
    // TODO: This should be output of the solver, not a local variable.
    double * snapshotTime;
    std::size_t numberOfSnapshots;
};

template<std::size_t DIM>
struct TimeIntegrationParameters
{
    /// \brief The integration method to use.
    typename DGMaxTimeIntegration<DIM>::IntegrationMethod method;
    /// \brief The stability parameter to use
    double stab;

    /// \brief The size of each time step in the time integration.
    double timeStepSize;
    /// \brief The number of time steps for the total time integration.
    std::size_t numberOfSteps;

    /// \brief Stride of the time steps where a snapshot should be taken
    std::size_t snapshotStride;

    /// \brief Configure the time stepping parameters in the traditional style.
    ///
    /// This is intended as compatibility method with the previous style of
    /// DGMax, where the time step size was automatically computed based on
    /// several of the problem properties.
    ///
    /// Note This does not configure the stab parameter.
    ///
    /// \param method The method to use
    /// \param endTime  The end time
    /// \param polynomialOrder The order of the elements used.
    /// \param numberOfElements The number per dimension of the unit cube.
    void configureTraditional(typename DGMaxTimeIntegration<DIM>::IntegrationMethod method,
            double endTime, std::size_t polynomialOrder, std::size_t numberOfElements)
    {
        this->method = method;
        switch (method)
        {
            case DGMaxTimeIntegration<DIM>::CO2:
                timeStepSize = 0.05 / ((2 * polynomialOrder + 1) * numberOfElements);
                break;
            case DGMaxTimeIntegration<DIM>::CO4:
                timeStepSize = 0.3 / ((2 * polynomialOrder + 1) * numberOfElements);
                break;
            default:
                logger.assert_debug(false, "unknown time integration method %", method);
        }
        numberOfSteps = (std::size_t) std::ceil(endTime/timeStepSize);
    }

    std::size_t numberOfSnapshots()
    {
        // 1+ for the snapshot at t=0. Possibly 2+ if we also need an extra one
        // for the last step.
        if (numberOfSteps % snapshotStride == 0)
        {
            return 1 + numberOfSteps / snapshotStride;
        }
        else
        {
            return 2 + numberOfSteps / snapshotStride;
        }
    }
};

#endif //ALGORITHMS_DGMAXTIMEINTEGRATION_h
