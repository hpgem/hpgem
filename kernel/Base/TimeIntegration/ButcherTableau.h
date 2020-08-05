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

#ifndef HPGEM_KERNEL_BUTCHERTABLEAU_H
#define HPGEM_KERNEL_BUTCHERTABLEAU_H

#include <vector>
#include "Logger.h"

namespace hpgem {

namespace TimeIntegration {
/** \brief Basic container for Butcher's tableaux.
 *
 *  Butcher tableaux are a tool to write down some numerical methods for solving
 *  ordinary differential equations compactly.
 *  Every tableau has a matrix a, vector b and vector c containing its
 * coefficients. Furthermore, the order and the number of stages of the
 * numerical method corresponding to the Butcher's tableau is stored.
 */
class ButcherTableau {
   public:
    virtual std::size_t getOrder() const = 0;
    virtual std::size_t getNumberOfStages() const = 0;
    std::size_t getNumStages() { return getNumberOfStages(); }
    virtual bool getTotalVariationDiminishing() const = 0;
    virtual double getA(std::size_t i, std::size_t j) const = 0;
    virtual double getB(std::size_t i) const = 0;
    virtual double getC(std::size_t i) const = 0;
    virtual bool hasErrorEstimate() const { return false; }
    /// returns an expansion coefficient ('b'), execpt this one is used to
    /// construct the error in the time derivative, rather than the time
    /// derivative
    virtual double getErrorCoefficient(std::size_t i) const {
        logger.assert_debug(!hasErrorEstimate(),
                            "This butcher tableau promises to implement "
                            "getErrorCoefficient, but fails to do so.");
        logger(ERROR,
               "This Butcher tableau does not know how to estimate errors");
        return 0.;
    }
    virtual ~ButcherTableau() = default;
};
}  // namespace TimeIntegration

namespace Base {
/// \deprecated reinsert namespace TimeIntegration into Base to conform with
/// previous location
using namespace TimeIntegration;
}  // namespace Base

}  // namespace hpgem

#endif  // HPGEM_KERNEL_BUTCHERTABLEAU_H
