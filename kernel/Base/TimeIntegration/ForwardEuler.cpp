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

#include "ForwardEuler.h"
#include "Logger.h"

namespace hpgem {

namespace TimeIntegration {
ForwardEuler::ForwardEuler() {
    order_ = 1;
    numberOfStages_ = 1;
    totalVariationDiminishing_ = true;

    // make a_
    std::vector<double> aRow;
    a_.push_back(aRow);

    // make b_ and c_
    b_ = {1.0};
    c_ = {0.0};
}

std::size_t ForwardEuler::getOrder() const { return order_; }

std::size_t ForwardEuler::getNumberOfStages() const { return numberOfStages_; }

bool ForwardEuler::getTotalVariationDiminishing() const {
    return totalVariationDiminishing_;
}

double ForwardEuler::getA(std::size_t i, std::size_t j) const {
    logger.assert_debug(i < getNumberOfStages(),
                        "Asked for stage %, but there are only % stages", i,
                        getNumberOfStages());
    logger.assert_debug(j < i,
                        "Asked for implicit coefficient %, but this is an "
                        "explicit butcher tableau",
                        j);
    return a_[i][j];
}

double ForwardEuler::getB(std::size_t i) const {
    logger.assert_debug(i < getNumberOfStages(),
                        "Asked for stage %, but there are only % stages", i,
                        getNumberOfStages());
    return b_[i];
}

double ForwardEuler::getC(std::size_t i) const {
    logger.assert_debug(i < getNumberOfStages(),
                        "Asked for stage %, but there are only % stages", i,
                        getNumberOfStages());
    return c_[i];
}
}  // namespace TimeIntegration

}  // namespace hpgem
