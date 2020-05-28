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

#include "AllTimeIntegrators.h"
#include "ButcherTableau.h"
#include "RK4Methods.h"
#include "ForwardEuler.h"
#include "DormandPrince.h"
#include "BogackiShampine.h"
#include "MidPoint.h"
#include "RK2TVD.h"
#include "RK3TVD.h"
#include "Logger.h"

namespace TimeIntegration {
AllTimeIntegrators::AllTimeIntegrators() {
    vecOfIntegrators_.push_back(&ForwardEuler::instance());
    vecOfIntegrators_.push_back(&MidPoint::instance());
    vecOfIntegrators_.push_back(&RK2TVD::instance());
    vecOfIntegrators_.push_back(&RK3TVD::instance());
    vecOfIntegrators_.push_back(&RK4_4::instance());
    vecOfIntegrators_.push_back(&DormandPrince::instance());
    vecOfIntegrators_.push_back(&BogackiShampine::instance());
}

AllTimeIntegrators& AllTimeIntegrators::Instance() {
    static AllTimeIntegrators theInstance;
    return theInstance;
}

ButcherTableau* AllTimeIntegrators::getRule(std::size_t order,
                                            std::size_t numberOfStages,
                                            bool totalVariationDiminishing) {

    for (ButcherTableau* rule : vecOfIntegrators_) {
        // Check if the exact rule asked exists.
        if (rule->getOrder() == order &&
            rule->getNumberOfStages() == numberOfStages &&
            rule->getTotalVariationDiminishing() == totalVariationDiminishing) {
            return rule;
        }
    }
    for (ButcherTableau* rule : vecOfIntegrators_) {

        // Relax the given TVD condition to find a rule that satisfies the order
        // and stages
        if (rule->getOrder() == order &&
            rule->getNumberOfStages() == numberOfStages &&
            rule->getTotalVariationDiminishing() != totalVariationDiminishing) {
            logger(WARN,
                   "Warning: The TVD specification for your rule does not "
                   "exist. Using an existing rule with the same order and "
                   "stage instead. ");
            return rule;
        }
    }

    logger(ERROR, "Could not find the Runge Kutta method you're looking for.");
    return nullptr;
}
}  // namespace TimeIntegration
