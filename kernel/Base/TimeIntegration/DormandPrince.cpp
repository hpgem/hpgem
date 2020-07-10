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

#include "DormandPrince.h"
#include "Logger.h"

namespace hpgem {

namespace TimeIntegration {
DormandPrince::DormandPrince() {
    order_ = 5;
    numberOfStages_ = 7;
    totalVariationDiminishing_ = false;

    // make a_
    std::vector<double> aRow;
    a_.push_back(aRow);
    aRow = {1. / 5.};
    a_.push_back(aRow);
    aRow = {3. / 40., 9. / 40.};
    a_.push_back(aRow);
    aRow = {44. / 45., -56. / 15., 32. / 9.};
    a_.push_back(aRow);
    aRow = {19372. / 6561., -25360. / 2187., 64448. / 6561., -212. / 729.};
    a_.push_back(aRow);
    aRow = {9017. / 3168., -355. / 33., 46732. / 5247., 49. / 176.,
            -5103. / 18656.};
    a_.push_back(aRow);
    aRow = {35. / 384.,     0.,       500. / 1113., 125. / 192.,
            -2187. / 6784., 11. / 84.};
    a_.push_back(aRow);

    // make b_ and c_
    b_ = {35. / 384., 0., 500. / 1113., 125. / 192., -2187. / 6784.,
          11. / 84.,  0.};
    // fourth-order coefficients
    //{   5179./57600., 0., 7571./16695., 393./640., -92097./339200.,
    // 187./2100., 1./40.};
    c_ = {0., 1. / 5., 3. / 10., 4. / 5., 8. / 9., 1., 1.};
    error_ = {35. / 384. - 5179. / 57600.,
              0. - 0.,
              500. / 1113. - 7571. / 16695.,
              125. / 192. - 393. / 640.,
              -2187. / 6784. - -92097. / 339200.,
              11. / 84. - 187. / 2100.,
              0. - 1. / 40.};
}

std::size_t DormandPrince::getOrder() const { return order_; }

std::size_t DormandPrince::getNumberOfStages() const { return numberOfStages_; }

bool DormandPrince::getTotalVariationDiminishing() const {
    return totalVariationDiminishing_;
}

double DormandPrince::getA(std::size_t i, std::size_t j) const {
    logger.assert_debug(i < getNumberOfStages(),
                        "Asked for stage %, but there are only % stages", i,
                        getNumberOfStages());
    logger.assert_debug(j < i,
                        "Asked for implicit coefficient %, but this is an "
                        "explicit butcher tableau",
                        j);
    return a_[i][j];
}

double DormandPrince::getB(std::size_t i) const {
    logger.assert_debug(i < getNumberOfStages(),
                        "Asked for stage %, but there are only % stages", i,
                        getNumberOfStages());
    return b_[i];
}

double DormandPrince::getC(std::size_t i) const {
    logger.assert_debug(i < getNumberOfStages(),
                        "Asked for stage %, but there are only % stages", i,
                        getNumberOfStages());
    return c_[i];
}

double DormandPrince::getErrorCoefficient(std::size_t i) const {
    logger.assert_debug(i < getNumberOfStages(),
                        "Asked for stage %, but there are only % stages", i,
                        getNumberOfStages());
    return error_[i];
}
}  // namespace TimeIntegration

}  // namespace hpgem
