/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, University of Twente
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "BogackiShampine.h"

TimeIntegration::BogackiShampine::BogackiShampine()
{
    order_ = 5;
    numberOfStages_ = 7;
    totalVariationDiminishing_ = false;

    //make a_
    std::vector<double> aRow;
    a_.push_back(aRow);
    aRow =
            {   1./2.};
    a_.push_back(aRow);
    aRow =
            {   0., 3./4.};
    a_.push_back(aRow);
    aRow =
            {   2./9., 1./3., 4./9.};
    a_.push_back(aRow);

    //make b_ and c_
    b_ =
            {   2./9., 1./3., 4./9., 0.};
    //second-order coefficients
    //{   7./24., 1./4., 1./3., 1./8.};
    c_ =
            {   0., 1./5., 3./10., 4./5., 8./9., 1., 1.};
    error_ = {   2./9. - 7./24.,
                 1./3. - 1./4.,
                 4./9. - 1./3.,
                 0.    - 1./8.,};
}

std::size_t TimeIntegration::BogackiShampine::getOrder() const
{
    return 3;
}

std::size_t TimeIntegration::BogackiShampine::getNumberOfStages() const
{
    return 4;
}

bool TimeIntegration::BogackiShampine::getTotalVariationDiminishing() const
{
    return false;
}

double TimeIntegration::BogackiShampine::getA(std::size_t i, std::size_t j) const
{
    return a_[i][j];
}

double TimeIntegration::BogackiShampine::getB(std::size_t i) const
{
    return b_[i];
}

double TimeIntegration::BogackiShampine::getC(std::size_t i) const
{
    return c_[i];
}

double TimeIntegration::BogackiShampine::getErrorCoefficient(std::size_t i) const
{
    return error_[i];
}

