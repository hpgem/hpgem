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

#ifndef HPGEM_KERNEL_BASISFUNCTIONS1DH1CONFORMINGLINE_H
#define HPGEM_KERNEL_BASISFUNCTIONS1DH1CONFORMINGLINE_H

#include "BaseBasisFunction.h"
#include <vector>

namespace hpgem {

namespace Geometry {
template <std::size_t DIM>
class PointReference;
}

namespace FE {

class BasisFunctionSet;
class OrientedBasisFunctionSet;

class BasisFunction1DVertexLine : public BaseBasisFunction {
   public:
    BasisFunction1DVertexLine(std::size_t node)
        : nodePosition_(2 * static_cast<int>(node) - 1) {
        logger.assert_debug(node < 2, "A line only has 2 nodes");
    }

    double eval(const Geometry::PointReference<1>& p) const override;

    double evalDeriv0(const Geometry::PointReference<1>& p) const override;

   private:
    int nodePosition_;
};

class BasisFunction1DInteriorLine : public BaseBasisFunction {
   public:
    BasisFunction1DInteriorLine(std::size_t polynomialOrder)
        : polynomialOrder_(polynomialOrder) {}

    double eval(const Geometry::PointReference<1>& p) const override;

    double evalDeriv0(const Geometry::PointReference<1>& p) const override;

   private:
    std::size_t polynomialOrder_;
};

BasisFunctionSet* createDGBasisFunctionSet1DH1Line(std::size_t polynomialOrder);

std::vector<BaseBasisFunction*> createDGBasisFunctions1DH1Line(
    std::size_t polynomialOrder);

BasisFunctionSet* createInteriorBasisFunctionSet1DH1Line(
    std::size_t polynomialOrder);

std::vector<const OrientedBasisFunctionSet*>
    createVertexBasisFunctionSet1DH1Line(std::size_t polynomialOrder);

}  // namespace FE

}  // namespace hpgem

#endif  // HPGEM_KERNEL_BASISFUNCTIONS1DH1CONFORMINGLINE_H
