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

#ifndef BASISFUNCTIONS3DH1CONFORMINGPRISM_HPP_
#define BASISFUNCTIONS3DH1CONFORMINGPRISM_HPP_

#include "Base/BaseBasisFunction.h"
#include <vector>

namespace Base {
class BasisFunctionSet;
class OrientedBasisFunctionSet;
}  // namespace Base

namespace Geometry {
template <std::size_t DIM>
class PointReference;
}

namespace Utilities {

class BasisFunction3DVertexPrism : public Base::BaseBasisFunction {
   public:
    BasisFunction3DVertexPrism(std::size_t node);

    double eval(const Geometry::PointReference<3>& p) const;

    double evalDeriv0(const Geometry::PointReference<3>& p) const;

    double evalDeriv1(const Geometry::PointReference<3>& p) const;

    double evalDeriv2(const Geometry::PointReference<3>& p) const;

   private:
    int nodePosition_;
    std::size_t node_;  // node is number inside triangle
};

class BasisFunction3DEdgePrism_0 : public Base::BaseBasisFunction {
   public:
    BasisFunction3DEdgePrism_0(std::size_t node0, std::size_t node1,
                               std::size_t polynomialOrder)
        : edgePosition_((static_cast<int>(node0) / 3) * 2 - 1),
          node0_(node0 % 3),
          node1_(node1 % 3),
          polynomialOrder_(polynomialOrder) {
        logger.assert_debug(node0 < 6, "A triangular prism only has 6 nodes");
        logger.assert_debug(node1 < 6, "A triangular prism only has 6 nodes");
        logger.assert_debug(
            node0 / 3 == node1 / 3,
            "Nodes % and % do not form an edge connected to a triangular face");
    }

    double eval(const Geometry::PointReference<3>& p) const;

    double evalDeriv0(const Geometry::PointReference<3>& p) const;

    double evalDeriv1(const Geometry::PointReference<3>& p) const;

    double evalDeriv2(const Geometry::PointReference<3>& p) const;

   private:
    int edgePosition_;
    std::size_t node0_, node1_, polynomialOrder_;
};

class BasisFunction3DEdgePrism_1 : public Base::BaseBasisFunction {
   public:
    BasisFunction3DEdgePrism_1(std::size_t node0, std::size_t node1,
                               std::size_t polynomialOrder);

    double eval(const Geometry::PointReference<3>& p) const;

    double evalDeriv0(const Geometry::PointReference<3>& p) const;

    double evalDeriv1(const Geometry::PointReference<3>& p) const;

    double evalDeriv2(const Geometry::PointReference<3>& p) const;

   private:
    int mirroring_;
    std::size_t node_, polynomialOrder_;
};

class BasisFunction3DFacePrism_0 : public Base::BaseBasisFunction {
   public:
    BasisFunction3DFacePrism_0(std::size_t node0, std::size_t node1,
                               std::size_t node2, std::size_t polynomialOrder0,
                               std::size_t polynomialOrder1);

    double eval(const Geometry::PointReference<3>& p) const;

    double evalDeriv0(const Geometry::PointReference<3>& p) const;

    double evalDeriv1(const Geometry::PointReference<3>& p) const;

    double evalDeriv2(const Geometry::PointReference<3>& p) const;

   private:
    int facePosition_;
    std::size_t polynomialOrder0_, polynomialOrder1_, node0_, node1_, node2_;
};

class BasisFunction3DFacePrism_1 : public Base::BaseBasisFunction {
   public:
    BasisFunction3DFacePrism_1(std::size_t node0, std::size_t node1,
                               std::size_t node2, std::size_t polynomialOrder0,
                               std::size_t polynomialOrder1);

    double eval(const Geometry::PointReference<3>& p) const;

    double evalDeriv0(const Geometry::PointReference<3>& p) const;

    double evalDeriv1(const Geometry::PointReference<3>& p) const;

    double evalDeriv2(const Geometry::PointReference<3>& p) const;

   private:
    int mirroring_;
    std::size_t node0_, node1_, polynomialOrder0_, polynomialOrder1_;
};

class BasisFunction3DInteriorPrism : public Base::BaseBasisFunction {
   public:
    BasisFunction3DInteriorPrism(std::size_t polynomialOrder0,
                                 std::size_t polynomialOrder1,
                                 std::size_t polynomialOrder2)
        : polnomialOrder0_(polynomialOrder0),
          polynomialOrder1_(polynomialOrder1),
          polynomialOrder2_(polynomialOrder2) {}

    double eval(const Geometry::PointReference<3>& p) const;

    double evalDeriv0(const Geometry::PointReference<3>& p) const;

    double evalDeriv1(const Geometry::PointReference<3>& p) const;

    double evalDeriv2(const Geometry::PointReference<3>& p) const;

   private:
    std::size_t polnomialOrder0_, polynomialOrder1_, polynomialOrder2_;
};

Base::BasisFunctionSet* createDGBasisFunctionSet3DH1ConformingPrism(
    std::size_t order);

Base::BasisFunctionSet* createInteriorBasisFunctionSet3DH1ConformingPrism(
    std::size_t order);

std::vector<const Base::BasisFunctionSet*>
    createVertexBasisFunctionSet3DH1ConformingPrism(std::size_t order);

std::vector<const Base::OrientedBasisFunctionSet*>
    createEdgeBasisFunctionSet3DH1ConformingPrism(std::size_t order);

std::vector<const Base::OrientedBasisFunctionSet*>
    createFaceBasisFunctionSet3DH1ConformingPrism(std::size_t order);

}  // namespace Utilities

#endif /* BASISFUNCTIONS3DH1CONFORMINGPRISM_HPP_ */
