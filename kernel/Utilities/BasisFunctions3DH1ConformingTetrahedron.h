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

#ifndef BASISFUNCTIONS3DH1CONFORMINGTETRAHEDRON_HPP_
#define BASISFUNCTIONS3DH1CONFORMINGTETRAHEDRON_HPP_

#include "Base/BaseBasisFunction.h"
#include <vector>

namespace Base
{
    class BasisFunctionSet;
    class OrientedBasisFunctionSet;
}

namespace Geometry
{
    template<std::size_t DIM>
    class PointReference;
}

namespace Utilities
{
    
    class BasisFunction3DVertexTetrahedron : public Base::BaseBasisFunction
    {
    public:
        BasisFunction3DVertexTetrahedron(std::size_t node)
                : node_(node)
        {
            logger.assert_debug(node < 4, "A tetrahedron only has 4 nodes");
        }
        
        double eval(const Geometry::PointReference<3>& p) const;

        double evalDeriv0(const Geometry::PointReference<3>& p) const;

        double evalDeriv1(const Geometry::PointReference<3>& p) const;

        double evalDeriv2(const Geometry::PointReference<3>& p) const;

    private:
        std::size_t node_;
    };
    
    class BasisFunction3DEdgeTetrahedron : public Base::BaseBasisFunction
    {
    public:
        BasisFunction3DEdgeTetrahedron(std::size_t node0, std::size_t node1, std::size_t polynomialOrder)
                : node0_(node0), node1_(node1), polynomialOrder_(polynomialOrder)
        {
            logger.assert_debug(node0 < 4, "A tetrahedron only has 4 nodes");
            logger.assert_debug(node1 < 4, "A tetrahedron only has 4 nodes");
        }
        
        double eval(const Geometry::PointReference<3>& p) const;

        double evalDeriv0(const Geometry::PointReference<3>& p) const;

        double evalDeriv1(const Geometry::PointReference<3>& p) const;

        double evalDeriv2(const Geometry::PointReference<3>& p) const;

    private:
        std::size_t node0_, node1_, polynomialOrder_;
    };
    
    class BasisFunction3DFaceTetrahedron : public Base::BaseBasisFunction
    {
    public:
        BasisFunction3DFaceTetrahedron(std::size_t node0, std::size_t node1, std::size_t node2, std::size_t polynomialOrder0, std::size_t polynomialOrder1)
                : node0_(node0), node1_(node1), node2_(node2), polynomialOrder0_(polynomialOrder0), polynomialOrder1_(polynomialOrder1)
        {
            logger.assert_debug(node0 < 4, "A tetrahedron only has 4 nodes");
            logger.assert_debug(node1 < 4, "A tetrahedron only has 4 nodes");
            logger.assert_debug(node2 < 4, "A tetrahedron only has 4 nodes");
        }
        
        double eval(const Geometry::PointReference<3>& p) const;

        double evalDeriv0(const Geometry::PointReference<3>& p) const;

        double evalDeriv1(const Geometry::PointReference<3>& p) const;

        double evalDeriv2(const Geometry::PointReference<3>& p) const;

    private:
        std::size_t node0_, node1_, node2_, polynomialOrder0_, polynomialOrder1_;
    };
    
    class BasisFunction3DInteriorTetrahedron : public Base::BaseBasisFunction
    {
    public:
        BasisFunction3DInteriorTetrahedron(std::size_t polynomialOrder0, std::size_t polynomialOrder1, std::size_t polynomialOrder2)
                : polynomialOrder0_(polynomialOrder0), polynomialOrder1_(polynomialOrder1), polynomialOrder2_(polynomialOrder2)
        {
        }
        
        double eval(const Geometry::PointReference<3>& p) const;

        double evalDeriv0(const Geometry::PointReference<3>& p) const;

        double evalDeriv1(const Geometry::PointReference<3>& p) const;

        double evalDeriv2(const Geometry::PointReference<3>& p) const;

    private:
        std::size_t polynomialOrder0_, polynomialOrder1_, polynomialOrder2_;
    };
    
    Base::BasisFunctionSet* createDGBasisFunctionSet3DH1Tetrahedron(std::size_t order);
    
    Base::BasisFunctionSet* createInteriorBasisFunctionSet3DH1Tetrahedron(std::size_t order);
    
    std::vector<const Base::BasisFunctionSet*> createVertexBasisFunctionSet3DH1Tetrahedron(std::size_t order);
    
    std::vector<const Base::OrientedBasisFunctionSet*> createEdgeBasisFunctionSet3DH1Tetrahedron(std::size_t order);
    
    std::vector<const Base::OrientedBasisFunctionSet*> createFaceBasisFunctionSet3DH1Tetrahedron(std::size_t order);

}

#endif
