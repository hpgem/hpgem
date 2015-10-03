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

#include "BasisFunctions2DH1ConformingSquare.h"
#include "helperFunctions.h"
#include "Base/BasisFunctionSet.h"
#include "Base/OrientedBasisFunctionSet.h"
#include "Geometry/ReferenceSquare.h"
#include "Geometry/PointReference.h"
#include "Logger.h"

//only uses the constant basis functions
#include "BasisFunctionsCollection_A.h"

namespace Utilities
{
    
    double BasisFunction2DVertexSquare::eval(const Geometry::PointReference<2>& p) const
    {
        return (1 + nodePosition0_ * p[0]) * (1 + nodePosition1_ * p[1]) / 4.;
    }
    
    double BasisFunction2DVertexSquare::evalDeriv0(const Geometry::PointReference<2>& p) const
    {
        return nodePosition0_ * (1 + nodePosition1_ * p[1]) / 4.;
    }
    
    double BasisFunction2DVertexSquare::evalDeriv1(const Geometry::PointReference<2>& p) const
    {
        return nodePosition1_ * (1 + nodePosition0_ * p[0]) / 4.;
    }
    
    BasisFunction2DFaceSquare_0::BasisFunction2DFaceSquare_0(std::size_t node0, std::size_t node1, std::size_t polynomialOrder)
            : polynomialOrder_(polynomialOrder)
    {
        logger.assert(node0 < 4, "A square only has 4 nodes");
        logger.assert(node1 < 4, "A square only has 4 nodes");
        logger.assert((node0 + node1) % 2 == 1, "please use BasisFunction2DFaceSquare_1 for edges that are aligned vertically");
        mirroring_ = (node0 > node1) ? -1 : 1;
        edgePosition_ = (node0 + node1 < 3) ? -1 : 1;
    }
    
    double BasisFunction2DFaceSquare_0::eval(const Geometry::PointReference<2>& p) const
    {
        return (1 + edgePosition_ * p[1]) * (1 - p[0]) * (1 + p[0]) * LobattoPolynomial(polynomialOrder_, mirroring_ * p[0]) / 8.;
    }
    
    double BasisFunction2DFaceSquare_0::evalDeriv0(const Geometry::PointReference<2>& p) const
    {
        return (1 + edgePosition_ * p[1]) * (-p[0] * LobattoPolynomial(polynomialOrder_, mirroring_ * p[0]) + (1 - p[0]) * (1 + p[0]) * LobattoPolynomialDerivative(polynomialOrder_, mirroring_ * p[0]) * mirroring_ / 2.) / 4.;
    }
    
    double BasisFunction2DFaceSquare_0::evalDeriv1(const Geometry::PointReference<2>& p) const
    {
        return edgePosition_ * (1 - p[0]) * (1 + p[0]) * LobattoPolynomial(polynomialOrder_, mirroring_ * p[0]) / 8.;
    }
    
    BasisFunction2DFaceSquare_1::BasisFunction2DFaceSquare_1(std::size_t node0, std::size_t node1, std::size_t polynomialOrder)
            : polynomialOrder_(polynomialOrder)
    {
        logger.assert(node0 < 4, "A square only has 4 nodes");
        logger.assert(node1 < 4, "A square only has 4 nodes");
        logger.assert((node0 + node1) % 2 == 0, "please use BasisFunction2DFaceSquare_0 for edges that are aligned horizontally");
        mirroring_ = (node0 > node1) ? -1 : 1;
        edgePosition_ = (node0 + node1 < 3) ? -1 : 1;
    }
    
    double BasisFunction2DFaceSquare_1::eval(const Geometry::PointReference<2>& p) const
    {
        return (1 + edgePosition_ * p[0]) * (1 - p[1]) * (1 + p[1]) * LobattoPolynomial(polynomialOrder_, mirroring_ * p[1]) / 8.;
    }
    
    double BasisFunction2DFaceSquare_1::evalDeriv0(const Geometry::PointReference<2>& p) const
    {
        return edgePosition_ * (1 - p[1]) * (1 + p[1]) * LobattoPolynomial(polynomialOrder_, mirroring_ * p[1]) / 8.;
    }
    
    double BasisFunction2DFaceSquare_1::evalDeriv1(const Geometry::PointReference<2>& p) const
    {
        return (1 + edgePosition_ * p[0]) * (-p[1] * LobattoPolynomial(polynomialOrder_, mirroring_ * p[1]) + (1 - p[1]) * (1 + p[1]) * LobattoPolynomialDerivative(polynomialOrder_, mirroring_ * p[1]) * mirroring_ / 2.) / 4.;
    }
    
    double BasisFunction2DInteriorSquare::eval(const Geometry::PointReference<2>& p) const
    {
        return (1 - p[0]) * (1 + p[0]) * (1 - p[1]) * (1 + p[1]) * LobattoPolynomial(polynomialOrder0_, p[0]) * LobattoPolynomial(polynomialOrder1_, p[1]) / 16.;
    }
    
    double BasisFunction2DInteriorSquare::evalDeriv0(const Geometry::PointReference<2>& p) const
    {
        return LobattoPolynomial(polynomialOrder1_, p[1]) * (1 - p[1]) * (1 + p[1]) / 4. * (-p[0] * LobattoPolynomial(polynomialOrder0_, p[0]) / 2. + (1 - p[0]) * (1 + p[0]) * LobattoPolynomialDerivative(polynomialOrder0_, p[0]) / 4.);
    }
    
    double BasisFunction2DInteriorSquare::evalDeriv1(const Geometry::PointReference<2>& p) const
    {
        return LobattoPolynomial(polynomialOrder0_, p[0]) * (1 - p[0]) * (1 + p[0]) / 4. * (-p[1] * LobattoPolynomial(polynomialOrder1_, p[1]) / 2. + (1 - p[1]) * (1 + p[1]) * LobattoPolynomialDerivative(polynomialOrder1_, p[1]) / 4.);
    }
    
    Base::BasisFunctionSet* createDGBasisFunctionSet2DH1Square(std::size_t order)
    {
        Base::BasisFunctionSet* result(new Base::BasisFunctionSet(order));
        if(order > 0)
        {
            for (std::size_t i = 0; i < 4; ++i)
            {
                result->addBasisFunction(new BasisFunction2DVertexSquare(i));
            }
            for (std::size_t i = 0; i + 2 <= order; ++i)
            {
                result->addBasisFunction(new BasisFunction2DFaceSquare_0(0, 1, i));
                result->addBasisFunction(new BasisFunction2DFaceSquare_0(2, 3, i));
                result->addBasisFunction(new BasisFunction2DFaceSquare_1(0, 2, i));
                result->addBasisFunction(new BasisFunction2DFaceSquare_1(1, 3, i));
                for (std::size_t j = 0; j + 2 <= order; ++j)
                {
                    result->addBasisFunction(new BasisFunction2DInteriorSquare(i, j));
                }
            }
        }
        else
        {
            result->addBasisFunction(new Base::Basis_A0_2D);
        }
        return result;
    }
    
    Base::BasisFunctionSet* createInteriorBasisFunctionSet2DH1Square(std::size_t order)
    {
        logger.assert(order > 0, "Trying to create a conforming, constant basis function set, did you mean the constant solution?");
        Base::BasisFunctionSet* result(new Base::BasisFunctionSet(order));
        for (std::size_t i = 0; i + 2 <= order; ++i)
        {
            for (std::size_t j = 0; j + 2 <= order; ++j)
            {
                result->addBasisFunction(new BasisFunction2DInteriorSquare(i, j));
            }
        }
        return result;
    }
    
    std::vector<const Base::BasisFunctionSet*> createVertexBasisFunctionSet2DH1Square(std::size_t order)
    {
        logger.assert(order > 0, "Trying to create a conforming, constant basis function set, did you mean the constant solution?");
        std::vector<const Base::BasisFunctionSet*> result;
        Base::BasisFunctionSet* set;
        for (std::size_t i = 0; i < 4; ++i)
        {
            set = new Base::BasisFunctionSet(order);
            set->addBasisFunction(new BasisFunction2DVertexSquare(i));
            result.push_back(set);
        }
        return result;
    }
    
    std::vector<const Base::OrientedBasisFunctionSet*> createFaceBasisFunctionSet2DH1Square(std::size_t order)
    {
        logger.assert(order > 0, "Trying to create a conforming, constant basis function set, did you mean the constant solution?");
        std::vector<const Base::OrientedBasisFunctionSet*> result;
        Geometry::ReferenceSquare& square = Geometry::ReferenceSquare::Instance();
        Base::OrientedBasisFunctionSet* set;
        std::vector<std::size_t> vertexindices(2);
        for (std::size_t i = 0; i < 4; ++i)
        {
            set = new Base::OrientedBasisFunctionSet(order, 0, i);
            vertexindices = square.getCodim1EntityLocalIndices(i);
            for (std::size_t j = 0; j + 2 <= order; ++j)
            {
                if ((vertexindices[0] + vertexindices[1]) % 2 == 1)
                    set->addBasisFunction(new BasisFunction2DFaceSquare_0(vertexindices[0], vertexindices[1], j));
                else
                    set->addBasisFunction(new BasisFunction2DFaceSquare_1(vertexindices[0], vertexindices[1], j));
            }
            result.push_back(set);
            set = new Base::OrientedBasisFunctionSet(order, 1, i);
            for (std::size_t j = 0; j + 2 <= order; ++j)
            {
                if ((vertexindices[0] + vertexindices[1]) % 2 == 1)
                    set->addBasisFunction(new BasisFunction2DFaceSquare_0(vertexindices[1], vertexindices[0], j));
                else
                    set->addBasisFunction(new BasisFunction2DFaceSquare_1(vertexindices[1], vertexindices[0], j));
            }
            result.push_back(set);
        }
        return result;
    }

}
