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

#include "BasisFunctions3DH1ConformingPrism.h"
#include "helperFunctions.h"
#include "Base/BasisFunctionSet.h"
#include "Base/OrientedBasisFunctionSet.h"
#include "Geometry/ReferenceTriangularPrism.h"
#include "Geometry/PointReference.h"

namespace Utilities
{
    
    BasisFunction3DVertexPrism::BasisFunction3DVertexPrism(std::size_t node)
    {
        logger.assert(node < 6, "A triangular prism only has 6 nodes");
        nodePosition_ = (node / 3) * 2 - 1;
        node_ = node % 3;
    }
    
    double BasisFunction3DVertexPrism::eval(const Geometry::PointReference& p) const
    {
        return baricentric_2D(node_, p) * (1 + nodePosition_ * p[2]) / 2.;
    }
    
    double BasisFunction3DVertexPrism::evalDeriv0(const Geometry::PointReference& p) const
    {
        return baricentricDeriv(node_, 0) * (1 + nodePosition_ * p[2]) / 2.;
    }
    
    double BasisFunction3DVertexPrism::evalDeriv1(const Geometry::PointReference& p) const
    {
        return baricentricDeriv(node_, 1) * (1 + nodePosition_ * p[2]) / 2.;
    }
    
    double BasisFunction3DVertexPrism::evalDeriv2(const Geometry::PointReference& p) const
    {
        return baricentric_2D(node_, p) * nodePosition_ / 2.;
    }
    
    double BasisFunction3DEdgePrism_0::eval(const Geometry::PointReference& p) const
    {
        return (1 + edgePosition_ * p[2]) * baricentric_2D(node0_, p) * baricentric_2D(node1_, p) * LobattoPolynomial(polynomialOrder_, baricentric_2D(node0_, p) - baricentric_2D(node1_, p)) / 2.;
    }
    
    double BasisFunction3DEdgePrism_0::evalDeriv0(const Geometry::PointReference& p) const
    {
        return (1 + edgePosition_ * p[2]) * baricentricDeriv(node0_, 0) * (baricentric_2D(node1_, p) * LobattoPolynomial(polynomialOrder_, baricentric_2D(node0_, p) - baricentric_2D(node1_, p)) + baricentric_2D(node0_, p) * baricentric_2D(node1_, p) * LobattoPolynomialDerivative(polynomialOrder_, baricentric_2D(node0_, p) - baricentric_2D(node1_, p))) / 2. + (1 + edgePosition_ * p[2]) * baricentricDeriv(node1_, 0) * (baricentric_2D(node0_, p) * LobattoPolynomial(polynomialOrder_, baricentric_2D(node0_, p) - baricentric_2D(node1_, p)) - baricentric_2D(node0_, p) * baricentric_2D(node1_, p) * LobattoPolynomialDerivative(polynomialOrder_, baricentric_2D(node0_, p) - baricentric_2D(node1_, p))) / 2.;
    }
    
    double BasisFunction3DEdgePrism_0::evalDeriv1(const Geometry::PointReference& p) const
    {
        return (1 + edgePosition_ * p[2]) * baricentricDeriv(node0_, 1) * (baricentric_2D(node1_, p) * LobattoPolynomial(polynomialOrder_, baricentric_2D(node0_, p) - baricentric_2D(node1_, p)) + baricentric_2D(node0_, p) * baricentric_2D(node1_, p) * LobattoPolynomialDerivative(polynomialOrder_, baricentric_2D(node0_, p) - baricentric_2D(node1_, p))) / 2. + (1 + edgePosition_ * p[2]) * baricentricDeriv(node1_, 1) * (baricentric_2D(node0_, p) * LobattoPolynomial(polynomialOrder_, baricentric_2D(node0_, p) - baricentric_2D(node1_, p)) - baricentric_2D(node0_, p) * baricentric_2D(node1_, p) * LobattoPolynomialDerivative(polynomialOrder_, baricentric_2D(node0_, p) - baricentric_2D(node1_, p))) / 2.;
    }
    
    double BasisFunction3DEdgePrism_0::evalDeriv2(const Geometry::PointReference& p) const
    {
        return edgePosition_ * baricentric_2D(node0_, p) * baricentric_2D(node1_, p) * LobattoPolynomial(polynomialOrder_, baricentric_2D(node0_, p) - baricentric_2D(node1_, p)) / 2.;
    }
    
    BasisFunction3DEdgePrism_1::BasisFunction3DEdgePrism_1(std::size_t node0, std::size_t node1, std::size_t polynomialOrder)
            : node_(node0 % 3), polynomialOrder_(polynomialOrder)
    {
        logger.assert(node0 < 6, "A triangular prism only has 6 nodes");
        logger.assert(node1 < 6, "A triangular prism only has 6 nodes");
        logger.assert(node0 % 3 == node1 % 3, "This edge is aligned next to a triangular face");
        mirroring_ = node0 < node1 ? -1 : 1;
    }
    
    double BasisFunction3DEdgePrism_1::eval(const Geometry::PointReference& p) const
    {
        return baricentric_2D(node_, p) * (1 - p[2]) * (1 + p[2]) * LobattoPolynomial(polynomialOrder_, mirroring_ * p[2]) / 4.;
    }
    
    double BasisFunction3DEdgePrism_1::evalDeriv0(const Geometry::PointReference& p) const
    {
        return baricentricDeriv(node_, 0) * (1 - p[2]) * (1 + p[2]) * LobattoPolynomial(polynomialOrder_, mirroring_ * p[2]) / 4.;
    }
    
    double BasisFunction3DEdgePrism_1::evalDeriv1(const Geometry::PointReference& p) const
    {
        return baricentricDeriv(node_, 1) * (1 - p[2]) * (1 + p[2]) * LobattoPolynomial(polynomialOrder_, mirroring_ * p[2]) / 4.;
    }
    
    double BasisFunction3DEdgePrism_1::evalDeriv2(const Geometry::PointReference& p) const
    {
        return baricentric_2D(node_, p) * (-p[2] * LobattoPolynomial(polynomialOrder_, mirroring_ * p[2]) / 2. + (1 - p[2]) * (1 + p[2]) * LobattoPolynomialDerivative(polynomialOrder_, mirroring_ * p[2]) * mirroring_ / 4.);
    }
    
    BasisFunction3DFacePrism_0::BasisFunction3DFacePrism_0(std::size_t node0, std::size_t node1, std::size_t node2, std::size_t polynomialOrder0, std::size_t polynomialOrder1)
            : polynomialOrder0_(polynomialOrder0), polynomialOrder1_(polynomialOrder1), node0_(node0 % 3), node1_(node1 % 3), node2_(node2 % 3)
    {
        logger.assert(node0 < 6, "A triangular prism only has 6 nodes");
        logger.assert(node1 < 6, "A triangular prism only has 6 nodes");
        logger.assert(node2 < 6, "A triangular prism only has 6 nodes");
        logger.assert(node0 / 3 == node1 / 3, "This is not a triangular face");
        logger.assert(node0 / 3 == node2 / 3, "This is not a triangular face");
        facePosition_ = (node0 / 3) * 2 - 1;
    }
    
    double BasisFunction3DFacePrism_0::eval(const Geometry::PointReference& p) const
    {
        double x0(baricentric_2D(node0_, p) - baricentric_2D(node1_, p)), x1(baricentric_2D(node1_, p) - baricentric_2D(node2_, p));
        return (1 + facePosition_ * p[2]) * baricentric_2D(node0_, p) * baricentric_2D(node1_, p) * baricentric_2D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) / 2.;
    }
    
    double BasisFunction3DFacePrism_0::evalDeriv0(const Geometry::PointReference& p) const
    {
        double x0(baricentric_2D(node0_, p) - baricentric_2D(node1_, p)), x1(baricentric_2D(node1_, p) - baricentric_2D(node2_, p));
        return (1 + facePosition_ * p[2]) * baricentricDeriv(node0_, 0) * (baricentric_2D(node1_, p) * baricentric_2D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) + baricentric_2D(node0_, p) * baricentric_2D(node1_, p) * baricentric_2D(node2_, p) * LobattoPolynomialDerivative(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1)) / 2. + (1 + facePosition_ * p[2]) * baricentricDeriv(node1_, 0) * (baricentric_2D(node0_, p) * baricentric_2D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) - baricentric_2D(node0_, p) * baricentric_2D(node1_, p) * baricentric_2D(node2_, p) * LobattoPolynomialDerivative(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) + baricentric_2D(node0_, p) * baricentric_2D(node1_, p) * baricentric_2D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomialDerivative(polynomialOrder1_, x1)) / 2. + (1 + facePosition_ * p[2]) * baricentricDeriv(node2_, 0) * (baricentric_2D(node0_, p) * baricentric_2D(node1_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) - baricentric_2D(node0_, p) * baricentric_2D(node1_, p) * baricentric_2D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomialDerivative(polynomialOrder1_, x1)) / 2.;
    }
    
    double BasisFunction3DFacePrism_0::evalDeriv1(const Geometry::PointReference& p) const
    {
        double x0(baricentric_2D(node0_, p) - baricentric_2D(node1_, p)), x1(baricentric_2D(node1_, p) - baricentric_2D(node2_, p));
        return (1 + facePosition_ * p[2]) * baricentricDeriv(node0_, 1) * (baricentric_2D(node1_, p) * baricentric_2D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) + baricentric_2D(node0_, p) * baricentric_2D(node1_, p) * baricentric_2D(node2_, p) * LobattoPolynomialDerivative(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1)) / 2. + (1 + facePosition_ * p[2]) * baricentricDeriv(node1_, 1) * (baricentric_2D(node0_, p) * baricentric_2D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) - baricentric_2D(node0_, p) * baricentric_2D(node1_, p) * baricentric_2D(node2_, p) * LobattoPolynomialDerivative(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) + baricentric_2D(node0_, p) * baricentric_2D(node1_, p) * baricentric_2D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomialDerivative(polynomialOrder1_, x1)) / 2. + (1 + facePosition_ * p[2]) * baricentricDeriv(node2_, 1) * (baricentric_2D(node0_, p) * baricentric_2D(node1_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) - baricentric_2D(node0_, p) * baricentric_2D(node1_, p) * baricentric_2D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomialDerivative(polynomialOrder1_, x1)) / 2.;
    }
    
    double BasisFunction3DFacePrism_0::evalDeriv2(const Geometry::PointReference& p) const
    {
        double x0(baricentric_2D(node0_, p) - baricentric_2D(node1_, p)), x1(baricentric_2D(node1_, p) - baricentric_2D(node2_, p));
        return facePosition_ * baricentric_2D(node0_, p) * baricentric_2D(node1_, p) * baricentric_2D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) / 2.;
    }
    
    BasisFunction3DFacePrism_1::BasisFunction3DFacePrism_1(std::size_t node0, std::size_t node1, std::size_t node2, std::size_t polynomialOrder0, std::size_t polynomialOrder1)
            : node0_(node0 % 3)
    {
        logger.assert(node0 < 6, "A triangular prism only has 6 nodes");
        logger.assert(node1 < 6, "A triangular prism only has 6 nodes");
        logger.assert(node2 < 6, "A triangular prism only has 6 nodes");
        if (node1 % 3 == node0_)
        {
            mirroring_ = node0 < node1 ? -1 : 1;
            polynomialOrder0_ = polynomialOrder0;
            polynomialOrder1_ = polynomialOrder1;
            node1_ = node2 % 3;
        }
        else
        {
            mirroring_ = node0 < node2 ? -1 : 1;
            polynomialOrder0_ = polynomialOrder1;
            polynomialOrder1_ = polynomialOrder0;
            node1_ = node1 % 3;
        }
    }
    
    double BasisFunction3DFacePrism_1::eval(const Geometry::PointReference& p) const
    {
        return (1 + p[2]) * (1 - p[2]) * baricentric_2D(node0_, p) * baricentric_2D(node1_, p) * LobattoPolynomial(polynomialOrder0_, baricentric_2D(node0_, p) - baricentric_2D(node1_, p)) * LobattoPolynomial(polynomialOrder1_, mirroring_ * p[2]) / 4.;
    }
    
    double BasisFunction3DFacePrism_1::evalDeriv0(const Geometry::PointReference& p) const
    {
        return (1 + p[2]) * (1 - p[2]) * LobattoPolynomial(polynomialOrder1_, mirroring_ * p[2]) * baricentricDeriv(node0_, 0) * (baricentric_2D(node1_, p) * LobattoPolynomial(polynomialOrder0_, baricentric_2D(node0_, p) - baricentric_2D(node1_, p)) + baricentric_2D(node0_, p) * baricentric_2D(node1_, p) * LobattoPolynomialDerivative(polynomialOrder0_, baricentric_2D(node0_, p) - baricentric_2D(node1_, p))) / 4. + (1 + p[2]) * (1 - p[2]) * LobattoPolynomial(polynomialOrder1_, mirroring_ * p[2]) * baricentricDeriv(node1_, 0) * (baricentric_2D(node0_, p) * LobattoPolynomial(polynomialOrder0_, baricentric_2D(node0_, p) - baricentric_2D(node1_, p)) - baricentric_2D(node0_, p) * baricentric_2D(node1_, p) * LobattoPolynomialDerivative(polynomialOrder0_, baricentric_2D(node0_, p) - baricentric_2D(node1_, p))) / 4.;
    }
    
    double BasisFunction3DFacePrism_1::evalDeriv1(const Geometry::PointReference& p) const
    {
        return (1 + p[2]) * (1 - p[2]) * LobattoPolynomial(polynomialOrder1_, mirroring_ * p[2]) * baricentricDeriv(node0_, 1) * (baricentric_2D(node1_, p) * LobattoPolynomial(polynomialOrder0_, baricentric_2D(node0_, p) - baricentric_2D(node1_, p)) + baricentric_2D(node0_, p) * baricentric_2D(node1_, p) * LobattoPolynomialDerivative(polynomialOrder0_, baricentric_2D(node0_, p) - baricentric_2D(node1_, p))) / 4. + (1 + p[2]) * (1 - p[2]) * LobattoPolynomial(polynomialOrder1_, mirroring_ * p[2]) * baricentricDeriv(node1_, 1) * (baricentric_2D(node0_, p) * LobattoPolynomial(polynomialOrder0_, baricentric_2D(node0_, p) - baricentric_2D(node1_, p)) - baricentric_2D(node0_, p) * baricentric_2D(node1_, p) * LobattoPolynomialDerivative(polynomialOrder0_, baricentric_2D(node0_, p) - baricentric_2D(node1_, p))) / 4.;
    }
    
    double BasisFunction3DFacePrism_1::evalDeriv2(const Geometry::PointReference& p) const
    {
        return baricentric_2D(node0_, p) * baricentric_2D(node1_, p) * LobattoPolynomial(polynomialOrder0_, baricentric_2D(node0_, p) - baricentric_2D(node1_, p)) * (-p[2] * LobattoPolynomial(polynomialOrder1_, mirroring_ * p[2]) / 2. + (1 + p[2]) * (1 - p[2]) * LobattoPolynomialDerivative(polynomialOrder1_, mirroring_ * p[2]) * mirroring_ / 4.);
    }
    
    double BasisFunction3DInteriorPrism::eval(const Geometry::PointReference& p) const
    {
        double x0(baricentric_2D(0, p) - baricentric_2D(1, p)), x1(baricentric_2D(1, p) - baricentric_2D(2, p));
        return baricentric_2D(0, p) * baricentric_2D(1, p) * baricentric_2D(2, p) * (1 - p[2]) * (1 + p[2]) * LobattoPolynomial(polnomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) * LobattoPolynomial(polynomialOrder2_, p[2]) / 4.;
    }
    
    double BasisFunction3DInteriorPrism::evalDeriv0(const Geometry::PointReference& p) const
    {
        double x0(baricentric_2D(0, p) - baricentric_2D(1, p)), x1(baricentric_2D(1, p) - baricentric_2D(2, p));
        return (1 - p[2]) * (1 + p[2]) * LobattoPolynomial(polynomialOrder2_, p[2]) / 4. * (baricentricDeriv(0, 0) * (baricentric_2D(1, p) * baricentric_2D(2, p) * LobattoPolynomial(polnomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) + baricentric_2D(0, p) * baricentric_2D(1, p) * baricentric_2D(2, p) * LobattoPolynomialDerivative(polnomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1)) + baricentricDeriv(1, 0) * (baricentric_2D(0, p) * baricentric_2D(2, p) * LobattoPolynomial(polnomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) - baricentric_2D(0, p) * baricentric_2D(1, p) * baricentric_2D(2, p) * LobattoPolynomialDerivative(polnomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) + baricentric_2D(0, p) * baricentric_2D(1, p) * baricentric_2D(2, p) * LobattoPolynomial(polnomialOrder0_, x0) * LobattoPolynomialDerivative(polynomialOrder1_, x1)));
    }
    
    double BasisFunction3DInteriorPrism::evalDeriv1(const Geometry::PointReference& p) const
    {
        double x0(baricentric_2D(0, p) - baricentric_2D(1, p)), x1(baricentric_2D(1, p) - baricentric_2D(2, p));
        return (1 - p[2]) * (1 + p[2]) * LobattoPolynomial(polynomialOrder2_, p[2]) / 4. * (baricentricDeriv(0, 1) * (baricentric_2D(1, p) * baricentric_2D(2, p) * LobattoPolynomial(polnomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) + baricentric_2D(0, p) * baricentric_2D(1, p) * baricentric_2D(2, p) * LobattoPolynomialDerivative(polnomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1)) + baricentricDeriv(2, 1) * (baricentric_2D(0, p) * baricentric_2D(1, p) * LobattoPolynomial(polnomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) - baricentric_2D(0, p) * baricentric_2D(1, p) * baricentric_2D(2, p) * LobattoPolynomial(polnomialOrder0_, x0) * LobattoPolynomialDerivative(polynomialOrder1_, x1)));
    }
    
    double BasisFunction3DInteriorPrism::evalDeriv2(const Geometry::PointReference& p) const
    {
        double x0(baricentric_2D(0, p) - baricentric_2D(1, p)), x1(baricentric_2D(1, p) - baricentric_2D(2, p));
        return baricentric_2D(0, p) * baricentric_2D(1, p) * baricentric_2D(2, p) * LobattoPolynomial(polnomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) * (-p[2] * LobattoPolynomial(polynomialOrder2_, p[2]) / 2. + (1 - p[2]) * (1 + p[2]) * LobattoPolynomialDerivative(polynomialOrder2_, p[2]) / 4.);
    }
    
    Base::BasisFunctionSet* createDGBasisFunctionSet3DH1ConformingPrism(std::size_t order)
    {
        Base::BasisFunctionSet* result = new Base::BasisFunctionSet(order);
        Geometry::ReferenceTriangularPrism& prism = Geometry::ReferenceTriangularPrism::Instance();
        std::vector<std::size_t> vectorOfPointIndexes(4);
        for (std::size_t i = 0; i < 6; ++i)
        {
            result->addBasisFunction(new BasisFunction3DVertexPrism(i));
        }
        for (std::size_t j = 0; j + 2 <= order; ++j)
        {
            for (std::size_t i = 0; i < 6; ++i)
            {
                vectorOfPointIndexes = prism.getCodim2EntityLocalIndices(i);
                result->addBasisFunction(new BasisFunction3DEdgePrism_0(vectorOfPointIndexes[0], vectorOfPointIndexes[1], j));
            }
            for (std::size_t i = 6; i < 9; ++i)
            {
                vectorOfPointIndexes = prism.getCodim2EntityLocalIndices(i);
                result->addBasisFunction(new BasisFunction3DEdgePrism_1(vectorOfPointIndexes[0], vectorOfPointIndexes[1], j));
            }
            for (std::size_t i = 0; i < 2; ++i)
            {
                vectorOfPointIndexes = prism.getCodim1EntityLocalIndices(i);
                if (j > 0)
                {
                    for (std::size_t k = 0; k < j; ++k)
                    {
                        result->addBasisFunction(new BasisFunction3DFacePrism_0(vectorOfPointIndexes[0], vectorOfPointIndexes[1], vectorOfPointIndexes[2], j - k - 1, k));
                    }
                }
            }
            for (std::size_t i = 2; i < 5; ++i)
            {
                vectorOfPointIndexes = prism.getCodim1EntityLocalIndices(i);
                result->addBasisFunction(new BasisFunction3DFacePrism_1(vectorOfPointIndexes[0], vectorOfPointIndexes[1], vectorOfPointIndexes[2], j, j));
                for (std::size_t k = 0; k < j; ++k)
                {
                    result->addBasisFunction(new BasisFunction3DFacePrism_1(vectorOfPointIndexes[0], vectorOfPointIndexes[1], vectorOfPointIndexes[2], j, k));
                    result->addBasisFunction(new BasisFunction3DFacePrism_1(vectorOfPointIndexes[0], vectorOfPointIndexes[1], vectorOfPointIndexes[2], k, j));
                }
            }
            for (std::size_t i = 0; i < j; ++i)
            {
                for (std::size_t k = 0; i + k < j; ++k)
                {
                    result->addBasisFunction(new BasisFunction3DInteriorPrism(i, j, k));
                    result->addBasisFunction(new BasisFunction3DInteriorPrism(j, i, k));
                    result->addBasisFunction(new BasisFunction3DInteriorPrism(i, k, j));
                }
            }
        }
        return result;
    }
    
    Base::BasisFunctionSet* createInteriorBasisFunctionSet3DH1ConformingPrism(std::size_t order)
    {
        Base::BasisFunctionSet* result = new Base::BasisFunctionSet(order);
        for (std::size_t i = 0; i + 3 <= order; ++i)
        {
            for (int j = 0; i + j + 3 <= order; ++j)
            {
                for (int k = 0; k + 2 <= order; ++k)
                {
                    result->addBasisFunction(new BasisFunction3DInteriorPrism(i, j, k));
                }
            }
        }
        return result;
    }
    
    std::vector<const Base::BasisFunctionSet*> createVertexBasisFunctionSet3DH1ConformingPrism(std::size_t order)
    {
        std::vector<const Base::BasisFunctionSet*> result;
        Base::BasisFunctionSet* set;
        for (std::size_t i = 0; i < 6; ++i)
        {
            set = new Base::BasisFunctionSet(order);
            set->addBasisFunction(new BasisFunction3DVertexPrism(i));
            result.push_back(set);
        }
        return result;
    }
    
    std::vector<const Base::OrientedBasisFunctionSet*> createEdgeBasisFunctionSet3DH1ConformingPrism(std::size_t order)
    {
        std::vector<const Base::OrientedBasisFunctionSet*> result;
        Base::OrientedBasisFunctionSet* set;
        Geometry::ReferenceTriangularPrism& prism = Geometry::ReferenceTriangularPrism::Instance();
        std::vector<std::size_t> vectorOfPointIndexes(2);
        for (std::size_t i = 0; i < 6; ++i)
        {
            vectorOfPointIndexes = prism.getCodim2EntityLocalIndices(i);
            set = new Base::OrientedBasisFunctionSet(order, 0, i);
            for (std::size_t j = 0; j + 2 <= order; ++j)
            {
                set->addBasisFunction(new BasisFunction3DEdgePrism_0(vectorOfPointIndexes[0], vectorOfPointIndexes[1], j));
            }
            result.push_back(set);
            set = new Base::OrientedBasisFunctionSet(order, 1, i);
            for (std::size_t j = 0; j + 2 <= order; ++j)
            {
                set->addBasisFunction(new BasisFunction3DEdgePrism_0(vectorOfPointIndexes[1], vectorOfPointIndexes[0], j));
            }
            result.push_back(set);
        }
        for (int i = 6; i < 9; ++i)
        {
            vectorOfPointIndexes = prism.getCodim2EntityLocalIndices(i);
            set = new Base::OrientedBasisFunctionSet(order, 0, i);
            for (std::size_t j = 0; j + 2 <= order; ++j)
            {
                set->addBasisFunction(new BasisFunction3DEdgePrism_1(vectorOfPointIndexes[0], vectorOfPointIndexes[1], j));
            }
            result.push_back(set);
            set = new Base::OrientedBasisFunctionSet(order, 1, i);
            for (std::size_t j = 0; j + 2 <= order; ++j)
            {
                set->addBasisFunction(new BasisFunction3DEdgePrism_1(vectorOfPointIndexes[1], vectorOfPointIndexes[0], j));
            }
            result.push_back(set);
        }
        return result;
    }
    
    std::vector<const Base::OrientedBasisFunctionSet*> CreateFaceBasisFunctionSet3DH1ConformingPrism(std::size_t order)
    {
        std::vector<const Base::OrientedBasisFunctionSet*> result;
        Base::OrientedBasisFunctionSet* set;
        Geometry::ReferenceTriangularPrism& prism = Geometry::ReferenceTriangularPrism::Instance();
        std::vector<std::size_t> vectorOfPointIndexes(4);
        for (std::size_t i = 0; i < 2; ++i)
        {
            vectorOfPointIndexes = prism.getCodim1EntityLocalIndices(i);
            set = new Base::OrientedBasisFunctionSet(order, 0, i);
            for (std::size_t j = 0; j + 3 <= order; ++j)
            {
                for (int k = 0; j + k + 3 <= order; ++k)
                {
                    set->addBasisFunction(new BasisFunction3DFacePrism_0(vectorOfPointIndexes[0], vectorOfPointIndexes[1], vectorOfPointIndexes[2], j, k));
                }
            }
            result.push_back(set);
            set = new Base::OrientedBasisFunctionSet(order, 1, i);
            for (std::size_t j = 0; j + 3 <= order; ++j)
            {
                for (int k = 0; j + k + 3 <= order; ++k)
                {
                    set->addBasisFunction(new BasisFunction3DFacePrism_0(vectorOfPointIndexes[0], vectorOfPointIndexes[2], vectorOfPointIndexes[1], j, k));
                }
            }
            result.push_back(set);
            set = new Base::OrientedBasisFunctionSet(order, 2, i);
            for (std::size_t j = 0; j + 3 <= order; ++j)
            {
                for (int k = 0; j + k + 3 <= order; ++k)
                {
                    set->addBasisFunction(new BasisFunction3DFacePrism_0(vectorOfPointIndexes[1], vectorOfPointIndexes[2], vectorOfPointIndexes[0], j, k));
                }
            }
            result.push_back(set);
            set = new Base::OrientedBasisFunctionSet(order, 3, i);
            for (std::size_t j = 0; j + 3 <= order; ++j)
            {
                for (int k = 0; j + k + 3 <= order; ++k)
                {
                    set->addBasisFunction(new BasisFunction3DFacePrism_0(vectorOfPointIndexes[1], vectorOfPointIndexes[0], vectorOfPointIndexes[2], j, k));
                }
            }
            result.push_back(set);
            set = new Base::OrientedBasisFunctionSet(order, 4, i);
            for (std::size_t j = 0; j + 3 <= order; ++j)
            {
                for (int k = 0; j + k + 3 <= order; ++k)
                {
                    set->addBasisFunction(new BasisFunction3DFacePrism_0(vectorOfPointIndexes[2], vectorOfPointIndexes[1], vectorOfPointIndexes[0], j, k));
                }
            }
            result.push_back(set);
            set = new Base::OrientedBasisFunctionSet(order, 5, i);
            for (std::size_t j = 0; j + 3 <= order; ++j)
            {
                for (int k = 0; j + k + 3 <= order; ++k)
                {
                    set->addBasisFunction(new BasisFunction3DFacePrism_0(vectorOfPointIndexes[2], vectorOfPointIndexes[0], vectorOfPointIndexes[1], j, k));
                }
            }
            result.push_back(set);
        }
        for (std::size_t i = 2; i < 5; ++i)
        {
            vectorOfPointIndexes = prism.getCodim1EntityLocalIndices(i);
            set = new Base::OrientedBasisFunctionSet(order, 0, i);
            for (std::size_t j = 0; j + 2 <= order; ++j)
            {
                for (std::size_t k = 0; k + 2 <= order; ++k)
                {
                    set->addBasisFunction(new BasisFunction3DFacePrism_0(vectorOfPointIndexes[0], vectorOfPointIndexes[1], vectorOfPointIndexes[2], j, k));
                }
            }
            result.push_back(set);
            set = new Base::OrientedBasisFunctionSet(order, 1, i);
            for (std::size_t j = 0; j + 2 <= order; ++j)
            {
                for (std::size_t k = 0; k + 2 <= order; ++k)
                {
                    set->addBasisFunction(new BasisFunction3DFacePrism_0(vectorOfPointIndexes[1], vectorOfPointIndexes[2], vectorOfPointIndexes[3], j, k));
                }
            }
            result.push_back(set);
            set = new Base::OrientedBasisFunctionSet(order, 2, i);
            for (std::size_t j = 0; j + 2 <= order; ++j)
            {
                for (std::size_t k = 0; k + 2 <= order; ++k)
                {
                    set->addBasisFunction(new BasisFunction3DFacePrism_0(vectorOfPointIndexes[2], vectorOfPointIndexes[3], vectorOfPointIndexes[0], j, k));
                }
            }
            result.push_back(set);
            set = new Base::OrientedBasisFunctionSet(order, 3, i);
            for (std::size_t j = 0; j + 2 <= order; ++j)
            {
                for (std::size_t k = 0; k + 2 <= order; ++k)
                {
                    set->addBasisFunction(new BasisFunction3DFacePrism_0(vectorOfPointIndexes[3], vectorOfPointIndexes[0], vectorOfPointIndexes[1], j, k));
                }
            }
            result.push_back(set);
            set = new Base::OrientedBasisFunctionSet(order, 4, i);
            for (std::size_t j = 0; j + 2 <= order; ++j)
            {
                for (std::size_t k = 0; k + 2 <= order; ++k)
                {
                    set->addBasisFunction(new BasisFunction3DFacePrism_0(vectorOfPointIndexes[3], vectorOfPointIndexes[2], vectorOfPointIndexes[1], j, k));
                }
            }
            result.push_back(set);
            set = new Base::OrientedBasisFunctionSet(order, 5, i);
            for (std::size_t j = 0; j + 2 <= order; ++j)
            {
                for (std::size_t k = 0; k + 2 <= order; ++k)
                {
                    set->addBasisFunction(new BasisFunction3DFacePrism_0(vectorOfPointIndexes[1], vectorOfPointIndexes[0], vectorOfPointIndexes[3], j, k));
                }
            }
            result.push_back(set);
            set = new Base::OrientedBasisFunctionSet(order, 6, i);
            for (std::size_t j = 0; j + 2 <= order; ++j)
            {
                for (std::size_t k = 0; k + 2 <= order; ++k)
                {
                    set->addBasisFunction(new BasisFunction3DFacePrism_0(vectorOfPointIndexes[2], vectorOfPointIndexes[1], vectorOfPointIndexes[0], j, k));
                }
            }
            result.push_back(set);
            set = new Base::OrientedBasisFunctionSet(order, 7, i);
            for (std::size_t j = 0; j + 2 <= order; ++j)
            {
                for (std::size_t k = 0; k + 2 <= order; ++k)
                {
                    set->addBasisFunction(new BasisFunction3DFacePrism_0(vectorOfPointIndexes[0], vectorOfPointIndexes[3], vectorOfPointIndexes[2], j, k));
                }
            }
            result.push_back(set);
        }
        return result;
    }

}
