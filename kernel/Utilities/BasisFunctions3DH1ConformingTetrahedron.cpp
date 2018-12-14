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

#include "BasisFunctions3DH1ConformingTetrahedron.h"
#include "helperFunctions.h"
#include "Base/BasisFunctionSet.h"
#include "Base/OrientedBasisFunctionSet.h"
#include "Geometry/ReferenceTetrahedron.h"
#include "Geometry/PointReference.h"

//only uses the constant basis functions
#include "BasisFunctionsCollection_A.h"

namespace Utilities
{
    
    double BasisFunction3DVertexTetrahedron::eval(const Geometry::PointReference<3>& p) const
    {
        return baricentric_3D(node_, p);
    }
    
    double BasisFunction3DVertexTetrahedron::evalDeriv0(const Geometry::PointReference<3>& p) const
    {
        return baricentricDeriv(node_, 0);
    }
    
    double BasisFunction3DVertexTetrahedron::evalDeriv1(const Geometry::PointReference<3>& p) const
    {
        return baricentricDeriv(node_, 1);
    }
    
    double BasisFunction3DVertexTetrahedron::evalDeriv2(const Geometry::PointReference<3>& p) const
    {
        return baricentricDeriv(node_, 2);
    }
    
    double BasisFunction3DEdgeTetrahedron::eval(const Geometry::PointReference<3>& p) const
    {
        return baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * LobattoPolynomial(polynomialOrder_, baricentric_3D(node0_, p) - baricentric_3D(node1_, p));
    }
    
    double BasisFunction3DEdgeTetrahedron::evalDeriv0(const Geometry::PointReference<3>& p) const
    {
        return baricentricDeriv(node0_, 0) * (baricentric_3D(node1_, p) * LobattoPolynomial(polynomialOrder_, baricentric_3D(node0_, p) - baricentric_3D(node1_, p)) + baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * LobattoPolynomialDerivative(polynomialOrder_, baricentric_3D(node0_, p) - baricentric_3D(node1_, p))) + baricentricDeriv(node1_, 0) * (baricentric_3D(node0_, p) * LobattoPolynomial(polynomialOrder_, baricentric_3D(node0_, p) - baricentric_3D(node1_, p)) - baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * LobattoPolynomialDerivative(polynomialOrder_, baricentric_3D(node0_, p) - baricentric_3D(node1_, p)));
    }
    
    double BasisFunction3DEdgeTetrahedron::evalDeriv1(const Geometry::PointReference<3>& p) const
    {
        return baricentricDeriv(node0_, 1) * (baricentric_3D(node1_, p) * LobattoPolynomial(polynomialOrder_, baricentric_3D(node0_, p) - baricentric_3D(node1_, p)) + baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * LobattoPolynomialDerivative(polynomialOrder_, baricentric_3D(node0_, p) - baricentric_3D(node1_, p))) + baricentricDeriv(node1_, 1) * (baricentric_3D(node0_, p) * LobattoPolynomial(polynomialOrder_, baricentric_3D(node0_, p) - baricentric_3D(node1_, p)) - baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * LobattoPolynomialDerivative(polynomialOrder_, baricentric_3D(node0_, p) - baricentric_3D(node1_, p)));
    }
    
    double BasisFunction3DEdgeTetrahedron::evalDeriv2(const Geometry::PointReference<3>& p) const
    {
        return baricentricDeriv(node0_, 2) * (baricentric_3D(node1_, p) * LobattoPolynomial(polynomialOrder_, baricentric_3D(node0_, p) - baricentric_3D(node1_, p)) + baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * LobattoPolynomialDerivative(polynomialOrder_, baricentric_3D(node0_, p) - baricentric_3D(node1_, p))) + baricentricDeriv(node1_, 2) * (baricentric_3D(node0_, p) * LobattoPolynomial(polynomialOrder_, baricentric_3D(node0_, p) - baricentric_3D(node1_, p)) - baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * LobattoPolynomialDerivative(polynomialOrder_, baricentric_3D(node0_, p) - baricentric_3D(node1_, p)));
    }
    
    double BasisFunction3DFaceTetrahedron::eval(const Geometry::PointReference<3>& p) const
    {
        double x0(baricentric_3D(node0_, p) - baricentric_3D(node1_, p)), x1(baricentric_3D(node1_, p) - baricentric_3D(node2_, p));
        return baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * baricentric_3D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1);
    }
    
    double BasisFunction3DFaceTetrahedron::evalDeriv0(const Geometry::PointReference<3>& p) const
    {
        double x0(baricentric_3D(node0_, p) - baricentric_3D(node1_, p)), x1(baricentric_3D(node1_, p) - baricentric_3D(node2_, p));
        return baricentricDeriv(node0_, 0) * (baricentric_3D(node1_, p) * baricentric_3D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) + baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * baricentric_3D(node2_, p) * LobattoPolynomialDerivative(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1)) + baricentricDeriv(node1_, 0) * (baricentric_3D(node0_, p) * baricentric_3D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) - baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * baricentric_3D(node2_, p) * LobattoPolynomialDerivative(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) + baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * baricentric_3D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomialDerivative(polynomialOrder1_, x1)) + baricentricDeriv(node2_, 0) * (baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) - baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * baricentric_3D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomialDerivative(polynomialOrder1_, x1));
    }
    
    double BasisFunction3DFaceTetrahedron::evalDeriv1(const Geometry::PointReference<3>& p) const
    {
        double x0(baricentric_3D(node0_, p) - baricentric_3D(node1_, p)), x1(baricentric_3D(node1_, p) - baricentric_3D(node2_, p));
        return baricentricDeriv(node0_, 1) * (baricentric_3D(node1_, p) * baricentric_3D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) + baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * baricentric_3D(node2_, p) * LobattoPolynomialDerivative(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1)) + baricentricDeriv(node1_, 1) * (baricentric_3D(node0_, p) * baricentric_3D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) - baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * baricentric_3D(node2_, p) * LobattoPolynomialDerivative(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) + baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * baricentric_3D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomialDerivative(polynomialOrder1_, x1)) + baricentricDeriv(node2_, 1) * (baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) - baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * baricentric_3D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomialDerivative(polynomialOrder1_, x1));
    }
    
    double BasisFunction3DFaceTetrahedron::evalDeriv2(const Geometry::PointReference<3>& p) const
    {
        double x0(baricentric_3D(node0_, p) - baricentric_3D(node1_, p)), x1(baricentric_3D(node1_, p) - baricentric_3D(node2_, p));
        return baricentricDeriv(node0_, 2) * (baricentric_3D(node1_, p) * baricentric_3D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) + baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * baricentric_3D(node2_, p) * LobattoPolynomialDerivative(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1)) + baricentricDeriv(node1_, 2) * (baricentric_3D(node0_, p) * baricentric_3D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) - baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * baricentric_3D(node2_, p) * LobattoPolynomialDerivative(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) + baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * baricentric_3D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomialDerivative(polynomialOrder1_, x1)) + baricentricDeriv(node2_, 2) * (baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) - baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * baricentric_3D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomialDerivative(polynomialOrder1_, x1));
    }
    
    double BasisFunction3DInteriorTetrahedron::eval(const Geometry::PointReference<3>& p) const
    {
        double x0(baricentric_3D(0, p) - baricentric_3D(1, p)), x1(baricentric_3D(1, p) - baricentric_3D(2, p)), x2(baricentric_3D(2, p) - baricentric_3D(3, p));
        return baricentric_3D(0, p) * baricentric_3D(1, p) * baricentric_3D(2, p) * baricentric_3D(3, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) * LobattoPolynomial(polynomialOrder2_, x2);
    }
    
    double BasisFunction3DInteriorTetrahedron::evalDeriv0(const Geometry::PointReference<3>& p) const
    {
        double x0(baricentric_3D(0, p) - baricentric_3D(1, p)), x1(baricentric_3D(1, p) - baricentric_3D(2, p)), x2(baricentric_3D(2, p) - baricentric_3D(3, p));
        return baricentricDeriv(0, 0) * (baricentric_3D(1, p) * baricentric_3D(2, p) * baricentric_3D(3, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) * LobattoPolynomial(polynomialOrder2_, x2) + baricentric_3D(0, p) * baricentric_3D(1, p) * baricentric_3D(2, p) * baricentric_3D(3, p) * LobattoPolynomialDerivative(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) * LobattoPolynomial(polynomialOrder2_, x2)) + baricentricDeriv(1, 0) * (baricentric_3D(0, p) * baricentric_3D(2, p) * baricentric_3D(3, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) * LobattoPolynomial(polynomialOrder2_, x2) - baricentric_3D(0, p) * baricentric_3D(1, p) * baricentric_3D(2, p) * baricentric_3D(3, p) * LobattoPolynomialDerivative(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) * LobattoPolynomial(polynomialOrder2_, x2) + baricentric_3D(0, p) * baricentric_3D(1, p) * baricentric_3D(2, p) * baricentric_3D(3, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomialDerivative(polynomialOrder1_, x1) * LobattoPolynomial(polynomialOrder2_, x2));
    }
    
    double BasisFunction3DInteriorTetrahedron::evalDeriv1(const Geometry::PointReference<3>& p) const
    {
        double x0(baricentric_3D(0, p) - baricentric_3D(1, p)), x1(baricentric_3D(1, p) - baricentric_3D(2, p)), x2(baricentric_3D(2, p) - baricentric_3D(3, p));
        return baricentricDeriv(0, 1) * (baricentric_3D(1, p) * baricentric_3D(2, p) * baricentric_3D(3, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) * LobattoPolynomial(polynomialOrder2_, x2) + baricentric_3D(0, p) * baricentric_3D(1, p) * baricentric_3D(2, p) * baricentric_3D(3, p) * LobattoPolynomialDerivative(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) * LobattoPolynomial(polynomialOrder2_, x2)) + baricentricDeriv(2, 1) * (baricentric_3D(0, p) * baricentric_3D(1, p) * baricentric_3D(3, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) * LobattoPolynomial(polynomialOrder2_, x2) - baricentric_3D(0, p) * baricentric_3D(1, p) * baricentric_3D(2, p) * baricentric_3D(3, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomialDerivative(polynomialOrder1_, x1) * LobattoPolynomial(polynomialOrder2_, x2) + baricentric_3D(0, p) * baricentric_3D(1, p) * baricentric_3D(2, p) * baricentric_3D(3, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) * LobattoPolynomialDerivative(polynomialOrder2_, x2));
    }
    
    double BasisFunction3DInteriorTetrahedron::evalDeriv2(const Geometry::PointReference<3>& p) const
    {
        double x0(baricentric_3D(0, p) - baricentric_3D(1, p)), x1(baricentric_3D(1, p) - baricentric_3D(2, p)), x2(baricentric_3D(2, p) - baricentric_3D(3, p));
        return baricentricDeriv(0, 2) * (baricentric_3D(1, p) * baricentric_3D(2, p) * baricentric_3D(3, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) * LobattoPolynomial(polynomialOrder2_, x2) + baricentric_3D(0, p) * baricentric_3D(1, p) * baricentric_3D(2, p) * baricentric_3D(3, p) * LobattoPolynomialDerivative(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) * LobattoPolynomial(polynomialOrder2_, x2)) + baricentricDeriv(3, 2) * (baricentric_3D(0, p) * baricentric_3D(1, p) * baricentric_3D(2, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) * LobattoPolynomial(polynomialOrder2_, x2) - baricentric_3D(0, p) * baricentric_3D(1, p) * baricentric_3D(2, p) * baricentric_3D(3, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) * LobattoPolynomialDerivative(polynomialOrder2_, x2));
    }
    
    Base::BasisFunctionSet* createDGBasisFunctionSet3DH1Tetrahedron(std::size_t order)
    {
        Base::BasisFunctionSet* result = new Base::BasisFunctionSet(order);
        if(order > 0)
        {
            Geometry::ReferenceTetrahedron& tetrahedron = Geometry::ReferenceTetrahedron::Instance();
            std::vector<std::size_t> vectorOfPointIndexes(3);
            for (std::size_t i = 0; i < 4; ++i)
            {
                result->addBasisFunction(new BasisFunction3DVertexTetrahedron(i));
            }
            for (std::size_t j = 0; j + 2 <= order; ++j)
            {
                for (std::size_t i = 0; i < 6; ++i)
                {
                    vectorOfPointIndexes = tetrahedron.getCodim2EntityLocalIndices(i);
                    result->addBasisFunction(new BasisFunction3DEdgeTetrahedron(vectorOfPointIndexes[0], vectorOfPointIndexes[1], j));
                }
                if (j > 0)
                {
                    for (std::size_t i = 0; i < 4; ++i)
                    {
                        vectorOfPointIndexes = tetrahedron.getCodim1EntityLocalIndices(i);
                        for (std::size_t k = 0; k + 1 <= j; ++k)
                        {
                            result->addBasisFunction(new BasisFunction3DFaceTetrahedron(vectorOfPointIndexes[0], vectorOfPointIndexes[1], vectorOfPointIndexes[2], j - k - 1, k));
                        }
                    }
                }
                if (j > 1)
                {
                    for (std::size_t i = 0; i + 2 <= j; ++i)
                    {
                        for (std::size_t k = 0; (i + k) + 2 <= j; ++k)
                        {
                            result->addBasisFunction(new BasisFunction3DInteriorTetrahedron(i, j - i - k - 2, k));
                        }
                    }
                }
            }
        }
        else
        {
            result->addBasisFunction(new Base::Basis_A0_3D);
        }
        return result;
    }
    
    Base::BasisFunctionSet* createInteriorBasisFunctionSet3DH1Tetrahedron(std::size_t order)
    {
        logger.assert_debug(order > 0, "Trying to create a conforming, constant basis function set, did you mean the constant solution?");
        Base::BasisFunctionSet* result = new Base::BasisFunctionSet(order);
        for (std::size_t i = 0; i + 4 <= order; ++i)
        {
            for (std::size_t j = 0; i + j + 4 <= order; ++j)
            {
                for (std::size_t k = 0; i + j + k + 4 <= order; ++k)
                {
                    result->addBasisFunction(new BasisFunction3DInteriorTetrahedron(i, j, k));
                }
            }
        }
        return result;
    }
    
    std::vector<const Base::BasisFunctionSet*> createVertexBasisFunctionSet3DH1Tetrahedron(std::size_t order)
    {
        logger.assert_debug(order > 0, "Trying to create a conforming, constant basis function set, did you mean the constant solution?");
        std::vector<const Base::BasisFunctionSet*> result;
        Base::BasisFunctionSet* set;
        for (std::size_t i = 0; i < 4; ++i)
        {
            set = new Base::BasisFunctionSet(order);
            set->addBasisFunction(new BasisFunction3DVertexTetrahedron(i));
            result.push_back(set);
        }
        return result;
    }
    
    std::vector<const Base::OrientedBasisFunctionSet*> createEdgeBasisFunctionSet3DH1Tetrahedron(std::size_t order)
    {
        logger.assert_debug(order > 0, "Trying to create a conforming, constant basis function set, did you mean the constant solution?");
        std::vector<const Base::OrientedBasisFunctionSet*> result;
        Base::OrientedBasisFunctionSet* set;
        Geometry::ReferenceTetrahedron& tetrahedron = Geometry::ReferenceTetrahedron::Instance();
        std::vector<std::size_t> vectorOfPointIndexes(2);
        for (std::size_t i = 0; i < 6; ++i)
        {
            set = new Base::OrientedBasisFunctionSet(order, 0, i);
            vectorOfPointIndexes = tetrahedron.getCodim2EntityLocalIndices(i);
            for (std::size_t j = 0; j + 2 <= order; ++j)
            {
                set->addBasisFunction(new BasisFunction3DEdgeTetrahedron(vectorOfPointIndexes[0], vectorOfPointIndexes[1], j));
            }
            result.push_back(set);
            set = new Base::OrientedBasisFunctionSet(order, 1, i);
            for (std::size_t j = 0; j + 2 <= order; ++j)
            {
                set->addBasisFunction(new BasisFunction3DEdgeTetrahedron(vectorOfPointIndexes[1], vectorOfPointIndexes[0], j));
            }
            result.push_back(set);
        }
        return result;
    }
    
    std::vector<const Base::OrientedBasisFunctionSet*> createFaceBasisFunctionSet3DH1Tetrahedron(std::size_t order)
    {
        logger.assert_debug(order > 0, "Trying to create a conforming, constant basis function set, did you mean the constant solution?");
        std::vector<const Base::OrientedBasisFunctionSet*> result;
        Base::OrientedBasisFunctionSet* set;
        Geometry::ReferenceTetrahedron& tetrahedron = Geometry::ReferenceTetrahedron::Instance();
        std::vector<std::size_t> vectorOfPointIndexes(3);
        for (std::size_t i = 0; i < 4; ++i)
        {
            vectorOfPointIndexes = tetrahedron.getCodim1EntityLocalIndices(i);
            set = new Base::OrientedBasisFunctionSet(order, 0, i);
            for (std::size_t j = 0; j + 3 <= order; ++j)
            {
                for (std::size_t k = 0; j + k + 3 <= order; ++k)
                {
                    set->addBasisFunction(new BasisFunction3DFaceTetrahedron(vectorOfPointIndexes[0], vectorOfPointIndexes[1], vectorOfPointIndexes[2], j, k));
                }
            }
            result.push_back(set);
            set = new Base::OrientedBasisFunctionSet(order, 1, i);
            for (std::size_t j = 0; j + 3 <= order; ++j)
            {
                for (std::size_t k = 0; j + k + 3 <= order; ++k)
                {
                    set->addBasisFunction(new BasisFunction3DFaceTetrahedron(vectorOfPointIndexes[0], vectorOfPointIndexes[2], vectorOfPointIndexes[1], j, k));
                }
            }
            result.push_back(set);
            set = new Base::OrientedBasisFunctionSet(order, 2, i);
            for (std::size_t j = 0; j + 3 <= order; ++j)
            {
                for (std::size_t k = 0; j + k + 3 <= order; ++k)
                {
                    set->addBasisFunction(new BasisFunction3DFaceTetrahedron(vectorOfPointIndexes[1], vectorOfPointIndexes[2], vectorOfPointIndexes[0], j, k));
                }
            }
            result.push_back(set);
            set = new Base::OrientedBasisFunctionSet(order, 3, i);
            for (std::size_t j = 0; j + 3 <= order; ++j)
            {
                for (std::size_t k = 0; j + k + 3 <= order; ++k)
                {
                    set->addBasisFunction(new BasisFunction3DFaceTetrahedron(vectorOfPointIndexes[1], vectorOfPointIndexes[0], vectorOfPointIndexes[2], j, k));
                }
            }
            result.push_back(set);
            set = new Base::OrientedBasisFunctionSet(order, 4, i);
            for (std::size_t j = 0; j + 3 <= order; ++j)
            {
                for (std::size_t k = 0; j + k + 3 <= order; ++k)
                {
                    set->addBasisFunction(new BasisFunction3DFaceTetrahedron(vectorOfPointIndexes[2], vectorOfPointIndexes[1], vectorOfPointIndexes[0], j, k));
                }
            }
            result.push_back(set);
            set = new Base::OrientedBasisFunctionSet(order, 5, i);
            for (std::size_t j = 0; j + 3 <= order; ++j)
            {
                for (std::size_t k = 0; j + k + 3 <= order; ++k)
                {
                    set->addBasisFunction(new BasisFunction3DFaceTetrahedron(vectorOfPointIndexes[2], vectorOfPointIndexes[0], vectorOfPointIndexes[1], j, k));
                }
            }
            result.push_back(set);
        }
        return result;
    }

}
