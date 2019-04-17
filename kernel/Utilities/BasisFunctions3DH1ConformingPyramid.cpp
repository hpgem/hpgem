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

#include "BasisFunctions3DH1ConformingPyramid.h"
#include "helperFunctions.h"
#include "Base/BasisFunctionSet.h"
#include "Base/OrientedBasisFunctionSet.h"
#include "Geometry/ReferenceTriangularPrism.h"
#include "Geometry/PointReference.h"

//only uses the constant basis functions
#include "BasisFunctionsPiecewiseConstant.h"

namespace Utilities {

    BasisFunction3DVertexPyramid::BasisFunction3DVertexPyramid(std::size_t node)
        : node_(node)
    {
        logger.assert_debug(node < 5, "A pyramid only has 5 nodes");
    }

    double BasisFunction3DVertexPyramid::eval(const Geometry::PointReference<3> &p) const
    {
        switch(node_)
        {
        case 0: //(tip)
            return p[2];
            //the basisfunctions should be bilinear on the bottom (z=0)
            //and linear on the other faces dus to conformance requirements
            //so on the triangular faces one of the components of the tensor product should either be 0 or be divided away
            //the triangular faces are x=1-z; -x=1-z; y=1-z; -y=1-z
        case 1:
            return (1-p[2]-p[0])*(1-p[2]-p[1])/4/(1-p[2] + 1e-50);
        case 2:
            return (1-p[2]+p[0])*(1-p[2]-p[1])/4/(1-p[2] + 1e-50);
        case 3:
            return (1-p[2]-p[0])*(1-p[2]+p[1])/4/(1-p[2] + 1e-50);
        case 4:
            return (1-p[2]+p[0])*(1-p[2]+p[1])/4/(1-p[2] + 1e-50);
        }
        logger(ERROR, "Something magic happened, please run again with memory checking tools enabled");
        return 0;
    }

    double BasisFunction3DVertexPyramid::evalDeriv0(const Geometry::PointReference<3> &pR) const
    {
        auto p = pR;
        if(std::abs(p[0]) + std::abs(p[1]) + std::abs(p[2] - 1.) < 1e-14)
        {
            p[2] -= 1e-10;
        }
        switch(node_)
        {
        case 0:
            return 0;
        case 1:
            return -(1-p[2]-p[1])/4/(1-p[2]);
        case 2:
            return (1-p[2]-p[1])/4/(1-p[2]);
        case 3:
            return -(1-p[2]+p[1])/4/(1-p[2]);
        case 4:
            return (1-p[2]+p[1])/4/(1-p[2]);
        }
        logger(ERROR, "Something magic happened, please run again with memory checking tools enabled");
        return 0;
    }

    double BasisFunction3DVertexPyramid::evalDeriv1(const Geometry::PointReference<3> &pR) const
    {
        auto p = pR;
        if(std::abs(p[0]) + std::abs(p[1]) + std::abs(p[2] - 1.) < 1e-14)
        {
            p[2] -= 1e-10;
        }
        switch(node_)
        {
        case 0:
            return 0;
        case 1:
            return -(1-p[2]-p[0])/4/(1-p[2]);
        case 2:
            return -(1-p[2]+p[0])/4/(1-p[2]);
        case 3:
            return (1-p[2]-p[0])/4/(1-p[2]);
        case 4:
            return (1-p[2]+p[0])/4/(1-p[2]);
        }
        logger(ERROR, "Something magic happened, please run again with memory checking tools enabled");
        return 0;
    }

    double BasisFunction3DVertexPyramid::evalDeriv2(const Geometry::PointReference<3> &pR) const
    {
        auto p = pR;
        if(std::abs(p[0]) + std::abs(p[1]) + std::abs(p[2] - 1.) < 1e-14)
        {
            p[2] -= 1e-10;
        }
        switch(node_)
        {
        case 0:
            return 1;
        case 1:
            return (1-p[2]-p[0])*(1-p[2]-p[1])/4/(1-p[2])/(1-p[2])-(1-p[2]-p[1])/4/(1-p[2])-(1-p[2]-p[0])/4/(1-p[2]);
        case 2:
            return (1-p[2]+p[0])*(1-p[2]-p[1])/4/(1-p[2])/(1-p[2])-(1-p[2]-p[1])/4/(1-p[2])-(1-p[2]+p[0])/4/(1-p[2]);
        case 3:
            return (1-p[2]-p[0])*(1-p[2]+p[1])/4/(1-p[2])/(1-p[2])-(1-p[2]+p[1])/4/(1-p[2])-(1-p[2]-p[0])/4/(1-p[2]);
        case 4:
            return (1-p[2]+p[0])*(1-p[2]+p[1])/4/(1-p[2])/(1-p[2])-(1-p[2]+p[1])/4/(1-p[2])-(1-p[2]+p[0])/4/(1-p[2]);
        }
        logger(ERROR, "Something magic happened, please run again with memory checking tools enabled");
        return 0;
    }

    Base::BasisFunctionSet* createDGBasisFunctionSet3DH1ConformingPyramid(std::size_t order)
    {
        logger.assert_debug(order < 2, "Only linear basis functions have been implemented so far");
        auto result = new Base::BasisFunctionSet(order);
        if(order == 1)
        {
            result->addBasisFunction(new BasisFunction3DVertexPyramid(0));
            result->addBasisFunction(new BasisFunction3DVertexPyramid(1));
            result->addBasisFunction(new BasisFunction3DVertexPyramid(2));
            result->addBasisFunction(new BasisFunction3DVertexPyramid(3));
            result->addBasisFunction(new BasisFunction3DVertexPyramid(4));
        }
        else
        {
            addPiecewiseConstantBasisFunction3D(*result);
        }
        return result;
    }

    Base::BasisFunctionSet* createInteriorBasisFunctionSet3DH1ConformingPyramid(std::size_t order)
    {
        logger.assert_debug(order == 1, "Only linear basis functions have been implemented so far");
        return new Base::BasisFunctionSet(order);
    }

    std::vector<const Base::BasisFunctionSet*> createVertexBasisFunctionSet3DH1ConformingPyramid(std::size_t order)
    {
        logger.assert_debug(order == 1, "Only linear basis functions have been implemented so far");
        auto result = std::vector<const Base::BasisFunctionSet*>{};
        for(std::size_t i = 0; i < 5; ++i)
        {
            auto next = new Base::BasisFunctionSet(order);
            next->addBasisFunction(new BasisFunction3DVertexPyramid(i));
            result.push_back(next);
        }
        return result;
    }

    std::vector<const Base::OrientedBasisFunctionSet*> createEdgeBasisFunctionSet3DH1ConformingPyramid(std::size_t order)
    {
        logger.assert_debug(order == 1, "Only linear basis functions have been implemented so far");
        auto result = std::vector<const Base::OrientedBasisFunctionSet*>{};
        for(std::size_t i = 0; i < 8; ++i)
        {
            result.push_back(new Base::OrientedBasisFunctionSet(order, 0, i));
            result.push_back(new Base::OrientedBasisFunctionSet(order, 1, i));
        }
        return result;
    }

    std::vector<const Base::OrientedBasisFunctionSet*> createFaceBasisFunctionSet3DH1ConformingPyramid(std::size_t order)
    {
        logger.assert_debug(order == 1, "Only linear basis functions have been implemented so far");
        auto result = std::vector<const Base::OrientedBasisFunctionSet*>{};
        for(std::size_t i = 0; i < 8; ++i)
        {
            result.push_back(new Base::OrientedBasisFunctionSet(order, i, 0));
        }
        for(std::size_t i = 1; i < 5; ++i)
        {
            for(std::size_t j = 0; j < 6; ++j)
            {
                result.push_back(new Base::OrientedBasisFunctionSet(order, i, j));
            }

        }
        return result;
    }

}
