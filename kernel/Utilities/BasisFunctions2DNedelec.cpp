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

#include "BasisFunctions2DNedelec.h"
#include "helperFunctions.h"
#include "Base/BasisFunctionSet.h"
#include "Geometry/ReferenceTriangle.h"
#include "Geometry/PointReference.h"

//only uses the constant basis functions
#include "BasisFunctionsCollection_A.h"


namespace Utilities
{
    namespace
    {
        LinearAlgebra::SmallVector<2> baricentricDeriv(std::size_t node)
        {
            LinearAlgebra::SmallVector<2> ret;
            if (node == 0)
            {
                ret[0] = -1;
                ret[1] = -1;
            }
            else
            {
                //clear the return vector so we don't return trash
                ret[0] = 0;
                ret[1] = 0;
                ret[node - 1] = 1;
            }
            return ret;
        }
    }    

    BasisCurlEdgeNedelec2D::BasisCurlEdgeNedelec2D(std::size_t degree1, std::size_t degree2, std::size_t localFirstVertex, std::size_t localSecondVertex)
            : deg1(degree1), deg2(degree2), i(localFirstVertex), j(localSecondVertex)
    {
        logger.assert_debug(deg1 < 1 && deg2 < 1, "2D Nedelec basis functions only implemented for p = 1");
        logger.assert_debug(i < 3 && j < 3, "A triangle only has 3 nodes");
    }
    
    void BasisCurlEdgeNedelec2D::eval(const Geometry::PointReference<2>& p, LinearAlgebra::SmallVector<2>& ret) const 
    {
        LinearAlgebra::SmallVector<2> dummy;

        ret = baricentricDeriv(i);
        dummy = baricentricDeriv(j);

        double valI(baricentric_2D(i, p)),
	          valJ(baricentric_2D(j, p));

        ret*=valJ;
        dummy*=valI;
        ret-=dummy;
    }

    LinearAlgebra::SmallVector<2> BasisCurlEdgeNedelec2D::evalCurl(const Geometry::PointReference<2>& p) const 
    {
        LinearAlgebra::SmallVector<2> dummy,dummy2, ret;

        dummy = baricentricDeriv(i);
        dummy2 = baricentricDeriv(j);

        double valI(baricentric_2D(i,p)),
               valJ(baricentric_2D(j,p));

        ret[0] = -2*dummy2[1]*dummy[0] + 2*dummy[1]*dummy2[0]; 
        ret[1] = 0.0;
        return ret;
    }
    
    Base::BasisFunctionSet* createDGBasisFunctionSet2DNedelec(std::size_t order)
    {
        logger.assert_debug(order < 2, "2D Nedelec basis functions only implemented for p = 1");
        Base::BasisFunctionSet* bFset = new Base::BasisFunctionSet(order);
        bFset->addBasisFunction(new BasisCurlEdgeNedelec2D(0,0,0,1));
        bFset->addBasisFunction(new BasisCurlEdgeNedelec2D(0,0,0,2));
        bFset->addBasisFunction(new BasisCurlEdgeNedelec2D(0,0,1,2));
        return bFset;
    }
}
