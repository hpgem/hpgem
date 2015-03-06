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
#include "OutwardNormalVectorSign.hpp"
#include "MappingToRefLineToSquare.hpp"
#include "MappingToRefLineToTriangle.hpp"
#include "MappingToRefSquareToSquare.hpp"
#include "MappingToRefPointToLine.hpp"
#include "MappingToRefCubeToCube.hpp"
#include "MappingToRefSquareToCube.hpp"
#include "MappingToRefTriangleToTetrahedron.hpp"
    //#include "MappingToRefFaceToPyramid.hpp"
#include <typeinfo>

namespace Geometry
{

    // Added by M.T. Julianto, Feb 16, 2010
    double OutwardNormalVectorSign(const MappingReferenceToReference* const map)
    {
        if(typeid(*map)==typeid(const MappingToRefPointToLine0) ||
           typeid(*map)==typeid(const MappingToRefLineToTriangle0) ||
           typeid(*map)==typeid(const MappingToRefLineToTriangle2) ||
           typeid(*map)==typeid(const MappingToRefLineToSquare0) ||
           typeid(*map)==typeid(const MappingToRefLineToSquare2) ||
           typeid(*map)==typeid(const MappingToRefSquareToCube0) ||
           typeid(*map)==typeid(const MappingToRefSquareToCube2) ||
           typeid(*map)==typeid(const MappingToRefSquareToCube4))
        {
            return -1;
        }else if(typeid(*map)==typeid(const MappingToRefPointToLine1) ||
                 typeid(*map)==typeid(const MappingToRefLineToTriangle1) ||
                 typeid(*map)==typeid(const MappingToRefLineToSquare1) ||
                 typeid(*map)==typeid(const MappingToRefLineToSquare3) ||
                 typeid(*map)==typeid(const MappingToRefSquareToCube1) ||
                 typeid(*map)==typeid(const MappingToRefSquareToCube3) ||
                 typeid(*map)==typeid(const MappingToRefSquareToCube5) ||
                 typeid(*map)==typeid(const MappingToRefTriangleToTetrahedron0) ||
                 typeid(*map)==typeid(const MappingToRefTriangleToTetrahedron1) ||
                 typeid(*map)==typeid(const MappingToRefTriangleToTetrahedron2) ||
                 typeid(*map)==typeid(const MappingToRefTriangleToTetrahedron3))
        {
            return 1;
        }else{
                throw "outwardNormalVectorSign - don't know this face2el-map";
        }
        
        /*if (dynamic_cast<const MappingToRefPointToLine0*>(map) ||
        	dynamic_cast<const MappingToRefLineToTriangle0*>(map) ||
	        dynamic_cast<const MappingToRefLineToTriangle2*>(map) ||
	        dynamic_cast<const MappingToRefLineToSquare0*>(map)   ||
	        dynamic_cast<const MappingToRefLineToSquare2*>(map) ||
        	dynamic_cast<const MappingToRefSquareToCube0*>(map) ||
	        dynamic_cast<const MappingToRefSquareToCube2*>(map) ||
	        dynamic_cast<const MappingToRefSquareToCube4*>(map))
        {
            return -1.;
        }
        else
        {
            if (dynamic_cast<const MappingToRefPointToLine1*>(map) ||
            	dynamic_cast<const MappingToRefLineToTriangle1*>(map) ||
                dynamic_cast<const MappingToRefLineToSquare1*>(map)   ||
                dynamic_cast<const MappingToRefLineToSquare3*>(map) ||
                dynamic_cast<const MappingToRefSquareToCube1*>(map) ||
				dynamic_cast<const MappingToRefSquareToCube3*>(map) ||
				dynamic_cast<const MappingToRefSquareToCube5*>(map) ||
				dynamic_cast<const MappingToRefTriangleToTetrahedron0*>(map) ||
				dynamic_cast<const MappingToRefTriangleToTetrahedron1*>(map) ||
				dynamic_cast<const MappingToRefTriangleToTetrahedron2*>(map) ||
				dynamic_cast<const MappingToRefTriangleToTetrahedron3*>(map)
                //                dynamic_cast<const MappingToRefFaceToPyramid0*>(map) ||
                //                dynamic_cast<const MappingToRefFaceToPyramid1*>(map) ||
                //                dynamic_cast<const MappingToRefFaceToPyramid2*>(map) ||
                //                dynamic_cast<const MappingToRefFaceToPyramid3*>(map) ||
                //                dynamic_cast<const MappingToRefFaceToPyramid4*>(map) ||
                //                dynamic_cast<const MappingToRefFaceToTriangularPrism0*>(map) ||
                //                dynamic_cast<const MappingToRefFaceToTriangularPrism1*>(map) ||
                //                dynamic_cast<const MappingToRefFaceToTriangularPrism2*>(map) ||
                //                dynamic_cast<const MappingToRefFaceToTriangularPrism3*>(map) ||
                //                dynamic_cast<const MappingToRefFaceToTriangularPrism4*>(map) ||
                //                dynamic_cast<const MappingToRefFaceToTetrahedron0*>(map) ||
                //                dynamic_cast<const MappingToRefFaceToTetrahedron1*>(map) ||
                //                dynamic_cast<const MappingToRefFaceToTetrahedron2*>(map) ||
                //                dynamic_cast<const MappingToRefFaceToTetrahedron3*>(map)
                		         )
            {
                return +1.;
            }
            else
            {
                throw "outwardNormalVectorSign - don't know this face2el-map";
            }
        }*/
    }

	// Added by Vijaya for 4D Normal vector.
    /*template <>
    double OutwardNormalVectorSign<4>(
	const MappingReferenceToReference<3, 4>* const map)
    {
        if (dynamic_cast<const MappingCubeToHyperCube0*>(map) || // from ReferenceHyperCube
            dynamic_cast<const MappingCubeToHyperCube2*>(map) ||
            dynamic_cast<const MappingCubeToHyperCube4*>(map) ||
            dynamic_cast<const MappingCubeToHyperCube6*>(map))
        {
            return -1.;
        }
        else
        {
            if (dynamic_cast<const MappingCubeToHyperCube1*>(map) ||
                dynamic_cast<const MappingCubeToHyperCube3*>(map) ||
                dynamic_cast<const MappingCubeToHyperCube5*>(map) ||
                dynamic_cast<const MappingCubeToHyperCube7*>(map))
            {
                return +1.;
            }
            else
            {
                throw "outwardNormalVectorSign<3> - don't know this face2el-map";
            }
        }
        return 1;
    }*/
}
