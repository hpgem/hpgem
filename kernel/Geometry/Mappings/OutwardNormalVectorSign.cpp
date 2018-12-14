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
#include "OutwardNormalVectorSign.h"
#include "MappingToRefLineToSquare.h"
#include "MappingToRefLineToTriangle.h"
#include "MappingToRefSquareToSquare.h"
#include "MappingToRefPointToLine.h"
#include "MappingToRefCubeToCube.h"
#include "MappingToRefSquareToCube.h"
#include "MappingToRefTriangleToTetrahedron.h"
#include "MappingToRefFaceToTriangularPrism.h"
#include "MappingToRefFaceToPyramid.h"
#include <typeinfo>
#include <limits>

namespace Geometry
{
    double OutwardNormalVectorSign(const MappingReferenceToReference<1> *const map)
    {
        logger.assert_debug(map != nullptr, "Invalid mapping passed");
        if (typeid(*map) == typeid(const MappingToRefPointToLine0) ||
            typeid(*map) == typeid(const MappingToRefLineToTriangle0) ||
            typeid(*map) == typeid(const MappingToRefLineToTriangle2) ||
            typeid(*map) == typeid(const MappingToRefLineToSquare0) ||
            typeid(*map) == typeid(const MappingToRefLineToSquare2) ||
            typeid(*map) == typeid(const MappingToRefSquareToCube0) ||
            typeid(*map) == typeid(const MappingToRefSquareToCube2) ||
            typeid(*map) == typeid(const MappingToRefSquareToCube4))
        {
            return -1;
        }
        else if (typeid(*map) == typeid(const MappingToRefPointToLine1) ||
                 typeid(*map) == typeid(const MappingToRefLineToTriangle1) ||
                 typeid(*map) == typeid(const MappingToRefLineToSquare1) ||
                 typeid(*map) == typeid(const MappingToRefLineToSquare3) ||
                 typeid(*map) == typeid(const MappingToRefSquareToCube1) ||
                 typeid(*map) == typeid(const MappingToRefSquareToCube3) ||
                 typeid(*map) == typeid(const MappingToRefSquareToCube5) ||
                 typeid(*map) == typeid(const MappingToRefTriangleToTetrahedron0) ||
                 typeid(*map) == typeid(const MappingToRefTriangleToTetrahedron1) ||
                 typeid(*map) == typeid(const MappingToRefTriangleToTetrahedron2) ||
                 typeid(*map) == typeid(const MappingToRefTriangleToTetrahedron3) ||
                 typeid(*map) == typeid(const MappingToRefFaceToTriangularPrism0) ||
                 typeid(*map) == typeid(const MappingToRefFaceToTriangularPrism1) ||
                 typeid(*map) == typeid(const MappingToRefFaceToTriangularPrism2) ||
                 typeid(*map) == typeid(const MappingToRefFaceToTriangularPrism3) ||
                 typeid(*map) == typeid(const MappingToRefFaceToTriangularPrism4) ||
                 typeid(*map) == typeid(const MappingToRefFaceToPyramid0) ||
                 typeid(*map) == typeid(const MappingToRefFaceToPyramid1) ||
                 typeid(*map) == typeid(const MappingToRefFaceToPyramid2) ||
                 typeid(*map) == typeid(const MappingToRefFaceToPyramid3) ||
                 typeid(*map) == typeid(const MappingToRefFaceToPyramid4))
        {
            return 1;
        }
        else
        {
            logger(FATAL, "Face to Element mapping not known for given case. \n");
            return std::numeric_limits<double>::quiet_NaN();
        }
    }
}
