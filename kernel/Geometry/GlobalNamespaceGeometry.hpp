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

#ifndef _GlobalNamespaceGeometry_hpp
#define _GlobalNamespaceGeometry_hpp

#include <limits> 

namespace Geometry
{
    const int ZeroD  = 0;
    const int OneD   = 1;
    const int TwoD   = 2;
    const int ThreeD = 3;
    const int FourD  = 4;
    
    const int MaxInteger=std::numeric_limits<int>::max();
    const double SmallerDoubleThanMinimalSizeOfTheMesh=1.0e-10;
    
//    template <unsigned int DIM>
//    class TypedefsInGeometry
//    {
//    public:
//        typedef double                       CoordTypeT;
//        typedef unsigned int                 IndexT;
//        
//        typedef std::vector<NodeT>           VectorOfNodes;
//        typedef std::vector<NodeT* >         VectorOfPtrNodesT;
//        typedef std::vector<CoordTypeT >     VectorOfCoordsT;
//
//        typedef std::list<NodeT >            ListOfNodesT;
//        typedef std::list<NodeT* >           ListOfPtrNodesT;
//
//        typedef PhysicalGeometry<DIM>        PhysicalGeometryT;
//        typedef ReferenceGeometry<DIM>       ReferenceGeometryT;
//
//        typedef int MatrixT; // temp. until we get (if we get) linear algebra.
//    };
};

#endif///_GlobalNamespaceGeometry_hpp
