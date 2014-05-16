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

#ifndef PointReference_hpp
#define PointReference_hpp

#include "Point.hpp"

namespace Geometry
{
    class PointReference: public Point
    {

    public:

        typedef double                                  CoordTypeT;
        typedef Point                                   PointT;
        typedef PointReference                          PointReferenceT;
        typedef PointT::VectorOfCoordsT        VectorOfCoordsT;
        typedef PointT::IndexT                 IndexT;

    public:

        PointReference(unsigned int DIM):PointT(DIM) {}

        PointReference(const PointT& p): PointT(p) {}

        PointReference(CoordTypeT coords[],unsigned int DIM): Point(coords,DIM) {}

        PointReference(const VectorOfCoordsT& coord): PointT(coord){}
        
        PointReference operator *(double right);
        
        PointReference operator *(double right) const;
        
        PointReference operator +(const PointReferenceT& right);
        
        PointReference operator +(const PointReferenceT& right) const;
        
        PointReference operator -(const PointReferenceT& right);
        
        PointReference operator -(const PointReferenceT& right) const;
        
        PointReference& operator =(const PointReferenceT& rhs);
        
//        friend PointT operator*(const double& left, const PointReferenceT& right){return PointReferenceT(right.coordinates_*left)}
        
    };

}
;

#endif /* POINTREFERENCE_HPP_ */
