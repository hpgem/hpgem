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
#ifndef _POINT_HPP
#define _POINT_HPP

#include <iostream>

#include "LinearAlgebra/SmallVector.h"

namespace Geometry
{
    template<std::size_t DIM>
    class Point
    {
    public:
        /// Constructors.
        Point();
        
        ///\deprecated Because this Point ctr will truncate and take sizeof(dimension) points and will not give any warning. If possible, use Point(LinearAlgebra::SmallVector<DIM>&) instead.
        Point(double coords[]);

        Point(const Point& other);

        Point(const LinearAlgebra::SmallVector<DIM>& coord);

        void setCoordinate(std::size_t n, const double& coord);
        void setCoordinates(const LinearAlgebra::SmallVector<DIM>& coord);

        double& operator[](std::size_t n);
        const double& operator[](std::size_t n) const;

        Point& operator=(const Point& rhs);

        bool operator==(const Point& right) const;

        bool operator!=(const Point& right) const;

        bool operator<(const Point& right) const;

        Point& operator+=(const Point& right);

        Point& operator-=(const Point& right);

        Point& operator*=(double right);

        Point operator*(double right);

        double operator*(const Point& right) const
        {
            return coordinates_ * right.coordinates_;
        }

        Point operator*(double right) const;

        Point operator+(const Point& right);

        Point operator+(const Point& right) const;

        Point operator-(const Point& right);

        Point operator-(const Point& right) const;

        std::size_t size() const;

        std::size_t size();

        double getCoordinate(std::size_t n) const;
        const LinearAlgebra::SmallVector<DIM>& getCoordinates() const;

        friend Point operator-(const Point& right)
        {
            return right * -1.0;
        }

    protected:
        LinearAlgebra::SmallVector<DIM> coordinates_;
    };

    template<std::size_t DIM>
    Point<DIM> operator*(const double& left, const Point<DIM>& right);

}

#include "Point_Impl.h"
#endif /* defined(_POINT_HPP) */
