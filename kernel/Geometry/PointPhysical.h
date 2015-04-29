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

#ifndef POINTPHYSICAL_HPP_
#define POINTPHYSICAL_HPP_

#include "Point.h"
#include <complex>
namespace Geometry
{
    class PointPhysical : public Point
    {
        
    public:
        
        PointPhysical(std::size_t DIM)
                : Point(DIM)
        {
        }
        
        PointPhysical(const PointPhysical& p)
                : Point(p)
        {
        }

        explicit PointPhysical(const Point& p)
                : Point(p)
        {
        }
        
        PointPhysical(const VectorOfCoordsT& coord)
                : Point(coord)
        {
        }
        
        PointPhysical operator*(double right) const
        {
            return PointPhysical(Point::coordinates_ * right);
        }
        
        PointPhysical operator/(double right) const
        {
            return PointPhysical(Point::coordinates_ / right);
        }
        
        //please note that for type-safety this function cannot be removed in favor
        //of the Point::operator+
        PointPhysical operator+(const PointPhysical& right) const
        {
            logger.assert(size()==right.size(), "The sizes of the points do not match");
            return PointPhysical(Point::coordinates_ + right.coordinates_);
        }
        
        PointPhysical operator-(const PointPhysical& right) const
        {
            logger.assert(size()==right.size(), "The sizes of the points do not match");
            return PointPhysical(Point::coordinates_ - right.coordinates_);
        }
        
        PointPhysical operator-() const
        {
            return *this * -1.;
        }
        
        PointPhysical& operator=(const PointPhysical& right)
        {
            Point::coordinates_ = right.coordinates_;
            return *this;
        }
        
        void axpy(const double& alpha, const PointPhysical& x)
        {
            logger.assert(size()==x.size(), "The sizes of the points do not match");
            coordinates_.axpy(alpha, x.coordinates_);
        }
        
#ifdef HPGEM_USE_COMPLEX_PETSC
        
        const std::complex<double>* data() const
        {   
            static std::complex<double>* new_Data;

            for (std::size_t i = 0; i < PointT::coordinates_.size(); i++)
            {   
                new_Data[i] = PointT::coordinates_.data()[i];
            }
            return new_Data;
        }
#else
        
        const double* data() const
        {
            return coordinates_.data();
        }
        
#endif
    };
    
    PointPhysical operator*(double left, const PointPhysical& right);
}

#endif /* POINTPHYSICAL_HPP_ */
