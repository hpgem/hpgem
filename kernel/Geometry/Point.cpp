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
#include <iostream>
#include "Point.h"
using std::cout;
using std::endl;

namespace Geometry
{
    class Point;
    
    Point::Point(std::size_t DIM)
            : coordinates_(DIM)
    {
    }
    
    Point::Point(double coords[], std::size_t DIM)
            : coordinates_(coords, DIM)
    {
    }
    
    Point::Point(const Point& other)
            : coordinates_(other.coordinates_)
    {
        if (this->size() != other.size())
        { ///\bug this should be checked BEFORE data is copied
            std::cout << "ERROR HANDLER!!!!" << "Sizes do not coincide." << endl;
            //throw exception;
        }
    }
    
    Point::Point(const VectorOfCoordsT& coord)
            : coordinates_(coord)
    {
        if (coord.size() != this->size())
        { ///\bug this should be checked BEFORE data is copied
            std::cout << "ERROR HANDLER!!!!" << "Sizes do not coincide." << endl;
            //throw exception;
        }
    }
    
    bool Point::operator ==(const Point& right) const
    {
        return coordinates_ == right.coordinates_;
    }
    
    bool Point::operator!=(const Point& right) const
    {
        return !(*this == right);
    }
    
    bool Point::operator <(const Point& right) const
    {
        return coordinates_ < right.coordinates_;
    }
    
    Point& Point::operator +=(const Point& right)
    {
        coordinates_ += right.coordinates_;
        return *this;
    }
    
    Point& Point::operator -=(const Point& right)
    {
        coordinates_ -= right.coordinates_;
        return *this;
    }
    
    Point& Point::operator *=(double right)
    {
        coordinates_.operator *=(right);
        return *this;
    }
    
    Point Point::operator *(double right) const
    {
        return Point(coordinates_ * right);
    }
    
    Point Point::operator *(double right)
    {
        return Point(coordinates_ * right);
    }
    
    Point Point::operator +(const Point& right) const
    {
        return Point(coordinates_ + right.coordinates_);
    }
    
    Point Point::operator +(const Point& right)
    {
        return Point(coordinates_ + right.coordinates_);
    }
    
    Point Point::operator -(const Point& right) const
    {
        return Point(coordinates_ - right.coordinates_);
    }
    
    Point Point::operator -(const Point& right)
    {
        return Point(coordinates_ - right.coordinates_);
    }
    
    std::size_t Point::size()
    {
        return coordinates_.size();
    }
    
    std::size_t Point::size() const
    {
        return coordinates_.size();
    }
    
    double Point::getCoordinate(std::size_t n) const
    {
        logger.assert(n < size(), "In Point::getCoordinate, entry % is requested while the dimension is %", n, size());
        return coordinates_[n];
    }
    
    const typename Point::VectorOfCoordsT&
    Point::getCoordinates() const
    {
        return coordinates_;
    }
    
    void Point::setCoordinate(std::size_t n, const double& coord)
    {
        logger.assert(n < size(), "In Point::setCoordinate, trying to set entry % while the dimension is %", n, size());
        coordinates_[n] = coord;
    }
    
    void Point::setCoordinates(const VectorOfCoordsT& coord)
    {
        logger.assert(size() == coord.size(), "In Point::setCoordinates, the size of the given vector does not coincide with the dimension of the point.");
        coordinates_ = coord;        
    }
    
    double& Point::operator [](std::size_t n)
    { 
        logger.assert(n < size(), "In Point::operator[], entry % is requested while the dimension is %", n, size());
        return coordinates_[n];
    }
    
    const double& Point::operator [](std::size_t n) const
    {
        logger.assert(n < size(), "In Point::operator[] const, entry % is requested while the dimension is %", n, size());
        return coordinates_[n];
    }
    
    Point&
    Point::operator =(const Point& rhs)
    { ///\bug no size checking
        logger.assert(size() == rhs.size(), "In Point::operator=, the dimension of the given Point does not coincide with the dimension of the point.");
        coordinates_ = rhs.coordinates_;
        return *this;
    }
    
    std::ostream& operator <<(std::ostream& os, const Point& point)
    {
        os << "point={";
        for (std::size_t i = 0; i < point.coordinates_.size(); i++)
        {
            if (i < point.coordinates_.size() - 1)
                os << point.coordinates_[i] << ',';
            else
                os << point.coordinates_[i];
        }
        os << "} ";
        return os;
    }
    
    Point operator *(const double& left, const Point& right)
    {
        return Point(right.coordinates_ * left);
    }
}
;
