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
#include "Point.hpp"
using std::cout;
using std::endl;

namespace Geometry
{
    class Point;
    
    Point::Point(unsigned int DIM):
        coordinates_(DIM)
    {
    }
    
    Point::Point( CoordTypeT coords[],unsigned int DIM):
        coordinates_(coords,DIM)
    {
    }

    Point::Point(const PointT& other):
        coordinates_(other.coordinates_)
    {
            if (this->size() != other.size())
            {///\bug this should be checked BEFORE data is copied
                std::cout<<"ERROR HANDLER!!!!"<<"Sizes do not coincide."<<endl;
                    //throw exception;
            }
    }
    
    Point::Point(const VectorOfCoordsT& coord):
        coordinates_(coord)
    {
        if (coord.size()!=this->size())
        {///\bug this should be checked BEFORE data is copied
            std::cout<<"ERROR HANDLER!!!!"<<"Sizes do not coincide."<<endl;
                //throw exception;
        }
    }

    typename Point::CoordTypeT
    Point::getCoordinate(IndexT n)const
    {
        if (n<this->size())
        {
            return coordinates_[n];
        }
        else
        {
            std::cout<<"ERROR HANDLER!!!!"<<"Sizes do not coincide."<<endl;
                //throw exception;
            return 0.0;
        }
    }

    typename Point::VectorOfCoordsT
    Point::getCoordinates() const
    {
        //cout << "###################" << coordinates_;
        return coordinates_;
    }

    void
    Point::setCoordinate(IndexT n, const CoordTypeT& coord)
    {
        if (n<this->size())
        {
            coordinates_[n]=coord;
        }
        else
        {
            std::cout<<"ERROR HANDLER!!!!"<<"Sizes do not coincide."<<endl;
                //throw exception;
        }
    }

    void
    Point::setCoordinates(const VectorOfCoordsT& coord)
    {
        if (coord.size()==this->size())
        {
            coordinates_ = coord;
        }
        else
        {
            std::cout<<"ERROR HANDLER!!!!"<<"Sizes do not coincide."<<endl;
                //throw exception;
        }

    }


    typename Point::CoordTypeT&
    Point::operator [] (IndexT n)
    {///\bug no size checking
        return coordinates_[n];
    }
    
    const typename Point::CoordTypeT&
    Point::operator [] (IndexT n)const
    {///\bug no size checking
        return coordinates_[n];
    }
    
    Point&
    Point::operator = (const PointT& rhs)
    {///\bug no size checking
        coordinates_ = rhs.coordinates_;return *this;
    }
};
