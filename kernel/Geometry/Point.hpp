/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, Univesity of Twenete
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef _POINT_HPP
#define _POINT_HPP

#include <list>
#include <iostream>
#include "../LinearAlgebra/NumericalVector.hpp"
using std::cout;
using std::endl;
using namespace LinearAlgebra;

namespace Geometry
{
	class Point
	{
    public:
        
    public:
	    /// Typedefs.
	    typedef double                         CoordTypeT;
        typedef Geometry::Point                PointT;
        typedef NumericalVector                VectorOfCoordsT;
        typedef unsigned int                   IndexT;

	public:
        /// Constructors.
        Point(unsigned int DIM);
        /// Warning!!! This way Point ctr will truncate and take sizeof(dimension) points and will not give any warning. Be sure you took the right dimension.
        Point(CoordTypeT coords[],unsigned int DIM);
        
        Point(const Point& other);
        
        Point(const VectorOfCoordsT& coord);
        
        
        void                setCoordinate(IndexT n, const CoordTypeT& coord);
        void                setCoordinates(const VectorOfCoordsT& coord);
        
        CoordTypeT&         operator [] (IndexT n);
        const CoordTypeT&   operator [] (IndexT n)const;
        
//        CoordTypeT&         operator () (IndexT n);
//        const CoordTypeT&   operator () (IndexT n) const;
        PointT&             operator = (const Point& rhs);
        
        bool                operator== (const Point& right) const {return coordinates_==right.coordinates_;}

        bool                operator== (const Point& right){return coordinates_==right.coordinates_;}
        
        bool                operator< (const Point& right) const {return coordinates_<right.coordinates_;}
        
        PointT&             operator+= (const Point& right){coordinates_+=right.coordinates_; return *this;}

        PointT&             operator-= (const Point& right){coordinates_-=right.coordinates_; return *this;}
        
        PointT&             operator*= (double right){coordinates_.operator *=(right); return *this;}
        
        PointT              operator* (double right){return  PointT(coordinates_ * right);}

        PointT              operator* (double right)const{return  PointT(coordinates_ * right);}
        
        PointT              operator+ (const Point& right){return  PointT(coordinates_ + right.coordinates_);}
        
        PointT              operator+ (const Point& right)const{return  PointT(coordinates_ + right.coordinates_);}
        
        PointT              operator- (const Point& right){return  PointT(coordinates_ - right.coordinates_);}
        
        
        PointT              operator- (const Point& right)const{return  PointT(coordinates_ - right.coordinates_);}
        
        
        unsigned int        size()const{return coordinates_.size();}
        
        unsigned int        size(){return coordinates_.size();}
        

        CoordTypeT          getCoordinate(IndexT n)const;
        VectorOfCoordsT     getCoordinates()const;
                    friend PointT       operator-(const Point& right){return PointT(right * -1.0);}
        
        friend PointT       operator*(const double& left, const PointT& right){return PointT(right.coordinates_*left);}
        
            /// Output routine.
        friend std::ostream& operator<<(std::ostream& os, const PointT& point)
        {
            
                // cout << "Size in ostream="<< point.coordinates_.size()<<endl;
            os <<"point={";

            for (unsigned int i = 0; i < point.coordinates_.size(); i++)
            {
                if (i<point.coordinates_.size()-1)
                    os << point.coordinates_[i]<< ',';
                else
                    os << point.coordinates_[i];
            }
            os <<"} ";
            return os;
        }

	protected:
        VectorOfCoordsT coordinates_;
	};
};

#endif /* defined(_NODE_HPP) */
