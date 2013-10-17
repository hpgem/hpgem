
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
	template <unsigned int DIM>
	class Point
	{
    public:
        
    public:
	    /// Typedefs.
	    typedef double                         CoordTypeT;
        typedef Geometry::Point<DIM>           PointT;
        typedef NumericalVector                VectorOfCoordsT;
        typedef unsigned int                   IndexT;

	public:
        /// Constructors.
        Point();
        /// Warning!!! This way Point ctr will truncate and take sizeof(dimension) points and will not give any warning. Be sure you took the right dimension.
        Point(CoordTypeT coords[]);
        
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
#include "Point_Impl.hpp"
#endif /* defined(_NODE_HPP) */
