//
//  Point.cpp
//  
//
//  Created by Shavarsh Nurijanyan on 2/3/13.
//
//
#include <iostream>
#include "point.hpp"
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
