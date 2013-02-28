//
//  Point.cpp
//  
//
//  Created by Shavarsh Nurijanyan on 2/3/13.
//
//
#include <iostream>
using std::cout;
using std::endl;

namespace Geometry
{
    template<unsigned int DIM>
    class Point;
    
    template<unsigned int DIM>
    Point<DIM>::Point():
        coordinates_(DIM)
    {
    }
    
    template<unsigned int DIM>
    Point<DIM>::Point( CoordTypeT coords[]):
        coordinates_(coords, DIM)
    {
    }

    template<unsigned int DIM>
    Point<DIM>::Point(const PointT& other):
        coordinates_(other.coordinates_)
    {
            if (this->size() != other.size())
            {
                std::cout<<"ERROR HANDLER!!!!"<<"Sizes do not coincide."<<endl;
                    //throw exception;
            }
    }
    
    template<unsigned int DIM>
    Point<DIM>::Point(const VectorOfCoordsT& coord):
        coordinates_(coord)
    {
        if (coord.size()!=this->size())
        {
            std::cout<<"ERROR HANDLER!!!!"<<"Sizes do not coincide."<<endl;
                //throw exception;
        }
    }

    template<unsigned int DIM>
    typename Point<DIM>::CoordTypeT
    Point<DIM>::getCoordinate(IndexT n)const
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

    template<unsigned int DIM>
    typename Point<DIM>::VectorOfCoordsT
    Point<DIM>::getCoordinates() const
    {
        //cout << "###################" << coordinates_;
        return coordinates_;
    }

    template<unsigned int DIM>
    void
    Point<DIM>::setCoordinate(IndexT n, const CoordTypeT& coord)
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

    template<unsigned int DIM>
    void
    Point<DIM>::setCoordinates(const VectorOfCoordsT& coord)
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


    template<unsigned int DIM>
    typename Point<DIM>::CoordTypeT&
    Point<DIM>::operator [] (IndexT n)
    {
        return coordinates_[n];
    }
    
    template<unsigned int DIM>
    const typename Point<DIM>::CoordTypeT&
    Point<DIM>::operator [] (IndexT n)const
    {
        return coordinates_[n];
    }
    
    template<unsigned int DIM>
    Point<DIM>&
    Point<DIM>::operator = (const PointT& rhs)
    {
        coordinates_ = rhs.coordinates_;return *this;
    }
};
