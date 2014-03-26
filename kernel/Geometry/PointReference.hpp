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
        
        PointReferenceT  operator* (double right)
                {return PointReferenceT(PointT::coordinates_ * right);}
        
        PointReferenceT  operator* (double right) const
                {return PointReferenceT(PointT::coordinates_ * right);}
        
        PointReferenceT  operator+ (const PointReferenceT& right)
                {return PointReferenceT(PointT::coordinates_  + right.coordinates_);}
        
        PointReferenceT  operator+ (const PointReferenceT& right) const
                {return PointReferenceT(PointT::coordinates_  + right.coordinates_);}
        
        PointReferenceT  operator- (const PointReferenceT& right)
                {return PointReferenceT(PointT::coordinates_  - right.coordinates_);}
        
        PointReferenceT  operator- (const PointReferenceT& right) const
                {return PointReferenceT(PointT::coordinates_  - right.coordinates_);}
        
        PointReferenceT& operator = (const PointReferenceT& rhs)
        {
            this->coordinates_ = rhs.coordinates_;
            return *this;
        }
        
//        friend PointT operator*(const double& left, const PointReferenceT& right){return PointReferenceT(right.coordinates_*left)}
        
    };
};

#endif /* POINTREFERENCE_HPP_ */
