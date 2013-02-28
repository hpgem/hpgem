#ifndef PointReference_hpp
#define PointReference_hpp

#include "Point.hpp"

namespace Geometry
{
    template<unsigned int DIM>
    class PointReference: public Point<DIM>
    {

    public:

        typedef double                                  CoordTypeT;
        typedef Point<DIM>                              PointT;
        typedef PointReference<DIM>                     PointReferenceT;
        typedef typename PointT::VectorOfCoordsT        VectorOfCoordsT;
        typedef typename PointT::IndexT                 IndexT;

    public:

        PointReference() {}

        PointReference(const PointT& p): PointT(p) {}

        PointReference(CoordTypeT coords[]): Point<DIM>(coords) {}

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
