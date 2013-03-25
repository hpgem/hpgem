#ifndef BaseBasisFunction_hpp
#define BaseBasisFunction_hpp

#include "../Geometry/PointReference.hpp"

namespace Base
{
    template <unsigned int dim>
    class BaseBasisFunction;
    
    template <>
    class BaseBasisFunction<1>
    {
    public:
        typedef Geometry::PointReference<1> PointReferenceT;
    public:
        virtual ~BaseBasisFunction() {};

        virtual double eval(const PointReferenceT& p) const = 0;
        virtual double evalDeriv0(const PointReferenceT& p) const = 0;
    };

    template <>
    class BaseBasisFunction<2>
    {
    public:
        typedef Geometry::PointReference<2> PointReferenceT;
    public:
        virtual ~BaseBasisFunction() {};

        virtual double eval(const PointReferenceT& p) const = 0;
        virtual double evalDeriv0(const PointReferenceT& p) const = 0;
        virtual double evalDeriv1(const PointReferenceT& p) const = 0;
    };

    template <>
    class BaseBasisFunction<3>
    {
    public:
        typedef Geometry::PointReference<3> PointReferenceT;
    public:
        virtual ~BaseBasisFunction() {};
        
        virtual double eval(const PointReferenceT& p) const = 0;
        virtual double evalDeriv0(const PointReferenceT& p) const = 0;
        virtual double evalDeriv1(const PointReferenceT& p) const = 0;
        virtual double evalDeriv2(const PointReferenceT& p) const = 0;
    };

    template <>
    class BaseBasisFunction<4>
    {
    public:
        typedef Geometry::PointReference<4> PointReferenceT;
    public:
        virtual ~BaseBasisFunction() {};

        virtual double eval(const PointReferenceT& p) const = 0;
        virtual double evalDeriv0(const PointReferenceT& p) const = 0;
        virtual double evalDeriv1(const PointReferenceT& p) const = 0;
        virtual double evalDeriv2(const PointReferenceT& p) const = 0;
        virtual double evalDeriv3(const PointReferenceT& p) const = 0;
};
    
};

#endif // BaseBasisFunction_hpp
