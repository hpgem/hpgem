#ifndef BaseBasisFunction_hpp
#define BaseBasisFunction_hpp

#include "../Geometry/PointReference.hpp"

namespace Base
{
    /*template <unsigned int dim>
    class BaseBasisFunction;
    
    template <>
    class BaseBasisFunction<1>
    {
    public:
        typedef Geometry::PointReference<1> PointReferenceT;
    public:
        virtual ~BaseBasisFunction() {};

        virtual double eval(const PointReferenceT& p) const = 0;
        virtual void   eval(const PointReferenceT& p, NumericalVector& ret) const
        {
	    ret.resize(1);
	    ret[0]=eval(p);
	}
        virtual double evalDeriv0(const PointReferenceT& p) const = 0;
	virtual void   evalCurl(const PointReferenceT& p, NumericalVector& ret) const {throw "The curl of a scalar valued basisfunction is not implemented. Perhaps you meant evalDeriv?";}
    };

    template <>
    class BaseBasisFunction<2>
    {
    public:
        typedef Geometry::PointReference<2> PointReferenceT;
    public:
        virtual ~BaseBasisFunction() {};

        virtual double eval(const PointReferenceT& p) const = 0;
        virtual void   eval(const PointReferenceT& p, NumericalVector& ret) const
        {
	    ret.resize(1);
	    ret[0]=eval(p);
	}
        virtual double evalDeriv0(const PointReferenceT& p) const = 0;
        virtual double evalDeriv1(const PointReferenceT& p) const = 0;
	virtual void   evalCurl(const PointReferenceT& p, NumericalVector& ret) const {throw "The curl of a scalar valued basisfunction is not implemented. Perhaps you meant evalDeriv?";}
    };

    template <>
    class BaseBasisFunction<3>
    {
    public:
        typedef Geometry::PointReference<3> PointReferenceT;
    public:
        virtual ~BaseBasisFunction() {};
        
        virtual double eval(const PointReferenceT& p) const = 0;
        virtual void   eval(const PointReferenceT& p, NumericalVector& ret) const
        {
	    ret.resize(1);
	    ret[0]=eval(p);
	}
        virtual double evalDeriv0(const PointReferenceT& p) const = 0;
        virtual double evalDeriv1(const PointReferenceT& p) const = 0;
        virtual double evalDeriv2(const PointReferenceT& p) const = 0;
	virtual void   evalCurl(const PointReferenceT& p, NumericalVector& ret) const {throw "The curl of a scalar valued basisfunction is not implemented. Perhaps you meant evalDeriv?";}
    };*/

    class BaseBasisFunction
    {
    public:
        typedef Geometry::PointReference PointReferenceT;
    public:
        virtual ~BaseBasisFunction() {};

        virtual double eval(const PointReferenceT& p) const = 0;
        virtual void   eval(const PointReferenceT& p, NumericalVector& ret) const
        {
	    ret.resize(1);
	    ret[0]=eval(p);
	}
        virtual double evalDeriv0(const PointReferenceT& p) const {throw "The DIMension of your problem is too low to warrant taking a derivative in this direction";};
        virtual double evalDeriv1(const PointReferenceT& p) const {throw "The DIMension of your problem is too low to warrant taking a derivative in this direction";};
        virtual double evalDeriv2(const PointReferenceT& p) const {throw "The DIMension of your problem is too low to warrant taking a derivative in this direction";};
        virtual double evalDeriv3(const PointReferenceT& p) const {throw "The DIMension of your problem is too low to warrant taking a derivative in this direction";};
        virtual void   evalCurl(const PointReferenceT& p, NumericalVector& ret) const {throw "The curl of a scalar valued basisfunction is not implemented. Perhaps you meant evalDeriv?";}
    };
    
};

#endif // BaseBasisFunction_hpp
