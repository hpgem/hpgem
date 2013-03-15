#ifndef BasisFunctionsCollection_A_hpp
#define BasisFunctionsCollection_A_hpp

#include "BaseBasisFunction.hpp"
#include <vector>

namespace Base
{

    //---------------------------------------------------------------------------
    //------------------------------------------------ basis collection for dim=1
    
    //! Basis function on 1D: u(x) = 1
    struct Basis_A0_1D : public Base::BaseBasisFunction<1>
    {
        virtual double Eval(const PointReferenceT& p) const         { return 1; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 1D: u(x) = x
    struct Basis_A1_1D : public Base::BaseBasisFunction<1>
    {
        virtual double Eval(const PointReferenceT& p) const         { return p[0]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 1.; }
    };

    //! Basis function on 1D: u(x) = x^2
    struct Basis_A2_1D : public Base::BaseBasisFunction<1>
    {
        virtual double Eval(const PointReferenceT& p) const         { return p[0]*p[0]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 2.*p[0]; }
    };

    //! Basis function on 1D: u(x) = x^3
    struct Basis_A3_1D : public Base::BaseBasisFunction<1>
    {
        virtual double Eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[0]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 3.*p[0]*p[0]; }
    };

    //---------------------------------------------------------------------------
    //------------------------------------------------ basis collection for dim=2
    
    //! Basis function on 2D: u(x,y) = 1
    struct Basis_A0_2D : public Base::BaseBasisFunction<2>
    {
        virtual double Eval(const PointReferenceT& p) const         { return 1; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 2D: u(x,y) = x
    struct Basis_A1_2D : public Base::BaseBasisFunction<2>
    {
        virtual double Eval(const PointReferenceT& p) const         { return p[0]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 1.; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 2D: u(x,y) = y
    struct Basis_A2_2D : public Base::BaseBasisFunction<2>
    {
        virtual double Eval(const PointReferenceT& p) const         { return p[1]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 1.; }
    };

    //! Basis function on 2D: u(x,y) = x*y
    struct Basis_A3_2D : public Base::BaseBasisFunction<2>
    {
        virtual double Eval(const PointReferenceT& p) const         { return p[0]*p[1]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return p[1]; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return p[0]; }
    };

    //! Basis function on 2D: u(x,y) = x^2
    struct Basis_A4_2D : public Base::BaseBasisFunction<2>
    {
        virtual double Eval(const PointReferenceT& p) const         { return p[0]*p[0]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 2.*p[0]; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 2D: u(x,y) = y^2
    struct Basis_A5_2D : public Base::BaseBasisFunction<2>
    {
        virtual double Eval(const PointReferenceT& p) const         { return p[1]*p[1]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 2.*p[1]; }
    };

    //! Basis function on 2D: u(x,y) = y*x^2
    struct Basis_A6_2D : public Base::BaseBasisFunction<2>
    {
        virtual double Eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[1]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 2.*p[0]*p[1]; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return p[0]*p[0]; }
    };

    //! Basis function on 2D: u(x,y) = x*y^2
    struct Basis_A7_2D : public Base::BaseBasisFunction<2>
    {
        virtual double Eval(const PointReferenceT& p) const         { return p[0]*p[1]*p[1]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return p[1]*p[1]; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 2.*p[0]*p[1]; }
    };

    //! Basis function on 2D: u(x,y) = x^3
    struct Basis_A8_2D : public Base::BaseBasisFunction<2>
    {
        virtual double Eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[0]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 3.*p[0]*p[0]; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 2D: u(x,y) = y^3
    struct Basis_A9_2D : public Base::BaseBasisFunction<2>
    {
        virtual double Eval(const PointReferenceT& p) const         { return p[1]*p[1]*p[1]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 3.*p[1]*p[1]; }
    };


    //---------------------------------------------------------------------------
    //------------------------------------------------ basis collection for dim=3
    
    //! Basis function on 3D: u(x,y,z) = 1
    struct Basis_A0_3D : public Base::BaseBasisFunction<3>
    {
        virtual double Eval(const PointReferenceT& p) const         { return 1; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = x
    struct Basis_A1_3D : public Base::BaseBasisFunction<3>
    {
        virtual double Eval(const PointReferenceT& p) const         { return p[0]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 1.; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = y
    struct Basis_A2_3D : public Base::BaseBasisFunction<3>
    {
        virtual double Eval(const PointReferenceT& p) const         { return p[1]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 1.; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = z
    struct Basis_A3_3D : public Base::BaseBasisFunction<3>
    {
        virtual double Eval(const PointReferenceT& p) const         { return p[2]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return 1.; }
    };

    //! Basis function on 3D: u(x,y,z) = x*y
    struct Basis_A4_3D: public Base::BaseBasisFunction<3>
    {
        virtual double Eval(const PointReferenceT& p) const         { return p[0]*p[1]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return p[1]; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return p[0]; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = x*z
    struct Basis_A5_3D: public Base::BaseBasisFunction<3>
    {
        virtual double Eval(const PointReferenceT& p) const         { return p[0]*p[2]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return p[2]; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return p[0]; }
    };

    //! Basis function on 3D: u(x,y,z) = y*z
    struct Basis_A6_3D: public Base::BaseBasisFunction<3>
    {
        virtual double Eval(const PointReferenceT& p) const         { return p[1]*p[2]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return p[2]; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return p[1]; }
    };

    //! Basis function on 3D: u(x,y,z) = x^2
    struct Basis_A7_3D: public Base::BaseBasisFunction<3>
    {
        virtual double Eval(const PointReferenceT& p) const         { return p[0]*p[0]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 2.*p[0]; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = y^2
    struct Basis_A8_3D: public Base::BaseBasisFunction<3>
    {
        virtual double Eval(const PointReferenceT& p) const         { return p[1]*p[1]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 2.*p[1]; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = z^2
    struct Basis_A9_3D: public Base::BaseBasisFunction<3>
    {
        virtual double Eval(const PointReferenceT& p) const         { return p[2]*p[2]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return 2.*p[2]; }
    };

    //! Basis function on 3D: u(x,y,z) = x*y*z
    struct Basis_A10_3D: public Base::BaseBasisFunction<3>
    {
        virtual double Eval(const PointReferenceT& p) const         { return p[0]*p[1]*p[2]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return p[1]*p[2]; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return p[0]*p[2]; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return p[0]*p[1]; }
    };

    //! Basis function on 3D: u(x,y,z) = y*x^2
    struct Basis_A11_3D: public Base::BaseBasisFunction<3>
    {
        virtual double Eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[1]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 2.*p[0]*p[1]; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return p[0]*p[0]; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = z*x^2
    struct Basis_A12_3D: public Base::BaseBasisFunction<3>
    {
        virtual double Eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[2]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 2.*p[0]*p[2]; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return p[0]*p[0]; }
    };

    //! Basis function on 3D: u(x,y,z) = x*y^2
    struct Basis_A13_3D: public Base::BaseBasisFunction<3>
    {
        virtual double Eval(const PointReferenceT& p) const         { return p[0]*p[1]*p[1]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return p[1]*p[1]; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 2.*p[0]*p[1]; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = z*y^2
    struct Basis_A14_3D: public Base::BaseBasisFunction<3>
    {
        virtual double Eval(const PointReferenceT& p) const         { return p[2]*p[1]*p[1]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 2.*p[1]*p[2]; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return p[1]*p[1]; }
    };

    //! Basis function on 3D: u(x,y,z) = x*z^2
    struct Basis_A15_3D: public Base::BaseBasisFunction<3>
    {
        virtual double Eval(const PointReferenceT& p) const         { return p[0]*p[2]*p[2]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return p[2]*p[2]; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return 2.*p[0]*p[2]; }
    };

    //! Basis function on 3D: u(x,y,z) = y*z^2
    struct Basis_A16_3D: public Base::BaseBasisFunction<3>
    {
        virtual double Eval(const PointReferenceT& p) const         { return p[1]*p[2]*p[2]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return p[2]*p[2]; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return 2.*p[1]*p[2]; }
    };

    //! Basis function on 3D: u(x,y,z) = x^3
    struct Basis_A17_3D: public Base::BaseBasisFunction<3>
    {
        virtual double Eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[0]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 3.*p[0]*p[0]; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = y^3
    struct Basis_A18_3D: public Base::BaseBasisFunction<3>
    {
        virtual double Eval(const PointReferenceT& p) const         { return p[1]*p[1]*p[1]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 3.*p[1]*p[1]; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = z^3
    struct Basis_A19_3D: public Base::BaseBasisFunction<3>
    {
        virtual double Eval(const PointReferenceT& p) const         { return p[2]*p[2]*p[2]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return 3.*p[2]*p[2]; }
    };
    
};
#endif  // BasisFunctionsCollection_A_hpp
