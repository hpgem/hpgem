#ifndef BasisFunctionsCollection_A_hpp
#define BasisFunctionsCollection_A_hpp

#include "BaseBasisFunction.hpp"
#include <vector>

namespace Base
{

    //---------------------------------------------------------------------------
    //------------------------------------------------ basis collection for dim=1
    
    //! Basis function on 1D: u(x) = 1
    struct Basis_A0_1D : public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return 1; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 1D: u(x) = x
    struct Basis_A1_1D : public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 1.; }
    };

    //! Basis function on 1D: u(x) = x^2
    struct Basis_A2_1D : public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[0]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 2.*p[0]; }
    };

    //! Basis function on 1D: u(x) = x^3
    struct Basis_A3_1D : public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[0]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 3.*p[0]*p[0]; }
    };

    //! Basis function on 1D: u(x) = x^4
    struct Basis_A4_1D : public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[0]*p[0]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 4.*p[0]*p[0]*p[0]; }
    };

    //! Basis function on 1D: u(x) = x^5
    struct Basis_A5_1D : public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[0]*p[0]*p[0]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 5.*p[0]*p[0]*p[0]*p[0]; }
    };

    //---------------------------------------------------------------------------
    //------------------------------------------------ basis collection for dim=2
    
    //! Basis function on 2D: u(x,y) = 1
    struct Basis_A0_2D : public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return 1; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 2D: u(x,y) = x
    struct Basis_A1_2D : public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 1.; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 2D: u(x,y) = y
    struct Basis_A2_2D : public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[1]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 1.; }
    };

    //! Basis function on 2D: u(x,y) = x*y
    struct Basis_A3_2D : public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[1]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return p[1]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return p[0]; }
    };

    //! Basis function on 2D: u(x,y) = x^2
    struct Basis_A4_2D : public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[0]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 2.*p[0]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 2D: u(x,y) = y^2
    struct Basis_A5_2D : public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[1]*p[1]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 2.*p[1]; }
    };

    //! Basis function on 2D: u(x,y) = y*x^2
    struct Basis_A6_2D : public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[1]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 2.*p[0]*p[1]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return p[0]*p[0]; }
    };

    //! Basis function on 2D: u(x,y) = x*y^2
    struct Basis_A7_2D : public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[1]*p[1]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return p[1]*p[1]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 2.*p[0]*p[1]; }
    };

    //! Basis function on 2D: u(x,y) = x^3
    struct Basis_A8_2D : public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[0]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 3.*p[0]*p[0]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 2D: u(x,y) = y^3
    struct Basis_A9_2D : public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[1]*p[1]*p[1]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 3.*p[1]*p[1]; }
    };

    //! Basis function on 2D: u(x,y) = x^3*y
    struct Basis_A10_2D : public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[0]*p[1]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 3.*p[0]*p[0]*p[1]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return p[0]*p[0]*p[0]; }
    };

    //! Basis function on 2D: u(x,y) = x^2*y^2
    struct Basis_A11_2D : public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[1]*p[1]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 2.*p[0]*p[1]*p[1]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 2.*p[0]*p[0]*p[1]; }
    };

    //! Basis function on 2D: u(x,y) = x*y^3
    struct Basis_A12_2D : public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[1]*p[1]*p[1]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return p[1]*p[1]*p[1]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 3.*p[0]*p[1]*p[1]; }
    };

    //! Basis function on 2D: u(x,y) = x^4
    struct Basis_A13_2D : public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[0]*p[0]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 4.*p[0]*p[0]*p[0]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 2D: u(x,y) = y^4
    struct Basis_A14_2D : public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[1]*p[1]*p[1]*p[1]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 4.*p[1]*p[1]*p[1]; }
    };

    //! Basis function on 2D: u(x,y) = x^4*y
    struct Basis_A15_2D : public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[0]*p[0]*p[1]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 4.*p[0]*p[0]*p[0]*p[1]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return p[0]*p[0]*p[0]*p[0]; }
    };

    //! Basis function on 2D: u(x,y) = x^3*y^2
    struct Basis_A16_2D : public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[0]*p[1]*p[1]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 3.*p[0]*p[0]*p[1]*p[1]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 2.*p[0]*p[0]*p[0]*p[1]; }
    };

    //! Basis function on 2D: u(x,y) = x^2*y^3
    struct Basis_A17_2D : public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[1]*p[1]*p[1]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 2.*p[0]*p[1]*p[1]*p[1]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 3.*p[0]*p[0]*p[1]*p[1]; }
    };

    //! Basis function on 2D: u(x,y) = x*y^4
    struct Basis_A18_2D : public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[1]*p[1]*p[1]*p[1]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return p[1]*p[1]*p[1]*p[1]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 4.*p[0]*p[1]*p[1]*p[1]; }
    };

    //! Basis function on 2D: u(x,y) = x^5
    struct Basis_A19_2D : public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[0]*p[0]*p[0]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 5.*p[0]*p[0]*p[0]*p[0]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 2D: u(x,y) = y^5
    struct Basis_A20_2D : public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[1]*p[1]*p[1]*p[1]*p[1]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 5.*p[1]*p[1]*p[1]*p[1]; }
    };


    //---------------------------------------------------------------------------
    //------------------------------------------------ basis collection for dim=3
    
    //! Basis function on 3D: u(x,y,z) = 1
    struct Basis_A0_3D : public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return 1; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = x
    struct Basis_A1_3D : public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 1.; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = y
    struct Basis_A2_3D : public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[1]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 1.; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = z
    struct Basis_A3_3D : public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[2]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 1.; }
    };

    //! Basis function on 3D: u(x,y,z) = x*y
    struct Basis_A4_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[1]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return p[1]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return p[0]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = x*z
    struct Basis_A5_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[2]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return p[2]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return p[0]; }
    };

    //! Basis function on 3D: u(x,y,z) = y*z
    struct Basis_A6_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[1]*p[2]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return p[2]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return p[1]; }
    };
    //! Basis function on 3D: u(x,y,z) = x*y*z
    struct Basis_A7_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[1]*p[2]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return p[1]*p[2]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return p[0]*p[2]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return p[0]*p[1]; }
    };
    //! Basis function on 3D: u(x,y,z) = x^2
    struct Basis_A8_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[0]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 2.*p[0]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = y^2
    struct Basis_A9_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[1]*p[1]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 2.*p[1]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = z^2
    struct Basis_A10_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[2]*p[2]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 2.*p[2]; }
    };

    


    //! Basis function on 3D: u(x,y,z) = y*x^2
    struct Basis_A11_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[1]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 2.*p[0]*p[1]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return p[0]*p[0]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = z*x^2
    struct Basis_A12_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[2]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 2.*p[0]*p[2]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return p[0]*p[0]; }
    };

    //! Basis function on 3D: u(x,y,z) = x*y^2
    struct Basis_A13_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[1]*p[1]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return p[1]*p[1]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 2.*p[0]*p[1]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = z*y^2
    struct Basis_A14_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[2]*p[1]*p[1]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 2.*p[1]*p[2]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return p[1]*p[1]; }
    };

    //! Basis function on 3D: u(x,y,z) = x*z^2
    struct Basis_A15_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[2]*p[2]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return p[2]*p[2]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 2.*p[0]*p[2]; }
    };

    //! Basis function on 3D: u(x,y,z) = y*z^2
    struct Basis_A16_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[1]*p[2]*p[2]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return p[2]*p[2]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 2.*p[1]*p[2]; }
    };

    //! Basis function on 3D: u(x,y,z) = x^3
    struct Basis_A17_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[0]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 3.*p[0]*p[0]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = y^3
    struct Basis_A18_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[1]*p[1]*p[1]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 3.*p[1]*p[1]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = z^3
    struct Basis_A19_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[2]*p[2]*p[2]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 3.*p[2]*p[2]; }
    };

    //! Basis function on 3D: u(x,y,z) = x^2*y*z
    struct Basis_A20_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[1]*p[2]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 2.*p[0]*p[1]*p[2]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return p[0]*p[0]*p[2]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return p[0]*p[0]*p[1]; }
    };

    //! Basis function on 3D: u(x,y,z) = x*y^2*z
    struct Basis_A21_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[1]*p[1]*p[2]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return p[1]*p[1]*p[2]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 2.*p[0]*p[1]*p[2]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return p[0]*p[1]*p[1]; }
    };

    //! Basis function on 3D: u(x,y,z) = x*y*z^2
    struct Basis_A22_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[1]*p[2]*p[2]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return p[1]*p[2]*p[2]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return p[0]*p[2]*p[2]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 2.*p[0]*p[1]*p[2]; }
    };

    //! Basis function on 3D: u(x,y,z) = x^3*y
    struct Basis_A23_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[0]*p[1]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 3.*p[0]*p[0]*p[1]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return p[0]*p[0]*p[0]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = x^3*z
    struct Basis_A24_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[0]*p[2]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 3.*p[0]*p[0]*p[2]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return p[0]*p[0]*p[0]; }
    };

    //! Basis function on 3D: u(x,y,z) = y^3*x
    struct Basis_A25_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[1]*p[1]*p[1]*p[0]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return p[1]*p[1]*p[1]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 3.*p[1]*p[1]*p[0]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = y^3*z
    struct Basis_A26_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[1]*p[1]*p[1]*p[2]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 3.*p[1]*p[1]*p[2]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return p[1]*p[1]*p[1]; }
    };

    //! Basis function on 3D: u(x,y,z) = z^3*x
    struct Basis_A27_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[2]*p[2]*p[2]*p[0]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return p[2]*p[2]*p[2]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 3.*p[2]*p[2]*p[0]; }
    };

    //! Basis function on 3D: u(x,y,z) = z^3*y
    struct Basis_A28_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[2]*p[2]*p[2]*p[1]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return p[2]*p[2]*p[2]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 3.*p[2]*p[2]*p[1]; }
    };

    //! Basis function on 3D: u(x,y,z) = x^2*y^2
    struct Basis_A29_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[1]*p[1]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 2.*p[0]*p[1]*p[1]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 2.*p[0]*p[0]*p[1]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = y^2*z^2
    struct Basis_A30_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[1]*p[1]*p[2]*p[2]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 2.*p[1]*p[2]*p[2]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 2.*p[1]*p[1]*p[2]; }
    };

    //! Basis function on 3D: u(x,y,z) = z^2*x^2
    struct Basis_A31_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[2]*p[2]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 2.*p[0]*p[2]*p[2]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 2.*p[0]*p[0]*p[2]; }
    };

    //! Basis function on 3D: u(x,y,z) = x^4
    struct Basis_A32_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[0]*p[0]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 4.*p[0]*p[0]*p[0]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = y^4
    struct Basis_A33_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[1]*p[1]*p[1]*p[1]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 4.*p[1]*p[1]*p[1]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = z^4
    struct Basis_A34_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[2]*p[2]*p[2]*p[2]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 4.*p[2]*p[2]*p[2]; }
    };

    //! Basis function on 3D: u(x,y,z) = x^3*y*z
    struct Basis_A35_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[0]*p[1]*p[2]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 3.*p[0]*p[0]*p[1]*p[2]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return p[0]*p[0]*p[0]*p[2]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return p[0]*p[0]*p[0]*p[1]; }
    };

    //! Basis function on 3D: u(x,y,z) = x*y^3*z
    struct Basis_A36_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[1]*p[1]*p[1]*p[2]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return p[1]*p[1]*p[1]*p[2]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 3.*p[0]*p[1]*p[1]*p[2]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return p[0]*p[1]*p[1]*p[1]; }
    };

    //! Basis function on 3D: u(x,y,z) = x*y*z^3
    struct Basis_A37_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[1]*p[2]*p[2]*p[2]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return p[1]*p[2]*p[2]*p[2]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return p[0]*p[2]*p[2]*p[2]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 3.*p[0]*p[1]*p[2]*p[2]; }
    };

    //! Basis function on 3D: u(x,y,z) = x^2*y^2*z
    struct Basis_A38_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[1]*p[1]*p[2]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 2.*p[0]*p[1]*p[1]*p[2]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 2.*p[0]*p[0]*p[1]*p[2]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return p[0]*p[0]*p[1]*p[1]; }
    };

    //! Basis function on 3D: u(x,y,z) = x^2*y*z^2
    struct Basis_A39_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[1]*p[2]*p[2]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 2.*p[0]*p[1]*p[2]*p[2]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return p[0]*p[0]*p[2]*p[2]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 2.*p[0]*p[0]*p[1]*p[2]; }
    };

    //! Basis function on 3D: u(x,y,z) = x*y^2*z^2
    struct Basis_A40_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[1]*p[1]*p[2]*p[2]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return p[1]*p[1]*p[2]*p[2]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 2.*p[0]*p[1]*p[2]*p[2]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 2.*p[0]*p[1]*p[1]*p[2]; }
    };

    //! Basis function on 3D: u(x,y,z) = x^4*y
    struct Basis_A41_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[0]*p[0]*p[1]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 4.*p[0]*p[0]*p[0]*p[1]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return p[0]*p[0]*p[0]*p[0]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = x^4*z
    struct Basis_A42_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[0]*p[0]*p[2]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 4.*p[0]*p[0]*p[0]*p[2]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return p[0]*p[0]*p[0]*p[0]; }
    };

    //! Basis function on 3D: u(x,y,z) = y^4*x
    struct Basis_A43_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[1]*p[1]*p[1]*p[1]*p[0]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return p[1]*p[1]*p[1]*p[1]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 4.*p[1]*p[1]*p[1]*p[0]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = y^4*z
    struct Basis_A44_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[1]*p[1]*p[1]*p[1]*p[2]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 4.*p[1]*p[1]*p[1]*p[2]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return p[1]*p[1]*p[1]*p[1]; }
    };

    //! Basis function on 3D: u(x,y,z) = z^4*x
    struct Basis_A45_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[2]*p[2]*p[2]*p[2]*p[0]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return p[2]*p[2]*p[2]*p[2]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 4.*p[2]*p[2]*p[2]*p[0]; }
    };

    //! Basis function on 3D: u(x,y,z) = z^4*y
    struct Basis_A46_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[2]*p[2]*p[2]*p[2]*p[1]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return p[2]*p[2]*p[2]*p[2]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 4.*p[2]*p[2]*p[2]*p[1]; }
    };

    //! Basis function on 3D: u(x,y,z) = x^3*y^2
    struct Basis_A47_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[0]*p[1]*p[1]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 3.*p[0]*p[0]*p[1]*p[1]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 2.*p[0]*p[0]*p[0]*p[1]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = x^3*z^2
    struct Basis_A48_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[0]*p[2]*p[2]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 3.*p[0]*p[0]*p[2]*p[2]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 2.*p[0]*p[0]*p[0]*p[2]; }
    };

    //! Basis function on 3D: u(x,y,z) = y^3*x^2
    struct Basis_A49_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[1]*p[1]*p[1]*p[0]*p[0]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 2.*p[1]*p[1]*p[1]*p[0]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 3.*p[1]*p[1]*p[0]*p[0]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = y^3*z^2
    struct Basis_A50_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[1]*p[1]*p[1]*p[2]*p[2]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 3.*p[1]*p[1]*p[2]*p[2]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 2.*p[1]*p[1]*p[1]*p[2]; }
    };

    //! Basis function on 3D: u(x,y,z) = z^3*x^2
    struct Basis_A51_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[2]*p[2]*p[2]*p[0]*p[0]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 2.*p[2]*p[2]*p[2]*p[0]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 3.*p[2]*p[2]*p[0]*p[0]; }
    };

    //! Basis function on 3D: u(x,y,z) = z^3*y^2
    struct Basis_A52_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[2]*p[2]*p[2]*p[1]*p[1]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 2.*p[2]*p[2]*p[2]*p[1]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 3.*p[2]*p[2]*p[1]*p[1]; }
    };

    //! Basis function on 3D: u(x,y,z) = x^5
    struct Basis_A53_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[0]*p[0]*p[0]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 5.*p[0]*p[0]*p[0]*p[0]; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = y^5
    struct Basis_A54_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[1]*p[1]*p[1]*p[1]*p[1]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 5.*p[1]*p[1]*p[1]*p[1]; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = z^5
    struct Basis_A55_3D: public Base::BaseBasisFunction
    {
        virtual double eval(const PointReferenceT& p) const         { return p[2]*p[2]*p[2]*p[2]*p[2]; }
        virtual double evalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double evalDeriv2(const PointReferenceT& p) const   { return 5.*p[2]*p[2]*p[2]*p[2]; }
    };
    
};
#endif  // BasisFunctionsCollection_A_hpp
