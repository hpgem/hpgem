/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, University of Twente
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef BasisFunctionsCollection_B_hpp
#define BasisFunctionsCollection_B_hpp

#include "BaseBasisFunction.hpp"
#include <vector>

namespace Base
{

    //---------------------------------------------------------------------------
    //------------------------------------------------ basis collection for dim=1
    
    //! Basis function on 1D: u(x) = 1
    class Basis_B0_1D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return 1; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 1D: u(x) = x
    class Basis_B1_1D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return p[0]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 1.; }
    };

    //! Basis function on 1D: u(x) = x^2
    class Basis_B2_1D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return p[0]*p[0]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 2.*p[0]; }
    };

    //! Basis function on 1D: u(x) = x^3
    class Basis_B3_1D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[0]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 3.*p[0]*p[0]; }
    };

    //---------------------------------------------------------------------------
    //------------------------------------------------ basis collection for dim=2
    
    //! Basis function on 2D: u(x,y) = 1
    class Basis_B0_2D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return 1; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 2D: u(x,y) = x
    class Basis_B1_2D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return p[0]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 1.; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 2D: u(x,y) = y
    class Basis_B2_2D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return p[1]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 1.; }
    };

    //! Basis function on 2D: u(x,y) = x*y
    class Basis_B3_2D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return p[0]*p[1]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return p[1]; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return p[0]; }
    };

    //! Basis function on 2D: u(x,y) = x^2
    class Basis_B4_2D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return p[0]*p[0]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 2.*p[0]; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 2D: u(x,y) = y^2
    class Basis_B5_2D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return p[1]*p[1]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 2.*p[1]; }
    };

    //! Basis function on 2D: u(x,y) = y*x^2
    class Basis_B6_2D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[1]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 2.*p[0]*p[1]; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return p[0]*p[0]; }
    };

    //! Basis function on 2D: u(x,y) = x*y^2
    class Basis_B7_2D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return p[0]*p[1]*p[1]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return p[1]*p[1]; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 2.*p[0]*p[1]; }
    };

    //! Basis function on 2D: u(x,y) = x^3
    class Basis_B8_2D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[0]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 3.*p[0]*p[0]; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 2D: u(x,y) = y^3
    class Basis_B9_2D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return p[1]*p[1]*p[1]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 3.*p[1]*p[1]; }
    };


    //---------------------------------------------------------------------------
    //------------------------------------------------ basis collection for dim=3
    
    //! Basis function on 3D: u(x,y,z) = 1
    class Basis_B0_3D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return 1; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = x
    class Basis_B1_3D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return p[0]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 1.; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = y
    class Basis_B2_3D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return p[1]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 1.; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = z
    class Basis_B3_3D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return p[2]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return 1.; }
    };

    //! Basis function on 3D: u(x,y,z) = x*y
    class Basis_B4_3D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return p[0]*p[1]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return p[1]; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return p[0]; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = x*z
    class Basis_B5_3D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return p[0]*p[2]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return p[2]; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return p[0]; }
    };

    //! Basis function on 3D: u(x,y,z) = y*z
    class Basis_B6_3D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return p[1]*p[2]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return p[2]; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return p[1]; }
    };

    //! Basis function on 3D: u(x,y,z) = x^2
    class Basis_B7_3D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return p[0]*p[0]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 2.*p[0]; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = y^2
    class Basis_B8_3D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return p[1]*p[1]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 2.*p[1]; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = z^2
    class Basis_B9_3D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return p[2]*p[2]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return 2.*p[2]; }
    };

    //! Basis function on 3D: u(x,y,z) = x*y*z
    class Basis_B10_3D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return p[0]*p[1]*p[2]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return p[1]*p[2]; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return p[0]*p[2]; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return p[0]*p[1]; }
    };

    //! Basis function on 3D: u(x,y,z) = y*x^2
    class Basis_B11_3D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[1]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 2.*p[0]*p[1]; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return p[0]*p[0]; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = z*x^2
    class Basis_B12_3D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[2]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 2.*p[0]*p[2]; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return p[0]*p[0]; }
    };

    //! Basis function on 3D: u(x,y,z) = x*y^2
    class Basis_B13_3D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return p[0]*p[1]*p[1]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return p[1]*p[1]; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 2.*p[0]*p[1]; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = z*y^2
    class Basis_B14_3D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return p[2]*p[1]*p[1]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 2.*p[1]*p[2]; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return p[1]*p[1]; }
    };

    //! Basis function on 3D: u(x,y,z) = x*z^2
    class Basis_B15_3D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return p[0]*p[2]*p[2]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return p[2]*p[2]; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return 2.*p[0]*p[2]; }
    };

    //! Basis function on 3D: u(x,y,z) = y*z^2
    class Basis_B16_3D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return p[1]*p[2]*p[2]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return p[2]*p[2]; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return 2.*p[1]*p[2]; }
    };

    //! Basis function on 3D: u(x,y,z) = x^3
    class Basis_B17_3D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return p[0]*p[0]*p[0]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 3.*p[0]*p[0]; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = y^3
    class Basis_B18_3D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return p[1]*p[1]*p[1]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 3.*p[1]*p[1]; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return 0.; }
    };

    //! Basis function on 3D: u(x,y,z) = z^3
    class Basis_B19_3D : public Base::BaseBasisFunction
    {
    public:
        virtual double Eval(const PointReferenceT& p) const         { return p[2]*p[2]*p[2]; }
        virtual double EvalDeriv0(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv1(const PointReferenceT& p) const   { return 0.; }
        virtual double EvalDeriv2(const PointReferenceT& p) const   { return 3.*p[2]*p[2]; }
    };
};
#endif  // BasisFunctionsCollection_A_hpp
