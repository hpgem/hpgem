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
//

#include "PhysGradientOfBasisFunction.hpp"
#include "Geometry/PointReference.hpp"
#include "Geometry/Jacobian.hpp"
#include "Element.hpp"
#include "Geometry/ReferenceGeometry.hpp"
#include "ElementCacheData.hpp"

namespace Utilities {

	// Fixed for 1D case: M.T.Julianto Feb 14, 2010
	/*   template <>
	 struct PhysGradientOfBasisFunction<1, Base::Element<1> >
	 {
	 
	 typedef Geometry::PointReference<1>             PointReferenceT;
	 // typedef double (*DerivativeOfBasisFuncPtr)(const Geometry::PointReference<1>&);
	 typedef LinearAlgebra::NumericalVector RetType;
	 
	 PhysGradientOfBasisFunction(const Base::Element<1>* e, unsigned int bFuncNr):
	 myElement_(e),
	 myFunctionNumber_(bFuncNr)
	 {
	 
	 }
	 
	 //! Evaluation operator, also compatible with integration routines.
	 void operator()(const PointReferenceT& p, RetType& r) const
	 {
	 // Calculate the gradient with respect to ref space coords and
	 // store it in B:
	 double B = myElement_->basisFunctionDeriv(myFunctionNumber_,0, p);
	 //myDeriveFunctionPtr_[myFunctionNumber_*1](P);
	 
	 // get the jacobian at P:
	 Geometry::ElementGeometry<1>::JacobianT jac;
	 
	 myElement_->calcJacobian(p, jac);
	 
	 r[0]=1.0/jac(0,0) * B;
	 }
	 
	 private:
	 const Base::Element<1>*                   myElement_;
	 // const DerivativeOfBasisFuncPtr*     myDeriveFunctionPtr_; // a pointer to the derivatives
	 const unsigned int                  myFunctionNumber_;

	 };
	 
	 
	 template <>
	 struct PhysGradientOfBasisFunction<2, Base::Element<2> >
	 {
	 typedef Geometry::PointReference<2>             PointReferenceT;
	 // typedef double (*DerivativeOfBasisFuncPtr)(const Geometry::PointReference<2>&);
	 typedef LinearAlgebra::NumericalVector RetType;
	 
	 PhysGradientOfBasisFunction(const Base::Element<2>* e, unsigned int bFuncNr):
	 myElement_(e),
	 myFunctionNumber_(bFuncNr)
	 {
	 
	 }
	 
	 //! Evaluation operator, also compatible with integration routines.
	 void operator()(const PointReferenceT& p, RetType& r) const
	 {
	 
	 // Calculate the gradient with respect to ref space coords and
	 // store it in B:
	 LinearAlgebra::Matrix B(2, 1);
	 for (unsigned int d = 0; d < 2; ++d)
	 B(d, 0) = myElement_->basisFunctionDeriv(myFunctionNumber_,d, p);
	 //myDeriveFunctionPtr_[myFunctionNumber_*2+d](P);
	 
	 // get the jacobian at P:
	 Geometry::ElementGeometry<2>::JacobianT jac;
	 myElement_->calcJacobian(p, jac);
	 
	 //LinearAlgebra::Matrix A_T_Inv;
	 
	 double a=jac(0,0);
	 double b=jac(1,0);
	 double c=jac(0,1);
	 double d=jac(1,1);
	 
	 double det=1.0/(a*d -b*c);
	 
	 
	 r[0]=det * (d*B(0,0) - b*B(1,0));
	 
	 r[1]=det * (a*B(1,0) - c*B(0,0));
	 
	 }
	 
	 private:
	 const Base::Element<2>*                   myElement_;
	 // const DerivativeOfBasisFuncPtr*     myDeriveFunctionPtr_; // a pointer to the derivatives
	 const unsigned int                  myFunctionNumber_;
	 };
	 
	 template <>
	 struct PhysGradientOfBasisFunction<3, Base::Element<3> >
	 {
	 typedef Geometry::PointReference<3>             PointReferenceT;
	 //  typedef double (*DerivativeOfBasisFuncPtr)(const Geometry::PointReference<3>&);
	 typedef LinearAlgebra::NumericalVector RetType;
	 
	 
	 PhysGradientOfBasisFunction(const Base::Element<3>* e, unsigned int bFuncNr):
	 myElement_(e),
	 myFunctionNumber_(bFuncNr)
	 {
	 
	 }
	 
	 //! Evaluation operator, also compatible with integration routines.
	 void operator()(const PointReferenceT& p, RetType& r) const
	 {
	 
	 // Calculate the gradient with respect to ref space coords and
	 // store it in B:
	 LinearAlgebra::Matrix B(3, 1);
	 for (unsigned int d = 0; d < 3; ++d)
	 B(d, 0) =myElement_->basisFunctionDeriv(myFunctionNumber_,d, p);
	 //myDeriveFunctionPtr_[myFunctionNumber_*3+d](P);
	 
	 // get the jacobian at P:
	 Geometry::ElementGeometry<3>::JacobianT jac;
	 
	 myElement_->calcJacobian(p, jac);
	 
	 
	 double a=jac(0,0);
	 double b=jac(1,0);
	 double c=jac(2,0);
	 double d=jac(0,1);
	 double e=jac(1,1);
	 double f=jac(2,1);
	 double g=jac(0,2);
	 double h=jac(1,2);
	 double i=jac(2,2);
	 
	 double t1=e*i - f*h;
	 double t2=f*g - d*i;
	 double t3=d*h - e*g;
	 
	 double det=1.0/(a*t1 + b*t2 + c*t3);
	 
	 r[0]=det * (t1 * B(0,0) 
	 + (c*h - b*i) * B(1,0) // + (b*i - c*h) * B(1,0) 
	 + (b*f - c*e) * B(2,0));
	 
	 r[1]=det * (t2 * B(0,0) 
	 + (a*i - c*g) * B(1,0) 
	 + (c*d - a*f) * B(2,0));
	 
	 r[2]=det * (t3 * B(0,0) 
	 + (b*g - a*h) * B(1,0) 
	 + (a*e - b*d) * B(2,0));
	 
	 }
	 
	 private:
	 const Base::Element<3>*                    myElement_;
	 //const DerivativeOfBasisFuncPtr*     myDeriveFunctionPtr_; // a pointer to the derivatives
	 const unsigned int                  myFunctionNumber_;
	 };*/
	void PhysGradientOfBasisFunction::operator ()(const PointReferenceT& p, RetType& r) const {
		const unsigned int DIM = p.size();
		r.resize(DIM);
                r*=0;
		static RetType dummy(DIM);
                dummy.resize(DIM);
		static Geometry::Jacobian jac(DIM, DIM);
                jac.resize(DIM,DIM);
		myElement_->calcJacobian(p, jac);
		myElement_->getReferenceGeometry()->getBasisFunctionDerivative(myFunction_, p, dummy);
		jac = jac.inverse();
		//r*=jac;///\todo can someone who knows BLAS update the linAlg routines?
		for (int i = 0; i < DIM; ++i) {
                        //std::cout<<dummy[i]<<" ";
			for (int j = 0; j < DIM; ++j) {
				r[i] += dummy[j] * jac(j, i);
			}
		}
                //std::cout<<std::endl;
	}

}
