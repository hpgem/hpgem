//
//  PhysGradientOfBasisFunction.h
//  
//
//  Created by Shavarsh Nurijanyan on 7/25/13.
//
//

#ifndef ____PhysGradientOfBasisFunction__
#define ____PhysGradientOfBasisFunction__

#include <iostream>
#include "Geometry/PointReference.hpp"
#include "Base/Element.hpp"


namespace Utilities
{
    /*! For a basis function we also need its physical space gradient as opposed
     *  to the one in reference space (which can be harvested by evaluating the
     *  derivatives of basis functions). Hence this class computes the
     *  reference space gradient, transforms it with the Jacobian of the mapping
     *  and thus yields the physical space gradient. */
    template <unsigned int DIM, class EType>
    struct PhysGradientOfBasisFunction
    {
        typedef Geometry::PointReference<DIM>             PointReferenceT;
        
            //typedef double (*DerivativeOfBasisFuncPtr)(const PointReferenceT&);
        typedef typename LinearAlgebra::NumericalVector RetType;
        
        PhysGradientOfBasisFunction(const EType* e, unsigned int bFuncNr):
            myElement_(e),
            myFunctionNumber_(bFuncNr)
	    {
            
        }
        
            //! Evaluation operator, also compatible with integration routines.
        void operator()(const PointReferenceT* p, RetType& r) const
	    {
                ////SHOULD BE OPENED AFTER MATRIX FIX!!!!
            
            
                // Calculate the gradient with respect to ref space coords and
                // store it in B:
//            LinearAlgebra::Matrix  B(DIM, 1);
//            for (unsigned int d = 0; d < DIM; ++d)
//                B(d, 0) = myDeriveFunctionPtr_[myFunctionNumber_*DIM+d](P);
//                
//                    // get the jacobian at P:
//            typename ElementGeometry<DIM>::JacobianT jac;
//            myElement_->calcJacobian(P, jac);
//            
//            std::vector<int> ipiv (DIM);  // pivot vector
//            lapack::getrf(
//                          static_cast<ublas::matrix<double, ublas::column_major>&>
//                          (jac), ipiv);		// factorize
//            lapack::getrs(
//                          'T',
//                          static_cast<ublas::matrix<double, ublas::column_major>&>
//                          (jac), ipiv, B);	// solve for jac^T
//            
//            for(unsigned int i = 0; i < DIM; ++i)
//            {
//                r[i] = B(i, 0);
//            }
	    }
        
    private:
        const EType*                        myElement_;
        const unsigned int                  myFunctionNumber_;
    };
    
        // Fixed for 1D case: M.T.Julianto Feb 14, 2010
    template <>
    struct PhysGradientOfBasisFunction<1, Base::Element<1> >
    {
        
        typedef Geometry::PointReference<1>             PointReferenceT;
            // typedef double (*DerivativeOfBasisFuncPtr)(const Geometry::PointReference<1>&);
         typedef typename LinearAlgebra::NumericalVector RetType;
        
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
            typename Geometry::ElementGeometry<1>::JacobianT jac;
            
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
        typedef typename LinearAlgebra::NumericalVector RetType;
        
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
            typename  Geometry::ElementGeometry<2>::JacobianT jac;
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
          typedef typename LinearAlgebra::NumericalVector RetType;
        
        
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
            typename  Geometry::ElementGeometry<3>::JacobianT jac;
            
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
    };
    
} // namespace


#endif /* defined(____PhysGradientOfBasisFunction__) */
