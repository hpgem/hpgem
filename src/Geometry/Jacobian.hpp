/*
 * Jacobian.hpp
 *
 *  Created on: Feb 5, 2013
 *      Author: nicorivas
 */

#ifndef JACOBIAN_HPP_
#define JACOBIAN_HPP_

#include "../LinearAlgebra/Matrix.hpp"
#include "../LinearAlgebra/NumericalVector.hpp"
#include "PointPhysical.hpp"




namespace Geometry
{
    class Jacobian: public LinearAlgebra::Matrix
    {
        public:
        typedef Jacobian            JacobianT;
        typedef PointPhysical                PhysicalPointT;

        public:
            // Constructors.
        Jacobian(unsigned int dimTo,unsigned int dimFrom);
        Jacobian(const JacobianT& jacobian);
       
        double determinant()const;
        void   computeWedgeStuffVector(NumericalVector& p)const;
        
        virtual ~Jacobian(){}
        
//        void operator *=(const Jacobian& j2)
//            {
//                    //((LinearAlgebra::Matrix)(*this)) *= ((LinearAlgebra::Matrix)(j2));
//            }
     
        /*! (OC): ConcatenatedMapping has to be able to do a matrix product on the
         Jacobians of two (successively applied) mappings. Therefore we provide
         the function multiplyJacobiansInto. */
        
        void multiplyJacobiansInto(const Jacobian& jac2,
                                   Jacobian& jres)
        {
                // TODO: This is very inefficient, because of Anthony's code in LinearAlgebra.
                LinearAlgebra::Matrix& matThis= *this;
                LinearAlgebra::Matrix& matRes= jres;
                                             
                const LinearAlgebra::Matrix& mat2=jac2;
            
                matRes = matThis.operator*(mat2);
            
        }
        
    };
    
        /// \bug create as a part of Jacobian class
 
}
#endif /* JACOBIAN_HPP_ */
