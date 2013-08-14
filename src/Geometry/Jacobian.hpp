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
    template <unsigned int dimFrom, unsigned int dimTo>
    class Jacobian: public LinearAlgebra::Matrix
    {
        public:
        typedef Jacobian<dimFrom, dimTo>            JacobianT;
        typedef PointPhysical<dimTo>                PhysicalPointT;

        public:
            // Constructors.
        Jacobian();
        Jacobian(const JacobianT& jacobian);
       
        double determinant()const;
        void   computeWedgeStuffVector(PhysicalPointT& p)const;
        
        virtual ~Jacobian(){}
        
//        void operator *=(const Jacobian& j2)
//            {
//                    //((LinearAlgebra::Matrix)(*this)) *= ((LinearAlgebra::Matrix)(j2));
//            }
     
        /*! (OC): ConcatenatedMapping has to be able to do a matrix product on the
         Jacobians of two (successively applied) mappings. Therefore we provide
         the function multiplyJacobiansInto. */
        
        template <unsigned int dim3>
        void multiplyJacobiansInto(const Jacobian<dim3, dimFrom>& jac2,
                                   Jacobian<dim3, dimTo>& jres)
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
#include "Jacobian_Impl.hpp"
#endif /* JACOBIAN_HPP_ */
