

namespace Geometry
{
    template <unsigned int dimFrom, unsigned int dimTo>
    class Jacobian;
    
    template <unsigned int dimFrom, unsigned int dimTo>
    Jacobian<dimFrom, dimTo>::Jacobian():
        LinearAlgebra::Matrix(dimTo,dimFrom)
    {
    }
    template <unsigned int dimFrom, unsigned int dimTo>
    Jacobian<dimFrom, dimTo>::Jacobian(const JacobianT& jacobian):
        LinearAlgebra::Matrix(jacobian)
    {
    }

    template <unsigned int dimFrom, unsigned int dimTo>
    void
    Jacobian<dimFrom, dimTo>::computeWedgeStuffVector(PhysicalPointT& p)const
    {
        NumericalVector&      v((NumericalVector&)p);
        
        const LinearAlgebra::Matrix& jac=*this;

            //cout << "jacobian="<<jac<<endl;

        jac.computeWedgeStuffVector(v);
    }

        /// The computation of Jacobians are harcoded up until 4D, to make it faster.
    template <unsigned int dimFrom, unsigned int dimTo>
    double
    Jacobian<dimFrom, dimTo>::determinant()const
    {
        if (dimFrom!=dimTo)
        {
            cout<<"Jacobian should be square to have a determinant!"<<endl;
        }
        else
        {
            switch (dimFrom)
            {
                case 1:
                    return (*this)(0,0);
                break;
                case 2:
                    return (*this)(0,0) * (*this)(1,1) - (*this)(0,1) * (*this)(1,0);
                break;

                case 3:
                    return (*this)(0,0) * ((*this)(1,1) * (*this)(2,2) - (*this)(1,2) * (*this)(2,1))
                                      -(*this)(0,1) * ((*this)(1,0) * (*this)(2,2) - (*this)(2,0) * (*this)(1,2))
                                       +(*this)(0,2) * ((*this)(1,0) * (*this)(2,1) - (*this)(2,0) * (*this)(1,1));
                break;

                case 4:
                    return ((*this)(3,0)*(*this)(2,1)*(*this)(0,3)-(*this)(2,0)*(*this)(3,1)*(*this)(0,3))*(*this)(1,2)
                                         +(-(*this)(3,0)*(*this)(0,3)*(*this)(2,2)+(*this)(2,0)*(*this)(0,3)*(*this)(3,2))*(*this)(1,1)
                                         +((*this)(3,1)*(*this)(0,3)*(*this)(2,2)-(*this)(2,1)*(*this)(0,3)*(*this)(3,2))*(*this)(1,0)
                                         +(-(*this)(3,0)*(*this)(2,1)*(*this)(1,3)+(*this)(2,0)*(*this)(3,1)*(*this)(1,3)
                                                 +(-(*this)(2,0)*(*this)(3,3)+(*this)(3,0)*(*this)(2,3))*(*this)(1,1)
                                                 +((*this)(2,1)*(*this)(3,3)-(*this)(3,1)*(*this)(2,3))*(*this)(1,0))*(*this)(0,2)
                                                 +((*this)(3,0)*(*this)(1,3)*(*this)(2,2)-(*this)(2,0)*(*this)(1,3)*(*this)(3,2)
                                                         +((*this)(2,0)*(*this)(3,3)-(*this)(3,0)*(*this)(2,3))*(*this)(1,2)
                                                         +(-(*this)(2,2)*(*this)(3,3)+(*this)(2,3)*(*this)(3,2))*(*this)(1,0))*(*this)(0,1)
                                                         +(-(*this)(3,1)*(*this)(1,3)*(*this)(2,2)+(*this)(2,1)*(*this)(1,3)*(*this)(3,2)
                                                                 +((*this)(3,1)*(*this)(2,3)-(*this)(2,1)*(*this)(3,3))*(*this)(1,2)
                                                                 +(*this)(1,1)*((*this)(2,2)*(*this)(3,3)-(*this)(2,3)*(*this)(3,2)))*(*this)(0,0);
                             // ... says Maple; this can possibly be done more efficiently,
                             // maybe even with LU (with pivoting, though...)
                    break;
                default:
                    cout<<"This dimension is not implemented"<<endl;
                    return -1;
            }
        }
    }
};
