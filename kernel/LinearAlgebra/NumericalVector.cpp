/*! \file NumericalVector.cpp
\brief This has the Numerical Vector for hpGEM
Details.
*/

#include "NumericalVector.hpp"

namespace LinearAlgebra
{

    extern "C"

    {
    
    
    ///This is the gernal scale times vector + vector from blas, hence from blas level 1. Here we also use on a matrix by treating as a vector
    int daxpy_(unsigned int* N, double* DA, double* DX,unsigned int* INCX, double* DY, unsigned int* INCY);
    
    
    }
    
    
    NumericalVector::NumericalVector() : data_() {}
    
    NumericalVector::NumericalVector(int m) : data_(m){}
    
    NumericalVector::NumericalVector(const NumericalVector& other) : data_(other.data_){}
    
    NumericalVector::NumericalVector(const double array[], int size) : data_(array,size){}
    
    void NumericalVector::resize(unsigned int size) { data_.resize(size); }
    
    NumericalVector& NumericalVector::operator= (const NumericalVector& right){data_=right.data_; return *this;}
    
    NumericalVector NumericalVector::operator+ (const NumericalVector& right)
    {
        NumericalVector result(*this);
        result.data_+=right.data_;
        #ifndef IAMNIC0
        //  std::cout << "Nico you are using a slow operator be warned of the Numerical Vector, please try and rewrite to use += not +" << std::endl;
        #endif
        return result;
    }
    
   
    NumericalVector NumericalVector::operator+ (const NumericalVector& right) const
    {
        NumericalVector result(*this);
        result.data_+=right.data_;
#ifndef IAMNICO
        //   std::cout << "Nico you are using a slow operator be warned of the Numerical Vector, please try and rewrite to use += not +" << std::endl;
#endif
        return result;
    }
    
    NumericalVector NumericalVector::operator- (const NumericalVector& right)
    {
        NumericalVector result(*this);
        result.data_-=right.data_;
#ifndef IAMNICO
        //   std::cout << "Nico you are using a slow operator be warned of the Numerical Vector, please try and rewrite to use -= not -" << std::endl;
#endif
        return result;
    }
    
    NumericalVector NumericalVector::operator- (const NumericalVector& right) const
    {
        NumericalVector result(*this);
        result.data_-=right.data_;
#ifndef IAMNICO
        //  std::cout << "Nico you are using a slow operator be warned of the Numerical Vector, please try and rewrite to use -= not -" << std::endl;
#endif
        return result;
    }
    
    NumericalVector NumericalVector::operator* (const double& right)
    {
        NumericalVector result(*this);
        result.data_*=right;
#ifndef IAMNICO
        //   std::cout << "Nico You are using a slow operator be warned of the Numerical Vector, please try and rewrite to use *= not *" << std::endl;
#endif
        return result;
    }
    
    NumericalVector NumericalVector::operator* (const double& right) const
    {
        NumericalVector result(*this);
        result.data_*=right;
#ifndef IAMNICO
        //  std::cout << "Nico you are using a slow operator be warned of the Numerical Vector, please try and rewrite to use *= not *" << std::endl;
#endif
        return result;
    }
    
    double NumericalVector::operator* (const NumericalVector& right) const
    {
        ///\TODO replace with BLAS (I dont know where to find them)
        if(this->size()!=right.size())
            throw "vector \\cdot vector product only defined for vectors of the same sizes";
        double result(0);
        for(int i=0;i<right.size();++i){
            result+=data_[i]*right[i];
        }
        return result;
    }
    
    NumericalVector& NumericalVector::operator/= (const double& right)
    {
        data_/=right;
        return *this;
    }
    
    NumericalVector NumericalVector::operator/ (const double& right)
    {
        NumericalVector result(*this);
        return (result/=right);
        
    }
    
    void NumericalVector::axpy(double a, const NumericalVector& x)
    {
        
        unsigned int size=data_.size();
        
        unsigned int i_one=1;
        
        
        daxpy_(&size, &a, &((*(const_cast<NumericalVector *> (&x)))[0]), &i_one, &((*this)[0]) , &i_one);
        
        
    }
    
    bool NumericalVector::operator== (const NumericalVector& right)
    {
        for (int i = 0; i < data_.size(); ++i)
        {
            if (data_[i] != right.data_[i]) return false;
        }
        return true;
    }
    
    bool NumericalVector::operator== (const NumericalVector& right) const
    {
        for (int i = 0; i < data_.size(); ++i)
        {
            if (data_[i] != right.data_[i]) return false;
        }
        return true;
    }
    
    bool NumericalVector::operator< (const NumericalVector& right) const
    {
        for(int i=0;i<data_.size();++i)
        {
            if(data_[i]<right.data_[i]) return true;
            if(data_[i]>right.data_[i]) return false;
        }
	    return false;
    }
    
    NumericalVector& NumericalVector::operator+= (const NumericalVector& right){data_+=right.data_; return *this;}
    
    NumericalVector& NumericalVector::operator-= (const NumericalVector& right){data_-=right.data_; return *this;}
    
    NumericalVector& NumericalVector::operator*= (const double& right){data_*=right; return *this;}
    
    inline double& NumericalVector::operator[] (const unsigned int n) {return data_[n];}
    
    NumericalVector operator*(const double& left, const NumericalVector& right)
    {
        NumericalVector result(right);
        result.data_*=left;
        return result;
    }
    
    NumericalVector   operator-(const NumericalVector& right){return NumericalVector(right * -1.0);}
    
    
    ostream& operator<<(ostream& os, const NumericalVector& A)
    {
        os<< '[';
        for (int i=0;i<A.data_.size()-1;i++){os<< A(i)<<',';}
        os<< A(A.data_.size()-1)<<']';
        return os;
    }
    
    

    
    
    
}
