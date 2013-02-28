#ifndef NumericalVectorHPP
#define NumericalVectorHPP

#include "GlobalNamespaceLinearAlgebra.hpp"
#include <valarray>
#include <iostream>

#define IAMNICO

namespace LinearAlgebra
{
    extern "C"
    
    {
     
        
        ///This is the gernal scale times vector + vector from blas, hence from blas level 1. Here we also use on a matrix by treating as a vector
        int daxpy_(unsigned int* N, double* DA, double* DX,unsigned int* INCX, double* DY, unsigned int* INCY);
        
        
    }
    
    //This is dervied from valarray so import that information
    using std::valarray;
    using std::ostream;
    
    /// \class NumericalVector
    /// \brief This is a vector of doubles
    ///
    /// \details
    /// This implements a vector of doubles and all the standard operators for it.
    /// Note it is encapulating a valarray for its data storage.
    class NumericalVector
    {
        
    public:
        
        NumericalVector() : data_() {}
        
        NumericalVector(int m) : data_(m){}
        
        NumericalVector(const NumericalVector& other) : data_(other.data_){}
     
        
        //NumericalVector(double* array, int size) : data_(array,size){}
        
        NumericalVector(double array[], int size) : data_(array,size){}
        
        void resize(unsigned int size) { data_.resize(size); }

        NumericalVector& operator= (const NumericalVector& right){data_=right.data_; return *this;}
        
        NumericalVector operator+ (const NumericalVector& right)
        {
            NumericalVector result(*this);
            result.data_+=right.data_; 
            #ifndef IAMNIC0
              //  std::cout << "Nico you are using a slow operator be warned of the Numerical Vector, please try and rewrite to use += not +" << std::endl;
            #endif
            return result;
        }
        
        NumericalVector operator+ (const NumericalVector& right) const
        {
            NumericalVector result(*this);
            result.data_+=right.data_; 
            #ifndef IAMNICO
             //   std::cout << "Nico you are using a slow operator be warned of the Numerical Vector, please try and rewrite to use += not +" << std::endl;
            #endif
            return result;
        }
        
       NumericalVector operator- (const NumericalVector& right)
        {
            NumericalVector result(*this);
            result.data_-=right.data_;
            #ifndef IAMNICO
             //   std::cout << "Nico you are using a slow operator be warned of the Numerical Vector, please try and rewrite to use -= not -" << std::endl;
            #endif
            return result;
        }
        
        NumericalVector operator- (const NumericalVector& right) const
        {
            NumericalVector result(*this);
            result.data_-=right.data_;
            #ifndef IAMNICO
              //  std::cout << "Nico you are using a slow operator be warned of the Numerical Vector, please try and rewrite to use -= not -" << std::endl;
            #endif
            return result;
        }
        
        NumericalVector operator* (const double& right)
        {
            NumericalVector result(*this);
            result.data_*=right;
            #ifndef IAMNICO
             //   std::cout << "Nico You are using a slow operator be warned of the Numerical Vector, please try and rewrite to use *= not *" << std::endl;
            #endif
            return result;
        }
        
        NumericalVector operator* (const double& right) const
        {
            NumericalVector result(*this);
            result.data_*=right;
            #ifndef IAMNICO
              //  std::cout << "Nico you are using a slow operator be warned of the Numerical Vector, please try and rewrite to use *= not *" << std::endl;
            #endif
            return result;
        }
        
        NumericalVector& operator/= (const double& right)
        {
            data_/=right;
            return *this;
        }
        
        NumericalVector operator/ (const double& right)
        {
            NumericalVector result(*this);
            return (result/=right);
            
        }
        
        void axpy(double a, const NumericalVector& x)
        {
            
            unsigned int size=data_.size();
            
            unsigned int i_one=1;
            
            
            daxpy_(&size, &a, &((*(const_cast<NumericalVector *> (&x)))[0]), &i_one, &((*this)[0]) , &i_one);
            
            
        }
        
        bool operator== (const NumericalVector& right)
        {
            for (int i = 0; i < data_.size(); ++i)
            {
                if (data_[i] != right.data_[i]) return false;
            }
            return true;
        }

        bool operator== (const NumericalVector& right) const
        {
            for (int i = 0; i < data_.size(); ++i)
            {
                if (data_[i] != right.data_[i]) return false;
            }
            return true;
        }

        NumericalVector& operator+= (const NumericalVector& right){data_+=right.data_; return *this;}
        
        NumericalVector& operator-= (const NumericalVector& right){data_-=right.data_; return *this;}
        
        NumericalVector& operator*= (const double& right){data_*=right; return *this;}
        
        inline double& operator[] (const unsigned int n) {return data_[n];}
        
        inline const double&  operator[] (const unsigned int n) const {return data_[n];}
        
        inline double& operator() (const unsigned int n) {return data_[n];}
        
        inline const double&  operator() (const unsigned int n) const {return data_[n];}
        
        int size() const {return data_.size();}
        
        int size() {return data_.size();}
        
        friend NumericalVector operator*(const double& left, const NumericalVector& right)
        {
            NumericalVector result(right);
            result.data_*=left;
            return result;
        }
        
        friend NumericalVector   operator-(const NumericalVector& right){return NumericalVector(right * -1.0);}
 
        friend ostream& operator<<(ostream& os, const NumericalVector& A)
        {
            os<< '[';
            for (int i=0;i<A.data_.size()-1;i++){os<< A(i)<<',';}
            os<< A(A.data_.size()-1)<<']';
            return os;
        }
        
        
   
    private:

        valarray<double> data_;
        
    };
    
}
#endif
