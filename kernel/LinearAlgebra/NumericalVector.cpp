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
    
    #ifdef LA_STL_VECTOR
        /// \bug This constrcutor could be faster
        NumericalVector::NumericalVector(const double array[], int size)
        {
            data_.resize(size);
            for (int i=0;i<size;i++)
            {
                data_[i]=array[i];
            }
        }
    #else
        NumericalVector::NumericalVector(const double array[], int size) : data_(array,size){}
    #endif
    
    void NumericalVector::resize(unsigned int size) { if(size!=data_.size())data_.resize(size); }
    
    NumericalVector& NumericalVector::operator= (const NumericalVector& right){data_=right.data_; return *this;}
    
    NumericalVector NumericalVector::operator+ (const NumericalVector& right)
    {
        NumericalVector result(*this);
        #ifdef LA_STL_VECTOR
            /// \bug operator not implemented
        #else
            result.data_+=right.data_;
        #endif
        
        return result;
    }
    
   
    NumericalVector NumericalVector::operator+ (const NumericalVector& right) const
    {
        NumericalVector result(*this);
        #ifdef LA_STL_VECTOR
            /// \bug operator not implemented
        #else
            result.data_+=right.data_;
        #endif
        
        return result;
    }
    
    NumericalVector NumericalVector::operator- (const NumericalVector& right)
    {
        NumericalVector result(*this);
        #ifdef LA_STL_VECTOR
            /// \bug operator not implemeted
        #else
            result.data_-=right.data_;
        #endif
        return result;
    }
    
    NumericalVector NumericalVector::operator- (const NumericalVector& right) const
    {
        NumericalVector result(*this);
        #ifdef LA_STL_VECTOR
            /// \bug operator not implemeted
        #else
            result.data_-=right.data_;
        #endif
        return result;
    }
    
    NumericalVector NumericalVector::operator* (const double& right)
    {
        NumericalVector result(*this);
        #ifdef LA_STL_VECTOR
            /// \bug operator not implemented
        #else
            result.data_*=right;
        #endif
        
        return result;
    }
    
    NumericalVector NumericalVector::operator* (const double& right) const
    {
        NumericalVector result(*this);
        #ifdef LA_STL_VECTOR
            /// \bug operator not implemented
        #else
            result.data_*=right;
        #endif

        return result;
    }
    
    double NumericalVector::operator* (const NumericalVector& right) const
    {
        ///\TODO replace with BLAS (I dont know where to find them)
        #ifdef LA_STL_VECTOR
            /// \bug operator not implemented
            return 0;
        #else
            return (data_*right.data_).sum();
        #endif
    }
    
    NumericalVector& NumericalVector::operator/= (const double& right)
    {
        #ifdef LA_STL_VECTOR
            /// \bug operator not implemented
        #else
            data_/=right;
        #endif
        
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
    	if(data_.size()!=right.data_.size())
    		return false;
        for (int i = 0; i < data_.size(); ++i)
        {
            if (data_[i] != right.data_[i]) return false;
        }
        return true;
    }
    
    bool NumericalVector::operator== (const NumericalVector& right) const
    {
    	int n=data_.size();
    	if(n!=right.data_.size())
    		return false;
        for (int i = 0; i < n; ++i)
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
    

    NumericalVector& NumericalVector::operator+= (const NumericalVector& right)
        {
            #ifdef LA_STL_VECTOR
            /// \bug operator not implemented
            #else
                data_+=right.data_;
            #endif
            return *this;
        }
    
    NumericalVector& NumericalVector::operator-= (const NumericalVector& right)
        {
            #ifdef LA_STL_VECTOR
                /// \bug operator not implemented
            #else
                data_-=right.data_;
            #endif
            return *this;
        }
    
    NumericalVector& NumericalVector::operator*= (const double& right)
        {
            #ifdef LA_STL_VECTOR
                /// \bug operator not implemented
            #else
                data_*=right;
            #endif
            return *this;
        }
    
    double& NumericalVector::operator[] (const unsigned int n) {return data_[n];}
    
    NumericalVector operator*(const double& left, const NumericalVector& right)
    {
        NumericalVector result(right);
        #ifdef LA_STL_VECTOR
            /// \bug operator not implemented
        #else
            result.data_*=left;
        #endif
        return result;
    }
    
    NumericalVector   operator-(const NumericalVector& right){return NumericalVector(right * -1.0);}
    
    
    std::ostream& operator<<(std::ostream& os, const NumericalVector& A)
    {
        os<< '[';
        for (int i=0;i<A.data_.size()-1;i++){os<< A(i)<<',';}
        os<< A(A.data_.size()-1)<<']';
        return os;
    }
    
    

    
    
    
}
