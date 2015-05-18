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

#include "MiddleSizeVector.h"

#include "Logger.h"

/**
 * \bug part of the fortran interface is passes integers as unsigned int and part of the fortran interface passes integers as int
 */
namespace LinearAlgebra
{
    
    extern "C"

    {
        
        ///This is the general scalar times vector + vector from blas, hence from blas level 1. Here we also use on a matrix by treating as a vector
        int daxpy_(unsigned int* N, double* DA, double* DX, unsigned int* INCX, double* DY, unsigned int* INCY);
    
    }
    
    MiddleSizeVector::MiddleSizeVector()
            : data_()
    {
    }
    
    MiddleSizeVector::MiddleSizeVector(std::size_t m)
            : data_(m)
    {
    }
    
    MiddleSizeVector::MiddleSizeVector(std::initializer_list<double> l)
            : data_(l)
    {
    }

    MiddleSizeVector::MiddleSizeVector(const MiddleSizeVector& other)
            : data_(other.data_)
    {
    }
    
    MiddleSizeVector::MiddleSizeVector(MiddleSizeVector&& other)
            : data_(std::move(other.data_))
    {
    }
    
#ifdef LA_STL_VECTOR
    MiddleSizeVector::MiddleSizeVector(const double array[], std::size_t size)
            : data_(array, array + size)
    {
    }

#ifdef HPGEM_USE_COMPLEX_PETSC
    MiddleSizeVector::MiddleSizeVector(const std::complex<double> array[], int size)
    {   
        data_.resize(size);

        for(std::size_t i = 0; i<size; i++)
        {   
            if(std::imag(array[i]) !=0)
            {   
                logger(INFO, "Imaginary part nonzero in the entry of vector");
            }
            data_[i]=std::real(array[i]);
        }

    }
#endif
#else
    MiddleSizeVector::MiddleSizeVector(const double array[], std::size_t size)
    : data_(array, size)
    {}
#endif
    
    void MiddleSizeVector::resize(std::size_t size)
    {
        if (size != data_.size())
        {
            data_.resize(size);
        }
    }
    
    MiddleSizeVector& MiddleSizeVector::operator=(const MiddleSizeVector& right)
    {
        data_ = right.data_;
        return *this;
    }
    
    MiddleSizeVector& MiddleSizeVector::operator=(const std::initializer_list<double> l)
    {
    	data_ = l;
    	return *this;
    }

    MiddleSizeVector MiddleSizeVector::operator+(const MiddleSizeVector& right) const
    {
        MiddleSizeVector result(*this);
        logger.assert(data_.size() == right.data_.size(), "Vectors don't have the same size");
#ifdef LA_STL_VECTOR
        for (std::size_t i = 0; i < data_.size(); i++)
            result.data_[i] += right.data_[i];
#else
        result.data_+=right.data_;
#endif
        
        return result;
    }
    
    MiddleSizeVector MiddleSizeVector::operator-(const MiddleSizeVector& right) const
    {
        MiddleSizeVector result(*this);
        logger.assert(data_.size() == right.data_.size(), "Vectors don't have the same size");
#ifdef LA_STL_VECTOR
        for (std::size_t i = 0; i < data_.size(); i++)
            result.data_[i] -= right.data_[i];
#else
        result.data_-=right.data_;
#endif
        return result;
    }
    
    MiddleSizeVector MiddleSizeVector::operator*(const double& right) const
    {
        MiddleSizeVector result(*this);
#ifdef LA_STL_VECTOR
        for (double& d : result.data_)
            d *= right;
#else
        result.data_*=right;
#endif
        
        return result;
    }
    
    double MiddleSizeVector::operator*(const MiddleSizeVector& right) const
    {
        ///\TODO replace with BLAS 
        logger.assert(data_.size() == right.data_.size(), "Vectors don't have equal length.");
#ifdef LA_STL_VECTOR
        double sum = 0;
        for (std::size_t i = 0; i < data_.size(); i++)
            sum += data_[i] * right.data_[i];
        return sum;
#else
        return (data_*right.data_).sum();
#endif
    }
    
    MiddleSizeVector& MiddleSizeVector::operator/=(const double& right)
    {
#ifdef LA_STL_VECTOR
        for (double& d : data_)
            d /= right;
#else
        data_/=right;
#endif
        
        return *this;
    }
    
    MiddleSizeVector MiddleSizeVector::operator/(const double& right) const
    {
        MiddleSizeVector result(*this);
        return (result /= right);
        
    }
    
    void MiddleSizeVector::axpy(double a, const MiddleSizeVector& x)
    {
        logger.assert(x.size() == data_.size(), "Vectors dont have the same size");
        unsigned int size = data_.size();
        
        unsigned int i_one = 1;
        
        daxpy_(&size, &a, const_cast<MiddleSizeVector *>(&x)->data_.data(), &i_one, ((*this).data_.data()), &i_one);
        
    }
    
    bool MiddleSizeVector::operator==(const MiddleSizeVector& right) const
    {
        return (data_ == right.data_);
    }
    
    bool MiddleSizeVector::operator<(const MiddleSizeVector& right) const
    {
        for (std::size_t i = 0; i < data_.size() && i < right.data_.size(); ++i)
        {
            if (data_[i] < right.data_[i])
            {
                return true;
            }
            if (data_[i] > right.data_[i])
            {
                return false;
            }
        }
        return false;
    }
    
    MiddleSizeVector& MiddleSizeVector::operator+=(const MiddleSizeVector& right)
    {
        logger.assert(data_.size() == right.data_.size(), "Vectors don't have the same size");
#ifdef LA_STL_VECTOR
        for (std::size_t i = 0; i < data_.size(); i++)
            data_[i] += right.data_[i];
#else
        data_+=right.data_;
#endif
        return *this;
    }
    
    MiddleSizeVector& MiddleSizeVector::operator-=(const MiddleSizeVector& right)
    {
        logger.assert(data_.size() == right.data_.size(), "Vectors don't have the same size");
#ifdef LA_STL_VECTOR
        for (std::size_t i = 0; i < data_.size(); i++)
            data_[i] -= right.data_[i];
        
#else
        data_-=right.data_;
#endif
        return *this;
    }
    
    MiddleSizeVector& MiddleSizeVector::operator*=(const double& right)
    {
#ifdef LA_STL_VECTOR
        for (double& d : data_)
            d *= right;
#else
        data_*=right;
#endif
        return *this;
    }
    
    double& MiddleSizeVector::operator[](std::size_t n)
    {
        logger.assert(n < data_.size(), "Requested entry %, but there are only % entries", n, data_.size());
        return data_[n];
    }
    
    MiddleSizeVector operator*(const double& left, const MiddleSizeVector& right)
    {
        MiddleSizeVector result(right);
#ifdef LA_STL_VECTOR
        for (double& d : result.data_)
            d *= left;
#else
        result.data_*=left;
#endif
        return result;
    }
    
    MiddleSizeVector operator-(const MiddleSizeVector& right)
    {
        return MiddleSizeVector(right * -1.0);
    }
    
    std::ostream& operator<<(std::ostream& os, const MiddleSizeVector& A)
    {
        os << '[';
        for (std::size_t i = 0; i < A.data_.size() - 1; i++)
        {
            os << A(i) << ',';
        }
        os << A(A.data_.size() - 1) << ']';
        return os;
    }

#ifdef HPGEM_USE_COMPLEX_PETSC

const std::complex<double>* MiddleSizeVector::data() const
{   
    static std::vector<std::complex<double>> new_Data(data_.size());

    for (std::size_t i = 0; i < data_.size(); i++)
    {   
        new_Data[i] = data_[i];
    }

    return new_Data.data();
}
#endif

}
