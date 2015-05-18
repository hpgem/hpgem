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

#ifndef NumericalVectorHPP
#define NumericalVectorHPP

//This is derived from valarray or vector so import that information
#include "Logger.h"
#ifdef LA_STL_VECTOR
#include <vector>
#include <cmath>
#else
#include <valarray>
#endif
#include <iostream>
#include <complex>

namespace LinearAlgebra
{
    
#ifdef LA_STL_VECTOR
    using std::vector;
#else
    using std::valarray;
#endif
    
    /// \class NumericalVector
    /// \brief This is a vector of doubles
    ///
    /// \details
    /// This implements a vector of doubles and all the standard operators for it.
    /// Note it is encapulating a valarray for its data storage.
    class MiddleSizeVector
    {
        
    public:
        
        MiddleSizeVector();

        explicit MiddleSizeVector(std::size_t m);

        MiddleSizeVector(std::initializer_list<double> t);

        MiddleSizeVector(const MiddleSizeVector& other);
        
        MiddleSizeVector(MiddleSizeVector&& other);

        MiddleSizeVector(const double array[], std::size_t size);

        //Constructor to accommodate complex<double>
        MiddleSizeVector(const std::complex<double> array[], int size);

        void resize(std::size_t size);

        MiddleSizeVector& operator=(const MiddleSizeVector& right);

        MiddleSizeVector& operator=(const std::initializer_list<double> l);

        MiddleSizeVector operator+(const MiddleSizeVector& right) const;

        MiddleSizeVector operator-(const MiddleSizeVector& right) const;

        MiddleSizeVector operator*(const double& right) const;

        ///Computes inner product between two vectors.
        double operator*(const MiddleSizeVector& right) const;

        MiddleSizeVector& operator/=(const double& right);

        MiddleSizeVector operator/(const double& right) const;

        void axpy(double a, const MiddleSizeVector& x);

        /// This function is dangerous to use, since it compares doubles without 
        /// a tolerance interval to see if they are equal.
        /// Needs fixing if someone wants to use valarray.
        bool operator==(const MiddleSizeVector& right) const;

        /// This function is dangerous to use, since it compares doubles without
        /// a tolerance interval to see if they are equal.
        bool operator<(const MiddleSizeVector& right) const;

        MiddleSizeVector& operator+=(const MiddleSizeVector& right);

        MiddleSizeVector& operator-=(const MiddleSizeVector& right);

        MiddleSizeVector& operator*=(const double& right);

        double& operator[](std::size_t n);

        const double& operator[](std::size_t n) const
        {
            logger.assert(n < data_.size(), "Requested entry %, but there are only % entries", n, data_.size());
            return data_[n];
        }
        
        double& operator()(std::size_t n)
        {
            logger.assert(n < data_.size(), "Requested entry %, but there are only % entries", n, data_.size());
            return data_[n];
        }
        
        const double& operator()(std::size_t n) const
        {
            logger.assert(n < data_.size(), "Requested entry %, but there are only % entries", n, data_.size());
            return data_[n];
        }
        
        std::size_t size() const
        {
            return data_.size();
        }
#ifdef HPGEM_USE_COMPLEX_PETSC
        
        const std::complex<double>* data() const;

#else
        const double* data() const
        {
            return data_.data();
        }
        
        double* data()
        {
            return data_.data();
        }
        
#endif
        
        friend MiddleSizeVector operator*(const double& left, const MiddleSizeVector& right);

        friend MiddleSizeVector operator-(const MiddleSizeVector& right);

        friend std::ostream& operator<<(std::ostream& os, const MiddleSizeVector& A);

    private:
#ifdef LA_STL_VECTOR
        vector<double> data_;
#else
        valarray<double> data_;
#endif
        
    };

}
#endif
