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

#ifndef SMALLVECTOR_H_
#define SMALLVECTOR_H_

#include "Logger.h"
#include <array>
#include <cmath>
#include "MiddleSizeVector.h"
#include <algorithm>
#include <numeric>

namespace LinearAlgebra
{
    class MiddleSizeVector;
    /// \class NumericalVector
    /// \brief This is a vector of doubles
    ///
    /// \details
    /// This implements a vector of doubles and all the standard operators for it.
    /// Note it is encapulating a valarray for its data storage.
    template<std::size_t nRows>
    class SmallVector
    {

    public:

        SmallVector()
            : data_()
        {
        }

        /*SmallVector(std::array<double, nRows> t)
            : data_(t)
        {
        }*/

        SmallVector(const SmallVector& other)
            : data_(other.data_)
        {
        }

        SmallVector(const MiddleSizeVector& other)
            : data_()
        {
            logger.assert(other.size() == nRows, "Cannot construct a vector of size % from a vector of size %", nRows, other.size());
            std::copy(other.data(), other.data() + nRows, data_.begin());
        }

        SmallVector(SmallVector&& other)
            : data_(std::move(other.data_))
        {
        }

        SmallVector(const double array[])
            : data_()
        {
            std::copy(array, array + nRows, data_.begin());
        }

        SmallVector& operator=(const SmallVector& right)
        {
            std::copy(right.data_.begin(), right.data_.end(), data_.begin());
            return *this;
        }

        SmallVector& operator=(const std::array<double, nRows> l)
        {
            std::copy(l.begin(), l.end(), data_.begin());
            return *this;
        }

        SmallVector operator+(const SmallVector& right) const
        {
            SmallVector result;
            std::transform(data_.begin(), data_.end(), right.data_.begin(), result.data_.begin(), std::plus<double>());
            return result;
        }

        SmallVector operator-(const SmallVector& right) const
        {
            SmallVector result;
            std::transform(data_.begin(), data_.end(), right.data_.begin(), result.data_.begin(), std::minus<double>());
            return result;
        }

        SmallVector operator*(const double& right) const
        {
            SmallVector result;
            std::transform(data_.begin(), data_.end(), result.data_.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, right));
            return result;
        }

        ///Computes inner product between two vectors.
        double operator*(const SmallVector& right) const
        {
            return std::inner_product(right.data_.begin(), right.data_.end(), data_.begin(), 0.);
        }

        SmallVector& operator/=(const double& right)
        {
            std::transform(data_.begin(), data_.end(), data_.begin(), std::bind(std::divides<double>(), std::placeholders::_1, right));
            return *this;
        }

        SmallVector operator/(const double& right) const
        {
            SmallVector result;
            std::transform(data_.begin(), data_.end(), result.data_.begin(), std::bind(std::divides<double>(), std::placeholders::_1, right));
            return result;
        }

        void axpy(double a, const SmallVector& x)
        {
            for(std::size_t i = 0; i < nRows; ++i)
            {
                data_[i] += a * x[i];
            }
        }

        /// This function is dangerous to use, since it compares doubles without
        /// a tolerance interval to see if they are equal.
        bool operator==(const SmallVector& right) const
        {
            for(std::size_t i = 0; i < nRows; ++i)
            {
                if(data_[i] != right[i])
                {
                    return false;
                }
            }
            return true;
        }

        /// This function is dangerous to use, since it compares doubles without
        /// a tolerance interval to see if they are equal.
        bool operator<(const SmallVector& right) const
        {
            for(std::size_t i = 0; i < nRows; ++i)
            {
                if(data_[i] < right[i])
                {
                    return true;
                }
                if(data_[i] > right[i])
                {
                    return false;
                }
            }
            return false;
        }

        SmallVector& operator+=(const SmallVector& right)
        {
            std::transform(data_.begin(), data_.end(), right.data_.begin(), data_.begin(), std::plus<double>());
            return *this;
        }

        SmallVector& operator-=(const SmallVector& right)
        {
            std::transform(data_.begin(), data_.end(), right.data_.begin(), data_.begin(), std::minus<double>());
            return *this;
        }

        SmallVector& operator*=(const double& right)
        {
            std::transform(data_.begin(), data_.end(), data_.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, right));
            return *this;
        }

        double& operator[](std::size_t n)
        {
            logger.assert(n < nRows, "Requested entry %, but there are only % entries", n, nRows);
            return data_[n];
        }

        const double& operator[](std::size_t n) const
        {
            logger.assert(n < nRows, "Requested entry %, but there are only % entries", n, nRows);
            return data_[n];
        }

        double& operator()(std::size_t n)
        {
            logger.assert(n < nRows, "Requested entry %, but there are only % entries", n, nRows);
            return data_[n];
        }

        const double& operator()(std::size_t n) const
        {
            logger.assert(n < nRows, "Requested entry %, but there are only % entries", n, nRows);
            return data_[n];
        }

        std::size_t size() const
        {
            return nRows;
        }
        const double* data() const
        {
            return data_.data();
        }

        double* data()
        {
            return data_.data();
        }

        SmallVector operator-() const
        {
            return *this * -1.;
        }

    private:
        std::array<double, nRows> data_;

    };

    template<std::size_t nRows>
    SmallVector<nRows> operator*(const double& left, const SmallVector<nRows>& right)
    {
        return right * left;
    }

    template<std::size_t nRows>
    std::ostream& operator<<(std::ostream& os, const SmallVector<nRows>& A)
    {
        os << "(";
        for(std::size_t i = 0; i < nRows; ++i)
        {
            os << A[i] << " ";
        }
        os << ")";
        return os;
    }

    template<std::size_t nRows>
    MiddleSizeVector::MiddleSizeVector(const SmallVector<nRows>& other)
        : data_(other.data(), other.data() + nRows)
    {
        logger(WARN, "Constructing middle size vector from small vector, consider using small vectors everywhere for fixed length vectors of size <= 4");
    }

}
#endif
