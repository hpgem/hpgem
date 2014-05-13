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

#include "GlobalNamespaceLinearAlgebra.hpp"
#include <valarray>
#include <iostream>


#define IAMNICO

namespace LinearAlgebra
{
    
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
        
        NumericalVector();
        
        NumericalVector(int m);
        
        NumericalVector(const NumericalVector& other);
        
        NumericalVector(const double array[], int size);
        
        void resize(unsigned int size);

        NumericalVector& operator= (const NumericalVector& right);
        
        NumericalVector operator+ (const NumericalVector& right);
        
        NumericalVector operator+ (const NumericalVector& right) const;
        
        NumericalVector operator- (const NumericalVector& right);
        
        NumericalVector operator- (const NumericalVector& right) const;
        
        NumericalVector operator* (const double& right);
        
        NumericalVector operator* (const double& right) const;
        
        double operator* (const NumericalVector& right) const;

        NumericalVector& operator/= (const double& right);
        
        NumericalVector operator/ (const double& right);
        
        void axpy(double a, const NumericalVector& x);
        
        bool operator== (const NumericalVector& right);

        bool operator== (const NumericalVector& right) const;
        
        bool operator< (const NumericalVector& right) const;

        NumericalVector& operator+= (const NumericalVector& right);
        
        NumericalVector& operator-= (const NumericalVector& right);
        
        NumericalVector& operator*= (const double& right);
        
        double& operator[] (const unsigned int n);
        
        const double&  operator[] (const unsigned int n) const {return data_[n];}
        
        double& operator() (const unsigned int n) {return data_[n];}
        
        const double&  operator() (const unsigned int n) const {return data_[n];}
        
        int size() const {return data_.size();}
        
        int size() {return data_.size();}
        
        friend NumericalVector operator*(const double& left, const NumericalVector& right);
        
        friend NumericalVector   operator-(const NumericalVector& right);
 
        friend ostream& operator<<(ostream& os, const NumericalVector& A);
        
   
    private:

        valarray<double> data_;
        
    };
    
}
#endif
