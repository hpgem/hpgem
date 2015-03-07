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

#include "Jacobian.h"
#include "Logger.h"

namespace Geometry
{
    class Jacobian;
    
    Jacobian::Jacobian(std::size_t dimFrom, std::size_t dimTo)
            : LinearAlgebra::Matrix(dimFrom, dimTo)
    {
    }
    
    Jacobian::Jacobian(const JacobianT& jacobian)
            : LinearAlgebra::Matrix(jacobian)
    {
    }
    
    /*void
     Jacobian::computeWedgeStuffVector(NumericalVector& p)const
     {
     NumericalVector&      v((NumericalVector&)p);
     
     const LinearAlgebra::Matrix& jac=*this;

     //cout << "jacobian="<<jac<<endl;

     jac.computeWedgeStuffVector(v);
     }*/

    /// The computation of Jacobians are harcoded up until 4D, to make it faster.
    double Jacobian::determinant() const
    {
        std::size_t dimTo(getNRows()), dimFrom(getNCols());
        
        logger.assert(dimFrom == dimTo, "Jacobian should be square to have a determinant!");
        
        switch (dimFrom)
        {
            case 0:
                return 1;
                break;
            case 1:
                return (*this)(0, 0);
                break;
            case 2:
                return (*this)(0, 0) * (*this)(1, 1) - (*this)(0, 1) * (*this)(1, 0);
                break;
                
            case 3:
                return (*this)(0, 0) * ((*this)(1, 1) * (*this)(2, 2) - (*this)(1, 2) * (*this)(2, 1)) - (*this)(0, 1) * ((*this)(1, 0) * (*this)(2, 2) - (*this)(2, 0) * (*this)(1, 2)) + (*this)(0, 2) * ((*this)(1, 0) * (*this)(2, 1) - (*this)(2, 0) * (*this)(1, 1));
                break;
                
            case 4:
                return ((*this)(3, 0) * (*this)(2, 1) * (*this)(0, 3) - (*this)(2, 0) * (*this)(3, 1) * (*this)(0, 3)) * (*this)(1, 2) + (-(*this)(3, 0) * (*this)(0, 3) * (*this)(2, 2) + (*this)(2, 0) * (*this)(0, 3) * (*this)(3, 2)) * (*this)(1, 1) + ((*this)(3, 1) * (*this)(0, 3) * (*this)(2, 2) - (*this)(2, 1) * (*this)(0, 3) * (*this)(3, 2)) * (*this)(1, 0) + (-(*this)(3, 0) * (*this)(2, 1) * (*this)(1, 3) + (*this)(2, 0) * (*this)(3, 1) * (*this)(1, 3) + (-(*this)(2, 0) * (*this)(3, 3) + (*this)(3, 0) * (*this)(2, 3)) * (*this)(1, 1) + ((*this)(2, 1) * (*this)(3, 3) - (*this)(3, 1) * (*this)(2, 3)) * (*this)(1, 0)) * (*this)(0, 2) + ((*this)(3, 0) * (*this)(1, 3) * (*this)(2, 2) - (*this)(2, 0) * (*this)(1, 3) * (*this)(3, 2) + ((*this)(2, 0) * (*this)(3, 3) - (*this)(3, 0) * (*this)(2, 3)) * (*this)(1, 2) + (-(*this)(2, 2) * (*this)(3, 3) + (*this)(2, 3) * (*this)(3, 2)) * (*this)(1, 0)) * (*this)(0, 1) + (-(*this)(3, 1) * (*this)(1, 3) * (*this)(2, 2) + (*this)(2, 1) * (*this)(1, 3) * (*this)(3, 2) + ((*this)(3, 1) * (*this)(2, 3) - (*this)(2, 1) * (*this)(3, 3)) * (*this)(1, 2) + (*this)(1, 1) * ((*this)(2, 2) * (*this)(3, 3) - (*this)(2, 3) * (*this)(3, 2))) * (*this)(0, 0);
                // ... says Maple; this can possibly be done more efficiently,
                // maybe even with LU (with pivoting, though...)
                break;
            default:
                logger(ERROR, "Computing the Jacobian for dimension % is not implemented", dimFrom);
                break;
        }
        return 0;
    }
}
;
