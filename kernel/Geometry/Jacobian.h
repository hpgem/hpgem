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

#ifndef JACOBIAN_HPP_
#define JACOBIAN_HPP_

#include "LinearAlgebra/SmallMatrix.h"

namespace Geometry
{
    /**
     * Jacobian of a mapping from dimFrom to dimTo. The standardized layout for Jacobians means that this has dimTo rows and dimFrom columns.
     * This may mean the template parameters appear to be in reverse order if you like to treat a Jacobian as any other matrix
     */
    template<std::size_t dimFrom, std::size_t dimTo>
    class Jacobian : public LinearAlgebra::SmallMatrix<dimTo, dimFrom>
    {

    public:
        // Constructors.
        Jacobian();
        Jacobian(const Jacobian& jacobian) : LinearAlgebra::SmallMatrix<dimTo, dimFrom>(jacobian) { }

        Jacobian(const LinearAlgebra::SmallMatrix<dimTo, dimFrom>& matrix)
                : LinearAlgebra::SmallMatrix<dimTo, dimFrom>(matrix)
        {
        }
        
        double determinant() const;
        
        /*! ConcatenatedMapping has to be able to do a matrix product on the
         Jacobians of two (successively applied) mappings. Therefore we provide
         the function multiplyJacobiansInto. */
        //version for concatenatedMapping(RefFace->RefFace->RefElem)
        Jacobian multiplyJacobiansInto(const Jacobian<dimFrom, dimFrom>& jac2)
        {
            return this->operator*(jac2);            
        }
        
        //version for computeNormalVector (RefFace->RefElem->PhysElem)
        Jacobian<dimFrom - 1, dimFrom> multiplyJacobiansInto(const Jacobian<dimFrom - 1, dimFrom>& jac2)
        {
            logger.assert(dimFrom == dimTo, "One of the matrixes must be square");
            return this->operator*(jac2);
        }

    };

}

#include "Jacobian_Impl.h"

#endif /* JACOBIAN_HPP_ */
