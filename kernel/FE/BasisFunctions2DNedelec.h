/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2014, University of Twente
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors
 may be used to endorse or promote products derived from this software without
 specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef HPGEM_KERNEL_BASISFUNCTIONS2DNEDELEC_H
#define HPGEM_KERNEL_BASISFUNCTIONS2DNEDELEC_H

#include <vector>
#include <memory>

namespace hpgem {

namespace FE {

class BaseBasisFunction;
class BasisFunctionSet;

/**
 * Implementation of Nedelec functions based on the basis for Raviart-Thomas
 * (RT) basis functions from [1]. We use the relation between H(div) and H(curl)
 * vector fields in 2D. If F = (Fx, Fy) is a vector field in H(div) then
 * rotating the vectors by 90 degrees (Fy, -Fx) one gets a vector field in
 * H(curl). Similarly, rotating the fields of the RT basis functions, one gets a
 * basis for the Nedelec space.
 *
 * As noted in [1] we can split the basis functions in two categories, edge and
 * internal. With the slight difference that we use order = p = k+1 (paper) we
 * have:
 * - 3*p edge basis functions, one per edge
 * - (p-1)*p internal basis functions
 * For a total of p*(p+2) basis functions.
 *
 * [1] Computational bases for RT_k and BDM_k on triangles
 * VJ Ervin, Computers & Mathematics with Applications 64 (8) 2764--2774
 */

/**
 * Create BasisFunctionSet for DG (elementwise) 2D Nedelec functions.
 * @param order The order of the basis functions (highest polynomial degree)
 * @return The set of basis functions, user owns the pointer.
 */
BasisFunctionSet* createDGBasisFunctionSet2DNedelec(std::size_t order);

/**
 * Create BasisFunctionSet for DG (elementwise) 2D Nedelec functions.
 * @param order The order of the basis functions (highest polynomial degree)
 * @return A vector with the basis functions. The user owns all the pointers.
 */
std::vector<BaseBasisFunction*> createDGBasisFunctions2DNedelec(
    std::size_t order);
}  // namespace FE

}  // namespace hpgem

#endif  // HPGEM_KERNEL_BASISFUNCTIONS2DNEDELEC_H
