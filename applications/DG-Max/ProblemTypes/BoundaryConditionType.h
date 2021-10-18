/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2021, University of Twente
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
#ifndef HPGEM_BOUNDARYCONDITIONTYPE_H
#define HPGEM_BOUNDARYCONDITIONTYPE_H

#include <functional>
#include "DGMaxLogger.h"
#include "Base/Face.h"

namespace DGMax {

/**
 * Enumeration of the possible boundary conditions that could be used on a
 * boundary face.
 *
 * Note: Not all boundary conditions may be supported by all algorithms.
 */
enum class BoundaryConditionType {
    /**
     * Boundary condition of the form n x E = n x g_D, thus prescribing the
     * tangential part of the field on the boundary face.
     */
    DIRICHLET,
    /**
     * Boundary condition of the form n x (curl E) = n x g_N, thus prescribing
     * the tangential part of the curl on the boundary face.
     */
    NEUMANN,
    /**
     * Boundary condition of the form
     * n x (Curl E)  + i omega sqrt[epsilon_r mu_r] E_tangent = n x g_N,
     * corresponding to the Silver-Muller radiation condition applied to the
     * boundary (instead of at infinity). This is an Impedance boundary
     * condition with a fixed impedance. Mathematically it is of a Robin type.
     *
     * Note: The right hand side has a 'n x' to enforce that only the tangential
     * part of g_N is used.
     */
    SILVER_MULLER,
    /**
     * Dummy value for internal faces, should not be used for external faces.
     *
     * This is primarily intended so that we can assign each face a value
     */
    INTERNAL,
};

inline bool isNaturalBoundary(BoundaryConditionType type) {
    switch (type) {
        case BoundaryConditionType::DIRICHLET:
            return false;
        case BoundaryConditionType::NEUMANN:  // Fall through
        case BoundaryConditionType::SILVER_MULLER:
            return true;
        case BoundaryConditionType::INTERNAL:
            DGMaxLogger.assert_always(false,
                                      "Internal boundary to isNaturalBoundary");
            return false;
        default:
            DGMaxLogger.assert_always(false,
                                      "isNaturalBoundary not implemented for "
                                      "this boundary condition type");
            return false;
    }
}

using BoundaryConditionIndicator =
    typename std::function<BoundaryConditionType(const hpgem::Base::Face&)>;

}  // namespace DGMax

#endif  // HPGEM_BOUNDARYCONDITIONTYPE_H
