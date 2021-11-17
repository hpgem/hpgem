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
#ifndef HPGEM_PROBLEMMODE_H
#define HPGEM_PROBLEMMODE_H

namespace DGMax {

/**
 * \brief Specifies the field to use: E or H
 *
 * For the second order (time harmonic) Maxwell problem one has a problem in
 * either the electric field E or the H-field, where the other is eliminated.
 *
 * The resulting problems are (up to boundary conditions)
 *  Curl (mu^{-1} Curl E) - epsilon omega^2 E = J = - i omega j,
 *                              Div epsilon E = 0,
 * with Div J = 0, for the E-field and
 *  Curl (epsilon^{-1} Curl H) - mu omega^2 H = Curl j,
 *                                   Div mu H = 0,
 * for the H-field.
 * Notice that the differences between these formulations are minor:
 *  E is replaced by H,
 *  mu and epsilon switch,
 *  the source term changes,
 *  (boundary conditions may also change).
 * Thus the same discretization could be used for both problems (though the
 * analysis etc. may differ). This enumeration allows specifying which of the
 * two problems is used/assumed.
 */
enum class ProblemField { ELECTRIC_FIELD, MAGNETIC_FIELD };

}  // namespace DGMax

#endif  // HPGEM_PROBLEMMODE_H
