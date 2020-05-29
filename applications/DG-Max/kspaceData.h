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

#include <vector>
#include "LinearAlgebra/SmallVector.h"
#include "LinearAlgebra/MiddleSizeMatrix.h"
#include "LinearAlgebra/MiddleSizeVector.h"

/**
 * The tesselation in k-space does not fit the hpGEM-framework because of its
 * continuous nature this is a specialized class to store the relevant
 * information and provide the relevant computations
 */
class KspaceData {
    /*

    std::vector<LinearAlgebra::NumericalVector> kpoints_, deltak_;
    std::vector<std::vector<double> > omegaAtKpoints_;
    std::vector<std::vector<int> > elements_;
    std::vector<std::vector<LinearAlgebra::NumericalVector> >
functionValuesAtKpoints_; int current_; int minimumsize_;

public:


     * constructor takes a desired number of k-point per direction and
constructs a tesselation of small tetrahedra
     * that can be used for linear interpolation. the actual number of k-points
where the eigenvalue problem needs to be solved
     * is (n)(n+1)(n+2)/6. makes sure the first k-point is the origin so no
matrix updates are needed for the first pass

    KspaceData(int pointsPerDirection);


     * returns true as long as there is still a point that needs computing

    bool hasNextPoint();


     * iterator over the k-points. Call this when you are done for the current
k-point and it gets you a update-vector

    LinearAlgebra::NumericalVector& nextPoint();


     * sets the found values of omega at the current k-point. Make sure you set
the same number of omega each k-point and that they are sorted in the same way
(accounting for band crossings)

    void setOmega(std::vector<double>& omega);


     * optional: sets the value of f from the expression
int(f*delta(omega-omega_n))dk

    void setFunctionValues(std::vector<LinearAlgebra::NumericalVector>&
functionValues);


     * returns the density of states for a given value omega

    void getIntegral(double omega, LinearAlgebra::NumericalVector& result);

    */
};
