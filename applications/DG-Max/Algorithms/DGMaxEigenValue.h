/*
This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
http://hpgem.org/about-the-code/team

This code is distributed using BSD 3-Clause License. A copy of which can found below.


Copyright (c) 2018, Univesity of Twenete
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef ALGORITHMS_DGMAXEIGENVALUE_h
#define ALGORITHMS_DGMAXEIGENVALUE_h

#include "../ProblemTypes/EigenValueProblem.h"

#include "../BaseExtended.h"
#include "DGMaxDiscretization.h"

#include "Utilities/GlobalIndexing.h"

#include <slepceps.h>

//TODO: It might be better to call this differently
class DGMaxEigenValue
{

public:

    DGMaxEigenValue(hpGemUIExtentions& base);
    void solve(const EigenValueProblem<DIM>& input, double stab);
    // TODO: A nice wrapper of EPS that does RAII would be nicer
    EPS createEigenSolver();
    void destroyEigenSolver(EPS& eps);


private:
    void initializeMatrices(double stab);
    void makeShiftMatrix(const Base::MeshManipulator<DIM>& mesh, const Utilities::GlobalIndexing& indexing,
                         const LinearAlgebra::SmallVector<DIM>& direction, Vec& waveVecMatrix) const;

    void extractEigenValues(const EPS& solver, std::vector<PetscScalar>& result);

    std::vector<Base::Face*> findPeriodicBoundaryFaces() const;
    LinearAlgebra::SmallVector<DIM> boundaryFaceShift(const Base::Face *face) const;

    hpGemUIExtentions& base_;
    DGMaxDiscretization discretization_;
};


#endif //ALGORITHMS_DGMAXEIGENVALUE_h
