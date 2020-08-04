/*
This file forms part of hpGEM. This package has been developed over a number of
years by various people at the University of Twente and a full list of
contributors can be found at http://hpgem.org/about-the-code/team

This code is distributed using BSD 3-Clause License. A copy of which can found
below.


Copyright (c) 2018, Univesity of Twenete
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef HPGEM_APP_DIVDGMAXEIGENVALUE_H
#define HPGEM_APP_DIVDGMAXEIGENVALUE_H

#include <petscis.h>
#include <petscvec.h>
#include <slepceps.h>

#include "Utilities/GlobalIndexing.h"

#include "DivDGMaxDiscretization.h"
#include "ProblemTypes/AbstractEigenvalueResult.h"
#include "ProblemTypes/AbstractEigenvalueSolver.h"
#include "ProblemTypes/EigenvalueProblem.h"

using namespace hpgem;

template <std::size_t DIM>
class DivDGMaxEigenvalue : public AbstractEigenvalueSolver<DIM> {
   public:
    class Result : public AbstractEigenvalueResult<DIM> {
       public:
        Result(EigenvalueProblem<DIM> problem,
               std::vector<std::vector<PetscScalar>> eigenvalues);
        const EigenvalueProblem<DIM>& originalProblem() const final;
        const std::vector<double> frequencies(std::size_t point) const final;

       private:
        const EigenvalueProblem<DIM> problem_;
        const std::vector<std::vector<PetscScalar>> eigenvalues_;
    };

    DivDGMaxEigenvalue(Base::MeshManipulator<DIM>& mesh, std::size_t order,
                       typename DivDGMaxDiscretization<DIM>::Stab stab,
                       DivDGMaxDiscretizationBase::DivType divType =
                           DivDGMaxDiscretizationBase::DivType::DG);
    std::unique_ptr<AbstractEigenvalueResult<DIM>> solve(
        const EigenvalueProblem<DIM>& input) override;

   private:
    void extractEigenvalues(const EPS& solver,
                            std::vector<PetscScalar>& result) const;

    void makeShiftMatrix(LinearAlgebra::SmallVector<DIM>& direction,
                         const Utilities::GlobalIndexing& index,
                         Vec& waveVecMatrix);

    Base::MeshManipulator<DIM>& mesh_;
    typename DivDGMaxDiscretization<DIM>::Stab stab_;
    std::size_t order_;
    DivDGMaxDiscretizationBase::DivType divType_;
    DivDGMaxDiscretization<DIM> discretization;
};

#endif  // HPGEM_APP_DIVDGMAXEIGENVALUE_H
