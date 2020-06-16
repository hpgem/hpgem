/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2020, University of Twente
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

#include "DGMaxEVConvergenceTest.h"

#include "Base/ConfigurationData.h"
#include "Base/MeshManipulator.h"

#include "DGMaxLogger.h"
#include "DGMaxProgramUtils.h"

#include "Algorithms/DGMaxEigenvalue.h"

namespace DGMax {

template <std::size_t DIM>
std::unique_ptr<AbstractEigenvalueResult<DIM>>
    DGMaxEVConvergenceTest<DIM>::runInternal(std::size_t level) {

    logger.assert_always(level < meshFileNames_.size(), "No such mesh");

    std::size_t unknowns = solverConfig_.useProjector_ ? 2 : 1;
    Base::ConfigurationData configData(unknowns, 1);
    auto mesh = DGMax::readMesh<DIM>(
        meshFileNames_[level], &configData,
        [&](const Geometry::PointPhysical<DIM> &p) {
            return jelmerStructure(p, testCase_.getStructureId());
        });
    DGMaxLogger(INFO, "Loaded mesh % with % local elements.",
                meshFileNames_[level], mesh->getNumberOfElements());
    KSpacePath<DIM> path =
        KSpacePath<DIM>::singleStepPath(testCase_.getKPoint());
    EigenvalueProblem<DIM> input(path, testCase_.getNumberOfEigenvalues());

    DGMaxEigenvalue<DIM> solver(*mesh, this->order_, this->solverConfig_);
    return solver.solve(input);
}

// Template instantiation
template class DGMaxEVConvergenceTest<2>;
template class DGMaxEVConvergenceTest<3>;

};  // namespace DGMax
