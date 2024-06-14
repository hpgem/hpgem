/*
This file forms part of hpGEM. This package has been developed over a number of
years by various people at the University of Twente and a full list of
contributors can be found at http://hpgem.org/about-the-code/team

This code is distributed using BSD 3-Clause License. A copy of which can found
below.


Copyright (c) 2020, Univesity of Twenete
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
#ifndef HPGEM_APP_ABSTRACTEIGENVALUERESULT_H
#define HPGEM_APP_ABSTRACTEIGENVALUERESULT_H

#include <array>
#include <vector>
#include "LinearAlgebra/SmallVector.h"
#include "LinearAlgebra/MiddleSizeMatrix.h"
#include "Output/VTKSpecificTimeWriter.h"

#include "DGMaxLogger.h"

using namespace hpgem;

/**
 * Callback interface for an AbstractEigenvalueSolver to allow access to the
 * result of solving at a single K-point in a discretization independent way.
 *
 * @tparam DIM The dimension of the original problem.
 */
template <std::size_t DIM>
class AbstractEigenvalueResult {
   public:
    virtual ~AbstractEigenvalueResult() = default;

    /**
     * @return The computed eigenfrequencies.
     */
    virtual std::vector<double> getFrequencies() = 0;

    /**
     * @return Wavevector point of this result
     */
    virtual const LinearAlgebra::SmallVector<DIM>& getKPoint() const = 0;

    /**
     * Write the field (or fields) of a specific eigenvalue (if supported).
     *
     * When not supported should remain as is, printing a warning.
     *
     * @param eigenvalue The eigenvalue index
     * @param writer The writer to use for the output.
     */
    virtual void writeField(std::size_t eigenvalue,
                            Output::VTKSpecificTimeWriter<DIM>& writer) {
        DGMaxLogger(WARN, "Not supported for this result implementation.");
    };

    /**
     * @return Mesh file used for the computation.
     */
    virtual const Base::MeshManipulator<DIM>* getMesh() const = 0;

    /**
     * \brief Compute field overlap integrals
     *
     * Compute the field overlap integrals integral(E_prev epsilon E_cur)
     * between the fields E_cur of the current k-point and the field E_prev of
     * the previous k-point. The rows correspond to the fields of the current
     * k-point, while the columns correspond to the fields of the previous
     * k-point.
     *
     * By default outputs a matrix of size (0,0) if not supported.
     *
     * @return A matrix with the overlap values.
     */
    virtual LinearAlgebra::MiddleSizeMatrix computeFieldOverlap() const {
        DGMaxLogger(WARN, "Field overlap computation not supported");
        return LinearAlgebra::MiddleSizeMatrix(0, 0);
    };

    virtual bool supportsWaveVectorDerivatives() const { return false; }

    virtual std::array<LinearAlgebra::MiddleSizeMatrix, DIM>
        computeWaveVectorDerivatives() const {
        DGMaxLogger(WARN, "Wave vector derivatives not supported");
        // Empty matrices
        std::array<LinearAlgebra::MiddleSizeMatrix, DIM> result;
        for (std::size_t i = 0; i < DIM; ++i) {
            result[i] = LinearAlgebra::MiddleSizeMatrix(0, 0);
        }
        return result;
    }
};

#endif  // HPGEM_APP_ABSTRACTEIGENVALUERESULT_H
