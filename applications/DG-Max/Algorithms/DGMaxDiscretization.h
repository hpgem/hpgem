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
#ifndef HPGEM_APP_DGMAXDISCRETIZATION_H
#define HPGEM_APP_DGMAXDISCRETIZATION_H

#include <functional>
#include <map>
#include <memory>
#include <set>

#include "Logger.h"
#include "Base/Element.h"
#include "Base/Face.h"
#include "Base/PhysicalElement.h"
#include "Base/PhysicalFace.h"
#include "Base/MeshManipulator.h"
#include "Geometry/PointPhysical.h"
#include "Integration/ElementIntegral.h"
#include "Integration/FaceIntegral.h"
#include "ProblemTypes/BoundaryConditionType.h"

using namespace hpgem;

/// Dimension independent constants of DGMaxDiscretization
class DGMaxDiscretizationBase {
   public:
    /**
     * Processed mass matrix, contents depends on the rescaling approach.
     */
    static const std::size_t MASS_MATRIX_ID = 0;
    static const std::size_t STIFFNESS_MATRIX_ID = 1;
    static const std::size_t PROJECTOR_MATRIX_ID = 2;
    static const std::size_t INITIAL_CONDITION_VECTOR_ID = 0;
    static const std::size_t INITIAL_CONDITION_DERIVATIVE_VECTOR_ID = 1;
    static const std::size_t SOURCE_TERM_VECTOR_ID = 2;

    static const std::size_t FACE_MATRIX_ID = 0;
    static const std::size_t FACE_VECTOR_ID = 0;

    enum NormType { L2, HCurl, DG };

    /**
     * Several options for how matrices are rescaled to aid the using algorithm.
     */
    enum MassMatrixHandling {
        /**
         * Compute the mass matrix, no extra steps
         */
        NORMAL,
        /**
         * Compute the inverse of the mass matrix and store it as
         * MASS_MATRIX_ID.
         */
        INVERT,
        /**
         * Symmetrically rescale the matrices. The effect is the same as if the
         * basis functions were orthonormalized with respect to the L2-epsilon
         * inner product.
         *
         * The mass matrix M is factored as LL^H = M. The following
         * transformations are done:
         *  - Instead of the mass matrix its factor L is stored
         *  - The stiffness matrix S is rescaled to L^{-1} S L^{-H}
         *  - The projector B is rescaled to B L^{-H}
         *  - Any resulting vector y needs to be unscaled y = L^H x, to obtain
         * basis function coefficients x.
         */
        ORTHOGONALIZE
    };
};

template <std::size_t DIM>
class DGMaxDiscretization : public DGMaxDiscretizationBase {
   public:
    using PointPhysicalT = hpgem::Geometry::PointPhysical<DIM>;
    using InputFunction =
        std::function<LinearAlgebra::SmallVector<DIM>(const PointPhysicalT&)>;
    using FaceInputFunction = std::function<LinearAlgebra::SmallVector<DIM>(
        Base::PhysicalFace<DIM>&)>;
    using TimeFunction = std::function<LinearAlgebra::SmallVector<DIM>(
        const PointPhysicalT&, double)>;

    /// Computed fields of the solution
    struct Fields {
        /// Real part of the field
        LinearAlgebra::SmallVector<DIM> realEField;
        /// Imaginary part of the field
        LinearAlgebra::SmallVector<DIM> imagEField;
    };

    DGMaxDiscretization(bool includeProjector = false);

    void setMatrixHandling(MassMatrixHandling matrixHandling) {
        matrixHandling_ = matrixHandling;
    }

    void setBoundaryIndicator(DGMax::BoundaryConditionIndicator indicator) {
        boundaryIndicator_ = indicator;
    }

    void initializeBasisFunctions(Base::MeshManipulator<DIM>& mesh,
                                  std::size_t order);

    /**
     * Compute element matrices and vectors
     * @param mesh The mesh
     * @param massMatrix The possible transformation applied to the mass matrix
     * @param elementVectors The element vectors to compute as mapping from id
     * to the function to use.
     */
    void computeElementIntegrands(
        Base::MeshManipulator<DIM>& mesh,
        const std::map<std::size_t, InputFunction>& elementVectors);
    void computeFaceIntegrals(
        Base::MeshManipulator<DIM>& mesh,
        const std::map<std::size_t, FaceInputFunction>& boundaryVectors,
        double stab);

    static std::string normName(NormType norm) {
        switch (norm) {
            case L2:
                return "L2";
            case HCurl:
                return "HCurl";
            case DG:
                return "DG";
            default:
                logger.assert_debug(false, "Unknown norm type %", norm);
                return "Unknown";
        }
    }

    std::map<NormType, double> computeError(Base::MeshManipulator<DIM>& mesh,
                                            std::size_t timeVector,
                                            InputFunction electricField,
                                            InputFunction electricFieldCurl,
                                            std::set<NormType> norms);

    Fields computeFields(
        const Base::Element* element, const Geometry::PointReference<DIM>& p,
        const LinearAlgebra::MiddleSizeVector& coefficients) const;
    LinearAlgebra::SmallVector<DIM> computeField(
        const Base::Element* element, const Geometry::PointReference<DIM>& p,
        const LinearAlgebra::MiddleSizeVector& coefficients) const;
    LinearAlgebra::SmallVector<DIM> computeCurlField(
        const Base::Element* element, const Geometry::PointReference<DIM>& p,
        const LinearAlgebra::MiddleSizeVector& coefficients) const;

   private:
    /**
     * Compute the element local matrices
     */
    void computeElementMatrices(Base::Element* element);
    /**
     * Post process the element local matrices based on matrixHandling_
     */
    void postProcessElementMatrices(Base::Element* element) const;

    // Element vector integrand for the source term
    void elementInnerProduct(Base::PhysicalElement<DIM>& el,
                             const InputFunction& function,
                             LinearAlgebra::MiddleSizeVector& ret) const;

    void computeFaceMatrix(Base::Face* face, double stab);
    void postProcessFaceMatrices(Base::Face* face) const;

    // The face vector integrand.
    void faceVector(Base::PhysicalFace<DIM>& fa,
                    const FaceInputFunction& boundaryCondition,
                    LinearAlgebra::MiddleSizeVector& ret,
                    DGMax::BoundaryConditionType bct, double stab) const;

    // TODO: Replace this by a better type than SmallVector<2>.
    LinearAlgebra::SmallVector<2> elementErrorIntegrand(
        Base::PhysicalElement<DIM>& el, bool computeCurl,
        std::size_t timeVector, InputFunction exactValues,
        InputFunction curlValues) const;
    double faceErrorIntegrand(Base::PhysicalFace<DIM>& fa,
                              std::size_t timeVector,
                              InputFunction exactValue) const;

    const bool includeProjector_;
    MassMatrixHandling matrixHandling_;
    DGMax::BoundaryConditionIndicator boundaryIndicator_;

    std::vector<std::shared_ptr<Base::CoordinateTransformation<DIM>>>
        transforms_;
    Integration::ElementIntegral<DIM> elementIntegrator_;
    Integration::FaceIntegral<DIM> faceIntegrator_;
};

#endif  // HPGEM_APP_DGMAXDISCRETIZATION_H
