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

#include <functional>
#include <map>
#include <memory>
#include <set>

#include "Logger.h"

using namespace hpgem;


#ifndef HPGEM_APP_DGMAXDISCRETIZATION_H
#define HPGEM_APP_DGMAXDISCRETIZATION_H

// Forward definitions
namespace Base {
class Element;
template <std::size_t>
class MeshManipulator;
template <std::size_t DIM>
class PhysicalElement;
template <std::size_t DIM>
class PhysicalFace;
}  // namespace Base

namespace LinearAlgebra {
class MiddleSizeMatrix;
class MiddleSizeVector;

template <std::size_t DIM>
class SmallVector;
}  // namespace LinearAlgebra

namespace Geometry {
template <std::size_t DIM>
class PointReference;
template <std::size_t DIM>
class PointPhysical;
}  // namespace Geometry

template <std::size_t DIM>
class DGMaxDiscretization {
   public:
    static const std::size_t MASS_MATRIX_ID = 0;
    static const std::size_t STIFFNESS_MATRIX_ID = 1;
    static const std::size_t INITIAL_CONDITION_VECTOR_ID = 0;
    static const std::size_t INITIAL_CONDITION_DERIVATIVE_VECTOR_ID = 1;
    static const std::size_t SOURCE_TERM_VECTOR_ID = 2;

    static const std::size_t FACE_MATRIX_ID = 0;
    static const std::size_t FACE_VECTOR_ID = 0;

    using PointPhysicalT = Geometry::PointPhysical<DIM>;
    using InputFunction = std::function<void(const PointPhysicalT&,
                                             LinearAlgebra::SmallVector<DIM>&)>;
    using FaceInputFunction =
        std::function<void(const PointPhysicalT&, Base::PhysicalFace<DIM>&,
                           LinearAlgebra::SmallVector<DIM>&)>;
    using TimeFunction = std::function<void(const PointPhysicalT&, double,
                                            LinearAlgebra::SmallVector<DIM>&)>;

    void initializeBasisFunctions(Base::MeshManipulator<DIM>& mesh,
                                  std::size_t order);

    void computeElementIntegrands(
        Base::MeshManipulator<DIM>& mesh, bool invertMassMatrix,
        const InputFunction& sourceTerm, const InputFunction& initialCondition,
        const InputFunction& initialConditionDerivative) const;
    void computeFaceIntegrals(Base::MeshManipulator<DIM>& mesh,
                              const FaceInputFunction& boundaryCondition,
                              double stab) const;

    enum NormType { L2, HCurl, DG };
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
                                            std::set<NormType> norms) const;

    LinearAlgebra::SmallVector<DIM> computeField(
        const Base::Element* element, const Geometry::PointReference<DIM>& p,
        const LinearAlgebra::MiddleSizeVector& coefficients) const;
    LinearAlgebra::SmallVector<DIM> computeCurlField(
        const Base::Element* element, const Geometry::PointReference<DIM>& p,
        const LinearAlgebra::MiddleSizeVector& coefficients) const;

   private:
    // Mass matrix of the element
    void elementMassMatrix(Base::PhysicalElement<DIM>& el,
                           LinearAlgebra::MiddleSizeMatrix& ret) const;
    // Stiffness matrix, TODO: Is this the real stiffness or is this more
    // something like curl-stiffness?
    void elementStiffnessMatrix(Base::PhysicalElement<DIM>& el,
                                LinearAlgebra::MiddleSizeMatrix& ret) const;
    // Element vector integrand for the source term
    void elementInnerProduct(Base::PhysicalElement<DIM>& el,
                             const InputFunction& function,
                             LinearAlgebra::MiddleSizeVector& ret) const;

    // TODO: when this code can be tested, merge the two functions
    // The 'standard part' in the face matrix (jump-avg + avg-jump)
    void faceMatrix(Base::PhysicalFace<DIM>& fa,
                    LinearAlgebra::MiddleSizeMatrix& ret) const;
    // The penalty term in the face matrix
    void facePenaltyMatrix(Base::PhysicalFace<DIM>& fa,
                           LinearAlgebra::MiddleSizeMatrix& ret,
                           double stab) const;
    // The face vector integrand.
    void faceVector(Base::PhysicalFace<DIM>& fa,
                    const FaceInputFunction& boundaryCondition,
                    LinearAlgebra::MiddleSizeVector& ret, double stab) const;

    // TODO: Replace this by a better type than SmallVector<2>.
    LinearAlgebra::SmallVector<2> elementErrorIntegrand(
        Base::PhysicalElement<DIM>& el, bool computeCurl,
        std::size_t timeVector, InputFunction exactValues,
        InputFunction curlValues) const;
    double faceErrorIntegrand(Base::PhysicalFace<DIM>& fa,
                              std::size_t timeVector,
                              InputFunction exactValue) const;
};

#endif  // HPGEM_APP_DGMAXDISCRETIZATION_H
