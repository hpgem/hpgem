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
#ifndef HPGEM_ABSTRACTDISCRETIZATION_H
#define HPGEM_ABSTRACTDISCRETIZATION_H

#include <functional>
#include <map>

#include <Base/MeshManipulator.h>
#include <Base/PhysicalFace.h>
#include <Geometry/PointPhysical.h>
#include <Geometry/PointReference.h>
#include <LinearAlgebra/SmallVector.h>
#include <LinearAlgebra/MiddleSizeVector.h>
#include <Output/VTKSpecificTimeWriter.h>

#include "../ProblemTypes/FieldPattern.h"
#include "../ProblemTypes/BoundaryConditionType.h"

namespace DGMax {

class AbstractDiscretizationBase {
   public:
    /// Matrix corresponding to (epsilon E, v)
    static const constexpr std::size_t MASS_MATRIX_ID = 0;
    /// Matrix corresponding to (mu^{-1} Curl E, Curl v)
    static const constexpr std::size_t STIFFNESS_MATRIX_ID = 1;

    /// Face matrix from discretizing (mu^{-1} Curl E, Curl v)
    static const constexpr std::size_t FACE_STIFFNESS_MATRIX_ID = 0;
    /// Matrix for the impedance term(s)
    static const constexpr std::size_t FACE_IMPEDANCE_MATRIX_ID = 1;

    /// Element vector corresponding to (J, v)
    static const constexpr std::size_t ELEMENT_VECTOR_ID = 0;
    static const constexpr std::size_t INITIAL_CONDITION_VECTOR_ID = 1;
    static const constexpr std::size_t INITIAL_CONDITION_DERIVATIVE_VECTOR_ID =
        2;

    /// Face vector corresponding to the boundary conditions
    static const constexpr std::size_t FACE_VECTOR_ID = 0;

    enum class LocalIntegrals {
        /**
         * Compute both matrix and vector integrals
         */
        ALL,
        /**
         * Compute only the vector intergals
         */
        ONLY_VECTORS
    };

    virtual std::size_t getOrder() const = 0;
    virtual std::size_t getNumberOfUnknowns() const = 0;
    virtual std::size_t getNumberOfElementMatrices() const = 0;
    virtual std::size_t getNumberOfFaceMatrices() const = 0;
};

template <std::size_t dim>
class AbstractDiscretization : public AbstractDiscretizationBase {
   public:
    using PointPhysicalT = hpgem::Geometry::PointPhysical<dim>;
    using PointReferenceT = hpgem::Geometry::PointReference<dim>;
    using InputFunction = std::function<hpgem::LinearAlgebra::SmallVectorC<dim>(
        const Base::Element& element, const PointPhysicalT&)>;
    using FaceInputFunction =
        std::function<hpgem::LinearAlgebra::SmallVectorC<dim>(
            hpgem::Base::PhysicalFace<dim>&)>;
    using FieldT = hpgem::LinearAlgebra::SmallVectorC<dim>;

    virtual void initializeBasisFunctions(
        hpgem::Base::MeshManipulator<dim>& mesh) const = 0;

    void computeElementIntegrals(
        hpgem::Base::MeshManipulator<dim>& mesh,
        const std::map<std::size_t, InputFunction>& elementVectors,
        double dispersionOmega,
        LocalIntegrals integrals = LocalIntegrals::ALL) {
        computeElementIntegralsImpl(mesh, elementVectors, dispersionOmega,
                                    integrals);
    }

    void computeFaceIntegrals(
        hpgem::Base::MeshManipulator<dim>& mesh,
        const std::map<std::size_t, FaceInputFunction>& faceVectors,
        double dispersionOmega,
        BoundaryConditionIndicator boundaryIndicator =
            [](const hpgem::Base::Face&) {
                return BoundaryConditionType::DIRICHLET;
            },
        LocalIntegrals integrals = LocalIntegrals::ALL) {
        computeFaceIntegralsImpl(mesh, faceVectors, dispersionOmega,
                                 boundaryIndicator, integrals);
    }

    virtual double computeL2Error(hpgem::Base::MeshManipulator<dim>& mesh,
                                  std::size_t timeIntegrationVectorId,
                                  InputFunction electricField) = 0;

    virtual FieldT computeField(
        const hpgem::Base::Element* element, const PointReferenceT& p,
        double omega,
        const hpgem::LinearAlgebra::MiddleSizeVector& coefficients) const = 0;
    virtual FieldT computeCurlField(
        const hpgem::Base::Element* element, const PointReferenceT& p,
        double omega,
        const hpgem::LinearAlgebra::MiddleSizeVector& coefficients) const = 0;

    virtual void writeFields(hpgem::Output::VTKSpecificTimeWriter<dim>& writer,
                             double omega,
                             std::size_t timeIntegrationVectorId) const = 0;

    /**
     * Compute the energy flux through a face by integrating normal component of
     * the numerical Poynting vector.
     *
     * Note:
     *  - The direction of the normal is outward from the element on the given
     *    side.
     *  - When a background field is give, 4 fluxes will be computed
     *      0. The flux of the field itself
     *      1. The flux of the background field
     *      2. The flux from cross terms between the field and the background
     *      3. The flux from the combined field (should equal 0+1+2)
     *  - The implementation of the Poynting vector may include extra numerical
     *   terms from the discretization.
     *
     * @param face The face
     * @param side The side which is the inside
     * @param waveNumber The angular wavenumber omega used in the computation
     *    (Poynting vector scales as 1/wavenumber)
     * @param timeIntegrationVectorId id of the vector with coefficients
     * @param background Optional background field
     * @return The fluxes
     */
    virtual LinearAlgebra::SmallVector<4> computeEnergyFluxes(
        Base::Face& face, Base::Side side, double waveNumber,
        std::size_t timeIntegrationVectorId,
        const FieldPattern<dim>* background) {
        return {};
    }

   protected:
    virtual void computeElementIntegralsImpl(
        hpgem::Base::MeshManipulator<dim>& mesh,
        const std::map<std::size_t, InputFunction>& elementVectors,
        double dispersionOmega, LocalIntegrals integrals) = 0;

    virtual void computeFaceIntegralsImpl(
        hpgem::Base::MeshManipulator<dim>& mesh,
        const std::map<std::size_t, FaceInputFunction>& faceVectors,
        double dispersionOmega, BoundaryConditionIndicator boundaryIndicator,
        LocalIntegrals integrals) = 0;
};

}  // namespace DGMax

#endif  // HPGEM_ABSTRACTDISCRETIZATION_H
