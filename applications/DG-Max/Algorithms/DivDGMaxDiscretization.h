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

#ifndef HPGEM_APP_DIVDGMAXDISCRETIZATION_H
#define HPGEM_APP_DIVDGMAXDISCRETIZATION_H

#include <functional>
#include <memory>

#include "AbstractDiscretization.h"

#include "Base/HCurlConformingTransformation.h"
#include "Base/H1ConformingTransformation.h"
#include "Base/MeshManipulator.h"
#include "Integration/ElementIntegral.h"
#include "Integration/FaceIntegral.h"
#include "LinearAlgebra/MiddleSizeMatrix.h"
#include "LinearAlgebra/MiddleSizeVector.h"
#include "LinearAlgebra/SmallVector.h"
#include "Output/VTKSpecificTimeWriter.h"

#include "ProblemTypes/BoundaryConditionType.h"
#include "Material.h"

// Forward definitions
namespace hpgem {
namespace Base {
template <std::size_t>
class PhysicalElement;
template <std::size_t DIM>
class PhysicalFace;
}  // namespace Base

namespace Utilities {
class ElementLocalIndexing;
class FaceLocalIndexing;
}  // namespace Utilities

}  // namespace hpgem

/// Base class with all dimensionless constants
class DivDGMaxDiscretizationBase : public DGMax::AbstractDiscretizationBase {
   public:
    enum class FluxType { IP, BREZZI };

    struct Stab {
        /// \brief Stabilization parameter for the tangential part of u/v.
        double stab1;
        /// \brief Stabilization parameter for the normal part of u/v.
        double stab2;
        /// \brief Stabilization parameter for p/q.
        double stab3;
        /// Numerical flux and stabilization to use for the tangential part of
        /// u/v (vector part)
        FluxType fluxType1 = FluxType::IP;
        /// Numerical flux and stabilization to use for the normal part of u/v
        /// (vector part)
        FluxType fluxType2 = FluxType::IP;
        /// Numerical flux and stabilization to use for p/q (scalar part)
        FluxType fluxType3 = FluxType::IP;

        bool hasFlux(FluxType flux) {
            return fluxType1 == flux || fluxType2 == flux || fluxType3 == flux;
        }

        void setAllFluxeTypes(FluxType type) {
            fluxType1 = type;
            fluxType2 = type;
            fluxType3 = type;
        }
    };
};

using namespace hpgem;
/// \brief Discontinuous Galerkin discretization for Maxwell, where the
/// divergence constraint (div E = 0) is part of the discretization.
///
/// This implementation is based on chapter 5 of Devashish2017PhD, and similar
/// to for example Lu2016JSciComput. It decomposes E = u + grad p, forming a
/// mixed system (for eigenvalue problems) of the form
/// [ A   B ] [u] = λ [ M 0 ] [u]
/// [ BT -C ] [p]     [ 0 0 ] [p]
///
/// Where A corresponds to the curl-curl operator, B the coupling between
/// u and p, and C is a stabilization term. The matrix M is the mass
/// matrix, corresponding to the omega^2 E term in the timeharmonic
/// formulation.
template <std::size_t DIM>
class DivDGMaxDiscretization : public DGMax::AbstractDiscretization<DIM>,
                               public virtual DivDGMaxDiscretizationBase {
   public:
    /// Value class for the solution.
    struct Fields {

        Fields()
            : electricField(), electricFieldCurl(), potential(0), material(){};

        // The electric field
        LinearAlgebra::SmallVectorC<DIM> electricField;
        // Curl of the electric field, corresponds to the H field by
        // Curl E = i omega mu_r H (omega = the wave number)
        LinearAlgebra::SmallVectorC<DIM> electricFieldCurl;
        // Complex valued p scalar function;
        std::complex<double> potential;
        // Material constants
        DGMax::Material material;
    };

    // See notes in DGMaxDiscretization

    using typename DGMax::AbstractDiscretization<DIM>::PointPhysicalT;
    using typename DGMax::AbstractDiscretization<DIM>::PointReferenceT;
    using typename DGMax::AbstractDiscretization<DIM>::InputFunction;
    using typename DGMax::AbstractDiscretization<DIM>::FaceInputFunction;

    DivDGMaxDiscretization(std::size_t order, Stab stab);

    size_t getOrder() const override { return order_; }
    size_t getNumberOfUnknowns() const override { return 2; }
    size_t getNumberOfElementMatrices() const override { return 2; }
    size_t getNumberOfFaceMatrices() const override { return 2; }

    void initializeBasisFunctions(Base::MeshManipulator<DIM>& mesh) const final;

    // TODO: LJ include the same norms as in DGMaxDiscretization
    double computeL2Error(Base::MeshManipulator<DIM>& mesh,
                          std::size_t timeIntegrationVectorId,
                          InputFunction electricField) final;

    Fields computeFields(
        const Base::Element* element,
        const Geometry::PointReference<DIM>& point,
        const LinearAlgebra::MiddleSizeVector& coefficients) const;

    LinearAlgebra::SmallVectorC<DIM> computeField(
        const Base::Element* element,
        const Geometry::PointReference<DIM>& point,
        const LinearAlgebra::MiddleSizeVector& coefficients) const final;

    LinearAlgebra::SmallVectorC<DIM> computeCurlField(
        const hpgem::Base::Element* element, const PointReferenceT& p,
        const hpgem::LinearAlgebra::MiddleSizeVector& coefficients) const final;

    std::complex<double> computePotential(
        const Base::Element* element,
        const Geometry::PointReference<DIM>& point,
        const LinearAlgebra::MiddleSizeVector& coefficients) const;

    void writeFields(Output::VTKSpecificTimeWriter<DIM>& output,
                     std::size_t timeIntegrationVectorId) const final;

    double computeEnergyFlux(Base::Face& face, hpgem::Base::Side side,
                             double wavenumber,
                             std::size_t timeIntegrationVectorId) final;

   private:
    void computeElementIntegralsImpl(
        Base::MeshManipulator<DIM>& mesh,
        const std::map<std::size_t, InputFunction>& elementVectors,
        LocalIntegrals integrals) final;

    void computeFaceIntegralsImpl(
        Base::MeshManipulator<DIM>& mesh,
        const std::map<std::size_t, FaceInputFunction>& boundaryVectors,
        DGMax::BoundaryConditionIndicator indicator,
        LocalIntegrals integrals) final;

    /// Compute Mass and Stiffness matrix for the element
    void computeElementMatrices(Base::Element* element,
                                Utilities::ElementLocalIndexing& indexing);

    /// Element part of matrix B and B^T, with zeros around it (- grad p, eps v)
    void elementScalarVectorCoupling(
        Base::PhysicalElement<DIM>& el,
        LinearAlgebra::MiddleSizeMatrix& ret) const;

    /// The source vector for harmonic problems.
    void elementSourceVector(Base::PhysicalElement<DIM>& el,
                             const InputFunction& source,
                             LinearAlgebra::MiddleSizeVector& ret) const;

    /// Part 1 of face integrand for the stiffness matrix. Consists of the
    /// terms:
    ///
    ///  1.  -[[v]]_t {{mu^{-1} curl u}} - [[u]]_t {{mu^{-1} curl v}}
    ///  2. For IP-stab1: stab1/diameter * [[u]]_t [[v]]_t
    ///  3. For IP-stab2: stab2*diameter/espMax [[eps u]]_n . [[eps v]]_n
    void faceStiffnessMatrixFieldIntegrand(
        Base::PhysicalFace<DIM>& fa,
        const Utilities::FaceLocalIndexing& indexing,
        LinearAlgebra::MiddleSizeMatrix& ret) const;
    /// Part 2 of the face integrand for the stiffness matrix, consisting of the
    /// terms with the potential. This contributes two parts
    ///
    ///  1. The coupling scalar-vector: [[p]] {{eps v}}
    ///  2. For IP-stab3: stab3/diameter * epsMax [[p]] [[q]]
    void addFaceMatrixPotentialIntegrand(
        Base::PhysicalFace<DIM>& fa,
        const Utilities::FaceLocalIndexing& indexing,
        DGMax::BoundaryConditionType bct,
        LinearAlgebra::MiddleSizeMatrix& ret) const;

    LinearAlgebra::MiddleSizeMatrix brezziFluxBilinearTerm(
        Base::Face* face, DGMax::BoundaryConditionType bct);

    LinearAlgebra::MiddleSizeMatrix computeFaceImpedanceIntegrand(
        Base::PhysicalFace<DIM>& face, Utilities::FaceLocalIndexing& indexing,
        DGMax::BoundaryConditionType& bct) const;

    /// \brief Compute mass matrix for vector components on elements adjacent to
    /// a face
    ///
    /// Computes the mass matrix corresponding to the vector valued basis
    /// functions for the element(s) adjacent to the face. For an internal face
    /// the result is a block diagonal matrix with the two mass matrices for the
    /// separate elements.
    ///
    /// Note that the resulting matrix only contains the degrees of freedom for
    /// the vector valued basis functions.
    /// \param face The face to compute the local mass matrix for
    /// \return The mass matrix
    LinearAlgebra::MiddleSizeMatrix computeFaceVectorMassMatrix(
        Base::Face* face);

    /// \brief Compute mass matrix for scalar basis functions on elements
    /// adjacent to a face
    ///
    /// Same as computeFaceVectorMassMatrix but for the scalar basis functions.
    /// \param face The face to compute the matrix forr
    /// \return The mass matrix
    LinearAlgebra::MiddleSizeMatrix computeFaceScalarMassMatrix(
        Base::Face* face);

    /// \brief Compute projection matrix of the jump of the scalar basis
    /// functions
    ///     for the lift operator.
    ///
    /// Compute the projection matrix
    ///   R_{ij} = integral_F p_j n.{{u_i}} dS
    /// with p_j the scalar basis functions, u_i the vector valued ones and n
    /// the outward pointing normal to the element on which p_j has support.
    ///
    /// This term is needed for the implementation of the lifting operator for
    /// the Brezzi fluxes.
    /// \param face The face to compute it on
    /// \return The (dofs u) by (dofs p) projection matrix.
    LinearAlgebra::MiddleSizeMatrix computeScalarLiftProjector(
        Base::Face* face);

    /// \brief Compute the projection matrix of the tangential jump of the
    /// vector
    ///   basis functions for the lifting operator.
    ///
    /// Compute the projection matrix
    ///   R_{ij} = integral_F (n x u_j) {{u_i}} dS
    /// with u_i, u_j the vector basis functions and n the outward pointing
    /// normal
    ///   to the element on which u_j has support.
    /// \param face The face to compute it on
    /// \return The (dofs u)^2 projection matrix.
    LinearAlgebra::MiddleSizeMatrix computeVectorLiftProjector(
        Base::Face* face);

    /// \brief Compute the projection matrix for the normal part of the
    ///         vector basis functions in the lifting operators.
    ///
    /// Compute the projection matrix
    ///   S_{ij} = integral_F [[epsilon u_j]]_n {{p_i}} dS
    /// with u_j and p_i the basis functions for the vector part and scalar part
    /// and epsilon is the permittivity.
    ///
    /// \param face The face to compute it on
    /// \return The (dofs p) by (dofs u) projection matrix.
    LinearAlgebra::MiddleSizeMatrix computeVectorNormalLiftProjector(
        Base::Face* face);

    void faceBoundaryVector(Base::PhysicalFace<DIM>& fa,
                            const FaceInputFunction& boundaryValue,
                            LinearAlgebra::MiddleSizeVector& ret,
                            DGMax::BoundaryConditionType bct) const;

    /// Compute contribution of the brezzi flux to the face vector on the
    /// boundary
    LinearAlgebra::MiddleSizeVector brezziFluxBoundaryVector(
        Base::Face* face, const FaceInputFunction& boundaryValue);

    double elementErrorIntegrand(Base::PhysicalElement<DIM>& el,
                                 std::size_t timeVector,
                                 const InputFunction& exactValues) const;

    std::size_t order_;
    Stab stab_;

    /// Shared parts for computing integrals
    Integration::ElementIntegral<DIM> elementIntegrator_;
    Integration::FaceIntegral<DIM> faceIntegrator_;
    std::shared_ptr<Base::HCurlConformingTransformation<DIM>> fieldTransform_;
    std::shared_ptr<Base::H1ConformingTransformation<DIM>> potentialTransform_;
};

// TODO: Deduction fails for a templated variant, hence using explicit versions
// here

std::ostream& operator<<(std::ostream& os,
                         const DivDGMaxDiscretizationBase::Stab& stab);

#endif  // HPGEM_APP_DIVDGMAXDISCRETIZATION_H
