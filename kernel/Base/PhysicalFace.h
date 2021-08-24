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

#ifndef HPGEM_KERNEL_PHYSICALFACE_H
#define HPGEM_KERNEL_PHYSICALFACE_H

#include "PhysicalElement.h"

#include "FaceMatrix.h"
#include "Face.h"
#include "LazyCached.h"

namespace hpgem {
namespace Base {

template <std::size_t DIM>
class CoordinateTransformation;
template <std::size_t DIM>
class H1ConformingTransformation;

// class is final as a reminder that there is no virtual destructor
// note that none of the functions in here is marked const, because a
// PhysicalFace reserves the right to alter its internal state to optimize
// future repeated calls note that names in this class match the names in Face
// unless this makes no sense when you use a physical face in the kernel be
// careful to avoid infinite recursion
///\todo generalize implementation to support the cached data
template <std::size_t DIM>
class PhysicalFace final {
   public:
    PhysicalFace(bool forInternalFace)
        : left(),
          right(),
          isInternal_(forInternalFace),
          hasPointReference(false),
          hasFace(false),
          basisFunctionDeriv_(this,
                              &PhysicalFace<DIM>::computeBasisFunctionDeriv),
          basisFunctionCurl_(this,
                             &PhysicalFace<DIM>::computeBasisFunctionCurl),
          basisFunctionDiv_(this, &PhysicalFace<DIM>::computeBasisFunctionDiv)
    // other fields will be initialized when we have more
    // information
    {
        std::shared_ptr<Base::CoordinateTransformation<DIM>> transform(
            new H1ConformingTransformation<DIM>{});
        transform_.push_back(transform);
        hasFaceMatrix = false;
        hasFaceVector = false;
        hasLeftRightMatrix = false;
        hasRightLeftMatrix = false;
    }

    PhysicalFace(const PhysicalFace&) = delete;
    PhysicalFace(PhysicalFace&&) = delete;

    /// value of basis function i at the current reference point; indexing
    /// functions in the right element after functions in the left element
    double basisFunction(std::size_t i);
    double basisFunction(std::size_t i, std::size_t unknown);

    /// value of basis function i at the current reference point; indexing the
    /// left and the right element separately
    double basisFunction(Side side, std::size_t i);
    double basisFunction(Side side, std::size_t i, std::size_t unknown);

    /// derivative of basis function i at the current reference point; indexing
    /// functions in the right element after functions in the left element
    const LinearAlgebra::SmallVector<DIM>& basisFunctionDeriv(
        std::size_t i, std::size_t unknown = 0);

    /// derivative of basis function i at the current reference point; indexing
    /// the left and the right element separately
    const LinearAlgebra::SmallVector<DIM>& basisFunctionDeriv(
        Side side, std::size_t i, std::size_t unknown = 0);

    /// value of basis function i multiplied by the normal vector at the current
    /// reference point; indexing functions in the right element after functions
    /// in the left element
    /// \details this normal vector includes a scaling by the surface area of
    /// the face
    const LinearAlgebra::SmallVector<DIM>& basisFunctionNormal(std::size_t i);
    const LinearAlgebra::SmallVector<DIM>& basisFunctionNormal(
        std::size_t i, std::size_t unknown);

    /// value of basis function i multiplied by the normal vector at the current
    /// reference point; indexing the left and the right element separately
    /// \details this normal vector includes a scaling by the surface area of
    /// the face
    const LinearAlgebra::SmallVector<DIM>& basisFunctionNormal(Side side,
                                                               std::size_t i);
    const LinearAlgebra::SmallVector<DIM>& basisFunctionNormal(
        Side side, std::size_t i, std::size_t unknown);

    /// value of basis function i multiplied by the normal vector at the current
    /// reference point; indexing functions in the right element after functions
    /// in the left element
    /// \details does not do any scaling so you probably have to scale the
    /// integrand separately
    const LinearAlgebra::SmallVector<DIM>& basisFunctionUnitNormal(
        std::size_t i);
    const LinearAlgebra::SmallVector<DIM>& basisFunctionUnitNormal(
        std::size_t i, std::size_t unknown);

    /// value of basis function i multiplied by the normal vector at the current
    /// reference point; indexing the left and the right element separately
    /// \details does not do any scaling so you probably have to scale the
    /// integrand separately
    const LinearAlgebra::SmallVector<DIM>& basisFunctionUnitNormal(
        Side side, std::size_t i);
    const LinearAlgebra::SmallVector<DIM>& basisFunctionUnitNormal(
        Side side, std::size_t i, std::size_t unknown);

    /// value of basis function i at the current reference point; indexing
    /// functions in the right element after functions in the left element
    void basisFunction(std::size_t i, LinearAlgebra::SmallVector<DIM>& result);
    void basisFunction(std::size_t i, LinearAlgebra::SmallVector<DIM>& result,
                       std::size_t unknown);

    /// value of basis function i at the current reference point; indexing the
    /// left and the right element separately
    void basisFunction(Side side, std::size_t i,
                       LinearAlgebra::SmallVector<DIM>& result);
    void basisFunction(Side side, std::size_t i,
                       LinearAlgebra::SmallVector<DIM>& result,
                       std::size_t unknown);

    /// curl of basis function i at the current reference point; indexing
    /// functions in the right element after functions in the left element
    const LinearAlgebra::SmallVector<DIM>& basisFunctionCurl(
        std::size_t i, std::size_t unknown = 0);

    /// curl of basis function i at the current reference point; indexing the
    /// left and the right element separately
    const LinearAlgebra::SmallVector<DIM>& basisFunctionCurl(
        Side side, std::size_t i, std::size_t unknown = 0);

    /// divergence of basis function i at the current reference point; indexing
    /// functions in the right element after functions in the left element
    const double& basisFunctionDiv(std::size_t i, std::size_t unknown = 0);

    /// divergence of basis function i at the current reference point; indexing
    /// the left and the right element separately
    const double& basisFunctionDiv(Side side, std::size_t i,
                                   std::size_t unknown = 0);

    /// value of basis function i multiplied by the normal vector at the current
    /// reference point; indexing functions in the right element after functions
    /// in the left element
    /// \details this normal vector includes a scaling by the surface area of
    /// the face
    void basisFunctionNormalCross(std::size_t i,
                                  LinearAlgebra::SmallVector<DIM>& result);
    void basisFunctionNormalCross(std::size_t i,
                                  LinearAlgebra::SmallVector<DIM>& result,
                                  std::size_t unknown);

    /// value of basis function i multiplied by the normal vector at the current
    /// reference point; indexing the left and the right element separately
    /// \details this normal vector includes a scaling by the surface area of
    /// the face
    void basisFunctionNormalCross(Side side, std::size_t i,
                                  LinearAlgebra::SmallVector<DIM>& result);
    void basisFunctionNormalCross(Side side, std::size_t i,
                                  LinearAlgebra::SmallVector<DIM>& result,
                                  std::size_t unknown);

    /// value of basis function i multiplied by the normal vector at the current
    /// reference point; indexing functions in the right element after functions
    /// in the left element
    /// \details does not do any scaling so you probably have to scale the
    /// integrand separately
    void basisFunctionUnitNormalCross(std::size_t i,
                                      LinearAlgebra::SmallVector<DIM>& result);
    void basisFunctionUnitNormalCross(std::size_t i,
                                      LinearAlgebra::SmallVector<DIM>& result,
                                      std::size_t unknown);

    /// value of basis function i multiplied by the normal vector at the current
    /// reference point; indexing the left and the right element separately
    /// \details does not do any scaling so you probably have to scale the
    /// integrand separately
    void basisFunctionUnitNormalCross(Side side, std::size_t i,
                                      LinearAlgebra::SmallVector<DIM>& result);
    void basisFunctionUnitNormalCross(Side side, std::size_t i,
                                      LinearAlgebra::SmallVector<DIM>& result,
                                      std::size_t unknown);

    /// the trace of the solution
    const LinearAlgebra::MiddleSizeVector& getSolution(Side side);

    /// the trace of the solution
    void getSolution(Side side,
                     std::vector<LinearAlgebra::SmallVector<DIM>>& result);

    /// derivative of the trace of the solution
    const std::vector<LinearAlgebra::SmallVector<DIM>>& getSolutionDeriv(
        Side side);

    /// curl of the trace of the solution
    const std::vector<LinearAlgebra::SmallVector<DIM>>& getSolutionCurl(
        Side side);

    /// divergence of the trace of the solution
    const LinearAlgebra::MiddleSizeVector& getSolutionDiv(Side side);

    /// trace of the solution multiplied by the normal vector
    std::vector<LinearAlgebra::SmallVector<DIM>> getSolutionNormal(Side side);

    /// trace of the solution multiplied by the normal vector
    std::vector<LinearAlgebra::SmallVector<DIM>> getSolutionUnitNormal(
        Side side);

    /// trace of the solution multiplied by the normal vector
    void getSolutionNormal(Side side, LinearAlgebra::SmallVector<DIM>& result);

    /// trace of the solution multiplied by the normal vector
    void getSolutionUnitNormal(Side side,
                               LinearAlgebra::SmallVector<DIM>& result);

    /// the current reference point
    const Geometry::PointReference<DIM - 1>& getPointReference();

    /// the current physical point
    const Geometry::PointPhysical<DIM>& getPointPhysical();

    /// the normal vector
    const LinearAlgebra::SmallVector<DIM>& getNormalVector();

    /// the normal vector, scaled to have unit length
    const LinearAlgebra::SmallVector<DIM>& getUnitNormalVector();

    /// the length of the normal vector
    double getRelativeSurfaceArea();

    ///\deprecated assumes the area of the reference geometry is 1
    double getSurfaceArea() { return getRelativeSurfaceArea(); }

    /// note that this matrix and the side based result matrix have no implied
    /// coupling
    FaceMatrix& getResultMatrix();

    /// note that these 4 matrices and the full result matrix have no implied
    /// coupling
    LinearAlgebra::MiddleSizeMatrix& getResultMatrix(Side iSide, Side jSide);

    /// vector is initialized with space to store all expansion coefficients and
    /// has no implied coupling with the short result vectors
    LinearAlgebra::MiddleSizeVector& getResultVector();

    /// vector is initialized with space to store only expansion coefficients
    /// for one side and has no implied coupling with the long result vector
    LinearAlgebra::MiddleSizeVector& getResultVector(Side side);

    /// check if this PhysicalFace is an internal face or a face on a periodic
    /// boundary.
    bool isInternal();

    ///\deprecated Does not conform naming conventions, use
    /// getNumberOfBasisFunctions instead
    std::size_t getNumOfBasisFunctions() { return getNumberOfBasisFunctions(); }

    /// get the total number of basis functions that might be nonzero on the
    /// face
    std::size_t getNumberOfBasisFunctions() {
        return face_->getNumberOfBasisFunctions();
    }

    std::size_t getNumberOfBasisFunctions(std::size_t unknown) {
        return face_->getNumberOfBasisFunctions(unknown);
    }

    ///\deprecated Does not conform naming conventions, use getNumberOfUnknowns
    /// instead
    std::size_t getNumOfUnknowns() { return getNumberOfUnknowns(); }

    /// get the number of variables present in the problem
    std::size_t getNumberOfUnknowns() {
        return face_->getPtrElementLeft()->getNumberOfUnknowns();
    }

    /// convert a function index and a variable index to a single index in the
    /// contiguous numbering of the face
    std::size_t convertToSingleIndex(Base::Side side, std::size_t functionId,
                                     std::size_t variableId) {
        return face_->convertToSingleIndex(side, functionId, variableId);
    }

    /// provides access to the left and right physical elements in case you need
    /// it
    PhysicalElement<DIM>& getPhysicalElement(Base::Side side) {
        if (side == Base::Side::LEFT) {
            return left;
        }
        logger.assert_debug(isInternal(),
                            "This physical face is meant for boundaries "
                            "and can only see left elements");
        return right;
    }

    /// provides access to the element on the opposite side of iSide.
    PhysicalElement<DIM>& getPhysicalElementOpposite(Base::Side side) {
        if (side == Base::Side::LEFT) {
            logger.assert_debug(isInternal(),
                                "This physical face is meant for boundaries "
                                "and can only see left elements");
            return right;
        }
        return left;
    }

    /// get the index of the face
    std::size_t getID() { return face_->getID(); }

    /// direct access to the underlying face in case it is needed
    const Face* getFace();

    /// getTransform should only be needed internally
    const CoordinateTransformation<DIM>* getTransform();
    const CoordinateTransformation<DIM>* getTransform(std::size_t unknown);

    void setPointReference(const Geometry::PointReference<DIM - 1>& point);
    void setFace(const Face* face);
    void setTransform(std::shared_ptr<CoordinateTransformation<DIM>>& transform,
                      std::size_t unknown = 0);

    void setQuadratureRule(QuadratureRules::GaussQuadratureRule* rule);
    void setQuadraturePointIndex(std::size_t index);

   private:
    void resetLazyCaches(std::size_t currentUnknowns);
    void computeBasisFunctionDeriv(
        std::vector<LinearAlgebra::SmallVector<DIM>>& values,
        std::size_t unknown);
    void computeBasisFunctionCurl(
        std::vector<LinearAlgebra::SmallVector<DIM>>& values,
        std::size_t unknown);
    void computeBasisFunctionDiv(std::vector<double>& values,
                                 std::size_t unknown);

    /// Updates the transformation information between the left and right faces
    void updateLeftRightTransform();
    /// Compute a set of direction vectors on the left/right face from
    ///
    /// helper for #updateLeftRightTransform()
    static LinearAlgebra::SmallMatrix<DIM, DIM> computeDirectionVectors(
        const std::array<LinearAlgebra::SmallVector<DIM>, DIM>& points);

    PhysicalElement<DIM> left, right;
    std::vector<std::size_t> nLeftBasisFunctions;

    std::vector<std::vector<LinearAlgebra::SmallVector<DIM>>>
        basisFunctionNormal_;
    std::vector<std::vector<LinearAlgebra::SmallVector<DIM>>>
        vectorBasisFunctionNormal_;
    std::vector<std::vector<LinearAlgebra::SmallVector<DIM>>>
        basisFunctionUnitNormal_;
    std::vector<std::vector<LinearAlgebra::SmallVector<DIM>>>
        vectorBasisFunctionUnitNormal_;
    LazyVectorCached<std::vector<LinearAlgebra::SmallVector<DIM>>>
        basisFunctionDeriv_;
    LazyVectorCached<std::vector<LinearAlgebra::SmallVector<DIM>>>
        basisFunctionCurl_;
    LazyVectorCached<std::vector<double>> basisFunctionDiv_;

    std::vector<LinearAlgebra::SmallVector<DIM>> solutionNormal_;
    std::vector<LinearAlgebra::SmallVector<DIM>> vectorSolutionNormal_;
    std::vector<LinearAlgebra::SmallVector<DIM>> solutionUnitNormal_;
    std::vector<LinearAlgebra::SmallVector<DIM>> vectorSolutionUnitNormal_;

    Geometry::PointReference<DIM - 1> pointReference_;
    QuadratureRules::GaussQuadratureRule* quadratureRule_;
    const Face* face_;
    std::vector<std::shared_ptr<CoordinateTransformation<DIM>>> transform_;
    LinearAlgebra::SmallVector<DIM> normal;
    LinearAlgebra::SmallVector<DIM> unitNormal;
    double normalNorm;

    /**
     * Transformation matrix connecting the left and right faces (up to
     * translation). Given a reference point xref on the face we get two
     * coordinates x_left(xref) and x_right(xref). These two coordinates are
     * related by the affine transformation:
     * x_left = R x_right + x_off
     * For most internal faces the coordinates on the left and right are the
     * same and so R = I and x_off = 0. For periodically connected faces we will
     * have x_off != 0 and/or R != I.
     *
     * In case R != I we have that the faces are rotated[1] with respect to
     * each other. This poses a problem for all vector based quantities. With R
     * != I the left and right faces are not aligned. Hence a vector that is
     * normal with respect to the left face, may be tangent to the right face.
     *
     * For computing face based quantities (like face integrals) we need to
     * choose a fixed coordinate system for the vectors. The choice here is to
     * take that of the left element. All quantities from the right element are
     * transformed such that it looks like the right element was directly
     * stitched to the face.
     *
     * [1] In general R can be a general transformation matrix. However, as both
     * faces have the same shape and area we are restricted to a rotation
     * matrix. Moreover, we assume that the rotation is a proper rotation.
     */
    ValueCoordinateTransformationData<DIM> rightVectorTransform;
    /**
     * Whether a transformation of vectors is needed
     */
    bool requiresTransformation;

    // need to store this to keep it existing
    std::shared_ptr<const Geometry::MappingReferenceToReference<1>>
        mapToLeftElement, mapToRightElement;

    FaceMatrix resultMatrix;
    LinearAlgebra::MiddleSizeMatrix leftRightMatrix, rightLeftMatrix;
    LinearAlgebra::MiddleSizeVector resultVector;

    bool isInternal_, hasPointReference, hasFace;

    bool hasFaceMatrix, hasFaceVector, hasLeftRightMatrix, hasRightLeftMatrix;
    std::vector<bool> hasBasisFunctionNormal, hasBasisFunctionUnitNormal,
        hasVectorBasisFunctionNormal, hasVectorBasisFunctionUnitNormal;
    bool hasSolutionNormal, hasSolutionUnitNormal, hasVectorSolutionNormal,
        hasVectorSolutionUnitNormal;
    bool hasNormal, hasUnitNormal, hasNormalNorm;
};
}  // namespace Base
}  // namespace hpgem
#include "PhysicalFace_Impl.h"

#endif  // HPGEM_KERNEL_PHYSICALFACE_H
