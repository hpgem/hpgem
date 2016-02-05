/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, University of Twente
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef PHYSICALFACE_H_
#define PHYSICALFACE_H_

#include "PhysicalElement.h"
#include "FaceMatrix.h"
#include "Face.h"

namespace Base
{

    ///Forward declaration of all coordinate transformations.
    ///\todo check if these can be circumvented
    template<std::size_t DIM>
    class CoordinateTransformation;
    template<std::size_t DIM>
    class H1ConformingTransformation;
    template<std::size_t DIM>
    class HCurlConformingTransformation;
    template<std::size_t DIM>
    class HDivConformingTransformation;
    template<std::size_t DIM>
    class DoNotScaleIntegrands;
    template<std::size_t DIM>
    class IdentityTransformation;


    //class is final as a reminder that there is no virtual destructor
    //note that none of the functions in here is marked const, because a PhysicalFace reserves the right to alter its internal state to optimize future repeated calls
    //note that names in this class match the names in Face unless this makes no sense
    //when you use a physical face in the kernel be careful to avoid infinite recursion
    ///\todo generalize implementation to support the cached data
    template<std::size_t DIM>
    class PhysicalFace final
    {
    public:
        PhysicalFace(bool forInternalFace)
            : left(), right(), transform_((new H1ConformingTransformation<DIM>())), isInternal_(forInternalFace), hasPointReference(false), hasFace(false)  //other fields will be initialized when we have more information
        {
                hasFaceMatrix = false;
                hasFaceVector = false;
                hasLeftRightMatrix = false;
                hasRightLeftMatrix = false;
        }

        PhysicalFace(const PhysicalFace&) = delete;
        PhysicalFace(PhysicalFace&&) = delete;

        ///value of basis function i at the current reference point; indexing functions in the right element after functions in the left element
        double basisFunction(std::size_t i);

        ///value of basis function i at the current reference point; indexing the left and the right element separately
        double basisFunction(Side side, std::size_t i);

        ///derivative of basis function i at the current reference point; indexing functions in the right element after functions in the left element
        const LinearAlgebra::SmallVector<DIM>& basisFunctionDeriv(std::size_t i);

        ///derivative of basis function i at the current reference point; indexing the left and the right element separately
        const LinearAlgebra::SmallVector<DIM>& basisFunctionDeriv(Side side, std::size_t i);

        ///value of basis function i multiplied by the normal vector at the current reference point; indexing functions in the right element after functions in the left element
        /// \details this normal vector includes a scaling by the surface area of the face
        const LinearAlgebra::SmallVector<DIM>& basisFunctionNormal(std::size_t i);

        ///value of basis function i multiplied by the normal vector at the current reference point; indexing the left and the right element separately
        /// \details this normal vector includes a scaling by the surface area of the face
        const LinearAlgebra::SmallVector<DIM>& basisFunctionNormal(Side side, std::size_t i);

        ///value of basis function i multiplied by the normal vector at the current reference point; indexing functions in the right element after functions in the left element
        /// \details does not do any scaling so you probably have to scale the integrand separately
        const LinearAlgebra::SmallVector<DIM>& basisFunctionUnitNormal(std::size_t i);

        ///value of basis function i multiplied by the normal vector at the current reference point; indexing the left and the right element separately
        /// \details does not do any scaling so you probably have to scale the integrand separately
        const LinearAlgebra::SmallVector<DIM>& basisFunctionUnitNormal(Side side, std::size_t i);


        ///value of basis function i at the current reference point; indexing functions in the right element after functions in the left element
        void basisFunction(std::size_t i, LinearAlgebra::SmallVector<DIM>& result);

        ///value of basis function i at the current reference point; indexing the left and the right element separately
        void basisFunction(Side side, std::size_t i, LinearAlgebra::SmallVector<DIM>& result);

        ///curl of basis function i at the current reference point; indexing functions in the right element after functions in the left element
        const LinearAlgebra::SmallVector<DIM>& basisFunctionCurl(std::size_t i);

        ///curl of basis function i at the current reference point; indexing the left and the right element separately
        const LinearAlgebra::SmallVector<DIM>& basisFunctionCurl(Side side, std::size_t i);

        ///divergence of basis function i at the current reference point; indexing functions in the right element after functions in the left element
        const double& basisFunctionDiv(std::size_t i);

        ///divergence of basis function i at the current reference point; indexing the left and the right element separately
        const double& basisFunctionDiv(Side side, std::size_t i);

        ///value of basis function i multiplied by the normal vector at the current reference point; indexing functions in the right element after functions in the left element
        /// \details this normal vector includes a scaling by the surface area of the face
        void basisFunctionNormal(std::size_t i, LinearAlgebra::SmallVector<DIM>& result);

        ///value of basis function i multiplied by the normal vector at the current reference point; indexing the left and the right element separately
        /// \details this normal vector includes a scaling by the surface area of the face
        void basisFunctionNormal(Side side, std::size_t i, LinearAlgebra::SmallVector<DIM>& result);

        ///value of basis function i multiplied by the normal vector at the current reference point; indexing functions in the right element after functions in the left element
        /// \details does not do any scaling so you probably have to scale the integrand separately
        void basisFunctionUnitNormal(std::size_t i, LinearAlgebra::SmallVector<DIM>& result);

        ///value of basis function i multiplied by the normal vector at the current reference point; indexing the left and the right element separately
        /// \details does not do any scaling so you probably have to scale the integrand separately
        void basisFunctionUnitNormal(Side side, std::size_t i, LinearAlgebra::SmallVector<DIM>& result);


        ///the trace of the solution
        const LinearAlgebra::MiddleSizeVector& getSolution(Side side);

        ///the trace of the solution
        void getSolution(Side side, std::vector<LinearAlgebra::SmallVector<DIM> >& result);

        ///derivative of the trace of the solution
        const std::vector<LinearAlgebra::SmallVector<DIM> >& getSolutionDeriv(Side side);

        ///curl of the trace of the solution
        const std::vector<LinearAlgebra::SmallVector<DIM> >& getSolutionCurl(Side side);

        ///divergence of the trace of the solution
        const LinearAlgebra::MiddleSizeVector& getSolutionDiv(Side side);

        ///trace of the solution multiplied by the normal vector
        std::vector<LinearAlgebra::SmallVector<DIM> > getSolutionNormal(Side side);

        ///trace of the solution multiplied by the normal vector
        std::vector<LinearAlgebra::SmallVector<DIM> > getSolutionUnitNormal(Side side);

        ///trace of the solution multiplied by the normal vector
        void getSolutionNormal(Side side, LinearAlgebra::SmallVector<DIM>& result);

        ///trace of the solution multiplied by the normal vector
        void getSolutionUnitNormal(Side side, LinearAlgebra::SmallVector<DIM>& result);

        ///the current reference point
        const Geometry::PointReference<DIM - 1>& getPointReference();

        ///the current physical point
        const Geometry::PointPhysical<DIM>& getPointPhysical();

        ///the normal vector
        const LinearAlgebra::SmallVector<DIM>& getNormalVector();

        ///the normal vector, scaled to have unit length
        const LinearAlgebra::SmallVector<DIM>& getUnitNormalVector();

        ///the length of the normal vector
        double getRelativeSurfaceArea();

        ///\deprecated assumes the area of the reference geometry is 1
        double getSurfaceArea()
        {
            return getRelativeSurfaceArea();
        }

        ///note that this matrix and the side based result matrix have no implied coupling
        FaceMatrix& getResultMatrix();

        ///note that these 4 matrices and the full result matrix have no implied coupling
        LinearAlgebra::MiddleSizeMatrix& getResultMatrix(Side iSide, Side jSide);

        ///vector is initialized with space to store all expansion coefficients and has no implied coupling with the short result vectors
        LinearAlgebra::MiddleSizeVector& getResultVector();

        ///vector is initialized with space to store only expansion coefficients for one side and has no implied coupling with the long result vector
        LinearAlgebra::MiddleSizeVector& getResultVector(Side side);

        ///check if this PhysicalFace is an internal face or a face on a periodic boundary.
        bool isInternal();
        
        ///\deprecated Does not conform naming conventions, use getNumberOfBasisFunctions instead
        std::size_t getNumOfBasisFunctions()
        {
            return getNumberOfBasisFunctions();
        }

        ///get the total number of basis functions that might be nonzero on the face
        std::size_t getNumberOfBasisFunctions()
        {
            return face_->getNumberOfBasisFunctions();
        }
        
        ///\deprecated Does not conform naming conventions, use getNumberOfUnknowns instead
        std::size_t getNumOfUnknowns()
        {
            return getNumberOfUnknowns();
        }

        /// get the number of variables present in the problem
        std::size_t getNumberOfUnknowns()
        {
            return face_->getPtrElementLeft()->getNumberOfUnknowns();
        }

        ///convert a function index and a variable index to a single index in the contiguous numbering of the face
        std::size_t convertToSingleIndex(Base::Side side, std::size_t functionId, std::size_t variableId)
        {
            return face_->convertToSingleIndex(side, functionId, variableId);
        }

        ///provides access to the left and right physical elements in case you need it
        PhysicalElement<DIM>& getPhysicalElement(Base::Side side)
        {
            if(side == Base::Side::LEFT)
            {
                return left;
            }
            else
            {
                logger.assert(isInternal(), "This physical face is meant for boundaries and can only see left elements");
                return right;
            }
        }

        ///provides access to the element on the opposite side of iSide.
        PhysicalElement<DIM>& getPhysicalElementOpposite(Base::Side side)
        {
            if(side == Base::Side::LEFT)
            {
                logger.assert(isInternal(), "This physical face is meant for boundaries and can only see left elements");
                return right;
            }
            else
            {
                return left;
            }
        }

        ///get the index of the face
        std::size_t getID()
        {
            return face_->getID();
        }

        ///direct access to the underlying face in case it is needed
        const Face* getFace();

        ///getTransform should only be needed internally
        const CoordinateTransformation<DIM>* getTransform();

        void setPointReference(const Geometry::PointReference<DIM - 1>& point);
        void setFace(const Face* face);
        void setTransform(std::shared_ptr<CoordinateTransformation<DIM> >& transform);

        template <class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ///Boost only allows to serialize types it knows about. Usually this is not a problem, but we have never written a subclass of
            ///CoordinateTransformation before. It is also possible to register new types without writing them. If you want to use a new
            ///coordinate transformation from the kernel with PhysicalElement you should also register it here. If you have an application-
            ///specific coordinate transformation in your application, you should probably register it over there instead.
            ///\todo think of a better solution
            //ar.template register_type<H1ConformingTransformation<DIM>>();
            ar.template register_type<HCurlConformingTransformation<DIM>>();
            ar.template register_type<HDivConformingTransformation<DIM>>();
            ar.template register_type<DoNotScaleIntegrands<DIM>>();
            ar.template register_type<IdentityTransformation<DIM>>();

            ar & transform_;
            ///\todo check if other data should be saved as well
        }

        void setQuadratureRule(QuadratureRules::GaussQuadratureRule *rule);
        void setQuadraturePointIndex(std::size_t index);

    private:
        PhysicalElement<DIM> left, right;
        std::size_t nLeftBasisFunctions;

        std::vector<LinearAlgebra::SmallVector<DIM> > basisFunctionNormal_;
        std::vector<LinearAlgebra::SmallVector<DIM> > vectorBasisFunctionNormal_;
        std::vector<LinearAlgebra::SmallVector<DIM> > basisFunctionUnitNormal_;
        std::vector<LinearAlgebra::SmallVector<DIM> > vectorBasisFunctionUnitNormal_;
        std::vector<LinearAlgebra::SmallVector<DIM> > solutionNormal_;
        std::vector<LinearAlgebra::SmallVector<DIM> > vectorSolutionNormal_;
        std::vector<LinearAlgebra::SmallVector<DIM> > solutionUnitNormal_;
        std::vector<LinearAlgebra::SmallVector<DIM> > vectorSolutionUnitNormal_;

        Geometry::PointReference<DIM - 1> pointReference_;
        QuadratureRules::GaussQuadratureRule* quadratureRule_;
        const Face* face_;
        std::shared_ptr<CoordinateTransformation<DIM> > transform_;
        LinearAlgebra::SmallVector<DIM> normal;
        LinearAlgebra::SmallVector<DIM> unitNormal;
        double normalNorm;

        //need to store this to keep it existing
        std::shared_ptr<const Geometry::MappingReferenceToReference<1>> mapToLeftElement, mapToRightElement;

        FaceMatrix resultMatrix;
        LinearAlgebra::MiddleSizeMatrix leftRightMatrix, rightLeftMatrix;
        LinearAlgebra::MiddleSizeVector resultVector;

        bool isInternal_, hasPointReference, hasFace;

        bool hasFaceMatrix, hasFaceVector, hasLeftRightMatrix, hasRightLeftMatrix;
        bool hasBasisFunctionNormal, hasBasisFunctionUnitNormal, hasVectorBasisFunctionNormal, hasVectorBasisFunctionUnitNormal;
        bool hasSolutionNormal, hasSolutionUnitNormal, hasVectorSolutionNormal, hasVectorSolutionUnitNormal;
        bool hasNormal, hasUnitNormal, hasNormalNorm;
    };
}

#include "PhysicalFace_Impl.h"
    
#endif /* PHYSICALFACE_H_ */
