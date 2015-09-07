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

    template<std::size_t DIM>
    class CoordinateTransformation;
    template<std::size_t DIM>
    class H1ConformingTransformation;

    //class is final as a reminder that there is no virtual destructor
    //note that none of the functions in here is marked const, because a PhysicalFace reserves the right to alter its internal state to optimize future repeated calls
    //note that names in this class match the names in Face unless this makes no sense
    //when you use a physical face in the kernel be careful to avoid infinite recursion
    ///\TODO generalize implementation to support the cached data
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

        ///get the total number of basis functions that might be nonzero on the face
        std::size_t getNumOfBasisFunctions()
        {
            return face_->getNumberOfBasisFunctions();
        }

        /// get the number of variables present in the problem
        std::size_t getNumOfUnknowns()
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

        const Geometry::PointReference<DIM - 1>* pointReference_;
        const Face* face_;
        std::shared_ptr<CoordinateTransformation<DIM> > transform_;
        LinearAlgebra::SmallVector<DIM> normal;
        LinearAlgebra::SmallVector<DIM> unitNormal;
        double normalNorm;

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
#include "CoordinateTransformation.h"
    
    template<std::size_t DIM>
    inline double Base::PhysicalFace<DIM>::basisFunction(std::size_t i)
    {
        logger.assert(hasPointReference && hasFace, "Need a location to evaluate the data");
        if(i < nLeftBasisFunctions)
        {
            return left.basisFunction(i);
        }
        else
        {
            logger.assert(isInternal_, "basis function index out of bounds");
            return right.basisFunction(i - nLeftBasisFunctions);
        }
    }
    
    template<std::size_t DIM>
    inline double Base::PhysicalFace<DIM>::basisFunction(Side side, std::size_t i)
    {
        logger.assert(hasPointReference && hasFace, "Need a location to evaluate the data");
        if(side == Side::LEFT)
        {
            return left.basisFunction(i);
        }
        else
        {
            logger.assert(isInternal_, "cannot find the right element for a boundary face");
            return right.basisFunction(i);
        }
    }
    
    template<std::size_t DIM>
    inline const LinearAlgebra::SmallVector<DIM>& Base::PhysicalFace<DIM>::basisFunctionDeriv(std::size_t i)
    {
        logger.assert(hasPointReference && hasFace, "Need a location to evaluate the data");
        if(i < nLeftBasisFunctions)
        {
            return left.basisFunctionDeriv(i);
        }
        else
        {
            logger.assert(isInternal_, "basis function index out of bounds");
            return right.basisFunctionDeriv(i - nLeftBasisFunctions);
        }
    }
    
    template<std::size_t DIM>
    inline const LinearAlgebra::SmallVector<DIM>& Base::PhysicalFace<DIM>::basisFunctionDeriv(Side side, std::size_t i)
    {
        logger.assert(hasPointReference && hasFace, "Need a location to evaluate the data");
        if(side == Side::LEFT)
        {
            return left.basisFunctionDeriv(i);
        }
        else
        {
            logger.assert(isInternal_, "cannot find the right element for a boundary face");
            return right.basisFunctionDeriv(i);
        }
    }
    
    template<std::size_t DIM>
    inline const LinearAlgebra::SmallVector<DIM>& Base::PhysicalFace<DIM>::basisFunctionNormal(std::size_t i)
    {
        logger.assert(hasPointReference && hasFace, "Need a location to evaluate the data");
        if(hasBasisFunctionNormal)
        {
            return basisFunctionNormal_[i];
        }
        else
        {
            hasBasisFunctionNormal = true;
            for(std::size_t j = 0; j < face_->getNumberOfBasisFunctions(); ++j)
            {
                basisFunctionNormal_[j] = getNormalVector() * basisFunction(j);
                if(j >= nLeftBasisFunctions)
                {
                    basisFunctionNormal_[j] *= -1.;
                }
            }
            return basisFunctionNormal_[i];
        }
    }
    
    template<std::size_t DIM>
    inline const LinearAlgebra::SmallVector<DIM>& Base::PhysicalFace<DIM>::basisFunctionNormal(Side side, std::size_t i)
    {
        logger.assert(hasPointReference && hasFace, "Need a location to evaluate the data");
        if(side == Side::LEFT)
        {
            return basisFunctionNormal(i);
        }
        else
        {
            logger.assert(isInternal_, "cannot find the right element for a boundary face");
            return basisFunctionNormal(i + nLeftBasisFunctions);
        }
    }
    
    template<std::size_t DIM>
    inline const LinearAlgebra::SmallVector<DIM>& Base::PhysicalFace<DIM>::basisFunctionUnitNormal(std::size_t i)
    {
        logger.assert(hasPointReference && hasFace, "Need a location to evaluate the data");
        if(hasBasisFunctionUnitNormal)
        {
            return basisFunctionUnitNormal_[i];
        }
        else
        {
            hasBasisFunctionUnitNormal = true;
            for(std::size_t j = 0; j < face_->getNumberOfBasisFunctions(); ++j)
            {
                basisFunctionUnitNormal_[j] = getUnitNormalVector() * basisFunction(j);
                if(j >= nLeftBasisFunctions)
                {
                    basisFunctionUnitNormal_[j] *= -1.;
                }
            }
            return basisFunctionUnitNormal_[i];
        }
    }
    
    template<std::size_t DIM>
    inline const LinearAlgebra::SmallVector<DIM>& Base::PhysicalFace<DIM>::basisFunctionUnitNormal(Side side, std::size_t i)
    {
        logger.assert(hasPointReference && hasFace, "Need a location to evaluate the data");
        if(side == Side::LEFT)
        {
            return basisFunctionUnitNormal(i);
        }
        else
        {
            logger.assert(isInternal_, "cannot find the right element for a boundary face");
            return basisFunctionUnitNormal(i + nLeftBasisFunctions);
        }
    }
    
    template<std::size_t DIM>
    inline void Base::PhysicalFace<DIM>::basisFunction(std::size_t i, LinearAlgebra::SmallVector<DIM>& result)
    {
        logger.assert(hasPointReference && hasFace, "Need a location to evaluate the data");
        if(i < nLeftBasisFunctions)
        {
            left.basisFunction(i, result);
            return;
        }
        else
        {
            logger.assert(isInternal_, "basis function index out of bounds");
            right.basisFunction(i - nLeftBasisFunctions, result);
            return;
        }
    }
    
    template<std::size_t DIM>
    inline void Base::PhysicalFace<DIM>::basisFunction(Side side, std::size_t i, LinearAlgebra::SmallVector<DIM>& result)
    {
        logger.assert(hasPointReference && hasFace, "Need a location to evaluate the data");
        if(side == Side::LEFT)
        {
            left.basisFunction(i, result);
        }
        else
        {
            logger.assert(isInternal_, "cannot find the right element for a boundary face");
            right.basisFunction(i, result);
        }
    }
    
    template<std::size_t DIM>
    inline const LinearAlgebra::SmallVector<DIM>& Base::PhysicalFace<DIM>::basisFunctionCurl(std::size_t i)
    {
        logger.assert(hasPointReference && hasFace, "Need a location to evaluate the data");
        if(i < nLeftBasisFunctions)
        {
            return left.basisFunctionCurl(i);
        }
        else
        {
            logger.assert(isInternal_, "basis function index out of bounds");
            return right.basisFunctionCurl(i - nLeftBasisFunctions);
        }
    }
    
    template<std::size_t DIM>
    inline const LinearAlgebra::SmallVector<DIM>& Base::PhysicalFace<DIM>::basisFunctionCurl(Side side, std::size_t i)
    {
        logger.assert(hasPointReference && hasFace, "Need a location to evaluate the data");
        if(side == Side::LEFT)
        {
            return left.basisFunctionCurl(i);
        }
        else
        {
            logger.assert(isInternal_, "cannot find the right element for a boundary face");
            return right.basisFunctionCurl(i);
        }
    }
    
    template<std::size_t DIM>
    inline void Base::PhysicalFace<DIM>::basisFunctionNormal(std::size_t i, LinearAlgebra::SmallVector<DIM>& result)
    {
        logger.assert(hasPointReference && hasFace, "Need a location to evaluate the data");
        if(hasVectorBasisFunctionNormal)
        {
            result = vectorBasisFunctionNormal_[i];
            return;
        }
        else
        {
            hasVectorBasisFunctionNormal = true;
            for(std::size_t j = 0; j < face_->getNumberOfBasisFunctions(); ++j)
            {
                basisFunction(i, result);
                vectorBasisFunctionNormal_[j] = LinearAlgebra::SmallMatrix<DIM, DIM - 1>{{getNormalVector(), result}}.computeWedgeStuffVector();
                if(j >= nLeftBasisFunctions)
                {
                    vectorBasisFunctionNormal_[j] *= -1.;
                }
            }
            result = vectorBasisFunctionNormal_[i];
            return;
        }
    }
    
    template<std::size_t DIM>
    inline void Base::PhysicalFace<DIM>::basisFunctionNormal(Side side, std::size_t i, LinearAlgebra::SmallVector<DIM>& result)
    {
        logger.assert(hasPointReference && hasFace, "Need a location to evaluate the data");
        if(side == Side::LEFT)
        {
            return basisFunctionNormal(i, result);
        }
        else
        {
            logger.assert(isInternal_, "cannot find the right element for a boundary face");
            return basisFunctionNormal(i + nLeftBasisFunctions, result);
        }
    }
    
    template<std::size_t DIM>
    inline void Base::PhysicalFace<DIM>::basisFunctionUnitNormal(std::size_t i, LinearAlgebra::SmallVector<DIM>& result)
    {
        logger.assert(hasPointReference && hasFace, "Need a location to evaluate the data");
        if(hasVectorBasisFunctionUnitNormal)
        {
            result = vectorBasisFunctionUnitNormal_[i];
            return;
        }
        else
        {
            hasVectorBasisFunctionUnitNormal = true;
            for(std::size_t j = 0; j < face_->getNumberOfBasisFunctions(); ++j)
            {
                basisFunction(i, result);
                vectorBasisFunctionUnitNormal_[j] = LinearAlgebra::SmallMatrix<DIM, DIM - 1>{{getUnitNormalVector(), result}}.computeWedgeStuffVector();
                if(j >= nLeftBasisFunctions)
                {
                    vectorBasisFunctionUnitNormal_[j] *= -1.;
                }
            }
            result = vectorBasisFunctionUnitNormal_[i];
            return;
        }
    }
    
    template<std::size_t DIM>
    inline void Base::PhysicalFace<DIM>::basisFunctionUnitNormal(Side side, std::size_t i, LinearAlgebra::SmallVector<DIM>& result)
    {
        logger.assert(hasPointReference && hasFace, "Need a location to evaluate the data");
        if(side == Side::LEFT)
        {
            return basisFunctionUnitNormal(i, result);
        }
        else
        {
            logger.assert(isInternal_, "cannot find the right element for a boundary face");
            return basisFunctionUnitNormal(i + nLeftBasisFunctions, result);
        }
    }
    
    template<std::size_t DIM>
    inline const LinearAlgebra::MiddleSizeVector& Base::PhysicalFace<DIM>::getSolution(Side side)
    {
        logger.assert(hasPointReference && hasFace, "Need a location to evaluate the data");
        if(side == Side::LEFT)
        {
            return left.getSolution();
        }
        else
        {
            logger.assert(isInternal_, "cannot find the right element for a boundary face");
            return right.getSolution();
        }
    }
    
    template<std::size_t DIM>
    inline void Base::PhysicalFace<DIM>::getSolution(Side side, std::vector<LinearAlgebra::SmallVector<DIM> >& result)
    {
        logger.assert(hasPointReference && hasFace, "Need a location to evaluate the data");
        if(side == Side::LEFT)
        {
            return left.getSolution(result);
        }
        else
        {
            logger.assert(isInternal_, "cannot find the right element for a boundary face");
            return right.getSolution(result);
        }
    }
    
    template<std::size_t DIM>
    inline const std::vector<LinearAlgebra::SmallVector<DIM> >& Base::PhysicalFace<DIM>::getSolutionDeriv(Side side)
    {
        logger.assert(hasPointReference && hasFace, "Need a location to evaluate the data");
        if(side == Side::LEFT)
        {
            return left.getSolutionDeriv();
        }
        else
        {
            logger.assert(isInternal_, "cannot find the right element for a boundary face");
            return right.getSolutionDeriv();
        }
    }
    
    template<std::size_t DIM>
    inline const std::vector<LinearAlgebra::SmallVector<DIM> >& Base::PhysicalFace<DIM>::getSolutionCurl(Side side)
    {
        logger.assert(hasPointReference && hasFace, "Need a location to evaluate the data");
        if(side == Side::LEFT)
        {
            return left.getSolutionCurl();
        }
        else
        {
            logger.assert(isInternal_, "cannot find the right element for a boundary face");
            return right.getSolutionCurl();
        }
    }
    
    template<std::size_t DIM>
    inline std::vector<LinearAlgebra::SmallVector<DIM> > Base::PhysicalFace<DIM>::getSolutionNormal(Side side)
    {
        logger.assert(hasPointReference && hasFace, "Need a location to evaluate the data");
        if(side == Side::LEFT)
        {
            return getNormalVector() * left.getSolution();
        }
        else
        {
            logger.assert(isInternal_, "cannot find the right element for a boundary face");
            return -getNormalVector() * right.getSolution();
        }
    }
    
    template<std::size_t DIM>
    inline std::vector<LinearAlgebra::SmallVector<DIM> > Base::PhysicalFace<DIM>::getSolutionUnitNormal(Side side)
    {
        logger.assert(hasPointReference && hasFace, "Need a location to evaluate the data");
        if(side == Side::LEFT)
        {
            return getUnitNormalVector() * left.getSolution();
        }
        else
        {
            logger.assert(isInternal_, "cannot find the right element for a boundary face");
            return -getUnitNormalVector() * right.getSolution();
        }
    }
    
    template<std::size_t DIM>
    inline void Base::PhysicalFace<DIM>::getSolutionNormal(Side side, LinearAlgebra::SmallVector<DIM>& result)
    {
        logger(ERROR, "Not supported by Element yet");
    }
    
    template<std::size_t DIM>
    inline void Base::PhysicalFace<DIM>::getSolutionUnitNormal(Side side, LinearAlgebra::SmallVector<DIM>& result)
    {
        logger(ERROR, "Not supported by Element yet");
    }
    
    template<std::size_t DIM>
    inline const Geometry::PointReference<DIM - 1>& Base::PhysicalFace<DIM>::getPointReference()
    {
        logger.assert(hasPointReference, "Need a location to evaluate the data");
        return *pointReference_;
    }
    
    template<std::size_t DIM>
    inline const Geometry::PointPhysical<DIM>& Base::PhysicalFace<DIM>::getPointPhysical()
    {
        logger.assert(hasPointReference && hasFace, "Need a location to evaluate the data");
        return left.getPointPhysical();
    }
    
    template<std::size_t DIM>
    inline const LinearAlgebra::SmallVector<DIM>& Base::PhysicalFace<DIM>::getNormalVector()
    {
        logger.assert(hasPointReference && hasFace, "Need a location to evaluate the data");
        if(hasNormal)
        {
            return normal;
        }
        else
        {
            hasNormal = true;
            normal = face_->getNormalVector(*pointReference_);
            return normal;
        }
    }
    
    template<std::size_t DIM>
    inline const LinearAlgebra::SmallVector<DIM>& Base::PhysicalFace<DIM>::getUnitNormalVector()
    {
        logger.assert(hasPointReference && hasFace, "Need a location to evaluate the data");
        if(hasUnitNormal)
        {
            return unitNormal;
        }
        else
        {
            hasUnitNormal = true;
            unitNormal = getNormalVector() / getRelativeSurfaceArea();
            return unitNormal;
        }
    }
    
    template<std::size_t DIM>
    inline double Base::PhysicalFace<DIM>::getRelativeSurfaceArea()
    {
        logger.assert(hasPointReference && hasFace, "Need a location to evaluate the data");
        if(hasNormalNorm)
        {
            return normalNorm;
        }
        else
        {
            hasNormalNorm = true;
            normalNorm = L2Norm(getNormalVector());
            return normalNorm;
        }
    }
    
    template<std::size_t DIM>
    inline Base::FaceMatrix& Base::PhysicalFace<DIM>::getResultMatrix()
    {
        logger.assert(hasPointReference && hasFace, "Need a location to evaluate the data");
        logger.assert(hasFaceMatrix, "Matrix has already been requested for this face/point combination");
        hasFaceMatrix = false;
        return resultMatrix;
    }
    
    template<std::size_t DIM>
    inline LinearAlgebra::MiddleSizeMatrix& Base::PhysicalFace<DIM>::getResultMatrix(Side iSide, Side jSide)
    {
        logger.assert(hasPointReference && hasFace, "Need a location to evaluate the data");
        if(iSide == Side::LEFT)
        {
            if(jSide == Side::LEFT)
            {
                return left.getResultMatrix();
            }
            else
            {
                logger.assert(isInternal_, "cannot find the right element for a boundary face");
                logger.assert(hasLeftRightMatrix, "Matrix has already been requested for this face/point combination");
                hasLeftRightMatrix = false;
                return leftRightMatrix;
            }
        }
        else
        {
            logger.assert(isInternal_, "cannot find the right element for a boundary face");
            if(jSide == Side::LEFT)
            {
                logger.assert(hasRightLeftMatrix, "Matrix has already been requested for this face/point combination");
                hasRightLeftMatrix = false;
                return rightLeftMatrix;
            }
            else
            {
                return right.getResultMatrix();
            }
        }
    }
    
    template<std::size_t DIM>
    inline LinearAlgebra::MiddleSizeVector& Base::PhysicalFace<DIM>::getResultVector()
    {
        logger.assert(hasPointReference && hasFace, "Need a location to evaluate the data");
        logger.assert(hasFaceVector, "Vector has already been requested for this face/point combination");
        hasFaceVector = false;
        return resultVector;
    }
    
    template<std::size_t DIM>
    inline LinearAlgebra::MiddleSizeVector& Base::PhysicalFace<DIM>::getResultVector(Side side)
    {
        logger.assert(hasPointReference && hasFace, "Need a location to evaluate the data");
        if(side == Side::LEFT)
        {
            return left.getResultVector();
        }
        else
        {
            logger.assert(isInternal_, "cannot find the right element for a boundary face");
            return right.getResultVector();
        }
    }
    
    template<std::size_t DIM>
    inline bool Base::PhysicalFace<DIM>::isInternal()
    {
        return isInternal_;
    }
    
    template<std::size_t DIM>
    inline const Base::Face* Base::PhysicalFace<DIM>::getFace()
    {
        logger.assert(hasFace, "Need a location to evaluate the data");
        return face_;
    }
    
    template<std::size_t DIM>
    inline const Base::CoordinateTransformation<DIM>* Base::PhysicalFace<DIM>::getTransform()
    {
        return transform_.get();
    }
    
    template<std::size_t DIM>
    inline void Base::PhysicalFace<DIM>::setPointReference(const Geometry::PointReference<DIM - 1>& point)
    {
        pointReference_ = &point;
        hasPointReference = true;
        if(hasFace)
        {
            left.setPointReference(face_->mapRefFaceToRefElemL(point));
            if(isInternal_)
            {
                right.setPointReference(face_->mapRefFaceToRefElemR(point));
            }
        }
        //even if they are already computed, the information is now out of date
        hasBasisFunctionNormal = false;
        hasBasisFunctionUnitNormal = false;
        hasVectorBasisFunctionNormal = false;
        hasVectorBasisFunctionUnitNormal = false;
        hasSolutionNormal = false;
        hasSolutionUnitNormal = false;
        hasVectorSolutionNormal = false;
        hasVectorSolutionUnitNormal = false;
        hasNormal = false;
        hasUnitNormal = false;
        hasNormalNorm = false;
        if(!hasFaceMatrix)
        {
            resultMatrix *= 0;
        }
        if(!hasFaceVector)
        {
            resultVector *= 0;
        }
        if(!hasLeftRightMatrix)
        {
            leftRightMatrix *= 0;
        }
        if(!hasRightLeftMatrix)
        {
            rightLeftMatrix *= 0;
        }
        hasLeftRightMatrix = true;
        hasRightLeftMatrix = true;
        hasFaceMatrix = true;
        hasFaceVector = true;
    }
    
    template<std::size_t DIM>
    inline void Base::PhysicalFace<DIM>::setFace(const Face* face)
    {
        logger.assert(isInternal_ == face->isInternal(), "This face is not supported by this physical face");
        face_ = face;
        if(!hasFace)
        {
            std::size_t leftCoefficients = face->getPtrElementLeft()->getNumberOfUnknowns() * face->getPtrElementLeft()->getNumberOfBasisFunctions();
            nLeftBasisFunctions = face->getPtrElementLeft()->getNumberOfBasisFunctions();
            std::size_t rightCoefficients = 0;
            std::size_t basisFunctions = face->getPtrElementLeft()->getNumberOfBasisFunctions();
            if(isInternal_)
            {
                basisFunctions += face->getPtrElementRight()->getNumberOfBasisFunctions();
                rightCoefficients = face->getPtrElementRight()->getNumberOfUnknowns() * face->getPtrElementRight()->getNumberOfBasisFunctions();
            }
            resultMatrix.resize(leftCoefficients, rightCoefficients);
            leftRightMatrix.resize(leftCoefficients, rightCoefficients);
            rightLeftMatrix.resize(rightCoefficients, leftCoefficients);
            resultVector.resize(leftCoefficients + rightCoefficients);
            basisFunctionNormal_.resize(basisFunctions);
            vectorBasisFunctionNormal_.resize(basisFunctions);
            basisFunctionUnitNormal_.resize(basisFunctions);
            vectorBasisFunctionUnitNormal_.resize(basisFunctions);
        }
        hasFace = true;
        left.setElement(face->getPtrElementLeft());
        if(isInternal_)
        {
            right.setElement(face->getPtrElementRight());
        }
        if(hasPointReference)
        {
            left.setPointReference(face_->mapRefFaceToRefElemL(*pointReference_));
            if(isInternal_)
            {
                right.setPointReference(face_->mapRefFaceToRefElemR(*pointReference_));
            }
        }
        //even if they are already computed, the information is now out of date
        hasBasisFunctionNormal = false;
        hasBasisFunctionUnitNormal = false;
        hasVectorBasisFunctionNormal = false;
        hasVectorBasisFunctionUnitNormal = false;
        hasSolutionNormal = false;
        hasSolutionUnitNormal = false;
        hasVectorSolutionNormal = false;
        hasVectorSolutionUnitNormal = false;
        hasNormal = false;
        hasUnitNormal = false;
        hasNormalNorm = false;
        if(!hasFaceMatrix)
        {
            resultMatrix *= 0;
        }
        if(!hasFaceVector)
        {
            resultVector *= 0;
        }
        if(!hasLeftRightMatrix)
        {
            leftRightMatrix *= 0;
        }
        if(!hasRightLeftMatrix)
        {
            rightLeftMatrix *= 0;
        }
        hasLeftRightMatrix = true;
        hasRightLeftMatrix = true;
        hasFaceMatrix = true;
        hasFaceVector = true;
    }
    
    template<std::size_t DIM>
    inline void Base::PhysicalFace<DIM>::setTransform(std::shared_ptr<CoordinateTransformation<DIM> >& transform)
    {
        transform_ = transform;
        left.setTransformation(transform);
        if(isInternal_)
        {
            right.setTransformation(transform);
        }
        //even if they are already computed, the information is now out of date
        hasBasisFunctionNormal = false;
        hasBasisFunctionUnitNormal = false;
        hasVectorBasisFunctionNormal = false;
        hasVectorBasisFunctionUnitNormal = false;
        hasSolutionNormal = false;
        hasSolutionUnitNormal = false;
        hasVectorSolutionNormal = false;
        hasVectorSolutionUnitNormal = false;
        if(!hasFaceMatrix)
        {
            resultMatrix *= 0;
        }
        if(!hasFaceVector)
        {
            resultVector *= 0;
        }
        if(!hasLeftRightMatrix)
        {
            leftRightMatrix *= 0;
        }
        if(!hasRightLeftMatrix)
        {
            rightLeftMatrix *= 0;
        }
        hasLeftRightMatrix = true;
        hasRightLeftMatrix = true;
        hasFaceMatrix = true;
        hasFaceVector = true;
    }
    
#endif /* PHYSICALFACE_H_ */
