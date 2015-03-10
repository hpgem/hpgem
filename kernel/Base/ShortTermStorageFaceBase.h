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

#ifndef SHORTTERMSTORAGEFACEBASE_HPP_
#define SHORTTERMSTORAGEFACEBASE_HPP_

#include "Base/Face.h"

#include "Geometry/PointReference.h"
#include "LinearAlgebra/NumericalVector.h"
#include "Side.h"
#include "BasisFunctionSet.h"


namespace Base
{
    
    /**
     * A face that computes and stores all basis function values and derivatives.
     * This way user can just type ret(j,i)=basisFunctionDeriv(i,p)*basisFunctionDeriv(j,p)
     * without being bothered by unnecessary information about Jacobians or transformations
     * and without having to compute the Jacobian twice for every entry in the face matrix
     *
     * This class will automatically recompute all data whenever a new point is passed to a basisfunction
     *
     * This cannot be directly implemented in Face because then all Faces will store all this data
     * resulting in a massive storage overhead
     *
     * The actual transformations required depend on the function space you are working in, so the actual work
     * is delegated to a subclass that doesn't need to shield all the getters.
     *
     * You cannot use this class to modify faces
     *
     * Be VERY careful to not put this type of face in a mesh, the extra storage needed for this type of faces will likely crash your program
     * Once proper error checking/handling is implemented safeguards will be added to make this a bit more difficult
     */
    class ShortTermStorageFaceBase : public Face
    {
    public:
        //The user should be able to use this as if it were a Face, and it is
        //quicker in integration routines than Face, since it stores the values
        //of the transformed basisfunctions and derivatives of basisfunctions.
        //Note that in the constructor the face on which all operations in this
        //object are performed is a null-pointer until the operator= is called.
        //While this is suboptimal, in practice the user makes only one
        //ShortTermStorageFaceBase which gets assigned all face in a for-loop.
        //Therefore, if the user loops over all the faces anyway, it is not
        //necessary to assign any face in the constructor.
        ShortTermStorageFaceBase(std::size_t dimension, bool useCache = false)
                : Face(), face_(nullptr),
                currentPoint_(dimension - 1),
                normal_(dimension),
                recomputeCache_(true), useCache_(useCache), currentPointIndex_(-1)
        {
        }
        
        virtual Face& operator=(const Face& face);
        
        virtual void computeData();

        virtual ~ShortTermStorageFaceBase() override
        {
            //keep the face alive!
        }
        
        LinearAlgebra::NumericalVector getNormalVector(const ReferencePointT& pRefFace) const override;
        virtual LinearAlgebra::NumericalVector getNormalVector(const ReferencePointT& pRefFace);

        double basisFunction(std::size_t i, const Geometry::PointReference& p) const override
        {
            throw "No storage functionality was implemented! Are you working in a vector valued function space?";
        }
        
        virtual double basisFunction(std::size_t i, const Geometry::PointReference& p)
        {
            throw "No storage functionality was implemented! Are you working in a vector valued function space?";
        }
        
        void basisFunction(std::size_t i, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret) const override
        {
            throw "No storage functionality was implemented! Are you working in a scalar function space?";
        }
        
        virtual void basisFunction(std::size_t i, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret)
        {
            throw "No storage functionality was implemented! Are you working in a scalar function space?";
        }
        
        LinearAlgebra::NumericalVector basisFunctionNormal(std::size_t i, const LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p) const override
        {
            throw "No storage functionality was implemented! Are you working in an unusual function space?";
        }
        
        virtual LinearAlgebra::NumericalVector basisFunctionNormal(std::size_t i, const LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p)
        {
            throw "No storage functionality was implemented! Are you working in an unusual function space?";
        }
        
        double basisFunctionDeriv(std::size_t i, std::size_t jDir, const Geometry::PointReference& p) const override
        {
            return face_->basisFunctionDeriv(i, jDir, p);
        }
        
        LinearAlgebra::NumericalVector basisFunctionDeriv(std::size_t i, const Geometry::PointReference& p) const override
        {
            throw "No storage functionality was implemented! Did you mean basisFunctionCurl?";
        }
        
        virtual LinearAlgebra::NumericalVector basisFunctionDeriv(std::size_t i, const Geometry::PointReference& p)
        {
            throw "No storage functionality was implemented! Did you mean basisFunctionCurl?";
        }
        
        LinearAlgebra::NumericalVector basisFunctionCurl(std::size_t i, const Geometry::PointReference& p) const override
        {
            throw "No storage functionality was implemented! Did you mean basisFunctionDeriv?";
        }
        
        virtual LinearAlgebra::NumericalVector basisFunctionCurl(std::size_t i, const Geometry::PointReference& p)
        {
            throw "No storage functionality was implemented! Did you mean basisFunctionDeriv?";
        }
        
        //if this is needed a lot, also store this
        Geometry::PointPhysical referenceToPhysical(const Geometry::PointReference& pointReference) const override;

        //caching functionality
        
        //! \brief Start caching (geometry) information now.
        void cacheOn();

        //! \brief Stop using cache.
        void cacheOff();

        //! \brief Set recompute the cache ON.
        void recomputeCacheOn();

        //! \brief Set recompute the cache OFF.
        void recomputeCacheOff();

        //make sure all the other functions map to the other face
        
        const ElementT* getPtrElementLeft() const override
        {
            return face_->getPtrElementLeft();
        }
        
        const ElementT* getPtrElementRight() const override
        {
            return face_->getPtrElementRight();
        }
        
        const FaceQuadratureRule* getGaussQuadratureRule() const override
        {
            return face_->getGaussQuadratureRule();
        }
        
        bool isInternal() const override
        {
            return face_->isInternal();
        }

        std::size_t getNrOfBasisFunctions() const override
        {
            return face_->getNrOfBasisFunctions();
        }
        
        std::size_t getLocalNrOfBasisFunctions() const override
        {
            return face_->getLocalNrOfBasisFunctions();
        }
        
        std::size_t getID() const override
        {
            return face_->getID();
        }
        
        const ElementGeometryT* getElementGLeft() const override
        {
            return face_->getElementGLeft();
        }
        
        const ElementGeometryT* getPtrElementGRight() const override
        {
            return face_->getPtrElementGRight();
        }
        
        std::size_t localFaceNumberLeft() const override
        {
            return face_->localFaceNumberLeft();
        }
        
        std::size_t localFaceNumberRight() const override
        {
            return face_->localFaceNumberRight();
        }
        
        Geometry::FaceType getFaceType() const override
        {
            return face_->getFaceType();
        }
        
        std::size_t getFaceToFaceMapIndex() const override
        {
            return face_->getFaceToFaceMapIndex();
        }
        
        const ReferenceFaceGeometryT* getReferenceGeometry() const override
        {
            return face_->getReferenceGeometry();
        }
        
        Geometry::PointReference mapRefFaceToRefElemL(const ReferencePointT& pRefFace) const override
        {
            return face_->mapRefFaceToRefElemL(pRefFace);
        }
        
        Geometry::PointReference mapRefFaceToRefElemR(const ReferencePointT& pRefFace) const override
        {
            return face_->mapRefFaceToRefElemR(pRefFace);
        }
        
        Geometry::PointReference mapRefFaceToRefFace(const ReferencePointT& pIn) const override
        {
            return face_->mapRefFaceToRefFace(pIn);
        }
        
        RefFaceToRefElementMapping refFaceToRefElemMapL() const override
        {
            return face_->refFaceToRefElemMapL();
        }
        
        RefFaceToRefElementMapping refFaceToRefElemMapR() const override
        {
            return face_->refFaceToRefElemMapR();
        }
        
        LinearAlgebra::Matrix getFaceMatrixMatrix(std::size_t matrixID = 0) const override
        {
            return face_->getFaceMatrixMatrix(matrixID);
        }
        
        LinearAlgebra::NumericalVector getFaceVector(std::size_t vectorID = 0) const override
        {
            return face_->getFaceVector(vectorID);
        }
        
        const VecCacheT& getVecCacheData() const override
        {
            return face_->FaceData::getVecCacheData();
        }
        
        UserFaceData* getUserData() const override
        {
            return face_->getUserData();
        }
        
        const std::size_t convertToSingleIndex(Base::Side sideId, std::size_t basisFunctionId, std::size_t unknownId) const override
        {
            return face_->convertToSingleIndex(sideId, basisFunctionId, unknownId);
        }
        
    private:
        
        ShortTermStorageFaceBase(const ShortTermStorageFaceBase&): currentPoint_(0)
        {
            throw "you are already storing the data, no need to store it twice!";
        }
        
        ShortTermStorageFaceBase& operator=(const ShortTermStorageFaceBase&)
        {
            throw "you are already storing the data, no need to store it twice!";
        }
        
    protected:
        const Face* face_;

        Geometry::PointReference currentPoint_;
        LinearAlgebra::NumericalVector normal_;

        std::vector<LinearAlgebra::NumericalVector> basisFunctionValues_, basisFunctionsTimesNormal_, basisFunctionDerivatives_, basisFunctionCurlValues_;
    private:
        
        bool useCache_;
        bool recomputeCache_;
        //currentPointIndex_ is set to -1 in constructor, so can't be a std::size_t currently.
        int currentPointIndex_;
    };

}

#endif /* SHORTTERMSTORAGEFACEBASE_HPP_ */
