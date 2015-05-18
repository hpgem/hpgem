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
#include "LinearAlgebra/MiddleSizeVector.h"
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
     * \todo add the safeguards
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
                currentPoint_(Geometry::PointReferenceFactory::instance()->makePoint(dimension - 1)),
                normal_(dimension),
                useCache_(useCache), recomputeCache_(true), currentPointIndex_(-1)
        {
        }
                
        ShortTermStorageFaceBase(const ShortTermStorageFaceBase &other) = delete;
        
        virtual Face& operator=(Face& face);
        
        virtual void computeData();

        virtual ~ShortTermStorageFaceBase()
        {
            //keep the face alive!
        }
        
        LinearAlgebra::MiddleSizeVector getNormalVector(const ReferencePointT& pRefFace) const override;
        virtual LinearAlgebra::MiddleSizeVector getNormalVector(const ReferencePointT& pRefFace);

        double basisFunction(std::size_t i, const Geometry::PointReference& p) const override
        {
            logger(ERROR, "No storage functionality was implemented! Are you working in a vector valued function space?");
            return 0;
        }

        double basisFunction(Side iSide, std::size_t iBasisFunction, const Geometry::PointReference& p) const override
        {
            logger(ERROR, "No storage functionality was implemented! Are you working in a vector valued function space?");
            return 0;
        }
        
        virtual double basisFunction(std::size_t i, const Geometry::PointReference& p)
        {
            logger(ERROR, "No storage functionality was implemented! Are you working in a vector valued function space?");
            return 0;
        }
        
        virtual double basisFunction(Side iSide, std::size_t iBasisFunction, const Geometry::PointReference& p)
        {
            logger(ERROR, "No storage functionality was implemented! Are you working in a vector valued function space?");
            return 0;
        }
        
        void basisFunction(std::size_t i, const Geometry::PointReference& p, LinearAlgebra::MiddleSizeVector& ret) const override
        {
            logger(ERROR, "No storage functionality was implemented! Are you working in a scalar function space?");
        }
        
        virtual void basisFunction(std::size_t i, const Geometry::PointReference& p, LinearAlgebra::MiddleSizeVector& ret)
        {
            logger(ERROR, "No storage functionality was implemented! Are you working in a scalar function space?");
        }
        
        LinearAlgebra::MiddleSizeVector basisFunctionNormal(std::size_t i, const LinearAlgebra::MiddleSizeVector& normal, const Geometry::PointReference& p) const override
        {
            logger(ERROR, "No storage functionality was implemented! Are you working in an unusual function space?");
            return LinearAlgebra::MiddleSizeVector();
        }
        
        LinearAlgebra::MiddleSizeVector basisFunctionNormal(Side iSide, std::size_t i, const LinearAlgebra::MiddleSizeVector& normal, const Geometry::PointReference& p) const override
        {
            logger(ERROR, "No storage functionality was implemented! Are you working in an unusual function space?");
            return LinearAlgebra::MiddleSizeVector();
        }
        
        virtual LinearAlgebra::MiddleSizeVector basisFunctionNormal(std::size_t i, const LinearAlgebra::MiddleSizeVector& normal, const Geometry::PointReference& p)
        {
            logger(ERROR, "No storage functionality was implemented! Are you working in an unusual function space?");
            return LinearAlgebra::MiddleSizeVector();
        }
        
        virtual LinearAlgebra::MiddleSizeVector basisFunctionNormal(Side iSide, std::size_t i, const LinearAlgebra::MiddleSizeVector& normal, const Geometry::PointReference& p)
        {
            logger(ERROR, "No storage functionality was implemented! Are you working in an unusual function space?");
            return LinearAlgebra::MiddleSizeVector();
        }
        
        double basisFunctionDeriv(std::size_t i, std::size_t jDir, const Geometry::PointReference& p) const override
        {
            return face_->basisFunctionDeriv(i, jDir, p);
        }
        
        LinearAlgebra::MiddleSizeVector basisFunctionDeriv(std::size_t i, const Geometry::PointReference& p) const override
        {
            logger(ERROR, "No storage functionality was implemented! Did you mean basisFunctionCurl?");
            return LinearAlgebra::MiddleSizeVector();
        }
        
        LinearAlgebra::MiddleSizeVector basisFunctionDeriv(Side iSide, std::size_t i, const Geometry::PointReference& p) const override
        {
            logger(ERROR, "No storage functionality was implemented! Did you mean basisFunctionCurl?");
            return LinearAlgebra::MiddleSizeVector();
        }
        
        virtual LinearAlgebra::MiddleSizeVector basisFunctionDeriv(std::size_t i, const Geometry::PointReference& p)
        {
            logger(ERROR, "No storage functionality was implemented! Did you mean basisFunctionCurl?");
            return LinearAlgebra::MiddleSizeVector();
        }
        
        virtual LinearAlgebra::MiddleSizeVector basisFunctionDeriv(Side iSide, std::size_t i, const Geometry::PointReference& p)
        {
            logger(ERROR, "No storage functionality was implemented! Did you mean basisFunctionCurl?");
            return LinearAlgebra::MiddleSizeVector();
        }
        
        LinearAlgebra::MiddleSizeVector basisFunctionCurl(std::size_t i, const Geometry::PointReference& p) const override
        {
            logger(ERROR, "No storage functionality was implemented! Did you mean basisFunctionDeriv?");
            return LinearAlgebra::MiddleSizeVector();
        }
        
        virtual LinearAlgebra::MiddleSizeVector basisFunctionCurl(std::size_t i, const Geometry::PointReference& p)
        {
            logger(ERROR, "No storage functionality was implemented! Did you mean basisFunctionDeriv?");
            return LinearAlgebra::MiddleSizeVector();
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
        
        const Element* getPtrElementLeft() const override final
        {
            return face_->getPtrElementLeft();
        }
        
        Element* getPtrElementLeft() override final
        {
            return face_->getPtrElementLeft();
        }
        
        const Element* getPtrElementRight() const override final
        {
            return face_->getPtrElementRight();
        }
        
        Element* getPtrElementRight() override final
        {
            return face_->getPtrElementRight();
        }
        
        const Element* getPtrElement(Side side) const override final
        {
            return face_->getPtrElement(side);
        }

        const FaceQuadratureRule* getGaussQuadratureRule() const override final
        {
            return face_->getGaussQuadratureRule();
        }
        
        bool isInternal() const override final
        {
            return face_->isInternal();
        }

        std::size_t getNrOfBasisFunctions() const override final
        {
            return face_->getNrOfBasisFunctions();
        }
        
        std::size_t getLocalNrOfBasisFunctions() const override final
        {
            return face_->getLocalNrOfBasisFunctions();
        }
        
        std::size_t getID() const override final
        {
            return face_->getID();
        }
        
        const ElementGeometryT* getElementGLeft() const override final
        {
            return face_->getElementGLeft();
        }
        
        const ElementGeometryT* getPtrElementGRight() const override final
        {
            return face_->getPtrElementGRight();
        }
        
        std::size_t localFaceNumberLeft() const override final
        {
            return face_->localFaceNumberLeft();
        }
        
        std::size_t localFaceNumberRight() const override final
        {
            return face_->localFaceNumberRight();
        }
        
        Geometry::FaceType getFaceType() const override final
        {
            return face_->getFaceType();
        }
        
        std::size_t getFaceToFaceMapIndex() const override final
        {
            return face_->getFaceToFaceMapIndex();
        }
        
        const ReferenceFaceGeometryT* getReferenceGeometry() const override final
        {
            return face_->getReferenceGeometry();
        }
        
        const Geometry::PointReference& mapRefFaceToRefElemL(const ReferencePointT& pRefFace) const override final
        {
            return face_->mapRefFaceToRefElemL(pRefFace);
        }
        
        const Geometry::PointReference& mapRefFaceToRefElemR(const ReferencePointT& pRefFace) const override final
        {
            return face_->mapRefFaceToRefElemR(pRefFace);
        }
        
        const Geometry::PointReference& mapRefFaceToRefFace(const ReferencePointT& pIn) const override final
        {
            return face_->mapRefFaceToRefFace(pIn);
        }
        
        RefFaceToRefElementMappingPtr refFaceToRefElemMapL() const override final
        {
            return face_->refFaceToRefElemMapL();
        }
        
        RefFaceToRefElementMappingPtr refFaceToRefElemMapR() const override final
        {
            return face_->refFaceToRefElemMapR();
        }
        
        LinearAlgebra::MiddleSizeVector getTimeLevelData(std::size_t timeLevel, std::size_t unknown = 0) const override final
        {
            return face_->getTimeLevelData(timeLevel, unknown);
        }

        LinearAlgebra::MiddleSizeMatrix getFaceMatrixMatrix(std::size_t matrixID = 0) const override final
        {
            return face_->getFaceMatrixMatrix(matrixID);
        }
        
        LinearAlgebra::MiddleSizeVector getFaceVector(std::size_t vectorID = 0) const override final
        {
            return face_->getFaceVector(vectorID);
        }
        
        UserFaceData* getUserData() const override final
        {
            return face_->getUserData();
        }
        
        std::size_t convertToSingleIndex(Side side, std::size_t scalarBasisFunctionId, std::size_t varId) const override final
        {
            return face_->convertToSingleIndex(side, scalarBasisFunctionId, varId);
        }
        
        Side getSide(std::size_t faceBasisFunctionId) const override final
        {
            return face_->getSide(faceBasisFunctionId);
        }
        
        std::size_t getElementBasisFunctionId(std::size_t faceBasisFunctionId) const override final
        {
            return face_->getElementBasisFunctionId(faceBasisFunctionId);
        }
        
    private:
        
        ShortTermStorageFaceBase& operator=(const ShortTermStorageFaceBase&)
        {
            logger(ERROR, "you are already storing the data, no need to store it twice!");
            return *this;
        }
        
    protected:
        Face* face_;

        const Geometry::PointReference* currentPoint_;
        LinearAlgebra::MiddleSizeVector normal_;

        std::vector<LinearAlgebra::MiddleSizeVector> basisFunctionValues_, basisFunctionsTimesNormal_, basisFunctionDerivatives_, basisFunctionCurlValues_;
    private:
        
        bool useCache_;
        bool recomputeCache_;
        //currentPointIndex_ is set to -1 in constructor, so can't be a std::size_t currently.
        int currentPointIndex_;
    };

}

#endif /* SHORTTERMSTORAGEFACEBASE_HPP_ */
