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

#include "Base/Face.hpp"

///\bug resolves field has incomplete type
#include "Geometry/PointReference.hpp"
#include "LinearAlgebra/NumericalVector.hpp"

namespace Base {

	/**
	 * An face that computes and stores all basis function values and derivatives.
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
	class ShortTermStorageFaceBase: public Face {
	public:

		ShortTermStorageFaceBase(unsigned int dimension, bool useCache=false) :
			Face(),face_(NULL),               //I dont like that face_ is not defined before operator= is called at least once
			currentPoint_(dimension-1),       //but I want to give users the ability to pass alternative wrappers to the integrators
			normal_(dimension),               //without forcing them to pick a random face that is going to be discarded anyway
			recomputeCache_(true),
			useCache_(useCache),
			currentPointIndex_(-1){}

		virtual Face& operator=(const Face& face){//todo check that &face and this are different things (logger)
			face_=&face;
			if(currentPoint_.size()==0){
				computeData();
			}else{
                 /// \bug This should go back to NAN at some point. Again to fix problems with math and STL::vector
				currentPoint_[0]=1./0.;
			}
			currentPointIndex_=-1;
			return *this;
		}

		virtual void computeData();

		virtual ~ShortTermStorageFaceBase() {
			//keep the face alive!
		}

		virtual void getNormalVector(const ReferencePointOnTheFaceT& pRefFace, LinearAlgebra::NumericalVector& v) const;
		virtual void getNormalVector(const ReferencePointOnTheFaceT& pRefFace, LinearAlgebra::NumericalVector& v);

		virtual double basisFunction(unsigned int i, const Geometry::PointReference& p) const {throw "No storage functionality was implemented! Are you working in a vector valued function space?";}
		virtual double basisFunction(unsigned int i, const Geometry::PointReference& p) {throw "No storage functionality was implemented! Are you working in a vector valued function space?";}

		virtual void basisFunction(unsigned int i, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret) const {throw "No storage functionality was implemented! Are you working in a scalar function space?";}
		virtual void basisFunction(unsigned int i, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret) {throw "No storage functionality was implemented! Are you working in a scalar function space?";}

		virtual void basisFunctionNormal(unsigned int i, const LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret) const {throw "No storage functionality was implemented! Are you working in an unusual function space?";}
		virtual void basisFunctionNormal(unsigned int i, const LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret) {throw "No storage functionality was implemented! Are you working in an unusual function space?";}

		virtual double basisFunctionDeriv(unsigned int i, unsigned int jDir, const Geometry::PointReference& p) const {return face_->basisFunctionDeriv(i,jDir,p);}

		virtual void basisFunctionDeriv(unsigned int i, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret) const {throw "No storage functionality was implemented! Did you mean basisFunctionCurl?";}
		virtual void basisFunctionDeriv(unsigned int i, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret) {throw "No storage functionality was implemented! Did you mean basisFunctionCurl?";}

		virtual void basisFunctionCurl(unsigned int i, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret) const {throw "No storage functionality was implemented! Did you mean basisFunctionDeriv?";}
		virtual void basisFunctionCurl(unsigned int i, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret) {throw "No storage functionality was implemented! Did you mean basisFunctionDeriv?";}

		//if this is needed a lot, also store this
		virtual void referenceToPhysical(const Geometry::PointReference& pointReference, PointPhysicalT& pointPhysical) const;

        //caching functionality

        //! \brief Start caching (geometry) information now.
		virtual void    cacheOn();

		//! \brief Stop using cache.
		virtual void    cacheOff();

		//! \brief Set recompute the cache ON.
		virtual void    recomputeCacheOn();

		//! \brief Set recompute the cache OFF.
		virtual void    recomputeCacheOff();

		//make sure all the other functions map to the other face

		virtual const ElementT* getPtrElementLeft() const {
			return face_->getPtrElementLeft();
		}

		virtual const ElementT* getPtrElementRight() const {
			return face_->getPtrElementRight();
		}

		virtual const FaceQuadratureRule* getGaussQuadratureRule() const {
			return face_->getGaussQuadratureRule();
		}

		virtual bool isInternal() const {
			return face_->isInternal();
		}

		//virtual VecCacheT&       getVecCacheData() { return vecCacheData_; } not sure if ugly or non-const for a reason

		virtual int getNrOfBasisFunctions() const {
			return face_->getNrOfBasisFunctions();
		}

		virtual int getLocalNrOfBasisFunctions() const {
			return face_->getLocalNrOfBasisFunctions();
		}

		virtual int getID() const {
			return face_->getID();
		}

		virtual const ElementGeometryT* getElementGLeft() const {
			return face_->getElementGLeft();
		}

		virtual const ElementGeometryT* getPtrElementGRight() const {
			return face_->getPtrElementGRight();
		}

		virtual unsigned int localFaceNumberLeft() const {
			return face_->localFaceNumberLeft();
		}

		virtual unsigned int localFaceNumberRight() const {
			return face_->localFaceNumberRight();
		}

		virtual Geometry::FaceType getFaceType() const {
			return face_->getFaceType();
		}

		virtual int getFaceToFaceMapIndex() const {
			return face_->getFaceToFaceMapIndex();
		}

		virtual const ReferenceFaceGeometryT* getReferenceGeometry() const {
			return face_->getReferenceGeometry();
		}

		virtual void mapRefFaceToRefElemL(const ReferencePointOnTheFaceT& pRefFace, ReferencePointT& pRefEl) const {
			face_->mapRefFaceToRefElemL(pRefFace, pRefEl);
		}

		virtual void mapRefFaceToRefElemR(const ReferencePointOnTheFaceT& pRefFace, ReferencePointT& pRefEl) const {
			face_->mapRefFaceToRefElemR(pRefFace, pRefEl);
		}

		virtual void mapRefFaceToRefFace(const ReferencePointOnTheFaceT& pIn, ReferencePointOnTheFaceT& pOut) const {
			face_->mapRefFaceToRefFace(pIn, pOut);
		}

		virtual RefFaceToRefElementMapping refFaceToRefElemMapL() const {
			return face_->refFaceToRefElemMapL();
		}

		virtual RefFaceToRefElementMapping refFaceToRefElemMapR() const {
			return face_->refFaceToRefElemMapR();
		}

		virtual void getFaceMatrix(LinearAlgebra::Matrix& matrix, unsigned int matrixID = 0) const {return face_->getFaceMatrix(matrix,matrixID);}

		virtual void getFaceVector(LinearAlgebra::NumericalVector& vector, unsigned int vectorID = 0) const {return face_->getFaceVector(vector,vectorID);}

		virtual const VecCacheT& getVecCacheData() const {
			return face_->FaceData::getVecCacheData();
		}

		virtual UserFaceData* getUserData() const {
			return face_->getUserData();
		}

	private:

		//ShortTermStorageFaceBase(const ShortTermStorageFaceBase&){throw "you are already storing the data, no need to store it twice!";}
		ShortTermStorageFaceBase& operator=(const ShortTermStorageFaceBase&){throw "you are already storing the data, no need to store it twice!";}

	protected:
		const Face* face_;

		Geometry::PointReference currentPoint_;
		LinearAlgebra::NumericalVector normal_;

		std::vector<LinearAlgebra::NumericalVector> basisFunctionValues_,basisFunctionsTimesNormal_,basisFunctionDerivatives_;
	private:

        bool useCache_;
        bool recomputeCache_;
        int currentPointIndex_;
	};

}

#endif /* SHORTTERMSTORAGEFACEBASE_HPP_ */
