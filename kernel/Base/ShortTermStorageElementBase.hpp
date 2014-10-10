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


#ifndef SHORTTERMSTORAGEELEMENT_HPP_
#define SHORTTERMSTORAGEELEMENT_HPP_

#include "Base/Element.hpp"

///\BUG resolves field has incomplete type
#include "Geometry/PointReference.hpp"
#include "Geometry/Jacobian.hpp"

namespace Base{

	class Element;

	/**
	 * An element that computes and stores all basis function values and derivatives.
	 * This way user can just type ret(j,i)=basisFunctionDeriv(i,p)*basisFunctionDeriv(j,p)
	 * without being bothered by unnecessary information about Jacobians or transformations
	 * and without having to compute the Jacobian twice for every entry in the element matrix
	 *
	 * This class will automatically recompute all data whenever a new point is passed to a basisfunction
	 *
	 * This cannot be directly implemented in Element because then all Element will store all this data
	 * resulting in a massive storage overhead
	 *
	 * The actual transformations required depend on the function space you are working in, so the actual work
	 * is delegated to a subclass that doesn't need to shield all the getters.
	 *
	 * You cannot use this class to modify elements
	 *
	 * Be VERY careful to not put this type of element in a mesh, the extra storage needed for this type of elements will likely crash your program
	 * Once proper error checking/handling is implemented safeguards will be added to make this a bit more difficult
	 */
	class ShortTermStorageElementBase:public Element{
	public:

		ShortTermStorageElementBase(unsigned int dimension, bool useCache=false):
			Element(),                        //the superclass is not meant for actual use
			element_(NULL),                   //I dont like that face_ is not defined before operator= is called at least once
			currentPoint_(dimension),         //but I want to give users the ability to pass alternative wrappers to the integrators
			jac_(dimension,dimension),        //without forcing them to pick a random face that is going to be discarded anyway
			recomputeCache_(true),
			currentPointIndex_(-1),
			useCache_(useCache){}

		///recomputes the jacobian, the physical point, functionvalues and derivatives of functions based on the current point
		virtual void computeData();

		Element& operator=(const Element& element){//todo check that &element and this are different things (errorChecker)
			element_=&element;
			currentPoint_[0]=1./0.;
			currentPointIndex_=-1;
			return *this;
		}

		ShortTermStorageElementBase(const ShortTermStorageElementBase& copy):element_(copy.element_),currentPoint_(copy.currentPoint_),jac_(copy.jac_),useCache_(copy.useCache_),currentPointIndex_(copy.currentPointIndex_),recomputeCache_(copy.recomputeCache_){}

		~ShortTermStorageElementBase(){
			//keep the element alive!
		}

		virtual double                          basisFunction(unsigned int i, const PointReferenceT& p)  {throw "No storage functionality was implemented! Are you working in a vector valued function space?";}

		virtual void                            basisFunction(unsigned int i, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret) {throw "No storage functionality was implemented! Are you working in a scalar function space?";}

		virtual double                          basisFunctionDeriv(unsigned int i, unsigned int jDir, const PointReferenceT& p) const {return element_->basisFunctionDeriv(i,jDir,p);}

		virtual void                            basisFunctionDeriv(unsigned int i,const PointReferenceT& p, LinearAlgebra::NumericalVector& ret,const Element* =NULL) {throw "No storage functionality was implemented! Did you mean basisFunctionCurl?";}

		virtual void                            basisFunctionCurl(unsigned int i, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret) {throw "No storage functionality was implemented! Did you mean basisFunctionDeriv?";}

		virtual double basisFunction(unsigned int i, const PointReferenceT& p) const {throw "No storage functionality was implemented! Are you working in a vector valued function space?";}
		virtual void   basisFunction(unsigned int i, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret) const {throw "No storage functionality was implemented! Are you working in a scalar function space?";}
		virtual void   basisFunctionDeriv(unsigned int i, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret,const Element* =NULL) const {throw "No storage functionality was implemented! Did you mean basisFunctionCurl?";}
		virtual void   basisFunctionCurl (unsigned int i, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret) const {throw "No storage functionality was implemented! Did you mean basisFunctionDeriv?";}

        virtual void                                        calcJacobian(const PointReferenceT& pointReference, JacobianT& jacobian);

        virtual void                                        calcJacobian(const PointReferenceT& pointReference, JacobianT& jacobian) const;

		//if this is needed a lot, also store this
        virtual void                                        referenceToPhysical(const PointReferenceT& pointReference, PointPhysicalT& pointPhysical)const;

        //caching functionality

        //! \brief Start caching (geometry) information now.
		virtual void    cacheOn();

		//! \brief Stop using cache.
		virtual void    cacheOff();

		//! \brief Set recompute the cache ON.
		virtual void    recomputeCacheOn();

		//! \brief Set recompute the cache OFF.
		virtual void    recomputeCacheOff();

        //make sure all the other function map to the current element

		virtual unsigned int                    getID()const{return element_->getID();}

		virtual unsigned int                    getID(){return element_->getID();}

		virtual const GaussQuadratureRuleT*     getGaussQuadratureRule() const{return element_->getGaussQuadratureRule();}

		//virtual VecCacheT&                      getVecCacheData(){return element_->getVecCacheData();}//Not sure if ugly or non-const for a reason

		virtual void                            getSolution(unsigned int timeLevel, const PointReferenceT& p, SolutionVector& solution) const{return element_->getSolution(timeLevel,p,solution);}

		virtual int                             getLocalNrOfBasisFunctions() const{return element_->getLocalNrOfBasisFunctions();}

		virtual int                             getLocalNrOfBasisFunctionsVertex() const{return element_->getLocalNrOfBasisFunctionsVertex();}

		virtual const Face*                     getFace(int localFaceNr)const {return element_->getFace(localFaceNr);}

		virtual const Edge*                     getEdge(int localEdgeNr)const {return element_->getEdge(localEdgeNr);}

		virtual int                             getNrOfEdges() const{return element_->getNrOfEdges();}

#ifndef NDEBUG
		virtual const Base::BaseBasisFunction*  getBasisFunction(int i)const {return element_->getBasisFunction(i);}
#endif

        virtual void getElementMatrix(LinearAlgebra::Matrix& mat, int matrixID=0) const {element_->getElementMatrix(mat,matrixID);}

        virtual void getElementVector(LinearAlgebra::NumericalVector& vec, int vectorID=0) const{element_->getElementVector(vec,vectorID);}

        virtual const LinearAlgebra::Matrix&    getTimeLevelData(unsigned int timeLevel) const{return element_->getTimeLevelData(timeLevel);}

        virtual double                    getData(unsigned int timeLevel, unsigned int unknown, unsigned int basisFunction) const{return element_->getData(timeLevel,unknown,basisFunction);}

        virtual int                       getNrOfUnknows() const {return element_->getNrOfUnknows();}

        virtual int                       getNrOfBasisFunctions() const {return element_->getNrOfBasisFunctions();}

        virtual const VectorOfDoubles&          getResidue() const {return element_->getResidue();}

        virtual UserElementData*          getUserData() const {return element_->getUserData();}

        virtual const MappingReferenceToPhysicalT* const    getReferenceToPhysicalMap() const{return element_->getReferenceToPhysicalMap();}

        virtual const PhysicalGeometryT* const              getPhysicalGeometry() const{return element_->getPhysicalGeometry();}

        virtual unsigned int                                getNrOfNodes() const{return element_->getNrOfNodes();}

        virtual const ReferenceGeometryT* const             getReferenceGeometry() const{return element_->getReferenceGeometry();}

        virtual const RefinementGeometryT*                  getRefinementGeometry() const{return element_->getRefinementGeometry();}

	private:

		ShortTermStorageElementBase& operator=(const ShortTermStorageElementBase&){throw "you are already storing the data, no need to store it twice!";}

	protected:
		const Element* element_;
		Geometry::PointReference currentPoint_;

		Geometry::Jacobian jac_;

		std::vector<LinearAlgebra::NumericalVector> basisFunctionValues_,basisFunctionDerivatives_;
	private:

        bool useCache_;
        bool recomputeCache_;
        int currentPointIndex_;
	};
}



#endif /* SHORTTERMSTORAGEELEMENT_HPP_ */
