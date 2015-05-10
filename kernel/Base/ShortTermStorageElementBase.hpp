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
#include "Node.hpp"

namespace Base
{

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
    class ShortTermStorageElementBase : public Element
    {
    public:

        ShortTermStorageElementBase(std::size_t dimension, bool useCache = false) :
        Element(), //the superclass is not meant for actual use
        element_(nullptr), //I dont like that face_ is not defined before operator= is called at least once
        currentPoint_(dimension), //but I want to give users the ability to pass alternative wrappers to the integrators
        jac_(dimension, dimension), //without forcing them to pick a random face that is going to be discarded anyway
        recomputeCache_(true),
        currentPointIndex_(-1),
        useCache_(useCache) { }


        ///recomputes the jacobian, the physical point, functionvalues and derivatives of functions based on the current point
        virtual void computeData();



        Element& operator=(const Element& element)
        {//todo check that &element and this are different things (errorChecker)
            element_ = &element;
            /// \bug This should go back to NAN at some point. Again to fix problems with math and STL::vector
            ///\bug placing 0/0 here breaks one of the PETSc based self tests
            currentPoint_[0] = 1. / 0.;
            currentPointIndex_ = -1;
            return *this;
        }

        ShortTermStorageElementBase(const ShortTermStorageElementBase& copy) : element_(copy.element_), currentPoint_(copy.currentPoint_), jac_(copy.jac_), useCache_(copy.useCache_), currentPointIndex_(copy.currentPointIndex_), recomputeCache_(copy.recomputeCache_) { }

        ~ShortTermStorageElementBase() {
            //keep the element alive!
        }

        virtual double basisFunction(std::size_t i, const PointReferenceT& p)
        {
            throw "No storage functionality was implemented! Are you working in a vector valued function space?";
        }

        virtual void basisFunction(std::size_t i, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret)
        {
            throw "No storage functionality was implemented! Are you working in a scalar function space?";
        }

        double basisFunctionDeriv(std::size_t i, std::size_t jDir, const PointReferenceT& p) const override
        {
            return element_->basisFunctionDeriv(i, jDir, p);
        }

        virtual LinearAlgebra::NumericalVector basisFunctionDeriv(std::size_t i, const PointReferenceT& p, const Element* = nullptr)
        {
            throw "No storage functionality was implemented! Did you mean basisFunctionCurl?";
        }

        virtual LinearAlgebra::NumericalVector basisFunctionCurl(std::size_t i, const PointReferenceT& p)
        {
            throw "No storage functionality was implemented! Did you mean basisFunctionDeriv?";
        }

        double basisFunction(std::size_t i, const PointReferenceT& p) const override
        {
            throw "No storage functionality was implemented! Are you working in a vector valued function space?";
        }

        void basisFunction(std::size_t i, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret) const override
        {
            throw "No storage functionality was implemented! Are you working in a scalar function space?";
        }

        LinearAlgebra::NumericalVector basisFunctionDeriv(std::size_t i, const PointReferenceT& p, const Element* = nullptr) const override
        {
            throw "No storage functionality was implemented! Did you mean basisFunctionCurl?";
        }

        LinearAlgebra::NumericalVector basisFunctionCurl(std::size_t i, const PointReferenceT& p) const override
        {
            throw "No storage functionality was implemented! Did you mean basisFunctionDeriv?";
        }

        virtual Geometry::Jacobian calcJacobian(const PointReferenceT& pointReference);

        Geometry::Jacobian calcJacobian(const PointReferenceT& pointReference) const override;

        //if this is needed a lot, also store this
        Geometry::PointPhysical referenceToPhysical(const PointReferenceT& pointReference)const override;

        //caching functionality

        //! \brief Start caching (geometry) information now.
        void cacheOn();

        //! \brief Stop using cache.
        void cacheOff();

        //! \brief Set recompute the cache ON.
        void recomputeCacheOn();

        //! \brief Set recompute the cache OFF.
        void recomputeCacheOff();

        //make sure all the other function map to the current element

        std::size_t getID()const override
        {
            return element_->getID();
        }

        std::size_t getID() override
        {
            return element_->getID();
        }

        const GaussQuadratureRuleT* getGaussQuadratureRule() const override
        {
            return element_->getGaussQuadratureRule();
        }

        Element::SolutionVector getSolution(std::size_t timeLevel, const PointReferenceT& p) const override
        {
            return element_->getSolution(timeLevel, p);
        }

        std::size_t getLocalNrOfBasisFunctions() const override
        {
            return element_->getLocalNrOfBasisFunctions();
        }

        const Face* getFace(std::size_t localFaceNr)const override
        {
            return element_->getFace(localFaceNr);
        }

        const Edge* getEdge(std::size_t localEdgeNr)const override
        {
            return element_->getEdge(localEdgeNr);
        }

        const Node* getNode(std::size_t localNodeNr)const override
        {
            return element_->getNode(localNodeNr);
        }

        std::size_t getNrOfFaces() const override
        {
            return element_->getNrOfFaces();
        }

        std::size_t getNrOfEdges() const override
        {
            return element_->getNrOfEdges();
        }

        std::size_t getNrOfNodes() const override
        {
            return element_->getNrOfNodes();
        }

#ifndef NDEBUG

        const Base::BaseBasisFunction* getBasisFunction(std::size_t i)const override
        {
            return element_->getBasisFunction(i);
        }
#endif

        const LinearAlgebra::Matrix & getElementMatrix(std::size_t matrixID = 0) const override
        {
            return element_->getElementMatrix(matrixID);
        }

        LinearAlgebra::NumericalVector getElementVector(std::size_t vectorID = 0) const override
        {
            return element_->getElementVector(vectorID);
        }

        const LinearAlgebra::NumericalVector getTimeLevelData(std::size_t timeLevel, std::size_t unknown = 0) const override
        {
            return element_->getTimeLevelData(timeLevel, unknown);
        }

        double getData(std::size_t timeLevel, std::size_t unknown, std::size_t basisFunction) const override
        {
            return element_->getData(timeLevel, unknown, basisFunction);
        }

        std::size_t getNrOfUnknows() const override
        {
            return element_->getNrOfUnknows();
        }

        std::size_t getNrOfBasisFunctions() const override
        {
            return element_->getNrOfBasisFunctions();
        }

        const LinearAlgebra::NumericalVector& getResidue() const override
        {
            return element_->getResidue();
        }

        UserElementData* getUserData() const override
        {
            return element_->getUserData();
        }

        const MappingReferenceToPhysicalT * const getReferenceToPhysicalMap() const override
        {
            return element_->getReferenceToPhysicalMap();
        }

        const PhysicalGeometryT * const getPhysicalGeometry() const override
        {
            return element_->getPhysicalGeometry();
        }

        const ReferenceGeometryT * const getReferenceGeometry() const override
        {
            return element_->getReferenceGeometry();
        }

        const RefinementGeometryT* getRefinementGeometry() const override
        {
            return element_->getRefinementGeometry();
        }

        const std::size_t convertToSingleIndex(std::size_t basisFunctionId, std::size_t unknownId) const override
        {
            return element_->convertToSingleIndex(basisFunctionId, unknownId);
        }

    private:

        ShortTermStorageElementBase& operator=(const ShortTermStorageElementBase&)
        {
            throw "you are already storing the data, no need to store it twice!";
        }

    protected:
        const Element* element_;
        Geometry::PointReference currentPoint_;

        Geometry::Jacobian jac_;

        std::vector<LinearAlgebra::NumericalVector> basisFunctionValues_, basisFunctionDerivatives_, basisFunctionCurlValues_;
    private:

        bool useCache_;
        bool recomputeCache_;
        //currentPointIndex is set to -1 in constructor, so can't be a std::size_t currently.
        int currentPointIndex_;
    };
}



#endif /* SHORTTERMSTORAGEELEMENT_HPP_ */
