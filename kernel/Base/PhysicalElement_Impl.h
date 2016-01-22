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

namespace Base
{
    template<std::size_t DIM>
    inline double PhysicalElement<DIM>::basisFunction(std::size_t i)
    {
        logger.assert(hasPointReference && hasElement, "Need a location to evaluate the data");
        if (hasFunctionValue)
        {
            return basisFunctionValue[i];
        }
        else
        {
            hasFunctionValue = true;
            for (std::size_t j = 0; j < theElement_->getNumberOfBasisFunctions(); ++j)
            {
                if (hasQuadratureRule)
                {
                    if (doesMapQuadraturePointFromFace)
                    {
                        basisFunctionValue[j] = transform_->transform(theElement_->basisFunction(j, quadratureRule_, quadraturePointIndex_, faceToElementMap_), *this);
                    }
                    else
                    {
                        basisFunctionValue[j] = transform_->transform(theElement_->basisFunction(j, quadratureRule_, quadraturePointIndex_), *this);
                    }
                }
                else
                {
                    basisFunctionValue[j] = transform_->transform(theElement_->basisFunction(j, pointReference_), *this);
                }
            }
            return basisFunctionValue[i];
        }
    }

    template<std::size_t DIM>
    inline const LinearAlgebra::SmallVector<DIM>& PhysicalElement<DIM>::basisFunctionDeriv(std::size_t i)
    {
        logger.assert(hasPointReference && hasElement, "Need a location to evaluate the data");
        if (hasFunctionDeriv)
        {
            return basisFunctionDeriv_[i];
        }
        else
        {
            hasFunctionDeriv = true;
            for (std::size_t j = 0; j < theElement_->getNumberOfBasisFunctions(); ++j)
            {
                if (hasQuadratureRule)
                {
                    if (doesMapQuadraturePointFromFace)
                    {
                        basisFunctionDeriv_[j] = transform_->transformDeriv(theElement_->basisFunctionDeriv<DIM>(j, quadratureRule_, quadraturePointIndex_, faceToElementMap_), *this);
                    }
                    else
                    {
                        basisFunctionDeriv_[j] = transform_->transformDeriv(theElement_->basisFunctionDeriv<DIM>(j, quadratureRule_, quadraturePointIndex_), *this);
                    }
                }
                else
                {
                    basisFunctionDeriv_[j] = transform_->transformDeriv(theElement_->basisFunctionDeriv(j, pointReference_), *this);
                }
            }
            return basisFunctionDeriv_[i];
        }
    }

    template<std::size_t DIM>
    inline void PhysicalElement<DIM>::basisFunction(std::size_t i, LinearAlgebra::SmallVector<DIM>& result)
    {
        logger.assert(hasPointReference && hasElement, "Need a location to evaluate the data");
        if (hasVectorFunctionValue)
        {
            result = vectorBasisFunctionValue[i];
        }
        else
        {
            hasVectorFunctionValue = true;
            for (std::size_t j = 0; j < theElement_->getNumberOfBasisFunctions(); ++j)
            {
                if(hasQuadratureRule)
                {
                    if (doesMapQuadraturePointFromFace)
                    {
                        theElement_->basisFunction(j, quadratureRule_, quadraturePointIndex_, faceToElementMap_, vectorBasisFunctionValue[j]);
                        vectorBasisFunctionValue[j] = transform_->transform(vectorBasisFunctionValue[j], *this);
                    }
                    else
                    {
                        theElement_->basisFunction(j, quadratureRule_, quadraturePointIndex_, vectorBasisFunctionValue[j]);
                        vectorBasisFunctionValue[j] = transform_->transform(vectorBasisFunctionValue[j], *this);
                    }
                }
                else
                {
                    theElement_->basisFunction(j, pointReference_, vectorBasisFunctionValue[j]);
                    vectorBasisFunctionValue[j] = transform_->transform(vectorBasisFunctionValue[j], *this);
                }
            }
            result = vectorBasisFunctionValue[i];
        }
    }

    template<std::size_t DIM>
    inline const LinearAlgebra::SmallVector<DIM>& PhysicalElement<DIM>::basisFunctionCurl(std::size_t i)
    {
        logger.assert(hasPointReference && hasElement, "Need a location to evaluate the data");
        if (hasFunctionCurl)
        {
            return basisFunctionCurl_[i];
        }
        else
        {
            hasFunctionCurl = true;
            for (std::size_t j = 0; j < theElement_->getNumberOfBasisFunctions(); ++j)
            {
                if(hasQuadratureRule)
                {
                    if (doesMapQuadraturePointFromFace)
                    {
                        basisFunctionCurl_[j] = transform_->transformCurl(theElement_->basisFunctionCurl<DIM>(j, quadratureRule_, quadraturePointIndex_, faceToElementMap_), *this);
                    }
                    else
                    {
                        basisFunctionCurl_[j] = transform_->transformCurl(theElement_->basisFunctionCurl<DIM>(j, quadratureRule_, quadraturePointIndex_), *this);
                    }
                }
                else
                {
                    basisFunctionCurl_[j] = transform_->transformCurl(theElement_->basisFunctionCurl(j, pointReference_), *this);
                }
            }
            return basisFunctionCurl_[i];
        }
    }

    template<std::size_t DIM>
    inline const double& PhysicalElement<DIM>::basisFunctionDiv(std::size_t i)
    {
        logger.assert(hasPointReference && hasElement, "Need a location to evaluate the data");
        if (hasFunctionDiv)
        {
            return basisFunctionDiv_[i];
        }
        else
        {
            hasFunctionDiv = true;
            for (std::size_t j = 0; j < theElement_->getNumberOfBasisFunctions(); ++j)
            {
                if(hasQuadratureRule)
                {
                    if (doesMapQuadraturePointFromFace)
                    {
                        basisFunctionDiv_[j] = transform_->transformDiv(theElement_->basisFunctionDiv(j, quadratureRule_, quadraturePointIndex_, faceToElementMap_), *this);
                    }
                    else
                    {
                        basisFunctionDiv_[j] = transform_->transformDiv(theElement_->basisFunctionDiv(j, quadratureRule_, quadraturePointIndex_), *this);
                    }
                }
                else
                {
                    basisFunctionDiv_[j] = transform_->transformDiv(theElement_->basisFunctionDiv(j, pointReference_), *this);
                }
            }
            return basisFunctionDiv_[i];
        }
    }

    template<std::size_t DIM>
    inline const LinearAlgebra::MiddleSizeVector& PhysicalElement<DIM>::getSolution()
    {
        logger.assert(hasPointReference && hasElement, "Need a location to evaluate the data");
        if (hasSolution)
        {
            return solution;
        }
        else
        {
            hasSolution = true;
            solution = theElement_->getSolution(0, *this);
            return solution;
        }
    }

    template<std::size_t DIM>
    inline const std::vector<LinearAlgebra::SmallVector<DIM> >& PhysicalElement<DIM>::getSolutionDeriv()
    {
        logger.assert(hasPointReference && hasElement, "Need a location to evaluate the data");
        if (hasSolutionDeriv)
        {
            return solutionDeriv;
        }
        else
        {
            hasSolutionDeriv = true;
            solutionDeriv = theElement_->getSolutionGradient(0, *this);
            return solutionDeriv;
        }
    }

    template<std::size_t DIM>
    inline void PhysicalElement<DIM>::getSolution(std::vector<LinearAlgebra::SmallVector<DIM> >& result)
    {
        logger(ERROR, "not supported by element yet");
    }

    template<std::size_t DIM>
    inline const std::vector<LinearAlgebra::SmallVector<DIM> >& PhysicalElement<DIM>::getSolutionCurl()
    {
        logger(ERROR, "not supported by element yet");
        return std::vector<LinearAlgebra::SmallVector<DIM> >();
    }

    template<std::size_t DIM>
    inline const LinearAlgebra::MiddleSizeVector& PhysicalElement<DIM>::getSolutionDiv()
    {
        logger(ERROR, "not supported by element yet");
        //just inlining a default-constructed vector complains about stack return despite being dead code
        static LinearAlgebra::MiddleSizeVector dummy;
        return dummy;
    }

    template<std::size_t DIM>
    inline const Geometry::PointReference<DIM>& PhysicalElement<DIM>::getPointReference()
    {
        logger.assert(hasPointReference, "Need a location to evaluate the data");
        return pointReference_;
    }

    template<std::size_t DIM>
    inline const Geometry::PointPhysical<DIM>& PhysicalElement<DIM>::getPointPhysical()
    {
        logger.assert(hasPointReference && hasElement, "Need a location to evaluate the data");
        if (hasPointPhysical)
        {
            return pointPhysical;
        }
        else
        {
            hasPointPhysical = true;
            pointPhysical = theElement_->referenceToPhysical(pointReference_);
            return pointPhysical;
        }
    }

    template<std::size_t DIM>
    inline const Geometry::Jacobian<DIM, DIM>& PhysicalElement<DIM>::getJacobian()
    {
        logger.assert(hasPointReference && hasElement, "Need a location to evaluate the data");
        if (hasJacobian)
        {
            return jacobian;
        }
        else
        {
            hasJacobian = true;
            jacobian = theElement_->calcJacobian(pointReference_);
            return jacobian;
        }
    }

    template<std::size_t DIM>
    inline const Geometry::Jacobian<DIM, DIM>& PhysicalElement<DIM>::getInverseTransposeJacobian()
    {
        logger.assert(hasPointReference && hasElement, "Need a location to evaluate the data");
        if (hasInverseTransposeJacobian)
        {
            return inverseTransposeJacobian;
        }
        else
        {
            hasInverseTransposeJacobian = true;
            inverseTransposeJacobian = getTransposeJacobian().inverse();
            return inverseTransposeJacobian;
        }
    }

    template<std::size_t DIM>
    inline const Geometry::Jacobian<DIM, DIM>& PhysicalElement<DIM>::getTransposeJacobian()
    {
        logger.assert(hasPointReference && hasElement, "Need a location to evaluate the data");
        if (hasTransposeJacobian)
        {
            return transposeJacobian;
        }
        else
        {
            hasTransposeJacobian = true;
            transposeJacobian = getJacobian().transpose();
            return transposeJacobian;
        }
    }

    template<std::size_t DIM>
    inline double PhysicalElement<DIM>::getJacobianAbsDet()
    {
        logger.assert(hasPointReference && hasElement, "Need a location to evaluate the data");
        if (hasJacobianAbsDet)
        {
            return jacobianAbsDet;
        }
        else
        {
            hasJacobianAbsDet = true;
            jacobianAbsDet = std::abs(getJacobianDet());
            return jacobianAbsDet;
        }
    }

    template<std::size_t DIM>
    inline double PhysicalElement<DIM>::getJacobianDet()
    {
        logger.assert(hasPointReference && hasElement, "Need a location to evaluate the data");
        if (hasJacobianDet)
        {
            return jacobianDet;
        }
        else
        {
            hasJacobianDet = true;
            jacobianDet = getJacobian().determinant();
            return jacobianDet;
        }
    }

    template<std::size_t DIM>
    inline LinearAlgebra::MiddleSizeMatrix& PhysicalElement<DIM>::getResultMatrix()
    {
        logger.assert(hasPointReference && hasElement, "Need a location to evaluate the data");
        logger.assert(hasElementMatrix, "Can only provide the matrix once per coordinate");
        hasElementMatrix = false;
        return resultMatrix;
    }

    template<std::size_t DIM>
    inline LinearAlgebra::MiddleSizeVector& PhysicalElement<DIM>::getResultVector()
    {
        logger.assert(hasPointReference && hasElement, "Need a location to evaluate the data");
        logger.assert(hasElementVector, "Can only provide the vector once per coordinate");
        hasElementVector = false;
        return resultVector;
    }

    template<std::size_t DIM>
    inline const Base::Element* PhysicalElement<DIM>::getElement()
    {
        logger.assert(hasElement, "Need a location to evaluate the data");
        return theElement_;
    }

    template<std::size_t DIM>
    inline const Base::CoordinateTransformation<DIM>* PhysicalElement<DIM>::getTransformation()
    {
        return transform_.get();
    }

    template<std::size_t DIM>
    inline void PhysicalElement<DIM>::setPointReference(const Geometry::PointReference<DIM>& point)
    {
        pointReference_ = point;
        hasQuadratureRule = false;
        hasPointReference = true;
        //even if they are already computed, the information is now out of date
        hasFunctionValue = false;
        hasVectorFunctionValue = false;
        hasFunctionDeriv = false;
        hasFunctionCurl = false;
        hasFunctionDiv = false;
        hasSolution = false;
        hasVectorSolution = false;
        hasSolutionDeriv = false;
        hasSolutionCurl = false;
        hasSolutionDiv = false;
        hasPointPhysical = false;
        hasJacobian = false;
        hasTransposeJacobian = false;
        hasInverseTransposeJacobian = false;
        hasJacobianDet = false;
        hasJacobianAbsDet = false;
        if (!hasElementMatrix)
        {
            resultMatrix *= 0;
        }
        if (!hasElementVector)
        {
            resultVector *= 0;
        }
        hasElementMatrix = true;
        hasElementVector = true;
    }

    template<std::size_t DIM>
    inline void PhysicalElement<DIM>::setElement(const Element* element)
    {
        theElement_ = element;
        if (!hasElement)
        {
            std::size_t numberOfEntries = theElement_->getNumberOfBasisFunctions() * theElement_->getNumberOfUnknowns();
            resultMatrix.resize(numberOfEntries, numberOfEntries);
            resultVector.resize(numberOfEntries);
            basisFunctionValue.resize(theElement_->getNumberOfBasisFunctions());
            vectorBasisFunctionValue.resize(theElement_->getNumberOfBasisFunctions());
            basisFunctionDeriv_.resize(theElement_->getNumberOfBasisFunctions());
            basisFunctionCurl_.resize(theElement_->getNumberOfBasisFunctions());
            basisFunctionDiv_.resize(theElement_->getNumberOfBasisFunctions());
        }
        hasElement = true;
        //even if they are already computed, the information is now out of date
        hasFunctionValue = false;
        hasVectorFunctionValue = false;
        hasFunctionDeriv = false;
        hasFunctionCurl = false;
        hasFunctionDiv = false;
        hasSolution = false;
        hasVectorSolution = false;
        hasSolutionDeriv = false;
        hasSolutionCurl = false;
        hasSolutionDiv = false;
        hasPointPhysical = false;
        hasJacobian = false;
        hasTransposeJacobian = false;
        hasInverseTransposeJacobian = false;
        hasJacobianDet = false;
        hasJacobianAbsDet = false;
        if (!hasElementMatrix)
        {
            resultMatrix *= 0;
        }
        if (!hasElementVector)
        {
            resultVector *= 0;
        }
        hasElementMatrix = true;
        hasElementVector = true;
    }

    template<std::size_t DIM>
    inline void PhysicalElement<DIM>::setTransformation(std::shared_ptr<Base::CoordinateTransformation<DIM> >& transform)
    {
        transform_ = transform;
        //even if they are already computed, the information is now out of date
        hasFunctionValue = false;
        hasVectorFunctionValue = false;
        hasFunctionDeriv = false;
        hasFunctionCurl = false;
        hasFunctionDiv = false;
        hasSolution = false;
        hasVectorSolution = false;
        hasSolutionDeriv = false;
        hasSolutionCurl = false;
        hasSolutionDiv = false;
        if (!hasElementMatrix)
        {
            resultMatrix *= 0;
        }
        if (!hasElementVector)
        {
            resultVector *= 0;
        }
        hasElementMatrix = true;
        hasElementVector = true;
    }

    template<std::size_t DIM>
    inline void PhysicalElement<DIM>::setQuadratureRule(QuadratureRules::GaussQuadratureRule* rule)
    {
        quadratureRule_ = rule;
        hasQuadratureRule = true;
        doesMapQuadraturePointFromFace = false;
        setQuadraturePointIndex(0);
    }

    template<std::size_t DIM>
    inline void PhysicalElement<DIM>::setQuadratureRule(QuadratureRules::GaussQuadratureRule* rule, const Geometry::MappingReferenceToReference<1> *map)
    {
        quadratureRule_ = rule;
        faceToElementMap_ = map;
        hasQuadratureRule = true;
        doesMapQuadraturePointFromFace = true;
        setQuadraturePointIndex(0);
    }

    template<std::size_t DIM>
    inline void PhysicalElement<DIM>::setQuadraturePointIndex(std::size_t index)
    {
        logger.assert(hasQuadratureRule, "It only makes sense to index into a quadrature rule if you set one first");
        quadraturePointIndex_ = index;
        if(doesMapQuadraturePointFromFace)
        {
            Geometry::PointReference<DIM - 1> pointOnFace = quadratureRule_->getPoint(index);
            setPointReference(faceToElementMap_->transform(pointOnFace));
        }
        else
        {
            setPointReference(quadratureRule_->getPoint(index));
        }
        hasQuadratureRule = true;
    }
}
