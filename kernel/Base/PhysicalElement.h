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

#ifndef PHYSYCALELEMENT_H_
#define PHYSYCALELEMENT_H_

#include <cstdlib>
#include "Geometry/PointReference.h"
#include "Geometry/PointPhysical.h"
#include "Geometry/Jacobian.h"
#include "CoordinateTransformation.h"
#include "H1ConformingTransformation.h"

namespace Base
{
    template<std::size_t DIM>
    class PhysicalFace;

    //class is final as a reminder that there is no virtual destructor
    //note that none of the functions in here is marked const, because a PhysicalElement reserves the right to alter its internal state to optimize future repeated calls
    //note that names in this class match the names in Element unless this makes no sense
    //when you use a physical element in the kernel be careful to avoid infinite recursion
    ///\TODO generalize implementation to support the cached data
    template<std::size_t DIM>
    class PhysicalElement final
    {
    public:
        PhysicalElement()
            :  transform_((new H1ConformingTransformation<DIM>())), hasPointReference(false), hasElement(false) // other data will get initialized when we have more info
        {
        }

        PhysicalElement(const PhysicalElement& other) = delete;
        PhysicalElement(PhysicalElement&& other) = delete;

        ///value of basis function i at the current reference point
        double basisFunction(std::size_t i);

        ///derivative of basis function i at the current reference point
        const LinearAlgebra::SmallVector<DIM>& basisFunctionDeriv(std::size_t i);


        ///value of basis function i at the current reference point
        void basisFunction(std::size_t i, LinearAlgebra::SmallVector<DIM>& result);

        ///curl of basis function i at the current reference point
        const LinearAlgebra::SmallVector<DIM>& basisFunctionCurl(std::size_t i);


        ///value of the solution at the current reference point at time level 0
        const LinearAlgebra::MiddleSizeVector& getSolution();

        ///derivative of the solution at the current reference point at time level 0
        const std::vector<LinearAlgebra::SmallVector<DIM> >& getSolutionDeriv();


        ///value of the solution at the current reference point at time level 0
        void getSolution(std::vector<LinearAlgebra::SmallVector<DIM> >& result);

        ///curl of the solution at the current reference point at time level 0
        const std::vector<LinearAlgebra::SmallVector<DIM> >& getSolutionCurl();


        ///the current reference point
        const Geometry::PointReference<DIM>& getPointReference();

        ///the current physical point
        const Geometry::PointPhysical<DIM>& getPointPhysical();


        ///the Jacobian of the coordinate transformation
        const Geometry::Jacobian<DIM, DIM>& getJacobian();

        ///the transpose of the inverse of the Jacobian of the coordinate transformation
        const Geometry::Jacobian<DIM, DIM>& getInverseTransposeJacobian();

        ///the transpose of the Jacobian of the coordinate transformation
        const Geometry::Jacobian<DIM, DIM>& getTransposeJacobian();

        ///the absolute value of the determinant of the Jacobian of the coordinate transformation
        double getJacobianAbsDet();


        ///a middle size square matrix of size nBasisFunctions x nUnknowns
        ///\details this gets zeroed out every time the reference point is changed and is only resized by the physical element upon construction, so this could also be used for matrixes of different size
        LinearAlgebra::MiddleSizeMatrix& getResultMatrix();

        ///a middle size vector of size nBasisFunctions x nUnknowns
        ///\details this gets zeroed out every time the reference point is changed and is only resized by the physical element upon construction, so this could also be used for vectors of different size
        LinearAlgebra::MiddleSizeVector& getResultVector();

        ///the element (elements have extra functions for users that need them)
        const Base::Element* getElement();

        ///the transformation that was used to get from the reference element to the physical element
        const Base::CoordinateTransformation<DIM>* getTransformation();

        ///setters should only be needed internally
        void setPointReference(const Geometry::PointReference<DIM>& point);
        ///setters should only be needed internally
        void setElement(const Element* element);
        ///setters should only be needed internally
        void setTransformation(std::shared_ptr<Base::CoordinateTransformation<DIM> >& transform);

    private:

        const Base::Element* theElement_;
        const Geometry::PointReference<DIM>* pointReference_;
        std::shared_ptr<Base::CoordinateTransformation<DIM> > transform_;

        std::vector<double> basisFunctionValue;
        std::vector<LinearAlgebra::SmallVector<DIM> > vectorBasisFunctionValue;
        std::vector<LinearAlgebra::SmallVector<DIM> > basisFunctionDeriv_;
        std::vector<LinearAlgebra::SmallVector<DIM> > basisFunctionCurl_;

        LinearAlgebra::MiddleSizeVector solution;
        std::vector<LinearAlgebra::SmallVector<DIM> > vectorSolution;
        std::vector<LinearAlgebra::SmallVector<DIM> > solutionDeriv;
        std::vector<LinearAlgebra::SmallVector<DIM> > solutionCurl;

        Geometry::PointPhysical<DIM> pointPhysical;
        Geometry::Jacobian<DIM, DIM> jacobian, transposeJacobian, inverseTransposeJacobian;
        double jacobianAbsDet;

        LinearAlgebra::MiddleSizeMatrix resultMatrix;
        LinearAlgebra::MiddleSizeVector resultVector;

        bool hasPointReference, hasElement;

        bool hasElementMatrix, hasElementVector;
        //did we already compute this?
        bool hasFunctionValue, hasVectorFunctionValue, hasFunctionDeriv, hasFunctionCurl;
        bool hasSolution, hasVectorSolution, hasSolutionDeriv, hasSolutionCurl;
        bool hasPointPhysical, hasJacobian, hasTransposeJacobian, hasInverseTransposeJacobian, hasJacobianDet;
    };

    template<std::size_t DIM>
    inline double PhysicalElement<DIM>::basisFunction(std::size_t i)
    {
        logger.assert(hasPointReference && hasElement, "Need a location to evaluate the data");
        if(hasFunctionValue)
        {
            return basisFunctionValue[i];
        }
        else
        {
            hasFunctionValue = true;
            for(std::size_t j = 0; j < theElement_->getNrOfBasisFunctions(); ++j)
            {
                basisFunctionValue[j] = transform_->transform(theElement_->basisFunction(j, *pointReference_), *this);
            }
            return basisFunctionValue[i];
        }
    }

}
    
    template<std::size_t DIM>
    inline const LinearAlgebra::SmallVector<DIM>& Base::PhysicalElement<DIM>::basisFunctionDeriv(std::size_t i)
    {
        logger.assert(hasPointReference && hasElement, "Need a location to evaluate the data");
        if(hasFunctionDeriv)
        {
            return basisFunctionDeriv_[i];
        }
        else
        {
            hasFunctionDeriv = true;
            for(std::size_t j = 0; j < theElement_->getNrOfBasisFunctions(); ++j)
            {
                basisFunctionDeriv_[j] = transform_->transformDeriv(theElement_->basisFunctionDeriv(j, *pointReference_), *this);
            }
            return basisFunctionDeriv_[i];
        }
    }
    
    template<std::size_t DIM>
    inline void Base::PhysicalElement<DIM>::basisFunction(std::size_t i, LinearAlgebra::SmallVector<DIM>& result)
    {
        logger.assert(hasPointReference && hasElement, "Need a location to evaluate the data");
        if(hasVectorFunctionValue)
        {
            result = vectorBasisFunctionValue[i];
            return;
        }
        else
        {
            hasVectorFunctionValue = true;
            for(std::size_t j = 0; j < theElement_->getNrOfBasisFunctions(); ++j)
            {
                theElement_->basisFunction(j, *pointReference_, vectorBasisFunctionValue[j]);
                vectorBasisFunctionValue[j] = transform_->transform(vectorBasisFunctionValue[j], *this);
            }
            result = vectorBasisFunctionValue[i];
            return;
        }
    }
    
    template<std::size_t DIM>
    inline const LinearAlgebra::SmallVector<DIM>& Base::PhysicalElement<DIM>::basisFunctionCurl(std::size_t i)
    {
        logger.assert(hasPointReference && hasElement, "Need a location to evaluate the data");
        if(hasFunctionCurl)
        {
            return basisFunctionCurl_[i];
        }
        else
        {
            hasFunctionCurl = true;
            for(std::size_t j = 0; j < theElement_->getNrOfBasisFunctions(); ++j)
            {
                basisFunctionCurl_[j] = transform_->transformCurl(theElement_->basisFunctionCurl(j, *pointReference_), *this);
            }
            return basisFunctionCurl_[i];
        }
    }
    
    template<std::size_t DIM>
    inline const LinearAlgebra::MiddleSizeVector& Base::PhysicalElement<DIM>::getSolution()
    {
        logger.assert(hasPointReference && hasElement, "Need a location to evaluate the data");
        if(hasSolution)
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
    inline const std::vector<LinearAlgebra::SmallVector<DIM> >& Base::PhysicalElement<DIM>::getSolutionDeriv()
    {
        logger.assert(hasPointReference && hasElement, "Need a location to evaluate the data");
        if(hasSolutionDeriv)
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
    inline void Base::PhysicalElement<DIM>::getSolution(std::vector<LinearAlgebra::SmallVector<DIM> >& result)
    {
        logger(ERROR, "not supported by element yet");
    }
    
    template<std::size_t DIM>
    inline const std::vector<LinearAlgebra::SmallVector<DIM> >& Base::PhysicalElement<DIM>::getSolutionCurl()
    {
        logger(ERROR, "not supported by element yet");
        return std::vector<LinearAlgebra::SmallVector<DIM> >();
    }
    
    template<std::size_t DIM>
    inline const Geometry::PointReference<DIM>& Base::PhysicalElement<DIM>::getPointReference()
    {
        logger.assert(hasPointReference, "Need a location to evaluate the data");
        return *pointReference_;
    }
    
    template<std::size_t DIM>
    inline const Geometry::PointPhysical<DIM>& Base::PhysicalElement<DIM>::getPointPhysical()
    {
        logger.assert(hasPointReference && hasElement, "Need a location to evaluate the data");
        if(hasPointPhysical)
        {
            return pointPhysical;
        }
        else
        {
            hasPointPhysical = true;
            pointPhysical = theElement_->referenceToPhysical(*pointReference_);
            return pointPhysical;
        }
    }
    
    template<std::size_t DIM>
    inline const Geometry::Jacobian<DIM, DIM>& Base::PhysicalElement<DIM>::getJacobian()
    {
        logger.assert(hasPointReference && hasElement, "Need a location to evaluate the data");
        if(hasJacobian)
        {
            return jacobian;
        }
        else
        {
            hasJacobian = true;
            jacobian = theElement_->calcJacobian(*pointReference_);
            return jacobian;
        }
    }
    
    template<std::size_t DIM>
    inline const Geometry::Jacobian<DIM, DIM>& Base::PhysicalElement<DIM>::getInverseTransposeJacobian()
    {
        logger.assert(hasPointReference && hasElement, "Need a location to evaluate the data");
        if(hasInverseTransposeJacobian)
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
    inline const Geometry::Jacobian<DIM, DIM>& Base::PhysicalElement<DIM>::getTransposeJacobian()
    {
        logger.assert(hasPointReference && hasElement, "Need a location to evaluate the data");
        if(hasTransposeJacobian)
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
    inline double Base::PhysicalElement<DIM>::getJacobianAbsDet()
    {
        logger.assert(hasPointReference && hasElement, "Need a location to evaluate the data");
        if(hasJacobianDet)
        {
            return jacobianAbsDet;
        }
        else
        {
            hasJacobianDet = true;
            jacobianAbsDet = std::abs(getJacobian().determinant());
            return jacobianAbsDet;
        }
    }
    
    template<std::size_t DIM>
    inline LinearAlgebra::MiddleSizeMatrix& Base::PhysicalElement<DIM>::getResultMatrix()
    {
        logger.assert(hasPointReference && hasElement, "Need a location to evaluate the data");
        return resultMatrix;
    }
    
    template<std::size_t DIM>
    inline LinearAlgebra::MiddleSizeVector& Base::PhysicalElement<DIM>::getResultVector()
    {
        logger.assert(hasPointReference && hasElement, "Need a location to evaluate the data");
        return resultVector;
    }
    
    template<std::size_t DIM>
    inline const Base::Element* Base::PhysicalElement<DIM>::getElement()
    {
        logger.assert(hasElement, "Need a location to evaluate the data");
        return theElement_;
    }
    
    template<std::size_t DIM>
    inline const Base::CoordinateTransformation<DIM>* Base::PhysicalElement<DIM>::getTransformation()
    {
        return transform_.get();
    }
    
    template<std::size_t DIM>
    inline void Base::PhysicalElement<DIM>::setPointReference(const Geometry::PointReference<DIM>& point)
    {
        pointReference_ = &point;
        hasPointReference = true;
        //even if they are already computed, the information is now out of date
        hasFunctionValue = false;
        hasVectorFunctionValue = false;
        hasFunctionDeriv = false;
        hasFunctionCurl = false;
        hasSolution = false;
        hasVectorSolution = false;
        hasSolutionDeriv = false;
        hasSolutionCurl = false;
        hasPointPhysical = false;
        hasJacobian = false;
        hasTransposeJacobian = false;
        hasInverseTransposeJacobian = false;
        hasJacobianDet = false;
        if(!hasElementMatrix)
        {
            resultMatrix *= 0;
        }
        if(!hasElementVector)
        {
            resultVector *= 0;
        }
        hasElementMatrix = true;
        hasElementVector = true;
    }
    
    template<std::size_t DIM>
    inline void Base::PhysicalElement<DIM>::setElement(const Element* element)
    {
        theElement_ = element;
        if(!hasElement)
        {
            std::size_t numEntries = theElement_->getNrOfBasisFunctions() * theElement_->getNrOfUnknows();
            resultMatrix.resize(numEntries, numEntries);
            resultVector.resize(numEntries);
            basisFunctionValue.resize(theElement_->getNrOfBasisFunctions());
            vectorBasisFunctionValue.resize(theElement_->getNrOfBasisFunctions());
            basisFunctionDeriv_.resize(theElement_->getNrOfBasisFunctions());
            basisFunctionCurl_.resize(theElement_->getNrOfBasisFunctions());
        }
        hasElement = true;
        //even if they are already computed, the information is now out of date
        hasFunctionValue = false;
        hasVectorFunctionValue = false;
        hasFunctionDeriv = false;
        hasFunctionCurl = false;
        hasSolution = false;
        hasVectorSolution = false;
        hasSolutionDeriv = false;
        hasSolutionCurl = false;
        hasPointPhysical = false;
        hasJacobian = false;
        hasTransposeJacobian = false;
        hasInverseTransposeJacobian = false;
        hasJacobianDet = false;
        if(!hasElementMatrix)
        {
            resultMatrix *= 0;
        }
        if(!hasElementVector)
        {
            resultVector *= 0;
        }
        hasElementMatrix = true;
        hasElementVector = true;
    }
    
    template<std::size_t DIM>
    inline void Base::PhysicalElement<DIM>::setTransformation(std::shared_ptr<Base::CoordinateTransformation<DIM> >& transform)
    {
        transform_ = transform;
        //even if they are already computed, the information is now out of date
        hasFunctionValue = false;
        hasVectorFunctionValue = false;
        hasFunctionDeriv = false;
        hasFunctionCurl = false;
        hasSolution = false;
        hasVectorSolution = false;
        hasSolutionDeriv = false;
        hasSolutionCurl = false;
        if(!hasElementMatrix)
        {
            resultMatrix *= 0;
        }
        if(!hasElementVector)
        {
            resultVector *= 0;
        }
        hasElementMatrix = true;
        hasElementVector = true;
    }
    
#endif /* PHYSYCALELEMENT_H_ */
