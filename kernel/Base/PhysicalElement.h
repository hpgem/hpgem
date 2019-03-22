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
#include "Element.h"

namespace Base
{
    template<std::size_t DIM>
    class PhysicalFace;
    template<std::size_t DIM>
    class CoordinateTransformation;

    //class is final as a reminder that there is no virtual destructor
    //note that none of the functions in here is marked const, because a PhysicalElement reserves the right to alter its internal state to optimize future repeated calls
    //note that names in this class match the names in Element unless this makes no sense
    //when you use a physical element in the kernel be careful to avoid infinite recursion
    ///\todo generalize implementation to support the cached data
    template<std::size_t DIM>
    class PhysicalElement final
    {
    public:
        PhysicalElement();

        PhysicalElement(const PhysicalElement &other) = delete;

        PhysicalElement(PhysicalElement &&other) = delete;

        ///value of basis function i at the current reference point
        double basisFunction(std::size_t i);
        double basisFunction(std::size_t i, std::size_t unknown);

        ///derivative of basis function i at the current reference point
        const LinearAlgebra::SmallVector<DIM> &basisFunctionDeriv(std::size_t i);
        const LinearAlgebra::SmallVector<DIM> &basisFunctionDeriv(std::size_t i, std::size_t unknown);


        ///value of basis function i at the current reference point
        void basisFunction(std::size_t i, LinearAlgebra::SmallVector<DIM> &result);
        void basisFunction(std::size_t i, LinearAlgebra::SmallVector<DIM> &result, std::size_t unknown);

        ///curl of basis function i at the current reference point
        const LinearAlgebra::SmallVector<DIM> &basisFunctionCurl(std::size_t i);
        const LinearAlgebra::SmallVector<DIM> &basisFunctionCurl(std::size_t i, std::size_t unknown);

        ///divergence of basis function i at the current reference point
        const double &basisFunctionDiv(std::size_t i);
        const double &basisFunctionDiv(std::size_t i, std::size_t unknown);

        
        ///value of the solution at the current reference point at time level 0
        const LinearAlgebra::MiddleSizeVector &getSolution();

        ///derivative of the solution at the current reference point at time level 0
        const std::vector<LinearAlgebra::SmallVector<DIM> > &getSolutionDeriv();


        ///value of the solution at the current reference point at time level 0
        void getSolution(std::vector<LinearAlgebra::SmallVector<DIM> > &result);

        ///curl of the solution at the current reference point at time level 0
        const std::vector<LinearAlgebra::SmallVector<DIM> > &getSolutionCurl();

        ///divergence of the solution at the current reference point at time level 0
        const LinearAlgebra::MiddleSizeVector &getSolutionDiv();


        ///the current reference point
        const Geometry::PointReference<DIM> &getPointReference();

        ///the current physical point
        const Geometry::PointPhysical<DIM> &getPointPhysical();


        ///the Jacobian of the coordinate transformation
        const Geometry::Jacobian<DIM, DIM> &getJacobian();

        ///the transpose of the inverse of the Jacobian of the coordinate transformation
        const Geometry::Jacobian<DIM, DIM> &getInverseTransposeJacobian();

        ///the transpose of the Jacobian of the coordinate transformation
        const Geometry::Jacobian<DIM, DIM> &getTransposeJacobian();

        ///the absolute value of the determinant of the Jacobian of the coordinate transformation
        double getJacobianAbsDet();

        ///the determinant of the Jacobian of the coordinate transformation
        double getJacobianDet();

        ///a middle size square matrix of size nBasisFunctions x nUnknowns
        ///\details this gets zeroed out every time the reference point is changed and is only resized by the physical element upon construction, so this could also be used for matrixes of different size
        LinearAlgebra::MiddleSizeMatrix &getResultMatrix();

        ///a middle size vector of size nBasisFunctions x nUnknowns
        ///\details this gets zeroed out every time the reference point is changed and is only resized by the physical element upon construction, so this could also be used for vectors of different size
        LinearAlgebra::MiddleSizeVector &getResultVector();

        ///the element (elements have extra functions for users that need them)
        const Base::Element *getElement();

        ///\deprecated Does not conform naming conventions, use getNumberOfBasisFunctions instead
        std::size_t getNumOfBasisFunctions()
        {
            return getNumberOfBasisFunctions();
        }

        ///the number of basis functions that are nonzero in the element
        std::size_t getNumberOfBasisFunctions()
        {
            return theElement_->getNumberOfBasisFunctions();
        }
        
        std::size_t getNumberOfBasisFunctions(std::size_t unknown)
        {
            return theElement_->getNumberOfBasisFunctions(unknown);
        }
        
        std::size_t getTotalNumberOfBasisFunctions()
        {
            return theElement_->getTotalNumberOfBasisFunctions();
        }

        ///\deprecated Does not conform naming conventions, use getNumberOfUnknowns instead
        std::size_t getNumOfUnknowns()
        {
            return getNumberOfUnknowns();
        }

        ///the number of unknowns present in the problem
        std::size_t getNumberOfUnknowns()
        {
            return theElement_->getNumberOfUnknowns();
        }

        ///the id of the element
        std::size_t getID()
        {
            return theElement_->getID();
        }

        ///combine a function index and a variable index to a single index that can be used for indexing matrices or vectors
        std::size_t convertToSingleIndex(std::size_t functionId, std::size_t variableId)
        {
            //currently calling the function is too fast to be worth storing the variable
            return theElement_->convertToSingleIndex(functionId, variableId);
        }

        ///the transformation that was used to get from the reference element to the physical element (should only be needed internally)
        const Base::CoordinateTransformation<DIM> *getTransformation();
        
        const Base::CoordinateTransformation<DIM> *getTransformation(std::size_t unknown);

        ///setters should only be needed internally
        void setPointReference(const Geometry::PointReference<DIM> &point);

        ///setters should only be needed internally
        void setElement(const Element *element);

        ///setters should only be needed internally
        void setTransformation(std::shared_ptr<Base::CoordinateTransformation<DIM> > &transform, std::size_t unknown = 0);
        
        ///setters should only be needed internally
        ///sets a quadrature rule that provides points in the element
        void setQuadratureRule(QuadratureRules::GaussQuadratureRule *rule);

        ///setters should only be needed internally
        ///sets a quadrature rule that provides points on the face and a mapping that maps then into the element
        void setQuadratureRule(QuadratureRules::GaussQuadratureRule *rule, const Geometry::MappingReferenceToReference<1> *map);

        ///setters should only be needed internally
        void setQuadraturePointIndex(std::size_t index);

    private:

        const Base::Element *theElement_;
        Geometry::PointReference<DIM> pointReference_;
        QuadratureRules::GaussQuadratureRule *quadratureRule_;
        const Geometry::MappingReferenceToReference<1> *faceToElementMap_;
        std::size_t quadraturePointIndex_;
        std::vector<std::shared_ptr<Base::CoordinateTransformation<DIM> > > transform_;

        std::vector<std::vector<double> > basisFunctionValue;
        std::vector<std::vector<LinearAlgebra::SmallVector<DIM> > > vectorBasisFunctionValue;
        std::vector<std::vector<LinearAlgebra::SmallVector<DIM> > > basisFunctionDeriv_;
        std::vector<std::vector<LinearAlgebra::SmallVector<DIM> > > basisFunctionCurl_;
        std::vector<std::vector<double> > basisFunctionDiv_;

        LinearAlgebra::MiddleSizeVector solution;
        std::vector<LinearAlgebra::SmallVector<DIM> > vectorSolution;
        std::vector<LinearAlgebra::SmallVector<DIM> > solutionDeriv;
        std::vector<LinearAlgebra::SmallVector<DIM> > solutionCurl;
        LinearAlgebra::MiddleSizeVector solutionDiv;

        Geometry::PointPhysical<DIM> pointPhysical;
        Geometry::Jacobian<DIM, DIM> jacobian, transposeJacobian, inverseTransposeJacobian;
        double jacobianAbsDet;
        double jacobianDet;

        LinearAlgebra::MiddleSizeMatrix resultMatrix;
        LinearAlgebra::MiddleSizeVector resultVector;

        bool hasPointReference, hasElement, hasQuadratureRule, doesMapQuadraturePointFromFace;

        bool hasElementMatrix, hasElementVector;
        //did we already compute this?
        std::vector<bool> hasFunctionValue;
        std::vector<bool> hasVectorFunctionValue;
        std::vector<bool> hasFunctionDeriv;
        std::vector<bool> hasFunctionCurl;
        std::vector<bool> hasFunctionDiv;        
        bool hasSolution, hasVectorSolution, hasSolutionDeriv, hasSolutionCurl, hasSolutionDiv;
        bool hasPointPhysical, hasJacobian, hasTransposeJacobian, hasInverseTransposeJacobian, hasJacobianDet, hasJacobianAbsDet;
    };

}

#include "PhysicalElement_Impl.h"
#endif /* PHYSYCALELEMENT_H_ */
