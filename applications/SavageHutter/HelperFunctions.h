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
#ifndef HELPERFUNCTIONS_H
#define	HELPERFUNCTIONS_H
#include <cmath>
#include "Base/Element.h"
#include "LinearAlgebra/MiddleSizeVector.h"
#include "Geometry/PointReference.h"
#include "Integration/ElementIntegral.h"
#include "Base/UserData.h"

namespace Helpers
{

    struct DryFlag : public UserElementData
    {
        bool isDry;
    };

    int sign(const double x);

    template <std::size_t DIM>
    LinearAlgebra::MiddleSizeVector getSolution( Base::PhysicalElement<DIM> &element, const LinearAlgebra::MiddleSizeVector &solutionCoefficients, const std::size_t numVariables)
    {
        const std::size_t numBasisFunctions = element.getNumberOfBasisFunctions();
        LinearAlgebra::MiddleSizeVector solution(numVariables);
        for (std::size_t iFun = 0; iFun < numBasisFunctions; ++ iFun)
        {
            const double basisFunctionValue = element.basisFunction(iFun);
            for (std::size_t iVar = 0; iVar < numVariables; ++ iVar)
            {
                std::size_t iVB = element.convertToSingleIndex(iFun, iVar);
                solution(iVar) += basisFunctionValue * solutionCoefficients(iVB);
            }
        }
        return solution;
    }

    ///Compute the average of the height and discharge in the given element
    template <std::size_t DIM>
    LinearAlgebra::MiddleSizeVector computeAverageOfSolution(Base::Element* element, const LinearAlgebra::MiddleSizeVector &solutionCoefficients, Integration::ElementIntegral<DIM>& elementIntegrator)
    {
        std::size_t numOfVariables = element->getNumberOfUnknowns();
        const std::function < LinearAlgebra::MiddleSizeVector(Base::PhysicalElement<DIM>&) > integrandFunction =
            [ = ](Base::PhysicalElement<DIM>& elt) -> LinearAlgebra::MiddleSizeVector {
                LinearAlgebra::MiddleSizeVector solution = Helpers::getSolution<DIM>(elt, solutionCoefficients, numOfVariables);
                logger(DEBUG, "Solution at quadrature point: %", solution);
                return solution;
            };
        LinearAlgebra::MiddleSizeVector average = (elementIntegrator.integrate(element, integrandFunction, element->getGaussQuadratureRule()));
        //\todo: generalise to other than rectangular grid
        Geometry::PointPhysical<DIM> p0 = element->getPhysicalGeometry()->getLocalNodeCoordinates(0);
        Geometry::PointPhysical<DIM> p1 = element->getPhysicalGeometry()->getLocalNodeCoordinates(1);
        average /= Base::L2Norm(p1 - p0);
        if (DIM == 2)
        {            
            Geometry::PointPhysical<DIM> p2 = element->getPhysicalGeometry()->getLocalNodeCoordinates(2);
            average /= Base::L2Norm(p2-p0);
        }

        logger(DEBUG, "Average over element %: %", element->getID(), average);
        logger.assert_always(average(0) > - 1e-16, "Average height negative on "
            "element %! (%), u: %", element->getID(), average, average(1) / average(0));
        return average;
    }

    template <std::size_t DIM>
    LinearAlgebra::MiddleSizeVector projectOnBasisFuns(Base::Element *elt, std::function<double(const Geometry::PointReference<DIM>&) > myFun, Integration::ElementIntegral<DIM>& elementIntegrator)
    {
        const std::size_t numBasisFuns = elt->getNumberOfBasisFunctions();
        LinearAlgebra::MiddleSizeVector projection(numBasisFuns);
        for (std::size_t i = 0; i < numBasisFuns; ++ i)
        {
            const std::function < double(Base::PhysicalElement<DIM>&) > integrandFunction = [ = ](Base::PhysicalElement<DIM>& element) -> double 
            {
                return myFun(element.getPointReference()) * element.basisFunction(i);
            };
            double val = elementIntegrator.integrate(elt, integrandFunction, elt->getGaussQuadratureRule());
            projection[i] = val;
        }
        LinearAlgebra::MiddleSizeMatrix massMatrix(numBasisFuns, numBasisFuns);
        for (std::size_t i = 0; i < numBasisFuns; ++ i)
        {
            for (std::size_t j = 0; j < numBasisFuns; ++ j)
            {
                const std::function < double(Base::PhysicalElement<DIM>&) > massFun = [ = ](Base::PhysicalElement<DIM>& element) -> double 
                {
                    return element.basisFunction(j) * element.basisFunction(i);
                };
                massMatrix(i, j) = elementIntegrator.integrate(elt, massFun);
            }
        }
        massMatrix.solve(projection);
        return projection;
    }
}
#endif	/* HELPERFUNCTIONS_H */

