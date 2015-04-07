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

#ifndef AcousticWaveH
#define AcousticWaveH

#include <fstream>
#include <iomanip> 

#include "Base/CommandLineOptions.h"
#include "Base/ConfigurationData.h"
#include "Base/Element.h"
#include "Base/Face.h"
#include "Base/HpgemAPISimplified.h"
#include "Base/MpiContainer.h"
#include "Base/RectangularMeshDescriptor.h"
#include "Base/TimeIntegration/AllTimeIntegrators.h"
#include "Geometry/PointReference.h"
#include "Integration/ElementIntegral.h"
#include "Integration/FaceIntegral.h"
#include "Output/TecplotDiscontinuousSolutionWriter.h"
#include "Output/TecplotSingleElementWriter.h"
#include "Utilities/BasisFunctions1DH1ConformingLine.h"
#include "Utilities/BasisFunctions2DH1ConformingSquare.h"
#include "Utilities/BasisFunctions2DH1ConformingTriangle.h"
#include "Utilities/BasisFunctions3DH1ConformingCube.h"
#include "Utilities/BasisFunctions3DH1ConformingTetrahedron.h"
#include "SavageHutterRightHandSideComputer.h"

#include "Logger.h"

//todo: make the functions override final, but at the moment my parser does not 
//understand the override and final keywords, which makes development harder
class SavageHutter : public Base::HpgemAPISimplified
{
public:
    SavageHutter
    (
     const std::size_t dimension,
     const std::size_t numOfVariables,
     const std::size_t polynomialOrder,
     const Base::ButcherTableau * const ptrButcherTableau
     );
        
    /// \brief Create a domain
    Base::RectangularMeshDescriptor createMeshDescription(const std::size_t numOfElementPerDirection);

    /// \brief Compute the real solution at a given point in space and time.
    LinearAlgebra::NumericalVector getExactSolution(const PointPhysicalT &pPhys, const double &time, const std::size_t orderTimeDerivative = 0);

    /// \brief Compute the initial solution at a given point in space and time.
    LinearAlgebra::NumericalVector getInitialSolution(const PointPhysicalT &pPhys, const double &startTime, const std::size_t orderTimeDerivative = 0);

    /// \brief Compute the integrand for the reference element for obtaining the initial solution.
    LinearAlgebra::NumericalVector integrandInitialSolutionOnRefElement(const Base::Element *ptrElement, const double &startTime, const Geometry::PointReference &pRef);

    /// \brief Integrate the initial solution for a single element.
    LinearAlgebra::NumericalVector integrateInitialSolutionAtElement(Base::Element * ptrElement, const double startTime, const std::size_t orderTimeDerivative);
    
    /// \brief Show the progress of the time integration.
    void showProgress(const double time, const std::size_t timeStepID)
    {
        if (timeStepID % 100 == 0)
        {
            logger(INFO, "% time steps computed.", timeStepID);
        }
        std::string fileName = "data" + std::to_string(timeStepID);
        std::ofstream myFile(fileName);
        for (Base::Element* element : meshes_[0]->getElementsList())
        {
            Geometry::PointReference pRef(1);
            pRef.setCoordinate(0, -1);
            double pPhys  = (double) element->getID() / meshes_[0]->getElementsList().size();
            myFile << std::setw(5) << pPhys << '\t' << std::setprecision(12) << element->getSolution(0, pRef)(0) << std::endl;  
            pRef.setCoordinate(0, 1);
            pPhys  = (double) (element->getID() + 1) / meshes_[0]->getElementsList().size();
            myFile << std::setw(5) << pPhys << '\t'<< std::setprecision(12) << element->getSolution(0, pRef)(0) << std::endl;  
        }
    }
    
    LinearAlgebra::NumericalVector computeRightHandSideAtElement(Base::Element *ptrElement, LinearAlgebra::NumericalVector &solutionCoefficients, const double time);
    
    /// \brief Compute the right-hand side corresponding to a boundary face
    LinearAlgebra::NumericalVector computeRightHandSideAtFace
    (
     Base::Face *ptrFace,
     LinearAlgebra::NumericalVector &solutionCoefficients,
     const double time
     );
    
    /// \brief Compute the right-hand side corresponding to an internal face
    LinearAlgebra::NumericalVector computeRightHandSideAtFace
    (
     Base::Face *ptrFace,
     const Base::Side side,
     LinearAlgebra::NumericalVector &solutionCoefficientsLeft,
     LinearAlgebra::NumericalVector &solutionCoefficientsRight,
     const double time
     );

private: 
    /// Dimension of the domain
    const std::size_t DIM_;

    /// Number of variables
    const std::size_t numOfVariables_;
    
    SavageHutterRightHandSideComputer rhsComputer_;
        
};

#endif
