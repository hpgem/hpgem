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

#ifndef SavageHutterH
#define SavageHutterH

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


struct SHConstructorStruct
{
    std::size_t numOfVariables;
    std::size_t polyOrder;
    std::size_t numElements;
    Base::MeshType meshType;
    Base::ButcherTableau * ptrButcherTableau;
};

//todo: make the functions override final, but at the moment my parser does not 
//understand the override and final keywords, which makes development harder
class SavageHutter : public Base::HpgemAPISimplified<DIM>
{
public:
    SavageHutter(const std::size_t numOfVariables,
            const std::size_t polynomialOrder,
            const Base::ButcherTableau * const ptrButcherTableau);
    
    ///Alternative constructor with less input parameters. Furthermore, this 
    ///constructor also constructs the mesh and couples an object LimiterData to
    ///each element.
    SavageHutter(const SHConstructorStruct& inputValues);

    /// \brief Create a domain
    Base::RectangularMeshDescriptor<DIM> createMeshDescription(const std::size_t numOfElementPerDirection);

    /// \brief Compute the initial solution at a given point in space and time.
    LinearAlgebra::MiddleSizeVector getInitialSolution(const PointPhysicalT &pPhys, const double &startTime, const std::size_t orderTimeDerivative = 0);
    
    /// \brief Show the progress of the time integration.
    void showProgress(const double time, const std::size_t timeStepID)
    {
        if (timeStepID % 1 == 0)
        {
            logger(INFO, "% time steps computed.", timeStepID);
        }
    }

    LinearAlgebra::MiddleSizeVector computeRightHandSideAtElement(Base::Element *ptrElement, LinearAlgebra::MiddleSizeVector &solutionCoefficients, const double time);

    /// \brief Compute the right-hand side corresponding to an internal face
    LinearAlgebra::MiddleSizeVector computeRightHandSideAtFace(Base::Face *ptrFace,
            const Base::Side side,
            LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeft,
            LinearAlgebra::MiddleSizeVector &solutionCoefficientsRight,
            const double time);
    
    LinearAlgebra::MiddleSizeVector computeRightHandSideAtFace
        (
         Base::Face *ptrFace,
         LinearAlgebra::MiddleSizeVector &solutionCoefficients,
         const double time
         );

    ///At the beginning of each time step, it will be checked if a limiter should
    ///be used for this element. If so, it is saved in the LimiterData struct.
    void useLimiterForElement(Base::Element *element);
    
    void computeOneTimeStep(double &time, const double dt);
    void limitSolution();
    
    ///Auxiliary function for checking if a limiter should be used.
    LinearAlgebra::MiddleSizeVector computeVelocity(LinearAlgebra::MiddleSizeVector numericalSolution);
    
    ///Compute the average of the height and discharge in the given element
    LinearAlgebra::MiddleSizeVector computeAverageOfSolution(Base::Element *element);
    
    ///Compute the minimum of the height in the given element
    double getMinimumHeight(const Base::Element *element);
    
    ///If a limiter should be used, use the min-mod limiter. Save the values of 
    ///the left side and right side in the struct LimiterData.
    void limitWithMinMod(Base::Element *element, const std::size_t iVar);
    
    void tasksBeforeSolving() override final
    {
        //todo: for one face integral you used referenceFaceIntegral (which does not scale with the magnitude of the normal) and for the other you used integrate (which does scale)
        //so it is not clear to me whether or not you need scaling. Please fix as needed
        this->faceIntegrator_.setTransformation(std::shared_ptr<Base::CoordinateTransformation<DIM> >(new Base::DoNotScaleIntegrands<DIM>(new Base::H1ConformingTransformation<DIM>())));
        Base::HpgemAPISimplified<DIM>::tasksBeforeSolving();
    }

    int sign(const double x)
    {
        if (std::abs(x) < 1e-10)
            return 0;
        return ((x < 0)? -1 : 1);
    }
    
    LinearAlgebra::MiddleSizeVector projectOnBasisFuns(Base::Element *elt, std::function<double(const PointReferenceT&)> myFun);
    
    void changeHeight( Base::Element* element, double minimum);
    
    ///If all values are approximately 0, set all coefficients to 0.
    void zeroMaker (Base::Element* element);
    
private:

    /// Number of variables
    const std::size_t numOfVariables_;

    SavageHutterRightHandSideComputer rhsComputer_;
    
    friend class SavageHutterRightHandSideComputer;

};

#endif
