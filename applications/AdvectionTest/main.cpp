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

#include <cmath>
#include <functional>
#include <chrono>
#include "Base/CommandLineOptions.h"
#include "Base/ConfigurationData.h"
#include "Base/Element.h"
#include "Base/Face.h"
#include "Base/HpgemAPILinear.h"
#include "Base/RectangularMeshDescriptor.h"
#include "Integration/ElementIntegral.h"
#include "Integration/FaceIntegral.h"
#include "Output/TecplotDiscontinuousSolutionWriter.h"
#include "Output/TecplotSingleElementWriter.h"
#include "Utilities/BasisFunctions2DH1ConformingTriangle.h"

#include "Logger.h"

/// Linear advection equation du/dt + a[0] du/dx + a[1] du/dy = 0.
/// This class is meant for testing purposes.
//  Please verify that nobody went and tested a broken feature
//  Please keep the problem modelled here reasonably close to the linear advection problem
///\todo Write self-test
class AdvectionTest : public Base::HpgemAPILinear
{
public:
    ///Constructor. Assign all private variables.
    AdvectionTest(int p) :
    HpgemAPILinear(DIM_, 1, p)
    {
        //Choose the "direction" of the advection.
        //This cannot be implemented with iterators, and since the dimension is
        //not always 2, this is the most generic way to write it.
        a.resize(DIM_);
        for (std::size_t i = 0; i < DIM_; ++i)
        {
            a[i] = 0.1 + 0.1 * i;
        }
    }
    
    /// Create a mesh description
    Base::RectangularMeshDescriptor createMeshDescription(const std::size_t numOfElementPerDirection) override final
    {
        //describes a rectangular domain
        Base::RectangularMeshDescriptor description(DIM_);
        
        //this demo will use the square [0,1]^2
        for (std::size_t i = 0; i < DIM_; ++i)
        {
            description.bottomLeft_[i] = 0;
            description.topRight_[i] = 1;
            //Define elements in each direction.
            description.numElementsInDIM_[i] = numOfElementPerDirection;
            
            //Choose whether you want periodic boundary conditions or other (solid wall)
            //boundary conditions.
            description.boundaryConditions_[i] = Base::BoundaryType::PERIODIC;
        }
        
        return description;
    }
    
    ///Compute phi_i*(a.grad(phi_j)) on a reference point on an element for all
    ///basisfunctions phi_i and phi_j.
    ///You pass the reference point to the basisfunctions. Internally the basisfunctions will be mapped to the physical element
    ///so you wont have to do any transformations yourself
    LinearAlgebra::Matrix computeIntegrandStiffnessMatrixAtElement(const Base::Element *element, const PointReferenceT &point) override final
    {
        std::size_t numBasisFuncs = element->getNrOfBasisFunctions();
        LinearAlgebra::Matrix  result(numBasisFuncs, numBasisFuncs, 0);
        for (std::size_t i = 0; i < numBasisFuncs; ++i)
        {
            for (std::size_t j = 0; j < numBasisFuncs; ++j)
            {
                result(j, i) = element->basisFunction(i, point) * (a * element->basisFunctionDeriv(j, point));
            }
        }
        
        return result;
    }
    
    /// \brief Compute the integrals of the left-hand side associated with faces.
    ///
    ///For every internal face, we want to compute the integral of the flux
    ///for all basisfunctions phi_i and phi_j that are non-zero on that face.
    ///For boundary faces, similar expressions can be obtained depending of the type of boundary condition.
    ///This function will compute these integrands for all basisfunctions phi_i and phi_j
    ///on a certain face at a reference point p. Then the integral can later be computed with appropriate (Gauss-)quadrature rules.
    ///The resulting matrix of values is then given in the matrix integrandVal, to which we passed a reference when calling it.
    ///Please note that you pass a reference point to the basisfunctions and the
    ///transformations are done internally. The class FaceMatrix consists of four element matrices for internal faces and one element matrix for faces on the boundary. Each element matrix corresponds to a pair of two adjacent elements of the face.
    Base::FaceMatrix computeIntegrandStiffnessMatrixAtFace(const Base::Face *face, const LinearAlgebra::NumericalVector &normal, const PointReferenceT &point) override final
    {
        //Get the number of basis functions, first of both sides of the face and
        //then only the basis functions associated with the left and right element.
        std::size_t numBasisFuncs = face->getNrOfBasisFunctions();
        std::size_t nLeft = face->getPtrElementLeft()->getNrOfBasisFunctions();
        std::size_t nRight = 0;
        if(face->isInternal())
        {
            nRight = face->getPtrElementLeft()->getNrOfBasisFunctions();
        }
        
        //Resize the result to the correct size and set all elements to 0.
        Base::FaceMatrix integrandVal(nLeft, nRight);
        integrandVal *= 0;
        
        //Check if the normal is in the same direction as the advection.
        //Note that normal does not have length 1!
        const double A = (a * normal) / Base::L2Norm(normal);
        
        //Compute all entries of the integrand at this point:
        for (std::size_t i = 0; i < numBasisFuncs; ++i)
        {
            Base::Side sideBasisFunction = face->getSide(i);
            for (std::size_t j = 0; j < numBasisFuncs; ++j)
            {
                //Give the terms of the upwind flux.
                //Advection in the same direction as outward normal of the left element:
                if ((A > 1e-12) && (sideBasisFunction == Base::Side::LEFT))
                {
                    integrandVal(j, i) = -(a * face->basisFunctionNormal(j, normal, point)) * face->basisFunction(i, point);
                }
                //Advection in the same direction as outward normal of right element:
                else if ((A < -1e-12) && (sideBasisFunction == Base::Side::RIGHT))
                {
                    integrandVal(j, i) = -(a * face->basisFunctionNormal(j, normal, point)) * face->basisFunction(i, point);
                }
                //Advection orthogonal to normal:
                else if (std::abs(A) < 1e-12)
                {
                    integrandVal(j, i) = -(a * face->basisFunctionNormal(j, normal, point)) * face->basisFunction(i, point) / 2.0;
                }
            }
        }
        
        return integrandVal;
    }
    
    /// Define a solution at time zero.
    double getSolutionAtTimeZero(const PointPhysicalT& point)
    {
        return (std::sin(2 * M_PI * point[0]) * std::sin(2 * M_PI * point[1]));
    }
    
    /// Define the exact solution. In this case that is \f$ u_0(\vec{x}-\vec{a}t) \f$, where \f$ u_0 \f$ is the solution at time zero.
    LinearAlgebra::NumericalVector getExactSolution(const PointPhysicalT& point, const double &time, const std::size_t orderTimeDerivative) override final
    {
        LinearAlgebra::NumericalVector exactSolution(1);
        if(orderTimeDerivative == 0)
        {
            PointPhysicalT displacement(-a*time);
            exactSolution(0) = getSolutionAtTimeZero(point + displacement);
            return exactSolution;
        }
        else
        {
            logger(ERROR, "No exact solution for order time derivative % implemented", orderTimeDerivative);
            exactSolution(0) = 0;
            return exactSolution;
        }
    }
    
    /// Define the initial conditions. In this case it is just the exact solution at the start time.
    LinearAlgebra::NumericalVector getInitialSolution(const PointPhysicalT& point, const double &startTime, const std::size_t orderTimeDerivative) override final
    {
        return getExactSolution(point, startTime, orderTimeDerivative);
    }
    
private:
    
    //Dimension of the problem
    static const std::size_t DIM_;
    
    ///Advective vector
    LinearAlgebra::NumericalVector a;
};

const std::size_t AdvectionTest::DIM_(2);

auto& n = Base::register_argument<std::size_t>('n', "numelems", "Number of Elements", true);
auto& p = Base::register_argument<std::size_t>('p', "poly", "Polynomial order", true);

///Make the problem and solve it.
int main(int argc, char **argv)
{
    Base::parse_options(argc, argv);
    try
    {
        // Choose a mesh type (e.g. TRIANGULAR, RECTANGULAR).
        const Base::MeshType meshType = Base::MeshType::RECTANGULAR;
        
        // Choose variable name(s). Since we have a scalar function, we only need to choose one name.
        std::vector<std::string> variableNames;
        variableNames.push_back("u");
        
        //Construct our problem with n elements in every direction and polynomial order p
        AdvectionTest test(p.getValue());
        
        //Create the mesh
        test.createMesh(n.getValue(), meshType);
        
        // Set the names for the output file
        test.setOutputNames("output", "AdvectionTest", "AdvectionTest", variableNames);
        
        //Run the simulation and write the solution
        
        auto startTime = std::chrono::steady_clock::now();
        
        test.solve(Base::startTime.getValue(), Base::endTime.getValue(), Base::dt.getValue(), Base::numberOfSnapshots.getValue(), true);
        
        auto endTime = std::chrono::steady_clock::now();
        
        logger(INFO, "Simulation took %ms.", std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count());
        
        return 0;
    }
    //If something went wrong, print the error message and return -1.
    catch (const char* e)
    {
        std::cerr << e;
    }
    return -1;
}
