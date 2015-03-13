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
#include "Base/CommandLineOptions.h"
#include "Base/ConfigurationData.h"
#include "Base/Element.h"
#include "Base/Face.h"
#include "Base/HpgemUISimplified.h"
#include "Base/ProblemDescriptorTimeDependent.h"
#include "Base/RectangularMeshDescriptor.h"
#include "Base/TimeIntegration/AllTimeIntegrators.h"
#include "Base/TimeIntegration/ButcherTableau.h"
#include "Integration/ElementIntegral.h"
#include "Integration/FaceIntegral.h"
#include "Output/TecplotDiscontinuousSolutionWriter.h"
#include "Output/TecplotSingleElementWriter.h"
#include "Utilities/BasisFunctions2DH1ConformingTriangle.h"
#include "Utilities/BasisFunctions3DH1ConformingTetrahedron.h"

struct errorContainer
{
    errorContainer() = default;
    
    errorContainer(double pointError)
            : infError(std::abs(pointError)), L2Error(pointError * pointError)
    {
    }
    
    //scaling is done as part of the integration and should not be squared
    void axpy(double a, errorContainer x)
    {
        infError = std::max(x.infError, infError);
        L2Error += x.L2Error * a;
    }
    
    errorContainer& operator*=(double a)
    {
        L2Error *= a;
        return *this;
    }
    
    double getInfError()
    {
        return infError;
    }
    
    double getL2Error()
    {
        return std::sqrt(L2Error);
    }
    
private:
    double infError;
    double L2Error;
};

///Linear advection equation du/dt + a[0] du/dx + a[1] du/dy = 0.
///The first self-contained (no PETSc) program to make it into the SVN
//This is a scratch-pad application, please verify that nobody went and tested a broken feature
//please keep the problem modelled here reasonably close to the linear advection problem

class AdvectionTest : public Base::ProblemDescriptorTimeDependent
{
public:
    ///Constructor. Assign all private variables.
    
    AdvectionTest(std::size_t n, std::size_t p, const Base::ButcherTableau* integrator)
            : ProblemDescriptorTimeDependent(DIM_, p, integrator), numElements_(n), polyOrder_(p)
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
    
    ///set up the mesh
    
    bool initialise()
    {
        //describes a rectangular domain
        RectangularMeshDescriptorT description(DIM_);
        
        //this demo will use the square [0,1]^2
        for (std::size_t i = 0; i < DIM_; ++i)
        {
            description.bottomLeft_[i] = 0;
            description.topRight_[i] = 1;
            //Define elements in each direction.
            description.numElementsInDIM_[i] = numElements_;
            
            //Choose whether you want periodic boundary conditions or other (solid wall)
            //boundary conditions.
            description.boundaryConditions_[i] = Base::Boundary::PERIODIC;
        }
        
        //create a triangular mesh. The magic two and three magic ones that are passed to this function
        //specify the number of element matrices, the number of element vectors,
        //the number of face matrices and the number of face vectors (in that order).
        addMesh(description, Base::TRIANGULAR, 2, 1, 1, 1);
        
        //tell hpGEM to use basis functions that are discontinuous and are designed for triangles
        //this is likely to get automated by hpGEM at some point in the future
        meshes_[0]->setDefaultBasisFunctionSet(Utilities::createDGBasisFunctionSet3DH1Tetrahedron(polyOrder_));
        return true;
    }
    
    ///Compute phi_i*(a.grad(phi_j)) on a reference point on an element for all 
    ///basisfunctions phi_i and phi_j.
    ///You pass the reference point to the basisfunctions. Internally the basisfunctions will be mapped to the physical element
    ///so you wont have to do any transformations yourself
    
    void elementIntegrand(const ElementT *element, const PointReferenceT& point, LinearAlgebra::Matrix& result)
    {
        std::size_t numBasisFuncs = element->getNrOfBasisFunctions();
        result.resize(numBasisFuncs, numBasisFuncs);
        LinearAlgebra::NumericalVector phiDerivJ(DIM_);
        for (std::size_t i = 0; i < numBasisFuncs; ++i)
        {
            for (std::size_t j = 0; j < numBasisFuncs; ++j)
            {
                phiDerivJ = element->basisFunctionDeriv(j, point);
                result(j, i) = element->basisFunction(i, point) * (a * phiDerivJ);
            }
        }
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
    ///transformations are done internally. If you expect 4 matrices here, 
    ///you can assume that integrandVal is block structured with 4 blocks in total such
    ///that basisfunctions belonging to the left element are on the left and top.
    
    void faceIntegrand(const FaceT *face, const LinearAlgebra::NumericalVector& normal, const PointReferenceT& point, LinearAlgebra::Matrix& integrandVal)
    {
        //Get the number of basis functions, first of both sides of the face and
        //then only the basis functions associated with the left element.
        std::size_t numBasisFuncs = face->getNrOfBasisFunctions();
        std::size_t nLeft = face->getPtrElementLeft()->getNrOfBasisFunctions();
        
        //Resize the result to the correct size and set all elements to 0.
        integrandVal.resize(numBasisFuncs, numBasisFuncs);
        integrandVal *= 0;
        
        //Check if the normal is in the same direction as the advection.
        //Note that normal does not have length 1!
        const double A = (a * normal) / Base::L2Norm(normal);
        
        LinearAlgebra::NumericalVector phiNormalJ(DIM_);
        
        //Compute all entries of the integrand at this point:
        for (std::size_t i = 0; i < numBasisFuncs; ++i)
        {
            for (std::size_t j = 0; j < numBasisFuncs; ++j)
            {
                //Get phi_j normal_j at this point.
                phiNormalJ = face->basisFunctionNormal(j, normal, point);
                
                //Give the terms of the upwind flux.
                //Advection in the same direction as outward normal of the left element:
                if ((A > 1e-12) && (i < nLeft))
                {
                    integrandVal(j, i) = -(a * phiNormalJ) * face->basisFunction(i, point);
                }
                //Advection in the same direction as outward normal of right element:
                else if ((A < -1e-12) && (i >= nLeft))
                {
                    integrandVal(j, i) = -(a * phiNormalJ) * face->basisFunction(i, point);
                }
                //Advection orthogonal to normal:
                else if (std::abs(A) < 1e-12)
                {
                    integrandVal(j, i) = -(a * phiNormalJ) * face->basisFunction(i, point) / 2.0;
                }
            }
        }
    }
    
    ///The vector edition of the face integrand is meant for implementation of boundary conditions
    ///This is a periodic problem, so it just return 0
    
    void faceIntegrand(const FaceT *face, const LinearAlgebra::NumericalVector& normal, const PointReferenceT& point, LinearAlgebra::NumericalVector& result)
    {
        std::size_t numBasisFuncs = face->getNrOfBasisFunctions();
        result.resize(numBasisFuncs);
        result *= 0;
    }
    
    ///Define the initial conditions, in this case sin(2pi x)* sin(2pi y).
    
    double initialConditions(const PointPhysicalT& point)
    {
        return (std::sin(2 * M_PI * point[0]) * std::sin(2 * M_PI * point[1]));
    }
    
    ///interpolates the initial conditions
    
    void elementIntegrand(const ElementT *element, const PointReferenceT& point, LinearAlgebra::NumericalVector& integrandVal)
    {
        std::size_t numBasisFuncs = element->getNrOfBasisFunctions();
        //Compute the physical coordinates of the reference point
        PointPhysicalT pPhys = element->referenceToPhysical(point);
        
        //Resize the vector and compute the value of the basis function times the
        //value of the initial conditions in this point.
        integrandVal.resize(numBasisFuncs);
        for (std::size_t i = 0; i < numBasisFuncs; ++i)
        {
            integrandVal[i] = element->basisFunction(i, point) * initialConditions(pPhys);
        }
    }
    
    ///provide information about your solution that you want to use for visualisation
    
    void writeToTecplotFile(const ElementT *element, const PointReferenceT& point, std::ostream& out)
    {
        LinearAlgebra::NumericalVector value(1);
        value = element->getSolution(0, point);
        out << value[0];
    }
    
    /// For every element, compute the right hand side of the system of equation,
    /// namely rhs = S * u
    
    void computeRhsLocal()
    {
        LinearAlgebra::Matrix stiffness;
        LinearAlgebra::NumericalVector rhs, oldData;
        for (Base::Element *element : meshes_[0]->getElementsList())
        {
            //collect data
            std::size_t numBasisFuncs = element->getNrOfBasisFunctions();
            stiffness.resize(numBasisFuncs, numBasisFuncs);
            stiffness = element->getElementMatrix(0);
            
            //Get the data of the current time step
            oldData = element->getCurrentData();
            
            //compute rhs = M*u + dt*S*u
            rhs = stiffness * oldData;
            
            //save the rhs in the element
            element->setResidue(rhs);
        }
    }
    
    ///For every face, compute the contribution to the right hand side of the face,
    ///namely rhs = FaceMat * rhs
    
    void computeRhsFaces()
    {
        for (Base::Face *face : meshes_[0]->getFacesList())
        {
            std::size_t numBasisFuncs = face->getNrOfBasisFunctions();
            
            LinearAlgebra::Matrix faceMatrix(numBasisFuncs, numBasisFuncs);
            faceMatrix = face->getFaceMatrixMatrix();
            LinearAlgebra::NumericalVector rhs = face->getCurrentData();
            
            //compute the flux
            rhs = faceMatrix * rhs;
            face->setResidue(rhs);
        }
    }
    
    /// Compute the analytical solution. Hard coded for the initial condition
    /// sin(2pi x)* sin(2pi y).
    
    double analyticalSolution(const PointPhysicalT& point, double time)
    {
        return (std::sin(2 * M_PI * (point[0] - a[0] * time)) * std::sin(2 * M_PI * (point[1] - a[1] * time)));
    }
    
    /// \brief Compute the L^2 error. 
    ///
    /// To compute the error, integrate the square of the difference between the
    /// analytical solution and approximation over every element.
    /// This was originally done by hand and hence very ugly. Automated version still does not look very nice. Feel free to improve.
    /// Furthermore, the maximal difference between the analytical solution and
    /// approximation on Gauss points is computed and printed as Linf error.
    
    double computeError()
    {
        errorContainer total;
        //LinearAlgebra::NumericalVector temp;
        //double analytical, approx;
        double totalError = 0;
        double errorInf = 0;
        Geometry::Jacobian jac(DIM_, DIM_);
        Integration::ElementIntegral integral(false);
        //compute the values of approximation and analytical solution
        for (Base::Element *element : meshes_[0]->getElementsList())
        {
            errorContainer next;
            next = integral.integrate<errorContainer>(element, {[ & ](const Base::Element* element, const Geometry::PointReference& pRef)->errorContainer
            {   
                LinearAlgebra::NumericalVector temp = element->getSolution(0, pRef);
                Geometry::PointPhysical pPhys = element->referenceToPhysical(pRef);
                //the error is the difference between the numerical solution and the analytical solution
                    return analyticalSolution(pPhys, endTime_) - temp[0];
                }});
            total.axpy(1., next);
            //original code in case you are interested
            /*
             const QuadratureRules::GaussQuadratureRule *rule = element->getGaussQuadratureRule();

             for (std::size_t point = 0; point < rule->nrOfPoints(); ++point)
             {
             PointReferenceT pRef(2);
             rule->getPoint(point, pRef);
             element->calcJacobian(pRef, jac);
             element->getSolution(0, pRef, temp);
             approx = temp[0];
             PointPhysicalT pPhys(2);
             element->referenceToPhysical(pRef, pPhys);
             analytical= analyticalSolution(pPhys, endTime_);
             totalError += (approx - analytical) * (approx - analytical) * rule->weight(point) * std::abs(jac.determinant());
             errorInf = std::max(std::abs(approx - analytical), errorInf);
             }
             */
        }
        totalError = total.getL2Error();
        errorInf = total.getInfError();
        std::cout << "Linf error: " << errorInf << std::endl;
        return std::sqrt(totalError);
    }
    
private:
    
    //number of elements per cardinal direction
    std::size_t numElements_;

    //polynomial order of the approximation
    std::size_t polyOrder_;

    //Dimension of the problem
    static const std::size_t DIM_;

    ///Advective vector
    LinearAlgebra::NumericalVector a;
};

const std::size_t AdvectionTest::DIM_(3);

auto& n = Base::register_argument<std::size_t>('n', "numelems", "Number of Elements", true);
auto& p = Base::register_argument<std::size_t>('p', "poly", "Polynomial order", true);

///Make the problem and solve it.

int main(int argc, char **argv)
{
    Base::parse_options(argc, argv);
    try
    {
        //Define the time integrator (ODE method) that we're going to use
        const Base::ButcherTableau *integrator = Base::AllTimeIntegrators::Instance().getRule(4, 4);
        
        //Construct our problem with n elements in every direction and polynomial order p
        AdvectionTest test(n.getValue(), p.getValue(), integrator);
        
        //Define how we want the solution to be written in the VTK files
        test.registerVTKWriteFunction([](Base::Element* element, const Geometry::PointReference& point, std::size_t timelevel) -> double
        {   
            LinearAlgebra::NumericalVector solution(1);
            solution = element->getSolution(timelevel, point);
            return solution[0];
        }, "value");
        
        //Run the simulation and write the solution
        test.solve();
        double error = test.computeError();
        
        std::cout << "L2-error: " << error << std::endl;
        
        return 0;
    }
    //If something went wrong, print the error message and return -1.
    catch (const char* e)
    {
        std::cout << e;
    }
    return -1;
}
