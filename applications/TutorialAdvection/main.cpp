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

#include "Base/HpgemUISimplified.hpp"
#include "Output/TecplotSingleElementWriter.hpp"
#include "Base/ElementCacheData.hpp"
#include "Base/FaceCacheData.hpp"
#include "Output/TecplotDiscontinuousSolutionWriter.hpp"
#include "Utilities/BasisFunctions2DH1ConformingTriangle.hpp"
#include "Base/Element.hpp"
#include "Integration/ReturnTrait1.hpp"
#include "Base/Element.hpp"
#include "Integration/FaceIntegral.hpp"
#include "Integration/ElementIntegral.hpp"
#include "Base/ConfigurationData.hpp"
#include "Base/RectangularMeshDescriptor.hpp"
#include "Base/ElementCacheData.hpp"
#include "Base/FaceCacheData.hpp"
#include "cassert"
#include "Base/CommandLineOptions.hpp"
#include "Output/TecplotDiscontinuousSolutionWriter.hpp"
#include <cmath>
#include <functional>



///Linear advection equation du/dt + a[0] du/dx + a[1] du/dy = 0.
///The first self-contained (no PETSc) program to make it into the SVN
///\todo Polish into example code or tutorial
///\todo Write self-test
class Advection : public Base::HpgemUISimplified
{
public:
    ///Constructor. Assign all private variables.
    Advection(int n, int p) : HpgemUISimplified(DIM_, p), numElements_(n), polyOrder_(p)
    {
        //Choose the "direction" of the advection.
        //This cannot be implemented with iterators, and since the dimension is
        //not always 2, this is the most generic way to write it.
        a.resize(DIM_);
        for (size_t i = 0; i < DIM_; ++i)
        {
            a[i] = 0.1 + 0.1*i;
        }
        
    }

    ///set up the mesh
    bool initialise() 
    {
        //describes a rectangular domain
        RectangularMeshDescriptorT description(DIM_);

        //this demo will use the square [0,1]^2
        for (size_t i = 0; i < DIM_; ++i)
        {
            description.bottomLeft_[i] = 0;
            description.topRight_[i] = 1;
            //Define elements in each direction.
            description.numElementsInDIM_[i] = numElements_;

            //Choose whether you want periodic boundary conditions or other (solid wall)
            //boundary conditions.
            description.boundaryConditions_[i] = RectangularMeshDescriptorT::PERIODIC;
        }

        //create a triangular mesh. The magic two and three magic ones that are passed to this function
        //specify the number of element matrices, the number of element vectors,
        //the number of face matrices and the number of face vectors (in that order).
        addMesh(description, Base::TRIANGULAR, 2, 1, 1, 1);

        //tell hpGEM to use basis functions that are discontinuous and are designed for triangles
        //this is likely to get automated by hpGEM at some point in the future
        meshes_[0]->setDefaultBasisFunctionSet(Utilities::createDGBasisFunctionSet2DH1Triangle(polyOrder_));
        return true;
    }

    ///Compute phi_i*(a.grad(phi_j)) on a reference point on an element for all 
    ///basisfunctions phi_i and phi_j.
    ///You pass the reference point to the basisfunctions. Internally the basisfunctions will be mapped to the physical element
    ///so you wont have to do any transformations yourself
    void elementIntegrand(const ElementT* element, const PointReferenceT& point, LinearAlgebra::Matrix& result)
    {
        size_t numBasisFuncs = element->getNrOfBasisFunctions();
        result.resize(numBasisFuncs, numBasisFuncs);
        LinearAlgebra::NumericalVector phiDerivJ(DIM_);
        for (size_t i = 0; i < numBasisFuncs; ++i)
        {
            for (size_t j = 0; j < numBasisFuncs; ++j)
            {
                element->basisFunctionDeriv(j, point, phiDerivJ);
                result(j, i) = element->basisFunction(i, point)*(a * phiDerivJ);
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
    void faceIntegrand(const FaceT* face, const LinearAlgebra::NumericalVector& normal, const PointReferenceOnTheFaceT& point, LinearAlgebra::Matrix& integrandVal)
    {
        //Get the number of basis functions, first of both sides of the face and
        //then only the basis functions associated with the left element.
        size_t numBasisFuncs = face->getNrOfBasisFunctions();
        size_t nLeft = face->getPtrElementLeft()->getNrOfBasisFunctions();
        
        //Resize the result to the correct size and set all elements to 0.
        integrandVal.resize(numBasisFuncs, numBasisFuncs);
        integrandVal *= 0;

        //Check if the normal is in the same direction as the advection.
        //Note that normal does not have length 1!
        double A = (a * normal) / Base::L2Norm(normal);
        
        LinearAlgebra::NumericalVector phiNormalJ(DIM_);
        
        //Compute all entries of the integrand at this point:
        for (size_t i = 0; i < numBasisFuncs; ++i)
        {
            for (size_t j = 0; j < numBasisFuncs; ++j)
            {
                //Get phi_j normal_j at this point.
                face->basisFunctionNormal(j, normal, point, phiNormalJ);

                //Give the terms of the upwind flux.
                //Advection in the same direction as outward normal of the left element:
                if ((A > 1e-12)&&(i < nLeft)) 
                {
                    integrandVal(j, i) = -(a * phiNormalJ) * face->basisFunction(i, point);
                }
                //Advection in the same direction as outward normal of right element:
                else if ((A<-1e-12)&&(i >= nLeft))
                {
                    integrandVal(j, i) = -(a * phiNormalJ) * face->basisFunction(i, point);
                }
                //Advection orthogonal to normal:
                else if (std::abs(A) < 1e-12)
                {
                    integrandVal(j, i) = -(a * phiNormalJ) * face->basisFunction(i, point) / 2.;
                }
            }
        }
    }

    ///The vector edition of the face integrand is meant for implementation of boundary conditions
    ///This is a periodic problem, so it just return 0
    void faceIntegrand(const FaceT* face, const LinearAlgebra::NumericalVector& normal, const PointReferenceOnTheFaceT& point, LinearAlgebra::NumericalVector& result)
    {
        size_t n = face->getNrOfBasisFunctions();
        result.resize(n);
        result *= 0;
    }
    
    ///Define the initial conditions, in this case sin(2pi x)* sin(2pi y).
    double initialConditions(const PointPhysicalT& point)
    {
        return std::sin(2 * M_PI * point[0]) * std::sin(2 * M_PI * point[1]);
    }

    ///interpolates the initial conditions
    void elementIntegrand(const ElementT* element, const PointReferenceT& point, LinearAlgebra::NumericalVector& integrandVal)
    {
        //Compute the physical coordinates of the reference point
        PointPhysicalT pPhys(DIM_);
        element->referenceToPhysical(point, pPhys);
        
        //Resize the vector and compute the value of the basis function times the
        //value of the initial conditions in this point.
        integrandVal.resize(element->getNrOfBasisFunctions());
        for (size_t i = 0; i < element->getNrOfBasisFunctions(); ++i)
        {
            integrandVal[i] = element->basisFunction(i, point) * initialConditions(pPhys);
        }
    }

    ///provide information about your solution that you want to use for visualisation
    void writeToTecplotFile(const ElementT* element, const PointReferenceT& point, std::ostream& out)
    {
        LinearAlgebra::NumericalVector value(1);
        element->getSolution(0, point, value);
        out << value[0];
    }

    ///TODO put this in HpgemUISimplified
    // solve Mx=`residue`
    void interpolate()
    {
        LinearAlgebra::NumericalVector solution;
        for (Base::Element* element : meshes_[0]->getElementsList())
        {
            size_t numBasisFuncs = element->getNrOfBasisFunctions();
            LinearAlgebra::Matrix mass = element->getMassMatrix();
            solution = element->getResidue();
            mass.solve(solution);
            element->setTimeLevelData(0, solution);
        }
    }    
    
    void computeRhsLocal()
    {
        LinearAlgebra::Matrix stiffness;
        LinearAlgebra::NumericalVector residual, oldData;
        for (Base::Element* element : meshes_[0]->getElementsList())
        {
            //collect data
            //first number of basis functions, then the mass matrix and finally
            //the stiffness matrix
            size_t numBasisFuncs = element->getNrOfBasisFunctions();
            stiffness.resize(numBasisFuncs, numBasisFuncs);
            element->getElementMatrix(stiffness, 0);
            LinearAlgebra::Matrix mass = element->getMassMatrix();
            
            //Get the data of the initial time
            //At the moment, oldData has 1 row, numBasisFuncs columns
            oldData = element->getTimeLevelData(0);
            
            //compute residual=M*u+dt*S*u
            residual = (mass + (dt_ * stiffness)) * oldData;
            
            element->setResidue(residual);
        }
    }
    
    void computeRhsFaces()
    {
        LinearAlgebra::Matrix faceMatrix;
        for (Base::Face* face : meshes_[0]->getFacesList())
        {
            size_t numBasisFuncs = face->getNrOfBasisFunctions();
            size_t numBasisFuncsLeft = face->getPtrElementLeft()->getNrOfBasisFunctions();
            faceMatrix.resize(numBasisFuncs, numBasisFuncs);
            face->getFaceMatrix(faceMatrix);            
            LinearAlgebra::NumericalVector residue = face->getTimeLevelData(0);

            //compute the flux
            residue = faceMatrix*residue;
            residue *= dt_;
            face->setResidue(residue);
        }
    }
    
    ///This function contains everything that has to be done before time-integration
    ///can start. In this case, nothing.
    void beforeTimeIntegration()
    {
    }


private:

    //number of elements per cardinal direction
    int numElements_;

    //polynomial order of the approximation
    int polyOrder_;

    //Dimension of the problem
    static const unsigned int DIM_;
    
    ///advective terms.
    LinearAlgebra::NumericalVector a;
};

const unsigned int Advection::DIM_(2);

auto& n = Base::register_argument<std::size_t>('n', "numelems", "Number of Elements", true);
auto& p = Base::register_argument<std::size_t>('p', "poly", "Polynomial order", true);
int main(int argc, char **argv)
{
    Base::parse_options(argc, argv);
    try
    {

        Advection test(n.getValue(), p.getValue());
        test.registerVTKWriteFunction([](Base::Element* element, const Geometry::PointReference& point, size_t timelevel)->double
        {
            LinearAlgebra::NumericalVector solution(1);
            element->getSolution(timelevel, point, solution);
            return solution[0];
        }, "value");
        test.solve();
        return 0;
    }
    catch (const char* e)
    {
        std::cout << e;
    }
}

