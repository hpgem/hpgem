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

//This macro is no longer needed, but it may make your life easier if you want to change things that involve the link between PETSc and hpGEM in an IDE 
//#define hpGEM_INCLUDE_PETSC_SUPPORT

#include "Base/MpiContainer.hpp"
#include "Base/HpgemUISimplified.hpp"
#include "Base/Norm2.hpp"
#include "Utilities/GlobalMatrix.hpp"
#include "Utilities/GlobalVector.hpp"
#include "petscksp.h"
#include "Output/TecplotSingleElementWriter.hpp"
#include <fstream>
#include "Geometry/ReferenceTetrahedron.hpp"
#include "Utilities/BasisFunctions1DH1ConformingLine.hpp"
#include "Utilities/BasisFunctions2DH1ConformingSquare.hpp"
#include "Utilities/BasisFunctions2DH1ConformingTriangle.hpp"
#include "Utilities/BasisFunctions3DH1ConformingCube.hpp"
#include "Utilities/BasisFunctions3DH1ConformingTetrahedron.hpp"
#include "Base/ElementCacheData.hpp"
#include "Base/FaceCacheData.hpp"
#include "Output/TecplotDiscontinuousSolutionWriter.hpp"
#include "Base/Element.hpp"
#include "Base/RectangularMeshDescriptor.hpp"
#include "Integration/ElementIntegral.hpp"
#include "Base/CommandLineOptions.hpp"
#include "Output/VTKSpecificTimeWriter.hpp"

class Laplace : public Base::HpgemUISimplified
{
public:
    ///constructor: set the dimension of the problem, start the API of hpGEM and initialise the other fields
    //initialisation order is not fixed so DIM_ has to be passed to hpGEMUISimplified in a hardcoded way
    Laplace(int n, int p) : HpgemUISimplified(2, p), DIM_(2), n_(n), p_(p)
    {
        penaltyParameter_ = 3 * n_ * p_ * (p_ + DIM_ - 1) + 1;
    }

    ///set up the mesh
    bool virtual initialise()
    {
        //describes a rectangular domain
        RectangularMeshDescriptorT description(DIM_);

        //this demo will use a cube
        for (int i = 0; i < DIM_; ++i)
        {
            description.bottomLeft_[i] = 0;
            description.topRight_[i] = 1;
            description.numElementsInDIM_[i] = n_;

            //at the moment your options are SOLID_WALL and PERIODIC
            //once we decide what names of boundary conditions to support
            //it will become possible to appropriate boundary conditions
            //for your problem here
            description.boundaryConditions_[i] = RectangularMeshDescriptorT::SOLID_WALL;
        }

        //create a triangular mesh. The four magic ones that are passed to this function
        //specify the number of element matrices, the number of element vectors,
        //the number of face matrices and the number of face vectors (in that order)
        addMesh(description, Base::TRIANGULAR, 1, 1, 1, 1);

        //tell hpGEM to use basis functions that are discontinuous and are designed for triangles
        //this is likely to get automated by hpGEM at some point in the future
        meshes_[0]->setDefaultBasisFunctionSet(Utilities::createDGBasisFunctionSet2DH1Triangle(p_));
        return true;
    }

    //using conforming basis functions is more work
    //this will definitely be automated at some point in the future
    //this routine is only required if you use conforming basis functions
    void useConformingBasisFunctions()
    {
        //if you want to use conforming basis functions you have to create the bubble functions first (even if there are none)
        //take extra care here that the shape of your elements matches the shape the basis functions 'expect'
        meshes_[0]->setDefaultBasisFunctionSet(Utilities::createInteriorBasisFunctionSet2DH1Triangle(p_));
        //then you can create and assign the other basis functions in any order.
        std::vector<const Base::BasisFunctionSet*> bFsets;
        Utilities::createVertexBasisFunctionSet2DH1Triangle(p_, bFsets);
        meshes_[0]->addVertexBasisFunctionSet(bFsets);
        std::vector<const Base::OrientedBasisFunctionSet*> oBFsets;
        Utilities::createFaceBasisFunctionSet2DH1Triangle(p_, oBFsets);
        meshes_[0]->addFaceBasisFunctionSet(oBFsets);
    }

    ///You pass the reference point to the basis functions. Internally the basis functions will be mapped to the physical element
    ///so you wont have to do any transformations yourself
    virtual void elementIntegrand(const ElementT* element, const PointReferenceT& point, LinearAlgebra::Matrix& result)
    {
        //be careful that n changes meaning inside integration routines to represent the number of basis functions involved
        int n = element->getNrOfBasisFunctions();
        LinearAlgebra::NumericalVector phiDerivI(DIM_), phiDerivJ(DIM_);
        result.resize(n, n);
        for (int i = 0; i < n; ++i)
        {
            element->basisFunctionDeriv(i, point, phiDerivI);
            for (int j = 0; j < n; ++j)
            {
                element->basisFunctionDeriv(j, point, phiDerivJ);

                //the value to compute here is taken from the weak formulation
                result(j, i) = phiDerivI * phiDerivJ;
            }
        }
    }

    ///You pass the reference point to the basis functions. Internally the basis functions will be mapped to the physical element
    ///so you wont have to do any transformations yourself. If you expect 4 matrices here, you can assume that ret is block structured such
    ///that basis functions belonging to the left element are indexed first
    ///note that using a consistent flux has no effect if you also use conforming basis functions
    //this routine is only needed if you use discontinuous basis functions
    virtual void faceIntegrand(const FaceT* face, const LinearAlgebra::NumericalVector& normal, const PointReferenceOnTheFaceT& point, LinearAlgebra::Matrix& result)
    {
        int n = face->getNrOfBasisFunctions();
        result.resize(n, n);
        LinearAlgebra::NumericalVector phiNormalI(DIM_), phiNormalJ(DIM_), phiDerivI(DIM_), phiDerivJ(DIM_);
        PointPhysicalT pPhys(DIM_);
        face->referenceToPhysical(point, pPhys);
        for (int i = 0; i < n; ++i)
        {
            face->basisFunctionNormal(i, normal, point, phiNormalI);
            face->basisFunctionDeriv(i, point, phiDerivI);
            for (int j = 0; j < n; ++j)
            {
                face->basisFunctionNormal(j, normal, point, phiNormalJ);
                face->basisFunctionDeriv(j, point, phiDerivJ);
                if (face->isInternal())
                {
                    result(j, i) = -(phiNormalI * phiDerivJ + phiNormalJ * phiDerivI) / 2
                        + penaltyParameter_ * phiNormalI * phiNormalJ;

                    //for the moment you can figure out at what part of the boundary you are by looking at
                    //the point in physical space you are (pPhys) and implement a boundary condition accordingly
                    //once there are more labels for boundary conditions available you can use syntax that looks like
                    //face->getFaceType() to see what kind of boundary condition you are dealing with
                }
                else if (std::abs(pPhys[0]) < 1e-9 || std::abs(pPhys[0] - 1.) < 1e-9) //Dirichlet
                {
                    result(j, i) = -(phiNormalI * phiDerivJ + phiNormalJ * phiDerivI)
                        + penaltyParameter_ * phiNormalI * phiNormalJ * 2;
                    //}else if(std::abs(pPhys[0]-1)<1e-9){//Robin
                    //	ret(j,i)=-phiNormalI*phiNormalJ*0;
                }
                else
                { //homogeneous Neumann
                    result(j, i) = 0;
                }
            }
        }
    }

    ///The vector edition of the face integrand is meant for implementation of the boundary conditions
    //for conforming problems this functions only deals with non-homogeneous Neumann and Robin boundary conditions 
    virtual void faceIntegrand(const FaceT* face, const LinearAlgebra::NumericalVector& normal, const PointReferenceOnTheFaceT& point, LinearAlgebra::NumericalVector& result)
    {
        int n = face->getNrOfBasisFunctions();
        result.resize(n);
        PointPhysicalT pPhys(DIM_);
        face->referenceToPhysical(point, pPhys);
        if (std::abs(pPhys[0]) < 1e-9 || std::abs(pPhys[0] - 1) < 1e-9)
        { //Dirichlet
            LinearAlgebra::NumericalVector phiDeriv(DIM_);
            for (int i = 0; i < n; ++i)
            {
                face->basisFunctionDeriv(i, point, phiDeriv);
                result[i] = (-normal * phiDeriv / Utilities::norm2(normal)
                             + penaltyParameter_ * face->basisFunction(i, point)) * 0.;
            }
            /*}else if(std::abs(pPhys[1]-1)<1e-9){//Neumann and robin
             for(int i=0;i<n;++i){//be careful in 1D; this boundary condition always concerns df\dn
             ret[i]=face->basisFunction(i,p)*-1;
             }*/
        }
        else
        {
            for (int i = 0; i < n; ++i)
            {
                result[i] = 0;
            }
        }
    }

    //hpGEMUISimplified is originally designed for hyperbolic problems
    //those need initial conditions. Laplace and Poisson equations dont need
    //initial conditions, so just provide a dummy implementation
    virtual double initialConditions(const PointPhysicalT& point)
    {
        // initial conditions are not needed for a steady-state problem
        return 0;
    }
    double sourceTerm(const PointPhysicalT& point)
    {
        return sin(2 * M_PI * point[0]) * (4 * M_PI * M_PI) * cos(2 * M_PI * point[1]); //*cos(2*M_PI*p[2])*3;
    }

    ///interpolates the source term
    void elementIntegrand(const ElementT* element, const PointReferenceT& point, LinearAlgebra::NumericalVector& result)
    {
        PointPhysicalT pPhys(DIM_);
        element->referenceToPhysical(point, pPhys);
        result.resize(element->getNrOfBasisFunctions());
        for (int i = 0; i < element->getNrOfBasisFunctions(); ++i)
        {
            result[i] = element->basisFunction(i, point) * sourceTerm(pPhys);
        }
    }

    ///provide information about your solution that you want to use for visualisation
    void writeToTecplotFile(const ElementT* element, const PointReferenceT& point, std::ostream& out)
    {
        LinearAlgebra::NumericalVector value(1);
        element->getSolution(0, point, value);
        out << value[0];
    }

    ///guarantees the linear system keeps the Dirichlet boundary conditions in place. Assumes x already contains expansion coefficients for
    ///the boundaries in question and noise in all other entries
    //this routine in only needed if you use conforming basis functions
    void insertDirichletBoundary(Utilities::GlobalPetscMatrix& A, Utilities::GlobalPetscVector& b, Utilities::GlobalPetscVector& x)
    {
        int numberOfRows(0);
        std::vector<int> rows(0);
        Geometry::PointPhysical pPhys(DIM_);
        Geometry::PointReference centre(DIM_ - 1);
        for (Base::Face* face : meshes_[0]->getFacesList())
        {
            face->getReferenceGeometry()->getCenter(centre);
            face->referenceToPhysical(centre, pPhys);
            if (std::abs(pPhys[0]) < 1e-9 || std::abs(pPhys[0] - 1) < 1e-9)
            {
                A.getMatrixBCEntries(face, numberOfRows, rows);
            }
        }
        int ierr = MatZeroRows(A, numberOfRows, &rows[0], 1.0, x, b);
        CHKERRV(ierr);
    }
    
    bool solve()
    {
        doAllElementIntegration();
        doAllFaceIntegration();

        //create the global matrix using the first element matrix (first zero)
        //and the first face matrix (second zero)
        Utilities::GlobalPetscMatrix A(HpgemUI::meshes_[0], 0, 0);

        //do the same for the global vectors
        Utilities::GlobalPetscVector b(HpgemUI::meshes_[0], 0, 0), x(HpgemUI::meshes_[0]);

        //you explicitly have to tell the global vectors to assemble, because there are also
        //global vectors that have no initial data. An example is x, which will be used to store the solution
        //of the linear problem. No use case is seen for empty global matrixes so they assemble automatically.
        b.assemble();

        //KSP is a PETSc shorthand for [K]rylov [S]ubspace [P]roblem
        //we create one and use it to solve the linear system
        KSP ksp;
        KSPCreate(PETSC_COMM_WORLD, &ksp);
        KSPSetOperators(ksp, A, A);
        KSPSetFromOptions(ksp);
        KSPSetTolerances(ksp, 1e-12, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
        KSPSolve(ksp, b, x);
        KSPConvergedReason conferge;
        KSPGetConvergedReason(ksp, &conferge);
        int iterations;
        KSPGetIterationNumber(ksp, &iterations);
        std::cout << "KSP solver ended because of " << KSPConvergedReasons[conferge] <<
            " in " << iterations << " iterations." << std::endl;

        //once PETSc is done, feed the data back into the elements
        x.writeTimeLevelData(0);

        //so it can be used for post-processing
        std::ofstream outFile("output.dat." + std::to_string(Base::MPIContainer::Instance().getProcessorID()));
        //write tecplot data
        Output::TecplotDiscontinuousSolutionWriter writeFunc(outFile, "test", "01", "value");
        writeFunc.write(meshes_[0], "discontinuous solution", false, this);
        //AND paraview data
        Output::VTKSpecificTimeWriter paraWrite("output", meshes_[0]);
        paraWrite.write([](Base::Element* element, const Geometry::PointReference& point,size_t timelevel)->double
        {
            LinearAlgebra::NumericalVector value(1);
            element->getSolution(timelevel,point,value);
            return value[0];
        }
        ,"value");
        return true;
    }

private:

    //number of elements per cardinal direction
    int n_;

    //polynomial order of the approximation
    int p_;

    //the dimension of your problem
    const unsigned int DIM_;

    //this is a number that pops up in the derivation of the weak formulation when you do DG
    double penaltyParameter_;
};

auto& n = Base::register_argument<int>('n', "numElems", "number of elements per dimension", true);
auto& p = Base::register_argument<int>('p', "order", "polynomial order of the solution", true);
int main(int argc, char **argv)
{
    Base::parse_options(argc, argv);
    try
    {
        //create ...
        Laplace demo(n.getValue(), p.getValue());
        demo.initialise();

        //... and solve the problem
        demo.solve();

        return 0;
    }
    catch (const char* e)
    {
        std::cout << e;
    }
}

