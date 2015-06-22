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

#include "Base/MpiContainer.h"
#include "Base/HpgemUISimplified.h"
#include "Base/L2Norm.h"
#include "Utilities/GlobalMatrix.h"
#include "Utilities/GlobalVector.h"
#include "petscksp.h"
#include "Output/TecplotSingleElementWriter.h"
#include <fstream>
#include "Geometry/ReferenceTetrahedron.h"
#include "Utilities/BasisFunctions1DH1ConformingLine.h"
#include "Utilities/BasisFunctions2DH1ConformingSquare.h"
#include "Utilities/BasisFunctions2DH1ConformingTriangle.h"
#include "Utilities/BasisFunctions3DH1ConformingCube.h"
#include "Utilities/BasisFunctions3DH1ConformingTetrahedron.h"
#include "Base/ElementCacheData.h"
#include "Base/FaceCacheData.h"
#include "Output/TecplotDiscontinuousSolutionWriter.h"
#include "Base/Element.h"
#include "Base/RectangularMeshDescriptor.h"
#include "Integration/ElementIntegral.h"
#include "Base/CommandLineOptions.h"
#include "Output/VTKSpecificTimeWriter.h"
#include <chrono>

class Laplace : public Base::HpgemUISimplified
{
public:
    ///constructor: set the dimension of the problem, start the API of hpGEM and initialise the other fields
    //initialisation order is not fixed so DIM_ has to be passed to hpGEMUISimplified in a hardcoded way
    Laplace(int n, int p)
            : HpgemUISimplified(3, p), n_(n), p_(p), DIM_(3)
    {
        penaltyParameter_ = 3 * n_ * p_ * (p_ + DIM_ - 1) + 1;
    }
    
    ///set up the mesh
    bool initialise() override final
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
            description.boundaryConditions_[i] = Base::BoundaryType::SOLID_WALL;
        }
        
        //create a triangular mesh. The four magic ones that are passed to this function
        //specify the number of element matrices, the number of element vectors,
        //the number of face matrices and the number of face vectors (in that order)
        addMesh(description, Base::MeshType::TRIANGULAR, 1, 1, 1, 1);
        
        //tell hpGEM to use basis functions that are discontinuous and are designed for triangles
        //this is likely to get automated by hpGEM at some point in the future
        meshes_[0]->setDefaultBasisFunctionSet(Utilities::createDGBasisFunctionSet3DH1Tetrahedron(p_));
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
        bFsets = Utilities::createVertexBasisFunctionSet2DH1Triangle(p_);
        meshes_[0]->addVertexBasisFunctionSet(bFsets);
        std::vector<const Base::OrientedBasisFunctionSet*> oBFsets;
        oBFsets = Utilities::createFaceBasisFunctionSet2DH1Triangle(p_);
        meshes_[0]->addFaceBasisFunctionSet(oBFsets);
    }
    
    ///You pass the reference point to the basis functions. Internally the basis functions will be mapped to the physical element
    ///so you wont have to do any transformations yourself
    void elementIntegrand(const ElementT* element, const PointReferenceT& point, LinearAlgebra::Matrix& result) override final
    {
        //be careful that n changes meaning inside integration routines to represent the number of basis functions involved
        int n = element->getNrOfBasisFunctions();
        LinearAlgebra::NumericalVector phiDerivI(DIM_), phiDerivJ(DIM_);
        result.resize(n, n);
        for (int i = 0; i < n; ++i)
        {
            phiDerivI = element->basisFunctionDeriv(i, point);
            for (int j = 0; j < n; ++j)
            {
                phiDerivJ = element->basisFunctionDeriv(j, point);
                
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
    void faceIntegrand(const FaceT* face, const LinearAlgebra::NumericalVector& normal, const PointReferenceT& point, LinearAlgebra::Matrix& result) override final
    {
        int n = face->getNrOfBasisFunctions();
        result.resize(n, n);
        LinearAlgebra::NumericalVector phiNormalI(DIM_), phiNormalJ(DIM_), phiDerivI(DIM_), phiDerivJ(DIM_);
        PointPhysicalT pPhys = face->referenceToPhysical(point);
        for (int i = 0; i < n; ++i)
        {
            phiNormalI = face->basisFunctionNormal(i, normal, point);
            phiDerivI = face->basisFunctionDeriv(i, point);
            for (int j = 0; j < n; ++j)
            {
                phiNormalJ = face->basisFunctionNormal(j, normal, point);
                phiDerivJ = face->basisFunctionDeriv(j, point);
                if (face->isInternal())
                {
                    result(j, i) = -(phiNormalI * phiDerivJ + phiNormalJ * phiDerivI) / 2 + penaltyParameter_ * phiNormalI * phiNormalJ;
                    
                    //for the moment you can figure out at what part of the boundary you are by looking at
                    //the point in physical space you are (pPhys) and implement a boundary condition accordingly
                    //once there are more labels for boundary conditions available you can use syntax that looks like
                    //face->getFaceType() to see what kind of boundary condition you are dealing with
                }
                else if (std::abs(pPhys[0]) < 1e-9 || std::abs(pPhys[0] - 1.) < 1e-9) //Dirichlet
                {
                    result(j, i) = -(phiNormalI * phiDerivJ + phiNormalJ * phiDerivI) + penaltyParameter_ * phiNormalI * phiNormalJ * 2;
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
    void faceIntegrand(const FaceT* face, const LinearAlgebra::NumericalVector& normal, const PointReferenceT& point, LinearAlgebra::NumericalVector& result) override final
    {
        int n = face->getNrOfBasisFunctions();
        result.resize(n);
        PointPhysicalT pPhys = face->referenceToPhysical(point);
        if (std::abs(pPhys[0]) < 1e-9 || std::abs(pPhys[0] - 1) < 1e-9)
        { //Dirichlet
            LinearAlgebra::NumericalVector phiDeriv(DIM_);
            for (int i = 0; i < n; ++i)
            {
                phiDeriv = face->basisFunctionDeriv(i, point);
                result[i] = (-normal * phiDeriv / Base::L2Norm(normal) + penaltyParameter_ * face->basisFunction(i, point)) * 0.;
            }
            /*}else if(std::abs(pPhys[1]-1)<1e-9){//Neumann and robin
             for(int i=0;i<n;++i){//be careful in 1D; this boundary condition always concerns df/dn
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
    double initialConditions(const PointPhysicalT& point) override final
    {
        // initial conditions are not needed for a steady-state problem
        return 0;
    }
    double sourceTerm(const PointPhysicalT& point)
    {
        return std::sin(2 * M_PI * point[0]) * (4 * M_PI * M_PI) * std::cos(2 * M_PI * point[1]) * std::cos(2 * M_PI * point[2]) * 3;
    }
    
    ///interpolates the source term
    void elementIntegrand(const ElementT* element, const PointReferenceT& point, LinearAlgebra::NumericalVector& result) override final
    {
        PointPhysicalT pPhys = element->referenceToPhysical(point);
        result.resize(element->getNrOfBasisFunctions());
        for (int i = 0; i < element->getNrOfBasisFunctions(); ++i)
        {
            result[i] = element->basisFunction(i, point) * sourceTerm(pPhys);
        }
    }
    
    ///provide information about your solution that you want to use for visualisation
    void writeToTecplotFile(const ElementT* element, const PointReferenceT& point, std::ostream& out) override final
    {
        LinearAlgebra::NumericalVector value(1);
        value = element->getSolution(0, point);
        out << value[0];
    }
    
    ///guarantees the linear system keeps the Dirichlet boundary conditions in place. Assumes x already contains expansion coefficients for
    ///the boundaries in question and noise in all other entries
    //this routine in only needed if you use conforming basis functions
    void insertDirichletBoundary(Utilities::GlobalPetscMatrix& A, Utilities::GlobalPetscVector& b, Utilities::GlobalPetscVector& x)
    {
        std::size_t numberOfRows(0);
        std::vector<int> rows(0);
        Geometry::PointPhysical pPhys(DIM_);
        Geometry::PointReference centre(DIM_ - 1);
        for (Base::Face* face : meshes_[0]->getFacesList())
        {
            centre = face->getReferenceGeometry()->getCenter();
            pPhys = face->referenceToPhysical(centre);
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
        Utilities::GlobalPetscMatrix A(HpgemAPIBase::meshes_[0], 0, 0);
        
        //do the same for the global vectors
        Utilities::GlobalPetscVector b(HpgemAPIBase::meshes_[0], 0, 0), x(HpgemAPIBase::meshes_[0]);
        
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
        std::cout << "KSP solver ended because of " << KSPConvergedReasons[conferge] << " in " << iterations << " iterations." << std::endl;
        
        //once PETSc is done, feed the data back into the elements
        x.writeTimeLevelData(0);
        
        //so it can be used for post-processing
        std::ofstream outFile("output." + std::to_string(Base::MPIContainer::Instance().getProcessorID()) + ".dat");
        //write tecplot data
        Output::TecplotDiscontinuousSolutionWriter writeFunc(outFile, "test", "012", "value");
        writeFunc.write(meshes_[0], "discontinuous solution", false, this);
        //AND paraview data
        Output::VTKSpecificTimeWriter paraWrite("output", meshes_[0]);
        paraWrite.write([](Base::Element* element, const Geometry::PointReference& point, std::size_t timelevel)->double
        {   
            LinearAlgebra::NumericalVector value(1);
            value = element->getSolution(timelevel,point);
            return value[0];
        }, "value");
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

auto& numBasisFuns = Base::register_argument<int>('n', "numElems", "number of elements per dimension", true);
auto& p = Base::register_argument<int>('p', "order", "polynomial order of the solution", true);
int main(int argc, char **argv)
{
    Base::parse_options(argc, argv);
    try
    {
        auto start = std::chrono::high_resolution_clock::now();
        
        //create ...
        Laplace demo(numBasisFuns.getValue(), p.getValue());
        demo.initialise();
        //... and solve the problem
        demo.solve();
        
        auto end = std::chrono::high_resolution_clock::now();
        
        std::cout << "this simulation took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
        
        return 0;
    }
    catch (const char* e)
    {
        std::cout << e;
    }
}

