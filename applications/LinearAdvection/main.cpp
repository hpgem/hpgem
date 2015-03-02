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
#include "Utilities/BasisFunctions2DH1ConformingSquare.hpp"
#include <cmath>



///Linear advection equation
///The first self-contained (no PETSc) program to make it into the SVN
///If someone is bored, this should be polished into a demo application and/or a tutorial
///If someone is bored, this would be an excellent basis for a self test that tests all essential features of hpGEM at once

class Advection : public Base::HpgemUISimplified {
public:

    Advection(int n, int p) : HpgemUISimplified(DIM_, p), n_(n), p_(p) {
    }

    ///set up the mesh

    bool initialise() override {

        //describes a rectangular domain
        RectangularMeshDescriptorT description(DIM_);

        //this demo will use a cube
        for (int i = 0; i < DIM_; ++i) {
            description.bottomLeft_[i] = 0;
            description.topRight_[i] = 1;
            description.numElementsInDIM_[i] = n_;

            //at the moment your options are SOLID_WALL and PERIODIC
            //once we decide what names of boundary conditions to support
            //it will become possible to appropriate boundary conditions
            //for your problem here
            description.boundaryConditions_[i] = RectangularMeshDescriptorT::PERIODIC;
        }

        //create a triangular mesh. The four magic ones that are passed to this function
        //specify the number of element matrices, the number of element vectors,
        //the number of face matrices and the number of face vectors (in that order)
        addMesh(description, Base::TRIANGULAR, 2, 1, 1, 1);

        //tell hpGEM to use basis functions that are discontinuous and are designed for triangles
        //this is likely to get automated by hpGEM at some point in the future
        meshes_[0]->setDefaultBasisFunctionSet(Utilities::createDGBasisFunctionSet2DH1Triangle(p_));
        return true;
    }

    ///You pass the reference point to the basisfunctions. Internally the basisfunctions will be mapped to the physical element
    ///so you wont have to do any transformations yourself (constructs the mass matrix)

    void elementIntegrand(const ElementT* element, const PointReferenceT& point, LinearAlgebra::Matrix& result) override {
        int n = element->getNrOfBasisFunctions();
        result.resize(n, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                result(i, j) = element->basisFunction(i, point) * element->basisFunction(j, point);
            }
        }
    }

    ///Alternative way to define integrands; note that you will have to do the integration yourself if you use this way

    class advectiveTerm : public Integration::ElementIntegrandBase<LinearAlgebra::Matrix> {
    public:

        void elementIntegrand(const Base::Element* element, const Geometry::PointReference& point, LinearAlgebra::Matrix& result) override {
            int n = element->getNrOfBasisFunctions();
            result.resize(n, n);
            LinearAlgebra::NumericalVector b(DIM_);

            //choose the direction of advection
            for (double i = 0; i < DIM_; ++i) {
                b[i] = i / 10. + 0.1;
            }
            LinearAlgebra::NumericalVector phiDerivJ(DIM_);
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    element->basisFunctionDeriv(j, point, phiDerivJ);
                    result(j, i) = element->basisFunction(i, point)*(b * phiDerivJ);
                }
            }
        }
    } advection;

    ///You pass the reference point to the basisfunctions. Internally the basisfunctions will be mapped to the physical element
    ///so you wont have to do any transformations yourself. If you expect 4 matrices here, you can assume that ret is block structured such
    ///that basisfunctions belonging to the left element are indexed first

    void faceIntegrand(const FaceT* face, const LinearAlgebra::NumericalVector& normal, const PointReferenceT& point, LinearAlgebra::Matrix& result) override {
        int n = face->getNrOfBasisFunctions(), nLeft = face->getPtrElementLeft()->getNrOfBasisFunctions();
        result.resize(n, n);
        result *= 0;
        LinearAlgebra::NumericalVector b(DIM_), phiNormalJ(DIM_);

        //choose the direction of advection
        for (double i = 0; i < DIM_; ++i) {
            b[i] = i / 10. + 0.1;
        }
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                double A = (b * normal) / Base::L2Norm(normal);
                face->basisFunctionNormal(j, normal, point, phiNormalJ);

                //upwind flux
                if ((A > 1e-12)&&(i < nLeft)) {
                    result(j, i) = -(b * phiNormalJ) * face->basisFunction(i, point);
                } else if ((A<-1e-12)&&(i >= nLeft)) {
                    result(j, i) = -(b * phiNormalJ) * face->basisFunction(i, point);
                } else if (std::abs(A) < 1e-12) {
                    result(j, i) = -(b * phiNormalJ) * face->basisFunction(i, point) / 2.;
                }
            }
        }
    }

    ///The vector edition of the face integrand is meant for implementation of boundary conditions
    ///This is a periodic problem, so it just return 0

    void faceIntegrand(const FaceT* face, const LinearAlgebra::NumericalVector& normal, const PointReferenceT& point, LinearAlgebra::NumericalVector& result) override {
        int n = face->getNrOfBasisFunctions();
        result.resize(n);
        result *= 0;
    }

    virtual double initialConditions(const PointPhysicalT& point) {
        return std::sin(2 * M_PI * point[0]) * std::sin(2 * M_PI * point[1]); //*std::sin(2*M_PI*p[2]);
        //return p[0]+p[1];//+p[2];
    }

    ///interpolates the initial conditions

    void elementIntegrand(const ElementT* element, const PointReferenceT& point, LinearAlgebra::NumericalVector& result) {
        PointPhysicalT pPhys(DIM_);
        element->referenceToPhysical(point, pPhys);
        result.resize(element->getNrOfBasisFunctions());
        for (int i = 0; i < element->getNrOfBasisFunctions(); ++i) {
            result[i] = element->basisFunction(i, point) * initialConditions(pPhys);
        }
    }

    ///provide information about your solution that you want to use for visualisation

    void writeToTecplotFile(const ElementT* element, const PointReferenceT& point, std::ostream& out) {
        LinearAlgebra::NumericalVector value(1);
        element->getSolution(0, point, value);
        out << value[0];
    }
    
    ///TODO this cannot be automated because I dont know where the mass matrix is
    // solve Mx=`residue`
    void interpolate() override {
        LinearAlgebra::Matrix mass;
        LinearAlgebra::NumericalVector solution;
        for (Base::Element* element : meshes_[0]->getElementsList()) {
            int n(element->getNrOfBasisFunctions());
            mass.resize(n, n);
            mass = element->getElementMatrix(0);
            solution = element->getResidue();
            solution.resize(n);
            mass.solve(solution);
            element->setTimeLevelData(0, solution);
        }
    }


    void computeRhsLocal() override
    {
        LinearAlgebra::Matrix mass, stiffness;
        LinearAlgebra::NumericalVector oldData, residual;
        for (Base::Element* element : meshes_[0]->getElementsList()) {
            //collect data
            int n = element->getNrOfBasisFunctions();
            mass.resize(n, n);
            stiffness.resize(n, n);
            mass = element->getElementMatrix(0);
            stiffness = element->getElementMatrix(1);
            oldData = element->getTimeLevelData(0);
            //oldData.resize(n);
            //compute residual=M*u+dt*S*u
            residual = mass*oldData;
            residual.axpy(dt_, stiffness * oldData);
            element->setResidue(residual);
        }
    }


    void computeRhsFaces() override
    {
        LinearAlgebra::Matrix stiffness;
        LinearAlgebra::NumericalVector residue;
        for (Base::Face* face : meshes_[0]->getFacesList()) {
            int n(face->getNrOfBasisFunctions()), nLeft(face->getPtrElementLeft()->getNrOfBasisFunctions());
            stiffness.resize(n, n);
            face->getFaceMatrix(stiffness);
            residue.resize(n);

            ///TODO implement face->getData()
            //for now concatenate left and right data
            for (int i = 0; i < n; ++i) {
                if (i < nLeft) {
                    residue[i] = face->getPtrElementLeft()->getData(0, 0, i);
                } else {
                    residue[i] = face->getPtrElementRight()->getData(0, 0, i - nLeft);
                }
            }

            //compute the flux
            residue = stiffness*residue;
            residue *= dt_;
            face->setResidue(residue);
        }
    }


    void beforeTimeIntegration() override {
        
        //manual integration example
        Integration::ElementIntegral elIntegral(false);
        LinearAlgebra::Matrix stifness;
        for (Base::Element* element : meshes_[0]->getElementsList()) {
            int n(element->getNrOfBasisFunctions());
            stifness = elIntegral.integrate(element, &advection);
            element->setElementMatrix(stifness, 1);
        }
    }

private:

    //number of elements per cardinal direction
    int n_;

    //polynomial order of the approximation
    int p_;

    static const unsigned int DIM_;
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
        test.registerVTKWriteFunction([](Base::Element* element, const Geometry::PointReference& point, std::size_t timelevel)->double
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
        std::cout << e << std::endl;
        std::exit(3);
    }
}

