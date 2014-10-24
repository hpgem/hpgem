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
#include <cmath>

///Linear advection equation
///The first self-contained (no PETSc) program to make it into the SVN
///If someone is bored, this should be polished into a demo application and/or a tutorial
///If someone is bored, this would be an excellent basis for a self test that tests all essential features of hpGEM at once

class Advection : public Base::HpgemUISimplified, Output::TecplotSingleElementWriter {
public:

    Advection(int n, int p) : HpgemUISimplified(DIM_, p), n_(n), p_(p) {
    }

    ///set up the mesh

    bool virtual initialise() {

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

    virtual void elementIntegrand(const ElementT* element, const PointReferenceT& point, LinearAlgebra::Matrix& result) {
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

        virtual void elementIntegrand(const Base::Element* element, const Geometry::PointReference& point, LinearAlgebra::Matrix& result) {
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

    virtual void faceIntegrand(const FaceT* face, const LinearAlgebra::NumericalVector& normal, const PointReferenceOnTheFaceT& point, LinearAlgebra::Matrix& result) {
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

    virtual void faceIntegrand(const FaceT* face, const LinearAlgebra::NumericalVector& normal, const PointReferenceOnTheFaceT& point, LinearAlgebra::NumericalVector& result) {
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

    bool solve() {
        //do the integration
        doAllElementIntegration();
        doAllFaceIntegration();

        //manual integration example
        Integration::ElementIntegral elIntegral(false);
        LinearAlgebra::Matrix stifness;
        for (Base::Element* element : meshes_[0]->getElementsList()) {
            int n(element->getNrOfBasisFunctions());
            stifness.resize(n, n);
            elIntegral.integrate(element, &advection, stifness);
            element->setElementMatrix(stifness, 1);
        }

        //prepare tecplot output
        std::ofstream outFile("output.dat");
        Output::TecplotDiscontinuousSolutionWriter out(outFile, "simple advective movement", "01", "u");

        //set the final few variables
        double dt(0.0002), dtplot(0.05), tend(5.), t(0), tplot(t); //always plot the initial data
        LinearAlgebra::Matrix mass, solution, leftResidual, rightResidual;
        bool first(true);

        //finalise interpolation
        for (Base::Element* element : meshes_[0]->getElementsList()) {
            int n(element->getNrOfBasisFunctions());
            mass.resize(n, n);
            element->getElementMatrix(mass, 0);
            LinearAlgebra::NumericalVector initialCondition(n);
            element->getElementVector(initialCondition);
            solution.resize(n, 1);
            for (int i = 0; i < n; ++i) {///\BUG is it useful to have automated conversion from vectors to n by 1 matrices and back?
                solution[i] = initialCondition[i];
            }
            mass.solve(solution);
            solution.resize(1, n);
            element->setTimeLevelData(0, solution);
        }

        //start the time loop
        while (t <= tend + 1e-10) {
            if (t >= tplot - 1e-10) {
                out.write(meshes_[0], "t=" + std::to_string(t), !first, this);
                tplot += dtplot;
            }
            t += dt;

            //construct the RHS
            for (Base::Element* element : meshes_[0]->getElementsList()) {
                //collect data
                int n = element->getNrOfBasisFunctions();
                mass.resize(n, n);
                stifness.resize(n, n);
                leftResidual.resize(n, 1);
                element->getElementMatrix(mass, 0);
                element->getElementMatrix(stifness, 1);
                solution = element->getTimeLevelData(0);
                solution.resize(n, 1);
                //compute residual=M*u+dt*S*u
                leftResidual = mass*solution;
                leftResidual.axpy(dt, stifness * solution);
                leftResidual.resize(1, n);
                element->setResidue(leftResidual);
            }

            //this bit could do with some interface improvements
            for (Base::Face* face : meshes_[0]->getFacesList()) {
                int n(face->getNrOfBasisFunctions()), nLeft(face->getPtrElementLeft()->getNrOfBasisFunctions());
                stifness.resize(n, n);
                face->getFaceMatrix(stifness);
                solution.resize(n, 1);

                //concatenate left and right data
                for (int i = 0; i < n; ++i) {
                    if (i < nLeft) {
                        solution[i] = face->getPtrElementLeft()->getData(0, 0, i);
                    } else {
                        solution[i] = face->getPtrElementRight()->getData(0, 0, i - nLeft);
                    }
                }

                //compute the flux
                solution = stifness*solution;
                solution *= dt;
                leftResidual.resize(1, nLeft);
                rightResidual.resize(1, n - nLeft);

                //unconcatenate left and right data
                for (int i = 0; i < n; ++i) {
                    if (i < nLeft) {
                        leftResidual[i] = solution[i];
                    } else {
                        rightResidual[i - nLeft] = solution[i];
                    }
                }

                //add the flux to the residual
                leftResidual.axpy(1., face->getPtrElementLeft()->getResidue());
                face->getPtrElementLeft()->setResidue(leftResidual);
                if (nLeft < n) {
                    rightResidual.axpy(1., face->getPtrElementRight()->getResidue());
                    face->getPtrElementRight()->setResidue(rightResidual);
                }
            }

            //compute the solution in the next time step based on the residual
            for (Base::Element* element : meshes_[0]->getElementsList()) {
                int n(element->getNrOfBasisFunctions());
                mass.resize(n, n);
                solution.resize(n, 1);
                element->getElementMatrix(mass, 0);
                solution = element->getResidue();
                solution.resize(n, 1);
                mass.solve(solution);
                solution.resize(1, n);
                element->setTimeLevelData(0, solution);
            }
            first = false;
        }

        return true;
    }

private:

    //number of elements per cardinal direction
    int n_;

    //polynomial order of the approximation
    int p_;

    static const unsigned int DIM_;
};

const unsigned int Advection::DIM_(2);

int main(int argc, char **argv) {
    try {
        int n, p;
        if (argc > 2) {
            n = std::atoi(argv[1]);
            p = std::atoi(argv[2]);
        } else {
            throw "usage: LinearAdvection.out n p";
        }
        Advection test(n, p);
        test.initialise();
        test.solve();
        return 0;
    } catch (const char* e) {
        std::cout << e;
    }
}

