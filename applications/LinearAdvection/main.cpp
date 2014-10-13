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
#include <cmath>

///Linear advection equation
///The first self-contained (no PETSc) program to make it into the SVN
///If someone is bored, this should be polished into a demo application
///If someone is bored, this would be an excellent basis for a self test that test all essential features of hpGEM at once

const unsigned int DIM = 2;

class Advection: public Base::HpgemUISimplified, Output::TecplotSingleElementWriter {

public:
	Advection(int n, int p) : HpgemUISimplified(DIM, p), n_(n), p_(p) {
	}

	///set up the mesh
	bool virtual initialise() {
		RectangularMeshDescriptorT description(DIM);
		for (int i = 0; i < DIM; ++i) {
			description.bottomLeft_[i] = 0;
			description.topRight_[i] = 1;
			description.numElementsInDIM_[i] = n_;
			description.boundaryConditions_[i] = RectangularMeshDescriptorT::PERIODIC;
		}
		addMesh(description, Base::TRIANGULAR, 2, 1, 1, 1);
		return true;
	}

	///You pass the reference point to the basisfunctions. Internally the basisfunctions will be mapped to the physical element
	///so you wont have to do any transformations yourself (constructs the mass matrix)
	virtual void elementIntegrand(const ElementT* element, const PointReferenceT& p, LinearAlgebra::Matrix& ret) {
		int n = element->getNrOfBasisFunctions();
		ret.resize(n, n);
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				ret(i, j) = element->basisFunction(i,p) * element->basisFunction(j,p);
			}
		}
	}
        
        ///Alternative way to define integrands; note that you will have to do the integration yourself if you use this way
        class advectiveTerm:public Integration::ElementIntegrandBase<LinearAlgebra::Matrix>{
        public:
        virtual void elementIntegrand(const Base::Element* element, const Geometry::PointReference& p, LinearAlgebra::Matrix& ret){
                int n=element->getNrOfBasisFunctions();
                ret.resize(n,n);
                LinearAlgebra::NumericalVector b(DIM);
                for(double i=0;i<DIM;++i){
                    b[i]=i/10.+0.1;
                }
                LinearAlgebra::NumericalVector phiDerivJ(DIM);
                for(int i=0;i<n;++i){
                    for(int j=0;j<n;++j){
                        element->basisFunctionDeriv(j,p,phiDerivJ);
                        ret(j,i)=element->basisFunction(i,p)*(b*phiDerivJ);
                    }
                }
            }
        }advection;

	///You pass the reference point to the basisfunctions. Internally the basisfunctions will be mapped to the physical element
	///so you wont have to do any transformations yourself. If you expect 4 matrices here, you can assume that ret is block structured such
	///that basisfunctions belonging to the left element are indexed first
	virtual void faceIntegrand(const FaceT* face, const LinearAlgebra::NumericalVector& normal, const PointReferenceOnTheFaceT& p, LinearAlgebra::Matrix& ret) {
		int n = face->getNrOfBasisFunctions(),nLeft=face->getPtrElementLeft()->getNrOfBasisFunctions();
		ret.resize(n, n);
                ret*=0;
                LinearAlgebra::NumericalVector b(DIM),phiNormalJ(DIM);
                for(double i=0;i<DIM;++i){
                    b[i]=i/10.+0.1;
                }
                for(int i=0;i<n;++i){
                    for(int j=0;j<n;++j){
                        double A=(b*normal)/Base::L2Norm(normal);
                        face->basisFunctionNormal(j,normal,p,phiNormalJ);
                        if((A>1e-12)&&(i<nLeft)){
                            ret(j,i)=-(b*phiNormalJ)*face->basisFunction(i,p);
                        }else if ((A<-1e-12)&&(i>=nLeft)){
                            ret(j,i)=-(b*phiNormalJ)*face->basisFunction(i,p);
                        }else if (std::fabs(A)<1e-12){
                            ret(j,i)=-(b*phiNormalJ)*face->basisFunction(i,p)/2.;
                        }
                    }
                }
	}

	///The vector edition of the face integrand is meant for implementation of boundary conditions
        ///This is a periodic problem, so it just returns 0
	virtual void faceIntegrand(const FaceT* face, const LinearAlgebra::NumericalVector& normal, const PointReferenceOnTheFaceT& p, LinearAlgebra::NumericalVector& ret) {
		int n = face->getNrOfBasisFunctions();
		ret.resize(n);
                ret*=0;
	}

	virtual double initialConditions(const PointPhysicalT& p) {
            return std::sin(2*M_PI*p[0])*std::sin(2*M_PI*p[1]);//*std::sin(2*M_PI*p[2]);
            //return p[0]+p[1];//+p[2];
	}

	///interpolates the initial conditions
	void elementIntegrand(const ElementT* element, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret) {
		PointPhysicalT pPhys(DIM);
		element->referenceToPhysical(p, pPhys);
		ret.resize(element->getNrOfBasisFunctions());
		for (int i = 0; i < element->getNrOfBasisFunctions(); ++i) {
			ret[i] = element->basisFunction(i, p) * initialConditions(pPhys);
		}
	}

	void writeToTecplotFile(const ElementT* element, const PointReferenceT& p, std::ostream& out) {
		LinearAlgebra::NumericalVector value(1);
		element->getSolution(0, p, value);
		out << value[0];
	}

	bool solve() {
            //do the integration
            doAllElementIntegration();
            doAllFaceIntegration();
            
            //manual integration example
            Integration::ElementIntegral elIntegral(false);
            LinearAlgebra::Matrix mass,stifness,solution,leftResidual,rightResidual;
            for(Base::Element* element:meshes_[0]->getElementsList()){
                int n(element->getNrOfBasisFunctions());
                stifness.resize(n,n);
                elIntegral.integrate(element,&advection,stifness);
                element->setElementMatrix(stifness,1);
            }
            
            //prepare tecplot output
            std::ofstream outFile("output.dat");
            Output::TecplotDiscontinuousSolutionWriter out(outFile,"simple advective movement","01","u");
            
            //set the final few variables
            double dt(0.0002),dtplot(0.05),tend(5.),t(0),tplot(t);//always plot the initial data
            bool first(true);
            
            for(Base::Element* element:meshes_[0]->getElementsList()){
                int n(element->getNrOfBasisFunctions());
                mass.resize(n,n);
                element->getElementMatrix(mass,0);
                LinearAlgebra::NumericalVector initialCondition(n);
                element->getElementVector(initialCondition);
                solution.resize(n,1);
                for(int i=0;i<n;++i){///\BUG is it useful to have automated conversion from vectors to n by 1 matrices and back?
                    solution[i]=initialCondition[i];
                }
                mass.solve(solution);
                solution.resize(1,n);
                element->setTimeLevelData(0,solution);
            }
            
            //start the time loop
            while(t<=tend+1e-10){
                if(t>=tplot-1e-10){
                    out.write(meshes_[0],"t="+std::to_string(t),!first,this);
                    tplot+=dtplot;
                }
                t+=dt;
                
                //construct the RHS
                for(Base::Element* element:meshes_[0]->getElementsList()){
                    int n=element->getNrOfBasisFunctions();
                    mass.resize(n,n);
                    stifness.resize(n,n);
                    leftResidual.resize(n,1);
                    element->getElementMatrix(mass,0);
                    element->getElementMatrix(stifness,1);
                    solution=element->getTimeLevelData(0);
                    solution.resize(n,1);
                    leftResidual=mass*solution;
                    leftResidual.axpy(dt,stifness*solution);
                    leftResidual.resize(1,n);
                    element->setResidue(leftResidual);
                }
                
                //this bit could do with some interface improvements
                for(Base::Face* face:meshes_[0]->getFacesList()){
                    int n(face->getNrOfBasisFunctions()),nLeft(face->getPtrElementLeft()->getNrOfBasisFunctions());
                    stifness.resize(n,n);
                    face->getFaceMatrix(stifness);
                    solution.resize(n,1);
                    for(int i=0;i<n;++i){//concatenate left and right data
                        if(i<nLeft){
                            solution[i]=face->getPtrElementLeft()->getData(0,0,i);
                        }else{
                            solution[i]=face->getPtrElementRight()->getData(0,0,i-nLeft);
                        }
                    }
                    solution=stifness*solution;
                    solution*=dt;
                    leftResidual.resize(1,nLeft);
                    rightResidual.resize(1,n-nLeft);
                    for(int i=0;i<n;++i){//unconcatenate left and right data
                        if(i<nLeft){
                            leftResidual[i]=solution[i];
                        }else{
                            rightResidual[i-nLeft]=solution[i];
                        }
                    }
                    leftResidual.axpy(1.,face->getPtrElementLeft()->getResidue());
                    face->getPtrElementLeft()->setResidue(leftResidual);
                    if(nLeft<n){
                        rightResidual.axpy(1.,face->getPtrElementRight()->getResidue());
                        face->getPtrElementRight()->setResidue(rightResidual);
                    }
                }
                
                //solve the system
                for(Base::Element* element:meshes_[0]->getElementsList()){
                    int n(element->getNrOfBasisFunctions());
                    mass.resize(n,n);
                    solution.resize(n,1);
                    element->getElementMatrix(mass,0);
                    solution=element->getResidue();
                    solution.resize(n,1);
                    mass.solve(solution);
                    solution.resize(1,n);
                    element->setTimeLevelData(0,solution);
                }
                first=false;
            }
            
            return true;
	}

private:

	//number of elements per cardinal direction
	int n_;

	//polynomial order of the approximation
	int p_;
};

int main(int argc, char **argv) {
	try {
		int n, p;
		if (argc > 2) {
			n = std::atoi(argv[1]);
			p = std::atoi(argv[2]);
			argv[2] = argv[0];
			argc -= 2;
			argv += 2;
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

