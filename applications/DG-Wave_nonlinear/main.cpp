/*
 * main.cpp
 *
 *  Created on: Feb 28, 2014
 *      Author: brinkf
 */

#include "Base/HpgemUI.hpp"
#include "Base/Norm2.hpp"
#include "Output/TecplotSingleElementWriter.hpp"
#include <fstream>
#include "Utilities/BasisFunctions2DH1ConformingSquare.hpp"
#include "Integration/ElementIntegrandBase.hpp"
#include "Integration/FaceIntegrandBase.hpp"
#include "Integration/FaceIntegral.hpp"
#include "Integration/ElementIntegral.hpp"
#include "Base/GlobalData.hpp"
#include "Base/RectangularMeshDescriptor.hpp"
#include "Base/ConfigurationData.hpp"
#include "Base/Element.hpp"
#include "Base/Face.hpp"
#include "Output/TecplotDiscontinuousSolutionWriter.hpp"
#include "Base/ElementCacheData.hpp"
#include "Base/FaceCacheData.hpp"
#include "Geometry/PhysicalGeometry.hpp"
#include "Base/Edge.hpp"
#include "Geometry/Mappings/MappingReferenceToPhysical.hpp"
#include "Integration/QuadratureRules/AllGaussQuadratureRules.hpp"
#include "Geometry/ReferenceLine.hpp"
#include "Utilities/BasisFunctions3DH1ConformingCube.hpp"
#include <cmath>

//include headers that depend on c libraries last to minimize the risk of conflict...
#include "Utilities/GlobalMatrix.hpp"
#include "Utilities/GlobalVector.hpp"
#include "petscksp.h"//WARNING: secretly includes math.h and some other c libraries...

const unsigned int DIM = 2;
//const double L=4.963586542006828;
const double L = 12.;

class Fentonwave {
public:

    Fentonwave() {
        std::ifstream rfdata("rfwavedata.dat");
        rfdata.precision(14);
        rfdata>>n_;
        eta_ = new long double[n_ + 1];
        B_ = new long double[n_ + 1];
        for (int i = 0; i < n_ + 1; ++i) {
            rfdata >> eta_[i];
        }
        for (int i = 0; i < n_ + 1; ++i) {
            rfdata >> B_[i];
        }
        rfdata>>c_;
        rfdata>>k_;
        rfdata>>Q_;
        rfdata>>R_;
        rfdata.close();
    }

    void operator()(const double t, const Geometry::PointPhysical& p, LinearAlgebra::NumericalVector& ret) {
        ret[0] = 1;
        long double eps = 100.;
        long double x = p[0] - c_*t;
        while (eps > 1e-12) {
            long double u(B_[0]), v(0), du(0), dv(0);
            for (int i = 1; i < n_ + 1; ++i) {
                u += k_ * i * B_[i] * std::cos(i * k_ * x) * std::cosh(i * k_ * ret[0]) / std::cosh(i * k_);
                v += k_ * i * B_[i] * std::sin(i * k_ * x) * std::sinh(i * k_ * ret[0]) / std::cosh(i * k_);
                du += k_ * k_ * i * i * B_[i] * std::cos(i * k_ * x) * std::sinh(i * k_ * ret[0]) / std::cosh(i * k_);
                dv += k_ * k_ * i * i * B_[i] * std::sin(i * k_ * x) * std::cosh(i * k_ * ret[0]) / std::cosh(i * k_);
            }
            long double f = 0.5 * (u * u + v * v) + ret[0] - R_;
            long double df = (u * du + v * dv) + 1;
            ret[0] -= f / df;
            eps = std::abs(f);
        }
        long double cE = R_ - 0.5 * c_*c_;
        ret[1] = B_[0] * x + c_ * x - cE*t;
        for (int i = 1; i < n_ + 1; ++i) {
            ret[1] += B_[i] * std::sin(i * k_ * x) * std::cosh(i * k_ * (p[DIM - 1] + 1)) / std::cosh(i * k_);
        }
        ret[0] -= 1;
    }

    ~Fentonwave() {
        delete[] eta_;
        delete[] B_;
    }

private:
    int n_;
    long double *eta_;
    long double *B_;
    long double c_;
    long double k_;
    long double Q_;
    long double R_;
};

class waveMakerFileReader {
public:

    waveMakerFileReader() {
        std::ifstream position("PistonMotion.dat");
        std::ifstream derivative("PistonVelocity.dat");
        std::ifstream secondDerivative("PistonAcceleration.dat");
        double time, s, u, temp;
        while (!position.eof()) {
            position >> time;
            position>>s;
            assert(!derivative.eof());
            derivative>>temp;
            derivative>>u;
            assert(time == temp);
            t.push_back(time);
            x.push_back(s);
            dx.push_back(u);
            assert(!secondDerivative.eof());
            secondDerivative>>temp;
            secondDerivative>>u;
            assert(time == temp);
            ddx.push_back(u);
        }
        //for(int i=0;i<t.size();++i){
        //    t[i]=i/50.;
        //}
        //for(int i=1;i<x.size();++i){
        //    dx[i]=dx[i-1]+(ddx[i-1]+ddx[i])/2*(t[i]-t[i-1]);
        //    x[i]=x[i-1]+(dx[i-1]+dx[i])/2*(t[i]-t[i-1]);
        //}
    }

    std::pair<double, double> operator()(double time) {
        static bool output(true);
        static int i(0);
        //while(!((t[i]<=time)&&(t[i+1]>time))){//I assume the vectors are iterated through in a reasonably ordered manner and the time vector is sorted in a non-decreasing way
        //    output=true;
        //    if(t[i]>time)i--;else i++;
        //}
        i = time * 50.; //deliberate reproduction of a bug in the Gagarina code
        std::pair<double, double>ret;
        ret.first = (time - t[i]) * x[i + 1]+(t[i + 1] - time) * x[i];
        ret.second = (time - t[i]) * dx[i + 1]+(t[i + 1] - time) * dx[i];
        //ret.second=x[i+1]-x[i];
        if (output == true) {
            output = false;
            //std::cout<<"t: "<<t[i]<<" derivatives: "<<ret.second<<" "<<(x[i+1]-x[i])<<std::endl;
        }
        ret.first /= t[i + 1] - t[i];
        ret.second /= t[i + 1] - t[i];
        return ret;
    }

private:

    std::vector<double> t, x, dx, ddx;
} piston;

class DGWave : public Base::HpgemUI, public Output::TecplotSingleElementWriter {
public:

    DGWave(int n, int p) : HpgemUI(new Base::GlobalData(), new Base::ConfigurationData(DIM, 2, p, 1)), p_(p) {
        n_ = n;
    }

    ///set up the mesh

    bool virtual initialise() {
        RectangularMeshDescriptorT description(DIM);
        for (int i = 0; i < DIM; ++i) {
            description.bottomLeft_[i] = 0;
            description.topRight_[i] = L;
            description.numElementsInDIM_[i] = n_;
            description.boundaryConditions_[i] = RectangularMeshDescriptorT::PERIODIC;
        }
        description.numElementsInDIM_[1] = 4;
        description.bottomLeft_[DIM - 1] = 0.6;
        description.topRight_[DIM - 1] = 0;
        //description.numElementsInDIM_[DIM - 1] = n_ / 8;
        description.numElementsInDIM_[DIM - 1] = 10;
        description.boundaryConditions_[DIM - 1] = RectangularMeshDescriptorT::SOLID_WALL;
        description.boundaryConditions_[0] = RectangularMeshDescriptorT::SOLID_WALL;
        addMesh(description, Base::RECTANGULAR, 2, 0, 1, 1);
        addMesh(description, Base::RECTANGULAR, 0, 0, 3, 0); //construct surface matrices that dont have to account for surface tilt
        addMesh(description, Base::RECTANGULAR, 0, 0, 0, 0); //keeps original node coordinates
        meshes_[0]->setDefaultBasisFunctionSet(Utilities::createInteriorBasisFunctionSet2DH1Square(p_));
        meshes_[1]->setDefaultBasisFunctionSet(Utilities::createInteriorBasisFunctionSet2DH1Square(p_));
        std::vector<const Base::BasisFunctionSet*> bFsets;
        Utilities::createVertexBasisFunctionSet2DH1Square(p_, bFsets);
        meshes_[0]->addVertexBasisFunctionSet(bFsets);
        meshes_[1]->addVertexBasisFunctionSet(bFsets);
        std::vector<const Base::OrientedBasisFunctionSet*> oBFsets;
        Utilities::createFaceBasisFunctionSet2DH1Square(p_, oBFsets);
        meshes_[0]->addFaceBasisFunctionSet(oBFsets);
        meshes_[1]->addFaceBasisFunctionSet(oBFsets);
        //oBFsets.clear();
        //Utilities::createEdgeBasisFunctionSet3DH1Cube(p_,oBFsets);
        //meshes_[0]->addEdgeBasisFunctionSet(oBFsets);

        const std::vector<Geometry::PointPhysical>& nodes = meshes_[0]->getNodes();

        for (int i = 0; i < meshes_[0]->getNumberOfNodes(); ++i) {
            if (std::fabs(nodes[i][DIM - 1]) < 1e-9) {
                nodesBelowASurfaceNode_[i].push_back(const_cast<Geometry::PointPhysical*> (&nodes[i]));
            }
        }

        for (std::pair<const int, std::vector<Geometry::PointPhysical*> >& node : nodesBelowASurfaceNode_) {
            for (int i = 0; i < meshes_[0]->getNumberOfNodes(); ++i) {
                bool same = true;
                for (int j = 0; j < DIM - 1; ++j) {
                    if (nodes[i][j] != (*node.second[0])[j]) {
                        same = false;
                    }
                }
                if (same) {
                    surfaceNodesAboveANode_[i] = node.first;
                    if (i != node.first) {
                        node.second.push_back(const_cast<Geometry::PointPhysical*> (&nodes[i]));
                    }
                }
            }
        }

        for (Base::Face* face : meshes_[0]->getFacesList()) {
            Geometry::PointReference p(DIM - 1);
            Geometry::PointPhysical pPhys(DIM);
            face->getReferenceGeometry()->getCenter(p);
            face->referenceToPhysical(p, pPhys);
            if (std::fabs(pPhys[DIM - 1]) < 1e-9) {
                face->setFaceType(Geometry::OPEN_BC);
            }
        }
        for (Base::Face* face : meshes_[1]->getFacesList()) {
            Geometry::PointReference p(DIM - 1);
            Geometry::PointPhysical pPhys(DIM);
            face->getReferenceGeometry()->getCenter(p);
            face->referenceToPhysical(p, pPhys);
            if (std::fabs(pPhys[DIM - 1]) < 1e-9) {
                face->setFaceType(Geometry::OPEN_BC);
            }
        }

        return true;
    }

    class anonymous1 : public Integration::ElementIntegrandBase<LinearAlgebra::Matrix> {
    public:

        void elementIntegrand(const Base::Element* element, const Geometry::PointReference& p, LinearAlgebra::Matrix& ret) {
            int n = element->getNrOfBasisFunctions();
            ret.resize(n, n);
            LinearAlgebra::NumericalVector phiDerivI(DIM), phiDerivJ(DIM);
            for (int i = 0; i < n; ++i) {
                //element->getPhysicalGeometry()->getNodeCoordinates(i,pPhys);
                element->basisFunctionDeriv(i, p, phiDerivI);
                for (int j = 0; j <= i; ++j) {
                    element->basisFunctionDeriv(j, p, phiDerivJ);
                    ret(i, j) = ret(j, i) = phiDerivI*phiDerivJ;
                }
            }
        }
    } stifnessIntegrand;

    class anonymous2 : public Integration::FaceIntegrandBase<LinearAlgebra::Matrix> {
    public:

        void faceIntegrand(const Base::Face* face, const LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p, LinearAlgebra::Matrix& ret) {
            int n = face->getNrOfBasisFunctions();
            ret.resize(n, n);
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {//the basis functions belonging to internal parameters are 0 on the free surface anyway.
                    ret(i, j) = face->basisFunction(i, p) * face->basisFunction(j, p);
                }
            }
        }
    } massIntegrand;

    static void initialConditions(const PointPhysicalT& p, LinearAlgebra::NumericalVector& ret) {
        ret.resize(2);
        //double k=2.*M_PI;
        ret[0] = cos(M_PI/6.*p[0])*0; //wave maker
        ret[1] = 0;
        //double k=2.*M_PI;
        //ret[0]=cosh(k*(0+1))*sin(-k*p[0])*sqrt(k*tanh(k))/cosh(k)*0.03;//moving wave
        //ret[1]=cosh(k*(p[DIM-1]+1))*cos(-k*p[0])/cosh(k)*0.03;
        //exactSolution(0,p,ret);//fenton waves or generic case
    }

    static void exactSolution(const double t, const PointPhysicalT& p, LinearAlgebra::NumericalVector& ret) {
        ret.resize(2);
        //static Fentonwave data;
        //data(t,p,ret);
        //double k=2.*M_PI;
        ret[0] = 0;
        ret[1] = 0;
        //ret[0]=(cosh(k*(1))*sin(sqrt(k*tanh(k))*t-k*p[0])/cosh(k))*0.00001+1;//moving wave
        //ret[1]=(cosh(k*(1))*cos(sqrt(k*tanh(k))*t-k*p[0])/cosh(k))*0.00001/sqrt(k*tanh(k));  
    }

    static double waveMakerPosition() {
        return piston(t).first;
        //double omega = std::sqrt(M_PI/2*std::tanh(M_PI/2.));
        //return (1. - cos(omega * t)) * 6.;
        //return 0;
    }

    static double waveMakerDerivative() {
        return piston(t).second;
        //double omega = std::sqrt(M_PI/2*std::tanh(M_PI/2.));
        //return sin(omega * t) * omega * 6.;
        //return 0;
    }

    class anonymous3 : public Integration::FaceIntegrandBase<LinearAlgebra::Matrix> {
    public:

        void faceIntegrand(const FaceT* element, const LinearAlgebra::NumericalVector&, const PointReferenceT& p, LinearAlgebra::Matrix& ret) {
            PointPhysicalT pPhys(DIM);
            element->referenceToPhysical(p, pPhys);
            ret.resize(2, element->getNrOfBasisFunctions());
            LinearAlgebra::NumericalVector data;
            initialConditions(pPhys, data);
            for (int i = 0; i < element->getNrOfBasisFunctions(); ++i) {
                ret(0, i) = element->basisFunction(i, p) * data[0];
                ret(1, i) = element->basisFunction(i, p) * data[1];
            }
        }
    } interpolator;

    class anonymous4 : public Integration::FaceIntegrandBase<LinearAlgebra::NumericalVector> {
    public:

        virtual void faceIntegrand(const FaceT* face, const LinearAlgebra::NumericalVector& normal, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret) {
            PointPhysicalT pPhys(DIM);
            PointReferenceT pElement(DIM);
            face->mapRefFaceToRefElemL(p, pElement);
            face->getPtrElementLeft()->referenceToPhysical(pElement, pPhys);
            ret.resize(2);
            static LinearAlgebra::NumericalVector exact(2);
            if (face->getFaceType() == Geometry::OPEN_BC) {
                face->getPtrElementLeft()->getSolution(0, pElement, ret);
                exactSolution(t, pPhys, exact);
                ret -= exact;
                ret[0] *= ret[0];
                ret[1] *= ret[1];
            }
        }
    } error;

    class anonymous5 : public Integration::FaceIntegrandBase<LinearAlgebra::NumericalVector> {
    public:

        virtual void faceIntegrand(const FaceT* face, const LinearAlgebra::NumericalVector& normal, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret) {
            ret.resize(1);
            static LinearAlgebra::NumericalVector dummySolution(2);
            const LinearAlgebra::Matrix& expansioncoefficients = face->getPtrElementLeft()->getTimeLevelData(0);
            if (face->getFaceType() == Geometry::OPEN_BC) {
                dummySolution[0] = 0;
                for (int i = 0; i < face->getNrOfBasisFunctions(); ++i) {
                    dummySolution[0] += face->basisFunction(i, p) * expansioncoefficients(0, i);
                }
                ret[0] = dummySolution[0] * dummySolution[0]*9.8;
                ret[0] /= 2; //assumes g=9.8
            }
        }
    } faceEnergy;

    class anonymous6 : public Integration::ElementIntegrandBase<LinearAlgebra::NumericalVector> {
    public:

        void elementIntegrand(const Base::Element* element, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret) {
            int n = element->getNrOfBasisFunctions();
            ret.resize(1);
            LinearAlgebra::NumericalVector gradPhi(DIM), temp(DIM);
            const LinearAlgebra::Matrix& expansioncoefficients = element->getTimeLevelData(0);
            for (int i = 0; i < n; ++i) {
                element->basisFunctionDeriv(i, p, temp);
                temp *= expansioncoefficients(1, i);
                gradPhi += temp;
            }
            ret[0] = (gradPhi * gradPhi) / 2;
        }
    } elementEnergy;

    class anonymous7 : public Integration::ElementIntegrandBase<LinearAlgebra::Matrix> {
    public:

        void elementIntegrand(const Base::Element* element, const Geometry::PointReference& p, LinearAlgebra::Matrix& ret) {
            const LinearAlgebra::Matrix& coefficients = element->getTimeLevelData(0);
            static Geometry::Jacobian jac(DIM, DIM);
            jac.resize(DIM, DIM); //can only solve one dimension at a time
            element->calcJacobian(p, jac);
            double det = jac.determinant();
            jac /= det*det; //prevent excessive amounts of dividing later on
            const int n(element->getNrOfBasisFunctions());
            ret.resize(n, n);
            ret *= 0;
            static LinearAlgebra::NumericalVector phiDerivI(DIM), phiDeriv(DIM);
            phiDerivI.resize(DIM);
            phiDeriv.resize(DIM);
            static std::vector<double> dx(n), dy(n);
            double dxl(0), dyl(0);
            dx.resize(n);
            dy.resize(n);
            dx.clear();
            dy.clear();
            phiDeriv *= 0;
            for (int i = 0; i < n; ++i) {
                dx[i] = element->basisFunctionDeriv(i, 0, p);
                dy[i] = element->basisFunctionDeriv(i, 1, p);
                
                //use phiDerivI to collect contributions toward phiDeriv
                element->basisFunctionDeriv(i, p, phiDerivI);
                phiDerivI *= coefficients(1, i);
                phiDeriv += phiDerivI;
                dxl += dx[i] * coefficients(1, i);
                dyl += dy[i] * coefficients(1, i);
            }
            
            for (int i = 0; i < n; ++i) {
                
                //and set phiDerivI to its correct value here
                element->basisFunctionDeriv(i, p, phiDerivI);
                for (int j = 0; j < n; ++j) {
                    ret(i, j) += (2 * jac(1, 1) * dy[j] * dx[i] * dxl
                            - (jac(1, 0) * dy[j] + jac(1, 1) * dx[j])*(dy[i] * dxl + dx[i] * dyl)
                            + 2 * jac(1, 0) * dx[j] * dy[i] * dyl) // det / det
                            - (phiDerivI * phiDeriv)*(jac(0, 0) * dy[j] - jac(0, 1) * dx[j]) * std::fabs(det);
                }
            }
        }
    } nonlinearMatrix; //actually C_{i,l}(\eta_j)\phi_l (still requires summing different contributions in the same mesh column)

    /*class anonymous8:public Integration::ElementIntegrandBase<LinearAlgebra::Matrix>{
    public:
        void elementIntegrand(const Base::Element* element, const Geometry::PointReference& p, LinearAlgebra::Matrix& ret){
            const LinearAlgebra::Matrix& coefficients=element->getTimeLevelData(0);
            Geometry::Jacobian jac(DIM,DIM);
            element->calcJacobian(p,jac);
            double det=jac.determinant();
            const int n(element->getNrOfBasisFunctions());
            ret.resize(n,n);
            ret*=0;
            int nodes=element->getReferenceGeometry()->getCodim1ReferenceGeometry(0)->getNumberOfNodes();
            LinearAlgebra::NumericalVector phiDerivI(DIM),phiDerivL(DIM);
            std::vector<double> dx(n),dy(n);
            for(int i=0;i<n;++i){
                dx[i]=element->basisFunctionDeriv(i,0,p);
                dy[i]=element->basisFunctionDeriv(i,1,p);
            }
            for(int i=0;i<n;++i){
                element->basisFunctionDeriv(i,p,phiDerivI);
                for(int j=0;j<n;++j){
                    for(int l=0;l<n;++l){
                        if((element->getFace(3)->getFaceType()==Geometry::OPEN_BC)&&(i/nodes!=l/nodes)){
                            //the transpose (considered from phi) of the other nonlinearMatrix
                            element->basisFunctionDeriv(l,p,phiDerivL);
                            ret(l,j)+=(2*jac(1,1)*dy[j]*dx[i]*dx[l]
                                    -(jac(1,0)*dy[j]+jac(1,1)*dx[j])*(dy[i]*dx[l]+dx[i]*dy[l])
                                    +2*jac(1,0)*dx[j]*dy[i]*dy[l])*coefficients(1,i)/det/det;
                            ret(l,j)-=(phiDerivI*phiDerivL)*(jac(0,0)*dy[j]-jac(0,1)*dx[j])*coefficients(1,i)/std::fabs(det);
                        }
                    }
                }
            }
        }
    }nonlinearMatrix2;*/

    class anonymous9: public Integration::FaceIntegrandBase<LinearAlgebra::Matrix>{
    public:
        void faceIntegrand(const Base::Face* face, const LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p, LinearAlgebra::Matrix& ret){
            Geometry::PointPhysical physicalNode(DIM);
            int n=face->getNrOfBasisFunctions();
            ret.resize(n,n);
            ret*=0;
            std::vector<unsigned int> faceNodes;
            double halfL(10.);
            face->getPtrElementLeft()->getPhysicalGeometry()->getLocalFaceNodeIndices(face->localFaceNumberLeft(),faceNodes);
            for(int i=0;i<n;++i){
                for(int j:faceNodes){
                    for(int k:faceNodes){
                        face->getPtrElementLeft()->getPhysicalGeometry()->getNodeCoordinates(k, physicalNode);
                        if(physicalNode[0]<=halfL){
                            double x = (halfL - physicalNode[0]) / (halfL - waveMakerPosition());
                            ret(i,j) -= waveMakerDerivative() * face->basisFunction(i,p) * face->basisFunction(j,p) * x * face->basisFunctionDeriv(k,0,p) / Base::L2Norm(normal);
                        }
                    }
                }
            }
        }
    }massDerivative;

    class anonymous10 : public Integration::FaceIntegrandBase<LinearAlgebra::Matrix> {
    public:

        void faceIntegrand(const Base::Face* face, const LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p, LinearAlgebra::Matrix& ret) {
            Geometry::PointPhysical physicalNode(DIM);
            LinearAlgebra::NumericalVector phiDerivI(DIM);
            int n = face->getNrOfBasisFunctions();
            ret.resize(n, n);
            ret *= 0;
            std::vector<unsigned int> faceNodes;
            face->getPtrElementLeft()->getPhysicalGeometry()->getLocalFaceNodeIndices(face->localFaceNumberLeft(), faceNodes);
            double halfL(10.);
            for (int i = 0; i < n; ++i) {
                face->basisFunctionDeriv(i, p, phiDerivI);
                for (int k : faceNodes) {
                    face->getPtrElementLeft()->getPhysicalGeometry()->getNodeCoordinates(k, physicalNode);
                    for (int j = 0; j < n; ++j) {
                        if (physicalNode[0] <= halfL) {//artificial flow contribution towards the wave make (to compensate for node movement)
                            double x = (halfL - physicalNode[0]) / (halfL - waveMakerPosition());
                            ret(i, j) -= face->basisFunction(j, p) * phiDerivI[0] * waveMakerDerivative() * x * face->basisFunction(k, p);
                        }
                    }
                }
            }
        }
    } massDerivative2;

    class anonymous11 : public Integration::FaceIntegrandBase<LinearAlgebra::NumericalVector> {
    public:

        void faceIntegrand(const Base::Face* face, const LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret) {
            Geometry::PointPhysical pPhys(DIM);
            face->referenceToPhysical(p, pPhys);
            if (face->getFaceType() == Geometry::WALL_BC && std::fabs(pPhys[0] - waveMakerPosition()) < 1e-9) {
                int n = face->getNrOfBasisFunctions();
                ret.resize(n);
                for (int i = 0; i < n; ++i) {//mass conservation (for exact derivatives)
                    ret(i) = waveMakerDerivative() * face->basisFunction(i, p);
                }
            } else {
                //if(std::fabs(pPhys[0])-waveMakerPosition()<1e-5){
                //    std::cout<<pPhys<<" "<<waveMakerPosition()<<std::endl;
                //}
                int n = face->getNrOfBasisFunctions();
                ret.resize(n);
                ret *= 0;
            }
        }
    } waveMakerTerm;

    class anonymous12 : public Integration::FaceIntegrandBase<LinearAlgebra::Matrix> {
    public:

        void faceIntegrand(const Base::Face* face, const LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p, LinearAlgebra::Matrix& ret) {
            Geometry::PointPhysical pPhys(DIM);
            face->referenceToPhysical(p, pPhys);
            if (face->getFaceType() == Geometry::WALL_BC && std::fabs(pPhys[0] - waveMakerPosition()) < 1e-9) {
                int n = face->getNrOfBasisFunctions();
                std::vector<unsigned int> nodes;
                face->getPtrElementLeft()->getReferenceGeometry()->getCodim1EntityLocalIndices(face->localFaceNumberLeft(), nodes);
                ret.resize(n, n);
                ret *= 0;
                LinearAlgebra::NumericalVector phiDerivJ(DIM);
                for (int i = 0; i < n; ++i) {
                    for (int j : nodes) {
                        face->getPtrElementLeft()->getPhysicalGeometry()->getLocalNodeCoordinates(j,pPhys);
                        face->basisFunctionDeriv(j, p, phiDerivJ);
                        ret(j, i) = - waveMakerDerivative() * face->basisFunction(i, p) * face->basisFunctionDeriv(j, DIM - 1, p) / Base::L2Norm(normal);
                    }
                }
            } else {
                //if(std::fabs(pPhys[0])-waveMakerPosition()<1e-5){
                //    std::cout<<pPhys<<" "<<waveMakerPosition()<<std::endl;
                //}
                int n = face->getNrOfBasisFunctions();
                ret.resize(n, n);
                ret *= 0;
            }
        }
    } waveMakerDerivativeIntegrand;

    void writeToTecplotFile(const ElementT* element, const PointReferenceT& p, std::ostream& out) {
        LinearAlgebra::NumericalVector value(2);
        element->getSolution(0, p, value);
        out << value[0] << " " << value[1];

        LinearAlgebra::NumericalVector derivative(DIM), functionDeriv(DIM);
        const LinearAlgebra::Matrix& data = element->getTimeLevelData(0);
        for (int i = 0; i < element->getNrOfBasisFunctions(); ++i) {
            element->basisFunctionDeriv(i, p, functionDeriv);
            functionDeriv *= data(1, i);
            derivative += functionDeriv;
        }
        out << " " << derivative[0] << " " << derivative[1];

        Geometry::PointPhysical pPhys(DIM);
        element->referenceToPhysical(p, pPhys);
        exactSolution(t, pPhys, value);
        out << " " << value[0] << " " << value[1] + t;
    }

    //computes the mass matrix at the free surface

    void constructM() {
        Integration::FaceIntegral integral(false);
        LinearAlgebra::Matrix result;
        for (Base::Face* face : meshes_[1]->getFacesList()) {
            result.resize(face->getNrOfBasisFunctions(), face->getNrOfBasisFunctions());
            if (face->getFaceType() == Geometry::OPEN_BC) {
                integral.integrate(face, &massIntegrand, result);
                face->setFaceMatrix(result, 0);
                integral.integrate(face, &massDerivative2, result);
                face->setFaceMatrix(result, 1);
                integral.integrate(face, &massDerivative, result);
                face->setFaceMatrix(result, 2);
            } else {
                result *= 0;
                face->setFaceMatrix(result, 0);
                face->setFaceMatrix(result, 1);
                face->setFaceMatrix(result, 2);
            }
        }
    }

    void readInitialConditions() {
        Geometry::PointReference p(DIM - 1);
        Geometry::PointPhysical pPhys(DIM);
        LinearAlgebra::Matrix initialconditions(2, std::pow(2, DIM));
        LinearAlgebra::NumericalVector functionvalue(2);
        for (Base::Face* face : meshes_[0]->getFacesList()) {
            if (face->getFaceType() == Geometry::OPEN_BC) {
                int n = face->getReferenceGeometry()->getNumberOfNodes();
                for (int i = 0; i < n; ++i) {
                    face->getReferenceGeometry()->getNode(i, p);
                    face->referenceToPhysical(p, pPhys);
                    initialConditions(pPhys, functionvalue);
                    initialconditions(0, i + n) = functionvalue[0]; // double(n);
                    initialconditions(1, i + n) = functionvalue[1]; // double(n);
                    //if (std::fabs(pPhys[0]-waveMakerPosition()) < 1e-9 || std::fabs(pPhys[0] - L) < 1e-9) {
                    //    initialconditions(0, i + n) *= 2;
                    //    initialconditions(1, i + n) *= 2;
                    //}
                }
                face->getPtrElementLeft()->setTimeLevelData(0, initialconditions);
            }
        }
    }

    void constructS() {
        Integration::ElementIntegral elIntegral(false);
        Integration::FaceIntegral faIntegral(false);
        LinearAlgebra::Matrix result;
        LinearAlgebra::NumericalVector resultVec;
        for (Base::Element* element : meshes_[0]->getElementsList()) {
            result.resize(element->getNrOfBasisFunctions(), element->getNrOfBasisFunctions());
            elIntegral.integrate(element, &stifnessIntegrand, result);
            element->setElementMatrix(result, 0);
        }
        for (Base::Face* face : meshes_[0]->getFacesList()) {
            resultVec.resize(face->getNrOfBasisFunctions());
            if (face->getFaceType() == Geometry::WALL_BC) {
                faIntegral.integrate(face, &waveMakerTerm, resultVec);
                face->setFaceVector(resultVec, 0);
            } else {
                resultVec *= 0;
                face->setFaceVector(resultVec, 0);
            }
        }
    }

    //constructs a square edition of C. Another matrix will link the column entries to eta

    void constructC() {
        Integration::ElementIntegral elIntegral(false);
        Integration::FaceIntegral faIntegral(false);
        LinearAlgebra::Matrix result, missingMultiplier;
        for (Base::Element* element : meshes_[0]->getElementsList()) {
            result.resize(element->getNrOfBasisFunctions(), element->getNrOfBasisFunctions());
            //it is difficult to multiply by gamma inside the integral
            elIntegral.integrate(element, &nonlinearMatrix, result);

            //so do it here
            missingMultiplier.resize(element->getNrOfBasisFunctions(), element->getNrOfBasisFunctions());
            for (int i = 0; i < element->getNrOfBasisFunctions(); ++i) {
                missingMultiplier(i, i) = 1 - meshes_[2]->getNodes()[element->getPhysicalGeometry()->getNodeIndex(i)][DIM - 1] / 0.6;
            }
            //std::cout<<result<<std::endl<<missingMultiplier<<std::endl;
            result = result*missingMultiplier; //hpGEM matrixes appear to be transposed when compared to PETSc matrices
            //std::cout<<result<<std::endl;
            element->setElementMatrix(result, 1);
        }
        for (Base::Face* face : meshes_[0]->getFacesList()) {
            result.resize(face->getNrOfBasisFunctions(), face->getNrOfBasisFunctions());
            if (face->getFaceType() == Geometry::WALL_BC) {
                faIntegral.integrate(face, &waveMakerDerivativeIntegrand, result);
                const Base::Element * left(face->getPtrElementLeft());
                for (int i = 0; i < left->getNrOfBasisFunctions(); ++i) {
                    missingMultiplier(i, i) = 1 - meshes_[2]->getNodes()[left->getPhysicalGeometry()->getNodeIndex(i)][DIM - 1] / 0.6;
                }
                result = result * missingMultiplier;
            } else {
                result *= 0;
            }
            face->setFaceMatrix(result, 0);
        }
    }

    void moveNodes() {
        const Geometry::PointPhysical* firstNode, *referenceFirstNode;
        firstNode = &meshes_[0]->getNodes()[0];
        referenceFirstNode = &meshes_[2]->getNodes()[0];
        Geometry::PointReference p(DIM - 1), pElement(DIM);
        Geometry::PointPhysical pPhys(DIM);
        std::vector<unsigned int> faceNodes;
        LinearAlgebra::NumericalVector solution(2);
        static bool first(true);
        for (Base::Face* face : meshes_[0]->getFacesList()) {
            if (face->getFaceType() == Geometry::OPEN_BC) {
                face->getPtrElementLeft()->getPhysicalGeometry()->getGlobalFaceNodeIndices(face->localFaceNumberLeft(), faceNodes);
                for (int i = 0; i < faceNodes.size(); ++i) {
                    face->getReferenceGeometry()->getNode(i, p);
                    if (t == 0) {
                        if (first) {
                            solution[0] = 1;
                        } else {
                            face->referenceToPhysical(p, pPhys);
                            initialConditions(pPhys, solution);
                        }
                    } else {
                        face->mapRefFaceToRefElemL(p, pElement);
                        face->getPtrElementLeft()->getSolution(0, pElement, solution);
                    }
                    for (Geometry::PointPhysical* point : nodesBelowASurfaceNode_[faceNodes[i]]) {
                        double zOld = (*(point - firstNode + referenceFirstNode))[DIM - 1]; //there is probably a cleaner way to do this
                        //(*point)[DIM - 1] = -zOld + (1. - zOld)*(solution[0]);
                        //std::cout<<*point<<std::endl;
                        (*point)[DIM - 1] = (-zOld + (0.6 - zOld)*(solution[0]) / 0.6);
                    }
                }
            }
        }
        first = false;
        static double rOld(0.);
        double rNew(waveMakerPosition()), halfL(10.);
        int differenceBetweenMoves(0);
        for (const Geometry::PointPhysical& node : meshes_[0]->getNodes()) {
            if ((node)[0] < halfL) {
                differenceBetweenMoves++;
                double xOld = (*(&node - firstNode + referenceFirstNode))[0];
                const_cast<Geometry::PointPhysical&> (node)[0] = (halfL - xOld) / halfL * rNew + xOld; //(halfL - rNew)*((node)[0] - halfL) / (halfL - rOld) + halfL;
            }
        }
        firstNode = &meshes_[1]->getNodes()[0];
        for (const Geometry::PointPhysical& node : meshes_[1]->getNodes()) {
            if ((node)[0] < halfL) {
                differenceBetweenMoves--;
                double xOld = (*(&node - firstNode + referenceFirstNode))[0];
                const_cast<Geometry::PointPhysical&> (node)[0] = (halfL - xOld) / halfL * rNew + xOld; //(halfL - rNew)*((node)[0] - halfL) / (halfL - rOld) + halfL;
            }
        }
        assert(differenceBetweenMoves == 0);
        for (Base::Element* el : meshes_[0]->getElementsList()) {
            const_cast<Geometry::MappingReferenceToPhysical*> (el->getReferenceToPhysicalMap())->reinit(el->getPhysicalGeometry());
        }
        for (Base::Element* el : meshes_[1]->getElementsList()) {
            const_cast<Geometry::MappingReferenceToPhysical*> (el->getReferenceToPhysicalMap())->reinit(el->getPhysicalGeometry());
        }
        rOld = rNew;
    }

    void getSurfaceIS(Utilities::GlobalPetscMatrix& S, IS* surface, IS* rest) {
        int n(0);
        std::vector<int> facePositions;
        for (const Base::Face* face : meshes_[0]->getFacesList()) {
            if (face->getFaceType() == Geometry::OPEN_BC) {
                S.getMatrixBCEntries(face, n, facePositions);
            }
        }
        ISCreateGeneral(MPI_COMM_WORLD, n, &facePositions[0], PETSC_COPY_VALUES, surface);
        MatGetSize(S, &n, NULL);
        ISSort(*surface);
        ISComplement(*surface, 0, n, rest);
        ISDestroy(surface);
        ISComplement(*rest, 0, n, surface); //the complement of the complement does not contain duplicates
    }

    void printError() {
        Integration::FaceIntegral integral(false);
        LinearAlgebra::NumericalVector totalError(2), contribution(2);
        for (Base::Face* face : meshes_[0]->getFacesList()) {
            integral.integrate(face, &error, contribution);
            totalError += contribution;
            contribution[0] = 0;
            contribution[1] = 0;
        }
        std::cout << "t: " << t << " eta: " << sqrt(totalError[0]) << " phi: " << sqrt(totalError[1]) << std::endl;
    }

    void computeEnergy() {
        Integration::ElementIntegral elIntegral(false);
        Integration::FaceIntegral faIntegral(false);
        LinearAlgebra::NumericalVector totalEnergy(1), contribution(1);
        for (Base::Face* face : meshes_[0]->getFacesList()) {
            faIntegral.integrate(face, &faceEnergy, contribution);
            totalEnergy += contribution;
            contribution[0] = 0;
        }
        for (Base::Element* element : meshes_[0]->getElementsList()) {
            elIntegral.integrate(element, &elementEnergy, contribution);
            totalEnergy += contribution;
            contribution[0] = 0;
        }
        std::cout << "Energy: " << totalEnergy[0] << std::endl;
    }

    bool solveStormerVerlet() {
        bool first(true);
        std::cout.precision(14);
        moveNodes();
        constructM();
        readInitialConditions();

        Utilities::GlobalPetscVector eta(meshes_[0]), phi(meshes_[0]);
        Utilities::GlobalPetscMatrix M(meshes_[1], -1, 0), D(meshes_[1], -1, 1), dMdt(meshes_[1],-1,2);
        Vec phiS, phiOther, etaActually, interiorRHS, surfaceRHS, surfaceOld, surfaceUpdate, interiorUpdate, testThings;
        Mat surfaceMass, interiorStifness, surfaceStifness, mixStifness, backStifness, surfaceDMass, surfaceDMdt;
        Mat surfaceC, interiorC, surfaceB, interiorB;

        VecDuplicate(phi, &testThings);

        IS isSurface, isRest;
        getSurfaceIS(M, &isSurface, &isRest);

        std::ofstream outFile("output.dat");
        Output::TecplotDiscontinuousSolutionWriter writeFunc(outFile, "test", "01", "eta_num,phi_num,phi_x,phi_y,eta_exact,phi_exact");

        //deal with initial conditions
        double g(9.8), dt(0.00125);
        eta.constructFromTimeLevelData(0, 0);
        phi.constructFromTimeLevelData(0, 1);

        int numberOfSnapshots(401);
        int numberOfTimeSteps(80);
        //int numberOfTimeSteps(16); //placeholder parameters
        //int numberOfSnapshots(1 + n_ / 32. * 320);
        //double dt(0.01551125);

        MatGetSubMatrix(M, isSurface, isSurface, MAT_INITIAL_MATRIX, &surfaceMass);
        MatGetSubMatrix(D, isSurface, isSurface, MAT_INITIAL_MATRIX, &surfaceDMass);
        MatGetSubMatrix(dMdt, isSurface, isSurface, MAT_INITIAL_MATRIX, &surfaceDMdt);

        KSP interior, surface;
        KSPCreate(MPI_COMM_WORLD, &surface);
        KSPSetOperators(surface, surfaceMass, surfaceMass);
        KSPSetFromOptions(surface);
        KSPConvergedReason conferge;
        int iterations;

        VecGetSubVector(phi, isSurface, &phiS);
        VecGetSubVector(phi, isRest, &phiOther);
        VecGetSubVector(eta, isSurface, &etaActually);
        VecDuplicate(phiS, &surfaceRHS);
        VecDuplicate(phiS, &surfaceOld);
        VecDuplicate(phiS, &surfaceUpdate);
        VecDuplicate(phiOther, &interiorRHS);
        VecDuplicate(phiOther, &interiorUpdate);
        VecCopy(etaActually, surfaceRHS);
        //KSPSolve(surface,surfaceRHS,etaActually);
        //KSPGetConvergedReason(surface,&conferge);
        //KSPGetIterationNumber(surface,&iterations);
        //std::cout<<"Finalizing interpolation (1): KSP solver ended because of "<<KSPConvergedReasons[conferge]<<" in "<<iterations<<" iterations."<<std::endl;
        VecRestoreSubVector(eta, isSurface, &etaActually);
        VecRestoreSubVector(phi, isSurface, &phiS);
        VecRestoreSubVector(phi, isRest, &phiOther);


        eta.writeTimeLevelData(0, 0);
        moveNodes();
        readInitialConditions();
        eta.constructFromTimeLevelData(0, 0);
        CHKERRQ(VecGetSubVector(eta, isSurface, &etaActually));
        VecCopy(etaActually, surfaceRHS);
        KSPSetOperators(surface, surfaceMass, surfaceMass);
        //KSPSolve(surface,surfaceRHS,etaActually);
        //KSPGetConvergedReason(surface,&conferge);
        //KSPGetIterationNumber(surface,&iterations);
        //std::cout<<"Finalizing interpolation (1): KSP solver ended because of "<<KSPConvergedReasons[conferge]<<" in "<<iterations<<" iterations."<<std::endl;
        VecRestoreSubVector(eta, isSurface, &etaActually);
        eta.writeTimeLevelData(0, 0);

        //accurate construction of phi
        phi.constructFromTimeLevelData(0, 1);
        CHKERRQ(VecGetSubVector(phi, isSurface, &phiS));
        VecGetSubVector(eta, isSurface, &etaActually);
        VecCopy(phiS, surfaceRHS);
        //KSPSolve(surface,surfaceRHS,phiS);
        //KSPGetConvergedReason(surface,&conferge);
        //KSPGetIterationNumber(surface,&iterations);
        //std::cout<<"Finalizing interpolation (2): KSP solver ended because of "<<KSPConvergedReasons[conferge]<<" in "<<iterations<<" iterations."<<std::endl;
        VecRestoreSubVector(phi, isSurface, &phiS);
        phi.writeTimeLevelData(0, 1);

        constructM();
        M.reset();
        D.reset();
        dMdt.reset();
        MatGetSubMatrix(M, isSurface, isSurface, MAT_REUSE_MATRIX, &surfaceMass);
        MatGetSubMatrix(D, isSurface, isSurface, MAT_REUSE_MATRIX, &surfaceDMass);
        MatGetSubMatrix(dMdt, isSurface, isSurface, MAT_REUSE_MATRIX, &surfaceDMdt);

        constructS();
        constructC();

        Vec scaleVec;
        VecDuplicate(phi, &scaleVec);
        VecSet(scaleVec, 1.);
        VecGetSubVector(scaleVec, isSurface, &surfaceUpdate);
        VecSet(surfaceUpdate, dt / 2.);
        VecRestoreSubVector(scaleVec, isSurface, &surfaceUpdate);

        Utilities::GlobalPetscMatrix S(meshes_[0], 0, -1), SquareC(meshes_[0], 1, -1), squareB(meshes_[0], -1, 0);
        Utilities::GlobalPetscVector W(meshes_[0], -1, 0);
        W.assemble();
        Mat C, CollectCulumns, B;
        Vec interiorW, surfaceW;
        MatCreate(MPI_COMM_WORLD, &CollectCulumns);
        int n;
        MatGetSize(SquareC, &n, NULL); //extra rows hardly matter because they are empty anyway
        MatSetSizes(CollectCulumns, PETSC_DECIDE, PETSC_DECIDE, n, n);
        MatSetUp(CollectCulumns);
        for (std::pair<int, int> pair : surfaceNodesAboveANode_) {
            //if((pair.second+1)/(n_+1)*(n_+1)==(pair.second+1)){
            //MatSetValue(CollectCulumns,0,pair.first-n_-pair.first/(n_+1),1.,INSERT_VALUES);
            //}else{
            int first(pair.first/*-pair.first/(n_+1)*/), second(pair.second/*-pair.second/(n_+1)*/);
            if (DIM == 3) {
                if ((second + 1) / 5 * 5 == second + 1) {

                } else {
                    MatSetValue(CollectCulumns, second - second / 5, first - first / 5, 1., INSERT_VALUES);
                }
            } else {
                MatSetValue(CollectCulumns, second, first, 1., INSERT_VALUES);
            }
            //}
        }


        MatAssemblyBegin(CollectCulumns, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(CollectCulumns, MAT_FINAL_ASSEMBLY);


        MatMatMult(CollectCulumns, SquareC, MAT_INITIAL_MATRIX, 1., &C);
        MatMatMult(CollectCulumns, squareB, MAT_INITIAL_MATRIX, 1., &B);


        MatGetSubMatrix(C, isSurface, isSurface, MAT_INITIAL_MATRIX, &surfaceC);
        MatGetSubMatrix(C, isSurface, isRest, MAT_INITIAL_MATRIX, &interiorC);
        MatGetSubMatrix(B, isSurface, isSurface, MAT_INITIAL_MATRIX, &surfaceB);
        MatGetSubMatrix(B, isSurface, isRest, MAT_INITIAL_MATRIX, &interiorB);
        MatGetSubMatrix(S, isRest, isRest, MAT_INITIAL_MATRIX, &interiorStifness);
        MatGetSubMatrix(S, isSurface, isSurface, MAT_INITIAL_MATRIX, &surfaceStifness);
        MatGetSubMatrix(S, isRest, isSurface, MAT_INITIAL_MATRIX, &mixStifness);
        MatGetSubMatrix(S, isSurface, isRest, MAT_INITIAL_MATRIX, &backStifness);
        VecGetSubVector(W, isSurface, &surfaceW);
        VecGetSubVector(W, isRest, &interiorW);

        Mat BigMatrix;
        MatDuplicate(S, MAT_DO_NOT_COPY_VALUES, &BigMatrix);
        Vec BigVector, BigRHS;
        VecDuplicate(phi, &BigVector);
        VecDuplicate(phi, &BigRHS);
        VecGetSubVector(BigRHS, isSurface, &surfaceRHS);
        VecGetSubVector(BigRHS, isRest, &interiorRHS);
        KSPCreate(MPI_COMM_WORLD, &interior);
        KSPSetOperators(interior, BigMatrix, BigMatrix);
        MatSetOption(BigMatrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
        MatSetOption(BigMatrix, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
        KSPSetFromOptions(interior);

        double accuracy;

        KSPGetTolerances(interior, &accuracy, NULL, NULL, NULL);
        KSPSetTolerances(interior, 1e-12, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);

        writeFunc.write(meshes_[0], "solution at time 0", false, this);
        printError();
        //computeEnergy();//does not make sense since phi (internal) is not calculated yet
        for (int i = 1; i < numberOfSnapshots; ++i) {
            VecGetSubVector(phi, isSurface, &phiS);
            VecGetSubVector(phi, isRest, &phiOther);
            VecGetSubVector(eta, isSurface, &etaActually);
            for (int j = 0; j < numberOfTimeSteps; ++j) {

                MatMult(surfaceMass, phiS, surfaceOld);
                VecScale(surfaceOld, -2. / g / dt);

                t += dt / 2.;
                VecRestoreSubVector(phi, isSurface, &phiS);
                VecRestoreSubVector(phi, isRest, &phiOther);
                VecRestoreSubVector(W, isSurface, &surfaceW);
                VecRestoreSubVector(W, isRest, &interiorW);
                phi.writeTimeLevelData(0, 1);
                moveNodes();
                constructM();
                M.reset();//discrepancy with the Elena code
                D.reset();
                dMdt.reset();
                //MatGetSubMatrix(M, isSurface, isSurface, MAT_REUSE_MATRIX, &surfaceMass);
                MatGetSubMatrix(D, isSurface, isSurface, MAT_REUSE_MATRIX, &surfaceDMass);
                MatGetSubMatrix(dMdt, isSurface, isSurface, MAT_REUSE_MATRIX, &surfaceDMdt);
                MatAXPY(M,dt/2.,dMdt,SAME_NONZERO_PATTERN);




                constructC();
                constructS();
                W.assemble();
                SquareC.reset();
                squareB.reset();
                S.reset();
                MatGetSubMatrix(S, isRest, isRest, MAT_REUSE_MATRIX, &interiorStifness);
                MatGetSubMatrix(S, isSurface, isSurface, MAT_REUSE_MATRIX, &surfaceStifness);
                MatGetSubMatrix(S, isRest, isSurface, MAT_REUSE_MATRIX, &mixStifness);
                MatGetSubMatrix(S, isSurface, isRest, MAT_REUSE_MATRIX, &backStifness);
                MatMatMult(CollectCulumns, SquareC, MAT_REUSE_MATRIX, 1., &C);
                MatMatMult(CollectCulumns, squareB, MAT_REUSE_MATRIX, 1., &B);

                MatGetSubMatrix(C, isSurface, isSurface, MAT_REUSE_MATRIX, &surfaceC);
                MatGetSubMatrix(C, isSurface, isRest, MAT_REUSE_MATRIX, &interiorC);
                MatGetSubMatrix(B, isSurface, isSurface, MAT_REUSE_MATRIX, &surfaceB);
                MatGetSubMatrix(B, isSurface, isRest, MAT_REUSE_MATRIX, &interiorB);
                VecGetSubVector(W, isSurface, &surfaceW);
                VecGetSubVector(W, isRest, &interiorW);
                VecGetSubVector(phi, isSurface, &phiS);
                VecGetSubVector(phi, isRest, &phiOther);

                MatMultAdd(surfaceMass, etaActually, surfaceOld, surfaceOld);
                VecScale(surfaceOld, -g * dt / 2.);
                //MatMultAdd(surfaceMass, phiS, surfaceOld, surfaceOld);//move to beginning of time loop (use old M)
                MatMult(surfaceDMass, phiS, surfaceRHS);
                MatMultAdd(surfaceB, phiS, surfaceRHS, surfaceRHS);
                MatMultAdd(interiorB, phiOther, surfaceRHS, surfaceRHS);
                VecScale(surfaceRHS, 2);
                MatMultAdd(surfaceC, phiS, surfaceRHS, surfaceRHS);
                MatMultAdd(interiorC, phiOther, surfaceRHS, surfaceRHS);
                VecScale(surfaceRHS, dt / 4.);
                MatMultAdd(surfaceMass, phiS, surfaceRHS, surfaceRHS);
                VecAYPX(surfaceRHS, -1., surfaceOld);
                MatMult(interiorStifness, phiOther, interiorRHS);
                MatMultAdd(mixStifness, phiS, interiorRHS, interiorRHS);
                VecAXPY(interiorRHS, 1., interiorW);
                double normSurface, normInterior;
                VecNorm(etaActually, NORM_2, &normSurface);
                VecNorm(phi, NORM_2, &normInterior);
                int numberOfNewtonSteps(0);
                double relativeTolerance = accuracy * (normInterior + normSurface);
                VecNorm(surfaceRHS, NORM_2, &normSurface);
                VecNorm(interiorRHS, NORM_2, &normInterior);
                relativeTolerance=std::max(relativeTolerance,1e-16);
                while (normInterior + normSurface > relativeTolerance) {
                    numberOfNewtonSteps++;
                    std::cout<<"In first Newton iteration with solution errors: "<<normInterior<<" and "<<normSurface<<std::endl;

                    MatCopy(S, BigMatrix, DIFFERENT_NONZERO_PATTERN);
                    MatZeroRowsIS(BigMatrix, isSurface, 0., NULL, NULL);

                    CHKERRABORT(MPI_COMM_WORLD, MatAXPY(BigMatrix, -1., M, DIFFERENT_NONZERO_PATTERN));
                    CHKERRABORT(MPI_COMM_WORLD, MatAXPY(BigMatrix, -dt / 2., D, DIFFERENT_NONZERO_PATTERN));
                    CHKERRABORT(MPI_COMM_WORLD, MatAXPY(BigMatrix, -dt / 2., C, DIFFERENT_NONZERO_PATTERN));
                    CHKERRABORT(MPI_COMM_WORLD, MatAXPY(BigMatrix, -dt / 2., B, DIFFERENT_NONZERO_PATTERN));

                    //MatView(C,PETSC_VIEWER_STDOUT_WORLD);

                    VecRestoreSubVector(BigRHS, isSurface, &surfaceRHS);
                    VecRestoreSubVector(BigRHS, isRest, &interiorRHS);
                    if (first)
                        KSPSetOperators(interior, BigMatrix, BigMatrix);
                    else
                        KSPSetOperators(interior, BigMatrix, BigMatrix);
                    KSPSolve(interior, BigRHS, BigVector);
                    KSPGetConvergedReason(interior, &conferge);
                    KSPGetIterationNumber(interior, &iterations);
                    //std::cout<<"KSP solver ended because of "<<KSPConvergedReasons[conferge]<<" in "<<iterations<<" iterations."<<std::endl;
                    VecGetSubVector(BigVector, isSurface, &surfaceUpdate);
                    VecGetSubVector(BigVector, isRest, &interiorUpdate);
                    VecGetSubVector(BigRHS, isSurface, &surfaceRHS);
                    VecGetSubVector(BigRHS, isRest, &interiorRHS);
                    VecAXPY(phiS, -1.0, surfaceUpdate);
                    VecAXPY(phiOther, -1.0, interiorUpdate);
                    VecRestoreSubVector(BigVector, isSurface, &surfaceUpdate);
                    VecRestoreSubVector(BigVector, isRest, &interiorUpdate);

                    VecRestoreSubVector(phi, isSurface, &phiS);
                    VecRestoreSubVector(phi, isRest, &phiOther);
                    phi.writeTimeLevelData(0, 1);
                    constructC();
                    SquareC.reset();
                    //squareB.reset();//not based on \phi
                    MatMatMult(CollectCulumns, SquareC, MAT_REUSE_MATRIX, 1., &C);
                    //MatMatMult(CollectCulumns, squareB, MAT_REUSE_MATRIX, 1., &B);
                    MatGetSubMatrix(C, isSurface, isSurface, MAT_REUSE_MATRIX, &surfaceC);
                    MatGetSubMatrix(C, isSurface, isRest, MAT_REUSE_MATRIX, &interiorC);
                    //MatGetSubMatrix(B, isSurface, isSurface, MAT_REUSE_MATRIX, &surfaceB);
                    //MatGetSubMatrix(B, isSurface, isRest, MAT_REUSE_MATRIX, &interiorB);
                    VecGetSubVector(phi, isSurface, &phiS);
                    VecGetSubVector(phi, isRest, &phiOther);

                    MatMult(surfaceDMass, phiS, surfaceRHS);
                    MatMultAdd(surfaceB, phiS, surfaceRHS, surfaceRHS);
                    MatMultAdd(interiorB, phiOther, surfaceRHS, surfaceRHS);
                    VecScale(surfaceRHS, 2.);
                    MatMultAdd(surfaceC, phiS, surfaceRHS, surfaceRHS);
                    MatMultAdd(interiorC, phiOther, surfaceRHS, surfaceRHS);
                    VecScale(surfaceRHS, dt / 4.);
                    MatMultAdd(surfaceMass, phiS, surfaceRHS, surfaceRHS);
                    VecAYPX(surfaceRHS, -1., surfaceOld);
                    MatMult(interiorStifness, phiOther, interiorRHS);
                    MatMultAdd(mixStifness, phiS, interiorRHS, interiorRHS);
                    VecAXPY(interiorRHS, 1., interiorW);
                    VecNorm(surfaceRHS, NORM_2, &normSurface);
                    VecNorm(interiorRHS, NORM_2, &normInterior);

                }
                std::cout << "First set of equations solved after " << numberOfNewtonSteps << " newton iterations" << std::endl;

                std::cout.precision(15);
                    VecRestoreSubVector(phi, isSurface, &phiS);
                    VecRestoreSubVector(phi, isRest, &phiOther);
                    
                    /*double *data;
                    int size;
                    VecGetArray(phi,&data);
                    VecGetSize(phi,&size);
                    std::cout<<"phi"<<std::endl;
                    for(int i=0;i<size;++i){
                        std::cout<<data[i]<<std::endl;
                    }
                    VecRestoreArray(phi,&data);
                    
                    std::cout<<"eta"<<std::endl;
                    VecGetArray(etaActually,&data);
                    VecGetSize(etaActually,&size);
                    for(int i=0;i<size;++i){
                        std::cout<<data[i]<<std::endl;
                    }
                    VecRestoreArray(etaActually,&data);
                    std::cout<<std::endl;*/
                    
                    VecGetSubVector(phi, isSurface, &phiS);
                    VecGetSubVector(phi, isRest, &phiOther);
                
                MatMult(surfaceStifness, phiS, surfaceOld);
                MatMultAdd(backStifness, phiOther, surfaceOld, surfaceOld);
                MatMultAdd(surfaceDMass, etaActually, surfaceOld, surfaceOld);
                VecAXPY(surfaceOld, 1., surfaceW);
                VecScale(surfaceOld, dt / 2.);
                MatMultAdd(surfaceMass, etaActually, surfaceOld, surfaceOld);

                MatMult(surfaceStifness, phiS, surfaceRHS);
                MatMultAdd(backStifness, phiOther, surfaceRHS, surfaceRHS);
                MatMultAdd(surfaceDMass, etaActually, surfaceRHS, surfaceRHS);
                VecAXPY(surfaceRHS, 1., surfaceW);
                VecScale(surfaceRHS, -dt / 2.);
                MatMultAdd(surfaceMass, etaActually, surfaceRHS, surfaceRHS);
                VecAYPX(surfaceRHS, -1., surfaceOld);
                MatMult(interiorStifness, phiOther, interiorRHS);
                MatMultAdd(mixStifness, phiS, interiorRHS, interiorRHS);
                VecAXPY(interiorRHS, 1., interiorW);
                VecNorm(etaActually, NORM_2, &normSurface);
                VecNorm(phi, NORM_2, &normInterior);
                numberOfNewtonSteps = 0;
                relativeTolerance = accuracy * (normInterior + normSurface);
                VecNorm(surfaceRHS, NORM_2, &normSurface);
                VecNorm(interiorRHS, NORM_2, &normInterior);
                relativeTolerance=std::max(relativeTolerance,1e-16);
                while (normInterior + normSurface > relativeTolerance) {
                    numberOfNewtonSteps++;
                    std::cout<<"In second Newton iteration with solution errors: "<<normInterior<<" and "<<normSurface<<std::endl;

                    MatCopy(S, BigMatrix, DIFFERENT_NONZERO_PATTERN);
                    MatZeroRowsIS(BigMatrix, isSurface, 0., NULL, NULL);

                    //MatView(B,PETSC_VIEWER_DRAW_WORLD);
                    //MatView(M,PETSC_VIEWER_DRAW_WORLD);
                    //MatView(SquareC,PETSC_VIEWER_STDOUT_WORLD);
                    //MatView(BigMatrix,PETSC_VIEWER_DRAW_WORLD);

                    MatTranspose(BigMatrix, MAT_REUSE_MATRIX, &BigMatrix);
                    MatAXPY(BigMatrix, 1., B, DIFFERENT_NONZERO_PATTERN);
                    MatTranspose(BigMatrix, MAT_REUSE_MATRIX, &BigMatrix);
                    MatAXPY(BigMatrix, 1., C, DIFFERENT_NONZERO_PATTERN);
                    MatDiagonalScale(BigMatrix, NULL, scaleVec);

                    //MatView(BigMatrix,PETSC_VIEWER_DRAW_WORLD);

                    MatAXPY(BigMatrix, -1., M, DIFFERENT_NONZERO_PATTERN);
                    MatTranspose(BigMatrix, MAT_REUSE_MATRIX, &BigMatrix);
                    MatAXPY(BigMatrix, dt / 2., D, DIFFERENT_NONZERO_PATTERN);

                    //MatView(BigMatrix,PETSC_VIEWER_DRAW_WORLD);

                    VecRestoreSubVector(BigRHS, isSurface, &surfaceRHS);
                    VecRestoreSubVector(BigRHS, isRest, &interiorRHS);
                    if (first)
                        KSPSetOperators(interior, BigMatrix, BigMatrix);
                    else
                        KSPSetOperators(interior, BigMatrix, BigMatrix);
                    KSPSolve(interior, BigRHS, BigVector);
                    KSPGetConvergedReason(interior, &conferge);
                    KSPGetIterationNumber(interior, &iterations);
                    //std::cout<<"KSP solver ended because of "<<KSPConvergedReasons[conferge]<<" in "<<iterations<<" iterations."<<std::endl;
                    VecGetSubVector(BigRHS, isSurface, &surfaceRHS);
                    VecGetSubVector(BigRHS, isRest, &interiorRHS);
                    VecGetSubVector(BigVector, isSurface, &surfaceUpdate);
                    VecGetSubVector(BigVector, isRest, &interiorUpdate);
                    VecAXPY(etaActually, -1.0, surfaceUpdate);
                    VecAXPY(phiOther, -1.0, interiorUpdate);
                    VecRestoreSubVector(BigVector, isSurface, &surfaceUpdate);
                    VecRestoreSubVector(BigVector, isRest, &interiorUpdate);

                    VecRestoreSubVector(phi, isSurface, &phiS);
                    VecRestoreSubVector(phi, isRest, &phiOther);
                    VecRestoreSubVector(eta, isSurface, &etaActually);
                    VecRestoreSubVector(W, isSurface, &surfaceW);
                    VecRestoreSubVector(W, isRest, &interiorW);
                    eta.writeTimeLevelData(0, 0);
                    phi.writeTimeLevelData(0, 1);
                    moveNodes();
                    constructC();
                    constructS();
                    W.assemble();
                    S.reset();
                    MatGetSubMatrix(S, isRest, isRest, MAT_REUSE_MATRIX, &interiorStifness);
                    MatGetSubMatrix(S, isSurface, isSurface, MAT_REUSE_MATRIX, &surfaceStifness);
                    MatGetSubMatrix(S, isRest, isSurface, MAT_REUSE_MATRIX, &mixStifness);
                    MatGetSubMatrix(S, isSurface, isRest, MAT_REUSE_MATRIX, &backStifness);
                    SquareC.reset();
                    squareB.reset();
                    MatMatMult(CollectCulumns, SquareC, MAT_REUSE_MATRIX, 1., &C);
                    MatMatMult(CollectCulumns, squareB, MAT_REUSE_MATRIX, 1., &B);
                    MatGetSubMatrix(C, isSurface, isSurface, MAT_REUSE_MATRIX, &surfaceC);
                    MatGetSubMatrix(C, isSurface, isRest, MAT_REUSE_MATRIX, &interiorC);
                    MatGetSubMatrix(B, isSurface, isSurface, MAT_REUSE_MATRIX, &surfaceB);
                    MatGetSubMatrix(B, isSurface, isRest, MAT_REUSE_MATRIX, &interiorB);
                    VecGetSubVector(phi, isSurface, &phiS);
                    VecGetSubVector(phi, isRest, &phiOther);
                    VecGetSubVector(W, isSurface, &surfaceW);
                    VecGetSubVector(W, isRest, &interiorW);
                    VecGetSubVector(eta, isSurface, &etaActually);

                    MatMult(surfaceStifness, phiS, surfaceRHS);
                    MatMultAdd(backStifness, phiOther, surfaceRHS, surfaceRHS);
                    MatMultAdd(surfaceDMass, etaActually, surfaceRHS, surfaceRHS);
                    VecAXPY(surfaceRHS, 1., surfaceW);
                    VecScale(surfaceRHS, -dt / 2.);
                    MatMultAdd(surfaceMass, etaActually, surfaceRHS, surfaceRHS);
                    VecAYPX(surfaceRHS, -1., surfaceOld);
                    MatMult(interiorStifness, phiOther, interiorRHS);
                    MatMultAdd(mixStifness, phiS, interiorRHS, interiorRHS);
                    VecAXPY(interiorRHS, 1., interiorW);
                    VecNorm(surfaceRHS, NORM_2, &normSurface);
                    VecNorm(interiorRHS, NORM_2, &normInterior);
                }
                std::cout << "Second set of equations solved after " << numberOfNewtonSteps << " newton iterations" << std::endl;


            VecRestoreSubVector(phi, isSurface, &phiS);
            VecRestoreSubVector(phi, isRest, &phiOther);
            
            //VecView(phi,PETSC_VIEWER_STDOUT_WORLD);
            //VecView(etaActually,PETSC_VIEWER_STDOUT_WORLD);
            
            VecRestoreSubVector(eta, isSurface, &etaActually);
            eta.writeTimeLevelData(0, 0);
            phi.writeTimeLevelData(0, 1);
            //writeFunc.write(meshes_[0], "solution at time " + std::to_string(t), false, this);
            VecGetSubVector(phi, isSurface, &phiS);
            VecGetSubVector(phi, isRest, &phiOther);
            VecGetSubVector(eta, isSurface, &etaActually);


                MatMult(surfaceMass, etaActually, surfaceRHS);
                VecScale(surfaceRHS, g);
                MatMultAdd(surfaceDMass, phiS, surfaceRHS, surfaceRHS);
                MatMultAdd(surfaceB, phiS, surfaceRHS, surfaceRHS);
                MatMultAdd(interiorB, phiOther, surfaceRHS, surfaceRHS);
                VecScale(surfaceRHS, 2);
                MatMultAdd(interiorC, phiOther, surfaceRHS, surfaceRHS);
                MatMultAdd(surfaceC, phiS, surfaceRHS, surfaceRHS);
                VecScale(surfaceRHS, -dt / 4.);
                MatMultAdd(surfaceMass, phiS, surfaceRHS, surfaceRHS);

                t += dt / 2.;
                moveNodes();
                constructM();
                //dMdt.reset();
                M.reset();
                MatGetSubMatrix(M, isSurface, isSurface, MAT_REUSE_MATRIX, &surfaceMass);
                //MatGetSubMatrix(dMdt, isSurface, isSurface, MAT_REUSE_MATRIX, &surfaceDMdt);
                //MatAXPY(M,dt/2.,dMdt,SAME_NONZERO_PATTERN);
                
                KSPSetOperators(surface, surfaceMass, surfaceMass);
                KSPSolve(surface, surfaceRHS, phiS);
                KSPGetConvergedReason(surface, &conferge);
                KSPGetIterationNumber(surface, &iterations);
                std::cout << "updating \\phi: KSP solver ended because of " << KSPConvergedReasons[conferge] << " in " << iterations << " iterations." << std::endl;
                first = false;
            }
            VecRestoreSubVector(phi, isSurface, &phiS);
            VecRestoreSubVector(phi, isRest, &phiOther);
            VecRestoreSubVector(eta, isSurface, &etaActually);
            eta.writeTimeLevelData(0, 0);
            phi.writeTimeLevelData(0, 1);
            writeFunc.write(meshes_[0], "solution at time " + std::to_string(t), false, this);
            printError();
            computeEnergy();
        }
        return true;
    }

    bool solveSymplecticEuler() {
        bool first(true);
        std::cout.precision(14);
        moveNodes();
        constructM();
        readInitialConditions();

        Utilities::GlobalPetscVector eta(meshes_[0]), phi(meshes_[0]);
        Utilities::GlobalPetscMatrix M(meshes_[1], -1, 0), D(meshes_[1], -1, 1);
        Vec phiS, phiOther, etaActually, interiorRHS, surfaceRHS, surfaceOld, surfaceUpdate, interiorUpdate;
        Mat surfaceMass, interiorStifness, surfaceStifness, mixStifness, backStifness, surfaceDMass;
        Mat surfaceC, interiorC, surfaceB, interiorB;

        IS isSurface, isRest;
        getSurfaceIS(M, &isSurface, &isRest);

        std::ofstream outFile("output.dat");
        Output::TecplotDiscontinuousSolutionWriter writeFunc(outFile, "test", "01", "eta_num,phi_num,phi_x,phi_y,eta_exact,phi_exact");

        //deal with initial conditions
        double g(1), dt(0.00125);
        eta.constructFromTimeLevelData(0, 0);
        phi.constructFromTimeLevelData(0, 1);

        int numberOfSnapshots(1001);
        int numberOfTimeSteps(80);
        //int numberOfTimeSteps(16); //placeholder parameters
        //int numberOfSnapshots(1 + n_ / 32. * 320);
        //double dt(0.01551125);

        MatGetSubMatrix(M, isSurface, isSurface, MAT_INITIAL_MATRIX, &surfaceMass);
        MatGetSubMatrix(D, isSurface, isSurface, MAT_INITIAL_MATRIX, &surfaceDMass);

        KSP interior, surface;
        KSPCreate(MPI_COMM_WORLD, &surface);
        KSPSetOperators(surface, surfaceMass, surfaceMass);
        KSPSetFromOptions(surface);
        KSPConvergedReason conferge;
        int iterations;

        VecGetSubVector(phi, isSurface, &phiS);
        VecGetSubVector(phi, isRest, &phiOther);
        VecGetSubVector(eta, isSurface, &etaActually);
        VecDuplicate(phiS, &surfaceRHS);
        VecDuplicate(phiS, &surfaceOld);
        VecDuplicate(phiS, &surfaceUpdate);
        VecDuplicate(phiOther, &interiorRHS);
        VecDuplicate(phiOther, &interiorUpdate);
        VecCopy(etaActually, surfaceRHS);
        //KSPSolve(surface,surfaceRHS,etaActually);
        //KSPGetConvergedReason(surface,&conferge);
        //KSPGetIterationNumber(surface,&iterations);
        //std::cout<<"Finalizing interpolation (1): KSP solver ended because of "<<KSPConvergedReasons[conferge]<<" in "<<iterations<<" iterations."<<std::endl;
        VecRestoreSubVector(eta, isSurface, &etaActually);
        VecRestoreSubVector(phi, isSurface, &phiS);
        VecRestoreSubVector(phi, isRest, &phiOther);


        eta.writeTimeLevelData(0, 0);
        moveNodes();
        readInitialConditions();
        eta.constructFromTimeLevelData(0, 0);
        CHKERRQ(VecGetSubVector(eta, isSurface, &etaActually));
        VecCopy(etaActually, surfaceRHS);
        KSPSetOperators(surface, surfaceMass, surfaceMass);
        //KSPSolve(surface,surfaceRHS,etaActually);
        //KSPGetConvergedReason(surface,&conferge);
        //KSPGetIterationNumber(surface,&iterations);
        //std::cout<<"Finalizing interpolation (1): KSP solver ended because of "<<KSPConvergedReasons[conferge]<<" in "<<iterations<<" iterations."<<std::endl;
        VecRestoreSubVector(eta, isSurface, &etaActually);
        eta.writeTimeLevelData(0, 0);

        //accurate construction of phi
        phi.constructFromTimeLevelData(0, 1);
        CHKERRQ(VecGetSubVector(phi, isSurface, &phiS));
        VecGetSubVector(eta, isSurface, &etaActually);
        VecCopy(phiS, surfaceRHS);
        //KSPSolve(surface,surfaceRHS,phiS);
        //KSPGetConvergedReason(surface,&conferge);
        //KSPGetIterationNumber(surface,&iterations);
        //std::cout<<"Finalizing interpolation (2): KSP solver ended because of "<<KSPConvergedReasons[conferge]<<" in "<<iterations<<" iterations."<<std::endl;
        VecRestoreSubVector(phi, isSurface, &phiS);
        phi.writeTimeLevelData(0, 1);

        constructM();
        M.reset();
        D.reset();
        MatGetSubMatrix(M, isSurface, isSurface, MAT_REUSE_MATRIX, &surfaceMass);
        MatGetSubMatrix(D, isSurface, isSurface, MAT_REUSE_MATRIX, &surfaceDMass);

        constructS();
        constructC();

        Vec scaleVec;
        VecDuplicate(phi, &scaleVec);
        VecSet(scaleVec, 1.);
        VecGetSubVector(scaleVec, isSurface, &surfaceUpdate);
        VecSet(surfaceUpdate, dt / 2.);
        VecRestoreSubVector(scaleVec, isSurface, &surfaceUpdate);

        Utilities::GlobalPetscMatrix S(meshes_[0], 0, -1), SquareC(meshes_[0], 1, -1), squareB(meshes_[0], -1, 0);
        Utilities::GlobalPetscVector W(meshes_[0], -1, 0);
        W.assemble();
        Mat C, CollectCulumns, B;
        Vec interiorW, surfaceW;
        MatCreate(MPI_COMM_WORLD, &CollectCulumns);
        int n;
        MatGetSize(SquareC, &n, NULL); //extra rows hardly matter because they are empty anyway
        MatSetSizes(CollectCulumns, PETSC_DECIDE, PETSC_DECIDE, n, n);
        MatSetUp(CollectCulumns);
        for (std::pair<int, int> pair : surfaceNodesAboveANode_) {
            if ((pair.second + 1) / (n_ + 1)*(n_ + 1) == (pair.second + 1)) {
                //MatSetValue(CollectCulumns,0,pair.first-n_-pair.first/(n_+1),1.,INSERT_VALUES);
            } else {
                int first(pair.first - pair.first / (n_ + 1)), second(pair.second - pair.second / (n_ + 1));
                if (DIM == 3) {
                    if ((second + 1) / 5 * 5 == second + 1) {

                    } else {
                        MatSetValue(CollectCulumns, second - second / 5, first - first / 5, 1., INSERT_VALUES);
                    }
                } else {
                    MatSetValue(CollectCulumns, second, first, 1., INSERT_VALUES);
                }
            }
        }


        MatAssemblyBegin(CollectCulumns, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(CollectCulumns, MAT_FINAL_ASSEMBLY);
        MatMatMult(CollectCulumns, SquareC, MAT_INITIAL_MATRIX, 1., &C);
        MatMatMult(CollectCulumns, squareB, MAT_INITIAL_MATRIX, 1., &B);

        MatGetSubMatrix(C, isSurface, isSurface, MAT_INITIAL_MATRIX, &surfaceC);
        MatGetSubMatrix(C, isSurface, isRest, MAT_INITIAL_MATRIX, &interiorC);
        MatGetSubMatrix(B, isSurface, isSurface, MAT_INITIAL_MATRIX, &surfaceB);
        MatGetSubMatrix(B, isSurface, isRest, MAT_INITIAL_MATRIX, &interiorB);
        MatGetSubMatrix(S, isRest, isRest, MAT_INITIAL_MATRIX, &interiorStifness);
        MatGetSubMatrix(S, isSurface, isSurface, MAT_INITIAL_MATRIX, &surfaceStifness);
        MatGetSubMatrix(S, isRest, isSurface, MAT_INITIAL_MATRIX, &mixStifness);
        MatGetSubMatrix(S, isSurface, isRest, MAT_INITIAL_MATRIX, &backStifness);
        VecGetSubVector(W, isSurface, &surfaceW);
        VecGetSubVector(W, isRest, &interiorW);

        Mat BigMatrix;
        MatDuplicate(S, MAT_DO_NOT_COPY_VALUES, &BigMatrix);
        Vec BigVector, BigRHS;
        VecDuplicate(phi, &BigVector);
        VecDuplicate(phi, &BigRHS);
        VecGetSubVector(BigRHS, isSurface, &surfaceRHS);
        VecGetSubVector(BigRHS, isRest, &interiorRHS);
        KSPCreate(MPI_COMM_WORLD, &interior);
        KSPSetOperators(interior, BigMatrix, BigMatrix);
        MatSetOption(BigMatrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
        MatSetOption(BigMatrix, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
        KSPSetFromOptions(interior);

        double accuracy;

        KSPGetTolerances(interior, &accuracy, NULL, NULL, NULL);
        KSPSetTolerances(interior, 1e-12, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);

        writeFunc.write(meshes_[0], "solution at time 0", false, this);
        printError();
        //computeEnergy();//does not make sense since phi (internal) is not calculated yet
        for (int i = 1; i < numberOfSnapshots; ++i) {
            VecGetSubVector(phi, isSurface, &phiS);
            VecGetSubVector(phi, isRest, &phiOther);
            VecGetSubVector(eta, isSurface, &etaActually);
            for (int j = 0; j < numberOfTimeSteps; ++j) {

                t += dt;
                moveNodes();
                constructM();
                M.reset();
                MatGetSubMatrix(M, isSurface, isSurface, MAT_REUSE_MATRIX, &surfaceMass);
                VecRestoreSubVector(phi, isSurface, &phiS);
                VecRestoreSubVector(phi, isRest, &phiOther);
                phi.writeTimeLevelData(0, 1);
                constructC();
                constructS();
                SquareC.reset();
                S.reset();
                MatGetSubMatrix(S, isRest, isRest, MAT_REUSE_MATRIX, &interiorStifness);
                MatGetSubMatrix(S, isSurface, isSurface, MAT_REUSE_MATRIX, &surfaceStifness);
                MatGetSubMatrix(S, isRest, isSurface, MAT_REUSE_MATRIX, &mixStifness);
                MatGetSubMatrix(S, isSurface, isRest, MAT_REUSE_MATRIX, &backStifness);
                MatMatMult(CollectCulumns, SquareC, MAT_REUSE_MATRIX, 1., &C);
                MatGetSubMatrix(C, isSurface, isSurface, MAT_REUSE_MATRIX, &surfaceC);
                MatGetSubMatrix(C, isSurface, isRest, MAT_REUSE_MATRIX, &interiorC);
                VecGetSubVector(phi, isSurface, &phiS);
                VecGetSubVector(phi, isRest, &phiOther);

                MatMult(surfaceMass, etaActually, surfaceOld);
                VecScale(surfaceOld, -g * dt);
                MatMultAdd(surfaceMass, phiS, surfaceOld, surfaceOld); //move to beginning of time loop (use old M)

                MatMult(surfaceC, phiS, surfaceRHS);
                MatMultAdd(interiorC, phiOther, surfaceRHS, surfaceRHS);
                VecScale(surfaceRHS, dt / 2.);
                MatMultAdd(surfaceMass, phiS, surfaceRHS, surfaceRHS);
                VecAYPX(surfaceRHS, -1., surfaceOld);
                MatMult(interiorStifness, phiOther, interiorRHS);
                MatMultAdd(mixStifness, phiS, interiorRHS, interiorRHS);
                double normSurface, normInterior;
                VecNorm(surfaceRHS, NORM_2, &normSurface);
                VecNorm(interiorRHS, NORM_2, &normInterior);
                int numberOfNewtonSteps(0);
                while (normInterior + normSurface > accuracy) {
                    numberOfNewtonSteps++;
                    //std::cout<<"In Newton iteration with solution errors: "<<normInterior<<" and "<<normSurface<<std::endl;

                    MatCopy(S, BigMatrix, DIFFERENT_NONZERO_PATTERN);
                    MatZeroRowsIS(BigMatrix, isSurface, 0., NULL, NULL);

                    CHKERRABORT(MPI_COMM_WORLD, MatAXPY(BigMatrix, -1., M, DIFFERENT_NONZERO_PATTERN));
                    CHKERRABORT(MPI_COMM_WORLD, MatAXPY(BigMatrix, -dt, C, DIFFERENT_NONZERO_PATTERN));

                    //MatView(BigMatrix,PETSC_VIEWER_DRAW_WORLD);

                    VecRestoreSubVector(BigRHS, isSurface, &surfaceRHS);
                    VecRestoreSubVector(BigRHS, isRest, &interiorRHS);
                    if (first)
                        KSPSetOperators(interior, BigMatrix, BigMatrix);
                    else
                        KSPSetOperators(interior, BigMatrix, BigMatrix);
                    KSPSolve(interior, BigRHS, BigVector);
                    KSPGetConvergedReason(interior, &conferge);
                    KSPGetIterationNumber(interior, &iterations);
                    //std::cout<<"KSP solver ended because of "<<KSPConvergedReasons[conferge]<<" in "<<iterations<<" iterations."<<std::endl;
                    VecGetSubVector(BigVector, isSurface, &surfaceUpdate);
                    VecGetSubVector(BigVector, isRest, &interiorUpdate);
                    VecGetSubVector(BigRHS, isSurface, &surfaceRHS);
                    VecGetSubVector(BigRHS, isRest, &interiorRHS);
                    VecAXPY(phiS, -1.0, surfaceUpdate);
                    VecAXPY(phiOther, -1.0, interiorUpdate);
                    VecRestoreSubVector(BigVector, isSurface, &surfaceUpdate);
                    VecRestoreSubVector(BigVector, isRest, &interiorUpdate);

                    VecRestoreSubVector(phi, isSurface, &phiS);
                    VecRestoreSubVector(phi, isRest, &phiOther);
                    phi.writeTimeLevelData(0, 1);
                    constructC();
                    SquareC.reset();
                    MatMatMult(CollectCulumns, SquareC, MAT_REUSE_MATRIX, 1., &C);
                    MatGetSubMatrix(C, isSurface, isSurface, MAT_REUSE_MATRIX, &surfaceC);
                    MatGetSubMatrix(C, isSurface, isRest, MAT_REUSE_MATRIX, &interiorC);
                    VecGetSubVector(phi, isSurface, &phiS);
                    VecGetSubVector(phi, isRest, &phiOther);

                    MatMult(surfaceC, phiS, surfaceRHS);
                    MatMultAdd(interiorC, phiOther, surfaceRHS, surfaceRHS);
                    VecScale(surfaceRHS, dt / 2.);
                    MatMultAdd(surfaceMass, phiS, surfaceRHS, surfaceRHS);
                    VecAYPX(surfaceRHS, -1., surfaceOld);
                    MatMult(interiorStifness, phiOther, interiorRHS);
                    MatMultAdd(mixStifness, phiS, interiorRHS, interiorRHS);
                    VecNorm(surfaceRHS, NORM_2, &normSurface);
                    VecNorm(interiorRHS, NORM_2, &normInterior);
                }
                std::cout << "First set of equations solved after " << numberOfNewtonSteps << " newton iterations" << std::endl;

                MatMult(surfaceStifness, phiS, surfaceRHS);
                MatMultAdd(backStifness, phiOther, surfaceRHS, surfaceRHS);
                VecScale(surfaceRHS, dt / g);
                MatMultAdd(surfaceMass, etaActually, surfaceRHS, surfaceRHS);
                VecScale(surfaceRHS, g);

                KSPSetOperators(surface, surfaceMass, surfaceMass);
                KSPSolve(surface, surfaceRHS, etaActually);
                KSPGetConvergedReason(surface, &conferge);
                KSPGetIterationNumber(surface, &iterations);
                std::cout << "updating \\eta: KSP solver ended because of " << KSPConvergedReasons[conferge] << " in " << iterations << " iterations." << std::endl;
                first = false;
                VecRestoreSubVector(phi, isSurface, &phiS);
                VecRestoreSubVector(phi, isRest, &phiOther);
                VecRestoreSubVector(eta, isSurface, &etaActually);
                eta.writeTimeLevelData(0, 0);
                phi.writeTimeLevelData(0, 1);
                VecGetSubVector(phi, isSurface, &phiS);
                VecGetSubVector(phi, isRest, &phiOther);
                VecGetSubVector(eta, isSurface, &etaActually);
            }
            VecRestoreSubVector(phi, isSurface, &phiS);
            VecRestoreSubVector(phi, isRest, &phiOther);
            VecRestoreSubVector(eta, isSurface, &etaActually);
            eta.writeTimeLevelData(0, 0);
            phi.writeTimeLevelData(0, 1);
            writeFunc.write(meshes_[0], "solution at time t", false, this);
            printError();
            computeEnergy();
        }
        return true;
    }

    void testMatrices() {
        int iterations;
        KSPConvergedReason conferge;
        readInitialConditions();
        Utilities::GlobalPetscVector eta(meshes_[0]), phi(meshes_[0]);
        for (int i = 0; i < 3; ++i) {
            readInitialConditions();
            eta.constructFromTimeLevelData(0, 0);
            eta.writeTimeLevelData(0, 0);
            moveNodes();
        }
        Vec RHS;
        phi.constructFromTimeLevelData(0, 1);
        VecDuplicate(phi, &RHS);
        IS isSurface, isRest;
        constructS();
        Utilities::GlobalPetscMatrix S(meshes_[0], 0, -1);
        getSurfaceIS(S, &isSurface, &isRest);
        MatZeroRowsIS(S, isSurface, 1.0, phi, RHS);
        KSP laplace;
        KSPCreate(MPI_COMM_WORLD, &laplace);
        KSPSetOperators(laplace, S, S);
        KSPSetFromOptions(laplace);
        KSPSolve(laplace, RHS, phi);
        KSPGetConvergedReason(laplace, &conferge);
        KSPGetIterationNumber(laplace, &iterations);
        std::cout << "updating \\phi: KSP solver ended because of " << KSPConvergedReasons[conferge] << " in " << iterations << " iterations." << std::endl;
        phi.writeTimeLevelData(0, 1);
        std::ofstream outFile("output.dat");
        Output::TecplotDiscontinuousSolutionWriter writeFunc(outFile, "test", "01", "eta_num,phi_num,eta_exact,phi_exact");
        writeFunc.write(meshes_[0], "solution at time 0", false, this);
        t = 1e-12;
        Vec etaSurface, phiSurface, phiRest, surfaceDuplicate, restDuplicate;
        Mat mass, dMass, interiorS, dInteriorS, surfaceS, dSurfaceS, mixS, dMixS, backS, dBackS;
        Mat testM2, testInteriorC, testSurfaceC, testMixC, testBackC, testB;

        //construct analytical derivatives
        constructS();
        constructC();
        constructM();
        Utilities::GlobalPetscMatrix M(meshes_[1], -1, 0), M2(meshes_[1], -1, 1), SquareC(meshes_[0], 1, -1), SquareC2(meshes_[0], 2, -1), squareB(meshes_[0], -1, 0);
        Utilities::GlobalPetscVector W(meshes_[0], -1, 0);
        Mat C, C2, B, CollectCulumns;
        Mat surfaceC, interiorC, mixC, backC, surfaceM2;
        Vec dW;
        MatCreate(MPI_COMM_WORLD, &CollectCulumns);
        int n;
        MatGetSize(SquareC, &n, NULL); //extra rows hardly matter because they are empty anyway
        MatSetSizes(CollectCulumns, PETSC_DECIDE, PETSC_DECIDE, n, n);
        MatSetUp(CollectCulumns);

        t -= 1e-6;
        moveNodes();
        constructM();
        M.reset();
        //MatView(M,PETSC_VIEWER_DRAW_WORLD);
        MatDuplicate(M, MAT_COPY_VALUES, &mass);
        t += 2e-6;
        moveNodes();
        constructM();
        M.reset();
        MatAYPX(mass, -1., M, DIFFERENT_NONZERO_PATTERN);
        t -= 1e-6;
        MatScale(mass, 5e5);
        MatView(mass, PETSC_VIEWER_STDOUT_WORLD);
        MatView(M2, PETSC_VIEWER_STDOUT_WORLD);

        MatAXPY(mass, -1., M2, DIFFERENT_NONZERO_PATTERN);
        MatChop(mass, 1e-6);
        MatView(mass, PETSC_VIEWER_STDOUT_WORLD);
        MatView(mass, PETSC_VIEWER_DRAW_WORLD);



        for (std::pair<int, int> pair : surfaceNodesAboveANode_) {
            //if((pair.second+1)/(n_+1)*(n_+1)==(pair.second+1)){
            //MatSetValue(CollectCulumns,0,pair.first-n_-pair.first/(n_+1),1.,INSERT_VALUES);
            //}else{
            MatSetValue(CollectCulumns, pair.second/*-pair.second/(n_+1)*/, pair.first/*-pair.first/(n_+1)*/, 1., INSERT_VALUES);
            //}
        }
        MatAssemblyBegin(CollectCulumns, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(CollectCulumns, MAT_FINAL_ASSEMBLY);
        MatMatMult(CollectCulumns, SquareC, MAT_INITIAL_MATRIX, 1., &C);
        MatMatMult(CollectCulumns, SquareC2, MAT_INITIAL_MATRIX, 1., &C2);
        MatMatMult(CollectCulumns, squareB, MAT_INITIAL_MATRIX, 1., &B);
        MatGetSubMatrix(M2, isSurface, isSurface, MAT_INITIAL_MATRIX, &surfaceM2);
        MatGetSubMatrix(C, isSurface, isSurface, MAT_INITIAL_MATRIX, &surfaceC);
        MatGetSubMatrix(C, isSurface, isRest, MAT_INITIAL_MATRIX, &interiorC);
        MatGetSubMatrix(C2, isSurface, isSurface, MAT_INITIAL_MATRIX, &mixC);
        MatGetSubMatrix(C2, isSurface, isRest, MAT_INITIAL_MATRIX, &backC);

        MatDuplicate(surfaceM2, MAT_DO_NOT_COPY_VALUES, &testM2);
        MatDuplicate(surfaceC, MAT_DO_NOT_COPY_VALUES, &testSurfaceC);
        MatDuplicate(interiorC, MAT_DO_NOT_COPY_VALUES, &testInteriorC);
        MatDuplicate(mixC, MAT_DO_NOT_COPY_VALUES, &testBackC);
        MatDuplicate(backC, MAT_DO_NOT_COPY_VALUES, &testMixC);
        MatDuplicate(B, MAT_DO_NOT_COPY_VALUES, &testB);
        MatSetOption(testM2, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
        MatSetOption(testInteriorC, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
        MatSetOption(testSurfaceC, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
        MatSetOption(testMixC, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
        MatSetOption(testBackC, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
        MatSetOption(testB, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
        VecGetSubVector(eta, isSurface, &etaSurface);
        VecGetSubVector(phi, isSurface, &phiSurface);
        VecGetSubVector(phi, isRest, &phiRest);
        VecDuplicate(etaSurface, &surfaceDuplicate);
        VecDuplicate(phiRest, &restDuplicate);
        //construct numerical derivatives
        VecDuplicate(W, &dW);
        for (int i = 0; i < n_; ++i) {
            //high value
            VecSetValue(etaSurface, i, 1e-6, ADD_VALUES);
            VecAssemblyBegin(etaSurface);
            VecAssemblyEnd(etaSurface);
            VecRestoreSubVector(eta, isSurface, &etaSurface);
            eta.writeTimeLevelData(0, 0);
            moveNodes();
            writeFunc.write(meshes_[0], "solution at time 0", false, this);
            VecGetSubVector(eta, isSurface, &etaSurface);
            constructM();
            constructS();
            M.reset();
            S.reAssemble();
            W.assemble();
            VecCopy(W, dW);
            MatGetSubMatrix(M, isSurface, isSurface, MAT_INITIAL_MATRIX, &dMass);
            MatGetSubMatrix(S, isSurface, isSurface, MAT_INITIAL_MATRIX, &dSurfaceS);
            MatGetSubMatrix(S, isRest, isRest, MAT_INITIAL_MATRIX, &dInteriorS);
            MatGetSubMatrix(S, isRest, isSurface, MAT_INITIAL_MATRIX, &dMixS);
            MatGetSubMatrix(S, isSurface, isRest, MAT_INITIAL_MATRIX, &dBackS);

            //low value
            VecSetValue(etaSurface, i, -2e-6, ADD_VALUES);
            VecAssemblyBegin(etaSurface);
            VecAssemblyEnd(etaSurface);
            VecRestoreSubVector(eta, isSurface, &etaSurface);
            eta.writeTimeLevelData(0, 0);
            moveNodes();
            writeFunc.write(meshes_[0], "solution at time 0", false, this);
            VecGetSubVector(eta, isSurface, &etaSurface);
            constructM();
            constructS();
            M.reset();
            S.reAssemble();
            W.assemble();
            MatGetSubMatrix(M, isSurface, isSurface, MAT_INITIAL_MATRIX, &mass);
            MatGetSubMatrix(S, isSurface, isSurface, MAT_INITIAL_MATRIX, &surfaceS);
            MatGetSubMatrix(S, isRest, isRest, MAT_INITIAL_MATRIX, &interiorS);
            MatGetSubMatrix(S, isRest, isSurface, MAT_INITIAL_MATRIX, &mixS);
            MatGetSubMatrix(S, isSurface, isRest, MAT_INITIAL_MATRIX, &backS);
            //MatView(mass,PETSC_VIEWER_DRAW_WORLD);
            //MatView(dMass,PETSC_VIEWER_DRAW_WORLD);
            MatAXPY(dMass, -1., mass, SAME_NONZERO_PATTERN);
            MatAXPY(dSurfaceS, -1., surfaceS, SAME_NONZERO_PATTERN);
            //MatView(dSurfaceS,PETSC_VIEWER_DRAW_WORLD);
            MatAXPY(dInteriorS, -1., interiorS, SAME_NONZERO_PATTERN);
            MatAXPY(dMixS, -1., mixS, SAME_NONZERO_PATTERN);
            MatAXPY(dBackS, -1., backS, SAME_NONZERO_PATTERN);
            VecAXPY(dW, -1, W);
            MatScale(dMass, 5e5);
            MatScale(dSurfaceS, 5e5);
            MatScale(dInteriorS, 5e5);
            MatScale(dMixS, 5e5);
            MatScale(dBackS, 5e5);
            VecScale(dW, 5e5);


            //multiply by extra vector and store
            VecSetValue(etaSurface, i, 1e-6, ADD_VALUES);
            VecAssemblyBegin(etaSurface);
            VecAssemblyEnd(etaSurface);
            MatMult(dMass, etaSurface, surfaceDuplicate);
            double *vecData;
            VecGetArray(surfaceDuplicate, &vecData);
            for (int j = 0; j < n_; ++j) {
                if (std::fabs(vecData[j]) > 1e-9) {
                    MatSetValue(testM2, i, j, vecData[j], INSERT_VALUES);
                }
            }
            VecRestoreArray(surfaceDuplicate, &vecData);
            MatMult(dSurfaceS, phiSurface, surfaceDuplicate);
            VecGetArray(surfaceDuplicate, &vecData);
            for (int j = 0; j < n_; ++j) {
                if (std::fabs(vecData[j]) > 1e-9) {
                    MatSetValue(testSurfaceC, i, j, vecData[j], INSERT_VALUES);
                }
            }
            VecRestoreArray(surfaceDuplicate, &vecData);
            MatMult(dBackS, phiRest, surfaceDuplicate);
            VecGetArray(surfaceDuplicate, &vecData);
            for (int j = 0; j < n_; ++j) {
                if (std::fabs(vecData[j]) > 1e-9) {
                    MatSetValue(testBackC, i, j, vecData[j], INSERT_VALUES);
                }
            }
            VecRestoreArray(surfaceDuplicate, &vecData);
            MatMult(dMixS, phiSurface, restDuplicate);
            VecGetArray(restDuplicate, &vecData);
            for (int j = 0; j < n_ * n_ / 8; ++j) {
                if (std::fabs(vecData[j]) > 1e-9) {
                    MatSetValue(testMixC, i, j, vecData[j], INSERT_VALUES);
                }
            }
            VecRestoreArray(restDuplicate, &vecData);
            MatMult(dInteriorS, phiRest, restDuplicate);
            VecGetArray(restDuplicate, &vecData);
            for (int j = 0; j < n_ * n_ / 8; ++j) {
                if (std::fabs(vecData[j]) > 1e-9) {
                    MatSetValue(testInteriorC, i, j, vecData[j], INSERT_VALUES);
                }
            }
            VecRestoreArray(restDuplicate, &vecData);
            VecGetArray(dW, &vecData);
            for (int j = 0; j < n_ * (n_ / 8 + 1); ++j) {
                if (std::fabs(vecData[j]) > 1e-9) {
                    MatSetValue(testB, i, j, vecData[j], INSERT_VALUES);
                }
            }
            VecRestoreArray(dW, &vecData);
        }
        //check for correctness
        MatAssemblyBegin(testM2, MAT_FINAL_ASSEMBLY);
        MatAssemblyBegin(testBackC, MAT_FINAL_ASSEMBLY);
        MatAssemblyBegin(testInteriorC, MAT_FINAL_ASSEMBLY);
        MatAssemblyBegin(testMixC, MAT_FINAL_ASSEMBLY);
        MatAssemblyBegin(testSurfaceC, MAT_FINAL_ASSEMBLY);
        MatAssemblyBegin(testB, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(testM2, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(testBackC, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(testInteriorC, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(testMixC, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(testSurfaceC, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(testB, MAT_FINAL_ASSEMBLY);
        MatView(testB, PETSC_VIEWER_STDOUT_WORLD);

        MatAXPY(testM2, -1, surfaceM2, SAME_NONZERO_PATTERN);
        MatAXPY(testBackC, -1, mixC, SAME_NONZERO_PATTERN);
        MatAXPY(testInteriorC, -1, interiorC, SAME_NONZERO_PATTERN);
        MatAXPY(testMixC, -1, backC, SAME_NONZERO_PATTERN);
        MatAXPY(testSurfaceC, -1, surfaceC, SAME_NONZERO_PATTERN);
        MatAXPY(testB, -1, B, SAME_NONZERO_PATTERN);
        MatChop(testM2, 1e-6);
        MatChop(testBackC, 1e-6);
        MatChop(testInteriorC, 1e-6);
        MatChop(testMixC, 1e-6);
        MatChop(testSurfaceC, 1e-6);
        MatChop(testB, 1e-6);
        //MatView(testM2,PETSC_VIEWER_STDOUT_WORLD);
        //MatView(testSurfaceC,PETSC_VIEWER_STDOUT_WORLD);
        //MatView(testInteriorC,PETSC_VIEWER_STDOUT_WORLD);
        MatView(testB, PETSC_VIEWER_STDOUT_WORLD);
        //MatView(testBackC,PETSC_VIEWER_STDOUT_WORLD);
        MatView(B, PETSC_VIEWER_STDOUT_WORLD);
        MatView(testB, PETSC_VIEWER_DRAW_WORLD);
        MatView(testM2, PETSC_VIEWER_DRAW_WORLD);
        MatView(testSurfaceC, PETSC_VIEWER_DRAW_WORLD);
        MatView(testInteriorC, PETSC_VIEWER_DRAW_WORLD);
        MatView(testMixC, PETSC_VIEWER_DRAW_WORLD);
        MatView(testBackC, PETSC_VIEWER_DRAW_WORLD);
    }


private:

    static std::map<int, std::vector<Geometry::PointPhysical*> > nodesBelowASurfaceNode_;
    std::map<int, int> surfaceNodesAboveANode_;

    //number of elements per cardinal direction
    static int n_;

    //polynomial order of the approximation
    int p_;
    static double t;
};

double DGWave::t = 0;
int DGWave::n_;
std::map<int, std::vector<Geometry::PointPhysical*> > DGWave::nodesBelowASurfaceNode_;

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
            throw "usage: LinearPotentialFlow.out n p [petsc-args]";
        }
        PetscInitialize(&argc, &argv, NULL, NULL);
        DGWave test(n, p);
        test.initialise();
        test.solveStormerVerlet();
        PetscFinalize();
        return 0;
    } catch (const char* e) {
        std::cout << e << std::endl;
    }
}





