/*
 * main.cpp
 *
 *  Created on: Feb 28, 2014
 *      Author: brinkf
 */

#define hpGEM_INCLUDE_PETSC_SUPPORT

#include "Base/HpgemAPIBase.h"
#include "Base/L2Norm.h"
#include "Output/TecplotSingleElementWriter.h"
#include <fstream>
#include "Utilities/BasisFunctions2DH1ConformingSquare.h"
#include "Integration/ElementIntegrandBase.h"
#include "Integration/FaceIntegrandBase.h"
#include "Integration/FaceIntegral.h"
#include "Integration/ElementIntegral.h"
#include "Base/GlobalData.h"
#include "Base/RectangularMeshDescriptor.h"
#include "Base/ConfigurationData.h"
#include "Base/Element.h"
#include "Base/Face.h"
#include "Output/TecplotDiscontinuousSolutionWriter.h"
#include "Base/ElementCacheData.h"
#include "Base/FaceCacheData.h"
#include <cmath>

#include "Utilities/GlobalMatrix.h"
#include "Utilities/GlobalVector.h"
#include "petscksp.h"

const unsigned int DIM = 2;

class DGWave : public Base::HpgemAPIBase, public Output::TecplotSingleElementWriter
{
public:
    DGWave(int n, int p)
            : HpgemAPIBase(new Base::GlobalData(), new Base::ConfigurationData(DIM, 2, p)), n_(n), p_(p)
    {
    }
    
    ///set up the mesh
    bool initialise()
    {
        Base::RectangularMeshDescriptor description(DIM);
        for (int i = 0; i < DIM; ++i)
        {
            description.bottomLeft_[i] = 0;
            description.topRight_[i] = 1;
            description.numElementsInDIM_[i] = n_;
            description.boundaryConditions_[i] = Base::BoundaryType::PERIODIC;
        }
        description.topRight_[DIM - 1] = -1;
        description.numElementsInDIM_[DIM - 1] /= 8;
        description.boundaryConditions_[DIM - 1] = Base::BoundaryType::SOLID_WALL;
        addMesh(description, Base::MeshType::RECTANGULAR, 1, 1, 1, 1);
        meshes_[0]->setDefaultBasisFunctionSet(Utilities::createInteriorBasisFunctionSet2DH1Square(p_));
        std::vector<const Base::BasisFunctionSet*> bFsets;
        bFsets = Utilities::createVertexBasisFunctionSet2DH1Square(p_);
        meshes_[0]->addVertexBasisFunctionSet(bFsets);
        std::vector<const Base::OrientedBasisFunctionSet*> oBFsets;
        oBFsets = Utilities::createFaceBasisFunctionSet2DH1Square(p_);
        meshes_[0]->addFaceBasisFunctionSet(oBFsets);
        return true;
    }
    
    class : public Integration::ElementIntegrandBase<LinearAlgebra::Matrix>
    {
    public:
        void elementIntegrand(const Base::Element* element, const Geometry::PointReference& p, LinearAlgebra::Matrix& ret) override final
        {
            int numBasisFuns = element->getNrOfBasisFunctions();
            ret.resize(numBasisFuns, numBasisFuns);
            for (int i = 0; i < numBasisFuns; ++i)
            {
                for (int j = 0; j <= i; ++j)
                {
                    ret(i, j) = ret(j, i) = element->basisFunctionDeriv(i, p) * element->basisFunctionDeriv(j, p);
                }
            }
        }
    } stifnessIntegrand;

    class : public Integration::FaceIntegrandBase<LinearAlgebra::Matrix>
    {
    public:
        void faceIntegrand(const Base::Face* face, const LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p, LinearAlgebra::Matrix& ret) override final
        {
            int numBasisFuns = face->getNrOfBasisFunctions();
            ret.resize(numBasisFuns, numBasisFuns);
            for (int i = 0; i < numBasisFuns; ++i)
            {
                for (int j = 0; j < numBasisFuns; ++j)
                { //the basis functions belonging to internal parameters are 0 on the free surface anyway.
                    ret(i, j) = face->basisFunction(i, p) * face->basisFunction(j, p);
                }
            }
        }
    } massIntegrand;

    static void initialConditions(const PointPhysicalT& p, LinearAlgebra::NumericalVector& ret)
    {
        ret.resize(2);
        //ret[0]=cos(2*M_PI*p[0])*cosh(2*M_PI*(p[DIM-1]+1))/cosh(2*M_PI);//standing wave
        //ret[1]=0;
        double k = 2. * M_PI;
        ret[0] = cosh(k * (p[DIM - 1] + 1)) * sin(-k * p[0]) * sqrt(k * tanh(k)) / cosh(k) * 0.001; //moving wave
        ret[1] = cosh(k * (p[DIM - 1] + 1)) * cos(-k * p[0]) / cosh(k) * 0.001;
    }

    static void exactSolution(const double t, const PointPhysicalT& p, LinearAlgebra::NumericalVector& ret)
    {
        ret.resize(2);
        //ret[0]=cos(2*M_PI*p[0])*cos(sqrt(2*M_PI*tanh(2*M_PI))*t)*cosh(2*M_PI*(p[DIM-1]+1))/cosh(2*M_PI);//standing wave
        //ret[1]=cos(2*M_PI*p[0])*sin(-sqrt(2*M_PI*tanh(2*M_PI))*t)/sqrt(2*M_PI*tanh(2*M_PI))*cosh(2*M_PI*(p[DIM-1]+1))/cosh(2*M_PI);
        double k = 2. * M_PI;
        ret[0] = cosh(k * (p[DIM - 1] + 1)) * sin(sqrt(k * tanh(k)) * t - k * p[0]) * sqrt(k * tanh(k)) / cosh(k) * 0.001; //moving wave
        ret[1] = cosh(k * (p[DIM - 1] + 1)) * cos(sqrt(k * tanh(k)) * t - k * p[0]) / cosh(k) * 0.001;
        
    }
    
    class : public Integration::FaceIntegrandBase<LinearAlgebra::NumericalVector>
    {
    public:
        void faceIntegrand(const Base::Face* element, const LinearAlgebra::NumericalVector&, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret) override final
        {
            PointPhysicalT pPhys = element->referenceToPhysical(p);
            ret.resize(2 * element->getNrOfBasisFunctions());
            LinearAlgebra::NumericalVector data;
            initialConditions(pPhys, data);
            for (int i = 0; i < element->getNrOfBasisFunctions(); ++i)
            {
                ret(i) = element->basisFunction(i, p) * data[0];
                ret(i + element->getNrOfBasisFunctions()) = element->basisFunction(i, p) * data[1];
            }
        }
    } interpolator;

    class : public Integration::FaceIntegrandBase<LinearAlgebra::NumericalVector>
    {
    public:
        void faceIntegrand(const Base::Face* face, const LinearAlgebra::NumericalVector& normal, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret) override final
        {
            PointPhysicalT pPhys = face->referenceToPhysical(p);
            PointReferenceT pElement = face->mapRefFaceToRefElemL(p);
            ret.resize(2);
            static LinearAlgebra::NumericalVector exact(2);
            if (std::abs(pPhys[DIM - 1]) < 1e-9)
            {
                ret = face->getPtrElementLeft()->getSolution(0, pElement);
                exactSolution(t, pPhys, exact);
                ret -= exact;
                ret[0] *= ret[0];
                ret[1] *= ret[1];
            }
        }
    } error;

    class : public Integration::FaceIntegrandBase<LinearAlgebra::NumericalVector>
    {
    public:
        
        void faceIntegrand(const Base::Face* face, const LinearAlgebra::NumericalVector& normal, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret) override final
        {
            PointPhysicalT pPhys = face->referenceToPhysical(p);
            PointReferenceT pElement = face->mapRefFaceToRefElemL(p);
            ret.resize(1);
            static LinearAlgebra::NumericalVector dummySolution(2), gradPhi(DIM), temp(DIM);
            const LinearAlgebra::NumericalVector& expansioncoefficients = face->getPtrElementLeft()->getTimeLevelData(0, 0);
            if (std::abs(pPhys[DIM - 1]) < 1e-9)
            {
                dummySolution[0] = 0;
                for (int i = 0; i < face->getNrOfBasisFunctions(); ++i)
                {
                    dummySolution[0] += face->basisFunction(i, p) * expansioncoefficients(i);
                }
                ret[0] = dummySolution[0] * dummySolution[0];
                ret[0] /= 2; //assumes g=1
            }
        }
    } faceEnergy;

    class : public Integration::ElementIntegrandBase<LinearAlgebra::NumericalVector>
    {
    public:
        
        void elementIntegrand(const Base::Element* element, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret) override final
        {
            int numBasisFuns = element->getNrOfBasisFunctions();
            ret.resize(1);
            LinearAlgebra::NumericalVector gradPhi(DIM), temp(DIM);
            const LinearAlgebra::NumericalVector& expansioncoefficients = element->getTimeLevelData(0, 1);
            for (int i = 0; i < numBasisFuns; ++i)
            {
                temp = element->basisFunctionDeriv(i, p);
                temp *= expansioncoefficients(i);
                gradPhi += temp;
            }
            ret[0] = (gradPhi * gradPhi) / 2;
        }
    } elementEnergy;

    void writeToTecplotFile(const Base::Element* element, const PointReferenceT& p, std::ostream& out) override final
    {
        LinearAlgebra::NumericalVector value = element->getSolution(0, p);
        out << value[0] << " " << value[1];
        Geometry::PointPhysical pPhys = element->referenceToPhysical(p);
        exactSolution(t, pPhys, value);
        out << " " << value[0] << " " << value[1];
    }
    
    //assumes that 'all' means 'all relevant'; computes the mass matrix at the free surface
    void doAllFaceIntegration()
    {
        Integration::FaceIntegral integral(false);
        Geometry::PointReference p(DIM - 1);
        Geometry::PointPhysical pPhys(DIM);
        LinearAlgebra::Matrix result;
        LinearAlgebra::NumericalVector initialconditions;
        for (Base::Face* face : meshes_[0]->getFacesList())
        {
            p = face->getReferenceGeometry()->getCenter();
            pPhys = face->referenceToPhysical(p);
            if (std::abs(pPhys[DIM - 1]) < 1e-9)
            {
                result = integral.integrate(face, &massIntegrand);
                initialconditions = integral.integrate(face, &interpolator);
                face->getPtrElementLeft()->setTimeLevelData(0, initialconditions);
            }
            else
            {
                result.resize(face->getNrOfBasisFunctions(), face->getNrOfBasisFunctions());
                for (int i = 0; i < result.size(); ++i)
                {
                    result[i] = 0;
                }
            }
            face->setFaceMatrix(result);
        }
    }

    void doAllElementIntegration()
    {
        Integration::ElementIntegral integral(false);
        LinearAlgebra::Matrix result;
        for (Base::Element* element : meshes_[0]->getElementsList())
        {
            result = integral.integrate(element, &stifnessIntegrand);
            element->setElementMatrix(result);
        }
    }
    
    //messy routine that collects the row numbers of basisfunctions active on the free surface for use with PETSc
    void getSurfaceIS(Utilities::GlobalPetscMatrix& S, IS* surface, IS* rest)
    {
        Geometry::PointReference p(DIM);
        Geometry::PointPhysical pPhys(DIM);
        std::size_t numBasisFuns(0);
        std::vector<int> facePositions;
        for (const Base::Face* face : meshes_[0]->getFacesList())
        {
            p = face->getReferenceGeometry()->getCenter();
            pPhys = face->referenceToPhysical(p);
            if (std::abs(pPhys[DIM - 1]) < 1e-9)
            {
                S.getMatrixBCEntries(face, numBasisFuns, facePositions);
            }
        }
        ISCreateGeneral(PETSC_COMM_WORLD, numBasisFuns, &facePositions[0], PETSC_COPY_VALUES, surface);
        int totalFuncs;
        MatGetSize(S, &totalFuncs, PETSC_NULL);
        ISSort(*surface);
        ISComplement(*surface, 0, totalFuncs, rest);
        ISDestroy(surface);
        ISComplement(*rest, 0, totalFuncs, surface); //the complement of the complement does not contain duplicates
    }

    void printError()
    {
        Integration::FaceIntegral integral(false);
        LinearAlgebra::NumericalVector totalError(2), contribution(2);
        for (Base::Face* face : meshes_[0]->getFacesList())
        {
            contribution = integral.integrate(face, &error);
            totalError += contribution;
            contribution[0] = 0;
            contribution[1] = 0;
        }
        std::cout << "t: " << t << " eta: " << sqrt(totalError[0]) << " phi: " << sqrt(totalError[1]) << std::endl;
    }

    void computeEnergy()
    {
        Integration::ElementIntegral elIntegral(false);
        Integration::FaceIntegral faIntegral(false);
        LinearAlgebra::NumericalVector totalEnergy(1), contribution(1);
        for (Base::Face* face : meshes_[0]->getFacesList())
        {
            contribution = faIntegral.integrate(face, &faceEnergy);
            totalEnergy += contribution;
            contribution[0] = 0;
        }
        for (Base::Element* element : meshes_[0]->getElementsList())
        {
            contribution = elIntegral.integrate(element, &elementEnergy);
            totalEnergy += contribution;
            contribution[0] = 0;
        }
        std::cout << "Energy: " << totalEnergy[0] << std::endl;
    }

    bool solve()
    {
        std::cout.precision(14);
        doAllElementIntegration();
        doAllFaceIntegration();
        
        Utilities::GlobalPetscVector eta(meshes_[0]), phi(meshes_[0]);
        Utilities::GlobalPetscMatrix M(meshes_[0], -1, 0), S(meshes_[0], 0, -1);
        Vec phiS, phiOther, etaActually, interiorRHS, surfaceRHS;
        Mat surfaceMass, interiorStifness, surfaceStifness, mixStifness, backStiffness;
        
        IS isSurface, isRest;
        getSurfaceIS(S, &isSurface, &isRest);
        
        std::ofstream outFile("output.dat");
        Output::TecplotDiscontinuousSolutionWriter writeFunc(outFile, "test", "01", "eta_num,phi_num,eta_exact,phi_exact");
        
        //deal with initial conditions
        double g(1.), dt(M_PI / n_ / 2);
        eta.constructFromTimeLevelData(0, 0);
        phi.constructFromTimeLevelData(0, 1);
        
        int numberOfSnapshots(65); //placeholder parameters
        int numberOfTimeSteps(n_ / 8);
        
        MatGetSubMatrix(M, isSurface, isSurface, MAT_INITIAL_MATRIX, &surfaceMass);
        MatGetSubMatrix(S, isRest, isRest, MAT_INITIAL_MATRIX, &interiorStifness);
        MatGetSubMatrix(S, isSurface, isSurface, MAT_INITIAL_MATRIX, &surfaceStifness);
        MatGetSubMatrix(S, isRest, isSurface, MAT_INITIAL_MATRIX, &mixStifness);
        MatGetSubMatrix(S, isSurface, isRest, MAT_INITIAL_MATRIX, &backStiffness);
        
        KSP interior, surface;
        KSPCreate(PETSC_COMM_WORLD, &interior);
        KSPSetOperators(interior, interiorStifness, interiorStifness);
        KSPSetInitialGuessNonzero(interior, PETSC_TRUE);
        KSPSetFromOptions(interior);
        KSPCreate(PETSC_COMM_WORLD, &surface);
        KSPSetOperators(surface, surfaceMass, surfaceMass);
        KSPSetFromOptions(surface);
        KSPConvergedReason conferge;
        int iterations;
        
        VecGetSubVector(phi, isSurface, &phiS);
        VecGetSubVector(phi, isRest, &phiOther);
        VecGetSubVector(eta, isSurface, &etaActually);
        VecDuplicate(phiS, &surfaceRHS);
        VecDuplicate(phiOther, &interiorRHS);
        VecCopy(phiS, surfaceRHS);
        KSPSolve(surface, surfaceRHS, phiS);
        KSPGetConvergedReason(surface, &conferge);
        KSPGetIterationNumber(surface, &iterations);
        std::cout << "Finalizing interpolation (1): KSP solver ended because of " << KSPConvergedReasons[conferge] << " in " << iterations << " iterations." << std::endl;
        VecCopy(etaActually, surfaceRHS);
        KSPSolve(surface, surfaceRHS, etaActually);
        KSPGetConvergedReason(surface, &conferge);
        KSPGetIterationNumber(surface, &iterations);
        std::cout << "Finalizing interpolation (2): KSP solver ended because of " << KSPConvergedReasons[conferge] << " in " << iterations << " iterations." << std::endl;
        VecRestoreSubVector(eta, isSurface, &etaActually);
        VecRestoreSubVector(phi, isSurface, &phiS);
        VecRestoreSubVector(phi, isRest, &phiOther);
        
        eta.writeTimeLevelData(0, 0);
        phi.writeTimeLevelData(0, 1);
        
        writeFunc.write(meshes_[0], "solution at time 0", false, this);
        for (int i = 1; i < numberOfSnapshots; ++i)
        {
            VecGetSubVector(phi, isSurface, &phiS);
            VecGetSubVector(phi, isRest, &phiOther);
            VecGetSubVector(eta, isSurface, &etaActually);
            for (int j = 0; j < numberOfTimeSteps; ++j)
            {
                t += dt;
                VecAXPY(phiS, -g * dt / 2, etaActually); //(can in principle be combined with final step, but this makes output easier)
                MatMult(mixStifness, phiS, interiorRHS);
                VecScale(interiorRHS, -1);
                KSPSolve(interior, interiorRHS, phiOther);
                KSPGetConvergedReason(interior, &conferge);
                KSPGetIterationNumber(interior, &iterations);
                std::cout << "Laplace problem: KSP solver ended because of " << KSPConvergedReasons[conferge] << " in " << iterations << " iterations." << std::endl;
                
                MatMult(backStiffness, phiOther, surfaceRHS);
                MatMultAdd(surfaceStifness, phiS, surfaceRHS, surfaceRHS);
                VecScale(surfaceRHS, dt);
                MatMultAdd(surfaceMass, etaActually, surfaceRHS, surfaceRHS);
                KSPSolve(surface, surfaceRHS, etaActually);
                KSPGetConvergedReason(surface, &conferge);
                KSPGetIterationNumber(surface, &iterations);
                std::cout << "Updating \\eta: KSP solver ended because of " << KSPConvergedReasons[conferge] << " in " << iterations << " iterations." << std::endl;
                
                VecAXPY(phiS, -g * dt / 2, etaActually);
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
    
private:
    
    //number of elements per cardinal direction
    int n_;

    //polynomial order of the approximation
    int p_;
    static double t;
};

double DGWave::t = 0;

auto& n = Base::register_argument<std::size_t>('n', "numelems", "Number of Elements", true);
auto& p = Base::register_argument<std::size_t>('p', "poly", "Polynomial order", true);

int main(int argc, char **argv)
{
    Base::parse_options(argc, argv);
    try
    {
        DGWave test(n.getValue(), p.getValue());
        test.initialise();
        test.solve();
        return 0;
    }
    catch (const char* e)
    {
        std::cout << e;
    }
}

