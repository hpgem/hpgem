#ifndef HEULER_HH
#define HEULER_HH

#include <string>
using std::string;

//#include <petscmat.h>
#include <petscksp.h>


#include "Output/TecplotDiscontinuousSolutionWriter.hpp"
#include "Output/TecplotPhysicalGeometryIterator.hpp"
#include "LinearAlgebra/NumericalVector.hpp"
#include "Base/GlobalData.hpp"
#include "Base/ConfigurationData.hpp"
#include "Base/HpgemUI.hpp"
#include "Base/Norm2.hpp"
#include "Base/PhysGradientOfBasisFunction.hpp"

#include "InitialConditions.hpp"
#include "TecplotOutput.hpp"
#include "Integration/ElementIntegrandBase.hpp"
#include "Integration/FaceIntegrandBase.hpp"
using Base::RectangularMeshDescriptor;
using Base::HpgemUI;
using Base::GlobalData;
using Base::ConfigurationData;
using namespace Base;


const unsigned int  DIM = 3;

typedef std::vector<LinearAlgebra::Matrix>   VectorOfMatrices;
struct ElementIntegralData
{
        //optimize later!
    ElementIntegralData operator*= (const double& scalar){xGrad_*=scalar; yGrad_*=scalar; zGrad_*=scalar; return *this;}
    void axpy(double a, const ElementIntegralData& x){ xGrad_.axpy(a, x.xGrad_); yGrad_.axpy(a, x.yGrad_); zGrad_.axpy(a, x.zGrad_);}
    
    LinearAlgebra::Matrix xGrad_;
    LinearAlgebra::Matrix yGrad_;
    LinearAlgebra::Matrix zGrad_;
};

struct FluxData
{
    FluxData(unsigned int nb)
    {
        left_.resize(nb);
        right_.resize(nb);
        
        for (unsigned int i=0; i<nb; ++i)
        {
            LinearAlgebra::Matrix& left= left_[i];
            LinearAlgebra::Matrix& right= right_[i];
            
            left.resize(12,nb);
            right.resize(12,nb);
        }
    }
    void print()
    {
        cout <<"left="<<endl;
        for (unsigned int n=0; n < left_.size();++n)
        {
            cout <<"n="<<n<<", "<<left_[n]<<endl;
        }
        
        cout <<"right="<<endl;
        for (unsigned int n=0; n < right_.size();++n)
        {
            cout <<"n="<<n<<", "<<right_[n]<<endl;
        }
        
    }
    FluxData operator*= (const double& scalar)
    {
        for (unsigned int n=0; n < left_.size();++n)
        {
            left_[n]*=scalar;
            right_[n]*=scalar;
        }
        return *this;
    }
    void axpy(double a, const FluxData& x)
    {
        for (unsigned int n=0; n < left_.size();++n)
        {
            left_[n].axpy(a, x.left_[n]);
            right_[n].axpy(a, x.right_[n]);
        }
    }
    
    VectorOfMatrices    left_;
    VectorOfMatrices    right_;
};



struct HEulerElementData: public UserElementData
{
    HEulerElementData(unsigned int ndof):
        massMatrix_(ndof, ndof),
        invMassMatrix_(ndof, ndof)
    {
    }
    
    LinearAlgebra::Matrix   massMatrix_;
    LinearAlgebra::Matrix   invMassMatrix_;
};

struct HEulerGlobalVariables: public GlobalData
{
    unsigned int    nElements_;
    Mat             DivergenceFreeMatrix_;
    
    double          dt_;
};



struct HEulerConfigurationData: public ConfigurationData
{
    enum  SolutionType 		{INCOMPRESSIBLE_WALLS, INCOMPRESSIBLE_ONETHIRDPERIODIC, INCOMPRESSIBLE_PERIODIC, COMPRESSIBLE_WALLS, COMPRESSIBLE_PERIODIC};

    HEulerConfigurationData(unsigned int numberOfUnknowns, unsigned int numberOfBasisFunctions, unsigned int  numberOfTimeLevels=1, const string& fileName="in.txt", SolutionType type=COMPRESSIBLE_PERIODIC):
        ConfigurationData(3,numberOfUnknowns, 2, numberOfTimeLevels=1),
        solutionType_(type)
    {
            /// reading from a file
        theta_=0.5;
        numOfPeriods_=10;
        numOfTimeStepInOnePeriod_=100;
        numOfPeriodsInOnePlotStep_=10;
        onePeriod_=1;
        
    }
    
public:
    SolutionType    solutionType_;
    
    unsigned int    nx_;
    unsigned int    ny_;
    unsigned int    nz_;
    
    double          lx_;
    double          ly_;
    double          lz_;
    
    double          theta_;
    
    double          numOfPeriods_;
    double          numOfTimeStepInOnePeriod_;
    double          numOfPeriodsInOnePlotStep_;
    double          onePeriod_;
};

class HEuler: public HpgemUI,public Integration::ElementIntegrandBase<ElementIntegralData>,public Integration::FaceIntegrandBase<FluxData>,public Integration::ElementIntegrandBase<LinearAlgebra::Matrix>
{
public:
    typedef HpgemUI                        HpgemUIT;
    typedef Integration::ElementIntegral   ElementIntegralT;
    typedef Integration::FaceIntegral      FaceIntegralT;
    typedef ExactSolutionBase              ExactSolutionT;
    typedef PointReference               PointReferenceOnTheFaceT;
    
        //using HpgemUIT::ElementT;
        //    using HpgemUIT::PointReferenceT;
    
public:
    HEuler(HEulerGlobalVariables* global, const HEulerConfigurationData* config);

    ~HEuler();
public:
    
    
    void printFullMatrixInfo(Mat& matrix, const string& name);
    
    bool initialiseMesh();

    ///calculates mass matrix
    void elementIntegrand(const ElementT* element, const PointReferenceT& p, LinearAlgebra::Matrix& massMatrix);


    void calculateLocalEnergy(const ElementT& element, const PointReferenceT& p, double& returnValue);
    
    void elementIntegrand(const ElementT* element, const PointReferenceT& p, ElementIntegralData& returnObject);
    
    void faceIntegrand(const FaceT* face,          const NumericalVector& normal,
                       const PointReferenceOnTheFaceT& p,  FluxData& ret);
    
    void createCompressibleSystem();
    void createIncompressibleSystem();

    
        /// create Mass Matrices, store them as a User Defined Element Data
        /// calculate projection of the every unknown on the FEM spaces.
    void initialConditions();
    
    void solve();
    
    void output(double time=0.0);
    
private:///utilities
    void outputMatrix(Mat& matrix, const string& name);
    void outputMatrix(const Mat& matrix, const string& name)const;

    
    void outputVectorMatlab(Vec& vec, const string& name);
    void outputVectorMatlab(const Vec& vec, const string& name)const;
    
    
    void correctInitialProjectionOfVelocity(const Vec& UInit, Vec& UCorrected)const;
    
    void calculatePressure(const Mat& A, const Mat& Ah,const Vec& UCorrected);
    
    
private:
    ExactSolutionT*          exactSolution_;
    Mat                      P_;
    Mat                      Q_;
};
#endif
