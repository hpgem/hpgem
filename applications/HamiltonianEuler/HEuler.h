#ifndef HEULER_HH
#define HEULER_HH

#include <string>
using std::string;

#include <petscksp.h>

#include "Output/TecplotDiscontinuousSolutionWriter.h"
#include "Output/TecplotPhysicalGeometryIterator.h"
#include "LinearAlgebra/MiddleSizeVector.h"
#include "Base/GlobalData.h"
#include "Base/ConfigurationData.h"
#include "Base/HpgemAPIBase.h"
#include "Base/PhysGradientOfBasisFunction.h"
#include "Base/UserData.h"
#include "Integration/ElementIntegral.h"
#include "Integration/FaceIntegral.h"

#include "InitialConditions.h"
#include "TecplotOutput.h"
#include "Integration/ElementIntegrandBase.h"
#include "Integration/FaceIntegrandBase.h"
#include "Base/FaceCacheData.h"
#include "Base/ElementCacheData.h"
using Base::RectangularMeshDescriptor;
using Base::HpgemAPIBase;
using Base::GlobalData;
using Base::ConfigurationData;
using namespace Base;

const unsigned int DIM = 3;

using VectorOfMatrices = std::vector<LinearAlgebra::MiddleSizeMatrix>;
struct ElementIntegralData
{
    //optimize later!
    ElementIntegralData operator*=(const double& scalar)
    {
        xGrad_ *= scalar;
        yGrad_ *= scalar;
        zGrad_ *= scalar;
        return *this;
    }
    void axpy(double a, const ElementIntegralData& x)
    {
        xGrad_.axpy(a, x.xGrad_);
        yGrad_.axpy(a, x.yGrad_);
        zGrad_.axpy(a, x.zGrad_);
    }
    
    LinearAlgebra::MiddleSizeMatrix xGrad_;
    LinearAlgebra::MiddleSizeMatrix yGrad_;
    LinearAlgebra::MiddleSizeMatrix zGrad_;
};

struct FluxData
{
    void resize(std::size_t nb)
    {
        left_.resize(nb);
        right_.resize(nb);
        
        for (unsigned int i = 0; i < nb; ++i)
        {
            LinearAlgebra::MiddleSizeMatrix& left = left_[i];
            LinearAlgebra::MiddleSizeMatrix& right = right_[i];
            
            left.resize(12, nb);
            right.resize(12, nb);
        }
        
    }
    
    void print()
    {
        cout << "left=" << endl;
        for (unsigned int n = 0; n < left_.size(); ++n)
        {
            cout << "n=" << n << ", " << left_[n] << endl;
        }
        
        cout << "right=" << endl;
        for (unsigned int n = 0; n < right_.size(); ++n)
        {
            cout << "n=" << n << ", " << right_[n] << endl;
        }
        
    }
    FluxData operator*=(const double& scalar)
    {
        for (unsigned int n = 0; n < left_.size(); ++n)
        {
            left_[n] *= scalar;
            right_[n] *= scalar;
        }
        return *this;
    }
    void axpy(double a, const FluxData& x)
    {
        for (unsigned int n = 0; n < left_.size(); ++n)
        {
            left_[n].axpy(a, x.left_[n]);
            right_[n].axpy(a, x.right_[n]);
        }
    }
    
    VectorOfMatrices left_;
    VectorOfMatrices right_;
};

struct HEulerElementData : public UserElementData
{
    HEulerElementData(unsigned int ndof)
            : massMatrix_(ndof, ndof), invMassMatrix_(ndof, ndof)
    {
    }
    
    LinearAlgebra::MiddleSizeMatrix massMatrix_;
    LinearAlgebra::MiddleSizeMatrix invMassMatrix_;
};

struct HEulerGlobalVariables : public GlobalData
{
    unsigned int nElements_;
    Mat DivergenceFreeMatrix_;

    double dt_;
};

struct HEulerConfigurationData : public ConfigurationData
{
    enum SolutionType
    {
        INCOMPRESSIBLE_WALLS, INCOMPRESSIBLE_ONETHIRDPERIODIC, INCOMPRESSIBLE_PERIODIC, COMPRESSIBLE_WALLS, COMPRESSIBLE_PERIODIC
    };

    HEulerConfigurationData(unsigned int numberOfUnknowns, unsigned int numberOfBasisFunctions, unsigned int numberOfTimeLevels = 1, const string& fileName = "in.txt", SolutionType type = COMPRESSIBLE_PERIODIC)
            : ConfigurationData(3, numberOfUnknowns, 2, numberOfTimeLevels = 1), solutionType_(type)
    {
        /// reading from a file
        theta_ = 0.5;
        numOfPeriods_ = 10;
        numOfTimeStepInOnePeriod_ = 100;
        numOfPeriodsInOnePlotStep_ = 10;
        onePeriod_ = 1;
        
    }
    
public:
    SolutionType solutionType_;

    unsigned int nx_;
    unsigned int ny_;
    unsigned int nz_;

    double lx_;
    double ly_;
    double lz_;

    double theta_;

    double numOfPeriods_;
    double numOfTimeStepInOnePeriod_;
    double numOfPeriodsInOnePlotStep_;
    double onePeriod_;
};

class HEuler : public HpgemAPIBase<3>, public Integration::ElementIntegrandBase<ElementIntegralData, 3>, public Integration::FaceIntegrandBase<FluxData, 3>, public Integration::ElementIntegrandBase<LinearAlgebra::MiddleSizeMatrix, 3>
{
public:
    using ElementIntegralT = Integration::ElementIntegral<3>;
    using FaceIntegralT = Integration::FaceIntegral<3>;
    using ExactSolutionT = ExactSolutionBase;
    using PointReferenceOnTheFaceT = PointReference<2>;

    
public:
    HEuler(HEulerGlobalVariables* global, const HEulerConfigurationData* config);

    ~HEuler();
public:
    
    void printFullMatrixInfo(Mat& matrix, const string& name);

    bool initialiseMesh();

    ///calculates mass matrix
    void elementIntegrand(Base::PhysicalElement<3>& el, LinearAlgebra::MiddleSizeMatrix& massMatrix);

    void calculateLocalEnergy(Base::PhysicalElement<3>& el, double& returnValue);

    void elementIntegrand(Base::PhysicalElement<3>& el, ElementIntegralData& returnObject);

    void faceIntegrand(Base::PhysicalFace<3>& el, FluxData& ret);

    void createCompressibleSystem();
    void createIncompressibleSystem();

    /// create Mass Matrices, store them as a User Defined Element Data
    /// calculate projection of the every unknown on the FEM spaces.
    void initialConditions();

    void solve();

    void output(double time = 0.0);

private:
    ///utilities
    void outputMatrix(Mat& matrix, const string& name);
    void outputMatrix(const Mat& matrix, const string& name) const;

    void outputVectorMatlab(Vec& vec, const string& name);
    void outputVectorMatlab(const Vec& vec, const string& name) const;

    void correctInitialProjectionOfVelocity(const Vec& UInit, Vec& UCorrected) const;

    void calculatePressure(const Mat& A, const Mat& Ah, const Vec& UCorrected);

private:
    ExactSolutionT* exactSolution_;
    Mat P_;
    Mat Q_;
};
#endif
