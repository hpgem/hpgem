#include "cstring"
using std::string;


#include "Output/TecplotDiscontinuousSolutionWriter.hpp"
#include "Output/TecplotPhysicalGeometryIterator.hpp"
#include "LinearAlgebra/NumericalVector.hpp"
#include "Base/GlobalData.hpp"
#include "Base/ConfigurationData.hpp"
#include "Base/HpgemUI.hpp"
#include "Base/Norm2.hpp"
using Base::RectangularMeshDescriptor;
using Base::HpgemUI;
using Base::GlobalData;
using Base::ConfigurationData;



typedef std::vector<LinearAlgebra::Matrix>   VectorOfMatrices;
//struct FluxData
//{
//    VectorOfMatrices
//}
const unsigned int  DIM = 3;

struct HEulerGlobalVariables: public GlobalData
{
    
};



struct HEulerConfigurationData: public ConfigurationData
{
    enum  SolutionType 		{INCOMPRESSIBLE_WALLS, INCOMPRESSIBLE_ONETHIRDPERIODIC, INCOMPRESSIBLE_PERIODIC, COMPRESSIBLE_WALLS, COMPRESSIBLE_PERIODIC};
    
    HEulerConfigurationData(const string& fileName, SolutionType type):
        ConfigurationData(),
        solutionType_(type)
    {
            /// reading from a file
    }
    
    
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

class HEuler: public HpgemUI<DIM>
{
public:
    typedef HpgemUI<DIM>    HpgemUIT;
    
        //using HpgemUIT::ElementT;
        //    using HpgemUIT::PointReferenceT;
    
public:
    HEuler(HEulerGlobalVariables* global, const HEulerConfigurationData* config):
        HpgemUI(global, config)
    {
        
    }
    
public:
    
    bool initialiseMesh()
    {
        RectangularMeshDescriptor<DIM> rectangularMesh;
       
        rectangularMesh.bottomLeft_[0]       = 0;
        rectangularMesh.bottomLeft_[1]       = 0;
        rectangularMesh.bottomLeft_[2]       = 0;
        rectangularMesh.topLeft_[0]          = 1;
        rectangularMesh.topLeft_[1]          = 1;
        rectangularMesh.topLeft_[2]          = 1;
        rectangularMesh.numElementsInDIM_[0] = 8;
        rectangularMesh.numElementsInDIM_[1] = 8;
        rectangularMesh.numElementsInDIM_[2] = 8;
        
        addMesh(rectangularMesh);
        return true;
    }

    void calculateMassMatrix(const ElementT& element, const PointReferenceT& p, LinearAlgebra::Matrix& massMatrix)
    {
        unsigned int numOfDegreesOfFreedom = element.getNrOfBasisFunctions();
        
        massMatrix.resize(numOfDegreesOfFreedom, numOfDegreesOfFreedom);

        for (unsigned int i=0; i < numOfDegreesOfFreedom; ++i)
        {
            for (unsigned int j=0; j < numOfDegreesOfFreedom; ++j)
            {
                massMatrix(i,j) = element.basisFunction(j,p) * element.basisFunction(i,p);
            }
        }
    }
    
    void calculateLocalEnergy(const ElementT& element, const PointReferenceT& p, double& returnValue)
    {
        double extra=0;
        ElementT::SolutionVector  solution;
        
        element.getSolution(0, p, solution);
        
        HEulerConfigurationData::SolutionType sType = static_cast<const HEulerConfigurationData*>(configData_)->solutionType_;
        
        if (sType  == HEulerConfigurationData::COMPRESSIBLE_WALLS || sType  == HEulerConfigurationData::COMPRESSIBLE_PERIODIC)
        {
            extra =0.5*(solution[0]*solution[0]);
        }
        returnValue = 0.5*(solution[1]*solution[1] + solution[2]*solution[2] + solution[3]*solution[3])+extra;
    }
    
    void elementIntegrand(const ElementT& element, const PointReferenceT& p, VectorOfMatrices& returnObject)
    {
        returnObject.resize(DIM);
        
        unsigned int numberOfDegreesOfFreedom = element.getNrOfBasisFunctions();
        
        LinearAlgebra::Matrix& xDerReturnData = returnObject[0];
        LinearAlgebra::Matrix& yDerReturnData = returnObject[1];
        LinearAlgebra::Matrix& zDerReturnData = returnObject[2];

        xDerReturnData.resize(numberOfDegreesOfFreedom, numberOfDegreesOfFreedom);
        yDerReturnData.resize(numberOfDegreesOfFreedom, numberOfDegreesOfFreedom);
        zDerReturnData.resize(numberOfDegreesOfFreedom, numberOfDegreesOfFreedom);
        
        
        for (unsigned int i=0; i < numberOfDegreesOfFreedom; ++i)
        {
            for (unsigned int j=0; j < numberOfDegreesOfFreedom; ++j)
            {
                xDerReturnData(i,j) = element.basisFunction(j,p) * element.basisFunctionDeriv(i,0, p);
                yDerReturnData(i,j) = element.basisFunction(j,p) * element.basisFunctionDeriv(i,1, p);
                zDerReturnData(i,j) = element.basisFunction(j,p) * element.basisFunctionDeriv(i,2, p);
            }
        }
    }
    
    void faceIntegrand(const FaceT& face,          const PointPhysicalT& normal,
                       const PointReferenceT& p,  VectorOfMatrices& ret)
    {
        if (face.isInternal())
		{
            const double magn = Utilities::norm2<DIM>(normal);
            
            double theta = static_cast<const HEulerConfigurationData*>(configData_)->theta_;
            double nx    = normal[0]/magn;
            double ny    = normal[1]/magn;
            double nz    = normal[2]/magn;
            
            LinearAlgebra::Matrix& leftReturnData = ret[0];
            LinearAlgebra::Matrix& rightReturnData = ret[1];
            
//            BFevalL
//                /// provide numerical fluxes for 
//            leftReturnData(0,i) = 	BFevalL*nx*(1-theta)*bFLu;
//            leftReturnData(1,i) = 	BFevalL*ny*(1-theta)*bFLv;
//            leftReturnData(2,i) = 	BFevalL*nz*(1-theta)*bFL;
//            leftReturnData(3,i) = 	BFevalR*nx*(theta)*bFLu;
//            leftReturnData(4,i) = 	BFevalR*ny*(theta)*bFLv;
//            leftReturnData(5,i) = 	BFevalR*nz*(theta)*bFL;
//            
//            leftReturnData(6,i) = 	bFL*nx*(1-theta)*BFevalLu;
//            leftReturnData(7,i) = 	bFL*ny*(1-theta)*BFevalLv;
//            leftReturnData(8,i) = 	bFL*nz*(1-theta)*BFevalL;
//            leftReturnData(9,i) = 	bFL*nx*(theta)*BFevalRu;
//            leftReturnData(10,i) = 	bFL*ny*(theta)*BFevalRv;
//            leftReturnData(11,i) = 	bFL*nz*(theta)*BFevalR;
//            
//            
//            rightReturnData(0,i) = 	BFevalL*nx*(1-theta)*bFRu;
//            rightReturnData(1,i) = 	BFevalL*ny*(1-theta)*bFRv;
//            rightReturnData(2,i) = 	BFevalL*nz*(1-theta)*bFR;
//            rightReturnData(3,i) = 	BFevalR*nx*(theta)*bFRu;
//            rightReturnData(4,i) = 	BFevalR*ny*(theta)*bFRv;
//            rightReturnData(5,i) = 	BFevalR*nz*(theta)*bFR;
//            
//            rightReturnData(6,i) = 	bFR*nx*(1-theta)*BFevalLu;
//            rightReturnData(7,i) = 	bFR*ny*(1-theta)*BFevalLv;
//            rightReturnData(8,i) = 	bFR*nz*(1-theta)*BFevalL;
//            rightReturnData(9,i) = 	bFR*nx*(theta)*BFevalRu;
//            rightReturnData(10,i) = bFR*ny*(theta)*BFevalRv;
//            rightReturnData(11,i) = bFR*nz*(theta)*BFevalR;
        }
    }
    
    void initialConditions(const PointPhysicalT& p)
    {
            /// create Mass Matrices, store them as a User Defined Element Data
            /// calculate projection of the every unknown on the FEM spaces.
    }
    
    void output()
    {
        std::ofstream file3D;
        file3D.open ("out.dat");
        int dimensionsToWrite[2] = {0,1};
        Output::TecplotDiscontinuousSolutionWriter<DIM> out(file3D,"RectangularMesh",dimensionsToWrite,"xy");
        out.write(meshes_[0],"holi",false);
    }
};