#include "Base/MeshManipulator.hpp"
#include "Geometry/PointPhysical.hpp"
#include "Geometry/PointReference.hpp"
#include "Geometry/ReferenceGeometry.hpp"

#include "Integration/ElementIntegral.hpp"
#include "Integration/FaceIntegral.hpp"
#include "Integration/ElementIntegrandBase.hpp"
#include "Integration/FaceIntegrandBase.hpp"
#include "Integration/QuadratureRules/AllGaussQuadratureRules.hpp"

template <unsigned int dim>
class MyElementIntegrandType : public Integration::ElementIntegrandBase<dim>
{
public:

    void operator()(const Base::Element<dim>* elem, const Geometry::PointReference<dim>& p, LinearAlgebra::NumericalVector& ret)
    {
        ret[0] = p[0];
    }
};

template <unsigned int dim>
class MyFaceIntegrandType : public Integration::FaceIntegrandBase<dim>
{
public:

    void operator()(const Base::Face<dim>* face, const Geometry::PointPhysical<dim>& normal, const Geometry::PointReference<dim-1>& p, LinearAlgebra::NumericalVector& ret)
    {
        //ret[0] = p[0];
        //ret[1] = p[1];
    }
};

int main()
{
    const unsigned int dim=2;

    Geometry::PointPhysical<dim> bottomLeft, topLeft;
    std::vector<unsigned int> numElementsOneD(dim);

    bottomLeft[0] = 0;
    bottomLeft[1] = 0;

    topLeft[0] = 2;
    topLeft[1] = 2;

    numElementsOneD[0] = 2;
    numElementsOneD[1] = 2;

    Base::ConfigurationData config(1,1,1);
    
    config.numberOfUnknowns_       = 1;
    config.numberOfTimeLevels_     = 1;
    config.numberOfBasisFunctions_ = 1;
    
    Base::MeshManipulator<dim> myTwoDDemoMesh(&config, 1,1);

    myTwoDDemoMesh.createRectangularMesh(bottomLeft,topLeft,numElementsOneD);

    typedef std::list<Base::Element<dim>*> ListOfElementsT;
    typedef std::list<Base::Face<dim>*>    ListOfFacesT;

    ListOfElementsT& elements = myTwoDDemoMesh.getElementsList();
    ListOfFacesT& faces = myTwoDDemoMesh.getFacesList();


    //Create a 
    bool isUseCache(false);
    Integration::ElementIntegral<dim> 	elIntegral(isUseCache);
    Integration::FaceIntegral<dim>    	faIntegral(isUseCache);

    MyElementIntegrandType<dim>     	myElIntegrand;
    MyFaceIntegrandType<dim>        	myFaIntegrand;
    LinearAlgebra::NumericalVector  	result(1);

    std::cout << "Integral over all elements....\n";
    unsigned int reqQuadratureOrder = 6;

    for (ListOfElementsT::iterator el=elements.begin(); el!= elements.end(); ++el)
    {
        std::cout << (*el)->getReferenceGeometry()->getName() << std::endl;
        // Integrate using a quadrature rule of the element
        elIntegral.integrate((*el), myElIntegrand, result);
        cout << result;

        cout<< "#####################################END of ELEMENT######"<<endl;
	}

    std::cout << "Finished: Integral over all elements....\n";

    std::cout << "Integral over all faces....\n";

    for (ListOfFacesT::iterator fa=faces.begin(); fa!= faces.end(); ++fa)
    {
        std::cout << (*fa)->getReferenceGeometry()->getName() << std::endl;
        // Integrate using a quadrature rule of the element
        faIntegral.integrate(*fa, myFaIntegrand, result);
    }
    std::cout << "Finished: Integral over all faces....\n";

    return 0;
}
