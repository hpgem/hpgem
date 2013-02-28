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
class PhysicalFunction
{
    public:

    void operator()(const Geometry::PointPhysical<dim>& normal,
                    LinearAlgebra::NumericalVector& ret)
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

    Base::MeshManipulator<dim> myTwoDDemoMesh(2,2);

    myTwoDDemoMesh.createRectangularMesh(bottomLeft,topLeft,numElementsOneD);

    typedef std::list<Base::Element<dim>*> ListOfElementsT;
    typedef std::list<Base::Face<dim> >    ListOfFacesT;

    ListOfElementsT& elements = myTwoDDemoMesh.getElementsList();
    ListOfFacesT& faces = myTwoDDemoMesh.getFacesList();


    bool isUseCache(false);
    Integration::ElementIntegral<dim> elIntegral(isUseCache);
    Integration::FaceIntegral<dim>    faIntegral(isUseCache);
    
    LinearAlgebra::NumericalVector  result(1);

    PhysicalFunction<dim> physicalFunction;


    std::cout << "Integral over all elements....\n";
    unsigned int reqQuadratureOrder = 6;
    for (ListOfElementsT::iterator el=elements.begin(); el!= elements.end(); ++el)
    {
        std::cout << (*el)->getReferenceGeometry()->getName() << std::endl;
        // Integrate using a quadrature rule of the element
        (*el)->transformToReferenceElement<PhysicalFunction<dim> >(physicalFunction);
        elIntegral.integrate<PhysicalFunction<dim> >(*(*el), physicalFunction, result);
        cout << result;
        
        cout<< "#####################################END of ELEMENT######"<<endl;

    }
    std::cout << "Finished: Integral over all elements....\n";
    return 0;
}
