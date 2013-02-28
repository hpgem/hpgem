#include "Base/MeshManipulator.hpp"
#include "Geometry/PointPhysical.hpp"
#include "Geometry/PointReference.hpp"
#include "Integration/ElementIntegral.hpp"
#include "Integration/FaceIntegral.hpp"
#include "Integration/ElementIntegrandBase.hpp"
#include "Integration/FaceIntegrandBase.hpp"
#include "Integration/QuadratureRules/AllGaussQuadratureRules.hpp"
#include "Integration/QuadratureRules/QuadratureRuleSet.hpp"

// forward declarations
void CollectAllRules(QuadratureRules::QuadratureRuleSet& QRset);

template <unsigned int dim, typename T>
class myElementIntegrandType : public Integration::ElementIntegrandBase<dim, T>
{
  public:
    ~myElementIntegrandType() {}

    virtual void operator()(const Base::Element<dim>& el, const Geometry::PointReference<dim>& p, T& ret)
    {
      std::cout << "Integrand: " << el.getReferenceGeometry()->getName() << " Gauss ";
      std::cout << p << std::endl;
    }
};


template <unsigned int dim, typename T>
class myFaceIntegrandType : public Integration::FaceIntegrandBase<dim, T>
{
  public:
    ~myFaceIntegrandType() {}

    virtual void operator()(const Base::Face<dim>& fa, 
                    const Geometry::PointPhysical<dim>& normal, 
                    const Geometry::PointReference<dim-1>& p, 
                    T& ret)
    {
      std::cout << fa.getReferenceGeometry()->getName() << " Gauss ";
      std::cout << p << std::endl;
    }
};



int main()
{
    const unsigned int dim=2;

    Geometry::PointPhysical<dim> bottomLeft, topLeft;
    std::vector<unsigned int> numElementsOneD(dim);

    bottomLeft[0]=0;  bottomLeft[1]=0;
    topLeft[0]=2;     topLeft[1]=2;
    numElementsOneD[0]=2;  numElementsOneD[1]=2;

    Base::MeshManipulator<dim> myTwoDDemoMesh(1,1);

    myTwoDDemoMesh.createRectangularMesh(bottomLeft,topLeft,numElementsOneD);

    typedef std::list<Base::Element<dim> > ListOfElementsT;
    typedef std::list<Base::Face<dim> >    ListOfFacesT;

    ListOfElementsT& elements = myTwoDDemoMesh.getElementsList();
    ListOfFacesT& faces = myTwoDDemoMesh.getFacesList();

    bool isUseCache(false);
    Integration::ElementIntegral<dim> elIntegral(isUseCache);
    Integration::FaceIntegral<dim>    faIntegral(isUseCache);
    
    myElementIntegrandType<dim, LinearAlgebra::NumericalVector>  myElIntegrand;
    myFaceIntegrandType<dim, LinearAlgebra::NumericalVector>     myFaIntegrand;
    LinearAlgebra::NumericalVector result;

    QuadratureRules::QuadratureRuleSet QRset;
    CollectAllRules(QRset);
    
    unsigned int reqQuadratureOrder = 2;
    for (ListOfElementsT::iterator el=elements.begin(); el!= elements.end(); ++el)
    {
        std::cout << el->getReferenceGeometry()->getName() << std::endl;
//         std::cout << el->getRefinementGeometry()->getName() << std::endl;

        // Set a quadrature rule to this element; this should be set somewhere in MeshManipulator(?)
//         QuadratureRules::GaussQuadratureRule<dim>* qr = QRset.GetRule(el->getReferenceGeometry(), reqQuadratureOrder);
        QuadratureRules::GaussQuadratureRule<dim>* qr = QRset.GetRule(&Geometry::ReferenceSquare::Instance(), reqQuadratureOrder);
//         QuadratureRules::GaussQuadratureRule<dim>* qr =  &QuadratureRules::Cn2_3_4::Instance();
        el->setGaussQuadratureRule(qr);

        // Integrate using a quadrature rule of the element
        elIntegral.Integrate<myElementIntegrandType<dim, LinearAlgebra::NumericalVector> >(*el, myElIntegrand, result);

        // Integrate using a given quadrature rule
        elIntegral.Integrate<myElementIntegrandType<dim, LinearAlgebra::NumericalVector> >(*el, qr, myElIntegrand, result);
    }
    
    for (ListOfFacesT::iterator fa=faces.begin(); fa!= faces.end(); ++fa)
    {
        std::cout << fa->getReferenceGeometry()->getName() << std::endl;
        
        // Set a quadrature rule to this face; this should be set somewhere in MeshManipulator(?)
//         QuadratureRules::GaussQuadratureRule<dim-1>* qr = QRset.GetRule(fa->getReferenceGeometry(), reqQuadratureOrder);
        QuadratureRules::GaussQuadratureRule<dim-1>* qr = QRset.GetRule(&Geometry::ReferenceLine::Instance(), reqQuadratureOrder);
//         QuadratureRules::GaussQuadratureRule<dim-1>* qr =  &QuadratureRules::Cn1_5_9::Instance();
        fa->setGaussQuadratureRule(qr);
        
        // Integrate using a quadrature rule of the element
        faIntegral.Integrate<myFaceIntegrandType<dim, LinearAlgebra::NumericalVector> >(*fa, myFaIntegrand, result);

        // Integrate using a given quadrature rule
        faIntegral.Integrate<myFaceIntegrandType<dim, LinearAlgebra::NumericalVector> >(*fa, qr, myFaIntegrand, result);
    }

    return 0;
}


void CollectAllRules(QuadratureRules::QuadratureRuleSet& QRset)
{
    QRset.AddRule(&QuadratureRules::Cn1_1_1::Instance());
    QRset.AddRule(&QuadratureRules::Cn1_3_4::Instance());
    QRset.AddRule(&QuadratureRules::Cn1_5_9::Instance());
    QRset.AddRule(&QuadratureRules::C1_7_x::Instance());
    QRset.AddRule(&QuadratureRules::Cn2_1_1::Instance());
    QRset.AddRule(&QuadratureRules::Cn2_3_4::Instance());
    QRset.AddRule(&QuadratureRules::Cn2_5_9::Instance());
    QRset.AddRule(&QuadratureRules::C2_7_4::Instance());
    QRset.AddRule(&QuadratureRules::Tn2_1_1::Instance());
    QRset.AddRule(&QuadratureRules::Tn2_2_1::Instance());
    QRset.AddRule(&QuadratureRules::Tn2_3_1::Instance());
    QRset.AddRule(&QuadratureRules::Tn2_4_1::Instance());
    QRset.AddRule(&QuadratureRules::T2_5_1::Instance());
    QRset.AddRule(&QuadratureRules::T2_6_1::Instance());
    QRset.AddRule(&QuadratureRules::T2_7_1::Instance());
    QRset.AddRule(&QuadratureRules::T2_8_1::Instance());
    QRset.AddRule(&QuadratureRules::T2_9_1::Instance());
    QRset.AddRule(&QuadratureRules::T2_10_1::Instance());
}

