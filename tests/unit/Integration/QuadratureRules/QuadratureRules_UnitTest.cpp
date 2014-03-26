#include <iostream>
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"
#include "Integration/QuadratureRules/QuadratureRuleSet.hpp"
#include "Integration/QuadratureRules/AllGaussQuadratureRules.hpp"
#include "Geometry/PointReference.hpp"
#include "Base/TestErrorDebug.hpp"

// forward declarations
void DescribeRule(const QuadratureRules::GaussQuadratureRule* qr);
void CollectAllRules(QuadratureRules::QuadratureRuleSet& QRset);

int main()
{
  QuadratureRules::QuadratureRuleSet QRset;
  CollectAllRules(QRset);

  int requestedOrder;

  std::cout << "==================================================\n";
  std::cout << " Query for quadrature rules for ReferenceLine\n";
  std::cout << "--------------------------------------------------\n";
  for (unsigned int requestedOrder = 1; requestedOrder < 10; ++requestedOrder)
  {
    std::cout << "\n**** request for order-" << requestedOrder << std::endl;
    DescribeRule(QRset.GetRule(&Geometry::ReferenceLine::Instance(), requestedOrder));
  }
  std::cout << "--------------------------------------------------\n\n";
  
// TestErrorDebug((1<0),"Gotcha!");

  std::cout << "==================================================\n";
  std::cout << " Query for quadrature rules for ReferenceSquare\n";
  std::cout << "--------------------------------------------------\n";
  for (unsigned int requestedOrder = 1; requestedOrder < 10; ++requestedOrder)
  {
    std::cout << "\n**** request for order-" << requestedOrder << std::endl;
    DescribeRule(QRset.GetRule(&Geometry::ReferenceSquare::Instance(), requestedOrder));
  }
  std::cout << "--------------------------------------------------\n\n";
  
  std::cout << "==================================================\n";
  std::cout << " Query for quadrature rules for ReferenceTriangle\n";
  std::cout << "--------------------------------------------------\n";
  for (unsigned int requestedOrder = 1; requestedOrder < 12; ++requestedOrder)
  {
    std::cout << "\n**** request for order-" << requestedOrder << std::endl;
    DescribeRule(QRset.GetRule(&Geometry::ReferenceTriangle::Instance(), requestedOrder));
  }
  std::cout << "--------------------------------------------------\n\n";
/*  
  std::cout << "==================================================\n";
  std::cout << " Query for quadrature rules for ReferenceCube\n";
  std::cout << "--------------------------------------------------\n";
  for (unsigned int requestedOrder = 1; requestedOrder < 10; ++requestedOrder)
  {
    std::cout << "\n**** request for order-" << requestedOrder << std::endl;
    DescribeRule(QRset.GetRule(&Geometry::ReferenceCube::Instance(), requestedOrder));
  }
  std::cout << "--------------------------------------------------\n\n";
  
  std::cout << "====================================================\n";
  std::cout << " Query for quadrature rules for ReferenceTetrahedron\n";
  std::cout << "----------------------------------------------------\n";
  for (unsigned int requestedOrder = 1; requestedOrder < 12; ++requestedOrder)
  {
    std::cout << "\n**** request for order-" << requestedOrder << std::endl;
    DescribeRule(QRset.GetRule(&Geometry::ReferenceTetrahedron::Instance(), requestedOrder));
  }
  std::cout << "--------------------------------------------------\n\n";
  
  std::cout << "=========================================================\n";
  std::cout << " Query for quadrature rules for ReferenceTriangularPrism\n";
  std::cout << "---------------------------------------------------------\n";
  for (unsigned int requestedOrder = 1; requestedOrder < 10; ++requestedOrder)
  {
    std::cout << "\n**** request for order-" << requestedOrder << std::endl;
    DescribeRule(QRset.GetRule(&Geometry::ReferenceTriangularPrism::Instance(), requestedOrder));
  }
  std::cout << "--------------------------------------------------\n\n";
  
  std::cout << "====================================================\n";
  std::cout << " Query for quadrature rules for ReferencePyramid\n";
  std::cout << "----------------------------------------------------\n";
  for (unsigned int requestedOrder = 1; requestedOrder < 10; ++requestedOrder)
  {
    std::cout << "\n**** request for order-" << requestedOrder << std::endl;
    DescribeRule(QRset.GetRule(&Geometry::ReferencePyramid::Instance(), requestedOrder));
  }
  std::cout << "--------------------------------------------------\n\n";
  
  std::cout << "====================================================\n";
  std::cout << " Query for quadrature rules for ReferenceHypercube\n";
  std::cout << "----------------------------------------------------\n";
  for (unsigned int requestedOrder = 1; requestedOrder < 10; ++requestedOrder)
  {
    std::cout << "\n**** request for order-" << requestedOrder << std::endl;
    DescribeRule(QRset.GetRule(&Geometry::ReferenceHypercube::Instance(), requestedOrder));
  }
  std::cout << "--------------------------------------------------\n\n";
*/  
  return 0;
}

//! Detaily describe a quadrature.
void DescribeRule(const QuadratureRules::GaussQuadratureRule* qr)
{
  if (qr!=NULL)
  {
    unsigned int nrPoints = qr->nrOfPoints();
    std::cout 	
        << "You get...."
	<< "name: " << qr->getName() 
	<< ", dimension: " << qr->dimension()
	<< ", order: " << qr->order()
	<< ", #points: " << nrPoints
	<< "\nPoints & weights:\n";
    
    for (unsigned int i=0; i<nrPoints; ++i)
    {
      Geometry::PointReference p(qr->dimension());
      qr->getPoint(i, p);
      std::cout << "\t" << p << " weight=" << qr->weight(i) << "\n";
    }
  }
  else 
     std::cout << "That's an empty rule!\n";
}


//! Collect all quadrature rules available .
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
/*    QRset.AddRule(&QuadratureRules::Cn3_1_1::Instance());
    QRset.AddRule(&QuadratureRules::Cn3_3_4::Instance());
    QRset.AddRule(&QuadratureRules::Cn3_5_9::Instance());
    QRset.AddRule(&QuadratureRules::C3_7_2::Instance());
    QRset.AddRule(&QuadratureRules::Tn3_1_1::Instance());
    QRset.AddRule(&QuadratureRules::Tn3_2_1::Instance());
    QRset.AddRule(&QuadratureRules::Tn3_3_1::Instance());
    QRset.AddRule(&QuadratureRules::Tn3_4_1::Instance());
    QRset.AddRule(&QuadratureRules::T3_5_1::Instance());
    QRset.AddRule(&QuadratureRules::T3_6_1::Instance());
    QRset.AddRule(&QuadratureRules::T3_7_1::Instance());
    QRset.AddRule(&QuadratureRules::T3_8_1::Instance());
    QRset.AddRule(&QuadratureRules::T3_9_1::Instance());
    QRset.AddRule(&QuadratureRules::T3_10_1::Instance());
    QRset.AddRule(&QuadratureRules::TriPrism_1_1::Instance());
    QRset.AddRule(&QuadratureRules::TriPrism_3_1::Instance());
    QRset.AddRule(&QuadratureRules::TriPrism_5_1::Instance());
    QRset.AddRule(&QuadratureRules::TriPrism_7_1::Instance());
    QRset.AddRule(&QuadratureRules::Pyramid_1_1::Instance());
    QRset.AddRule(&QuadratureRules::Pyramid_3_1::Instance());
    QRset.AddRule(&QuadratureRules::Pyramid_5_1::Instance());
    QRset.AddRule(&QuadratureRules::Pyramid_7_1::Instance());
    QRset.AddRule(&QuadratureRules::Cn4_1_1::Instance());
    QRset.AddRule(&QuadratureRules::Cn4_3_4::Instance());*/
}
    
