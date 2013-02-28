#include <iostream>
#include <vector>

#include "Base/BasisFunctionsCollection_A.hpp"
#include "Base/AssembleBasisFunctionSet.hpp"
#include "Base/BasisFunctionSet.hpp"


int main()
{
  // order of a basis function set
  unsigned int order;
  
  // 1D basis functions------------------------------------------------
  Base::BaseBasisFunction<1>::PointReferenceT p1(1);
  p1[0] = .5;

  order = 1;
  Base::BasisFunctionSet<1> BFset1(order);
  Base::AssembleBasisFunctionSet_1D_Ord1_A0(BFset1);
  for (unsigned int i=0; i<BFset1.size(); ++i)
  {
      std::cout << BFset1.Eval(i,p1) << "\t";
  }
  std::cout << std::endl;
  
  
  // 2D basis functions------------------------------------------------
  Base::BaseBasisFunction<2>::PointReferenceT p2(2);
  p2[0] = .5;  p2[1] = .3;

  order = 1;
  Base::BasisFunctionSet<2> BFset2(order);
  Base::AssembleBasisFunctionSet_2D_Ord1_A1(BFset2);
  for (unsigned int i=0; i<BFset2.size(); ++i)
  {
      std::cout << BFset2.Eval(i,p2) << "\t";
  }
  std::cout << std::endl;
  
  order = 2;
  Base::BasisFunctionSet<2> BFset3(order);
  Base::AssembleBasisFunctionSet_2D_Ord2_A0(BFset3);
  for (unsigned int i=0; i<BFset3.size(); ++i)
  {
      std::cout << BFset3.Eval(i,p2) << "\t";
  }
  std::cout << std::endl;
  


  // 3D basis functions------------------------------------------------
  Base::BaseBasisFunction<3>::PointReferenceT p3(3);
  p3[0] = .5;  p3[1] = .3;  p3[1] = .75;

  order = 2;
  Base::BasisFunctionSet<3> BFset4(order);
  Base::AssembleBasisFunctionSet_3D_Ord2_A0(BFset4);
  for (unsigned int i=0; i<BFset4.size(); ++i)
  {
      std::cout << BFset4.Eval(i,p3) << "\t";
  }
  std::cout << std::endl;
  
  order = 2;
  Base::BasisFunctionSet<3> BFset5(order);
  Base::AssembleBasisFunctionSet_3D_Ord2_A1(BFset5);
  for (unsigned int i=0; i<BFset5.size(); ++i)
  {
      std::cout << BFset5.Eval(i,p3) << "\t";
  }
  std::cout << std::endl;
}
