#include <vector>

// #include "BaseBasisFunction.hpp"
#include "BasisFunctionsCollection_A.hpp"
#include "AssembleBasisFunctionSet.hpp"
#include "BasisFunctionSet.hpp"

namespace Base
{
    //! Assemble set of BasisFunctions on 1D: { 1, x }
    void AssembleBasisFunctionSet_1D_Ord1_A0(Base::BasisFunctionSet<1>& myBFSet)
    {
      myBFSet.AddBasisFunction(new Base::Basis_A0_1D);
      myBFSet.AddBasisFunction(new Base::Basis_A1_1D);
    }
    
    //! Assemble set of BasisFunctions on 1D: { 1, x, x^2 }
    void AssembleBasisFunctionSet_1D_Ord2_A0(Base::BasisFunctionSet<1>& myBFSet)
    {
      myBFSet.AddBasisFunction(new Base::Basis_A0_1D);
      myBFSet.AddBasisFunction(new Base::Basis_A1_1D);
      myBFSet.AddBasisFunction(new Base::Basis_A2_1D);
    }
    
    //! Assemble set of BasisFunctions on 1D: { 1, x, x^2, x^3 }
    void AssembleBasisFunctionSet_1D_Ord3_A0(Base::BasisFunctionSet<1>& myBFSet)
    {
      myBFSet.AddBasisFunction(new Base::Basis_A0_1D);
      myBFSet.AddBasisFunction(new Base::Basis_A1_1D);
      myBFSet.AddBasisFunction(new Base::Basis_A2_1D);
      myBFSet.AddBasisFunction(new Base::Basis_A3_1D);
    }
    
    //! Assemble set of BasisFunctions on 2D: { 1, x, y }
    void AssembleBasisFunctionSet_2D_Ord1_A0(Base::BasisFunctionSet<2>& myBFSet)
    {
      myBFSet.AddBasisFunction(new Base::Basis_A0_2D);
      myBFSet.AddBasisFunction(new Base::Basis_A1_2D);
      myBFSet.AddBasisFunction(new Base::Basis_A2_2D);
    }
    
    //! Assemble set of BasisFunctions on 2D: { 1, x, y, xy }
    void AssembleBasisFunctionSet_2D_Ord1_A1(Base::BasisFunctionSet<2>& myBFSet)
    {
      myBFSet.AddBasisFunction(new Base::Basis_A0_2D);
      myBFSet.AddBasisFunction(new Base::Basis_A1_2D);
      myBFSet.AddBasisFunction(new Base::Basis_A2_2D);
      myBFSet.AddBasisFunction(new Base::Basis_A3_2D);
    }
    
    //! Assemble set of BasisFunctions on 2D: { 1, x, y, x^2, y^2 }
    void AssembleBasisFunctionSet_2D_Ord2_A0(Base::BasisFunctionSet<2>& myBFSet)
    {
      myBFSet.AddBasisFunction(new Base::Basis_A0_2D);
      myBFSet.AddBasisFunction(new Base::Basis_A1_2D);
      myBFSet.AddBasisFunction(new Base::Basis_A2_2D);
      myBFSet.AddBasisFunction(new Base::Basis_A4_2D);
      myBFSet.AddBasisFunction(new Base::Basis_A5_2D);
    }
    
    //! Assemble set of BasisFunctions on 2D: { 1, x, y, x*y, x^2, y^2 }
    void AssembleBasisFunctionSet_2D_Ord2_A1(Base::BasisFunctionSet<2>& myBFSet)
    {
      myBFSet.AddBasisFunction(new Base::Basis_A0_2D);
      myBFSet.AddBasisFunction(new Base::Basis_A1_2D);
      myBFSet.AddBasisFunction(new Base::Basis_A2_2D);
      myBFSet.AddBasisFunction(new Base::Basis_A3_2D);
      myBFSet.AddBasisFunction(new Base::Basis_A4_2D);
      myBFSet.AddBasisFunction(new Base::Basis_A5_2D);
    }
    
    //! Assemble set of BasisFunctions on 3D: { 1, x, y, z, x^2, y^2, z^2 }
    void AssembleBasisFunctionSet_3D_Ord2_A0(Base::BasisFunctionSet<3>& myBFSet)
    {
      myBFSet.AddBasisFunction(new Base::Basis_A0_3D);
      myBFSet.AddBasisFunction(new Base::Basis_A1_3D);
      myBFSet.AddBasisFunction(new Base::Basis_A2_3D);
      myBFSet.AddBasisFunction(new Base::Basis_A3_3D);
      myBFSet.AddBasisFunction(new Base::Basis_A7_3D);
      myBFSet.AddBasisFunction(new Base::Basis_A8_3D);
      myBFSet.AddBasisFunction(new Base::Basis_A9_3D);
    }
    
    //! Assemble set of BasisFunctions on 3D: { 1, x, y, z, xy, xz, yz, x^2, y^2, z^2, xyz }
    void AssembleBasisFunctionSet_3D_Ord2_A1(Base::BasisFunctionSet<3>& myBFSet)
    {
      myBFSet.AddBasisFunction(new Base::Basis_A0_3D);
      myBFSet.AddBasisFunction(new Base::Basis_A1_3D);
      myBFSet.AddBasisFunction(new Base::Basis_A2_3D);
      myBFSet.AddBasisFunction(new Base::Basis_A3_3D);
      myBFSet.AddBasisFunction(new Base::Basis_A4_3D);
      myBFSet.AddBasisFunction(new Base::Basis_A5_3D);
      myBFSet.AddBasisFunction(new Base::Basis_A6_3D);
      myBFSet.AddBasisFunction(new Base::Basis_A7_3D);
      myBFSet.AddBasisFunction(new Base::Basis_A8_3D);
      myBFSet.AddBasisFunction(new Base::Basis_A9_3D);
      myBFSet.AddBasisFunction(new Base::Basis_A10_3D);
    }

};
