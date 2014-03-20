#include <vector>

// #include "BaseBasisFunction.hpp"
#include "BasisFunctionsCollection_A.hpp"
#include "AssembleBasisFunctionSet.hpp"
#include "BasisFunctionSet.hpp"

namespace Base
{
    //! Assemble set of BasisFunctions on 1D: { 1, x }
    void AssembleBasisFunctionSet_1D_Ord1_A0(Base::BasisFunctionSet& myBFSet)
    {
      myBFSet.addBasisFunction(new Base::Basis_A0_1D);
      myBFSet.addBasisFunction(new Base::Basis_A1_1D);
    }
    
    //! Assemble set of BasisFunctions on 1D: { 1, x, x^2 }
    void AssembleBasisFunctionSet_1D_Ord2_A0(Base::BasisFunctionSet& myBFSet)
    {
      myBFSet.addBasisFunction(new Base::Basis_A0_1D);
      myBFSet.addBasisFunction(new Base::Basis_A1_1D);
      myBFSet.addBasisFunction(new Base::Basis_A2_1D);
    }
    
    //! Assemble set of BasisFunctions on 1D: { 1, x, x^2, x^3 }
    void AssembleBasisFunctionSet_1D_Ord3_A0(Base::BasisFunctionSet& myBFSet)
    {
      myBFSet.addBasisFunction(new Base::Basis_A0_1D);
      myBFSet.addBasisFunction(new Base::Basis_A1_1D);
      myBFSet.addBasisFunction(new Base::Basis_A2_1D);
      myBFSet.addBasisFunction(new Base::Basis_A3_1D);
    }
    
    //! Assemble set of BasisFunctions on 1D: { 1, x, x^2, x^3 }
    void AssembleBasisFunctionSet_1D_Ord4_A0(Base::BasisFunctionSet& myBFSet)
    {
      myBFSet.addBasisFunction(new Base::Basis_A0_1D);
      myBFSet.addBasisFunction(new Base::Basis_A1_1D);
      myBFSet.addBasisFunction(new Base::Basis_A2_1D);
      myBFSet.addBasisFunction(new Base::Basis_A3_1D);
      myBFSet.addBasisFunction(new Base::Basis_A4_1D);
    }

    //! Assemble set of BasisFunctions on 1D: { 1, x, x^2, x^3 }
    void AssembleBasisFunctionSet_1D_Ord5_A0(Base::BasisFunctionSet& myBFSet)
    {
      myBFSet.addBasisFunction(new Base::Basis_A0_1D);
      myBFSet.addBasisFunction(new Base::Basis_A1_1D);
      myBFSet.addBasisFunction(new Base::Basis_A2_1D);
      myBFSet.addBasisFunction(new Base::Basis_A3_1D);
      myBFSet.addBasisFunction(new Base::Basis_A4_1D);
      myBFSet.addBasisFunction(new Base::Basis_A5_1D);
    }

    //! Assemble set of BasisFunctions on 2D: { 1, x, y }
    void AssembleBasisFunctionSet_2D_Ord1_A0(Base::BasisFunctionSet& myBFSet)
    {
      myBFSet.addBasisFunction(new Base::Basis_A0_2D);
      myBFSet.addBasisFunction(new Base::Basis_A1_2D);
      myBFSet.addBasisFunction(new Base::Basis_A2_2D);
    }
    
    //! Assemble set of BasisFunctions on 2D: { 1, x, y, xy }
    void AssembleBasisFunctionSet_2D_Ord1_A1(Base::BasisFunctionSet& myBFSet)
    {
      myBFSet.addBasisFunction(new Base::Basis_A0_2D);
      myBFSet.addBasisFunction(new Base::Basis_A1_2D);
      myBFSet.addBasisFunction(new Base::Basis_A2_2D);
      myBFSet.addBasisFunction(new Base::Basis_A3_2D);
    }
    
    //! Assemble set of BasisFunctions on 2D: { 1, x, y, x^2, y^2 }
    void AssembleBasisFunctionSet_2D_Ord2_A0(Base::BasisFunctionSet& myBFSet)
    {
      myBFSet.addBasisFunction(new Base::Basis_A0_2D);
      myBFSet.addBasisFunction(new Base::Basis_A1_2D);
      myBFSet.addBasisFunction(new Base::Basis_A2_2D);
      myBFSet.addBasisFunction(new Base::Basis_A4_2D);
      myBFSet.addBasisFunction(new Base::Basis_A5_2D);
    }
    
    //! Assemble set of BasisFunctions on 2D: { 1, x, y, x*y, x^2, y^2 }
    void AssembleBasisFunctionSet_2D_Ord2_A1(Base::BasisFunctionSet& myBFSet)
    {
      myBFSet.addBasisFunction(new Base::Basis_A0_2D);
      myBFSet.addBasisFunction(new Base::Basis_A1_2D);
      myBFSet.addBasisFunction(new Base::Basis_A2_2D);
      myBFSet.addBasisFunction(new Base::Basis_A3_2D);
      myBFSet.addBasisFunction(new Base::Basis_A4_2D);
      myBFSet.addBasisFunction(new Base::Basis_A5_2D);
    }
    
    //! Assemble set of BasisFunctions on 2D: { 1, x, y, x^2, y^2 }
    void AssembleBasisFunctionSet_2D_Ord3_A0(Base::BasisFunctionSet& myBFSet)
    {
      myBFSet.addBasisFunction(new Base::Basis_A0_2D);
      myBFSet.addBasisFunction(new Base::Basis_A1_2D);
      myBFSet.addBasisFunction(new Base::Basis_A2_2D);
      myBFSet.addBasisFunction(new Base::Basis_A4_2D);
      myBFSet.addBasisFunction(new Base::Basis_A5_2D);
      myBFSet.addBasisFunction(new Base::Basis_A8_2D);
      myBFSet.addBasisFunction(new Base::Basis_A9_2D);
    }

    //! Assemble set of BasisFunctions on 2D: { 1, x, y, x*y, x^2, y^2 }
    void AssembleBasisFunctionSet_2D_Ord3_A1(Base::BasisFunctionSet& myBFSet)
    {
      myBFSet.addBasisFunction(new Base::Basis_A0_2D);
      myBFSet.addBasisFunction(new Base::Basis_A1_2D);
      myBFSet.addBasisFunction(new Base::Basis_A2_2D);
      myBFSet.addBasisFunction(new Base::Basis_A3_2D);
      myBFSet.addBasisFunction(new Base::Basis_A4_2D);
      myBFSet.addBasisFunction(new Base::Basis_A5_2D);
      myBFSet.addBasisFunction(new Base::Basis_A6_2D);
      myBFSet.addBasisFunction(new Base::Basis_A7_2D);
      myBFSet.addBasisFunction(new Base::Basis_A8_2D);
      myBFSet.addBasisFunction(new Base::Basis_A9_2D);
    }

    //! Assemble set of BasisFunctions on 2D: { 1, x, y, x^2, y^2 }
    void AssembleBasisFunctionSet_2D_Ord4_A0(Base::BasisFunctionSet& myBFSet)
    {
      myBFSet.addBasisFunction(new Base::Basis_A0_2D);
      myBFSet.addBasisFunction(new Base::Basis_A1_2D);
      myBFSet.addBasisFunction(new Base::Basis_A2_2D);
      myBFSet.addBasisFunction(new Base::Basis_A4_2D);
      myBFSet.addBasisFunction(new Base::Basis_A5_2D);
      myBFSet.addBasisFunction(new Base::Basis_A8_2D);
      myBFSet.addBasisFunction(new Base::Basis_A9_2D);
      myBFSet.addBasisFunction(new Base::Basis_A13_2D);
      myBFSet.addBasisFunction(new Base::Basis_A14_2D);
    }

    //! Assemble set of BasisFunctions on 2D: { 1, x, y, x*y, x^2, y^2 }
    void AssembleBasisFunctionSet_2D_Ord4_A1(Base::BasisFunctionSet& myBFSet)
    {
      myBFSet.addBasisFunction(new Base::Basis_A0_2D);
      myBFSet.addBasisFunction(new Base::Basis_A1_2D);
      myBFSet.addBasisFunction(new Base::Basis_A2_2D);
      myBFSet.addBasisFunction(new Base::Basis_A3_2D);
      myBFSet.addBasisFunction(new Base::Basis_A4_2D);
      myBFSet.addBasisFunction(new Base::Basis_A5_2D);
      myBFSet.addBasisFunction(new Base::Basis_A6_2D);
      myBFSet.addBasisFunction(new Base::Basis_A7_2D);
      myBFSet.addBasisFunction(new Base::Basis_A8_2D);
      myBFSet.addBasisFunction(new Base::Basis_A9_2D);
      myBFSet.addBasisFunction(new Base::Basis_A10_2D);
      myBFSet.addBasisFunction(new Base::Basis_A11_2D);
      myBFSet.addBasisFunction(new Base::Basis_A12_2D);
      myBFSet.addBasisFunction(new Base::Basis_A13_2D);
      myBFSet.addBasisFunction(new Base::Basis_A14_2D);
    }

    //! Assemble set of BasisFunctions on 2D: { 1, x, y, x^2, y^2 }
    void AssembleBasisFunctionSet_2D_Ord5_A0(Base::BasisFunctionSet& myBFSet)
    {
      myBFSet.addBasisFunction(new Base::Basis_A0_2D);
      myBFSet.addBasisFunction(new Base::Basis_A1_2D);
      myBFSet.addBasisFunction(new Base::Basis_A2_2D);
      myBFSet.addBasisFunction(new Base::Basis_A4_2D);
      myBFSet.addBasisFunction(new Base::Basis_A5_2D);
      myBFSet.addBasisFunction(new Base::Basis_A8_2D);
      myBFSet.addBasisFunction(new Base::Basis_A9_2D);
      myBFSet.addBasisFunction(new Base::Basis_A13_2D);
      myBFSet.addBasisFunction(new Base::Basis_A14_2D);
      myBFSet.addBasisFunction(new Base::Basis_A19_2D);
      myBFSet.addBasisFunction(new Base::Basis_A20_2D);
    }

    //! Assemble set of BasisFunctions on 2D: { 1, x, y, x*y, x^2, y^2 }
    void AssembleBasisFunctionSet_2D_Ord5_A1(Base::BasisFunctionSet& myBFSet)
    {
      myBFSet.addBasisFunction(new Base::Basis_A0_2D);
      myBFSet.addBasisFunction(new Base::Basis_A1_2D);
      myBFSet.addBasisFunction(new Base::Basis_A2_2D);
      myBFSet.addBasisFunction(new Base::Basis_A3_2D);
      myBFSet.addBasisFunction(new Base::Basis_A4_2D);
      myBFSet.addBasisFunction(new Base::Basis_A5_2D);
      myBFSet.addBasisFunction(new Base::Basis_A6_2D);
      myBFSet.addBasisFunction(new Base::Basis_A7_2D);
      myBFSet.addBasisFunction(new Base::Basis_A8_2D);
      myBFSet.addBasisFunction(new Base::Basis_A9_2D);
      myBFSet.addBasisFunction(new Base::Basis_A10_2D);
      myBFSet.addBasisFunction(new Base::Basis_A11_2D);
      myBFSet.addBasisFunction(new Base::Basis_A12_2D);
      myBFSet.addBasisFunction(new Base::Basis_A13_2D);
      myBFSet.addBasisFunction(new Base::Basis_A14_2D);
      myBFSet.addBasisFunction(new Base::Basis_A15_2D);
      myBFSet.addBasisFunction(new Base::Basis_A16_2D);
      myBFSet.addBasisFunction(new Base::Basis_A17_2D);
      myBFSet.addBasisFunction(new Base::Basis_A18_2D);
      myBFSet.addBasisFunction(new Base::Basis_A19_2D);
      myBFSet.addBasisFunction(new Base::Basis_A20_2D);
    }

    //! Assemble set of BasisFunctions on 3D: { 1, x, y, z}
    void AssembleBasisFunctionSet_3D_Ord0_A0(Base::BasisFunctionSet& myBFSet)
    {
      myBFSet.addBasisFunction(new Base::Basis_A0_3D);
      myBFSet.addBasisFunction(new Base::Basis_A1_3D);
      myBFSet.addBasisFunction(new Base::Basis_A2_3D);
      myBFSet.addBasisFunction(new Base::Basis_A3_3D);
    }

    //! Assemble set of BasisFunctions on 3D: { 1, x, y, z, xy, xz, yz }
    void AssembleBasisFunctionSet_3D_Ord1_A1(Base::BasisFunctionSet& myBFSet)
    {
      myBFSet.addBasisFunction(new Base::Basis_A0_3D);
      myBFSet.addBasisFunction(new Base::Basis_A1_3D);
      myBFSet.addBasisFunction(new Base::Basis_A2_3D);
      myBFSet.addBasisFunction(new Base::Basis_A3_3D);
      myBFSet.addBasisFunction(new Base::Basis_A4_3D);
      myBFSet.addBasisFunction(new Base::Basis_A5_3D);
      myBFSet.addBasisFunction(new Base::Basis_A6_3D);
    }

    //! Assemble set of BasisFunctions on 3D: { 1, x, y, z, x^2, y^2, z^2 }
    void AssembleBasisFunctionSet_3D_Ord2_A0(Base::BasisFunctionSet& myBFSet)
    {
      myBFSet.addBasisFunction(new Base::Basis_A0_3D);
      myBFSet.addBasisFunction(new Base::Basis_A1_3D);
      myBFSet.addBasisFunction(new Base::Basis_A2_3D);
      myBFSet.addBasisFunction(new Base::Basis_A3_3D);
      myBFSet.addBasisFunction(new Base::Basis_A7_3D);
      myBFSet.addBasisFunction(new Base::Basis_A8_3D);
      myBFSet.addBasisFunction(new Base::Basis_A9_3D);
    }
    
    //! Assemble set of BasisFunctions on 3D: { 1, x, y, z, xy, xz, yz, x^2, y^2, z^2, xyz }
    void AssembleBasisFunctionSet_3D_Ord2_A1(Base::BasisFunctionSet& myBFSet)
    {
      myBFSet.addBasisFunction(new Base::Basis_A0_3D);
      myBFSet.addBasisFunction(new Base::Basis_A1_3D);
      myBFSet.addBasisFunction(new Base::Basis_A2_3D);
      myBFSet.addBasisFunction(new Base::Basis_A3_3D);
      myBFSet.addBasisFunction(new Base::Basis_A4_3D);
      myBFSet.addBasisFunction(new Base::Basis_A5_3D);
      myBFSet.addBasisFunction(new Base::Basis_A6_3D);
      myBFSet.addBasisFunction(new Base::Basis_A7_3D);
      myBFSet.addBasisFunction(new Base::Basis_A8_3D);
      myBFSet.addBasisFunction(new Base::Basis_A9_3D);
      myBFSet.addBasisFunction(new Base::Basis_A10_3D);
    }

    //! Assemble set of BasisFunctions on 3D: { 1, x, y, z, x^2, y^2, z^2 }
    void AssembleBasisFunctionSet_3D_Ord3_A0(Base::BasisFunctionSet& myBFSet)
    {
      myBFSet.addBasisFunction(new Base::Basis_A0_3D);
      myBFSet.addBasisFunction(new Base::Basis_A1_3D);
      myBFSet.addBasisFunction(new Base::Basis_A2_3D);
      myBFSet.addBasisFunction(new Base::Basis_A3_3D);
      myBFSet.addBasisFunction(new Base::Basis_A7_3D);
      myBFSet.addBasisFunction(new Base::Basis_A8_3D);
      myBFSet.addBasisFunction(new Base::Basis_A9_3D);
      myBFSet.addBasisFunction(new Base::Basis_A17_3D);
      myBFSet.addBasisFunction(new Base::Basis_A18_3D);
      myBFSet.addBasisFunction(new Base::Basis_A19_3D);
    }

    //! Assemble set of BasisFunctions on 3D: { 1, x, y, z, xy, xz, yz, x^2, y^2, z^2, xyz }
    void AssembleBasisFunctionSet_3D_Ord3_A1(Base::BasisFunctionSet& myBFSet)
    {
      myBFSet.addBasisFunction(new Base::Basis_A0_3D);
      myBFSet.addBasisFunction(new Base::Basis_A1_3D);
      myBFSet.addBasisFunction(new Base::Basis_A2_3D);
      myBFSet.addBasisFunction(new Base::Basis_A3_3D);
      myBFSet.addBasisFunction(new Base::Basis_A4_3D);
      myBFSet.addBasisFunction(new Base::Basis_A5_3D);
      myBFSet.addBasisFunction(new Base::Basis_A6_3D);
      myBFSet.addBasisFunction(new Base::Basis_A7_3D);
      myBFSet.addBasisFunction(new Base::Basis_A8_3D);
      myBFSet.addBasisFunction(new Base::Basis_A9_3D);
      myBFSet.addBasisFunction(new Base::Basis_A10_3D);
      myBFSet.addBasisFunction(new Base::Basis_A11_3D);
      myBFSet.addBasisFunction(new Base::Basis_A12_3D);
      myBFSet.addBasisFunction(new Base::Basis_A13_3D);
      myBFSet.addBasisFunction(new Base::Basis_A14_3D);
      myBFSet.addBasisFunction(new Base::Basis_A15_3D);
      myBFSet.addBasisFunction(new Base::Basis_A16_3D);
      myBFSet.addBasisFunction(new Base::Basis_A17_3D);
      myBFSet.addBasisFunction(new Base::Basis_A18_3D);
      myBFSet.addBasisFunction(new Base::Basis_A19_3D);
    }

    //! Assemble set of BasisFunctions on 3D: { 1, x, y, z, x^2, y^2, z^2 }
    void AssembleBasisFunctionSet_3D_Ord4_A0(Base::BasisFunctionSet& myBFSet)
    {
      myBFSet.addBasisFunction(new Base::Basis_A0_3D);
      myBFSet.addBasisFunction(new Base::Basis_A1_3D);
      myBFSet.addBasisFunction(new Base::Basis_A2_3D);
      myBFSet.addBasisFunction(new Base::Basis_A3_3D);
      myBFSet.addBasisFunction(new Base::Basis_A7_3D);
      myBFSet.addBasisFunction(new Base::Basis_A8_3D);
      myBFSet.addBasisFunction(new Base::Basis_A9_3D);
      myBFSet.addBasisFunction(new Base::Basis_A17_3D);
      myBFSet.addBasisFunction(new Base::Basis_A18_3D);
      myBFSet.addBasisFunction(new Base::Basis_A19_3D);
      myBFSet.addBasisFunction(new Base::Basis_A32_3D);
      myBFSet.addBasisFunction(new Base::Basis_A33_3D);
      myBFSet.addBasisFunction(new Base::Basis_A34_3D);
    }

    //! Assemble set of BasisFunctions on 3D: { 1, x, y, z, xy, xz, yz, x^2, y^2, z^2, xyz }
    void AssembleBasisFunctionSet_3D_Ord4_A1(Base::BasisFunctionSet& myBFSet)
    {
      myBFSet.addBasisFunction(new Base::Basis_A0_3D);
      myBFSet.addBasisFunction(new Base::Basis_A1_3D);
      myBFSet.addBasisFunction(new Base::Basis_A2_3D);
      myBFSet.addBasisFunction(new Base::Basis_A3_3D);
      myBFSet.addBasisFunction(new Base::Basis_A4_3D);
      myBFSet.addBasisFunction(new Base::Basis_A5_3D);
      myBFSet.addBasisFunction(new Base::Basis_A6_3D);
      myBFSet.addBasisFunction(new Base::Basis_A7_3D);
      myBFSet.addBasisFunction(new Base::Basis_A8_3D);
      myBFSet.addBasisFunction(new Base::Basis_A9_3D);
      myBFSet.addBasisFunction(new Base::Basis_A10_3D);
      myBFSet.addBasisFunction(new Base::Basis_A11_3D);
      myBFSet.addBasisFunction(new Base::Basis_A12_3D);
      myBFSet.addBasisFunction(new Base::Basis_A13_3D);
      myBFSet.addBasisFunction(new Base::Basis_A14_3D);
      myBFSet.addBasisFunction(new Base::Basis_A15_3D);
      myBFSet.addBasisFunction(new Base::Basis_A16_3D);
      myBFSet.addBasisFunction(new Base::Basis_A17_3D);
      myBFSet.addBasisFunction(new Base::Basis_A18_3D);
      myBFSet.addBasisFunction(new Base::Basis_A19_3D);
      myBFSet.addBasisFunction(new Base::Basis_A20_3D);
      myBFSet.addBasisFunction(new Base::Basis_A21_3D);
      myBFSet.addBasisFunction(new Base::Basis_A22_3D);
      myBFSet.addBasisFunction(new Base::Basis_A23_3D);
      myBFSet.addBasisFunction(new Base::Basis_A24_3D);
      myBFSet.addBasisFunction(new Base::Basis_A25_3D);
      myBFSet.addBasisFunction(new Base::Basis_A26_3D);
      myBFSet.addBasisFunction(new Base::Basis_A27_3D);
      myBFSet.addBasisFunction(new Base::Basis_A28_3D);
      myBFSet.addBasisFunction(new Base::Basis_A29_3D);
      myBFSet.addBasisFunction(new Base::Basis_A30_3D);
      myBFSet.addBasisFunction(new Base::Basis_A31_3D);
      myBFSet.addBasisFunction(new Base::Basis_A32_3D);
      myBFSet.addBasisFunction(new Base::Basis_A33_3D);
      myBFSet.addBasisFunction(new Base::Basis_A34_3D);
    }

    //! Assemble set of BasisFunctions on 3D: { 1, x, y, z, x^2, y^2, z^2 }
    void AssembleBasisFunctionSet_3D_Ord5_A0(Base::BasisFunctionSet& myBFSet)
    {
      myBFSet.addBasisFunction(new Base::Basis_A0_3D);
      myBFSet.addBasisFunction(new Base::Basis_A1_3D);
      myBFSet.addBasisFunction(new Base::Basis_A2_3D);
      myBFSet.addBasisFunction(new Base::Basis_A3_3D);
      myBFSet.addBasisFunction(new Base::Basis_A7_3D);
      myBFSet.addBasisFunction(new Base::Basis_A8_3D);
      myBFSet.addBasisFunction(new Base::Basis_A9_3D);
      myBFSet.addBasisFunction(new Base::Basis_A17_3D);
      myBFSet.addBasisFunction(new Base::Basis_A18_3D);
      myBFSet.addBasisFunction(new Base::Basis_A19_3D);
      myBFSet.addBasisFunction(new Base::Basis_A32_3D);
      myBFSet.addBasisFunction(new Base::Basis_A33_3D);
      myBFSet.addBasisFunction(new Base::Basis_A34_3D);
      myBFSet.addBasisFunction(new Base::Basis_A53_3D);
      myBFSet.addBasisFunction(new Base::Basis_A54_3D);
      myBFSet.addBasisFunction(new Base::Basis_A55_3D);
    }

    //! Assemble set of BasisFunctions on 3D: { 1, x, y, z, xy, xz, yz, x^2, y^2, z^2, xyz }
    void AssembleBasisFunctionSet_3D_Ord5_A1(Base::BasisFunctionSet& myBFSet)
    {
      myBFSet.addBasisFunction(new Base::Basis_A0_3D);
      myBFSet.addBasisFunction(new Base::Basis_A1_3D);
      myBFSet.addBasisFunction(new Base::Basis_A2_3D);
      myBFSet.addBasisFunction(new Base::Basis_A3_3D);
      myBFSet.addBasisFunction(new Base::Basis_A4_3D);
      myBFSet.addBasisFunction(new Base::Basis_A5_3D);
      myBFSet.addBasisFunction(new Base::Basis_A6_3D);
      myBFSet.addBasisFunction(new Base::Basis_A7_3D);
      myBFSet.addBasisFunction(new Base::Basis_A8_3D);
      myBFSet.addBasisFunction(new Base::Basis_A9_3D);
      myBFSet.addBasisFunction(new Base::Basis_A10_3D);
      myBFSet.addBasisFunction(new Base::Basis_A11_3D);
      myBFSet.addBasisFunction(new Base::Basis_A12_3D);
      myBFSet.addBasisFunction(new Base::Basis_A13_3D);
      myBFSet.addBasisFunction(new Base::Basis_A14_3D);
      myBFSet.addBasisFunction(new Base::Basis_A15_3D);
      myBFSet.addBasisFunction(new Base::Basis_A16_3D);
      myBFSet.addBasisFunction(new Base::Basis_A17_3D);
      myBFSet.addBasisFunction(new Base::Basis_A18_3D);
      myBFSet.addBasisFunction(new Base::Basis_A19_3D);
      myBFSet.addBasisFunction(new Base::Basis_A20_3D);
      myBFSet.addBasisFunction(new Base::Basis_A21_3D);
      myBFSet.addBasisFunction(new Base::Basis_A22_3D);
      myBFSet.addBasisFunction(new Base::Basis_A23_3D);
      myBFSet.addBasisFunction(new Base::Basis_A24_3D);
      myBFSet.addBasisFunction(new Base::Basis_A25_3D);
      myBFSet.addBasisFunction(new Base::Basis_A26_3D);
      myBFSet.addBasisFunction(new Base::Basis_A27_3D);
      myBFSet.addBasisFunction(new Base::Basis_A28_3D);
      myBFSet.addBasisFunction(new Base::Basis_A29_3D);
      myBFSet.addBasisFunction(new Base::Basis_A30_3D);
      myBFSet.addBasisFunction(new Base::Basis_A31_3D);
      myBFSet.addBasisFunction(new Base::Basis_A32_3D);
      myBFSet.addBasisFunction(new Base::Basis_A33_3D);
      myBFSet.addBasisFunction(new Base::Basis_A34_3D);
      myBFSet.addBasisFunction(new Base::Basis_A35_3D);
      myBFSet.addBasisFunction(new Base::Basis_A36_3D);
      myBFSet.addBasisFunction(new Base::Basis_A37_3D);
      myBFSet.addBasisFunction(new Base::Basis_A38_3D);
      myBFSet.addBasisFunction(new Base::Basis_A39_3D);
      myBFSet.addBasisFunction(new Base::Basis_A40_3D);
      myBFSet.addBasisFunction(new Base::Basis_A41_3D);
      myBFSet.addBasisFunction(new Base::Basis_A42_3D);
      myBFSet.addBasisFunction(new Base::Basis_A43_3D);
      myBFSet.addBasisFunction(new Base::Basis_A44_3D);
      myBFSet.addBasisFunction(new Base::Basis_A45_3D);
      myBFSet.addBasisFunction(new Base::Basis_A46_3D);
      myBFSet.addBasisFunction(new Base::Basis_A47_3D);
      myBFSet.addBasisFunction(new Base::Basis_A48_3D);
      myBFSet.addBasisFunction(new Base::Basis_A49_3D);
      myBFSet.addBasisFunction(new Base::Basis_A50_3D);
      myBFSet.addBasisFunction(new Base::Basis_A51_3D);
      myBFSet.addBasisFunction(new Base::Basis_A52_3D);
      myBFSet.addBasisFunction(new Base::Basis_A53_3D);
      myBFSet.addBasisFunction(new Base::Basis_A54_3D);
      myBFSet.addBasisFunction(new Base::Basis_A55_3D);
    }

};
