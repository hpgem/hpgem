#ifndef AssembleBasisFunctionSet_hpp
#define AssembleBasisFunctionSet_hpp

#include "BasisFunctionSet.hpp"

namespace Base
{

    //! Assemble set of BasisFunctions on 1D: { 1, x }
    void AssembleBasisFunctionSet_1D_Ord1_A0(Base::BasisFunctionSet& myBFSet);

    //! Assemble set of BasisFunctions on 1D: { 1, x, x^2 }
    void AssembleBasisFunctionSet_1D_Ord2_A0(Base::BasisFunctionSet& myBFSet);
    
    //! Assemble set of BasisFunctions on 1D: { 1, x, x^2, x^3 }
    void AssembleBasisFunctionSet_1D_Ord3_A0(Base::BasisFunctionSet& myBFSet);
    
    //! Assemble set of BasisFunctions on 1D: { 1, x, x^2, x^3,x^4 }
    void AssembleBasisFunctionSet_1D_Ord4_A0(Base::BasisFunctionSet& myBFSet);

    //! Assemble set of BasisFunctions on 1D: { 1, x, x^2, x^3,x^4,x^5 }
    void AssembleBasisFunctionSet_1D_Ord5_A0(Base::BasisFunctionSet& myBFSet);

    //! Assemble set of BasisFunctions on 2D: { 1, x, y }
    void AssembleBasisFunctionSet_2D_Ord1_A0(Base::BasisFunctionSet& myBFSet);
    
    //! Assemble set of BasisFunctions on 2D: { 1, x, y, xy }
    void AssembleBasisFunctionSet_2D_Ord1_A1(Base::BasisFunctionSet& myBFSet);
    
    //! Assemble set of BasisFunctions on 2D: { 1, x, y, x^2, y^2 }
    void AssembleBasisFunctionSet_2D_Ord2_A0(Base::BasisFunctionSet& myBFSet);
    
    //! Assemble set of BasisFunctions on 2D: { 1, x, y, x*y, x^2, y^2 }
    void AssembleBasisFunctionSet_2D_Ord2_A1(Base::BasisFunctionSet& myBFSet);
    
    //! Assemble set of BasisFunctions on 2D: { 1, x, y, x^2, y^2 and some others}
    void AssembleBasisFunctionSet_2D_Ord3_A0(Base::BasisFunctionSet& myBFSet);

    //! Assemble set of BasisFunctions on 2D: { 1, x, y, x*y, x^2, y^2 and some others}
    void AssembleBasisFunctionSet_2D_Ord3_A1(Base::BasisFunctionSet& myBFSet);

    //! Assemble set of BasisFunctions on 2D: { 1, x, y, x^2, y^2 and some others}
    void AssembleBasisFunctionSet_2D_Ord4_A0(Base::BasisFunctionSet& myBFSet);

    //! Assemble set of BasisFunctions on 2D: { 1, x, y, x*y, x^2, y^2 and some others}
    void AssembleBasisFunctionSet_2D_Ord4_A1(Base::BasisFunctionSet& myBFSet);

    //! Assemble set of BasisFunctions on 2D: { 1, x, y, x^2, y^2 and some others}
    void AssembleBasisFunctionSet_2D_Ord5_A0(Base::BasisFunctionSet& myBFSet);

    //! Assemble set of BasisFunctions on 2D: { 1, x, y, x*y, x^2, y^2 and some others}
    void AssembleBasisFunctionSet_2D_Ord5_A1(Base::BasisFunctionSet& myBFSet);

    //! Assemble set of BasisFunctions on 3D: { 1, x, y, z }
    void AssembleBasisFunctionSet_3D_Ord1_A0(Base::BasisFunctionSet& myBFSet);
    
    //! Assemble set of BasisFunctions on 3D: { 1, x, y, z, xy, xz, yz }
    void AssembleBasisFunctionSet_3D_Ord1_A1(Base::BasisFunctionSet& myBFSet);

    //! Assemble set of BasisFunctions on 3D: { 1, x, y, z, x^2, y^2, z^2 }
    void AssembleBasisFunctionSet_3D_Ord2_A0(Base::BasisFunctionSet& myBFSet);

    //! Assemble set of BasisFunctions on 3D: { 1, x, y, z, xy, xz, yz, x^2, y^2, z^2, xyz }
    void AssembleBasisFunctionSet_3D_Ord2_A1(Base::BasisFunctionSet& myBFSet);

    //! Assemble set of BasisFunctions on 3D: { 1, x, y, z, x^2, y^2, z^2 and some others}
    void AssembleBasisFunctionSet_3D_Ord3_A0(Base::BasisFunctionSet& myBFSet);

    //! Assemble set of BasisFunctions on 3D: { 1, x, y, z, xy, xz, yz, x^2, y^2, z^2, and some others}
    void AssembleBasisFunctionSet_3D_Ord3_A1(Base::BasisFunctionSet& myBFSet);

    //! Assemble set of BasisFunctions on 3D: { 1, x, y, z, x^2, y^2, z^2 and some others}
    void AssembleBasisFunctionSet_3D_Ord4_A0(Base::BasisFunctionSet& myBFSet);

    //! Assemble set of BasisFunctions on 3D: { 1, x, y, z, xy, xz, yz, x^2, y^2, z^2, and some others}
    void AssembleBasisFunctionSet_3D_Ord4_A1(Base::BasisFunctionSet& myBFSet);

    //! Assemble set of BasisFunctions on 3D: { 1, x, y, z, x^2, y^2, z^2 and some others}
    void AssembleBasisFunctionSet_3D_Ord5_A0(Base::BasisFunctionSet& myBFSet);

    //! Assemble set of BasisFunctions on 3D: { 1, x, y, z, xy, xz, yz, x^2, y^2, z^2, and some others}
    void AssembleBasisFunctionSet_3D_Ord5_A1(Base::BasisFunctionSet& myBFSet);

};

#endif
