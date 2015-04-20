/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, University of Twente
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef AssembleBasisFunctionSet_hpp
#define AssembleBasisFunctionSet_hpp

namespace Base
{
	class BasisFunctionSet;

    //! Assemble set of BasisFunctions on 1D: { 1 }
    void AssembleBasisFunctionSet_1D_Ord0_A0(Base::BasisFunctionSet& myBFSet);

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
    
    

    //! Assemble set of BasisFunctions on 2D: { 1 }
    void AssembleBasisFunctionSet_2D_Ord0_A0(Base::BasisFunctionSet& myBFSet);

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

    
    
    //! Assemble set of BasisFunctions on 3D: { 1 }
    void AssembleBasisFunctionSet_3D_Ord0_A0(Base::BasisFunctionSet& myBFSet);

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
