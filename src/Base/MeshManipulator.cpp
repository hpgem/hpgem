//
//  MeshManipulator.cpp
//  
//
//  Created by Shavarsh Nurijanyan on 9/4/13.
//
//

#include "MeshManipulator.hpp"

namespace Base
{

        /// \brief This is the function that checks if two halfFaces nned to swapped or not.
        /// \details
    /*!
     *  First checks the firstNode index and if this is the same, swap based on the secondNode index.
     *  Therefore in the end that will be sorted based on firstNode and secondly on secondNode.
     !*/
    bool compareHalfFace(HalfFaceDescription first, HalfFaceDescription second){
        assert(first.nodeList.size()==second.nodeList.size());
        for(int i=0;i<first.nodeList.size();++i){
	    if(first.nodeList[i]>second.nodeList[i]){
	        return true;
	    }else if(second.nodeList[i]>first.nodeList[i]){
	        return false;
	    }
	}
	//give a consistent return value
	return first.elementNum>second.elementNum;
    }
//     bool compareHalfFace(HalfFaceDescription first, HalfFaceDescription second)
//     {
//         if (first.firstNode<second.firstNode) return true;
//         
//         if ((first.secondNode<second.secondNode) && (first.firstNode==second.firstNode)) return true;
//         
//         return false;
//     }
        //***********************************************************************************************************************
        //***********************************************************************************************************************
        //***********************************************************************************************************************

    template<>
    void
    MeshManipulator<1>::createDefaultBasisFunctions(unsigned int order)
    {
        Base::BasisFunctionSet<1>* bFset1 = new Base::BasisFunctionSet<1>(order);
        
        Base::AssembleBasisFunctionSet_1D_Ord3_A0(*bFset1);
        
        defaultSetOfBasisFunctions_ = bFset1;
    }

    template<>
    void
    MeshManipulator<2>::createDefaultBasisFunctions(unsigned int order)
    {
        Base::BasisFunctionSet<2>* bFset1 = new Base::BasisFunctionSet<2>(order);
        
        Base::AssembleBasisFunctionSet_2D_Ord2_A1(*bFset1);
        
        defaultSetOfBasisFunctions_ = bFset1;
    }

    template<>
    void
    MeshManipulator<3>::createDefaultBasisFunctions(unsigned int order)
    {
        Base::BasisFunctionSet<3>* bFset1 = new Base::BasisFunctionSet<3>(order);
        
        Base::AssembleBasisFunctionSet_3D_Ord2_A1(*bFset1);
        
        defaultSetOfBasisFunctions_ = bFset1;
    }
}