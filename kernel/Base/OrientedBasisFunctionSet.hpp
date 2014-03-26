/*
 * OrientedBasisFunctionSet.hpp
 *
 *  Created on: Feb 21, 2014
 *      Author: brinkf
 */

#ifndef ORIENTEDBASISFUNCTIONSET_HPP_
#define ORIENTEDBASISFUNCTIONSET_HPP_

#include "BasisFunctionSet.hpp"

namespace Base {

/**
 * has all the functionality of a basisfunctionset, but is meant for cases
 * where the orientation of the basisfunctions matter
 *
 * it is assumed that this set contains only part of the basis (for example only the functions with their DOF linked to a face)
 */
class OrientedBasisFunctionSet: public Base::BasisFunctionSet {
public:
	OrientedBasisFunctionSet(int order, int orientation,int face):BasisFunctionSet(order),orientation_(orientation),face_(face){}
	virtual ~OrientedBasisFunctionSet(){};

	bool checkOrientation(int codim0mapIndex,int faceIndex)const{return codim0mapIndex==orientation_&&faceIndex==face_;}
private:
	int orientation_;
	int face_;
};

} /* namespace Base */

#endif /* ORIENTEDBASISFUNCTIONSET_HPP_ */
