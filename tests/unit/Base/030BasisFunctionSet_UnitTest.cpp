/*
 * 030BasisFunctionSet_UnitTest.cpp
 *
 *  Created on: Apr 15, 2014
 *      Author: brinkf
 */

//naming convention: <Digit><ClassName>_UnitTest.cpp where <Digit> is a number that will make sure
//the unit tests are ordered such that the first failing unit test indicate the culprit class and
//other 'unit' tests may assume correct execution of all prior unit tests

#include "Base/AssembleBasisFunctionSet.hpp"
#include "cassert"

int main(){

	Base::BasisFunctionSet all1DbasisFunctions(5);
	Base::AssembleBasisFunctionSet_1D_Ord5_A0(all1DbasisFunctions);
	Geometry::PointReference point1D(1);
	for(int i=0;i<all1DbasisFunctions.size();++i){
		const Base::BaseBasisFunction* test=all1DbasisFunctions[i];
		for(point1D[0]=-1.5;point1D[0]<1.51;point1D[0]+=0.1){
			assert(("eval",test->eval(point1D)==all1DbasisFunctions.eval(i,point1D)));
			assert(("derivative",test->evalDeriv0(point1D)==all1DbasisFunctions.evalDeriv(i,0,point1D)));
		}
	}

	Base::BasisFunctionSet all2DbasisFunctions(5);
	Base::AssembleBasisFunctionSet_2D_Ord5_A0(all2DbasisFunctions);
	Geometry::PointReference point2D(2);
	for(int i=0;i<all2DbasisFunctions.size();++i){
		const Base::BaseBasisFunction* test=all2DbasisFunctions[i];
		for(point2D[0]=-1.5;point2D[0]<1.51;point2D[0]+=0.1){
			for(point2D[1]=-1.5;point2D[1]<1.51;point2D[1]+=0.1){
				assert(("eval",test->eval(point2D)==all2DbasisFunctions.eval(i,point2D)));
				assert(("derivative",test->evalDeriv0(point2D)==all2DbasisFunctions.evalDeriv(i,0,point2D)));
				assert(("derivative",test->evalDeriv1(point2D)==all2DbasisFunctions.evalDeriv(i,1,point2D)));
			}
		}
	}

	Base::BasisFunctionSet all3DbasisFunctions(5);
	Base::AssembleBasisFunctionSet_3D_Ord5_A0(all3DbasisFunctions);
	Geometry::PointReference point3D(3);
	for(int i=0;i<all3DbasisFunctions.size();++i){
		const Base::BaseBasisFunction* test=all3DbasisFunctions[i];
		for(point3D[0]=-1.5;point3D[0]<1.51;point3D[0]+=0.1){
			for(point3D[1]=-1.5;point3D[1]<1.51;point3D[1]+=0.1){
				for(point3D[2]=-1.5;point3D[2]<1.51;point3D[2]+=0.1){
					assert(("eval",test->eval(point3D)==all3DbasisFunctions.eval(i,point3D)));
					assert(("derivative",test->evalDeriv0(point3D)==all3DbasisFunctions.evalDeriv(i,0,point3D)));
					assert(("derivative",test->evalDeriv1(point3D)==all3DbasisFunctions.evalDeriv(i,1,point3D)));
					assert(("derivative",test->evalDeriv2(point3D)==all3DbasisFunctions.evalDeriv(i,2,point3D)));
				}
			}
		}
	}
}

