/*
 * 010BasisFunctions_UnitTest.cpp
 *
 *  Created on: Apr 16, 2014
 *      Author: brinkf
 */

//naming convention: <Digit><ClassName>_UnitTest.cpp where <Digit> is a number that will make sure
//the unit tests are ordered such that the first failing unit test indicate the culprit class and
//other 'unit' tests may assume correct execution of all prior unit tests

//see also Base/???BasisFunction_UnitTest.cpp - running the same series of checks on different basisfunctions
#include "Utilities/BasisFunctions1DH1ConformingLine.hpp"
#include "Utilities/BasisFunctions2DH1ConformingSquare.hpp"
#include "Utilities/BasisFunctions2DH1ConformingTriangle.hpp"
#include "Utilities/BasisFunctions3DH1ConformingCube.hpp"
#include "Utilities/BasisFunctions3DH1ConformingTetrahedron.hpp"
#include "Utilities/BasisFunctions3DH1ConformingPrism.hpp"
#include "cassert"

#include "Base/BasisFunctionSet.hpp"
#include "Base/L2Norm.hpp"
using Base::L2Norm;

int main() {

	// 1D

	Base::BasisFunctionSet *all1DbasisFunctions = Utilities::createDGBasisFunctionSet1DH1Line(5);
	Geometry::PointReference point1D(1);
	LinearAlgebra::NumericalVector ret(1);
	for(int i=0;i<all1DbasisFunctions->size();++i){
		const Base::BaseBasisFunction* test = (*all1DbasisFunctions)[i];
		for(point1D[0]=-1.5;point1D[0]<1.51;point1D[0]+=0.1){
			test->eval(point1D,ret);
			assert(("eval",(test->eval(point1D)-ret[0])<1e-12));

			point1D[0]+=-1.e-8;
			double x0=test->eval(point1D);
			point1D[0]+=2.e-8;
			double x1=test->eval(point1D);

			point1D[0]+=-1e-8;
			test->evalDeriv(point1D,ret);
			assert(("derivative",fabs(ret[0]-5.e7*(x1-x0))<1e-5));
			assert(("derivative",fabs(test->evalDeriv0(point1D)-5.e7*(x1-x0))<1e-5));
		}
	}

	// 2D

	Base::BasisFunctionSet *all2DbasisFunctions = Utilities::createDGBasisFunctionSet2DH1Square(5);
	Geometry::PointReference point2D(2);
	for(int i=0;i<all2DbasisFunctions->size();++i){
		const Base::BaseBasisFunction* test = (*all2DbasisFunctions)[i];
		for(point2D[0]=-1.5;point2D[0]<1.51;point2D[0]+=0.1){
			for(point2D[1]=-1.5;point2D[1]<1.51;point2D[1]+=0.1){
				test->eval(point2D,ret);
				assert(("eval",(test->eval(point2D)-ret[0])<1e-12));

				point2D[0]+=-1.e-8;
				double x0=test->eval(point2D);
				point2D[0]+=2.e-8;
				double x1=test->eval(point2D);

				point2D[0]+=-1e-8;
				ret.resize(2);
				test->evalDeriv(point2D,ret);
				assert(("derivative",fabs(ret[0]-5.e7*(x1-x0))<1e-5));
				assert(("derivative",fabs(test->evalDeriv0(point2D)-5.e7*(x1-x0))<1e-5));

				point2D[1]+=-1.e-8;
				x0=test->eval(point2D);
				point2D[1]+=2.e-8;
				x1=test->eval(point2D);

				point2D[1]+=-1e-8;
				assert(("derivative",fabs(ret[1]-5.e7*(x1-x0))<1e-5));
				assert(("derivative",fabs(test->evalDeriv1(point2D)-5.e7*(x1-x0))<1e-5));

				ret.resize(1);
			}
		}
	}

	delete all2DbasisFunctions;

	all2DbasisFunctions = Utilities::createDGBasisFunctionSet2DH1Triangle(5);
	for(int i=0;i<all2DbasisFunctions->size();++i){
		const Base::BaseBasisFunction* test = (*all2DbasisFunctions)[i];
		for(point2D[0]=-1.5;point2D[0]<1.51;point2D[0]+=0.1){
			for(point2D[1]=-1.5;point2D[1]<1.51;point2D[1]+=0.1){
				test->eval(point2D,ret);
				assert(("eval",(test->eval(point2D)-ret[0])<1e-12));

				point2D[0]+=-1.e-8;
				double x0=test->eval(point2D);
				point2D[0]+=2.e-8;
				double x1=test->eval(point2D);

				point2D[0]+=-1e-8;
				ret.resize(2);
				test->evalDeriv(point2D,ret);//exact to within absolute OR relative tolerace of 1e-5
				assert(("derivative",fabs(ret[0]-5.e7*(x1-x0))<1e-5||(L2Norm(ret)>1&&fabs(ret[0]-5.e7*(x1-x0))<1e-5*L2Norm(ret))));
				assert(("derivative",fabs(test->evalDeriv0(point2D)-5.e7*(x1-x0))<1e-5||(L2Norm(ret)>1&&fabs(ret[0]-5.e7*(x1-x0))<1e-5*L2Norm(ret))));

				point2D[1]+=-1.e-8;
				x0=test->eval(point2D);
				point2D[1]+=2.e-8;
				x1=test->eval(point2D);

				point2D[1]+=-1e-8;
				assert(("derivative",fabs(ret[1]-5.e7*(x1-x0))<1e-5||(L2Norm(ret)>1&&fabs(ret[1]-5.e7*(x1-x0))<1e-5*L2Norm(ret))));
				assert(("derivative",fabs(test->evalDeriv1(point2D)-5.e7*(x1-x0))<1e-5||(L2Norm(ret)>1&&fabs(ret[1]-5.e7*(x1-x0))<1e-5*L2Norm(ret))));

				ret.resize(1);
			}
		}
	}

	delete all2DbasisFunctions;

	//3D

	Base::BasisFunctionSet *all3DbasisFunctions = Utilities::createDGBasisFunctionSet3DH1Cube(5);
	Geometry::PointReference point3D(3);
	for(int i=0;i<all3DbasisFunctions->size();++i){
		const Base::BaseBasisFunction* test = (*all3DbasisFunctions)[i];
		for(point3D[0]=-1.5;point3D[0]<1.51;point3D[0]+=0.25){
			for(point3D[1]=-1.5;point3D[1]<1.51;point3D[1]+=0.25){
				for(point3D[2]=-1.5;point3D[2]<1.51;point3D[2]+=0.25){
					test->eval(point3D,ret);
					assert(("eval",(test->eval(point3D)-ret[0])<1e-12));

					point3D[0]+=-1.e-8;
					double x0=test->eval(point3D);
					point3D[0]+=2.e-8;
					double x1=test->eval(point3D);

					point3D[0]+=-1e-8;
					ret.resize(3);
					test->evalDeriv(point3D,ret);
					assert(("derivative",fabs(ret[0]-5.e7*(x1-x0))<1e-5||(L2Norm(ret)>1&&fabs(ret[0]-5.e7*(x1-x0))<1e-5*L2Norm(ret))));
					assert(("derivative",fabs(test->evalDeriv0(point3D)-5.e7*(x1-x0))<1e-5||(L2Norm(ret)>1&&fabs(ret[0]-5.e7*(x1-x0))<1e-5*L2Norm(ret))));

					point3D[1]+=-1.e-8;
					x0=test->eval(point3D);
					point3D[1]+=2.e-8;
					x1=test->eval(point3D);

					point3D[1]+=-1e-8;
					assert(("derivative",fabs(ret[1]-5.e7*(x1-x0))<1e-5||(L2Norm(ret)>1&&fabs(ret[1]-5.e7*(x1-x0))<1e-5*L2Norm(ret))));
					assert(("derivative",fabs(test->evalDeriv1(point3D)-5.e7*(x1-x0))<1e-5||(L2Norm(ret)>1&&fabs(ret[1]-5.e7*(x1-x0))<1e-5*L2Norm(ret))));

					point3D[2]+=-1.e-8;
					x0=test->eval(point3D);
					point3D[2]+=2.e-8;
					x1=test->eval(point3D);

					point3D[2]+=-1e-8;
					assert(("derivative",fabs(ret[2]-5.e7*(x1-x0))<1e-5||(L2Norm(ret)>1&&fabs(ret[2]-5.e7*(x1-x0))<1e-5*L2Norm(ret))));
					assert(("derivative",fabs(test->evalDeriv2(point3D)-5.e7*(x1-x0))<1e-5||(L2Norm(ret)>1&&fabs(ret[2]-5.e7*(x1-x0))<1e-5*L2Norm(ret))));

					ret.resize(1);
				}
			}
		}
	}

	delete all3DbasisFunctions;

	all3DbasisFunctions = Utilities::createDGBasisFunctionSet3DH1Tetrahedron(5);
	for(int i=0;i<all3DbasisFunctions->size();++i){
		const Base::BaseBasisFunction* test = (*all3DbasisFunctions)[i];
		for(point3D[0]=-1.5;point3D[0]<1.51;point3D[0]+=0.25){
			for(point3D[1]=-1.5;point3D[1]<1.51;point3D[1]+=0.25){
				for(point3D[2]=-1.5;point3D[2]<1.51;point3D[2]+=0.25){
					test->eval(point3D,ret);
					assert(("eval",(test->eval(point3D)-ret[0])<1e-12));

					point3D[0]+=-1.e-8;
					double x0=test->eval(point3D);
					point3D[0]+=2.e-8;
					double x1=test->eval(point3D);

					point3D[0]+=-1e-8;
					ret.resize(3);
					test->evalDeriv(point3D,ret);
					assert(("derivative",fabs(ret[0]-5.e7*(x1-x0))<1e-5||(L2Norm(ret)>1&&fabs(ret[0]-5.e7*(x1-x0))<1e-5*L2Norm(ret))));
					assert(("derivative",fabs(test->evalDeriv0(point3D)-5.e7*(x1-x0))<1e-5||(L2Norm(ret)>1&&fabs(ret[0]-5.e7*(x1-x0))<1e-5*L2Norm(ret))));

					point3D[1]+=-1.e-8;
					x0=test->eval(point3D);
					point3D[1]+=2.e-8;
					x1=test->eval(point3D);

					point3D[1]+=-1e-8;
					assert(("derivative",fabs(ret[1]-5.e7*(x1-x0))<1e-5||(L2Norm(ret)>1&&fabs(ret[1]-5.e7*(x1-x0))<1e-5*L2Norm(ret))));
					assert(("derivative",fabs(test->evalDeriv1(point3D)-5.e7*(x1-x0))<1e-5||(L2Norm(ret)>1&&fabs(ret[1]-5.e7*(x1-x0))<1e-5*L2Norm(ret))));

					point3D[2]+=-1.e-8;
					x0=test->eval(point3D);
					point3D[2]+=2.e-8;
					x1=test->eval(point3D);

					point3D[2]+=-1e-8;
					assert(("derivative",fabs(ret[2]-5.e7*(x1-x0))<1e-5||(L2Norm(ret)>1&&fabs(ret[2]-5.e7*(x1-x0))<1e-5*L2Norm(ret))));
					assert(("derivative",fabs(test->evalDeriv2(point3D)-5.e7*(x1-x0))<1e-5||(L2Norm(ret)>1&&fabs(ret[2]-5.e7*(x1-x0))<1e-5*L2Norm(ret))));

					ret.resize(1);
				}
			}
		}
	}

	delete all3DbasisFunctions;

	all3DbasisFunctions = Utilities::createDGBasisFunctionSet3DH1ConformingPrism(5);
	for(int i=0;i<all3DbasisFunctions->size();++i){
		const Base::BaseBasisFunction* test = (*all3DbasisFunctions)[i];
		for(point3D[0]=-1.5;point3D[0]<1.51;point3D[0]+=0.25){
			for(point3D[1]=-1.5;point3D[1]<1.51;point3D[1]+=0.25){
				for(point3D[2]=-1.5;point3D[2]<1.51;point3D[2]+=0.25){
					test->eval(point3D,ret);
					assert(("eval",(test->eval(point3D)-ret[0])<1e-12));

					point3D[0]+=-1.e-8;
					double x0=test->eval(point3D);
					point3D[0]+=2.e-8;
					double x1=test->eval(point3D);

					point3D[0]+=-1e-8;
					ret.resize(3);
					test->evalDeriv(point3D,ret);
					assert(("derivative",fabs(ret[0]-5.e7*(x1-x0))<1e-5||(L2Norm(ret)>1&&fabs(ret[0]-5.e7*(x1-x0))<1e-5*L2Norm(ret))));
					assert(("derivative",fabs(test->evalDeriv0(point3D)-5.e7*(x1-x0))<1e-5||(L2Norm(ret)>1&&fabs(ret[0]-5.e7*(x1-x0))<1e-5*L2Norm(ret))));

					point3D[1]+=-1.e-8;
					x0=test->eval(point3D);
					point3D[1]+=2.e-8;
					x1=test->eval(point3D);

					point3D[1]+=-1e-8;
					assert(("derivative",fabs(ret[1]-5.e7*(x1-x0))<1e-5||(L2Norm(ret)>1&&fabs(ret[1]-5.e7*(x1-x0))<1e-5*L2Norm(ret))));
					assert(("derivative",fabs(test->evalDeriv1(point3D)-5.e7*(x1-x0))<1e-5||(L2Norm(ret)>1&&fabs(ret[1]-5.e7*(x1-x0))<1e-5*L2Norm(ret))));

					point3D[2]+=-1.e-8;
					x0=test->eval(point3D);
					point3D[2]+=2.e-8;
					x1=test->eval(point3D);

					point3D[2]+=-1e-8;
					assert(("derivative",fabs(ret[2]-5.e7*(x1-x0))<1e-5||(L2Norm(ret)>1&&fabs(ret[2]-5.e7*(x1-x0))<1e-5*L2Norm(ret))));
					assert(("derivative",fabs(test->evalDeriv2(point3D)-5.e7*(x1-x0))<1e-5||(L2Norm(ret)>1&&fabs(ret[2]-5.e7*(x1-x0))<1e-5*L2Norm(ret))));

					ret.resize(1);
				}
			}
		}
	}

	delete all3DbasisFunctions;

	return 0;
}



