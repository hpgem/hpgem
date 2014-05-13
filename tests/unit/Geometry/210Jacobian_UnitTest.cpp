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

//naming convention: <Digit><ClassName>_UnitTest.cpp where <Digit> is a number that will make sure
//the unit tests are ordered such that the first failing unit test indicate the culprit class and
//other 'unit' tests may assume correct execution of all prior unit tests

#include "Geometry/Jacobian.hpp"
#include "cassert"

using Geometry::Jacobian;

int main(){

	//square jacobians and determinants

	Jacobian dim1(1,1);

	dim1(0,0)=3.14;
	assert(("1D determinant",fabs(dim1.determinant()-3.14)<1e-12));

	dim1(0,0)=-2.81;
	assert(("1D determinant",fabs(dim1.determinant()+2.81)<1e-12));

	dim1(0,0)=0;
	assert(("1D determinant",fabs(dim1.determinant())<1e-12));

	Jacobian dim2(2,2);

	dim2(0,0)=1.;
	dim2(1,0)=2.;
	dim2(0,1)=3.;
	dim2(1,1)=4.;
	assert(("2D determinant",fabs(dim2.determinant()+2.)<1e-12));

	dim2(0,0)=-1.;
	dim2(1,0)=2.;
	dim2(0,1)=3.;
	dim2(1,1)=4.;
	assert(("2D determinant",fabs(dim2.determinant()+10.)<1e-12));

	Jacobian copy(dim2);

	dim2(0,0)=1.;
	dim2(1,0)=-2.;
	dim2(0,1)=3.;
	dim2(1,1)=4.;
	assert(("2D determinant",fabs(dim2.determinant()-10.)<1e-12));

	dim2(0,0)=1.;
	dim2(1,0)=2.;
	dim2(0,1)=-3.;
	dim2(1,1)=4.;
	assert(("2D determinant",fabs(dim2.determinant()-10.)<1e-12));

	dim2(0,0)=1.;
	dim2(1,0)=2.;
	dim2(0,1)=3.;
	dim2(1,1)=-4.;
	assert(("2D determinant",fabs(dim2.determinant()+10.)<1e-12));

	dim2(0,0)=-1.;
	dim2(1,0)=-2.;
	dim2(0,1)=3.;
	dim2(1,1)=4.;
	assert(("2D determinant",fabs(dim2.determinant()-2.)<1e-12));

	dim2(0,0)=-1.;
	dim2(1,0)=2.;
	dim2(0,1)=-3.;
	dim2(1,1)=4.;
	assert(("2D determinant",fabs(dim2.determinant()-2.)<1e-12));

	dim2(0,0)=-1.;
	dim2(1,0)=2.;
	dim2(0,1)=3.;
	dim2(1,1)=-4.;
	assert(("2D determinant",fabs(dim2.determinant()+2.)<1e-12));

	dim2(0,0)=1.;
	dim2(1,0)=2.;
	dim2(0,1)=2.;
	dim2(1,1)=4.;
	assert(("2D determinant",fabs(dim2.determinant())<1e-12));
	assert(("copy constructor",fabs(copy.determinant()+10.)<1e-12));

	Jacobian product(2,2);

	dim2.multiplyJacobiansInto(copy,product);
	assert(("multiply JacobiansInto - square matrixes",product.getNCols()==2&&product.getNRows()==2&&fabs(product.determinant())<1e-12));

	Jacobian dim3(3,3);
	dim3(0,0)=1.;
	dim3(1,0)=2.;
	dim3(2,0)=3.;
	dim3(0,1)=4.;
	dim3(1,1)=5.;
	dim3(2,1)=6.;
	dim3(0,2)=7.;
	dim3(1,2)=8.;
	dim3(2,2)=9.;
	assert(("3D determinant",fabs(dim3.determinant())<1e-12));

	dim3(0,0)=-1.;
	dim3(1,0)=2.;
	dim3(2,0)=3.;
	dim3(0,1)=4.;
	dim3(1,1)=5.;
	dim3(2,1)=6.;
	dim3(0,2)=7.;
	dim3(1,2)=8.;
	dim3(2,2)=9.;
	assert(("3D determinant",fabs(dim3.determinant()-6.)<1e-12));

	dim3(0,0)=1.;
	dim3(1,0)=-2.;
	dim3(2,0)=3.;
	dim3(0,1)=4.;
	dim3(1,1)=5.;
	dim3(2,1)=6.;
	dim3(0,2)=7.;
	dim3(1,2)=8.;
	dim3(2,2)=9.;
	assert(("3D determinant",fabs(dim3.determinant()+24.)<1e-12));

	dim3(0,0)=1.;
	dim3(1,0)=2.;
	dim3(2,0)=-3.;
	dim3(0,1)=4.;
	dim3(1,1)=5.;
	dim3(2,1)=6.;
	dim3(0,2)=7.;
	dim3(1,2)=8.;
	dim3(2,2)=9.;
	assert(("3D determinant",fabs(dim3.determinant()-18.)<1e-12));

	dim3(0,0)=1.;
	dim3(1,0)=2.;
	dim3(2,0)=3.;
	dim3(0,1)=-4.;
	dim3(1,1)=5.;
	dim3(2,1)=6.;
	dim3(0,2)=7.;
	dim3(1,2)=8.;
	dim3(2,2)=9.;
	assert(("3D determinant",fabs(dim3.determinant()+48.)<1e-12));

	dim3(0,0)=1.;
	dim3(1,0)=2.;
	dim3(2,0)=3.;
	dim3(0,1)=4.;
	dim3(1,1)=-5.;
	dim3(2,1)=6.;
	dim3(0,2)=7.;
	dim3(1,2)=8.;
	dim3(2,2)=9.;
	assert(("3D determinant",fabs(dim3.determinant()-120.)<1e-12));

	dim3(0,0)=1.;
	dim3(1,0)=2.;
	dim3(2,0)=3.;
	dim3(0,1)=4.;
	dim3(1,1)=5.;
	dim3(2,1)=-6.;
	dim3(0,2)=7.;
	dim3(1,2)=8.;
	dim3(2,2)=9.;
	assert(("3D determinant",fabs(dim3.determinant()+72.)<1e-12));

	dim3(0,0)=1.;
	dim3(1,0)=2.;
	dim3(2,0)=3.;
	dim3(0,1)=4.;
	dim3(1,1)=5.;
	dim3(2,1)=6.;
	dim3(0,2)=-7.;
	dim3(1,2)=8.;
	dim3(2,2)=9.;
	assert(("3D determinant",fabs(dim3.determinant()-42.)<1e-12));

	dim3(0,0)=1.;
	dim3(1,0)=2.;
	dim3(2,0)=3.;
	dim3(0,1)=4.;
	dim3(1,1)=5.;
	dim3(2,1)=6.;
	dim3(0,2)=7.;
	dim3(1,2)=-8.;
	dim3(2,2)=9.;
	assert(("3D determinant",fabs(dim3.determinant()+96.)<1e-12));

	dim3(0,0)=1.;
	dim3(1,0)=2.;
	dim3(2,0)=3.;
	dim3(0,1)=4.;
	dim3(1,1)=5.;
	dim3(2,1)=6.;
	dim3(0,2)=7.;
	dim3(1,2)=8.;
	dim3(2,2)=-9.;
	assert(("3D determinant",fabs(dim3.determinant()-54.)<1e-12));

	Jacobian rectangle21(2,1);
	Jacobian rectangle12(1,2);
	Jacobian rectangle32(3,2);
	Jacobian rectangle23(2,3);

	rectangle21(0,0)=1.;
	rectangle21(0,1)=2.;
	rectangle12(0,0)=3.;
	rectangle12(1,0)=4.;

	rectangle32(0,0)=1.;
	rectangle32(0,1)=2.;
	rectangle32(0,2)=3.;
	rectangle32(1,0)=4.;
	rectangle32(1,1)=5.;
	rectangle32(1,2)=6.;
	rectangle23(0,0)=7.;
	rectangle23(0,1)=8.;
	rectangle23(1,0)=9.;
	rectangle23(1,1)=10.;
	rectangle23(2,0)=11.;
	rectangle23(2,1)=12.;

	product.resize(2,1);//not sure if Jacobians are meant to be able to do this

	dim1.multiplyJacobiansInto(rectangle21,product);
	assert(("multiply JacobiansInto - square & rectangular matrixes",product.getNCols()==1&&product.getNRows()==2));
	rectangle21.multiplyJacobiansInto(dim2,product);
	assert(("multiply JacobiansInto - square & rectangular matrixes",product.getNCols()==1&&product.getNRows()==2));

	product.resize(1,2);

	dim2.multiplyJacobiansInto(rectangle12,product);
	assert(("multiply JacobiansInto - square & rectangular matrixes",product.getNCols()==2&&product.getNRows()==1));
	rectangle12.multiplyJacobiansInto(dim1,product);
	assert(("multiply JacobiansInto - square & rectangular matrixes",product.getNCols()==2&&product.getNRows()==1));

	product.resize(2,3);

	dim3.multiplyJacobiansInto(rectangle23,product);
	assert(("multiply JacobiansInto - square & rectangular matrixes",product.getNCols()==3&&product.getNRows()==2));
	rectangle23.multiplyJacobiansInto(dim2,product);
	assert(("multiply JacobiansInto - square & rectangular matrixes",product.getNCols()==3&&product.getNRows()==2));

	product.resize(3,2);

	dim2.multiplyJacobiansInto(rectangle32,product);
	assert(("multiply JacobiansInto - square & rectangular matrixes",product.getNCols()==2&&product.getNRows()==3));
	rectangle32.multiplyJacobiansInto(dim3,product);
	assert(("multiply JacobiansInto - square & rectangular matrixes",product.getNCols()==2&&product.getNRows()==3));

	product.resize(1,3);

	rectangle21.multiplyJacobiansInto(rectangle32,product);
	assert(("multiply JacobiansInto - rectangular matrixes",product.getNCols()==3&&product.getNRows()==1));

	product.resize(3,1);

	rectangle23.multiplyJacobiansInto(rectangle12,product);
	assert(("multiply JacobiansInto - rectangular matrixes",product.getNCols()==1&&product.getNRows()==3));

	product.resize(1,1);

	rectangle21.multiplyJacobiansInto(rectangle12,product);
	assert(("multiply JacobiansInto - rectangular matrixes",product.getNCols()==1&&product.getNRows()==1&&fabs(product.determinant()-11.)<1e-12));

	product.resize(2,2);

	rectangle12.multiplyJacobiansInto(rectangle21,product);
	assert(("multiply JacobiansInto - rectangular matrixes",product.getNCols()==2&&product.getNRows()==2&&fabs(product.determinant())<1e-12));
	rectangle32.multiplyJacobiansInto(rectangle23,product);
	assert(("multiply JacobiansInto - rectangular matrixes",product.getNCols()==2&&product.getNRows()==2&&fabs(product.determinant()-36.)<1e-12));

	product.resize(3,3);

	rectangle23.multiplyJacobiansInto(rectangle32,product);
	assert(("multiply JacobiansInto - rectangular matrixes",product.getNCols()==3&&product.getNRows()==3&&fabs(product.determinant())<1e-12));

	return 0;
}

