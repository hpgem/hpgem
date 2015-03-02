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

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "LinearAlgebra/Matrix.hpp"
#include "LinearAlgebra/NumericalVector.hpp"


using namespace::std;

int main(int argc, char* argv[])
{
  
    //First test the empty constructor
    LinearAlgebra::Matrix AA1;
    
    //Second test the constuctor which defines the size
    LinearAlgebra::Matrix AA2(2,3);
    
    //Third test the constuctor to initiate to a set value e.g. 2 by 3 matrix with all values initised to 1.5
    LinearAlgebra::Matrix AA3(2,3,1.5);
    
    //Check the copy constructor
    LinearAlgebra::Matrix AA4(AA3);
    
    //Check resizing the matrix
    AA4.resize(2,5);
    
    //Check the Matrix assigment
    AA2=AA3;
    
    //Now test Matrix times vector
    LinearAlgebra::Matrix BB1(2,2);
    BB1(0,0)=1.0;
    BB1(0,1)=2.0;
    BB1(1,0)=3.0;
    BB1(1,1)=4.0;
    
    
    cout << "This is BB1 \n" << BB1 << "\n"; 
    
    //now test divide
    
    
    BB1/=2.0;
    cout << "This is BB1 divided by 2 \n" << BB1 << "\n";
    
    cout << "and by 2.0 again but using inline divide \n" << BB1/2.0 << "\n";
    
    //Now test output
    
    cout << BB1 <<"\n";
    
    cout << AA3 <<endl;
   
    LinearAlgebra::NumericalVector B2(2);
    LinearAlgebra::NumericalVector B3;
    B2(0)=1.0;
    B2(1)=2.0;
    
    cout<<"axpy"<<endl;
    B3=BB1*B2;
    
    cout<<"axpy"<<endl;
    //Test wedge product both quick form if vector exist and create vector form.
    B2=BB1.computeWedgeStuffVector();
    B3=BB1.computeWedgeStuffVector();
    
    //Matrix assigment test
    AA4=2.0;
    
    //Matrix times matrix test
    LinearAlgebra::Matrix CC1(2,3,2);
    LinearAlgebra::Matrix CC2(3,2,2);
    LinearAlgebra::Matrix CC3(2,2);
    
    CC3=CC1*CC2;
    
    
    cout<<"axpy"<<endl;
    
    //CC1.axpy(2.0,CC2);
    
    
    cout << CC1 << std::endl;
    
    LinearAlgebra::Matrix CC2_fix(3,2, 4);
    std::cout << "CC2_fix PRE" << std::endl;
    std::cout << CC2_fix << std::endl;
    
    CC2.axpy(2.0, CC2_fix);
    
    std::cout << "CC2_fix POST" << std::endl;
    std::cout << CC2_fix << std::endl;
    
    
    std::cout << CC2 << std::endl;
    

    CC3=BB1;
    CC3=BB1.LUfactorisation();
    
    cout << "\n Now the LU factorisation \n";
    
    cout << "Before : \n";
    
    cout << BB1;
     
    cout << "After : \n";
    cout << CC3;
    
    cout << "\n Now the inverse"  <<endl;
    
    CC3(0,0)=0.8147;
    CC3(0,1)=0.1270;
    CC3(1,0)=0.9058;
    CC3(1,1)=0.9134;
    
    LinearAlgebra::Matrix CC4(2,2);
    
    CC4 = CC3.inverse();
    
    cout << CC4; 
    
    CC3=CC3*CC4;
    
    cout << CC3;
    
    
    cout << "\n Now test the solution of Ax=B \n ";
    
    CC3.solve(CC4);
    
    cout << CC4;
    
    cout << "\n Now test the inverse of a 3 by3 \n";
    
    LinearAlgebra::Matrix DD(3,3);
    
    
    DD(0,0)=0.6324;
    DD(0,1)=0.5469;
    DD(0,2)=0.1576;
    DD(1,0)=0.0975;
    DD(1,1)=0.9575;
    DD(1,2)=0.9706;
    DD(2,0)=0.2785;
    DD(2,1)=0.9649;
    DD(2,2)=0.9572;
    
    LinearAlgebra::Matrix ans(3,3);
    
    ans = DD.inverse();
    
    cout << ans << std::endl;
    
    
    return 0;

}
	
	
//	LinearAlgebra::Matrix<double> AA(3,2);
//	
//	LinearAlgebra::Matrix<double> BB(2,3);
//	
//	LinearAlgebra::Matrix<double> CC(2,4);
//	
//	LinearAlgebra::Matrix<double> DD(3,2);
//	
//	LinearAlgebra::Matrix<double> ANS;
//    
//    LinearAlgebra::NumericalVector x(2);
//    LinearAlgebra::NumericalVector y(2);
//    
//    LinearAlgebra::Matrix<double> A(2,2);
//    
//    
//	x(0)=1.0;
//    x(1)=2.0;
//    
//    A(0,0)=1.0;
//    A(0,1)=2.0;
//    A(1,0)=3.0;
//    A(1,1)=4.0;
//	
//	AA(0,0)=1;
//	AA(0,1)=2;
//	AA(1,0)=3;
//	AA(1,1)=4;
//	AA(2,0)=5;
//	AA(2,1)=6;
//    
//    cout << x;
//    
//    y=AA*x;
//    
//    cout << "Help" << endl;
//    cout << y;
//    cout << "and now " <<endl;
//	
//	BB(0,0)=7;
//	BB(0,1)=8;
//	BB(0,2)=9;
//	BB(1,0)=10;
//	BB(1,1)=11;
//	BB(1,2)=12;
//	
//	CC(0,0)=7;
//	CC(0,1)=8;
//	CC(0,2)=9;
//	CC(0,3)=10;
//	CC(1,0)=11;
//	CC(1,1)=12;
//	CC(1,2)=13;
//	CC(1,3)=14;
//	
//	DD(0,0)=6;
//	DD(0,1)=5;
//	DD(1,0)=4;
//	DD(1,1)=3;
//	DD(2,0)=2;
//	DD(2,1)=1;
//	
//	//AA*=BB;
//	cout << AA;
//	//AA*=CC;
//	cout << AA;
//	
//	ANS=AA+DD;
//	cout << ANS;
//	//cout << DD;
//	//cout << AA;
//	
//	
// 	;


