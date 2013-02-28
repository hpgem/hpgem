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
    
    B3=BB1*B2;
    
    //Test wedge product both quick form if vector exist and create vector form.
    BB1.computeWedgeStuffVector(B2);
    B3=BB1.computeWedgeStuffVector();
    
    //Matrix assigment test
    AA4=2.0;
    
    //Matrix times matrix test
    LinearAlgebra::Matrix CC1(2,3,2);
    LinearAlgebra::Matrix CC2(3,2,2);
    LinearAlgebra::Matrix CC3(2,2);
    
    CC3=CC1*CC2;
    
    
    
    
    CC1.axpy(2.0,CC2);
    
    cout << CC1;
    
    
    
    

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


