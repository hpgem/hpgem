#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "LinearAlgebra/NumericalVector.hpp"

          
using namespace::std;

int main(int argc, char* argv[])
{
    
    LinearAlgebra::NumericalVector a(3);
    LinearAlgebra::NumericalVector b(3);
    
    cout << a.size();
    
    
    a[0]=2.0;
    a[1]=3.0;
    a[2]=1.0;
    
    cout << a << endl;
    
    b=3.0*a;
    
    cout << b << endl;
    
    a.axpy(3.0,b);
    cout << a;
    
   /** /*LinearAlgebra::NumericalVector c;
    Geometry::Point<2> p1 ((double []){1.0,2.0});
    Geometry::Point<2> p2 ((double []){1.0,3.0});
    Geometry::Point<2> p3;
    
     
    p3=p1;
	*/
    /*
    	
	LinearAlgebra::NumericalVector AA(3);
	
	LinearAlgebra::NumericalVector BB(3);
    
    LinearAlgebra::NumericalVector CC(3);
	
    LinearAlgebra::NumericalVector EE;
	
    FunnyVector DD(3);
	
	
	AA(0)=1;
	AA(1)=2;
	AA(2)=3;
	
	BB(0)=7;
	BB(1)=8;
	BB(2)=9;
    
    cout << AA <<endl;

    DD(0)=1;
	DD(1)=2;
	DD(2)=3;
    
    DD=AA;
    
    
    AA+=BB;
    
    cout << AA<<endl;
   
    CC=AA+BB;
    
    cout << CC << endl;
    
    */
	
}

