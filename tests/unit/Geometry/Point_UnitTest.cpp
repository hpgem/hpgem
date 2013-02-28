//
//  NodeUnitTest.cpp
//  
//
//  Created by Shavarsh Nurijanyan on 2/3/13.
//
//
#include "Geometry/Point.hpp"

using namespace std;


int main()
{

    Geometry::Point<3> point0;
    cout << "Default ctr"<<point0<<endl;
    cout <<"Size of ="<<point0<< "is "<< point0.size()<<endl;

    
    cout << "Point<3>() vector constructor ((double []){1.2, 1.3, 1.4})" << endl;
    double a[]={1.2, 1.3, 1.4};
    
   
      cout << "Works :)" << endl;
    
    Geometry::Point<3> point(a);
    cout << point<<endl;
    cout << "Works :)" << endl;
    
    cout << "Point() copy constructor" << endl;
    Geometry::Point<3> pointCopyConstructor(point);
    cout << pointCopyConstructor;
    cout << "Works :)" << endl;
    
    cout << "Point<2>() vector constructor ((double []){1.2, 1.3, 1.4})" << endl;
    //Warning! This way Point class will truncate and take two points, because of the dimension
    Geometry::Point<2> point2D((double []){1.01, 1.02, 1.5});
    cout << "Works :)" << endl;
    
    cout << point2D<<endl;

    cout << "setCoordinate()" << endl;
    Geometry::Point<2> pointIndepSet;
    pointIndepSet.setCoordinate(0, 1.888);
    pointIndepSet.setCoordinate(1, 1.777);
    pointIndepSet.setCoordinate(2, 100000);
    cout << pointIndepSet;
    cout << "Works :)" << endl;

    cout << "getCoordinate()" << endl;
    Geometry::Point<2> pointIndepGet;
    cout << pointIndepGet.getCoordinate(0) << " "
         << pointIndepGet.getCoordinate(1) << " "
         << pointIndepGet.getCoordinate(2);
    cout << "Works :)" << endl;
    
    cout << "[] Operator" << endl;
    Geometry::Point<4> outerPoint4D;
    outerPoint4D[0]=1.0001;
    outerPoint4D[1]=1.00041;
    outerPoint4D[2]=1.00301;
    outerPoint4D[3]=1.02001;
    cout << outerPoint4D<<endl;
    cout << "Works :)" << endl;
    
    cout << "Operators:" << endl;
    Geometry::Point<3> p1((double []){1.0,2.0,3.0});
    Geometry::Point<3> p2((double []){1.0,2.0,3.0});
    Geometry::Point<3> p3;
    cout << p1<<endl;
    cout << p2<<endl;
//    cout << "Point + Point" << endl;
//    p3 = p1 + p2;
//    cout << p3;
//    cout << "Point - Point" << endl;
//    p3 = p1 - p2;
//    cout << p3;
    
    cout << p3<<endl;
    cout <<"______________________"<<endl;
    
    cout << p1 <<"-"<<p2 << endl;
    cout << (p1-p2)<<endl;
    
    cout << p1 <<"+="<<p2 << endl;
    p1 += p2;
    cout << p1<<endl;
    
    cout <<"______________________"<<endl;
    cout << p1<<"-=" <<p2<< endl;
    p1 -= p2;
    cout << p1<<endl;
   
//    cout << "Point = Point + double" << endl;
//    p3 = p1+2.0;
//    cout << p3;
//    cout << "Point = Point - double" << endl;
//    p3 = p1-2.0;
//    cout << p3<<endl;
    cout <<"______________________"<<endl;
    cout << p1 << " * "<<22.0 << endl;
    p3 = p1*22.0;
    cout << p3<<endl;
//    cout << "Point += double" << endl;
//    p3 += 2.0;
//    cout << p3;
//    cout << "Point -= double" << endl;
//    p3 -= 2.0;
//    cout << p3;
    cout <<"______________________"<<endl;
    cout << p3<< " *="<<2.0 << endl;
    p3 *= 2.0;
    cout << p3<<endl;


    return 0;
}
