#include "Geometry/PointReference.hpp"
#include "Geometry/Point.hpp"

#include <iostream>
using namespace std;

int main()
{
    const unsigned int dim = 2;
    Geometry::PointReference p(dim);
    Geometry::PointReference q(p);
    p[0] = 1.;  p[1] = 1.8;
    
    Geometry::PointReference r(p);
    Geometry::PointReference s(dim);
    s = p;      // assign PointT into PointT?


    Geometry::PointReference::VectorOfCoordsT t(2);
    t[0] = 9.;  t[1] = 1.;
    Geometry::PointReference u(t);
    cout<<"9, 1 "<<u<<endl;
    Geometry::PointReference pp(t);
    pp[0] = 100.;
    cout<<pp<<endl;
    
    pp = p;     // assign PointReferenceT into PointT?
  
    cout<<pp<<endl;
    cout<<"COUTS"<<endl;
    cout << p << q << r << s << u << pp << endl;
    
        
    
    cout << pp <<"+="<<q << endl;
    pp += q;
    cout << pp<<endl;
    cout <<"______________________"<<endl;
    cout << r<<"-=" <<s<< endl;
    r -= s;
    cout << r<<endl;
    
        //    cout << "Point = Point + double" << endl;
        //    p3 = p1+2.0;
        //    cout << p3;
        //    cout << "Point = Point - double" << endl;
        //    p3 = p1-2.0;
        //    cout << p3<<endl;
    cout <<"______________________"<<endl;
    cout << u<< " * "<<22.0 << endl;
    s = (Geometry::PointReference)u*22.0;
    cout << s<<endl;
        //    cout << "Point += double" << endl;
        //    p3 += 2.0;
        //    cout << p3;
        //    cout << "Point -= double" << endl;
        //    p3 -= 2.0;
        //    cout << p3;
    cout <<"______________________"<<endl;
    cout << s<< " *="<<2.0 << endl;
    s *= 2.0;
    cout << s<<endl;


  return 0;
}
