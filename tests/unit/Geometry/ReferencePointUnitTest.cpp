//
//  ReferencePointUnitTest.cpp
//  
//
//  Created by Shavarsh Nurijanyan on 2/4/13.
//
//

#include "Geometry/PointReference.hpp"
#include "Geometry/ReferencePoint.hpp"
#include <vector>

int main ()
{
    cout << "Geometry::ReferencePoint point;" << endl;

    Geometry::ReferencePoint point;

    cout << "point.getName() "<< point.getName() << endl;

    Geometry::PointReference p(0);
    cout << "point.isInternalPoint() "<< point.isInternalPoint(p) << endl;
    
    std::vector<unsigned int> gi1;
    gi1.resize(1);
    gi1[0] = 1;
    std::vector<unsigned int> gi2;
    gi2.resize(1);
    gi2[0] = 1;
    cout << "point.getCodim0MappingIndex() "<< point.getCodim0MappingIndex(gi1,gi2) << endl;
    
    unsigned int i = 0;
    cout << "point.getCodim0MappingPtr() "<< point.getCodim0MappingPtr(i) << endl;
    
    return 0;
}
