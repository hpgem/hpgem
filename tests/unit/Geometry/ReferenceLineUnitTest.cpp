//
//  ReferenceLineUnitTest.cpp
//  
//
//  Created by Shavarsh Nurijanyan on 2/4/13.
//
//

#include "Geometry/ReferencePoint.hpp"
#include "Geometry/ReferenceLine.hpp"
#include "Geometry/Mappings/MappingToRefLineToLine.hpp"
using namespace std;

int main ()
{
    cout << "Geometry::ReferenceLine line;" << endl;

    Geometry::ReferenceLine line;
    
    cout << "line.getName() "<< line.getName() << endl;

    cout << "cout << line: ";
    cout << line;

    Geometry::PointReference<1> pCenter;
    line.getCenter(pCenter);
    cout << "line.getCenter(Geometry::PointReference<1>): " << pCenter << endl;

    Geometry::PointReference<1> p;
    line.getNode(0, p);
    cout << "line.getNode(0, p): " << p << endl;
    
    line.getNode(1, p);
    cout << "line.getNode(1, p): " << p << endl;
    
    p[0] = 1.0;
    cout << "line.isInternalPoint({1}): " << line.isInternalPoint(p) << endl;

    p[0] = 0.0;
    cout << "line.isInternalPoint({0}): " << line.isInternalPoint(p) << endl;

    p[0] = -1.0;
    cout << "line.isInternalPoint({-1}): " << line.isInternalPoint(p) << endl;

    p[0] = -2.0;
    cout << "line.isInternalPoint({-2}): " << line.isInternalPoint(p) << endl;

    cout << "line.getNrOfCodim1Entities(): " << line.getNrOfCodim1Entities() << " nodes" << endl;

    int a1[2]; int a2[2]; int map;
    /* Lines:
     *
     * 0-------1 0-------1
     *
     * Should give identity.
     *
     */
    a1[0] = 0; a1[1] = 1;
    a2[0] = 0; a2[1] = 1;
    const std::vector<unsigned int> globalNodeIndexes1(a1, a1 + sizeof(a1)/sizeof(a1[0]));
    const std::vector<unsigned int> globalNodeIndexes2(a2, a2 + sizeof(a2)/sizeof(a2[0]));

    map = line.getCodim0MappingIndex(globalNodeIndexes1,globalNodeIndexes2);
    cout << "line.getCodim0MappingIndex() = " << map << endl;

    /* Lines:
     *
     * 0-------1 1-------0
     *
     * Mirroring (the only other line mapping)
     *
     */
    a1[0] = 0; a1[1] = 1;
    a2[0] = 1; a2[1] = 0;
    const std::vector<unsigned int> globalNodeIndexes3(a1, a1 + sizeof(a1)/sizeof(a1[0]));
    const std::vector<unsigned int> globalNodeIndexes4(a2, a2 + sizeof(a2)/sizeof(a2[0]));

    map = line.getCodim0MappingIndex(globalNodeIndexes3,globalNodeIndexes4);
    cout << "line.getCodim0MappingIndex() = " << map << endl;

    cout << "line.getCodim0MappingPtr(0)" <<  endl;
    const Geometry::MappingReferenceToReference<1,1>* mappingReferenceToReference11;
    mappingReferenceToReference11 = line.getCodim0MappingPtr(0);
    cout << "line.getCodim0MappingPtr(1)" <<  endl;
    mappingReferenceToReference11 = line.getCodim0MappingPtr(1);


    cout << "line.getCodim1MappingPtr(0)" <<  endl;
    const Geometry::MappingReferenceToReference<0,1>* mappingReferenceToReference01;
    mappingReferenceToReference01 = line.getCodim1MappingPtr(0);
    cout << "line.getCodim1MappingPtr(1)" <<  endl;
    mappingReferenceToReference01 = line.getCodim1MappingPtr(1);
}
