//
//  ReferenceSquareUnitTest.cpp
//  
//
//  Created by Shavarsh Nurijanyan on 2/5/13.
//
//
#include "Geometry/ReferencePoint.hpp"
#include "Geometry/ReferenceLine.hpp"
#include "Geometry/ReferenceSquare.hpp"
#include "Geometry/Mappings/MappingToRefSquareToSquare.hpp"
using namespace std;

int main ()
{

    cout << "Geometry::ReferenceSquare square;" << endl;

    Geometry::ReferenceSquare square;
    
    cout << "square.getName() "<< square.getName() << endl;

    cout << "cout << square: ";
    cout << square;

    Geometry::PointReference pCenter(2);
    square.getCenter(pCenter);
    cout << "square.getCenter(Geometry::Point<2>): " << pCenter;
    
    Geometry::PointReference p(2);
    square.getNode(0, p);
    cout << "square.getVertexPoint(0, p): " << p;
    
    square.getNode(1, p);
    cout << "square.getVertexPoint(1, p): " << p;
    
    square.getNode(2, p);
    cout << "square.getVertexPoint(2, p): " << p;
    
    square.getNode(3, p);
    cout << "square.getVertexPoint(3, p): " << p;
    
    cout << "square.isInternalPoint({0,0}): " << square.isInternalPoint(pCenter) << endl;
    
    Geometry::PointReference pOut(2);
    
    pOut[0] = -2.0;
    pOut[1] = -2.0;
    
    cout << "square.isInternalPoint({-2.0,-2.0}): " << square.isInternalPoint(pOut) << endl;

    cout << "square.getNrOfCodim1Entities(): " << square.getNrOfCodim1Entities() << " lines" << endl;

    cout << "square.getNrOfCodim2Entities(): " << square.getNrOfCodim2Entities() << " points" << endl;
    
    std::vector<unsigned int> indexes;

    for (int face=0; face < square.getNrOfCodim1Entities(); ++face)
    {
        square.getCodim1EntityLocalIndices(face, indexes);
        cout << "square.getCodim1EntityLocalIndices("<<face<<"): (";
        for (int i = 0; i < indexes.size(); ++i)
        {
             if (i < indexes.size()-1) cout << indexes[i] << ",";
             else cout << indexes[i];
        }
        cout << ")" << endl;
    }

    const Geometry::ReferenceGeometry* referenceGeometry;
    for (int face=0; face < square.getNrOfCodim1Entities(); ++face)
    {
        referenceGeometry = square.getCodim1ReferenceGeometry(face);
        cout << "square.getCodim1ReferenceGeometry("<<face<<"): " << referenceGeometry->getName() << endl;
    }

    const Geometry::MappingReferenceToReference* ref2ToRef2Mapping;
    for (int index=0; index < 8; ++index)
    {
        ref2ToRef2Mapping = square.getCodim0MappingPtr(index);
    }

    const Geometry::MappingReferenceToReference* ref1ToRef2Mapping;
    for (int face=0; face < square.getNrOfCodim1Entities(); ++face)
    {
        ref1ToRef2Mapping = square.getCodim1MappingPtr(face);
    }

    int a1[4]; int a2[4]; int map;
    /* Squares:
     *
     * 2-------3 2-------3
     * |       | |       |
     * |       | |       |
     * |       | |       |
     * 0-------1 0-------1
     *
     * Should give identity.
     *
     */
    a1[0] = 0; a1[1] = 1; a1[2] = 2; a1[3] = 3;
    a2[0] = 0; a2[1] = 1; a2[2] = 2; a2[3] = 3;
    const std::vector<unsigned int> globalNodeIndexes1(a1, a1 + sizeof(a1)/sizeof(a1[0]));
    const std::vector<unsigned int> globalNodeIndexes2(a2, a2 + sizeof(a2)/sizeof(a2[0]));

    map = square.getCodim0MappingIndex(globalNodeIndexes1,globalNodeIndexes2);
    cout << "square.getCodim0MappingIndex() = " << map << endl;

    /* Squares:
     *
     * 2-------3 1-------3
     * |       | |       |
     * |       | |       |
     * |       | |       |
     * 0-------1 0-------2
     *
     * Should give x=y reflection.
     *
     */
    a1[0] = 0; a1[1] = 1; a1[2] = 2; a1[3] = 3;
    a2[0] = 0; a2[1] = 2; a2[2] = 1; a2[3] = 3;
    const std::vector<unsigned int> globalNodeIndexes3(a1, a1 + sizeof(a1)/sizeof(a1[0]));
    const std::vector<unsigned int> globalNodeIndexes4(a2, a2 + sizeof(a2)/sizeof(a2[0]));

    map = square.getCodim0MappingIndex(globalNodeIndexes3,globalNodeIndexes4);
    cout << "square.getCodim0MappingIndex() = " << map << endl;

    /* Squares:
     *
     * 2-------3 3-------2
     * |       | |       |
     * |       | |       |
     * |       | |       |
     * 0-------1 1-------0
     *
     * Should give x=0 reflection.
     *
     */
    a1[0] = 0; a1[1] = 1; a1[2] = 2; a1[3] = 3;
    a2[0] = 1; a2[1] = 0; a2[2] = 3; a2[3] = 2;
    const std::vector<unsigned int> globalNodeIndexes5(a1, a1 + sizeof(a1)/sizeof(a1[0]));
    const std::vector<unsigned int> globalNodeIndexes6(a2, a2 + sizeof(a2)/sizeof(a2[0]));

    map = square.getCodim0MappingIndex(globalNodeIndexes5,globalNodeIndexes6);
    cout << "square.getCodim0MappingIndex() = " << map << endl;

    /* Squares:
     *
     * 2-------3 3-------1
     * |       | |       |
     * |       | |       |
     * |       | |       |
     * 0-------1 2-------0
     *
     * Should give 90¼ rotation (clockwise).
     *
     */
    a1[0] = 0; a1[1] = 1; a1[2] = 2; a1[3] = 3;
    a2[0] = 2; a2[1] = 0; a2[2] = 3; a2[3] = 1;
    const std::vector<unsigned int> globalNodeIndexes7(a1, a1 + sizeof(a1)/sizeof(a1[0]));
    const std::vector<unsigned int> globalNodeIndexes8(a2, a2 + sizeof(a2)/sizeof(a2[0]));

    map = square.getCodim0MappingIndex(globalNodeIndexes7,globalNodeIndexes8 );
    cout << "square.getCodim0MappingIndex() = " << map << endl;

    /* Squares:
     *
     * 2-------3 0-------1
     * |       | |       |
     * |       | |       |
     * |       | |       |
     * 0-------1 2-------3
     *
     * Should give y=0 reflection.
     *
     */
    a1[0] = 0; a1[1] = 1; a1[2] = 2; a1[3] = 3;
    a2[0] = 2; a2[1] = 3; a2[2] = 0; a2[3] = 1;
    const std::vector<unsigned int> globalNodeIndexes9(a1, a1 + sizeof(a1)/sizeof(a1[0]));
    const std::vector<unsigned int> globalNodeIndexes10(a2, a2 + sizeof(a2)/sizeof(a2[0]));

    map = square.getCodim0MappingIndex(globalNodeIndexes9,globalNodeIndexes10);
    cout << "square.getCodim0MappingIndex() = " << map << endl;

    /* Squares:
     *
     * 2-------3 0-------2
     * |       | |       |
     * |       | |       |
     * |       | |       |
     * 0-------1 1-------3
     *
     * Should give 270¼ rotation (clockwise).
     *
     */
    a1[0] = 0; a1[1] = 1; a1[2] = 2; a1[3] = 3;
    a2[0] = 1; a2[1] = 3; a2[2] = 0; a2[3] = 2;
    const std::vector<unsigned int> globalNodeIndexes11(a1, a1 + sizeof(a1)/sizeof(a1[0]));
    const std::vector<unsigned int> globalNodeIndexes12(a2, a2 + sizeof(a2)/sizeof(a2[0]));

    map = square.getCodim0MappingIndex(globalNodeIndexes11,globalNodeIndexes12);
    cout << "square.getCodim0MappingIndex() = " << map << endl;

    /* Squares:
     *
     * 2-------3 2-------0
     * |       | |       |
     * |       | |       |
     * |       | |       |
     * 0-------1 3-------1
     *
     * Should give a y=-x reflection.
     *
     */
    a1[0] = 0; a1[1] = 1; a1[2] = 2; a1[3] = 3;
    a2[0] = 3; a2[1] = 1; a2[2] = 2; a2[3] = 0;
    const std::vector<unsigned int> globalNodeIndexes13(a1, a1 + sizeof(a1)/sizeof(a1[0]));
    const std::vector<unsigned int> globalNodeIndexes14(a2, a2 + sizeof(a2)/sizeof(a2[0]));

    map = square.getCodim0MappingIndex(globalNodeIndexes13,globalNodeIndexes14);
    cout << "square.getCodim0MappingIndex() = " << map << endl;

    /* Squares:
     *
     * 2-------3 1-------0
     * |       | |       |
     * |       | |       |
     * |       | |       |
     * 0-------1 3-------2
     *
     * Should give a 180¼ rotation.
     *
     */
    a1[0] = 0; a1[1] = 1; a1[2] = 2; a1[3] = 3;
    a2[0] = 3; a2[1] = 2; a2[2] = 1; a2[3] = 0;
    const std::vector<unsigned int> globalNodeIndexes15(a1, a1 + sizeof(a1)/sizeof(a1[0]));
    const std::vector<unsigned int> globalNodeIndexes16(a2, a2 + sizeof(a2)/sizeof(a2[0]));

    map = square.getCodim0MappingIndex(globalNodeIndexes15,globalNodeIndexes16);
    cout << "square.getCodim0MappingIndex() = " << map << endl;


    return 0;
}

