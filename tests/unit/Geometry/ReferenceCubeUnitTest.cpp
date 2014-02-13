//
//  ReferenceCubeUnitTest.cpp
//  
//
//  Created by Shavarsh Nurijanyan on 2/16/13.
//
//
#include "Geometry/ReferencePoint.hpp"
#include "Geometry/ReferenceLine.hpp"
#include "Geometry/ReferenceSquare.hpp"
#include "Geometry/ReferenceCube.hpp"
#include "Geometry/Mappings/MappingToRefCubeToCube.hpp"
using namespace std;

int main ()
{
    cout << "Geometry::ReferenceCube cube;" << endl;
    Geometry::ReferenceCube cube;
    
    cout << "cube.getName() " << cube.getName() << endl;

    cout << "cout << cube: ";
    cout << cube;

    Geometry::PointReference pCenter(3);
    cube.getCenter(pCenter);
    cout << "cube.getCenter(Geometry::ReferencePoint<3>): " << pCenter << endl;
    
    Geometry::PointReference p(3);
    for (int i = 0; i < 8; i++)
    {
        cube.getNode(i, p);
        cout << "cube.getNode("<<i<<", p): " << p << endl;
    }

    cout << "cube.isInternalPoint({0,0,0}): " << cube.isInternalPoint(pCenter) << endl;
    
    Geometry::PointReference pOut(3);
    
    pOut[0] = -2.0;
    pOut[1] = -2.0;
    pOut[2] = -2.0;
    
    cout << "cube.isInternalPoint({-2.0,-2.0,-2.0}): " << cube.isInternalPoint(pOut) << endl;

    pOut[0] = -1.0;
    pOut[1] = -1.0;
    pOut[2] = -1.0;

    cout << "cube.isInternalPoint({-1.0,-1.0,-1.0}): " << cube.isInternalPoint(pOut) << endl;

    pOut[0] = 0.0;
    pOut[1] = 0.0;
    pOut[2] = 0.0;

    cout << "cube.isInternalPoint({0.0,0.0,0.0}): " << cube.isInternalPoint(pOut) << endl;

    cout << "cube.getNrOfCodim1Entities(): " << cube.getNrOfCodim1Entities() << " faces" << endl;

    cout << "cube.getNrOfCodim2Entities(): " << cube.getNrOfCodim2Entities() << " edges" << endl;

    cout << "cube.getNrOfCodim3Entities(): " << cube.getNrOfCodim3Entities() << " points" << endl;
    
    std::vector<unsigned int> indexes;

    for (int face=0; face < cube.getNrOfCodim1Entities(); ++face)
    {
        cube.getCodim1EntityLocalIndices(face, indexes);
        cout << "cube.getCodim1EntityLocalIndices("<<face<<"): (";
        for (int i = 0; i < indexes.size(); ++i)
        {
             if (i < indexes.size()-1) cout << indexes[i] << ",";
             else cout << indexes[i];
        }
        cout << ")" << endl;
    }

    for (int edge=0; edge < cube.getNrOfCodim2Entities(); ++edge)
    {
        cube.getCodim2EntityLocalIndices(edge, indexes);
        cout << "cube.getCodim2EntityLocalIndices("<<edge<<"): (";
        for (int i = 0; i < indexes.size(); ++i)
        {
             if (i < indexes.size()-1) cout << indexes[i] << ",";
             else cout << indexes[i];
        }
        cout << ")" << endl;
    }

    // ------------------------------------------------------------------------------------------ \\

    const Geometry::ReferenceGeometry* referenceGeometry2D;
    for (int face=0; face < cube.getNrOfCodim1Entities(); ++face)
    {
        referenceGeometry2D = cube.getCodim1ReferenceGeometry(face);
        cout << "cube.getCodim1ReferenceGeometry("<<face<<"): " << referenceGeometry2D->getName() << endl;
        //cout << "cube.getCodim1ReferenceGeometry("<<face<<"): " << *((Geometry::ReferenceSquare*)(referenceGeometry2D));
    }

    const Geometry::ReferenceGeometry* referenceGeometry1D;
    for (int face=0; face < cube.getNrOfCodim2Entities(); ++face)
    {
        referenceGeometry1D = cube.getCodim2ReferenceGeometry(face);
        cout << "cube.getCodim1ReferenceGeometry("<<face<<"): " << referenceGeometry1D->getName() << endl;
        //cout << "cube.getCodim1ReferenceGeometry("<<face<<"): " << *((Geometry::ReferenceLine*)(referenceGeometry1D));
    }

    // ------------------------------------------------------------------------------------------ \\

    const Geometry::MappingReferenceToReference* ref3ToRef3Mapping;
    for (int index=0; index < 8; ++index)
    {
        ref3ToRef3Mapping = cube.getCodim0MappingPtr(index);
    }

    const Geometry::MappingReferenceToReference* ref2ToRef3Mapping;
    for (int face=0; face < cube.getNrOfCodim1Entities(); ++face)
    {
        ref2ToRef3Mapping = cube.getCodim1MappingPtr(face);
    }

    const Geometry::MappingReferenceToReference* ref1ToRef3Mapping;
    for (int face=0; face < cube.getNrOfCodim1Entities(); ++face)
    {
        ref1ToRef3Mapping = cube.getCodim2MappingPtr(face);
    }

    // ------------------------------------------------------------------------------------------ \\

    int a1[8]; int a2[8]; int map;

    // In the following we use the square drawings... see the comment in MappingCubeToCube.
    //
    // Squares:
    //
    // 2-------3 2-------3
    // |       | |       |
    // |       | |       |
    // |       | |       |
    // 0-------1 0-------1
    //
    // Should give identity.

    a1[0] = 0; a1[1] = 1; a1[2] = 2; a1[3] = 3; a1[4] = 4; a1[5] = 5; a1[6] = 6; a1[7] = 7;
    a2[0] = 0; a2[1] = 1; a2[2] = 2; a2[3] = 3; a2[4] = 4; a2[5] = 5; a2[6] = 6; a2[7] = 7;
    const std::vector<unsigned int> globalNodeIndexes1(a1, a1 + sizeof(a1)/sizeof(a1[0]));
    const std::vector<unsigned int> globalNodeIndexes2(a2, a2 + sizeof(a2)/sizeof(a2[0]));

    map = cube.getCodim0MappingIndex(globalNodeIndexes1,globalNodeIndexes2);
    cout << "cube.getCodim0MappingIndex() = " << map << endl;


    // Squares:
    //
    // 2-------3 1-------3
    // |       | |       |
    // |       | |       |
    // |       | |       |
    // 0-------1 0-------2
    //
    // Should give x=y reflection.

    a1[0] = 0; a1[1] = 1; a1[2] = 2; a1[3] = 3; a1[4] = 4; a1[5] = 5; a1[6] = 6; a1[7] = 7;
    a2[0] = 0; a2[1] = 2; a2[2] = 1; a2[3] = 3; a2[4] = 4; a2[5] = 6; a2[6] = 5; a2[7] = 7;
    const std::vector<unsigned int> globalNodeIndexes3(a1, a1 + sizeof(a1)/sizeof(a1[0]));
    const std::vector<unsigned int> globalNodeIndexes4(a2, a2 + sizeof(a2)/sizeof(a2[0]));

    map = cube.getCodim0MappingIndex(globalNodeIndexes3,globalNodeIndexes4);
    cout << "cube.getCodim0MappingIndex() = " << map << endl;

    // Squares:
    //
    // 2-------3 3-------2
    // |       | |       |
    // |       | |       |
    // |       | |       |
    // 0-------1 1-------0
    //
    // Should give x=0 reflection.

    a1[0] = 0; a1[1] = 1; a1[2] = 2; a1[3] = 3; a1[4] = 4; a1[5] = 5; a1[6] = 6; a1[7] = 7;
    a2[0] = 1; a2[1] = 0; a2[2] = 3; a2[3] = 2; a2[4] = 5; a2[5] = 4; a2[6] = 7; a2[7] = 6;
    const std::vector<unsigned int> globalNodeIndexes5(a1, a1 + sizeof(a1)/sizeof(a1[0]));
    const std::vector<unsigned int> globalNodeIndexes6(a2, a2 + sizeof(a2)/sizeof(a2[0]));

    map = cube.getCodim0MappingIndex(globalNodeIndexes5,globalNodeIndexes6);
    cout << "cube.getCodim0MappingIndex() = " << map << endl;

    // Squares:
    //
    // 2-------3 3-------1
    // |       | |       |
    // |       | |       |
    // |       | |       |
    // 0-------1 2-------0
    //
    // Should give 90¼ rotation (clockwise).

    a1[0] = 0; a1[1] = 1; a1[2] = 2; a1[3] = 3; a1[4] = 4; a1[5] = 5; a1[6] = 6; a1[7] = 7;
    a2[0] = 2; a2[1] = 0; a2[2] = 3; a2[3] = 1; a2[4] = 6; a2[5] = 4; a2[6] = 7; a2[7] = 5;
    const std::vector<unsigned int> globalNodeIndexes7(a1, a1 + sizeof(a1)/sizeof(a1[0]));
    const std::vector<unsigned int> globalNodeIndexes8(a2, a2 + sizeof(a2)/sizeof(a2[0]));

    map = cube.getCodim0MappingIndex(globalNodeIndexes7,globalNodeIndexes8 );
    cout << "cube.getCodim0MappingIndex() = " << map << endl;

    // Squares:
    //
    // 2-------3 0-------1
    // |       | |       |
    // |       | |       |
    // |       | |       |
    // 0-------1 2-------3
    //
    // Should give y=0 reflection.

    a1[0] = 0; a1[1] = 1; a1[2] = 2; a1[3] = 3; a1[4] = 4; a1[5] = 5; a1[6] = 6; a1[7] = 7;
    a2[0] = 2; a2[1] = 3; a2[2] = 0; a2[3] = 1; a2[4] = 6; a2[5] = 7; a2[6] = 4; a2[7] = 5;
    const std::vector<unsigned int> globalNodeIndexes9(a1, a1 + sizeof(a1)/sizeof(a1[0]));
    const std::vector<unsigned int> globalNodeIndexes10(a2, a2 + sizeof(a2)/sizeof(a2[0]));

    map = cube.getCodim0MappingIndex(globalNodeIndexes9,globalNodeIndexes10);
    cout << "cube.getCodim0MappingIndex() = " << map << endl;

    // Squares:
    //
    // 2-------3 0-------2
    // |       | |       |
    // |       | |       |
    // |       | |       |
    // 0-------1 1-------3
    //
    // Should give 270¼ rotation (clockwise).

    a1[0] = 0; a1[1] = 1; a1[2] = 2; a1[3] = 3; a1[4] = 4; a1[5] = 5; a1[6] = 6; a1[7] = 7;
    a2[0] = 1; a2[1] = 3; a2[2] = 0; a2[3] = 2; a2[4] = 5; a2[5] = 7; a2[6] = 4; a2[7] = 6;
    const std::vector<unsigned int> globalNodeIndexes11(a1, a1 + sizeof(a1)/sizeof(a1[0]));
    const std::vector<unsigned int> globalNodeIndexes12(a2, a2 + sizeof(a2)/sizeof(a2[0]));

    map = cube.getCodim0MappingIndex(globalNodeIndexes11,globalNodeIndexes12);
    cout << "cube.getCodim0MappingIndex() = " << map << endl;

    // Squares:
    //
    // 2-------3 2-------0
    // |       | |       |
    // |       | |       |
    // |       | |       |
    // 0-------1 3-------1
    //
    // Should give a y=-x reflection.

    a1[0] = 0; a1[1] = 1; a1[2] = 2; a1[3] = 3; a1[4] = 4; a1[5] = 5; a1[6] = 6; a1[7] = 7;
    a2[0] = 3; a2[1] = 1; a2[2] = 2; a2[3] = 0; a2[4] = 7; a2[5] = 5; a2[6] = 6; a2[7] = 4;
    const std::vector<unsigned int> globalNodeIndexes13(a1, a1 + sizeof(a1)/sizeof(a1[0]));
    const std::vector<unsigned int> globalNodeIndexes14(a2, a2 + sizeof(a2)/sizeof(a2[0]));

    map = cube.getCodim0MappingIndex(globalNodeIndexes13,globalNodeIndexes14);
    cout << "cube.getCodim0MappingIndex() = " << map << endl;

    // Squares:
    //
    // 2-------3 1-------0
    // |       | |       |
    // |       | |       |
    // |       | |       |
    // 0-------1 3-------2
    //
    // Should give a 180¼ rotation.

    a1[0] = 0; a1[1] = 1; a1[2] = 2; a1[3] = 3; a1[4] = 4; a1[5] = 5; a1[6] = 6; a1[7] = 7;
    a2[0] = 3; a2[1] = 2; a2[2] = 1; a2[3] = 0; a2[4] = 7; a2[5] = 6; a2[6] = 5; a2[7] = 4;
    const std::vector<unsigned int> globalNodeIndexes15(a1, a1 + sizeof(a1)/sizeof(a1[0]));
    const std::vector<unsigned int> globalNodeIndexes16(a2, a2 + sizeof(a2)/sizeof(a2[0]));

    map = cube.getCodim0MappingIndex(globalNodeIndexes15,globalNodeIndexes16);
    cout << "cube.getCodim0MappingIndex() = " << map << endl;

    return 0;
}
