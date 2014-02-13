//
//  ReferenceTriangleUnitTest.cpp
//  
//
//  Created by Shavarsh Nurijanyan on 2/5/13.
//
//
#include "Geometry/ReferenceTriangle.hpp"
using namespace std;

int main ()
{

    cout << "Geometry::ReferenceTriangle triangle;" << endl;

    Geometry::ReferenceTriangle triangle;
    
    cout << "triangle.getName() " << triangle.getName() << endl;

    cout << "cout << triangle";
    cout << triangle;

    Geometry::PointReference pCenter(2);
    triangle.getCenter(pCenter);
    cout << "triangle.getCenter(Geometry::Point<2>): " << pCenter << endl;
    
    
    Geometry::PointReference p(2);
    triangle.getNode(0, p);
    cout << "triangle.getVertexPoint(0, p): " << p << endl;
    
    triangle.getNode(1, p);
    cout << "triangle.getVertexPoint(1, p): " << p << endl;

    triangle.getNode(2, p);
    cout << "triangle.getVertexPoint(2, p): " << p << endl;

    cout << "triangle.isInternalPoint({center}): " << triangle.isInternalPoint(pCenter) << endl;

    Geometry::PointReference pOut(2);
    
    pOut[0] = -1.0;
    pOut[1] = -1.0;
    
    cout << "triangle.isInternalPoint({-1.0,-1.0}): " << triangle.isInternalPoint(pOut) << endl;
    
    cout << "square.getNrOfCodim1Entities(): " << triangle.getNrOfCodim1Entities() << " lines" << endl;

    cout << "square.getNrOfCodim2Entities(): " << triangle.getNrOfCodim2Entities() << " points" << endl;
    
    std::vector<unsigned int> indexes;

    for (int face=0; face < triangle.getNrOfCodim1Entities(); ++face)
    {
        triangle.getCodim1EntityLocalIndices(face, indexes);
        cout << "square.getCodim1EntityLocalIndices("<<face<<"): (";
        for (int i = 0; i < indexes.size(); ++i)
        {
             if (i < indexes.size()-1) cout << indexes[i] << ",";
             else cout << indexes[i];
        }
        cout << ")" << endl;
    }

    const Geometry::ReferenceGeometry* referenceGeometry;
    for (int face=0; face < triangle.getNrOfCodim1Entities(); ++face)
    {
        referenceGeometry = triangle.getCodim1ReferenceGeometry(face);
        cout << "triangle.getCodim1ReferenceGeometry("<<face<<"): " << referenceGeometry->getName() << endl;
    }

    const Geometry::MappingReferenceToReference* ref2ToRef2Mapping;
    for (int index=0; index < triangle.getNrOfCodim2Entities(); ++index)
    {
        ref2ToRef2Mapping = triangle.getCodim0MappingPtr(index);
    }

    const Geometry::MappingReferenceToReference* ref1ToRef2Mapping;
    for (int face=0; face < triangle.getNrOfCodim1Entities(); ++face)
    {
        ref1ToRef2Mapping = triangle.getCodim1MappingPtr(face);
    }

    int a1[3]; int a2[3]; int map;
    /* Triangles:
     *
     * 2         2
     * | \       | \
     * |   \     |   \
     * |     \   |     \
     * 0-------1 0-------1
     *
     * Should give identity.
     *
     */
    a1[0] = 0; a1[1] = 1; a1[2] = 2;
    a2[0] = 0; a2[1] = 1; a2[2] = 2;
    const std::vector<unsigned int> globalNodeIndexes1(a1, a1 + sizeof(a1)/sizeof(a1[0]));
    const std::vector<unsigned int> globalNodeIndexes2(a2, a2 + sizeof(a2)/sizeof(a2[0]));

    map = triangle.getCodim0MappingIndex(globalNodeIndexes1,globalNodeIndexes2);
    cout << "square.getCodim0MappingIndex() = " << map << endl;

    /* Triangles:
     *
     * 2         1
     * | \       | \
     * |   \     |   \
     * |     \   |     \
     * 0-------1 0-------2
     *
     */
    a1[0] = 0; a1[1] = 1; a1[2] = 2;
    a2[0] = 0; a2[1] = 2; a2[2] = 1;
    const std::vector<unsigned int> globalNodeIndexes3(a1, a1 + sizeof(a1)/sizeof(a1[0]));
    const std::vector<unsigned int> globalNodeIndexes4(a2, a2 + sizeof(a2)/sizeof(a2[0]));

    map = triangle.getCodim0MappingIndex(globalNodeIndexes3,globalNodeIndexes4);
    cout << "square.getCodim0MappingIndex() = " << map << endl;

    /* Triangles:
     *
     * 2         0
     * | \       | \
     * |   \     |   \
     * |     \   |     \
     * 0-------1 1-------2
     *
     */
    a1[0] = 0; a1[1] = 1; a1[2] = 2;
    a2[0] = 1; a2[1] = 2; a2[2] = 0;
    const std::vector<unsigned int> globalNodeIndexes5(a1, a1 + sizeof(a1)/sizeof(a1[0]));
    const std::vector<unsigned int> globalNodeIndexes6(a2, a2 + sizeof(a2)/sizeof(a2[0]));

    map = triangle.getCodim0MappingIndex(globalNodeIndexes5,globalNodeIndexes6);
    cout << "square.getCodim0MappingIndex() = " << map << endl;

    /* Triangles:
     *
     * 2         2
     * | \       | \
     * |   \     |   \
     * |     \   |     \
     * 0-------1 1-------0
     *
     */
    a1[0] = 0; a1[1] = 1; a1[2] = 2;
    a2[0] = 1; a2[1] = 0; a2[2] = 2;
    const std::vector<unsigned int> globalNodeIndexes7(a1, a1 + sizeof(a1)/sizeof(a1[0]));
    const std::vector<unsigned int> globalNodeIndexes8(a2, a2 + sizeof(a2)/sizeof(a2[0]));

    map = triangle.getCodim0MappingIndex(globalNodeIndexes7,globalNodeIndexes8);
    cout << "square.getCodim0MappingIndex() = " << map << endl;

    /* Triangles:
     *
     * 2         0
     * | \       | \
     * |   \     |   \
     * |     \   |     \
     * 0-------1 2-------1
     *
     */
    a1[0] = 0; a1[1] = 1; a1[2] = 2;
    a2[0] = 2; a2[1] = 1; a2[2] = 0;
    const std::vector<unsigned int> globalNodeIndexes9(a1, a1 + sizeof(a1)/sizeof(a1[0]));
    const std::vector<unsigned int> globalNodeIndexes10(a2, a2 + sizeof(a2)/sizeof(a2[0]));

    map = triangle.getCodim0MappingIndex(globalNodeIndexes9,globalNodeIndexes10);
    cout << "square.getCodim0MappingIndex() = " << map << endl;

    /* Triangles:
     *
     * 2         1
     * | \       | \
     * |   \     |   \
     * |     \   |     \
     * 0-------1 2-------0
     *
     */
    a1[0] = 0; a1[1] = 1; a1[2] = 2;
    a2[0] = 2; a2[1] = 0; a2[2] = 1;
    const std::vector<unsigned int> globalNodeIndexes11(a1, a1 + sizeof(a1)/sizeof(a1[0]));
    const std::vector<unsigned int> globalNodeIndexes12(a2, a2 + sizeof(a2)/sizeof(a2[0]));

    map = triangle.getCodim0MappingIndex(globalNodeIndexes11,globalNodeIndexes12);
    cout << "square.getCodim0MappingIndex() = " << map << endl;

    return 0;
}

