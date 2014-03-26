//
//  ReferenceTetrahedronUnitTest.cpp
//  
//
//  Created by Shavarsh Nurijanyan on 2/7/13.
//
//

#include "../../../kernel/Geometry/ReferenceTetrahedron.hpp"
using namespace std;


int main ()
{
    
    cout << "Geometry::ReferenceTetrahedron tetrahedron;" << endl;

    Geometry::ReferenceTetrahedron tetrahedron;

    cout << "square.getName() "<< tetrahedron.getName() << endl;

    cout << "cout << tetrahedron: ";
    cout << tetrahedron;

    Geometry::PointReference pCenter(3);
    tetrahedron.getCenter(pCenter);
    cout << "tetrahedron.getCenter(Geometry::Point<3>): " << pCenter;

    Geometry::PointReference p(3);
    tetrahedron.getNode(0, p);
    cout << "tetrahedron.getNode(0, p): " << p;

    tetrahedron.getNode(1, p);
    cout << "tetrahedron.getNode(1, p): " << p;

    tetrahedron.getNode(2, p);
    cout << "tetrahedron.getNode(2, p): " << p;

    tetrahedron.getNode(3, p);
    cout << "tetrahedron.getNode(3, p): " << p;

    cout << "tetrahedron.isInternalPoint({0,0}): " << tetrahedron.isInternalPoint(pCenter) << endl;
    
    Geometry::PointReference pOut(3);

        pOut[0] = -2.0;
        pOut[1] = -2.0;
        pOut[1] = -3.0;

    cout << "tetrahedron.isInternalPoint({-2.0,-2.0,-2.0}): " << tetrahedron.isInternalPoint(pOut) << endl;

    cout << "tetrahedron.getNrOfCodim1Entities(): " << tetrahedron.getNrOfCodim1Entities() << " lines" << endl;

    cout << "tetrahedron.getNrOfCodim2Entities(): " << tetrahedron.getNrOfCodim2Entities() << " points" << endl;

    std::vector<unsigned int> indexes;

    for (int face=0; face < tetrahedron.getNrOfCodim1Entities(); ++face)
    {
        tetrahedron.getCodim1EntityLocalIndices(face, indexes);
        cout << "tetrahedron.getCodim1EntityLocalIndices("<<face<<"): (";
        for (int i = 0; i < indexes.size(); ++i)
        {
             if (i < indexes.size()-1) cout << indexes[i] << ",";
             else cout << indexes[i];
        }
        cout << ")" << endl;
    }

    for (int edge=0; edge < tetrahedron.getNrOfCodim2Entities(); ++edge)
    {
        tetrahedron.getCodim2EntityLocalIndices(edge, indexes);
        cout << "tetrahedron.getCodim2EntityLocalIndices("<<edge<<"): (";
        for (int i = 0; i < indexes.size(); ++i)
        {
             if (i < indexes.size()-1) cout << indexes[i] << ",";
             else cout << indexes[i];
        }
        cout << ")" << endl;
    }
    
    const Geometry::ReferenceGeometry* referenceGeometry2D;
    for (int face=0; face < tetrahedron.getNrOfCodim1Entities(); ++face)
    {
        referenceGeometry2D = tetrahedron.getCodim1ReferenceGeometry(face);
        cout << "tetrahedron.getCodim1ReferenceGeometry("<<face<<"): " << referenceGeometry2D->getName() << endl;
    }

    const Geometry::ReferenceGeometry* referenceGeometry1D;
    for (int edge=0; edge < tetrahedron.getNrOfCodim2Entities(); ++edge)
    {
        referenceGeometry1D = tetrahedron.getCodim2ReferenceGeometry(edge);
        cout << "tetrahedron.getCodim1ReferenceGeometry("<<edge<<"): " << referenceGeometry1D->getName() << endl;
    }

    /*
    // There are no tetrahedron to tetrahedron mappings.
    const Geometry::MappingReferenceToReference<3,3>* ref3ToRef3Mapping;
    for (int index=0; index < 8; ++index)
    {
        ref3ToRef3Mapping = tetrahedron.getCodim0MappingPtr(index);
    }
    */

    const Geometry::MappingReferenceToReference* ref2ToRef3Mapping;
    for (int face=0; face < tetrahedron.getNrOfCodim1Entities(); ++face)
    {
        ref2ToRef3Mapping = tetrahedron.getCodim1MappingPtr(face);
    }

    /*
    // There are no line to tetrahedron mappings.
    const Geometry::MappingReferenceToReference<1,3>* ref1ToRef3Mapping;
    for (int face=0; face < tetrahedron.getNrOfCodim2Entities(); ++face)
    {
        ref1ToRef3Mapping = tetrahedron.getCodim2MappingPtr(face);
    }
    */

    return 0;
}

