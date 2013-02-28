//
//  ReferenceTriangularPrismUnitTest.cpp
//  
//
//  Created by Shavarsh Nurijanyan on 2/7/13.
//
//

#include "../../../src/Geometry/ReferenceTriangularPrism.hpp"
using namespace std;


int main ()
{
    
    cout << "Geometry::ReferenceTriangularPrism triangularprism;" << endl;

    Geometry::ReferenceTriangularPrism triangularprism;

    cout << "triangularprism.getName() "<< triangularprism.getName() << endl;

    cout << "cout << triangularprism: ";
    cout << triangularprism;

    Geometry::PointReference<3> pCenter;
    triangularprism.getCenter(pCenter);
    cout << "triangularprism.getCenter(Geometry::Point<3>): " << pCenter; cout << endl;

    Geometry::PointReference<3> p;
    triangularprism.getNode(0, p);
    cout << "triangularprism.getNode(0, p): " << p; cout << endl;

    triangularprism.getNode(1, p);
    cout << "triangularprism.getNode(1, p): " << p; cout << endl;

    triangularprism.getNode(2, p);
    cout << "triangularprism.getNode(2, p): " << p; cout << endl;

    triangularprism.getNode(3, p);
    cout << "triangularprism.getNode(3, p): " << p; cout << endl;

    triangularprism.getNode(4, p);
    cout << "triangularprism.getNode(4, p): " << p; cout << endl;

    triangularprism.getNode(5, p);
    cout << "triangularprism.getNode(5, p): " << p; cout << endl;

    cout << "triangularprism.isInternalPoint({0,0}): " << triangularprism.isInternalPoint(pCenter) << endl;
    
    Geometry::PointReference<3> pOut;

        pOut[0] = -2.0;
        pOut[1] = -2.0;
        pOut[2] = -2.0;

    cout << "triangularprism.isInternalPoint({-2.0,-2.0,-2.0}): " << triangularprism.isInternalPoint(pOut) << endl;

    cout << "triangularprism.getNrOfCodim1Entities(): " << triangularprism.getNrOfCodim1Entities() << " faces" << endl;

    cout << "triangularprism.getNrOfCodim2Entities(): " << triangularprism.getNrOfCodim2Entities() << " edges" << endl;

    cout << "triangularprism.getNrOfCodim2Entities(): " << triangularprism.getNrOfCodim2Entities() << " nodes" << endl;

    std::vector<unsigned int> indexes;

    for (int face=0; face < triangularprism.getNrOfCodim1Entities(); ++face)
    {
        triangularprism.getCodim1EntityLocalIndices(face, indexes);
        cout << "triangularprism.getCodim1EntityLocalIndices("<<face<<"): (";
        for (int i = 0; i < indexes.size(); ++i)
        {
             if (i < indexes.size()-1) cout << indexes[i] << ",";
             else cout << indexes[i];
        }
        cout << ")" << endl;
    }

    for (int edge=0; edge < triangularprism.getNrOfCodim2Entities(); ++edge)
    {
        triangularprism.getCodim2EntityLocalIndices(edge, indexes);
        cout << "triangularprism.getCodim2EntityLocalIndices("<<edge<<"): (";
        for (int i = 0; i < indexes.size(); ++i)
        {
             if (i < indexes.size()-1) cout << indexes[i] << ",";
             else cout << indexes[i];
        }
        cout << ")" << endl;
    }

    for (int edge=0; edge < triangularprism.getNrOfCodim3Entities(); ++edge)
    {
        triangularprism.getCodim3EntityLocalIndices(edge, indexes);
        cout << "triangularprism.getCodim3EntityLocalIndices("<<edge<<"): (";
        for (int i = 0; i < indexes.size(); ++i)
        {
             if (i < indexes.size()-1) cout << indexes[i] << ",";
             else cout << indexes[i];
        }
        cout << ")" << endl;
    }
    
    const Geometry::ReferenceGeometry<2>* referenceGeometry2D;
    for (int face=0; face < triangularprism.getNrOfCodim1Entities(); ++face)
    {
        referenceGeometry2D = triangularprism.getCodim1ReferenceGeometry(face);
        cout << "triangularprism.getCodim1ReferenceGeometry("<<face<<"): " << referenceGeometry2D->getName() << endl;
    }

    const Geometry::ReferenceGeometry<1>* referenceGeometry1D;
    for (int edge=0; edge < triangularprism.getNrOfCodim2Entities(); ++edge)
    {
        referenceGeometry1D = triangularprism.getCodim2ReferenceGeometry(edge);
        cout << "triangularprism.getCodim2ReferenceGeometry("<<edge<<"): " << referenceGeometry1D->getName() << endl;
    }

    /*
    // There are no triangularprism to triangularprism mappings.
    const Geometry::MappingReferenceToReference<3,3>* ref3ToRef3Mapping;
    for (int index=0; index < 8; ++index)
    {
        ref3ToRef3Mapping = triangularprism.getCodim0MappingPtr(index);
    }
    */

    const Geometry::MappingReferenceToReference<2,3>* ref2ToRef3Mapping;
    for (int face=0; face < triangularprism.getNrOfCodim1Entities(); ++face)
    {
        ref2ToRef3Mapping = triangularprism.getCodim1MappingPtr(face);
    }

    /*
    // There are no line to triangularprism mappings.
    const Geometry::MappingReferenceToReference<1,3>* ref1ToRef3Mapping;
    for (int face=0; face < triangularprism.getNrOfCodim2Entities(); ++face)
    {
        ref1ToRef3Mapping = triangularprism.getCodim2MappingPtr(face);
    }
    */

    return 0;
}

