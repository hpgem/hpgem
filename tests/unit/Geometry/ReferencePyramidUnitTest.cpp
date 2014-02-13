//
//  ReferencePyramidUnitTest.cpp
//  
//
//  Created by Shavarsh Nurijanyan on 2/7/13.
//
//

#include "../../../src/Geometry/ReferencePyramid.hpp"
using namespace std;


int main ()
{
    
    cout << "Geometry::ReferencePyramid pyramid;" << endl;

    Geometry::ReferencePyramid pyramid;

    cout << "pyramid.getName() "<< pyramid.getName() << endl;

    cout << "cout << pyramid: ";
    cout << pyramid;

    Geometry::PointReference pCenter(3);
    pyramid.getCenter(pCenter);
    cout << "pyramid.getCenter(Geometry::Point<3>): " << pCenter; cout << endl;;

    Geometry::PointReference p(3);
    pyramid.getNode(0, p);
    cout << "pyramid.getNode(0, p): " << p; cout << endl;

    pyramid.getNode(1, p);
    cout << "pyramid.getNode(1, p): " << p; cout << endl;

    pyramid.getNode(2, p);
    cout << "pyramid.getNode(2, p): " << p; cout << endl;

    pyramid.getNode(3, p);
    cout << "pyramid.getNode(3, p): " << p; cout << endl;

    pyramid.getNode(4, p);
    cout << "pyramid.getNode(4, p): " << p; cout << endl;

    cout << "pyramid.isInternalPoint({0,0}): " << pyramid.isInternalPoint(pCenter) << endl;
    
    Geometry::PointReference pOut(3);

        pOut[0] = -2.0;
        pOut[1] = -2.0;
        pOut[1] = -3.0;

    cout << "pyramid.isInternalPoint({-2.0,-2.0,-2.0}): " << pyramid.isInternalPoint(pOut) << endl;

    cout << "pyramid.getNrOfCodim1Entities(): " << pyramid.getNrOfCodim1Entities() << " faces" << endl;

    cout << "pyramid.getNrOfCodim2Entities(): " << pyramid.getNrOfCodim2Entities() << " edges" << endl;

    cout << "pyramid.getNrOfCodim3Entities(): " << pyramid.getNrOfCodim2Entities() << " nodes" << endl;

    std::vector<unsigned int> indexes;

    for (int face=0; face < pyramid.getNrOfCodim1Entities(); ++face)
    {
        pyramid.getCodim1EntityLocalIndices(face, indexes);
        cout << "pyramid.getCodim1EntityLocalIndices("<<face<<"): (";
        for (int i = 0; i < indexes.size(); ++i)
        {
             if (i < indexes.size()-1) cout << indexes[i] << ",";
             else cout << indexes[i];
        }
        cout << ")" << endl;
    }

    for (int edge=0; edge < pyramid.getNrOfCodim2Entities(); ++edge)
    {
        pyramid.getCodim2EntityLocalIndices(edge, indexes);
        cout << "pyramid.getCodim2EntityLocalIndices("<<edge<<"): (";
        for (int i = 0; i < indexes.size(); ++i)
        {
             if (i < indexes.size()-1) cout << indexes[i] << ",";
             else cout << indexes[i];
        }
        cout << ")" << endl;
    }
    
    const Geometry::ReferenceGeometry* referenceGeometry2D;
    for (int face=0; face < pyramid.getNrOfCodim1Entities(); ++face)
    {
        referenceGeometry2D = pyramid.getCodim1ReferenceGeometry(face);
        cout << "pyramid.getCodim1ReferenceGeometry("<<face<<"): " << referenceGeometry2D->getName() << endl;
    }

    const Geometry::ReferenceGeometry* referenceGeometry1D;
    for (int edge=0; edge < pyramid.getNrOfCodim2Entities(); ++edge)
    {
        referenceGeometry1D = pyramid.getCodim2ReferenceGeometry(edge);
        cout << "pyramid.getCodim1ReferenceGeometry("<<edge<<"): " << referenceGeometry1D->getName() << endl;
    }

    /*
    // There are no pyramid to pyramid mappings.
    const Geometry::MappingReferenceToReference<3,3>* ref3ToRef3Mapping;
    for (int index=0; index < 8; ++index)
    {
        ref3ToRef3Mapping = pyramid.getCodim0MappingPtr(index);
    }
    */

    const Geometry::MappingReferenceToReference* ref2ToRef3Mapping;
    for (int face=0; face < pyramid.getNrOfCodim1Entities(); ++face)
    {
        ref2ToRef3Mapping = pyramid.getCodim1MappingPtr(face);
    }

    /*
    // There are no line to pyramid mappings.
    const Geometry::MappingReferenceToReference<1,3>* ref1ToRef3Mapping;
    for (int face=0; face < pyramid.getNrOfCodim2Entities(); ++face)
    {
        ref1ToRef3Mapping = pyramid.getCodim2MappingPtr(face);
    }
    */

    return 0;
}

