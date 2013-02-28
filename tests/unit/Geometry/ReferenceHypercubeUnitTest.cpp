//
//  ReferenceHypercubeUnitTest.cpp
//  
//
//  Created by Shavarsh Nurijanyan on 2/7/13.
//
//
#include "Geometry/ReferenceHypercube.hpp"
using namespace std;

using namespace std;


int main ()
{
    cout << "Geometry::ReferenceHypercube hypercube;" << endl;
    Geometry::ReferenceHypercube hypercube;
    
    cout << "hypercube.getName() " << hypercube.getName() << endl;

    cout << "cout << hypercube: ";
    cout << hypercube;

    Geometry::PointReference<4> pCenter;
    hypercube.getCenter(pCenter);
    cout << "hypercube.getCenter(Geometry::ReferencePoint<4>): " << pCenter; cout << endl;
    
    Geometry::PointReference<4> p;
    for (int i = 0; i < 16; ++i)
    {
        hypercube.getNode(i, p);
        cout << "hypercube.getNode("<<i<<", p): " << p; cout << endl;
    }
    
    cout << "hypercube.isInternalPoint({0,0,0,0}): " << hypercube.isInternalPoint(pCenter) << endl;
    
    Geometry::PointReference<4> pOut;
    
    pOut[0] = -2.0;
    pOut[1] = -2.0;
    pOut[2] = -2.0;
    pOut[3] = -2.0;
    
    cout << "hypercube.isInternalPoint({-2.0,-2.0,-2.0,-2.0}): " << hypercube.isInternalPoint(pOut) << endl;

    pOut[0] = -1.0;
    pOut[1] = -1.0;
    pOut[2] = -1.0;
    pOut[3] = -1.0;

    cout << "hypercube.isInternalPoint({-1.0,-1.0,-1.0,-1.0}): " << hypercube.isInternalPoint(pOut) << endl;

    pOut[0] = 0.0;
    pOut[1] = 0.0;
    pOut[2] = 0.0;
    pOut[3] = 0.0;

    cout << "hypercube.getNrOfCodim1Entities(): " << hypercube.getNrOfCodim1Entities() << " cubes" << endl;

    cout << "hypercube.getNrOfCodim2Entities(): " << hypercube.getNrOfCodim2Entities() << " faces" << endl;

    cout << "hypercube.getNrOfCodim3Entities(): " << hypercube.getNrOfCodim3Entities() << " edges" << endl;
    
    std::vector<unsigned int> indexes;

    for (int face=0; face < hypercube.getNrOfCodim1Entities(); ++face)
    {
        hypercube.getCodim1EntityLocalIndices(face, indexes);
        cout << "hypercube.getCodim1EntityLocalIndices("<<face<<"): (";
        for (int i = 0; i < indexes.size(); ++i)
        {
             if (i < indexes.size()-1) cout << indexes[i] << ",";
             else cout << indexes[i];
        }
        cout << ")" << endl;
    }
    /*
    for (int edge=0; edge < hypercube.getNrOfCodim2Entities(); ++edge)
    {
        hypercube.getCodim2EntityLocalIndices(edge, indexes);
        cout << "hypercube.getCodim2EntityLocalIndices("<<edge<<"): (";
        for (int i = 0; i < indexes.size(); ++i)
        {
             if (i < indexes.size()-1) cout << indexes[i] << ",";
             else cout << indexes[i];
        }
        cout << ")" << endl;
    }
    */
    /*
    for (int edge=0; edge < hypercube.getNrOfCodim3Entities(); ++edge)
    {
        hypercube.getCodim3EntityLocalIndices(edge, indexes);
        cout << "hypercube.getCodim3EntityLocalIndices("<<edge<<"): (";
        for (int i = 0; i < indexes.size(); ++i)
        {
             if (i < indexes.size()-1) cout << indexes[i] << ",";
             else cout << indexes[i];
        }
        cout << ")" << endl;
    }
    */
    const Geometry::ReferenceGeometry<3>* referenceGeometry3D;
    for (int face=0; face < hypercube.getNrOfCodim1Entities(); ++face)
    {
        referenceGeometry3D = hypercube.getCodim1ReferenceGeometry(face);
        cout << "cube.getCodim1ReferenceGeometry("<<face<<"): " << referenceGeometry3D->getName() << endl;
    }

    /*
    const Geometry::ReferenceGeometry<2>* referenceGeometry2D;
    for (int face=0; face < hypercube.getNrOfCodim2Entities(); ++face)
    {
        referenceGeometry2D = hypercube.getCodim2ReferenceGeometry(face);
        cout << "cube.getCodim1ReferenceGeometry("<<face<<"): " << referenceGeometry2D->getName() << endl;
    }
    */

    /*
    const Geometry::ReferenceGeometry<1>* referenceGeometry1D;
    for (int face=0; face < hypercube.getNrOfCodim2Entities(); ++face)
    {
        referenceGeometry1D = hypercube.getCodim2ReferenceGeometry(face);
        cout << "cube.getCodim1ReferenceGeometry("<<face<<"): " << referenceGeometry1D->getName() << endl;
    }
    */

    /*
    const Geometry::MappingReferenceToReference<3,3>* ref3ToRef3Mapping;
    for (int index=0; index < 8; ++index)
    {
        ref3ToRef3Mapping = cube.getCodim0MappingPtr(index);
    }
    */

    const Geometry::MappingReferenceToReference<3,4>* ref3ToRef4Mapping;
    for (int face=0; face < hypercube.getNrOfCodim1Entities(); ++face)
    {
        ref3ToRef4Mapping = hypercube.getCodim1MappingPtr(face);
    }

    /*
    const Geometry::MappingReferenceToReference<1,3>* ref1ToRef3Mapping;
    for (int face=0; face < cube.getNrOfCodim1Entities(); ++face)
    {
        ref1ToRef3Mapping = cube.getCodim2MappingPtr(face);
    }
    */

    return 0;
}
