#include "Geometry/PointPhysical.hpp"
#include "Geometry/ReferenceTriangle.hpp"
#include "Geometry/PhysicalTriangle.hpp"
#include "Geometry/ReferenceCube.hpp"
#include "Geometry/PhysicalTetrahedron.hpp"
#include "Geometry/Mappings/MappingToPhysSimplexLinear.hpp"

#include <vector>
#include <iostream>
#include <fstream>

int main()
{
    /* General test of simplex mappings.
     *
     * We create or own global vector of nodes for every dimension.
     *
     * The output to a file is done to be read by 'Mathematica',
     * to have a visual representation of the mapping.
     */

    cout << "2D ~~~~~~~~~ 2D ~~~~~~~~~ 2D ~~~~~~~~~ 2D" << std::endl;

    std::ofstream file2D;
    file2D.open ("2DSimplexMapping.txt");

    // The ordering of the vertices is fundamental.
    // See comment on PhysicalTriangle.cpp

    /*
    Geometry::PointPhysical<2> point2D_1((double []){2.0,2.0});
    Geometry::PointPhysical<2> point2D_2((double []){4.0,2.0});
    Geometry::PointPhysical<2> point2D_3((double []){2.0,4.0});
    */

    // Physical geometry creation
    /*double dummy1[2]= {1.5,3.0};
    Geometry::PointPhysical<2> point2D_1(dummy1);
    
    double dummy2[2]= {4.0,3.9};
    Geometry::PointPhysical<2> point2D_2(dummy2);
    
    double dummy3[2]= {2.0,6.0};
    Geometry::PointPhysical<2> point2D_3(dummy3);

    std::vector<Geometry::PointPhysical<2> > points2D;
    points2D.push_back(point2D_1);
    points2D.push_back(point2D_2);
    points2D.push_back(point2D_3);

    file2D << point2D_1;
    file2D << point2D_2;
    file2D << point2D_3;

    unsigned int indexes2D[3] = {0,1,2};
    const std::vector<unsigned int> globalIndexes2D(indexes2D,indexes2D + sizeof(indexes2D)/sizeof(indexes2D[0]));

    Geometry::ReferenceTriangle referenceTriangle;
    Geometry::PhysicalTriangle physicalTriangle(globalIndexes2D,points2D,&referenceTriangle);

    const Geometry::PhysicalTriangle*const lala = &physicalTriangle;

    cout << "Geometry::MappingSimplexLinear<2> mapping(PhysicalTriangle)" << std::endl;
    Geometry::MappingToPhysSimplexLinear<2> mapping2D(lala);

    cout << "mapping.transform(PointReference,PointPhysical)" << std::endl;
    Geometry::PointReference<2> pointReference2D_1;
    Geometry::PointPhysical<2> pointPhysical2D_1;

    pointReference2D_1[0] =  0.0;
    pointReference2D_1[1] =  0.0;
    cout << "pointReference: " << pointReference2D_1; cout << std::endl;
    mapping2D.transform(pointReference2D_1,pointPhysical2D_1);
    cout << "pointPhysical: " << pointPhysical2D_1; cout << std::endl;

    file2D << pointReference2D_1;
    file2D << pointPhysical2D_1;

    pointReference2D_1[0] = +1.0;
    pointReference2D_1[1] =  0.0;
    cout << "pointReference: " << pointReference2D_1; cout << std::endl;;
    mapping2D.transform(pointReference2D_1,pointPhysical2D_1);
    cout << "pointPhysical: " << pointPhysical2D_1; cout << std::endl;;

    file2D << pointReference2D_1;
    file2D << pointPhysical2D_1;

    pointReference2D_1[0] =  0.0;
    pointReference2D_1[1] =  1.0;
    cout << "pointReference: " << pointReference2D_1; cout << std::endl;
    mapping2D.transform(pointReference2D_1,pointPhysical2D_1);
    cout << "pointPhysical: " << pointPhysical2D_1; cout << std::endl;

    file2D << pointReference2D_1;
    file2D << pointPhysical2D_1;

    pointReference2D_1[0] =  0.5/sqrt(2);
    pointReference2D_1[1] =  0.5/sqrt(2);
    cout << "pointReference: " << pointReference2D_1; cout << std::endl;;
    mapping2D.transform(pointReference2D_1,pointPhysical2D_1);
    cout << "pointPhysical: " << pointPhysical2D_1; cout << std::endl;;

    file2D << pointReference2D_1;
    file2D << pointPhysical2D_1;

    file2D.close();

    cout << "3D ~~~~~~~~~ 3D ~~~~~~~~~ 3D ~~~~~~~~~ 3D" << std::endl;

    // The ordering of the vertices is fundamental.
    // See comment on PhysicalCube.cpp

    std::ofstream file3D;
    file3D.open ("3DSimplexMapping.txt");

    double dummy4[3]= {2.0,2.0,2.0};
    Geometry::PointPhysical<3> point3D_1(dummy4);

    double dummy5[3]= {4.0,2.0,2.0};
    Geometry::PointPhysical<3> point3D_2(dummy5);

    double dummy6[3]= {2.0,4.0,2.0};
    Geometry::PointPhysical<3> point3D_3(dummy6);

    double dummy7[3]= {2.0,2.0,4.0};
    Geometry::PointPhysical<3> point3D_4(dummy7);

    std::vector<Geometry::PointPhysical<3> > points3D;
    points3D.push_back(point3D_1);
    points3D.push_back(point3D_2);
    points3D.push_back(point3D_3);
    points3D.push_back(point3D_4);

    file3D << point3D_1;
    file3D << point3D_2;
    file3D << point3D_4;
    file3D << point3D_3;

    unsigned int indexes3D[4] = {0,1,2,3};
    const std::vector<unsigned int> globalIndexes3D(indexes3D,indexes3D + sizeof(indexes3D)/sizeof(indexes3D[0]));

    Geometry::ReferenceCube referenceCube;
    Geometry::PhysicalTetrahedron physicalTetrahedron(globalIndexes3D,points3D,&referenceCube);

    const Geometry::PhysicalTetrahedron*const lele = &physicalTetrahedron;

    cout << "Geometry::MappingSimplexLinear<3> mapping(physicalTetrahedron)" << std::endl;
    Geometry::MappingToPhysSimplexLinear<3> mapping3D(lele);

    cout << "mapping.transform(PointReference,PointPhysical)" << std::endl;
    Geometry::PointReference<3> pointReference3D_1;
    Geometry::PointPhysical<3> pointPhysical3D_1;

    pointReference3D_1[0] = 0.0;
    pointReference3D_1[1] = 0.0;
    pointReference3D_1[2] = 0.0;
    cout << "pointReference: " << pointReference3D_1; cout << std::endl;
    mapping3D.transform(pointReference3D_1,pointPhysical3D_1);
    cout << "pointPhysical: " << pointPhysical3D_1; cout << std::endl;;
    file3D << pointReference3D_1;
    file3D << pointPhysical3D_1;

    pointReference3D_1[0] = +1.0;
    pointReference3D_1[1] =  0.0;
    pointReference3D_1[2] =  0.0;
    cout << "pointReference: " << pointReference3D_1; cout << std::endl;;
    mapping3D.transform(pointReference3D_1,pointPhysical3D_1);
    cout << "pointPhysical: " << pointPhysical3D_1; cout << std::endl;;
    file3D << pointReference3D_1;
    file3D << pointPhysical3D_1;

    pointReference3D_1[0] =  0.0;
    pointReference3D_1[1] = +1.0;
    pointReference3D_1[2] =  0.0;
    cout << "pointReference: " << pointReference3D_1; cout << std::endl;;
    mapping3D.transform(pointReference3D_1,pointPhysical3D_1);
    cout << "pointPhysical: " << pointPhysical3D_1; cout << std::endl;;
    file3D << pointReference3D_1;
    file3D << pointPhysical3D_1;

    pointReference3D_1[0] =  0.0;
    pointReference3D_1[1] =  0.0;
    pointReference3D_1[2] = +1.0;
    cout << "pointReference: " << pointReference3D_1; cout << std::endl;;
    mapping3D.transform(pointReference3D_1,pointPhysical3D_1);
    cout << "pointPhysical: " << pointPhysical3D_1; cout << std::endl;;
    file3D << pointReference3D_1;
    file3D << pointPhysical3D_1;

    return 0;*/
}
