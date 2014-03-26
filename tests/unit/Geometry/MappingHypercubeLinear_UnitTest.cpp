#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

#include "../../../kernel/Geometry/PointPhysical.hpp"
#include "../../../kernel/Geometry/ReferenceLine.hpp"
#include "../../../kernel/Geometry/PhysicalLine.hpp"
#include "../../../kernel/Geometry/ReferenceSquare.hpp"
#include "../../../kernel/Geometry/PhysicalQuadrilateral.hpp"
#include "../../../kernel/Geometry/ReferenceCube.hpp"
#include "../../../kernel/Geometry/PhysicalHexahedron.hpp"
#include "../../../kernel/Geometry/Mappings/MappingToPhysHypercubeLinear.hpp"

int main()
{
    /* General test of cube mappings.
     *
     * We create or own global vector of nodes for every dimension.
     *
     * The output to a file is done to be read by 'Mathematica',
     * to have a visual representation of the mapping.
     *
     * TODO: Add 4D. I'm not sure about the mapping in 4D, as it was not completely implemented
     *       in the old version of hpGEM. See their comments.
     */

    /*cout << "1D ~~~~~~~~~ 1D ~~~~~~~~~ 1D ~~~~~~~~~ 1D ~~~~~~~~~ 1D ~~~~~~~~~ 1D" << std::endl;

    double dummy1[1]= {2.0};
    
    Geometry::PointPhysical point1D_1(dummy1);///\bug PointPhysical does not have an array based initialization!!(see point instead)

    double dummy2[1]= {4.0};

    Geometry::PointPhysical point1D_2(dummy2);
    std::vector<Geometry::PointPhysical > points1D;
    points1D.push_back(point1D_1);
    points1D.push_back(point1D_2);

    unsigned int indexes1D[2] = {0,1};
    const std::vector<unsigned int> globalIndexes1D(indexes1D,indexes1D + sizeof(indexes1D)/sizeof(indexes1D[0]));

    Geometry::ReferenceLine referenceLine;

    Geometry::PhysicalLine physicalLine(globalIndexes1D,points1D,&referenceLine);

    const Geometry::PhysicalLine*const lala1D = &physicalLine;

    cout << "Geometry::MappingSimpleHypercubeLinear<1> mapping(physicalLine)";
    Geometry::MappingToPhysHypercubeLinear<1> mapping1D(lala1D);

    cout << "mapping.transform(PointReference,PointPhysical)" << std::endl;
    Geometry::PointReference pointReference1D_1;
    Geometry::PointPhysical pointPhysical1D_1(1);

    pointReference1D_1[0] = -1.0;
    cout << "pointReference: " << pointReference1D_1; cout << std::endl;
    mapping1D.transform(pointReference1D_1,pointPhysical1D_1);
    cout << "pointPhysical: " << pointPhysical1D_1; cout << std::endl;

    pointReference1D_1[0] = 0.0;
    cout << "pointReference: " << pointReference1D_1; cout << std::endl;
    mapping1D.transform(pointReference1D_1,pointPhysical1D_1);
    cout << "pointPhysical: " << pointPhysical1D_1; cout << std::endl;

    pointReference1D_1[0] = 1.0;
    cout << "pointReference: " << pointReference1D_1; cout << std::endl;
    mapping1D.transform(pointReference1D_1,pointPhysical1D_1);
    cout << "pointPhysical: " << pointPhysical1D_1; cout << std::endl;

    cout << "2D ~~~~~~~~~ 2D ~~~~~~~~~ 2D ~~~~~~~~~ 2D ~~~~~~~~~ 2D ~~~~~~~~~ 2D" << std::endl;

    std::ofstream file2D;
    file2D.open ("2Dmapping.txt");

    // The ordering of the vertices is fundamental.
    // See comment on PhysicalQuadrilateral.cpp

    double dummy3[2]= {2.0,2.0};
    Geometry::PointPhysical point2D_1(dummy3);

    double dummy4[2]= {4.0,2.0};
    Geometry::PointPhysical point2D_2(dummy4);
   
    double dummy5[2]= {2.0,4.0};
    Geometry::PointPhysical point2D_3(dummy5);
    
    double dummy6[2]= {4.0,4.0};
    Geometry::PointPhysical point2D_4(dummy6);

    //Geometry::PointPhysical<2> point2D_1((double []){2.3,2.0});
    //Geometry::PointPhysical<2> point2D_2((double []){4.7,2.8});
    //Geometry::PointPhysical<2> point2D_3((double []){2.1,4.4});
    //Geometry::PointPhysical<2> point2D_4((double []){4.4,4.2});

    std::vector<Geometry::PointPhysical > points2D;
    points2D.push_back(point2D_1);
    points2D.push_back(point2D_2);
    points2D.push_back(point2D_3);
    points2D.push_back(point2D_4);

    file2D << point2D_1;
    file2D << point2D_2;
    file2D << point2D_4; // Different order so that 'Mathematica' makes the correct square.
    file2D << point2D_3;

    unsigned int indexes2D[4] = {0,1,2,3};
    const std::vector<unsigned int> globalIndexes2D(indexes2D,indexes2D + sizeof(indexes2D)/sizeof(indexes2D[0]));

    Geometry::ReferenceSquare referenceSquare;
    Geometry::PhysicalQuadrilateral physicalQuadrilateral(globalIndexes2D,points2D,&referenceSquare);

    const Geometry::PhysicalQuadrilateral*const lala = &physicalQuadrilateral;

    cout << "Geometry::MappingSimpleHypercubeLinear<2> mapping(PhysicalQuadrilateral)" << std::endl;
    Geometry::MappingToPhysHypercubeLinear<2> mapping2D(lala);

    cout << "mapping.transform(PointReference,PointPhysical)" << std::endl;
    Geometry::PointReference pointReference2D_1(2);
    Geometry::PointPhysical pointPhysical2D_1(2);

    pointReference2D_1[0] = -1.0;
    pointReference2D_1[1] = -1.0;
    cout << "pointReference:" << pointReference2D_1; cout << std::endl;
    mapping2D.transform(pointReference2D_1,pointPhysical2D_1);
    cout << "pointPhysical:" << pointPhysical2D_1; cout << std::endl;

    file2D << pointReference2D_1;
    file2D << pointPhysical2D_1;

    pointReference2D_1[0] = -1.0;
    pointReference2D_1[1] = +1.0;
    cout << "pointReference:" << pointReference2D_1; cout << std::endl;
    mapping2D.transform(pointReference2D_1,pointPhysical2D_1);
    cout << "pointPhysical:" << pointPhysical2D_1; cout << std::endl;

    file2D << pointReference2D_1;
    file2D << pointPhysical2D_1;

    pointReference2D_1[0] = +1.0;
    pointReference2D_1[1] = -1.0;
    cout << "pointReference:" << pointReference2D_1; cout << std::endl;
    mapping2D.transform(pointReference2D_1,pointPhysical2D_1);
    cout << "pointPhysical:" << pointPhysical2D_1; cout << std::endl;

    file2D << pointReference2D_1;
    file2D << pointPhysical2D_1;

    pointReference2D_1[0] = +1.0;
    pointReference2D_1[1] = +1.0;
    cout << "pointReference:" << pointReference2D_1; cout << std::endl;
    mapping2D.transform(pointReference2D_1,pointPhysical2D_1);
    cout << "pointPhysical:" << pointPhysical2D_1; cout << std::endl;

    file2D << pointReference2D_1;
    file2D << pointPhysical2D_1;

    pointReference2D_1[0] = 0.0;
    pointReference2D_1[1] = 0.0;
    cout << "pointReference:" << pointReference2D_1; cout << std::endl;
    mapping2D.transform(pointReference2D_1,pointPhysical2D_1);
    cout << "pointPhysical:" << pointPhysical2D_1; cout << std::endl;

    file2D << pointReference2D_1;
    file2D << pointPhysical2D_1;

    // The following are just for the file, to be read by Mathematica.

    pointReference2D_1[0] = 0.5;
    pointReference2D_1[1] = 0.5;
    mapping2D.transform(pointReference2D_1,pointPhysical2D_1);
    file2D << pointReference2D_1;
    file2D << pointPhysical2D_1;

    pointReference2D_1[0] = -0.5;
    pointReference2D_1[1] = 0.5;
    mapping2D.transform(pointReference2D_1,pointPhysical2D_1);
    file2D << pointReference2D_1;
    file2D << pointPhysical2D_1;

    pointReference2D_1[0] = 0.5;
    pointReference2D_1[1] = -0.5;
    mapping2D.transform(pointReference2D_1,pointPhysical2D_1);
    file2D << pointReference2D_1;
    file2D << pointPhysical2D_1;

    pointReference2D_1[0] = -0.5;
    pointReference2D_1[1] = -0.5;
    mapping2D.transform(pointReference2D_1,pointPhysical2D_1);
    file2D << pointReference2D_1;
    file2D << pointPhysical2D_1;

    file2D.close();

    cout << "3D ~~~~~~~~~ 3D ~~~~~~~~~ 3D ~~~~~~~~~ 3D ~~~~~~~~~ 3D ~~~~~~~~~ 3D" << std::endl;

    // The ordering of the vertices is fundamental.
    // See comment on PhysicalCube.cpp

    std::ofstream file3D;
    file3D.open ("3Dmapping.txt");

    // Geometry::PointPhysical<3> point3D_1((double []){2.0,2.0,2.0});
    // Geometry::PointPhysical<3> point3D_2((double []){4.0,2.0,2.0});
    // Geometry::PointPhysical<3> point3D_3((double []){2.0,4.0,2.0});
    // Geometry::PointPhysical<3> point3D_4((double []){4.0,4.0,2.0});
    // Geometry::PointPhysical<3> point3D_5((double []){2.0,2.0,4.0});
    // Geometry::PointPhysical<3> point3D_6((double []){4.0,2.0,4.0});
    // Geometry::PointPhysical<3> point3D_7((double []){2.0,4.0,4.0});
    // Geometry::PointPhysical<3> point3D_8((double []){4.0,4.0,4.0});

    double dummy7[3]= {1.0,2.5,2.0};
    Geometry::PointPhysical<3> point3D_1(dummy7);
    
    double dummy8[3]= {4.0,2.5,2.0};
    Geometry::PointPhysical<3> point3D_2(dummy8);
    
    double dummy9[3]= {1.0,4.2,2.0};
    Geometry::PointPhysical<3> point3D_3(dummy9);
    
    double dummy10[3]= {4.0,4.2,2.0};
    Geometry::PointPhysical<3> point3D_4(dummy10);
    
    double dummy11[3]= {2.0,2.5,4.8};
    Geometry::PointPhysical<3> point3D_5(dummy11);
    
    double dummy12[3]= {4.0,2.5,4.8};
    Geometry::PointPhysical<3> point3D_6(dummy12);
    
    double dummy13[3]= {1.0,4.2,4.0};
    Geometry::PointPhysical<3> point3D_7(dummy13);
    
    double dummy14[3]= {4.0,4.2,4.0};
    Geometry::PointPhysical<3> point3D_8(dummy14);

    std::vector<Geometry::PointPhysical<3> > points3D;
    points3D.push_back(point3D_1);
    points3D.push_back(point3D_2);
    points3D.push_back(point3D_3);
    points3D.push_back(point3D_4);
    points3D.push_back(point3D_5);
    points3D.push_back(point3D_6);
    points3D.push_back(point3D_7);
    points3D.push_back(point3D_8);

    file3D << point3D_1;
    file3D << point3D_2;
    file3D << point3D_4;
    file3D << point3D_3;
    file3D << point3D_5;
    file3D << point3D_6;
    file3D << point3D_8;
    file3D << point3D_7;

    unsigned int indexes3D[8] = {0,1,2,3,4,5,6,7};
    const std::vector<unsigned int> globalIndexes3D(indexes3D,indexes3D + sizeof(indexes3D)/sizeof(indexes3D[0]));

    Geometry::ReferenceCube referenceCube;
    Geometry::PhysicalHexahedron physicalHexahedron(globalIndexes3D,points3D,&referenceCube);

    const Geometry::PhysicalHexahedron*const lele = &physicalHexahedron;

    cout << "Geometry::MappingSimpleHypercubeLinear<3> mapping(physicalCube)" << std::endl;
    Geometry::MappingToPhysHypercubeLinear<3> mapping3D(lele);

    cout << "mapping.transform(PointReference,PointPhysical)" << std::endl;
    Geometry::PointReference<3> pointReference3D_1;
    Geometry::PointPhysical<3> pointPhysical3D_1;

    pointReference3D_1[0] = -1.0;
    pointReference3D_1[1] = -1.0;
    pointReference3D_1[2] = -1.0;
    cout << "pointReference:" << pointReference3D_1; cout << std::endl;
    mapping3D.transform(pointReference3D_1,pointPhysical3D_1);
    cout << "pointPhysical:" << pointPhysical3D_1; cout << std::endl;
    file3D << pointReference3D_1;
    file3D << pointPhysical3D_1;

    pointReference3D_1[0] = -1.0;
    pointReference3D_1[1] = -1.0;
    pointReference3D_1[2] = +1.0;
    cout << "pointReference:" << pointReference3D_1; cout << std::endl;
    mapping3D.transform(pointReference3D_1,pointPhysical3D_1);
    cout << "pointPhysical:" << pointPhysical3D_1; cout << std::endl;
    file3D << pointReference3D_1;
    file3D << pointPhysical3D_1;

    pointReference3D_1[0] = -1.0;
    pointReference3D_1[1] = +1.0;
    pointReference3D_1[2] = -1.0;
    cout << "pointReference:" << pointReference3D_1; cout << std::endl;
    mapping3D.transform(pointReference3D_1,pointPhysical3D_1);
    cout << "pointPhysical:" << pointPhysical3D_1; cout << std::endl;
    file3D << pointReference3D_1;
    file3D << pointPhysical3D_1;

    pointReference3D_1[0] = -1.0;
    pointReference3D_1[1] = +1.0;
    pointReference3D_1[2] = +1.0;
    cout << "pointReference:" << pointReference3D_1; cout << std::endl;
    mapping3D.transform(pointReference3D_1,pointPhysical3D_1);
    cout << "pointPhysical:" << pointPhysical3D_1; cout << std::endl;
    file3D << pointReference3D_1;
    file3D << pointPhysical3D_1;

    pointReference3D_1[0] = +1.0;
    pointReference3D_1[1] = -1.0;
    pointReference3D_1[2] = -1.0;
    cout << "pointReference:" << pointReference3D_1; cout << std::endl;
    mapping3D.transform(pointReference3D_1,pointPhysical3D_1);
    cout << "pointPhysical:" << pointPhysical3D_1; cout << std::endl;
    file3D << pointReference3D_1;
    file3D << pointPhysical3D_1;

    pointReference3D_1[0] = +1.0;
    pointReference3D_1[1] = -1.0;
    pointReference3D_1[2] = +1.0;
    cout << "pointReference:" << pointReference3D_1; cout << std::endl;
    mapping3D.transform(pointReference3D_1,pointPhysical3D_1);
    cout << "pointPhysical:" << pointPhysical3D_1; cout << std::endl;
    file3D << pointReference3D_1;
    file3D << pointPhysical3D_1;

    pointReference3D_1[0] = +1.0;
    pointReference3D_1[1] = +1.0;
    pointReference3D_1[2] = -1.0;
    cout << "pointReference:" << pointReference3D_1; cout << std::endl;
    mapping3D.transform(pointReference3D_1,pointPhysical3D_1);
    cout << "pointPhysical:" << pointPhysical3D_1; cout << std::endl;
    file3D << pointReference3D_1;
    file3D << pointPhysical3D_1;

    pointReference3D_1[0] = +1.0;
    pointReference3D_1[1] = +1.0;
    pointReference3D_1[2] = +1.0;
    cout << "pointReference:" << pointReference3D_1; cout << std::endl;
    mapping3D.transform(pointReference3D_1,pointPhysical3D_1);
    cout << "pointPhysical:" << pointPhysical3D_1; cout << std::endl;
    file3D << pointReference3D_1;
    file3D << pointPhysical3D_1;

    return 0;*/
}
