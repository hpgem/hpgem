#include "../../../kernel/Geometry/Point.hpp"
#include "../../../kernel/Geometry/PointPhysical.hpp"
#include "../../../kernel/Geometry/ElementGeometry.hpp"
#include "../../../kernel/Geometry/PhysicalGeometry.hpp"
#include "../../../kernel/Geometry/ReferenceSquare.hpp"
#include "../../../kernel/LinearAlgebra/NumericalVector.hpp"
#include "../../../kernel/Geometry/Mappings/MappingToPhysHypercubeLinear.hpp"
#include <vector>


using namespace std;

int main()
{
//    //Geometry::PhysicalGeometry<3> physicalGeometry;
//    std::vector<Geometry::PointPhysical<2>* > points(4);
//    LinearAlgebra::NumericalVector outerPoint(2);
//    outerPoint[0] = 1.0;
//    outerPoint[1] = 1.0;
//    points[0] = new Geometry::PointPhysical<2>(outerPoint);
//    outerPoint[0] = 2.0;
//    outerPoint[1] = 1.0;
//    points[1] = new Geometry::PointPhysical<2>(outerPoint);
//    outerPoint[0] = 2.0;
//    outerPoint[1] = 2.0;
//    points[2] = new Geometry::PointPhysical<2>(outerPoint);
//    outerPoint[0] = 2.0;
//    outerPoint[1] = 2.0;
//    points[3] = new Geometry::PointPhysical<2>(outerPoint);
//    cout << *points[0] << std::endl;
//    cout << *points[1] << std::endl;
//    cout << *points[2] << std::endl;
//    cout << *points[3] << std::endl;
//    //
////    cout << "Creating Physical Geometry:" << std::endl;
////    Geometry::PhysicalGeometry<2> physicalGeometry(points);
////    cout << physicalGeometry;
//
//    cout << "Creating Reference Geometry:" << std::endl;
//    Geometry::ReferenceSquare referenceSquare;
//    cout << referenceSquare;
//
//    cout << "Creating The Mapping:" << std::endl;
//    const Geometry::PhysicalGeometry<2>*const pG = &physicalGeometry;
//    Geometry::MappingSimpleCubeNLinear<2> mapping(pG);
//
//    cout << "ElementGeometry Constructor: ";
//    Geometry::ElementGeometry<2> elementGeometry(points,&referenceSquare,&mapping);
//    cout << "Works :)" << std::endl;
//
//    cout << "ElementGeometry setReferenceToPhysicalMap: ";
//    elementGeometry.setReferenceToPhysicalMap(&mapping);
//    cout << "Works :)" << std::endl;
//
//    cout << "ElementGeometry setReferenceGeometry: ";
//    Geometry::ReferenceGeometry<2>* referenceGeometryPtr = &referenceSquare;
//    elementGeometry.setReferenceGeometry(referenceGeometryPtr);
//    cout << "Works :)" << std::endl;
//
//    cout << "ElementGeometry getReferenceToPhysicalMap: ";
//    Geometry::MappingReferenceToPhysical<2>* mappingPtr;
//    mappingPtr = elementGeometry.getReferenceToPhysicalMap();
//    cout << "Works :)" << std::endl;
//
//    cout << "ElementGeometry getPhysicalGeometry: ";
//    Geometry::PhysicalGeometry<2> physicalGeometry2(points);
//    physicalGeometry2 = elementGeometry.getPhysicalGeometry();
//    cout << "Works :)" << std::endl;
//
//    cout << "ElementGeometry getReferenceGeometry: ";
//    Geometry::ReferenceGeometry<2>* referenceGeometryPtr2;
//    referenceGeometryPtr2 = elementGeometry.getReferenceGeometry();
//    cout << "Works :)" << std::endl;
//
//    cout << "ElementGeometry referenceToPhysical: ";
//    Geometry::PointPhysical<2> pointPhysical;
//    Geometry::PointReference<2> pointReference((double []){0.5, 0.5});
//    cout << pointReference;
//    elementGeometry.referenceToPhysical(pointReference,pointPhysical);
//    cout << "Works :)" << std::endl;

    return 0;
}
