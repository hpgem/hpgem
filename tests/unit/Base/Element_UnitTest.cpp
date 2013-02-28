#include "../Base/Element.hpp"
#include "../Geometry/Point.hpp"
#include "../Geometry/PointPhysical.hpp"
#include "../Geometry/ElementGeometry.hpp"
#include "../Geometry/PhysicalGeometry.hpp"
#include "../Geometry/ReferenceSquare.hpp"
#include "../Geometry/MappingSimpleCubeNLinear.hpp"
#include "../LinearAlgebra/NumericalVector.hpp"
#include <vector>

using namespace std;

int main()
{
    //Geometry::PhysicalGeometry<3> physicalGeometry;
    std::vector<int> PointsIndexes(4);
    for (i==0;i<4;i++) {PointsIndexes(i)=i;}
    
        
    
   // LinearAlgebra::NumericalVector outerPoint(2);
   // outerPoint[0] = 1.0;
   // outerPoint[1] = 1.0;
   // points[0] = new Geometry::PointPhysical<2>(outerPoint);
   // outerPoint[0] = 2.0;
   // outerPoint[1] = 1.0;
   // points[1] = new Geometry::PointPhysical<2>(outerPoint);
   // outerPoint[0] = 2.0;
   // outerPoint[1] = 2.0;
   // points[2] = new Geometry::PointPhysical<2>(outerPoint);
   // outerPoint[0] = 2.0;
   // outerPoint[1] = 2.0;
   // points[3] = new Geometry::PointPhysical<2>(outerPoint);
   // cout << *points[0] << std::endl;
   // cout << *points[1] << std::endl;
   // cout << *points[2] << std::endl;
   // cout << *points[3] << std::endl;
    
    

    Geometry::ReferenceSquare referenceSquare;

    cout << referenceSquare;

    Geometry::PhysicalGeometry<2> physicalGeometry(points);
    const Geometry::PhysicalGeometry<2>*const lala = &physicalGeometry;
    Geometry::MappingSimpleCubeNLinear<2> mapping(lala);

    Base::Element<2> element(points,&referenceSquare,&mapping);
    //Base::Element<1> element1D;
    return 0;
}
