#include "../../../src/Geometry/PointReference.hpp"
#include "../../../src/Geometry/Mappings/MappingToRefPointToLine.hpp"
#include <vector>
int main()
{
    /* Test of all four PointToLine mappings.
     *
     * The point is just 0.
     */

    const Geometry::MappingToRefPointToLine0& mappingPointToLine0 = Geometry::MappingToRefPointToLine0::Instance();
    const Geometry::MappingToRefPointToLine1& mappingPointToLine1 = Geometry::MappingToRefPointToLine1::Instance();

    Geometry::PointReference<0> pr1;
    Geometry::PointReference<1> pr2;

    cout << "Geometry::MappingPointToLine0";
    mappingPointToLine0.transform(pr1,pr2);
    cout << "\t" << pr1 << " -> "<< pr2 << "\n";

    cout << "Geometry::MappingPointToLine1";
    mappingPointToLine1.transform(pr1,pr2);
    cout << "\t" << pr1 << " -> "<< pr2 << "\n";

    return 0;
}
