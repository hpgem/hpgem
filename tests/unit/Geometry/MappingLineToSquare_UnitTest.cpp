#include "../../../src/Geometry/PointReference.hpp"
#include "../../../src/Geometry/Mappings/MappingToRefLineToSquare.hpp"
#include <vector>
int main()
{
    /* Test of all four LineToSquare mappings.
     *
     *
     */

    const Geometry::MappingToRefLineToSquare0& mappingLineToSquare0 = Geometry::MappingToRefLineToSquare0::Instance();
    const Geometry::MappingToRefLineToSquare1& mappingLineToSquare1 = Geometry::MappingToRefLineToSquare1::Instance();
    const Geometry::MappingToRefLineToSquare2& mappingLineToSquare2 = Geometry::MappingToRefLineToSquare2::Instance();
    const Geometry::MappingToRefLineToSquare3& mappingLineToSquare3 = Geometry::MappingToRefLineToSquare3::Instance();

    Geometry::PointReference<1> pr1;
    Geometry::PointReference<2> pr2;

    cout << "Geometry::MappingLineToSquare0";
    for (double x=-1.0; x<1.0; x+=0.1)
    {
        pr1[0] = x;
        mappingLineToSquare0.transform(pr1,pr2);
        cout << "\t" << pr1 << " -> "<< pr2 << "\n";
    }

    cout << "Geometry::MappingLineToSquare1";
    for (double x=-1.0; x<1.0; x+=0.1)
    {
        pr1[0] = x;
        mappingLineToSquare1.transform(pr1,pr2);
        cout << "\t" << pr1 << " -> "<< pr2 << "\n";
    }

    cout << "Geometry::MappingLineToSquare2";
    for (double x=-1.0; x<1.0; x+=0.1)
    {
        pr1[0] = x;
        mappingLineToSquare2.transform(pr1,pr2);
        cout << "\t" << pr1 << " -> "<< pr2 << "\n";
    }

    cout << "Geometry::MappingLineToSquare3";
    for (double x=-1.0; x<1.0; x+=0.1)
    {
        pr1[0] = x;
        mappingLineToSquare3.transform(pr1,pr2);
        cout << "\t" << pr1 << " -> "<< pr2 << "\n";
    }
    return 0;
}
