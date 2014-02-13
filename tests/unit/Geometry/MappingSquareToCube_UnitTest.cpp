#include "Geometry/PointReference.hpp"
#include "Geometry/Mappings/MappingToRefSquareToCube.hpp"
#include <vector>

int main()
{
    const Geometry::MappingToRefSquareToCube0& mappingSquareToCube0 = Geometry::MappingToRefSquareToCube0::Instance();
    const Geometry::MappingToRefSquareToCube1& mappingSquareToCube1 = Geometry::MappingToRefSquareToCube1::Instance();
    const Geometry::MappingToRefSquareToCube2& mappingSquareToCube2 = Geometry::MappingToRefSquareToCube2::Instance();
    const Geometry::MappingToRefSquareToCube3& mappingSquareToCube3 = Geometry::MappingToRefSquareToCube3::Instance();
    const Geometry::MappingToRefSquareToCube4& mappingSquareToCube4 = Geometry::MappingToRefSquareToCube4::Instance();
    const Geometry::MappingToRefSquareToCube5& mappingSquareToCube5 = Geometry::MappingToRefSquareToCube5::Instance();

    Geometry::PointReference pr1(2);
    Geometry::PointReference pr2(3);

    cout << "Geometry::MappingSquareToCube0 ~ o ~ Geometry::MappingSquareToCube0" << "\n";
    for (double x=-1.0; x<=1.0; x+=0.5)
    {
        for (double y=-1.0; y<=1.0; y+=0.5)
        {
            pr1[0] = x;
            pr1[1] = y;
            mappingSquareToCube0.transform(pr1,pr2);
            cout << pr1 << " -> "<< pr2 << "\n";
        }
        cout << "\n";
    }

    cout << "Geometry::MappingSquareToCube1 ~ o ~ Geometry::MappingSquareToCube1" << "\n";
    for (double x=-1.0; x<=1.0; x+=0.5)
    {
        for (double y=-1.0; y<=1.0; y+=0.5)
        {
            pr1[0] = x;
            pr1[1] = y;
            mappingSquareToCube1.transform(pr1,pr2);
            cout << pr1 << " -> "<< pr2 << "\n";
        }
        cout << "\n";
    }

    cout << "Geometry::MappingSquareToCube2 ~ o ~ Geometry::MappingSquareToCube2" << "\n";
    for (double x=-1.0; x<=1.0; x+=0.5)
    {
        for (double y=-1.0; y<=1.0; y+=0.5)
        {
            pr1[0] = x;
            pr1[1] = y;
            mappingSquareToCube2.transform(pr1,pr2);
            cout << pr1 << " -> "<< pr2 << "\n";
        }
        cout << "\n";
    }

    cout << "Geometry::MappingSquareToCube3 ~ o ~ Geometry::MappingSquareToCube3" << "\n";
    for (double x=-1.0; x<=1.0; x+=0.5)
    {
        for (double y=-1.0; y<=1.0; y+=0.5)
        {
            pr1[0] = x;
            pr1[1] = y;
            mappingSquareToCube3.transform(pr1,pr2);
            cout << pr1 << " -> "<< pr2 << "\n";
        }
        cout << "\n";
    }

    cout << "Geometry::MappingSquareToCube4 ~ o ~ Geometry::MappingSquareToCube4" << "\n";
    for (double x=-1.0; x<=1.0; x+=0.5)
    {
        for (double y=-1.0; y<=1.0; y+=0.5)
        {
            pr1[0] = x;
            pr1[1] = y;
            mappingSquareToCube4.transform(pr1,pr2);
            cout << pr1 << " -> "<< pr2 << "\n";
        }
        cout << "\n";
    }

    cout << "Geometry::MappingSquareToCube5 ~ o ~ Geometry::MappingSquareToCube5" << "\n";
    for (double x=-1.0; x<=1.0; x+=0.5)
    {
        for (double y=-1.0; y<=1.0; y+=0.5)
        {
            pr1[0] = x;
            pr1[1] = y;
            mappingSquareToCube5.transform(pr1,pr2);
            cout << pr1 << " -> "<< pr2 << "\n";
        }
        cout << "\n";
    }
    
    return 0;
}
