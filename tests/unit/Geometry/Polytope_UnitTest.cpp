#include "../../../src/Geometry/ReferenceCube.hpp"
int main()
{
    //Geometry::Polytope<3> polytope = new Geometry::Polytope<3>(3);
    Geometry::ReferenceCube cube;
    cout << "getName(): "<< cube.getName() << std::endl;

    cout << "output operator: "<< cube;

    cout << "\t" << "getNumberOfPoints()="<< cube.getNumberOfPoints() << std::endl;
    cout << "\t" << "getPoint():" << std::endl;
    Geometry::Point<3> p;

    for (int i = 0; i<cube.numberOfVertices();++i)
    {
        cube.getVertexPoint(i, p);
        
        cout << "Get the" << i<<"-th of the  "<<cube.getName()<<" "<< p<<endl;
    }

    cout << "\t" << "getPointCoordinates():" << std::endl;
    Geometry::ReferenceCube::PointT::VectorOfCoordsT n;
    cube.getPointCoordinates(0,n);
    cout << "\t" << "(" << n[0] << " " << n[1] << " " << n[2] << ")" << std::endl;
    cube.getPointCoordinates(1,n);
    cout << "\t" << "(" << n[0] << " " << n[1] << " " << n[2] << ")" << std::endl;
    cube.getPointCoordinates(2,n);
    cout << "\t" << "(" << n[0] << " " << n[1] << " " << n[2] << ")" << std::endl;

    cout << "\t" << "addPoints(point):" << std::endl;

    Geometry::Point<3> point((double[]){1.2, 1.3, 1.5});
    cube.addPoint(point);
    
    for (int i = 0; i<cube.numberOfVertices();++i)
    {
        cube.getVertexPoint(i, point);
        
        cout << "Get the" << i<<"-th of the  "<<cube.getName()<<" "<< point<<endl;
    }
    
    cout << "\t" << "setPoint(1,point):" << std::endl;
    cube.setPoint(1,point);
    for (int i = 0; i<cube.numberOfVertices();++i)
    {
        cube.getVertexPoint(i, p);
        
        cout << "Get the" << i<<"-th of the  "<<cube.getName()<<" "<< p<<endl;
    }

        // cout << "\t" << "calcElementDiameter(1,point)=" << cube.calcElementDiameter() << std::endl;

    return 0;
}
