namespace Base
{
    template<unsigned int DIM>
    class Base;

    template<unsigned int DIM>
    void
    Base<DIM>::output()
    {
        std::ofstream file2D;
        file2D.open ("out.dat");
        int dimensionsToWrite[2] = {0,1};
        Output::TecplotDiscontinuousSolutionWriter<2> out(file2D,"RectangularMesh",dimensionsToWrite,"xy");
        out.write(mesh_,"holi",false);
    }

    template<unsigned int DIM>
    bool
    Base<DIM>::solve()
    {
        // Do complex reimmann integrals of fifth order using guassian-levitsky spline
        // adaptative methods.
        // mesh_.move(); // just for testing
        return true;
    }

    template<unsigned int DIM>
    bool
    Base<DIM>::initializeMesh()
    {
        Geometry::PointPhysical<2> bottomLeft, topLeft;
        std::vector<unsigned int> numElementsOneD(2);
        bottomLeft[0] = 0;
        bottomLeft[1] = 0;
        topLeft[0] = 1;
        topLeft[1] = 1;
        numElementsOneD[0] = 8;
        numElementsOneD[1] = 8;
        mesh_.createRectangularMesh(bottomLeft, topLeft, numElementsOneD);
        return true;
    }

    template<unsigned int DIM>
    bool
    Base<DIM>::initializeMeshMover(MeshMoverBaseT* meshMoverBase)
    {
        mesh_.setMeshMover(meshMoverBase);
        return true;
    }
}
