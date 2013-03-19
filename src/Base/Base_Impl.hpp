namespace Base
{
    template<unsigned int DIM>
    class Base;

    template<unsigned int DIM>
    bool
    Base<DIM>::solve()
    {
        initialise();
        checkInitialisation();
        mesh_.move(); // just for testing
        return true;
    }


    template<unsigned int DIM>
    bool
    Base<DIM>::initialiseMeshMover(MeshMoverBaseT* meshMoverBase)
    {
        mesh_.setMeshMover(meshMoverBase);
        return true;
    }
    
    template<unsigned int DIM>
    bool
    Base<DIM>::checkInitialisation()
    {
    return true;
    }
}
