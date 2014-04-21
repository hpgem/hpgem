    //
    //  TypedefsInGeometry.hpp
    //
    //
    //  Created by Shavarsh Nurijanyan on 2/2/13.
    //
    //

#ifndef _GlobalNamespaceGeometry_hpp
#define _GlobalNamespaceGeometry_hpp

#include <limits> 

namespace Geometry
{
    const int ZeroD  = 0;
    const int OneD   = 1;
    const int TwoD   = 2;
    const int ThreeD = 3;
    const int FourD  = 4;
    
    const int MaxInteger=std::numeric_limits<int>::max();
    const double SmallerDoubleThanMinimalSizeOfTheMesh=1.0e-10;
    
//    template <unsigned int DIM>
//    class TypedefsInGeometry
//    {
//    public:
//        typedef double                       CoordTypeT;
//        typedef unsigned int                 IndexT;
//        
//        typedef std::vector<NodeT>           VectorOfNodes;
//        typedef std::vector<NodeT* >         VectorOfPtrNodesT;
//        typedef std::vector<CoordTypeT >     VectorOfCoordsT;
//
//        typedef std::list<NodeT >            ListOfNodesT;
//        typedef std::list<NodeT* >           ListOfPtrNodesT;
//
//        typedef PhysicalGeometry<DIM>        PhysicalGeometryT;
//        typedef ReferenceGeometry<DIM>       ReferenceGeometryT;
//
//        typedef int MatrixT; // temp. until we get (if we get) linear algebra.
//    };
};

#endif///_GlobalNamespaceGeometry_hpp
