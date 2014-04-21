#include "OutwardNormalVectorSign.hpp"
#include "MappingToRefLineToSquare.hpp"
#include "MappingToRefLineToTriangle.hpp"
#include "MappingToRefSquareToSquare.hpp"
#include "MappingToRefPointToLine.hpp"
#include "MappingToRefCubeToCube.hpp"
#include "MappingToRefSquareToCube.hpp"
#include "MappingToRefTriangleToTetrahedron.hpp"
    //#include "MappingToRefFaceToPyramid.hpp"
namespace Geometry
{

    // Added by M.T. Julianto, Feb 16, 2010
    double OutwardNormalVectorSign(const MappingReferenceToReference* const map)
    {
        if (dynamic_cast<const MappingToRefPointToLine0*>(map) ||
        	dynamic_cast<const MappingToRefLineToTriangle0*>(map) ||
	        dynamic_cast<const MappingToRefLineToTriangle2*>(map) ||
	        dynamic_cast<const MappingToRefLineToSquare0*>(map)   ||
	        dynamic_cast<const MappingToRefLineToSquare2*>(map) ||
        	dynamic_cast<const MappingToRefSquareToCube0*>(map) ||
	        dynamic_cast<const MappingToRefSquareToCube2*>(map) ||
	        dynamic_cast<const MappingToRefSquareToCube4*>(map))
        {
            return -1.;
        }
        else
        {
            if (dynamic_cast<const MappingToRefPointToLine1*>(map) ||
            	dynamic_cast<const MappingToRefLineToTriangle1*>(map) ||
                dynamic_cast<const MappingToRefLineToSquare1*>(map)   ||
                dynamic_cast<const MappingToRefLineToSquare3*>(map) ||
                dynamic_cast<const MappingToRefSquareToCube1*>(map) ||
				dynamic_cast<const MappingToRefSquareToCube3*>(map) ||
				dynamic_cast<const MappingToRefSquareToCube5*>(map) ||
				dynamic_cast<const MappingToRefTriangleToTetrahedron0*>(map) ||
				dynamic_cast<const MappingToRefTriangleToTetrahedron1*>(map) ||
				dynamic_cast<const MappingToRefTriangleToTetrahedron2*>(map) ||
				dynamic_cast<const MappingToRefTriangleToTetrahedron3*>(map)
                //                dynamic_cast<const MappingToRefFaceToPyramid0*>(map) ||
                //                dynamic_cast<const MappingToRefFaceToPyramid1*>(map) ||
                //                dynamic_cast<const MappingToRefFaceToPyramid2*>(map) ||
                //                dynamic_cast<const MappingToRefFaceToPyramid3*>(map) ||
                //                dynamic_cast<const MappingToRefFaceToPyramid4*>(map) ||
                //                dynamic_cast<const MappingToRefFaceToTriangularPrism0*>(map) ||
                //                dynamic_cast<const MappingToRefFaceToTriangularPrism1*>(map) ||
                //                dynamic_cast<const MappingToRefFaceToTriangularPrism2*>(map) ||
                //                dynamic_cast<const MappingToRefFaceToTriangularPrism3*>(map) ||
                //                dynamic_cast<const MappingToRefFaceToTriangularPrism4*>(map) ||
                //                dynamic_cast<const MappingToRefFaceToTetrahedron0*>(map) ||
                //                dynamic_cast<const MappingToRefFaceToTetrahedron1*>(map) ||
                //                dynamic_cast<const MappingToRefFaceToTetrahedron2*>(map) ||
                //                dynamic_cast<const MappingToRefFaceToTetrahedron3*>(map)
                		         )
            {
                return +1.;
            }
            else
            {
                throw "outwardNormalVectorSign - don't know this face2el-map";
            }
        }
    }

	// Added by Vijaya for 4D Normal vector.
    /*template <>
    double OutwardNormalVectorSign<4>(
	const MappingReferenceToReference<3, 4>* const map)
    {
        /*
        if (dynamic_cast<const MappingCubeToHyperCube0*>(map) || // from ReferenceHyperCube
            dynamic_cast<const MappingCubeToHyperCube2*>(map) ||
            dynamic_cast<const MappingCubeToHyperCube4*>(map) ||
            dynamic_cast<const MappingCubeToHyperCube6*>(map))
        {
            return -1.;
        }
        else
        {
            if (dynamic_cast<const MappingCubeToHyperCube1*>(map) ||
                dynamic_cast<const MappingCubeToHyperCube3*>(map) ||
                dynamic_cast<const MappingCubeToHyperCube5*>(map) ||
                dynamic_cast<const MappingCubeToHyperCube7*>(map))
            {
                return +1.;
            }
            else
            {
                throw "outwardNormalVectorSign<3> - don't know this face2el-map";
            }
        }
        *//*
        return 1;
    }*/
}
