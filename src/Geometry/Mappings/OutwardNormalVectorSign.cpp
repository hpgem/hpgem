#include "OutwardNormalVectorSign.hpp"
#include "MappingToRefLineToSquare.hpp"
#include "MappingToRefSquareToSquare.hpp"
#include "MappingToRefPointToLine.hpp"
namespace Geometry
{

    // Added by M.T. Julianto, Feb 16, 2010
    template <>
    double OutwardNormalVectorSign<1>(const MappingReferenceToReference<0, 1>* const map)
    {
        if (dynamic_cast<const MappingToRefPointToLine0*>(map))
        {
            return -1.;
        }
        else
        {
            if (dynamic_cast<const MappingToRefPointToLine1*>(map))
            {
                return +1.;
            }
            else
            {
                throw "outwardNormalVectorSign<1> - don't know this face2el-map";
            }
        }
    }
    
    template <>
    double OutwardNormalVectorSign<2>(
	const MappingReferenceToReference<1, 2>* const map)
    {
        if (//dynamic_cast<const Line2TriangleMapping0*>(map) ||
	        //dynamic_cast<const Line2TriangleMapping2*>(map) ||
	        dynamic_cast<const MappingToRefLineToSquare0*>(map)   ||
	        dynamic_cast<const MappingToRefLineToSquare2*>(map))
        {
            return -1.;
        }
        else
        {
            if (//dynamic_cast<const Line2TriangleMapping1*>(map) ||
                dynamic_cast<const MappingToRefLineToSquare1*>(map)   ||
                dynamic_cast<const MappingToRefLineToSquare3*>(map))
            {
                return +1.;
	        }
            else
	        {
	            throw "outwardNormalVectorSign<2> - don't know this face2el-map";
	        }
        }

    }

    template <>
    double OutwardNormalVectorSign<3>(
	const MappingReferenceToReference<2, 3>* const map)
    {
        /*
        if (dynamic_cast<const MappingSquareToCube0*>(map) ||
	        dynamic_cast<const MappingSquareToCube2*>(map) ||
	        dynamic_cast<const MappingSquareToCube4*>(map))
        {
            return -1.;
        }
        else
        {
            if (dynamic_cast<const MappingSquareToCube1*>(map) ||
                dynamic_cast<const MappingSquareToCube3*>(map) ||
                dynamic_cast<const MappingSquareToCube5*>(map) ||
                dynamic_cast<const MappingFaceToPyramid0*>(map) ||
                dynamic_cast<const MappingFaceToPyramid1*>(map) ||
                dynamic_cast<const MappingFaceToPyramid2*>(map) ||
                dynamic_cast<const MappingFaceToPyramid3*>(map) ||
                dynamic_cast<const MappingFaceToPyramid4*>(map) ||
                dynamic_cast<const MappingFaceToTriangularPrism0*>(map) ||
                dynamic_cast<const MappingFaceToTriangularPrism1*>(map) ||
                dynamic_cast<const MappingFaceToTriangularPrism2*>(map) ||
                dynamic_cast<const MappingFaceToTriangularPrism3*>(map) ||
                dynamic_cast<const MappingFaceToTriangularPrism4*>(map) ||
                dynamic_cast<const MappingFaceToTetrahedron0*>(map) ||
                dynamic_cast<const MappingFaceToTetrahedron1*>(map) ||
                dynamic_cast<const MappingFaceToTetrahedron2*>(map) ||
                dynamic_cast<const MappingFaceToTetrahedron3*>(map))
            {
                return +1.;
            }
            else
            {
                throw "outwardNormalVectorSign<3> - don't know this face2el-map";
            }
        }
        */
        return 0;
    }

	// Added by Vijaya for 4D Normal vector.
    template <>
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
        */
        return 0;
    }
}
