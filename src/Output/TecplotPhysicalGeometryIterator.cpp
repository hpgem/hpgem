/*
 * TecplotPhysicalGeometryIterator.cpp
 *
 *  Created on: Feb 16, 2013
 *      Author: nicorivas
 */


#include "TecplotPhysicalGeometryIterator.hpp"

namespace Output
{
    TecplotPhysicalGeometryIterator::TecplotPhysicalGeometryIterator()
    {
        // sort in the sequences in BACKWARD order!
        // see tecplot manual chapter 4.3

        // 4D FE hyper cube: ORDER IS GUESS
        /*
        hypercubeNodes.push_back(14);
        hypercubeNodes.push_back(15);
        hypercubeNodes.push_back(13);
        hypercubeNodes.push_back(12);
        hypercubeNodes.push_back(10);
        hypercubeNodes.push_back(11);
        hypercubeNodes.push_back(9);
        hypercubeNodes.push_back(8);

        hypercubeNodes.push_back(6);
        hypercubeNodes.push_back(7);
        hypercubeNodes.push_back(5);
        hypercubeNodes.push_back(4);
        hypercubeNodes.push_back(2);
        hypercubeNodes.push_back(3);
        hypercubeNodes.push_back(1);
        hypercubeNodes.push_back(0);

        // 3D FE volumes: (all as bricks)
        */
        hexahedronNodes.push_back(6);
        hexahedronNodes.push_back(7);
        hexahedronNodes.push_back(5);
        hexahedronNodes.push_back(4);
        hexahedronNodes.push_back(2);
        hexahedronNodes.push_back(3);
        hexahedronNodes.push_back(1);
        hexahedronNodes.push_back(0);

        triangularPrismNodes.push_back(5);
        triangularPrismNodes.push_back(5);
        triangularPrismNodes.push_back(4);
        triangularPrismNodes.push_back(3);
        triangularPrismNodes.push_back(2);
        triangularPrismNodes.push_back(2);
        triangularPrismNodes.push_back(1);
        triangularPrismNodes.push_back(0);

        pyramidNodes.push_back(0);
        pyramidNodes.push_back(0);
        pyramidNodes.push_back(0);
        pyramidNodes.push_back(0);
        pyramidNodes.push_back(3);
        pyramidNodes.push_back(4);
        pyramidNodes.push_back(2);
        pyramidNodes.push_back(1);

        tetrahedronNodes.push_back(3);
        tetrahedronNodes.push_back(3);
        tetrahedronNodes.push_back(3);
        tetrahedronNodes.push_back(3);
        tetrahedronNodes.push_back(2);
        tetrahedronNodes.push_back(2);
        tetrahedronNodes.push_back(1);
        tetrahedronNodes.push_back(0);

        // 2D FE-surfaces:
        quadrilateralNodes.push_back(2);
        quadrilateralNodes.push_back(3);
        quadrilateralNodes.push_back(1);
        quadrilateralNodes.push_back(0);

        triangleNodes.push_back(2);
        triangleNodes.push_back(1);
        triangleNodes.push_back(0);
        triangleNodes.push_back(0); // double so that we have four points

        // for the FE-QUADRILATERAL element type

        // 1D data - is this supported by tecplot?
        lineNodes.push_back(1);
        lineNodes.push_back(0);
    }

//    template<unsigned int DIM>
//    void TecplotPhysicalGeometryIterator::acceptG(const Geometry::PhysicalGeometry<DIM>* geo)
//    {
//        const Geometry::PhysicalLine* line= dynamic_cast<const Geometry::PhysicalLine*>(geo);
//        if (line)
//        {
//            acceptLineGeometry(line);
//        }
//        else
//        {
//            const Geometry::PhysicalTriangle* triangle= dynamic_cast<const Geometry::PhysicalTriangle*>(geo);
//            if (triangle)
//            {
//                acceptTriangleGeometry(triangle);
//            }
//            else
//            {
//                const Geometry::PhysicalQuadrilateral* quad= dynamic_cast<const Geometry::PhysicalQuadrilateral*>(geo);
//                if (quad)
//                {
//                    acceptQuadrilateralGeometry(quad);
//                }
//                else
//                {
//                    const Geometry::PhysicalTetrahedron* tetr= dynamic_cast<const Geometry::PhysicalTetrahedron*>(geo);
//                    if (tetr)
//                    {
//                        acceptTetrahedronGeometry(tetr);
//                    }
//                    else
//                    {
//                        const Geometry::PhysicalPyramid* pyr= dynamic_cast<const Geometry::PhysicalPyramid*>(geo);
//                        if (pyr)
//                        {
//                            acceptPyramidGeometry(pyr);
//                        }
//                        else
//                        {
//                            const Geometry::PhysicalTriangularPrism* trPrism= dynamic_cast<const Geometry::PhysicalTriangularPrism*>(geo);
//                            if (trPrism)
//                            {
//                                acceptTriangularPrismGeometry(trPrism);
//                            }
//                            else
//                            {
//                                const Geometry::PhysicalHexahedron* hex= dynamic_cast<const Geometry::PhysicalHexahedron*>(geo);
//                                if (hex)
//                                {
//                                    acceptHexahedronGeometry(hex);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }

    /*
    void TecplotPhysicalGeometryIterator::acceptHyperCubeGeometry(const Geometry::PhysicalHypercube&)
    {
        currentSequencePtr = &hypercubeNodes;
        currentNode = hypercubeNodes.size() - 1;
    }
    */

    
    void TecplotPhysicalGeometryIterator::acceptHexahedronGeometry(const Geometry::PhysicalHexahedron*)
    {
        currentSequencePtr = &hexahedronNodes;
        currentNode = hexahedronNodes.size() - 1;
    }

    void TecplotPhysicalGeometryIterator::acceptTriangularPrismGeometry(const Geometry::PhysicalTriangularPrism*)
    {
        currentSequencePtr = &triangularPrismNodes;
        currentNode = triangularPrismNodes.size() - 1;
    }

    void TecplotPhysicalGeometryIterator::acceptPyramidGeometry(const Geometry::PhysicalPyramid*)
    {
        currentSequencePtr = &pyramidNodes;
        currentNode = pyramidNodes.size() - 1;
    }

    void TecplotPhysicalGeometryIterator::acceptTetrahedronGeometry(const Geometry::PhysicalTetrahedron*)
    {
        currentSequencePtr = &tetrahedronNodes;
        currentNode = tetrahedronNodes.size() - 1;
    }

    void TecplotPhysicalGeometryIterator::acceptQuadrilateralGeometry(const Geometry::PhysicalQuadrilateral*)
    {
        currentSequencePtr = &quadrilateralNodes;
        currentNode = quadrilateralNodes.size() - 1; // inverse transverse
    }

    void TecplotPhysicalGeometryIterator::acceptTriangleGeometry(const Geometry::PhysicalTriangle*)
    {
        currentSequencePtr = &triangleNodes;
        currentNode = triangleNodes.size() - 1;
    }

    void TecplotPhysicalGeometryIterator::acceptLineGeometry(const Geometry::PhysicalLine*)
    {
        currentSequencePtr = &lineNodes;
        currentNode = lineNodes.size() - 1;
    }

    bool TecplotPhysicalGeometryIterator::more() const
    {
        // this test is why InternalIndexType has to be int !!!
        return (currentNode >= 0);
    }

    unsigned int TecplotPhysicalGeometryIterator::getNodeNr()
    {
        return currentSequencePtr->operator[](currentNode--);
    }
}
