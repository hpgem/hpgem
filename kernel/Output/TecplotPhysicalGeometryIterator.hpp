#ifndef TECPLOTPHYSICALGEOMETRYITERATOR_HH
#define TECPLOTPHYSICALGEOMETRYITERATOR_HH

#include <vector>
#include "PhysicalGeometryAcceptor.hpp"

namespace Output
{
    /// TODO: Implement other geometries. For that we need the physical geometries.
    class TecplotPhysicalGeometryIterator: public PhysicalGeometryAcceptor
    {

    public:

        /// return reference to singleton instance:
        static TecplotPhysicalGeometryIterator& Instance()
        {
            static TecplotPhysicalGeometryIterator theInstance;
            return theInstance;
        }

        void acceptG(const Geometry::PhysicalGeometry* geo)
        {
                //std::cout<<(*geo);
            const Geometry::PhysicalLine* line= dynamic_cast<const Geometry::PhysicalLine*>(geo);
            if (line)
            {
                acceptLineGeometry(line);
            }
            else
            {
                const Geometry::PhysicalTriangle* triangle= dynamic_cast<const Geometry::PhysicalTriangle*>(geo);
                if (triangle)
                {
                    acceptTriangleGeometry(triangle);
                }
                else
                {
                    const Geometry::PhysicalQuadrilateral* quad= dynamic_cast<const Geometry::PhysicalQuadrilateral*>(geo);
                    if (quad)
                    {
                        acceptQuadrilateralGeometry(quad);
                    }
                    else
                    {
                        const Geometry::PhysicalTetrahedron* tetr= dynamic_cast<const Geometry::PhysicalTetrahedron*>(geo);
                        if (tetr)
                        {
                            acceptTetrahedronGeometry(tetr);
                        }
                        else
                        {
                            const Geometry::PhysicalPyramid* pyr= dynamic_cast<const Geometry::PhysicalPyramid*>(geo);
                            if (pyr)
                            {
                                acceptPyramidGeometry(pyr);
                            }
                            else
                            {
                                const Geometry::PhysicalTriangularPrism* trPrism= dynamic_cast<const Geometry::PhysicalTriangularPrism*>(geo);
                                if (trPrism)
                                {
                                    acceptTriangularPrismGeometry(trPrism);
                                }
                                else
                                {
                                    const Geometry::PhysicalHexahedron* hex= dynamic_cast<const Geometry::PhysicalHexahedron*>(geo);
                                    if (hex)
                                    {
                                        acceptHexahedronGeometry(hex);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        /// \brief Choose node sequence for the hypercube.
        //virtual void acceptHyperCubeGeometry(const Geometry::PhysicalHypercube&);

        /// \brief Choose node sequence for the hexahedron.
        virtual void acceptHexahedronGeometry(const Geometry::PhysicalHexahedron*);

        /// \brief Choose node sequence for the prism.
        virtual void acceptTriangularPrismGeometry(const Geometry::PhysicalTriangularPrism*);

        /// \brief Choose node sequence for the pyramid.
        virtual void acceptPyramidGeometry(const Geometry::PhysicalPyramid*);

        /// \brief Choose node sequence for the tetrahedron.
        virtual void acceptTetrahedronGeometry(const Geometry::PhysicalTetrahedron*);

        /// \brief Choose node sequence for quadrilateral.
        virtual void acceptQuadrilateralGeometry(const Geometry::PhysicalQuadrilateral*);

        /// \brief Choose node sequence for triangle.
        virtual void acceptTriangleGeometry(const Geometry::PhysicalTriangle*);

        /// \brief Choose node sequence for line.
        virtual void acceptLineGeometry(const Geometry::PhysicalLine*);

        /// \brief Check whether all nodes of current shape have been used.
        bool more() const;

        /// \brief Return the current node number.
        unsigned int getNodeNr();

    private:
        TecplotPhysicalGeometryIterator();
        TecplotPhysicalGeometryIterator(const TecplotPhysicalGeometryIterator&);
        TecplotPhysicalGeometryIterator& operator=(const TecplotPhysicalGeometryIterator&);
        virtual ~TecplotPhysicalGeometryIterator() { }

        typedef int InternalIndexType;
        typedef std::vector<unsigned int> VectorOfNodeIndexes;
        VectorOfNodeIndexes hypercubeNodes;
        VectorOfNodeIndexes hexahedronNodes;
        VectorOfNodeIndexes cubeNodes;
        VectorOfNodeIndexes triangularPrismNodes;
        VectorOfNodeIndexes pyramidNodes;
        VectorOfNodeIndexes tetrahedronNodes;
        VectorOfNodeIndexes quadrilateralNodes;
        VectorOfNodeIndexes triangleNodes;
        VectorOfNodeIndexes lineNodes;

        VectorOfNodeIndexes* currentSequencePtr;
        InternalIndexType currentNode;
    };
}
#endif
