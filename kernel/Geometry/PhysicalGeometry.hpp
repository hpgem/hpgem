#ifndef PHYSICALGEOMETRY_H_
#define PHYSICALGEOMETRY_H_

#include "PointPhysical.hpp"
#include "GlobalNamespaceGeometry.hpp"
#include "ReferenceGeometry.hpp"
#include <vector>
//#include "Output/PhysicalGeometryAcceptor.hpp"

namespace Geometry
{


    class PhysicalGeometry
    {
        /*! \class PhysicalGeometry
         * \brief PhysicalGeometry describes an actual physical shape in real space.
         * \details
         * You shouldn't create a PhysicalGeometry (although it is possible), but a particular
         * Physical<Shape>.
         *
         * It contains only the global indexes of its points in globalNodeIndexes_. These global
         * indexes refer to the global node container, of which every PhysicalGeometry has a
         * reference: nodes_.
         *
         * It also contains a reference to the corresponding referenceGeometry.
         *
         * ~ Point is the name of a class.               ~
         * ~ Node is a point that belongs to a geometry. ~
         */
        public:
                   
            typedef unsigned int                    PointIndexT;
            typedef unsigned int                    LocalNodeIndexT;
            typedef double                          CoordT;
            typedef ReferenceGeometry          ReferenceGeometryT;
            typedef PointPhysical              PointPhysicalT;
            typedef std::vector<PointIndexT>        VectorOfPointIndexesT;
            typedef std::vector<CoordT >            VectorOfCoordsT;
            typedef std::vector<PointPhysicalT>     VectorOfPhysicalPointsT;
            typedef PhysicalGeometry           PhysicalGeometryT;
            typedef std::vector<unsigned int>       GlobalNodeSetOnTheFace;
        

        public:

            /// \brief Constructor gets indexes of the nodes, a reference to the node container, and a pointer to the corresponding reference geometry.
            PhysicalGeometry(const VectorOfPointIndexesT& globalNodeIndexes,
                             const VectorOfPhysicalPointsT& nodes,
                             const ReferenceGeometryT* const  refG):
                globalNodeIndexes_(globalNodeIndexes),
                nodes_(nodes),
                refGeometry_(refG)
        
            {}
                
            virtual ~PhysicalGeometry() {}

            /// \brief Returns a pointer to the container of the global node indexes.
            //VectorOfPointIndexesT&          getNodeIndexes() {return globalNodeIndexes_;}

            /// \brief Returns a constant pointer to the container of the global node indexes.
            const VectorOfPointIndexesT&    getNodeIndexes() const {return globalNodeIndexes_;}

            /// \brief Returns a pointer to the global container of nodes.
            const VectorOfPhysicalPointsT&        getNodes() {return nodes_;}

            /// \brief Returns a constant pointer of the global container of nodes
            const VectorOfPhysicalPointsT&  getNodes() const {return nodes_;}

            /// \brief Returns the name of the particular geometry.
            virtual std::string             getName() const = 0;

            /// \brief Given a local index relative to globalNodeIndexes_, return the global node index.
            PointIndexT                     getNodeIndex(int localIndex) const
                                            {return globalNodeIndexes_[localIndex];}

            /// \brief Given a global index, returns a pointer to the corresponding point.
            const PointPhysicalT*                 getNodePtr(const int globalIndex) const
                                            {return &(nodes_[globalIndex]);}

            /// \brief Returns the number of nodes of this geometry.
            unsigned int                    getNumberOfNodes() const
                                            {return globalNodeIndexes_.size();}

            /// \brief Given a local index, assigns the physical coordinates of the corresponding point.
            // MTJ: TODO: this should be renamed to getLocalNodeCoordinates.............
            ///\TODO remove duplicate code
            void                            getNodeCoordinates(const int localIndex, PointPhysicalT& coords) const
            {coords = (nodes_)[globalNodeIndexes_[localIndex]].getCoordinates();}

            /// \brief Given a local index, assigns the physical coordinates of the corresponding point.
            void                            getLocalNodeCoordinates(const int localIndex, PointPhysicalT& coords) const
            {coords = (nodes_)[globalNodeIndexes_[localIndex]].getCoordinates();}

            /// \brief Given a global index, assigns the physical coordinates of the corresponding point.
            void                            getGlobalNodeCoordinates(const int globalIndex, PointPhysicalT& coords) const
            {coords = (nodes_)[globalIndex].getCoordinates();}

            /// \brief Given a local face index, return the global indices of the entities contained on that face.
            virtual void                    getGlobalFaceNodeIndices(const PointIndexT, VectorOfPointIndexesT&) const = 0;

            /// \brief Given a local face index, return the local indices of the entities contained on that face.
            virtual void                    getLocalFaceNodeIndices(const PointIndexT, VectorOfPointIndexesT&) const = 0;

            /// \brief Returns the number of faces via a call to ReferenceGeometry->getNrOfCodim1Entities();
            virtual unsigned int            getNrOfFaces() const = 0;
        
            /// \brief Returns a reference to the corresponding reference geometry.
            const ReferenceGeometryT* const getRefGeometry() const {return refGeometry_;}
        
            /// \brief Output operator
            friend std::ostream& operator<<(std::ostream& os,
                    const PhysicalGeometry& physicalGeometry)
            {
                os << "PhysicalGeometry=( ";
                
                for (int i = 0; i < physicalGeometry.getNumberOfNodes(); i++)
                {   
                    os << physicalGeometry.getNodeIndex(i) << " ";
                }
                os << ')' << std::endl;

                return os;
            }

        protected:
            /// Reference to the global node container.
            const VectorOfPhysicalPointsT&  nodes_;

            /// Reference to the container of global indexes of the nodes, relative to nodes_.
            VectorOfPointIndexesT     globalNodeIndexes_;
        
            const ReferenceGeometryT* const refGeometry_;
    };
};
#endif /* PHYSICALGEOMETRY_H_ */
