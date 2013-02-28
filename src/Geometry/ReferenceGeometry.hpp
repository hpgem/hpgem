#ifndef REFERENCEGEOMETRY_HH
#define REFERENCEGEOMETRY_HH

#include "Geometry/PointReference.hpp"
#include "Geometry/Mappings/MappingCodimensions.hpp"
#include "Geometry/Mappings/RefinementMapping.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"

#include <iostream>
#include <vector>
using std::cout;
using std::endl;

namespace QuadratureRules
{
    // forward declaration
    template <unsigned int dim>
    class GaussQuadratureRule;
}

namespace Geometry
{
    
    enum TypeOfReferenceGeometry
    {
        POINT,
        LINE,
        TRIANGLE,
        SQUARE,
        TETRAHEDRON,
        PYRAMID,
        CUBE,
        TRIANGULARPRISM,
        HYPERCUBE
    };
    
    template <unsigned int DIM>
    class ReferenceGeometry :
            public RefinementMapping<DIM>,
            public MappingCodimensions<DIM>
    {
    /*! \class ReferenceGeometry
     * \brief ReferenceGeometry stores a the information of a unitary geometry where the integration is made.
     * \details
     * ReferenceGeometry stores the information of the corresponding reference object, which
     * are used for integration routines. It is a pure virtual class; an interface for every
     * particular reference geometry. The information is a container of the actual reference
     * points.
     *
     * Constructors are protected, to avoid the creation of two identical physical geometries.
     * Only one of every needed type is necessary.
     */
    public:
        /// \bug this is a workaround for a g++ bug. Should read using typenames;
        typedef std::string                                         String;
        typedef unsigned int                                        IndexT;
        typedef typename Geometry::PointReference<DIM>              PointReferenceT;
        typedef typename std::vector<PointReferenceT >              VectorOfReferencePointsT;
        typedef typename VectorOfReferencePointsT::iterator         iterator;
        typedef typename VectorOfReferencePointsT::const_iterator   const_iterator;
        typedef std::vector<IndexT>                                 ListOfIndexesT;

    public:

        virtual ~ReferenceGeometry(){};

        /// \brief Check whether a given point is within the ReferenceGeometry.
        virtual bool            isInternalPoint(const PointReferenceT& point) const = 0;

        /// \brief Each reference geometry knows its center of mass.
        virtual void            getCenter(PointReferenceT& point) const = 0;

        /// \brief Return number of nodes of the reference shape.
        virtual unsigned int    getNumberOfNodes() const {return points_.size();}
        TypeOfReferenceGeometry getGeometryType() const  {return geometryType_;}

        /// \brief Given a local index, return (assign to point) the corresponding node.
        virtual void            getNode(const IndexT& localIndex, PointReferenceT& node) const
                                {node = points_[localIndex];}

        virtual int             getLocalNodeIndex(int face, int node)const = 0;

        /// \brief For debugging and checkpointing: a human-readable name.
        virtual std::string     getName() const = 0;
        
        // ================================== Quadrature rules =====================================
        
        /// \brief Add a quadrature rule into the list of valid quadrature rules for this geometry.
        virtual void addGaussQuadratureRule(typename QuadratureRules::GaussQuadratureRule<DIM>* const qr) = 0;
        
        /// \brief Get a valid quadrature for this geometry.
        virtual typename QuadratureRules::GaussQuadratureRule<DIM>* const getGaussQuadratureRule(int order) const = 0;

    protected:
        ReferenceGeometry(const TypeOfReferenceGeometry& geoT);
        ReferenceGeometry(unsigned int numberOfNodes, const TypeOfReferenceGeometry& geoT);
        ReferenceGeometry(const ReferenceGeometry& other);
        
    protected:
        /// Container of the actual points (no reference).
        VectorOfReferencePointsT        points_;
        /// An identifier of the type of referenceGeometry, that some say shouldn't be used.
        const TypeOfReferenceGeometry   geometryType_;
        
    };
};
#include "ReferenceGeometry_Impl.hpp"
#endif
