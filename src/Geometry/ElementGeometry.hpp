#ifndef ELEMENTGEOMETRY_H_
#define ELEMENTGEOMETRY_H_

#include "ReferenceGeometry.hpp"
#include "PhysicalGeometry.hpp"
#include "RefinementGeometry.hpp"
#include "Base/PhysicalSpaceFunctor.hpp"
#include "Mappings/MappingReferenceToPhysical.hpp"

             
namespace Geometry
{
    template <unsigned int DIM>
    class ElementGeometry
    {
    public:
        typedef Point<DIM>                              PointT;
        typedef PointPhysical<DIM>                      PointPhysicalT;
        typedef PointReference<DIM>                     PointReferenceT;
        typedef PhysicalGeometry<DIM>                   PhysicalGeometryT;
        typedef ReferenceGeometry<DIM>                  ReferenceGeometryT;
        typedef RefinementGeometry<DIM>                 RefinementGeometryT;
        typedef MappingReferenceToPhysical<DIM, DIM>    MappingReferenceToPhysicalT;
        typedef LinearAlgebra::Matrix                   MatrixT;
        typedef Jacobian<DIM, DIM>                      JacobianT;
        typedef unsigned int                            PointIndexT;
        typedef std::vector<PointPhysicalT>             VectorOfPhysicalPointsT;
        typedef std::vector<PointIndexT>                VectorOfPointIndexesT;
    public:
    
        /// New style constructor with one less pass
        /// \bug This will only work for cubes
        ElementGeometry(const VectorOfPointIndexesT&   globalNodeIndexes, const VectorOfPhysicalPointsT& nodes);
        
        /// Copy constructor
        ElementGeometry(const ElementGeometry& other);
        
        virtual ~ElementGeometry();

        /// Returns a pointer to the referenceToPhysicalMapping
        const MappingReferenceToPhysicalT* const    getReferenceToPhysicalMap() const;

        /// Returns a pointer to the physicalGeometry object.
        const PhysicalGeometryT* const              getPhysicalGeometry() const;
            /// Returns a pointer to the physicalGeometry object.
        unsigned int                                getNrOfNodes() const;
        /// Returns a pointer to the referenceGeometry object.
        const ReferenceGeometryT* const             getReferenceGeometry() const;
        /// Returns a pointer to the refinementGeometry object.
        const RefinementGeometryT*                  getRefinementGeometry() const;
        /// This method gets a PointReference, which specifies a coordinate in the ReferenceGeometry,
        /// and returns a PointPhysical which is the corresponding point in the PhysicalGeometry,
        /// given the mapping.
        void                                        referenceToPhysical(const PointReferenceT& pointReference, PointPhysicalT& pointPhysical);
        
        void                                        referenceToPhysical(const PointReferenceT& pointReference, PointPhysicalT& pointPhysical)const;
        
        /// This method gets a PointReference and returns the corresponding jacobian of the
        /// referenceToPhysicalMapping.
        void                                        calcJacobian(const PointReferenceT& pointReference, JacobianT& jacobian) const;

        /*! \brief Transform a physical space function on the reference element.

            To be able to query a function that takes physical space coordinates
            at reference space points we wrap it with a class that offers the
            necessary operator(), which takes care of the necessary
            transformations. */
        
        template <class FType>
        Base::PhysicalSpaceFunctor<DIM, FType> transformToReferenceElement(const FType& f) const
        {
            return Base::PhysicalSpaceFunctor<DIM, FType>(this, f);
        }
        

        void                                        enableRefinement();
        
        
    public:
            /// Output operator.
        friend ostream& operator<<(ostream& os, const ElementGeometry& elementGeometry)
        {
            os << "PhysicalGeometry={";
            
            for (int i = 0; i < elementGeometry.physicalGeometry_->getNumberOfNodes(); i++)
            {   
                os << (elementGeometry.physicalGeometry_)->getNodeIndex(i) << " ";
            }
            os << '}' << std::endl;
            return os;
        }
    private:
        
        static const ReferenceGeometryT* const          createReferenceGeometry(unsigned int size);
        
        static const PhysicalGeometryT* const           createPhysicalGeometry(const VectorOfPointIndexesT&      globalNodeIndexes,
                                                                               const VectorOfPhysicalPointsT&    nodes,
                                                                               const ReferenceGeometryT* const   geo);
        
        static const MappingReferenceToPhysicalT* const createMappings(unsigned int size, const PhysicalGeometryT* const pGeo);
        
    
    protected:
            /// The corresponding referenceGeometry object, for integration.
        const ReferenceGeometryT* const              referenceGeometry_;
        
            /// The physicalGeometry object contains pointers to the actual physical points, and
            /// a container of global node indexes.
        const PhysicalGeometryT* const               physicalGeometry_;

        /// The referenceToPhysicalMapping relates the coordinates of the reference object to the
        /// physical object; basically a matrix transformation.
        const MappingReferenceToPhysicalT*const      referenceToPhysicalMapping_;

        /// The corresponding refinementGeometry object
        RefinementGeometryT*                         refinementGeometry_;
    };
}
#include "ElementGeometry_Impl.hpp"
#endif
