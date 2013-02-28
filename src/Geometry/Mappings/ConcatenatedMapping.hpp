//
//  ConcatenatedMapping.h
//  
//
//  Created by Shavarsh Nurijanyan on 2/13/13.
//
//

#ifndef ____ConcatenatedMapping__
#define ____ConcatenatedMapping__


#include "../PointReference.hpp"
#include "../ReferencePoint.hpp"
#include "../PhysicalGeometry.hpp"
#include "MappingReferenceToReference.hpp"
#include "MappingReferenceToPhysical.hpp"
#include "../Jacobian.hpp"
    //using Geometry::multiplyJacobiansInto;
    //------------------------------------------------------------------------------

namespace Geometry
{
    /*! ConcatenatedMapping allows to compose two mappings of type
     *  Ref2RefSpaceMapping and treat the result as a mapping itself. This
     *  functionality is needed for the face-to-reference element
     *  mappings. There, for the mapping to the right (R) element, the
     *  face-to-face mapping (from left to right side) and the right side's
     *  face-to-element mapping have to be combined. The resulting mapping can
     *  be handed out by the face, and that's where ConcatenatedMapping comes
     *  in. Note that it is necessary to allow different dimensions, since a
     *  face-to-element mapping is always (dim-1) -> dim.
     *
     *  Internally, ConcatenatedMapping keeps references only, so that it is
     *  clear that this class will not delete anything. Hence the lifetime of
     *  the ConcatenatedMapping must not exceed the one of the component
     *  mappings, but this is hardly possible for reasons in the mathematical
     *  usage. */
    template <unsigned int dFrom, unsigned int dIntermediate, unsigned int dTo = dIntermediate>
    class ConcatenatedMapping: public MappingReferenceToReference<dFrom, dTo>
    {
       public:
            //! Ctor gets two references to existing mappings.
        ConcatenatedMapping(const MappingReferenceToReference<dFrom, dIntermediate>& m1,
                            const MappingReferenceToReference<dIntermediate, dTo>& m2):
            map1_(m1),
            map2_(m2)
	    {
        }
        
            //! Transformation is simply via the intermediate space.
        virtual void transform(const PointReference<dFrom>& pIn, PointReference<dTo>& pOut) const
	    {
            PointReference<dIntermediate> pLoc;
            
            map1_.transform(pIn, pLoc);
            map2_.transform(pLoc, pOut);
	    }
        
            //! To compute the Jacobian, the two component ones have to multiplied.
        virtual void calcJacobian(const PointReference<dFrom>& p, Jacobian<dFrom, dTo>& jac) const
	    {
            PointReference<dIntermediate> pIntermediate;
            map1_.transform(p, pIntermediate);
            
                
            Jacobian<dFrom, dIntermediate> j1;
            Jacobian<dIntermediate, dTo> j2;
            
            map1_.calcJacobian(p, j1);
            map2_.calcJacobian(pIntermediate, j2);
            
            j2.multiplyJacobiansInto(j1, jac);
	    }
        
    private:
        const MappingReferenceToReference<dFrom, dIntermediate>&    map1_;
        const MappingReferenceToReference<dIntermediate, dTo>&      map2_;
    };
} // close namespace Geometry


#endif /* defined(____ConcatenatedMapping__) */
