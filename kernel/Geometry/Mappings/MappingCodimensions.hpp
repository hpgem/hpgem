/* MappingCodimensions.hpp
*
*  Created on: Feb 13, 2013
*      Author: nicorivas
*/

#ifndef MAPPINGCODIMENSIONS_HPP_
#define MAPPINGCODIMENSIONS_HPP_

//#include "MappingReferenceToReference.hpp"
#include <vector>

namespace Geometry
{
    class ReferenceGeometry;

    class MappingReferenceToReference;


    /*!
     *  TODO: Comment 3D and 4D cases.
     *
     *  ~OC~ CodimMaps can be used as base class for ReferenceGeometry
     *  implementations. It contains the methods necessary to organize the
     *  topological and geometrical information regarding lower-dimensional
     *  components, as well as the mappings of the reference geometry onto
     *  itself. The former items concern mappings of faces and edges onto the
     *  element, the latter is needed for the face-to-face mappings. The class
     *  is specialized for each dimension, mainly because a "recursive"
     *  derivation also becomes very hard to manage (relicts of that may still
     *  be found in old cvs versions).  */

    /*template <unsigned int dim>
    class MappingCodimensions;

    template <>
    class MappingCodimensions<0> // Only for ReferencePoint.
    {
        public:

        /// Returns 0.
        virtual int
        getCodim0MappingIndex(const std::vector<unsigned int>&, const std::vector<unsigned int>&)
        const = 0;

        /// Returns 0.
        virtual const MappingReferenceToReference<0, 0>* getCodim0MappingPtr(const unsigned int) const = 0;

        virtual ~MappingCodimensions() {}
    };


    template <>
    class MappingCodimensions<1> // Only for ReferenceLine.
    {
        public:

        //! Given two sequences of LOCAL node indices of lines, then determine the index of the
        //! correct MappingLineToLine.
        //! (I think the comment below is wrong).
        //!
        //! ~OC~
        //! Given two sequences of the global node indices of a line determine the index of the
        //! correct MappingLineToLine.
        virtual int
        getCodim0MappingIndex(const std::vector<unsigned int>&, const std::vector<unsigned int>&)
        const = 0;

        //! Given an index to a ReferenceSquare2ReferenceSquareMapping return it.
        //! ~OC~ If a self-mapping index has been checked previously, access directly
        virtual const MappingReferenceToReference<1, 1>*
        getCodim0MappingPtr(const unsigned int) const = 0;

        //! Given the two node permutations, return the mapping (through getCodim0MappingIndex,
        //! see up)
        virtual const MappingReferenceToReference<1, 1>*
        getCodim0MappingPtr(const std::vector<unsigned int>& n1, const std::vector<unsigned int>& n2) const
	    {
            return getCodim0MappingPtr(getCodim0MappingIndex(n1, n2));
	    }

        //! Get number of 'faces' (which are nodes), that is, 2.
        virtual unsigned int getNrOfCodim1Entities() const = 0;

        //! Given a local 'face' index return the local node numbers on that face.
        virtual void getCodim1EntityLocalIndices(
                const unsigned int,
                std::vector<unsigned int>&) const = 0;

        //! Given a local face index return the mapping from ReferenceLine to ReferenceSquare.
        //! ~OC~ In the 1d case we need the points to be integrated over as 'faces'!
        virtual const MappingReferenceToReference<0, 1>* getCodim1MappingPtr(const unsigned int) const = 0;

        //! Return the reference geometry of the corresponding codimension 1 entity
        //! ~OC~ In the 1d case we need the points to be integrated over as 'faces'!
        virtual const ReferenceGeometry<0>* getCodim1ReferenceGeometry(const unsigned int) const = 0;

     	virtual ~MappingCodimensions() {}
    };


    template <>
    class MappingCodimensions<2> // ReferenceSquare, ReferenceTriangle.
    {
        public:

        //! Given two sequences of LOCAL node indices of squares, the determine the index of the
        //! correct MappingSquareToSquare.
        //! (I think the comment below is wrong).
        //!
        //! ~OC~
        //! Given two sequences of the global node indices of a square determine the index of the
        //! correct MappingReferenceSquareToReferenceSquare.
        virtual int
        getCodim0MappingIndex(const std::vector<unsigned int>&, const std::vector<unsigned int>&)
        const = 0;

        //! Given an index to a ReferenceSquare2ReferenceSquareMapping return it.
        //! ~OC~ If a self-mapping index has been checked previously, access directly
        virtual const MappingReferenceToReference<2, 2>* getCodim0MappingPtr(const unsigned int) const = 0;

        //! Given the two node permutations, return the corresponding mapping pointer
        //! (through getCodim0MappingIndex, see up)
        virtual const MappingReferenceToReference<2, 2>*
        getCodim0MappingPtr(const std::vector<unsigned int>& n1, const std::vector<unsigned int>& n2) const
	    {
            return getCodim0MappingPtr(getCodim0MappingIndex(n1, n2));
	    }

        //! Get number of faces.
        virtual unsigned int getNrOfCodim1Entities() const = 0;

        //! Given a local face index return the local node numbers on that face.
        virtual void
        getCodim1EntityLocalIndices(const unsigned int, std::vector<unsigned int>&) const = 0;

        //! Given a local face index return the mapping from ReferenceLine to Reference<Shape>.
        virtual const MappingReferenceToReference<1, 2>* getCodim1MappingPtr(const unsigned int) const = 0;

        //! Return the reference geometry of the corresponding codimension 1 entity
        //! (I think in 2D always ReferenceLine).
        virtual const ReferenceGeometry<1>* getCodim1ReferenceGeometry(const unsigned int) const = 0;

        //! Return the number of vertex.
        virtual int getNrOfCodim2Entities() const = 0;

        //! Return vertex number (identity).
        virtual void getCodim2EntityLocalIndices(const int, std::vector<unsigned int>&) const = 0;

     	virtual ~MappingCodimensions() {}
    };

    template <>
    class MappingCodimensions<3> // ReferenceCube, ReferenceTetrahedron, ReferencePyramid, ReferenceTriangularPrism
    {
        public:

        virtual int
        getCodim0MappingIndex(const std::vector<unsigned int>&, const std::vector<unsigned int>&)
        const = 0;

        virtual const MappingReferenceToReference<3, 3>* getCodim0MappingPtr(const unsigned int) const = 0;

        virtual const MappingReferenceToReference<3, 3>*
        getCodim0MappingPtr(const std::vector<unsigned int>& n1, const std::vector<unsigned int>& n2) const
	    {
            return getCodim0MappingPtr(getCodim0MappingIndex(n1, n2));
	    }

        virtual unsigned int getNrOfCodim1Entities() const = 0;

        virtual void getCodim1EntityLocalIndices(const unsigned int, std::vector<unsigned int>&) const = 0;

        virtual const MappingReferenceToReference<2, 3>* getCodim1MappingPtr(const unsigned int) const = 0;

        virtual const ReferenceGeometry<2>* getCodim1ReferenceGeometry(const unsigned int) const = 0;

        virtual unsigned int getNrOfCodim2Entities() const = 0;

        virtual void getCodim2EntityLocalIndices(const unsigned int, std::vector<unsigned int>&) const = 0;

        virtual const MappingReferenceToReference<1, 3>* getCodim2MappingPtr(const unsigned int) const = 0;

        virtual const ReferenceGeometry<1>* getCodim2ReferenceGeometry(const unsigned int) const = 0;

        virtual unsigned int getNrOfCodim3Entities() const = 0;

        virtual void getCodim3EntityLocalIndices(const unsigned int, std::vector<unsigned int>&) const = 0;

    	virtual ~MappingCodimensions() {}
    };*/

    class MappingCodimensions
    {
        public:

        virtual int
        getCodim0MappingIndex(const std::vector<unsigned int>&, const std::vector<unsigned int>&)
        const = 0;

        virtual const MappingReferenceToReference* getCodim0MappingPtr(const unsigned int) const = 0;

        virtual const MappingReferenceToReference* getCodim0MappingPtr(
                const std::vector<unsigned int>& n1,
                const std::vector<unsigned int>& n2) const
	    {
            return getCodim0MappingPtr(getCodim0MappingIndex(n1, n2));
	    }

        virtual unsigned int getNrOfCodim1Entities() const {return 0;}

        virtual void getCodim1EntityLocalIndices(const unsigned int, std::vector<unsigned int>&) const {throw "The DIMension of this entity is too low to warrant maps of this codimension";}

        virtual const MappingReferenceToReference* getCodim1MappingPtr(const unsigned int) const {throw "The DIMension of this entity is too low to warrant maps of this codimension";}

        virtual const ReferenceGeometry* getCodim1ReferenceGeometry(const unsigned int) const {throw "The DIMension of this entity is too low to warrant maps of this codimension";}

        virtual unsigned int getNrOfCodim2Entities() const {return 0;}

        virtual void getCodim2EntityLocalIndices(const unsigned int, std::vector<unsigned int>&) const {throw "The DIMension of this entity is too low to warrant maps of this codimension";}

        virtual const MappingReferenceToReference* getCodim2MappingPtr(const unsigned int) const {throw "The DIMension of this entity is too low to warrant maps of this codimension";}

        virtual const ReferenceGeometry* getCodim2ReferenceGeometry(const unsigned int) const {throw "The DIMension of this entity is too low to warrant maps of this codimension";}

        virtual unsigned int getNrOfCodim3Entities() const {return 0;}

        virtual void getCodim3EntityLocalIndices(const unsigned int, std::vector<unsigned int>&) const {throw "The DIMension of this entity is too low to warrant maps of this codimension";}

    	virtual ~MappingCodimensions() {}
    };
}
#endif
