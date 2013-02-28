//----------------------------------------------------------------
#ifndef ElementData_hpp
#define ElementData_hpp
//----------------------------------------------------------------
#include <vector>
#include "Base/Element.hpp"

namespace Base
{
    template <unsigned int DIM>
	class ElementData
	{
        /*!
         * ElementData contains the data for every element in the vector elementData_.
         * This is a two dimensional vector where the first dimensions is the time level,
         * and the second the number of unknowns times the number of basis functions.
         */
	    public:
            typedef double NumberT;

	    public:
// MTJ: default constructor, just for a moment; to be removed later!  Soon!
ElementData() {}
	      
            ElementData(unsigned int timeLevels, unsigned nrOfUnkowns, unsigned int nrOfBasisFunctions);

            /// Specify a time level index, return a vector containing the data for that time level.
            std::vector<NumberT> getTimeLevelData(unsigned int timeLevel);

            /// Specify a time level index, and a data index, return the corresponding data (double)
            double getData(unsigned int timeLevel, unsigned int index);

            /// Specify a time level index, an unknown and a basis function nr, return data (double)
            double getData(unsigned int timeLevel, unsigned int unknown, unsigned int basisFunction);

            virtual ~ElementData() {;}
            
	    private:
		    unsigned int timeLevels_;
		    unsigned int nrOfUnkowns_;
		    unsigned int nrOfBasisFunctions_;
		    std::vector<std::vector<NumberT> > elementData_;
	};
};
#include "ElementData_Impl.hpp"
#endif