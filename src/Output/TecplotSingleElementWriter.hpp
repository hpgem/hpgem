#ifndef TECPLOTSINGLEELEMENTWRITER_HPP
#define TECPLOTSINGLEELEMENTWRITER_HPP

//forward declarations
namespace Base{
class Element;
}
namespace Geometry{
class PointReference;
}

namespace Output{

	/**
	 * If you want to write your data to tecplot it is likely you already have the function
	 * writeToTecplotFile(const Base::Element*, const Geometry::PointReference&, ostream&)
	 * already implemented in some class already, so that class can simply inherit from this class
	 * to signal the TecplotDiscontinuousSolutionWriter that it does so
	 */
	class TecplotSingleElementWriter{
	public:
		///prints the desired information for the given point to the ostream
		virtual void writeToTecplotFile(const Base::Element*, const Geometry::PointReference&, std::ostream&)=0;
	};
}

#endif
