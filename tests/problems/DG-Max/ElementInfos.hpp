//this file has a container data structure for everything you want to know on a per element basis

#ifndef Elementinfos_hpp
#define Elementinfos_hpp

#include "Base/UserElementData.hpp"
#include "Base/Element.hpp"
#include "Base/Face.hpp"
#include "BasisFunctionCollection_Curl.hpp"
typedef Geometry::PointReference PointElementReferenceT;
#if __cplusplus<=199711L
#include "map"
typedef std::map<PointElementReferenceT,std::vector<NumericalVector> > myMap;
#else
#include "unordered_map"
typedef std::unordered_map<PointElementReferenceT,std::vector<NumericalVector> > myMap;//the unordered_map trades functionality for speed
#endif

typedef Base::threeDBasisFunction basisFunctionT;
typedef Geometry::PointReference PointFaceReferenceT;

//one there is a default way in hpGEM to configure code without haveing to recompile this and some other things should be grouped in another file
/**
 * Stores some parameters that should be available everywhere. 
 * Attributes need not be constant since this struct will be stored in a constant field.
 */
struct MaxwellData: public Base::GlobalData{
    double Sigma_;//allows for lossy media - untested except maybe by Domokos
    double StabCoeff_;//a_F for IP or 1+eta_F/4 for BR
    int NumberOfIntervals_;
    int PolynomialOrder_;
    double StartTime_;
    double EndTime_;
    
    MaxwellData(int numberOfIntervals, int polynomialOrder);
};

/**
 * Computes the transpose of the inverse of a 3 by 3 Jacobian.
 * \param [in] orig The matrix to be inverted
 * \param [out] inverse The inverse of the transpose of orig. This may not be the same matrix as orig.
 */
void InvertAndTranspose(Geometry::Jacobian& orig, Geometry::Jacobian& inverse);

/**
 * store the basisfunction values for the reference element (i.e. without any transformations) for later reuse
 * this class will also compute basisfunction values if they are not yet available
 */
class FunctionCache{
  
    static myMap valueCache_;
    static myMap curlCache_;
  
public:
    
    /**
     * gets the function values of the function on the reference element
     */
    static void getFunctionValuesVector(const Base::Element* element, const PointElementReferenceT& point, std::vector<NumericalVector>& values);
    
    /**
     * gets the curl of the function on the reference element
     */
    static void getFunctionCurlsVector(const Base::Element* element, const PointElementReferenceT& point, std::vector<NumericalVector>& curls);
};

/**
 * store some usefull information that needs to be computed everytime at the beginning of an integrand
 * specialized for tetrahedra
 */
class ElementInfos:public UserElementData{
public:
    double determinant_;
    Geometry::Jacobian Jacobian_,inverse_;
    double epsilon_;
    
    ElementInfos(const Base::Element& element);
    
    /**
     * gets the basisfunction values on one mantissa from the reference element and transforms them based on the current element
     */
    void makeFunctionValuesVector(const Base::Element* element, const PointElementReferenceT& point, std::vector<NumericalVector>& values);
    
    /**
     * gets the curls of the baisfunction values on one mantissa from the reference element and transforms them based on the current element
     */
    void makeFunctionCurlsVector(const Base::Element* element, const PointElementReferenceT& point, std::vector<NumericalVector>& curls);
};


#endif  // Elementinfos_hpp
