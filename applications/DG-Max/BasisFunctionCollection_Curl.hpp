//this file contains the basis functions and some special vector manipulations that are not needed for more sensible bases

#ifndef BasisFunctionsCollection_Curl_hpp
#define BasisFunctionsCollection_Curl_hpp

#include "Base/BaseBasisFunction.hpp"
#include "Base/MeshManipulator.hpp"

//should probably be moved to the utilities namespace and another file at some point
/**
 * Computes the 3D outer product of 2 vectors
 * \param [in] a,b The 2 vectors that you want the outer product of.
 * \param [out] ret  a x b
 */
void OuterProduct(const NumericalVector& a,const NumericalVector& b,NumericalVector& ret);

namespace Base
{    
    /**
     * Computes legendre polynomials of arbitrary degree using a recursive defenition
     * \param [in] degree The degree of the legendre polynomial.
     * \param [in] x The point where you want the function evaluated.
     * \return The legendre polynomial of the requested degree, evaluated at x.
     * \pre -1<=x<=1
     * \pre 0<=degree
     */
    double LegendrePolynomial(int degree, double x);
    
    /**
     * Computes the derivative legendre polynomials of arbitrary degree using a recursive defenition
     * \param [in] degree The degree of the legendre polynomial. This means that the degree of the derivative will be one lower than this number.
     * \param [in] x The point where you want the function evaluated.
     * \return Tthe legendre polynomial of the requested degree, evaluated at x.
     * \pre -1<=x<=1
     * \pre 0<=degree
     */
     double LegendrePolynomialDerivative(int degree, double x);
    
    /**
     * Linear scalar nodal basis functions, with nodes on the element vertices. They are used to give a consistent 
     * coordinate system for the elements.
     */
    struct Basis_Curl_Bari : public BaseBasisFunction
    {
        int VertexNr;
        Basis_Curl_Bari(int vertex);
        
        //pure abstract functions derived from superclass
        virtual double eval(const PointReferenceT& p) const;      
        virtual double evalDeriv0(const PointReferenceT& p) const;   
        virtual double evalDeriv1(const PointReferenceT& p) const;  
        virtual double evalDeriv2(const PointReferenceT& p) const;  
        
        /**
	 * Combines \ref evalDeriv0 , \ref evalDeriv1 and \ref evalDeriv2 into one function.
	 * \param [in] p The point where the derivative should be evaluated.
	 * \param [out] ret a length 3 vector equal to {\ref evalDeriv0(p),\ref evalDeriv1(p),\ref evalDeriv2(p)}
	 */
        virtual void evalDeriv(const PointReferenceT& p, NumericalVector& ret);
        static Basis_Curl_Bari barycentricFunctions[];
    };
    

    /**
     * Base class for 3D vector valued functions
     */
    struct threeDBasisFunction : public BaseBasisFunction{
      
        //this class has to be derived from basebasisfunctions but we have no use for a scalar eval or the evaluation of a gradient
	/// Function not supported
	virtual double eval(const PointReferenceT& p) const;    
	/// Function not supported
	virtual double evalDeriv0(const PointReferenceT& p) const;  
	/// Function not supported
	virtual double evalDeriv1(const PointReferenceT& p) const;  
	/// Function not supported
	virtual double evalDeriv2(const PointReferenceT& p) const;  
	
	/**
	 * gives a physical point where the function 'is the most active'. A good implementation tries to make sure that
	 * integral(||exp(ik*x)phi||^2) is approximately integral(||exp(ik*node)phi||^2) for all k,
	 * but usually a decent guess is returned
	 */
	virtual void getReasonableNode(const Element& element, Geometry::PointPhysical node)=0;
	
	/**
	 * Computes the value of the basis function.
	 * \param [in] p The point where the derivative should be evaluated.
	 * \param [out] ret The value of the basis function
	 */
	virtual void eval(const PointReferenceT& p, NumericalVector& ret) const=0;
	
	/**
	 * Computes the curl of the basis functions. Implementations should prefer analytic derivation over numerical.
	 * \param [in] p The point where the curl should be evaluated
	 * \param [out] ret The curl of the basis function
	 */
	virtual void evalCurl(const PointReferenceT& p, NumericalVector& ret) const=0;
    };
    
    //! Curl conforming edge functions.
    struct Basis_Curl_Edge : public threeDBasisFunction
    {
        const int deg,o,i;
        Basis_Curl_Edge(int degree, int localFirstVertex, int localSecondVertex);
        
	virtual void eval(const PointReferenceT& p, NumericalVector& ret) const;

        virtual void evalCurl(const PointReferenceT& p,NumericalVector& ret) const;
	
	virtual void getReasonableNode(const Element& element, Geometry::PointPhysical node);
    };
    
    //! Curl conforming edge based face functions.
    struct Basis_Curl_Edge_Face : public threeDBasisFunction
    {
        int deg,a,b,c;
        Basis_Curl_Edge_Face(int degree, int localOpposingVertex, int localSpecialVertex);
        
        virtual void eval(const PointReferenceT& p, NumericalVector& ret) const;
        
        virtual void evalCurl(const PointReferenceT& p, NumericalVector& ret) const;
	
	virtual void getReasonableNode(const Element& element, Geometry::PointPhysical node);
    };
    
    //! Curl conforming face functions.
    struct Basis_Curl_Face : public threeDBasisFunction
    {
        int deg1,deg2,a,b,c;
        Basis_Curl_Face(int degree1,int degree2, int localOpposingVertex, int direction);
        
        virtual void eval(const PointReferenceT& p, NumericalVector& ret) const;

        virtual void evalCurl(const PointReferenceT& p, NumericalVector& ret) const;
	
	virtual void getReasonableNode(const Element& element, Geometry::PointPhysical node);
    };
    
    //! Curl conforming face based interior functions.
    struct Basis_Curl_Face_interior : public threeDBasisFunction
    {
        int deg1,deg2,a,b,c,d;
        Basis_Curl_Face_interior(int degree1,int degree2, int localOpposingVertex);
        
        virtual void eval(const PointReferenceT& p, NumericalVector& ret) const;

        virtual void evalCurl(const PointReferenceT& p, NumericalVector& ret) const;
	
	virtual void getReasonableNode(const Element& element, Geometry::PointPhysical node);
    };
    
    //! curl conforming interior functions
    struct Basis_Curl_interior : public threeDBasisFunction
    {
        const int deg1,deg2,deg3,direction;
        Basis_Curl_interior(int degree1,int degree2,int degree3, int direction);
                
        virtual void eval(const PointReferenceT& p, NumericalVector& ret) const;

        virtual void evalCurl(const PointReferenceT& p, NumericalVector& ret) const;
	
	virtual void getReasonableNode(const Element& element, Geometry::PointPhysical node);
    };
};

/**
 * A small extension to allow overriding \ref createBasisFunctions. (Will probably be redundand in the release version of hpGEM2)
 */
class MyMeshManipulator : public Base::MeshManipulator  {
    
public:
    MyMeshManipulator (const Base::ConfigurationData* data , int order, bool xPer = 0, bool yPer = 0, bool zPer = 0);
  
private:
    void createBasisFunctions(unsigned int order);
};

#endif  // BasisFunctionsCollection_Curl_hpp
