#ifndef FILLMATRICES_HPP
#define FILLMATRICES_HPP

#include "Integration/ElementIntegrandBase.hpp"
#include "Integration/FaceIntegrandBase.hpp"
#include "Base/Face.hpp"
#include "Base/Element.hpp"

class hpGemUIExtentions;

/**
 * this class provides the matrixes for the stiffness matrix, so a choise about using the IP-method or the Brezzi method can be made on a high level
 */
class matrixFiller: protected Integration::ElementIntegrandBase<LinearAlgebra::Matrix>, protected Integration::FaceIntegrandBase<LinearAlgebra::Matrix>,
					protected Integration::ElementIntegrandBase<LinearAlgebra::NumericalVector>, protected Integration::FaceIntegrandBase<LinearAlgebra::NumericalVector>{
public:

    /// the actual global assembly routine
    virtual void fillMatrixes(hpGemUIExtentions* matrixContainer)=0;

    /**
     * integrand for the filling of the element contibutions to the stiffness matrix S
     * \param [in] element the element that is currently being integrated on
     * \param [in] p the gauss point
     * \param [out] ret the contributions to the stifness matrix from this point. This should not yet be scaled down with the weight of this point!
     */
    void elementIntegrand(const Base::Element* element, const Geometry::PointReference& p, LinearAlgebra::Matrix& ret);

    /**
     * Computes the space-dependent part of the source term. This assumes that the source term can be uncoupled in a time dependent part and a space dependent part
     * \param [in] p The point in space where the source is to be evaluated.
     * \param [out] ret The value of the source term.
     */
    void initialExactSolution(const Geometry::PointPhysical& p, NumericalVector& ret);
    virtual void boundaryConditions(const Geometry::PointPhysical &p, NumericalVector &ret){initialExactSolution(p,ret);}
    virtual void sourceTerm(const Geometry::PointPhysical &p, NumericalVector &ret){initialExactSolution(p,ret);ret*=M_PI*M_PI*8-1;};

    /**
      * integrand used for the computation of the space dependent source term
      * maps the reference point and the element to physical coordinates and then uses this to call the user provided initial conditions
      */
    void elementIntegrand(const Base::Element* element, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret);
};

/**
 * this class provides a global assembly structure that is suitable for use with the IP-method
 */
class matrixFillerIP:public matrixFiller{
	double stabCoeff_;

public:
	matrixFillerIP(double stabCoeff):stabCoeff_(stabCoeff){}

    virtual void fillMatrixes(hpGemUIExtentions* matrixContainer);

    /**
     * integrand for the filling of the penalty contibutions to the stiffness matrix S, for the Interior Penalty method.
     * The penalty parameter will be multiplied in later
     * \param [in] element the element that is currently being integrated on
     * \param [in] p the gauss point
     * \param [out] ret the contributions to the stifness matrix from this point. This should not yet be scaled down with the weight of this point!
     * For internal faces the integration expects that this matrix contains first contributions associated with the left element and then with the right element
     */
    void faceIntegrand(const Base::Face* face, const NumericalVector& normal, const Geometry::PointReference& p, LinearAlgebra::Matrix& ret);

    /**
     * Integrand used for the computation of the boundary contributions to the RHS. This will be scaled by the same time dependent factor as in the source therm, just like in the original code by Domokos.
     * This version is specialized for the Interior Penalty method. It also computes the terms that are common to both the Interior Penalty method and the Brezzi method because it safes some loops over elements where no work needs to be done
     */
    void faceIntegrand(const Base::Face* face, const NumericalVector& normal, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret);
};

/**
 * this class provides a global assembly stucture that is suitable for the Brezzi-flux
 */
class matrixFillerBR:public matrixFiller{
	double stabCoeff_;

public:
	matrixFillerBR(double stabCoeff):stabCoeff_(stabCoeff){}
    virtual void fillMatrixes(hpGemUIExtentions* matrixContainer);

    /**
     * integrand for the filling of a matrix of the form D_ij=integral(phi_i*(n x phi_j)), to enable computation of the penalty contributions to the stiffness matrix S
     * \param [in] element the element that is currently being integrated on
     * \param [in] p the gauss point
     * \param [out] ret the contributions to the stifness matrix from this point. This should not yet be scaled down with the weight of this point!
     * For internal faces the integration expects that this matrix contains first contributions associated with the left element and then with the right element
     */
    void faceIntegrand(const Base::Face* face, const NumericalVector& normal, const Geometry::PointReference& p, LinearAlgebra::Matrix& ret);

    /**
     * Integrand used for the computation of the boundary contributions to the RHS. This will be scaled by the same time dependent factor as in the source therm, just like in the original code by Domokos.
     * This version is used in the computation of D_i=integral(phi_i*(n x u_0))
     */
    void faceIntegrand(const Base::Face* face, const NumericalVector& normal, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret);

};

#endif
