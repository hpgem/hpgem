/*
 * MappingToPhysSimplexLinear.hpp
 *
 *  Created on: Feb 5, 2013
 *      Author: nicorivas
 */

#ifndef MAPPINGSIMPLEXLINEAR_CPP_
#define MAPPINGSIMPLEXLINEAR_CPP_

/*! The mapping is linear in every reference space dimension, with the offset to
 *  the origin of the dim-simplex. */
template<unsigned int DIM>
void MappingToPhysSimplexLinear<DIM>::
transform(const PointReferenceT& pointReference, PointPhysicalT& pointPhysical) const
{
    pointPhysical = a[0];
    for (unsigned int i = 1; i <= DIM; ++i) pointPhysical += pointReference[i-1] * a[i];
}

/*! The Jacobian results from the transform function by symbolic derivation. */
template<unsigned int DIM>
void MappingToPhysSimplexLinear<DIM>::
calcJacobian(const PointReferenceT& pointReference, JacobianT& jacobian) const
{
    for (unsigned int i = 0; i < DIM; i++)
    for (unsigned int j = 0; j < DIM; j++)
    {
        jacobian(i,j)=a[j+1][i];
    }
}

/*! For the simplices it actually saves computations to save the coefficients.*/
template<unsigned int DIM>
void MappingToPhysSimplexLinear<DIM>::reinit(const PhysicalGeometry*const physicalGeometry)
{
    PointPhysicalT p0(DIM),pi(DIM);
    physicalGeometry->getNodeCoordinates(0, p0);
    a[0]=p0;
    for (unsigned int i=1; i <= DIM; ++i)
    {
        physicalGeometry->getNodeCoordinates((int) i, pi);
        a[i] = pi-p0;
    }
}
#endif /* MAPPINGSIMPLECUBENLINEAR_H_ */
