/*
This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
http://hpgem.org/about-the-code/team

This code is distributed using BSD 3-Clause License. A copy of which can found below.


Copyright (c) 2018, Univesity of Twenete
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef ALGORITHMS_DIVDGMAXDISCRETIZATION_h
#define ALGORITHMS_DIVDGMAXDISCRETIZATION_h

#include "../BaseExtended.h"

/// \brief Discontinuous Galerkin discretization for Maxwell, where the
/// divergence constraint (div E = 0) is part of the discretization.
///
// TODO: Write a better introduction on the discretization_ used
// Should include that we can split the problem matrices into four parts A, B,
// C and M. Which for example combine to form a eigen value problem
// [ A  B ] [u] = λ [ M 0 ] [u]
// [ BT C ] [p]     [ 0 0 ] [p]
class DivDGMaxDiscretization
{
public:
    // TODO: static const std::size_t matrix/vector ids

    static const std::size_t ELEMENT_MASS_MATRIX_ID = 0;
    static const std::size_t ELEMENT_STIFFNESS_MATRIX_ID = 1;
    // Note: Missing are the initial conditions, taking up positions 0, 1
    static const std::size_t ELEMENT_SOURCE_VECTOR_ID = 2;
    static const std::size_t FACE_STIFFNESS_MATRIX_ID = 0;
    static const std::size_t FACE_BOUNDARY_VECTOR_ID = 0;

    struct Stab
    {
        /// \brief Stabilization parameter for the tangential part of u/v.
        double stab1;
        /// \brief Stabilization parameter for the normal part of u/v.
        double stab2;
        /// \brief Stabilization parameter for p/q.
        double stab3;
    };

    // See notes in DGMaxDiscretization
    using InputFunction = std::function<void(const PointPhysicalT &, LinearAlgebra::SmallVector<DIM>&)>;
    using FaceInputFunction = std::function<void(const PointPhysicalT &, Base::PhysicalFace<DIM>&, LinearAlgebra::SmallVector<DIM>&)>;

    void initializeBasisFunctions(Base::MeshManipulator<DIM> &mesh);

    void computeElementIntegrands(Base::MeshManipulator<DIM>& mesh, bool invertMassMatrix,
                                  const InputFunction& sourceTerm,
                                  const InputFunction& initialCondition,
                                  const InputFunction& initialConditionDerivative) const;

    void computeFaceIntegrals(Base::MeshManipulator<DIM>& mesh, FaceInputFunction boundaryCondition,
            Stab stab) const;

    // TODO: LJ include the same norms as in DGMaxDiscretization
    double computeL2Error(Base::MeshManipulator<DIM>& mesh, std::size_t timeVector, const InputFunction& electricField) const;

    LinearAlgebra::SmallVector<DIM> computeField(
            const Base::Element* element, const Geometry::PointReference<DIM>& point,
            const LinearAlgebra::MiddleSizeVector& coefficients) const;

private:
    // Element part of matrix M, with zero matrices around it
    void elementMassMatrix(Base::PhysicalElement<DIM>& el , LinearAlgebra::MiddleSizeMatrix& ret) const;
    // Stiffness matrix, TODO: Is this the real stiffness or is this more something like curl-stiffness?
    // Element part of matrix A, with zeros around it
    void elementStiffnessMatrix(Base::PhysicalElement<DIM>& el , LinearAlgebra::MiddleSizeMatrix& ret) const;
    // Element part of matrix B and B^T, with zeros around it
    void elementScalarVectorCoupling(Base::PhysicalElement<DIM>& el , LinearAlgebra::MiddleSizeMatrix& ret) const;

    // The source vector for harmonic problems.
    void elementSourceVector(Base::PhysicalElement<DIM>& el,
            const InputFunction &source, LinearAlgebra::MiddleSizeVector& ret) const;

    // The part [[v]]_t {{mu^{-1} curl u}} + [[u]]_t {{mu^{-1} curl v}} of the stiffness integrand.
    void faceStiffnessMatrix1(Base::PhysicalFace<DIM>& fa, LinearAlgebra::MiddleSizeMatrix& ret) const;
    // The tangential stability term stab * [[u]]_T [[v]]_T part of the stiffness integrand
    void faceStiffnessMatrix2(Base::PhysicalFace<DIM>& fa, LinearAlgebra::MiddleSizeMatrix& ret, double stab1) const;
    // The normal stability term stab [[eps u]]_N [[eps v]]_N
    void faceStiffnessMatrix3(Base::PhysicalFace<DIM>& fa, LinearAlgebra::MiddleSizeMatrix& ret, double stab2) const;
    // The face part of B and B^T, [[p]]_N {{eps v}}
    void faceScalarVectorCoupling(Base::PhysicalFace<DIM>& fa, LinearAlgebra::MiddleSizeMatrix& ret) const;
    // Matrix C, stab [[p]]_N [[q]]_n
    void faceStiffnessScalarMatrix4(Base::PhysicalFace<DIM>& fa, LinearAlgebra::MiddleSizeMatrix& ret, double stab3) const;

    void faceBoundaryVector(Base::PhysicalFace<DIM>& fa, const FaceInputFunction &boundaryValue,
                            LinearAlgebra::MiddleSizeVector &ret, double stab1) const;

    double elementErrorIntegrand(Base::PhysicalElement<DIM> &el, std::size_t timeVector,
            const InputFunction& exactValues) const;
};


#endif //ALGORITHMS_DIVDGMAXDISCRETIZATION_h
