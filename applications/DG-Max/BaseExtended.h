/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, University of Twente
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef BASEEXTENDED_HPP
#define BASEEXTENDED_HPP

//This file contains all the fiddly bits about interpolation and filling matrices and stuff

#include "Base/MeshManipulator.h"
#include "Base/HpgemAPIBase.h"
#include "DGMaxDim.h"

//#include <Base/Norm2.h>
#include <fstream>
#include <sstream>
#include <ctime>
#include "ElementInfos.h"
#include "Base/Element.h"


using PointPhysicalT = Geometry::PointPhysical<DIM>;
using PointElementReferenceT = Geometry::PointReference<DIM>;
using PointFaceReferenceT = Geometry::PointReference<DIM-1>;

using ElementT = Base::Element;
using FaceT = Base::Face;

/**
 * This class provides some significant extentions to HpgemUI that may also be useful for other problems
 */
class hpGemUIExtentions : public Base::HpgemAPIBase<DIM>
{
    
public:

    //! gives read-only acces to the configuration data
    const Base::ConfigurationData* getConfigData();

    /**
     * Constructor allows PETSc to parse any PETSc-related input and does most of the initialisations.
     *
     * This constructor is extremely bare bones: it only guarantees PETSc dataTypes exists, not that they are
     * actually usefull for computations or storing data.
     */
    hpGemUIExtentions(Base::GlobalData * const globalConfig, Base::ConfigurationData* elementConfig);
    
    /**
     * Deconstructor cleans up again and logs performance statistics of PETSc
     */
    ~hpGemUIExtentions() override;
    
    /**
     * Wrapper for protected function in superclass
     */
    std::size_t addMesh(Base::MeshManipulator<DIM>* mesh);

    Base::MeshManipulator<DIM>* getMesh(std::size_t meshId);


// Removed because it is not used by anybody, porting it to the new orginazation
// of a single-class-per-solver would take some investment in time both for the
// refactoring and for the understanding of the code.
//    /*
//     * Call after setting up the mesh. This will arrange that the matrices are assembled and then find the Density of States using eigenfunction expansions
//     */
//    void solveDOS();

// See comment on solveDOS()
//    /**
//     * given an eigenvector, prepare an expression for f(x) in int(f(x)*delta(omega) dx)
//     */
//    void makeFunctionValue(Vec eigenVector, LinearAlgebra::MiddleSizeVector& result);
// See comment on solveDOS()
//    void LDOSIntegrand(Base::PhysicalElement<DIM>& el, double& ret);
};

#endif
