/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2021, University of Twente
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors
 may be used to endorse or promote products derived from this software without
 specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef HPGEM_APP_MESHPREDECLARATIONS_H
#define HPGEM_APP_MESHPREDECLARATIONS_H

#include "idtypes.h"
#include "elementShape.h"
#include "ElementShapes.h"
#include "utils/tag.h"
#include "utils/TemplateArray.h"

#include "LinearAlgebra/SmallVector.h"

#include <array>
#include <limits>
#include <memory>
#include <vector>

// FILE STRUCTURE //
////////////////////
//
// The mesh is a combination of several intertwined classes. This design is
// complicated by that most of these classes are templated by at least the
// dimension of the mesh, thus preventing the natural header and implementation
// file separation.
//
// To structure these classes they have been split into several separated header
// files located in the MeshImpl folder. This header takes these smaller parts
// and combines them into a definition of the mesh for use.
//
// The structure in this file is as follows:
//  - Includes needed for both definition and implementation
//  - Predeclarations of the classes to facilitate the definitions and as hint
//    of the interface provided.
//  - Includes of definition files
//  - Includes of the implementation files

namespace Preprocessor {

template <std::size_t dimension, std::size_t gridDimension>
class MeshEntity;

template <std::size_t dimension>
class Element;

template <std::size_t dimension>
class Mesh;

}  // namespace Preprocessor

// Definitions
#include "MeshImpl/MeshEntity.h"
#include "MeshImpl/Element.h"
#include "MeshImpl/Mesh.h"

// Implementations
#include "MeshImpl/MeshEntity_Impl.h"
#include "MeshImpl/Element_Impl.h"
#include "MeshImpl/Mesh_Impl.h"

#endif  // HPGEM_MESHPREDECLARATIONS_H
