<Fancy ascII logo goes here>

# hpGEM

hpGEM is a c++ partial differential equation solver. It is intended for Discontinuous Galerkin Finite Element Methods, but can also do normal (conforming) Finite Element Methods and Finite Volume.

## Installation

For a quick, basic installation do

> mkdir ../hpGEM_build
> cd ../hpGEM_build
> cmake ../trunk
> make

For detailled instructions and alternative options see http://www.hpgem.org/obtaining-the-code/downloading-and-installing-version-2.x

## Use

Doing a finite element simulation is a two step process. The first step is to create or convert a mesh to the native hpgem format and then use this mesh as input for your finite element simulation. This section will describe the preprocessor in more detail. For additional information on how to do the simulation please have a look at your preferred solver or the tutorial applications. At the moment the preprocessor can process centaur meshes and meshes in the native hpgem format. In addition it can mesh a rectangular domain based on a parametrized description of what the mesh should look like. We will first describe in words what this mesh description looks like and then provide some examples as illustration. 

The first line of the description file contains the text "structured 1" (without quotes) followed by the dimension d of the domain. The next line contains d integers denoting the number of elements per dimension. The next two lines contain d floating point numbers denoting the lower left and upper right corners of the domain respectively. The next line contains either the word "rectangular" or the word "triangular" (both without quotes) the denote if you want a rectanular or a triangular mesh. A triangular mesh is formed by first generating a rectangular mesh and then splitting each element in either 2 triangles or 5 tetrahedra in an alternating fasion, depending on the dimension. It is not allowed to ask for a triangular mesh in 1 dimension. The last line contains d booleans with are true if the domain is periodic in that dimension or false if it is not. For periodic domains the simulation domain is assumed to cover 1 periodic repetitions of the actual domain. If you create a periodic triangular mesh, the directions that are periodic must contain an even number of elements or else the triangles can't alternate properly.

For example if the input file contains the following:

> structured 1 3
> 4 5 6
> 0. 0. 0.
> 1. 2. 3.5
> rectangular
> false false false

describes a mesh of 4x5x6 hexahedral elements on the nonperiodic domain [0, 1]x[0, 2]x[0, 3.5]

> structured 1 2
> 7 8
> 0. -3.14159265359
> 1. 3.14159265359
> triangular
> false true

describes a mesh of 7x8 triangular elements on the boundary of a cylinder with height 1 and radius 1. Note that it is ok to have 7 elements in the height-wise direction, but we must have an even number of elements in the circular direction to make the mesh fit properly.

## Licence

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
 
## Contact
 
For questions, comments, requests and bug reports, send an e-mail to support@hpgem.org
 
## More information
 
More information is available on the website: http://www.hpgem.org

## FAQ

Q: How do I use hpGEM?
A: There are a few tutorials provided in the applications folder.

Q: Where is the documentation?
A: There is documentation on the website. You can also locally build the documentation if you enable hpGEM_BUILD_DOCUMENTATION in cmake. The command for doing so is 'make doc'. It then becomes available in docs/html/index.html.

Q: How do I create my own application?
A: Copy one of the existing applications, or copy TemplateBasicProblem if you want to start from an empty application. Then edit the copied CMakeLists.txt and change the two occurrences of 'ExampleProblem.out' to the desired name of your application

Q: I added an application, but it isn't detected!
A: CMake cannot automatically detect new applications. Please use generate in cmake-gui or ccmake or type cmake again.

Q: hpGEM is slow!
A: hpGEM is used by users of highly varying background and with highly varying programming skills. For this reason the (default) Debug configuration is extra paranoid about making sure all functions are used as intended. This safeguard drastically slows down computation. Setting CMAKE_BUILD_TYPE to Release disables these safeguards and enables compiler optimisations. It is recommended that you validate that your application is working as intended before doing this.
