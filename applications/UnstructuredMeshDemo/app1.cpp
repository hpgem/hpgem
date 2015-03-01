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

#include "Base/HpgemUISimplified.hpp"
#include "Base/RectangularMeshDescriptor.hpp"
#include "Output/TecplotDiscontinuousSolutionWriter.hpp"
#include "Output/TecplotPhysicalGeometryIterator.hpp"
#include "LinearAlgebra/NumericalVector.hpp"
#include "Base/PhysGradientOfBasisFunction.hpp"
#include "Base/Norm2.hpp"
#include "Output/TecplotSingleElementWriter.hpp"
#include "Base/ElementCacheData.hpp"
#include "Base/FaceCacheData.hpp"
#include "Base/Element.hpp"
#include "Base/ConfigurationData.hpp"
#include "Geometry/PhysicalGeometry.hpp"

#include <chrono>
#include <functional>
#include <array>

using Base::RectangularMeshDescriptor;
using Base::HpgemUISimplified;

const unsigned int DIM = 2;

//proof of concept code for the unstructured mesh generator
//the main challenge here is to handle the sliver elements properly
//the connectivity changes in the domain are not really a problem

//describe a domain with a rotating rotor
std::function<double(Geometry::PointPhysical, double) > meshDescription = { [](Geometry::PointPhysical point, double t) ->double
    {
        Geometry::PointPhysical rotatedPoint{DIM};
        rotatedPoint[0] = std::cos(t) * (point[0] - 0.5) - std::sin(t) * (point[1] - 0.5) + 0.5;
        rotatedPoint[1] = std::sin(t) * (point[0] - 0.5) + std::cos(t) * (point[1] - 0.5) + 0.5;
        return std::max(std::max(std::max(std::abs(point[0] - 0.5) - 0.5, std::abs(point[1] - 0.5) - 0.5),
                                 -std::max(std::abs(rotatedPoint[0] - 0.5) - 0.05, std::abs(rotatedPoint[1] - 0.5) - 0.5)),
                        -std::max(std::abs(rotatedPoint[0] - 0.5) - 0.5, std::abs(rotatedPoint[1] - 0.5) - 0.05));
    } };

//suitable for a boundary layer problem (automatically makes the elements larger away from the boundary)
std::function<double(Geometry::PointPhysical, std::function<double(Geometry::PointPhysical) >, double) > refinement = { [] (Geometry::PointPhysical point, std::function<double(Geometry::PointPhysical) > distance, double t)->double
    {
        //solid wall boundary
        //return (distance(point)<-0.01) ? std::numeric_limits<double>::quiet_NaN() : 1;
        
        //open boundary
        return (distance(point)<-0.01 || std::abs(point[1]) < 0.01 || std::abs(point[0]) < 0.01 || std::abs(point[0] - 1) < 0.01 || std::abs(point[1] - 1) < 0.01) ? std::numeric_limits<double>::quiet_NaN() : 1;
        
        //no refinement
        //return 1;
    } };

//Note: the intended use of the prototype classes is to merge Dummy with SimpleDemoProblem
class Dummy: public Output::TecplotSingleElementWriter
{
public:
    Dummy(){}
    void writeToTecplotFile(const Base::Element* el, const Geometry::PointReference& p, std::ostream& os)
    {
        //write something so tecplot doesnt get confused
        os << 1;
    }
};

class SimpleDemoProblem : public HpgemUISimplified
{

public:

    SimpleDemoProblem() : HpgemUISimplified(2, 1, 1), t(0) { }
    
    bool initialise()
    {
        RectangularMeshDescriptor rectangularMesh(2);
        
        //bottomLeft_ and topRight_ are used here to construct an estimated bounding box of the domain
        rectangularMesh.bottomLeft_[0] = 0;
        rectangularMesh.bottomLeft_[1] = 0;
        
        //I also provide the locations of all corners of the domain as fixed nodes, this means I get away with quite a bad estimate
        //To emphasise how bad the estimate can be, use the empty bounding box
        rectangularMesh.topRight_[0] = 0;
        rectangularMesh.topRight_[1] = 0;
        
        //numElementsInDIM is a lie, the total number of nodes will end up being the product of these contributions
        //it is still a useful way to make estimates for element sizes in terms of quantities used elsewhere in hpGEM
        rectangularMesh.numElementsInDIM_[0] = 32;
        rectangularMesh.numElementsInDIM_[1] = 16;
        rectangularMesh.boundaryConditions_[0] = RectangularMeshDescriptor::SOLID_WALL;
        rectangularMesh.boundaryConditions_[1] = RectangularMeshDescriptor::SOLID_WALL;

        std::vector<PointPhysicalT> corners;

        //use std::bind to pass additional arguments to your function before passing it to the algorithm
        std::function<double(Geometry::PointPhysical) > domain = std::bind(meshDescription, std::placeholders::_1, 0);   
        
        PointPhysicalT newPoint{DIM};
        std::array<double,4> cornerLocations = {0., 0.45, 0.55, 1.};
        for(double first : cornerLocations)
        {
            newPoint[0] = first;
            for(double second : cornerLocations)
            {
                newPoint[1] = second;
                //despite the nomenclature here, it is possible to fix the location of any or all nodes in the domain
                //if you fix all nodes, the algorithm will construct a mesh based on the Delaunay triangulation of these nodes
                //it is recommended to fix at least all concave corners. The algorithm also cannot deal robustly with very sharp corners
                //or very small features, you will be warned if it is recommended to fix nodes in the vicinity of these
                //alternatively small features might be smoothed away by the algorithm
                corners.push_back(newPoint);
            }
        }    

        //there is no automated addMash functionality yet, manually construct a mesh
        Base::MeshManipulator *theMesh = new MeshManipulatorT(configData_);

        auto start = std::chrono::high_resolution_clock::now();
        
        //construct the grid; the scaling factor generally should be lowered to the elements are smaller
        theMesh->createUnstructuredMesh(rectangularMesh.bottomLeft_, rectangularMesh.topRight_, rectangularMesh.numElementsInDIM_[0] * rectangularMesh.numElementsInDIM_[1],
                                        domain, corners, std::bind(refinement, std::placeholders::_1, domain, 0), 1.05);

        auto end = std::chrono::high_resolution_clock::now();

        std::cout << "mesh generation took " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << "s" << std::endl;

        theMesh->getElementsList();
        //and add it to the mesh manipulator
        meshes_.push_back(theMesh);
        
        return true;
    }

    void elementIntegrand(const ElementT* element, const PointReferenceT& p, LinearAlgebra::Matrix& ret)
    {
        //not implemented
    }

    void faceIntegrand(const FaceT* face, const LinearAlgebra::NumericalVector& normal, const PointReferenceT& p, LinearAlgebra::Matrix& ret) 
    {
        //not implemented
    }


    void faceIntegrand(const FaceT* face, const LinearAlgebra::NumericalVector& normal,
                       const PointReferenceT& p,  LinearAlgebra::NumericalVector& ret)
    {
        //not implemented
    }

    void elementIntegrand(const ElementT* element, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret) 
    {
        //not implemented
    }

    double initialConditions(const PointPhysicalT& p)
    {
        //not implemented
        return 0.;
    }
    
    //create an output file for a single time step
    void output()
    {
        std::ofstream file2D;
        file2D.open ("out.dat");
        int dimensionsToWrite[2] = {0,1};
        Output::TecplotDiscontinuousSolutionWriter out(file2D, "RectangularMesh", "01", "one");
        Dummy d;
        out.write(meshes_[0], "holi", false, &d, t);
    }

    //create a sequence of meshes for a rotating rotor and output
    void run()
    {
        std::ofstream file2D;
        file2D.open("out.dat");
        int dimensionsToWrite[2] = { 0, 1 };
        Output::TecplotDiscontinuousSolutionWriter out(file2D, "RectangularMesh", "01", "one");
        Dummy d;

        auto start = std::chrono::high_resolution_clock::now();
        auto end = std::chrono::high_resolution_clock::now();
        double dt = M_PI / 2000.;
        double PlotT = 0.;
        double dtPlot = M_PI / 2000.;
        auto cumulativeDifference = end - start;

        //track the indices of the corners that need to be fixed
        std::vector<std::size_t> cornerIndexes;

        for (; t < 2 * M_PI; t += dt)
        {
            std::cout << t << std::endl;
            std::function<double(PointPhysicalT) > domain = std::bind(meshDescription, std::placeholders::_1, -t);
            
            cornerIndexes.clear();
            cornerIndexes.reserve(16);
            for(std::size_t i=0;i<16;++i)
            {
                cornerIndexes.push_back(i);
            }
        
            PointPhysicalT newPoint{DIM};
            std::array<double,4> cornerLocations = {0., 0.45, 0.55, 1.};
            //rotate the corners of the rotor
            std::size_t nodeIndex=0;
            std::size_t extraIndex=16;
            for(double first : cornerLocations)
            {
                for(double second : cornerLocations)
                {
                    newPoint[0] = first;
                    newPoint[1] = second;
                    if(!((std::abs(newPoint[0]) < 1e-13 || std::abs(newPoint[0] - 1.) < 1e-13) && (std::abs(newPoint[1]) < 1e-13 || std::abs(newPoint[1] - 1.) < 1e-13)))
                    {
                        double temp = std::cos(t) * (newPoint[0] - 0.5) - std::sin(t) * (newPoint[1] - 0.5) + 0.5;
                        newPoint[1] = std::sin(t) * (newPoint[0] - 0.5) + std::cos(t) * (newPoint[1] - 0.5) + 0.5;
                        newPoint[0] = temp;
                        if(domain(newPoint) > 1e-10)
                        {
                            //some basic trigonometry to keep a point of the intersection of the rotor and the boundary instead of outside the domain
                            double offset = (1. - std::cos(t)) / 2. / std::sin(t);
                            if(newPoint[0]<0. || newPoint[0]>1.)
                            {
                                newPoint[0] = first;
                                newPoint[1] = 0.5 - offset * (std::signbit(second-0.5) ? 1 : -1);
                            }
                            else
                            {
                                newPoint[1] = second;
                                newPoint[0] = 0.5 - offset * (std::signbit(first-0.5) ? 1 : -1);
                            }
                            temp = std::cos(t) * (newPoint[0] - 0.5) - std::sin(t) * (newPoint[1] - 0.5) + 0.5;
                            newPoint[1] = std::sin(t) * (newPoint[0] - 0.5) + std::cos(t) * (newPoint[1] - 0.5) + 0.5;
                            newPoint[0] = temp;
                            //a corner of the rotor that is outside of the domain will make two corners instead of one
                            //claim some extra nodes to fix them both
                            meshes_[0]->getNodes()[extraIndex] = newPoint;
                            cornerIndexes.push_back(extraIndex);
                            extraIndex++;
                            newPoint[0] = first;
                            newPoint[1] = second;
                            //some more basic trigonometry for the other point
                            temp = std::cos(t) * (newPoint[0] - 0.5) - std::sin(t) * (newPoint[1] - 0.5) + 0.5;
                            newPoint[1] = std::sin(t) * (newPoint[0] - 0.5) + std::cos(t) * (newPoint[1] - 0.5) + 0.5;
                            newPoint[0] = temp;
                            offset = (1. - std::sin(t) / 10.) / 2. / std::cos(t) + 0.5;
                            if(newPoint[0]<0. || newPoint[0]>1.)
                            {
                                newPoint[0] = first + offset * (std::signbit(first-0.5) ? 1 : -1);
                                newPoint[1] = 1 - second;
                            }
                            else
                            {
                                newPoint[1] = second + offset * (std::signbit(second-0.5) ? 1 : -1);
                                newPoint[0] = 1 - first;
                            }
                            temp = std::cos(t) * (newPoint[0] - 0.5) - std::sin(t) * (newPoint[1] - 0.5) + 0.5;
                            newPoint[1] = std::sin(t) * (newPoint[0] - 0.5) + std::cos(t) * (newPoint[1] - 0.5) + 0.5;
                            newPoint[0] = temp;
                        }
                        meshes_[0]->getNodes()[nodeIndex] = newPoint;
                    }
                    ++nodeIndex;
                }
            }
            
            start = std::chrono::high_resolution_clock::now();

            //update the mesh according to the rotated corners
            meshes_[0]->updateMesh(domain, cornerIndexes, std::bind(refinement, std::placeholders::_1, domain, -t), 1.05);

            end = std::chrono::high_resolution_clock::now();

            cumulativeDifference += end - start;

            if(t>=PlotT)
            {
                PlotT += dtPlot;
                out.write(meshes_[0], "holi", false, &d, t);
            }

        }
        std::cout << "mesh movement took " << std::chrono::duration_cast<std::chrono::seconds>(cumulativeDifference).count() << "s" << std::endl;
    }

    double t;
};

int main(int argc, char **argv)
{
    Base::parse_options(argc, argv);
    SimpleDemoProblem problem;
//     
//     problem.solve();
    problem.initialise();
    //problem.output();
    problem.run();
}


