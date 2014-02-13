/*
 * TecplotOutputTest.cpp
 *
 *  Created on: Feb 16, 2013
 *      Author: nicorivas
 */

#include "Base/MeshManipulator.hpp"
#include "Geometry/PointPhysical.hpp"
#include "Output/TecplotDiscontinuousSolutionWriter.hpp"
#include "Output/TecplotSingleElementWriter.hpp"
#include <vector>
#include <iostream>
#include <fstream>

class Dummy:public Output::TecplotSingleElementWriter
{
public:
    Dummy(){}
    void writeToTecplotFile(const Base::Element* el, const Geometry::PointReference& p, ostream& os)
    {
    }
};

int main()
{
    Base::ConfigurationData config(1,1,1);
    
    config.numberOfUnknowns_       = 1;
    config.numberOfTimeLevels_     = 1;
    config.numberOfBasisFunctions_ = 1;

    // Create a mesh
    Geometry::PointPhysical bottomLeft(2), topLeft(2);
    std::vector<unsigned int> numElementsOneD(2);

    bottomLeft[0]=0;
    bottomLeft[1]=0;

    topLeft[0]=8;
    topLeft[1]=8;

    numElementsOneD[0]=8;
    numElementsOneD[1]=8;

    Base::MeshManipulator myTwoDDemoMesh(&config);

    myTwoDDemoMesh.createRectangularMesh(bottomLeft,topLeft,numElementsOneD);
    myTwoDDemoMesh.outputMesh(std::cout);

    // TecPlot Output
    std::ofstream file2D;
    file2D.open ("out.dat");

    int dimensionsToWrite[2] = {0,1};

    Output::TecplotDiscontinuousSolutionWriter out(file2D,"RectangularMesh","01","xy");
    Dummy d;
    out.write(&myTwoDDemoMesh,"holi",false, &d);

}
