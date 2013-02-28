/*
 * Problem.cpp
 *
 *  Created on: Feb 21, 2013
 *      Author: nicorivas
 */

#include "Base/Base.hpp"
#include "MeshMover.hpp"
#include "Output/TecplotDiscontinuousSolutionWriter.hpp"
#include "Output/TecplotPhysicalGeometryIterator.hpp"

int main(int argc, char **argv)
{
    Base::Base<2> problem;

    problem.initializeMesh();

    Base::MeshMover<2> meshMover;

    problem.initializeMeshMover(&meshMover);

    problem.output();

    problem.solve();
}
