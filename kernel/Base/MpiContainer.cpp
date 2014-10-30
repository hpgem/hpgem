/* 
 * File:   MpiContainer.cpp
 * Author: brinkf
 * 
 * Created on October 30, 2014, 4:31 PM
 */

#include "MpiContainer.hpp"
#ifdef HPGEM_USE_MPI
#include <mpi.h>
#endif

namespace Base {

MPIContainer::MPIContainer() {
#ifdef HPGEM_USE_MPI
    if(MPI::Is_initialized()){
        MPI::Group groupID=MPI::COMM_WORLD.Get_group();
        communicator_=MPI::COMM_WORLD.Create( groupID );
        processorID_=communicator_.Get_rank();
        numProcessors_=communicator_.Get_size();
    }else{
        throw "This happened way too early, please initialise MPI first";
    }
#else
    numProcessors_=1;
    processorID_=0;
#endif
}

MPIContainer::~MPIContainer() {
}

int MPIContainer::getNumProcessors(){
        return numProcessors_;
}
    
int MPIContainer::getProcessorID(){
        return processorID_;
}
    
}
