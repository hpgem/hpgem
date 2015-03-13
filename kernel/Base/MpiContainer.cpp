/* 
 * File:   MpiContainer.cpp
 * Author: brinkf
 * 
 * Created on October 30, 2014, 4:31 PM
 */

#include "MpiContainer.h"
#ifdef HPGEM_USE_MPI
#include <mpi.h>
#endif

namespace Base
{
    
    MPIContainer::MPIContainer()
    {
#ifdef HPGEM_USE_MPI
        logger.assert(MPI::Is_initialized(), "Please initialise MPI first before"
            " calling the constructor of MPIContainer");
        MPI::Group groupID=MPI::COMM_WORLD.Get_group();
        communicator_=MPI::COMM_WORLD.Create( groupID );
        processorID_=communicator_.Get_rank();
        numProcessors_=communicator_.Get_size();
#else
        numProcessors_ = 1;
        processorID_ = 0;
#endif
    }
    
    MPIContainer::~MPIContainer()
    {
    }
    
    int MPIContainer::getNumProcessors()
    {
        return numProcessors_;
    }
    
    int MPIContainer::getProcessorID()
    {
        return processorID_;
    }

#ifdef HPGEM_USE_MPI
MPI::Intracomm& MPIContainer::getComm()
{   
    return communicator_;
}
#endif

}
