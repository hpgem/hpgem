/* 
 * File:   MpiContainer.hpp
 * Author: brinkf
 *
 * Created on October 30, 2014, 4:31 PM
 */

#ifdef HPGEM_USE_MPI
#include <mpi.h>
#endif

#ifndef MPICONTAINER_HPP
#define	MPICONTAINER_HPP

namespace Base{

class MPIContainer {
public:
    static MPIContainer& Instance(){
        static MPIContainer theInstance;
        return theInstance;
    }
    
    int getProcessorID();
    int getNumProcessors();
#ifdef HPGEM_USE_MPI
    MPI::Intracomm& getComm();
    
    template<class T>
    void broadcast(T&);
    
    template<class T>
    void send(T&);
#endif
    
private:
    MPIContainer();
    MPIContainer(const MPIContainer& orig)=delete;
    virtual ~MPIContainer();
    
    int processorID_;
    int numProcessors_;
#ifdef HPGEM_USE_MPI
    MPI::Intracomm communicator_;
#endif

};
}

#endif	/* MPICONTAINER_HPP */

