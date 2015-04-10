/* 
 * File:   MpiContainer.h
 * Author: brinkf
 *
 * Created on October 30, 2014, 4:31 PM
 */

#ifndef MPICONTAINER_HPP
#define	MPICONTAINER_HPP

#ifdef HPGEM_USE_MPI
#include <mpi.h>
#endif

#include <type_traits>
#include <vector>

namespace Base
{
    
    namespace Detail
    {
#ifdef HPGEM_USE_MPI
//template<class T>
//typename std::enable_if<std::is_class<T>::value, MPI::Datatype>::type
// toMPIType(T& t)
//{
    //fancy code that ends with calls to MPI::Datatype::Create_struct 
    //and MPI::Datatype::commit here
    //static_assert(false, "Undefined Datatype");
//}
    
    template<class T>
    typename std::enable_if<std::is_integral<T>::value, MPI::Datatype>::type
    toMPIType(T t)
    {   
        return MPI::Datatype::Match_size(MPI_TYPECLASS_INTEGER,sizeof(T));
    }

    template<class T>
    typename std::enable_if<std::is_floating_point<T>::value, MPI::Datatype>::type
    toMPIType(T t)
    {   
        return MPI::Datatype::Match_size(MPI_TYPECLASS_REAL,sizeof(T));
    }
    
    
#endif // HPGEM_USE_MPI
}

class MPIContainer
{
public:
    static MPIContainer& Instance()
    {
        static MPIContainer theInstance;
        return theInstance;
    }
    
    std::size_t getProcessorID();
    std::size_t getNumProcessors();
#ifdef HPGEM_USE_MPI
    MPI::Intracomm& getComm();

    template<class T>
    void broadcast(T& t, int id)
    {   
        MPI::Datatype type;

        communicator_.Bcast(t.data(), t.size(), Detail::toMPIType(*t.data()), id);
    }

    template<class T>
    void send(T& t, int to, int tag)
    {   
        pending_.push_back(communicator_.Isend(t.data(), t.size(), Detail::toMPIType(*t.data()), to, tag ));
    }

    template<class T>
    void receive(T& t, int to, int tag)
    {   
        pending_.push_back(communicator_.Irecv((void *)t.data(), t.size(), Detail::toMPIType(*t.data()), to, tag ));
    }

    void sync()
    {   
        MPI::Request::Waitall(pending_.size(),pending_.data());
        pending_.clear();
        communicator_.Barrier();
    }
public:
#endif
    
private:
    MPIContainer();
    MPIContainer(const MPIContainer& orig) = delete;

    std::size_t processorID_;
    std::size_t numProcessors_;
#ifdef HPGEM_USE_MPI
    std::vector<MPI::Request> pending_;
    MPI::Intracomm communicator_;
#endif
    
};

}

#endif	/* MPICONTAINER_HPP */

