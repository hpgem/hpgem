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
#include <complex>

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
        //convert integral data to the corresponding MPI type
        template<typename T>
        typename std::enable_if<std::is_integral<T>::value, MPI::Datatype>::type
        toMPIType(T t)
        {
            return MPI::Datatype::Match_size(MPI_TYPECLASS_INTEGER,sizeof(T));
        }

        //convert floating point data to the corresponding MPI type
        template<typename T>
        typename std::enable_if<std::is_floating_point<T>::value, MPI::Datatype>::type
        toMPIType(T t)
        {
            return MPI::Datatype::Match_size(MPI_TYPECLASS_REAL,sizeof(T));
        }

        //convert complex data to the corresponding MPI type
        template<typename T>
        typename std::enable_if<std::is_floating_point<T>::value, MPI::Datatype>::type
        toMPIType(std::complex<T> t)
        {
            return MPI::Datatype::Match_size(MPI_TYPECLASS_COMPLEX,sizeof(std::complex<T>));
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

    template<typename T>
    typename std::enable_if<!std::is_scalar<T>::value, void>::type
    broadcast(T& t, int id)
    {   
        if(t.size() > 0)
        {
            communicator_.Bcast(t.data(), t.size(), Detail::toMPIType(*t.data()), id);
        }
    }

    template<typename T>
    typename std::enable_if<std::is_scalar<T>::value, void>::type
    broadcast(T& t, int id)
    {
        communicator_.Bcast(&t, 1, Detail::toMPIType(t), id);
    }

    template<typename T>
    typename std::enable_if<!std::is_scalar<T>::value, void>::type
    send(T& t, int to, int tag)
    {   
        if(t.size() > 0)
        {
            pending_.push_back(communicator_.Isend(t.data(), t.size(), Detail::toMPIType(*t.data()), to, tag ));
        }
    }

    template<typename T>
    typename std::enable_if<std::is_scalar<T>::value, void>::type
    send(T& t, int to, int tag)
    {
        pending_.push_back(communicator_.Isend(&t, 1, Detail::toMPIType(t), to, tag ));
    }

    template<typename T>
    typename std::enable_if<!std::is_scalar<T>::value, void>::type
    receive(T& t, int to, int tag)
    {   
        if(t.size() > 0)
        {
            pending_.push_back(communicator_.Irecv((void *)t.data(), t.size(), Detail::toMPIType(*t.data()), to, tag ));
        }
    }

    template<typename T>
    typename std::enable_if<std::is_scalar<T>::value, void>::type
    receive(T& t, int to, int tag)
    {
        pending_.push_back(communicator_.Irecv((void *)&t, 1, Detail::toMPIType(t), to, tag ));
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

