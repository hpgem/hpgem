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
    typename std::enable_if<std::is_scalar<T>::value, void>::type
    broadcast(std::complex<T>& t, int id)
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
    typename std::enable_if<std::is_scalar<T>::value, void>::type
    send(std::complex<T>& t, int to, int tag)
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

    template<typename T>
    typename std::enable_if<std::is_scalar<T>::value, void>::type
    receive(std::complex<T>& t, int to, int tag)
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

