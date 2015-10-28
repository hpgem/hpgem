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
#include "Logger.h"

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

    ///a class providing wrapper functions for most commonly used MPI routines.
    ///all functions where only one processor has to do something different
    ///will default to using processor 0 as the special one unless an argument is specified
class MPIContainer
{
public:
    ///fetch the instance to be used for communication
    static MPIContainer& Instance()
    {
        static MPIContainer theInstance;
        return theInstance;
    }
    
    ///get the unique identifier associated with this processor
    std::size_t getProcessorID();
    
    ///\deprecated Does not conform naming conventions, use getNumberOfProcessors instead
    std::size_t getNumProcessors()
    {
        return getNumberOfProcessors();
    }
    
    ///get the total number of processors participating in this simulation
    std::size_t getNumberOfProcessors();

    ///pass a function that only has to be executed on one processor (that is not processor 0)
    ///this function cannot return anything, but it is allowed to alter its arguments
    template<typename... Args>
    void onlyOnOneProcessor(int id, std::function<void(Args...)> f, Args... args)
    {
        //remember to add accolades when there is more than one statement in the if clause
        //(omitted to reduce the number of #ifdefs)
#ifdef HPGEM_USE_MPI
        if(getProcessorID() == id)
#endif
            f(args...);
    }

    ///pass a function that only has to be executed on one processor
    ///this function cannot return anything, but it is allowed to alter its arguments
    template<typename... Args>
    void onlyOnOneProcessor(std::function<void(Args...)> f, Args... args)
    {
        onlyOnOneProcessor(0, f, args...);
    }

    ///process all pending asynchronous communication requests before continuing
    void sync()
    {
#if HPGEM_USE_MPI
        MPI::Request::Waitall(pending_.size(),pending_.data());
        pending_.clear();
        communicator_.Barrier();
#endif
        //we are automatically synced if there is no MPI
    }

    ///make sure a vector of information is available on all processors participating in this simulation
    ///all processors must call this function. The processor specified by the identifier id will do the sending
    ///the rest will do the recieving
    template<typename T>
    typename std::enable_if<!std::is_scalar<T>::value, void>::type
    broadcast(T& t, int id = 0)
    {   
#if HPGEM_USE_MPI
        if(t.size() > 0)
        {
            communicator_.Bcast(t.data(), t.size(), Detail::toMPIType(*t.data()), id);
        }
#endif
        //we are the only processor if there is no MPI
    }

    ///make sure a scalar is available on all processors participating in this simulation
    ///all processors must call this function. The processor specified by the identifier id will do the sending
    ///the rest will do the recieving
    template<typename T>
    typename std::enable_if<std::is_scalar<T>::value, void>::type
    broadcast(T& t, int id = 0)
    {
#if HPGEM_USE_MPI
        communicator_.Bcast(&t, 1, Detail::toMPIType(t), id);
#endif
        //we are the only processor if there is no MPI
    }

    ///make sure a scalar is available on all processors participating in this simulation
    ///all processors must call this function. The processor specified by the identifier id will do the sending
    ///the rest will do the recieving
    template<typename T>
    typename std::enable_if<std::is_scalar<T>::value, void>::type
    broadcast(std::complex<T>& t, int id = 0)
    {
#if HPGEM_USE_MPI
        communicator_.Bcast(&t, 1, Detail::toMPIType(t), id);
#endif
        //we are the only processor if there is no MPI
    }

    ///asynchronously send a vector to some other processor. This process must call sync before the buffer can be used again.
    ///the recieving process must call recieve with appropriate parameters and then call sync before it can use the data
    ///\param t the data
    ///\param to the processor to recieve the information
    ///\param tag an identifier for this specific send request. This must be unique among all send requests between
    ///the previous synchronisation step and the next one. Exactly one recieve request must also provide this tag and
    ///it must be done on the process specified by the 'to' parameter
    template<typename T>
    typename std::enable_if<!std::is_scalar<T>::value, void>::type
    send(T& t, int to, int tag)
    {   
        logger.assert(to < getNumberOfProcessors(), "trying to send to processor %, but there are only % processors available", to, getNumberOfProcessors());
        if(to == getProcessorID())
        {
            logger(WARN, "sending data to self!");
        }
#if HPGEM_USE_MPI
        if(t.size() > 0)
        {
            pending_.push_back(communicator_.Isend(t.data(), t.size(), Detail::toMPIType(*t.data()), to, tag ));
        }
#endif
    }

    ///asynchronously send a scalar to some other processor. This process must call sync before the buffer can be used again.
    ///the recieving process must call recieve with appropriate parameters and then call sync before it can use the data
    ///\param t the data
    ///\param to the processor to recieve the information
    ///\param tag an identifier for this specific send request. This must be unique among all send requests between
    ///the previous synchronisation step and the next one. Exactly one recieve request must also provide this tag and
    ///it must be done on the process specified by the 'to' parameter
    template<typename T>
    typename std::enable_if<std::is_scalar<T>::value, void>::type
    send(T& t, int to, int tag)
    {
        logger.assert(to < getNumberOfProcessors(), "trying to send to processor %, but there are only % processors available", to, getNumberOfProcessors());
        if(to == getProcessorID())
        {
            logger(WARN, "sending data to self!");
        }
#if HPGEM_USE_MPI
        pending_.push_back(communicator_.Isend(&t, 1, Detail::toMPIType(t), to, tag ));
#endif
    }

    ///asynchronously send a scalar to some other processor. This process must call sync before the buffer can be used again.
    ///the recieving process must call recieve with appropriate parameters and then call sync before it can use the data
    ///\param t the data
    ///\param to the processor to recieve the information
    ///\param tag an identifier for this specific send request. This must be unique among all send requests between
    ///the previous synchronisation step and the next one. Exactly one recieve request must also provide this tag and
    ///it must be done on the process specified by the 'to' parameter
    template<typename T>
    typename std::enable_if<std::is_scalar<T>::value, void>::type
    send(std::complex<T>& t, int to, int tag)
    {
        logger.assert(to < getNumberOfProcessors(), "trying to send to processor %, but there are only % processors available", to, getNumberOfProcessors());
        if(to == getProcessorID())
        {
            logger(WARN, "sending data to self!");
        }
#if HPGEM_USE_MPI
        pending_.push_back(communicator_.Isend(&t, 1, Detail::toMPIType(t), to, tag ));
#endif
    }

    ///asynchronously receive a vector from some other processor. This process must call sync before it can use the data.
    ///the sending process must call send with appropriate parameters and then call sync to allow this process to finish the sync
    ///\param t the data
    ///\param from the processor that sends the information
    ///\param tag an identifier for this specific receive request. This must be unique among all receive requests between
    ///the previous synchronisation step and the next one. Exactly one send request must also provide this tag and
    ///it must be done on the process specified by the 'from' parameter
    template<typename T>
    typename std::enable_if<!std::is_scalar<T>::value, void>::type
    receive(T& t, int from, int tag)
    {   
        logger.assert(from < getNumberOfProcessors(), "trying to receive from processor %, but there are only % processors available", from, getNumberOfProcessors());
        if(from == getProcessorID())
        {
            logger(WARN, "getting data from self!");
        }
#if HPGEM_USE_MPI
        if(t.size() > 0)
        {
            pending_.push_back(communicator_.Irecv((void *)t.data(), t.size(), Detail::toMPIType(*t.data()), from, tag ));
        }
#endif
    }

    ///asynchronously receive a scalar from some other processor. This process must call sync before it can use the data.
    ///the sending process must call send with appropriate parameters and then call sync to allow this process to finish the sync
    ///\param t the data
    ///\param from the processor that sends the information
    ///\param tag an identifier for this specific receive request. This must be unique among all receive requests between
    ///the previous synchronisation step and the next one. Exactly one send request must also provide this tag and
    ///it must be done on the process specified by the 'from' parameter
    template<typename T>
    typename std::enable_if<std::is_scalar<T>::value, void>::type
    receive(T& t, int from, int tag)
    {
        logger.assert(from < getNumberOfProcessors(), "trying from receive from processor %, but there are only % processors available", from, getNumberOfProcessors());
        if(from == getProcessorID())
        {
            logger(WARN, "getting data from self!");
        }
#if HPGEM_USE_MPI
        pending_.push_back(communicator_.Irecv((void *)&t, 1, Detail::toMPIType(t), from, tag ));
#endif
    }

    ///asynchronously receive a scalar from some other processor. This process must call sync before it can use the data.
    ///the sending process must call send with appropriate parameters and then call sync to allow this process to finish the sync
    ///\param t the data
    ///\param from the processor that sends the information
    ///\param tag an identifier for this specific receive request. This must be unique among all receive requests between
    ///the previous synchronisation step and the next one. Exactly one send request must also provide this tag and
    ///it must be done on the process specified by the 'from' parameter
    template<typename T>
    typename std::enable_if<std::is_scalar<T>::value, void>::type
    receive(std::complex<T>& t, int from, int tag)
    {
        logger.assert(from < getNumberOfProcessors(), "trying to receive from processor %, but there are only % processors available", from, getNumberOfProcessors());
        if(from == getProcessorID())
        {
            logger(WARN, "getting data from self!");
        }
#if HPGEM_USE_MPI
        pending_.push_back(communicator_.Irecv((void *)&t, 1, Detail::toMPIType(t), from, tag ));
#endif
    }

#ifdef HPGEM_USE_MPI

    ///make sure a vector of data gets collected on one processor. This routine should be called by all processes.
    ///The process that does the recieving is specified by id, the others will do the sending
    ///the MPI::OP parameter specifies how the data is to be treated to make it fit in one slot
    ///popular options include MPI::MAX, MPI::MIN, MPI::SUM and MPI::PROD to keep the
    ///maximum, the minimum, take the sum or the product respectively. Note that operations are only defined where they
    ///make sense for the used data type
    template<typename T>
    typename std::enable_if<!std::is_scalar<T>::value, void>::type
    reduce(T& t, MPI::Op operation, int id = 0)
    {
        if(t.size() > 0)
        {
            if(id == getProcessorID())
            {
                communicator_.Reduce(MPI::IN_PLACE, t.data(), t.size(), Detail::toMPIType(*t.data()), operation, id);
            }
            else
            {
                communicator_.Reduce(t.data(), nullptr, t.size(), Detail::toMPIType(*t.data()), operation, id);
            }
        }
    }

    ///make sure a scalar gets collected on one processor. This routine should be called by all processes.
    ///The process that does the recieving is specified by id, the others will do the sending
    ///the MPI::OP parameter specifies how the data is to be treated to make it fit in one slot
    ///popular options include MPI::MAX, MPI::MIN, MPI::SUM and MPI::PROD to keep the
    ///maximum, the minimum, take the sum or the product respectively. Note that operations are only defined where they
    ///make sense for the used data type
    template<typename T>
    typename std::enable_if<std::is_scalar<T>::value, void>::type
    reduce(T& t, MPI::Op operation, int id = 0)
    {
        if(id == getProcessorID())
        {
            communicator_.Reduce(MPI::IN_PLACE, &t, 1, Detail::toMPIType(t), operation, id);
        }
        else
        {
            communicator_.Reduce(&t, nullptr, 1, Detail::toMPIType(t), operation, id);
        }
    }

    ///make sure a vector of data gets collected on one processor. This routine should be called by all processes.
    ///The process that does the recieving is specified by id, the others will do the sending
    ///the MPI::Op parameter specifies how the data is to be treated to make it fit in one slot
    ///popular options include MPI::SUM and MPI::PROD to take the
    ///sum or the product respectively. Note that operations are only defined where they
    ///make sense for the used data type
    template<typename T>
    typename std::enable_if<std::is_scalar<T>::value, void>::type
    reduce(std::complex<T>& t, MPI::Op operation, int id = 0)
    {
        if(id == getProcessorID())
        {
            communicator_.Reduce(MPI::IN_PLACE, &t, 1, Detail::toMPIType(t), operation, id);
        }
        else
        {
            communicator_.Reduce(&t, nullptr, 1, Detail::toMPIType(t), operation, id);
        }
    }

    ///retrieve a communicator for use in direct MPI calls. Please be aware that an unspecified amount of resources,
    ///such as asynchronous communication tags, is already in use by hpGEM
    MPI::Intracomm& getComm();
public:
#endif
    
private:
    MPIContainer();
    MPIContainer(const MPIContainer& orig) = delete;

    std::size_t processorID_;
    std::size_t numberOfProcessors_;
#ifdef HPGEM_USE_MPI
    std::vector<MPI::Request> pending_;
    MPI::Intracomm communicator_;
#endif
    
};

}

#endif	/* MPICONTAINER_HPP */

