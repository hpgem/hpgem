/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2014, University of Twente
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors
 may be used to endorse or promote products derived from this software without
 specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "MpiContainer.h"
#ifdef HPGEM_USE_MPI
#include <mpi.h>
#endif

#include "Logger.h"

namespace Base {
namespace Detail {
#ifdef HPGEM_USE_MPI
MPI_Datatype matchSizeInternal(int typeclass, int size) {
    MPI_Datatype type;
    (void)MPI_Type_match_size(typeclass, size, &type);
    return type;
}
#endif
}  // namespace Detail

MPIContainer::MPIContainer() {
#ifdef HPGEM_USE_MPI

    int flag;
    MPI_Initialized(&flag);  // No error
    logger.assert_debug(flag != 0,
                        "Please initialise MPI first before"
                        " calling the constructor of MPIContainer");
    // TODO: There are plenty of places where we use the world comm
    // this should be replaced by this comm.

    // Assume that we can clone the world comm without error.
    MPI_Comm_dup(MPI_COMM_WORLD, &communicator_);
    // Make sure that errors are fatal, thus absolving us from checking for
    // the return codes.
    MPI_Comm_set_errhandler(communicator_, MPI_ERRORS_ARE_FATAL);
    // Cache some useful information
    MPI_Comm_rank(communicator_, &processorID_);
    MPI_Comm_size(communicator_, &numberOfProcessors_);
#else
    numberOfProcessors_ = 1;
    processorID_ = 0;
#endif
}

int MPIContainer::getNumberOfProcessors() { return numberOfProcessors_; }

int MPIContainer::getProcessorID() { return processorID_; }

#ifdef HPGEM_USE_MPI

void MPIContainer::sendWrapper(const void *buf, int count,
                               MPI_Datatype datatype, int dest, int tag) {
    pending_.emplace_back();
    MPI_Request *request = &*(pending_.end() - 1);
    MPI_Isend(buf, count, datatype, dest, tag, communicator_, request);
}

void MPIContainer::receiveWrapper(void *buf, int count, MPI_Datatype datatype,
                                  int source, int tag) {
    pending_.emplace_back();
    MPI_Request *request = &*(pending_.end() - 1);
    MPI_Irecv(buf, count, datatype, source, tag, communicator_, request);
}

void MPIContainer::reduceWrapper(void *buffer, int count, MPI_Datatype datatype,
                                 MPI_Op op, int root) {
    void *sendBuffer, *receiveBuffer;
    if (root == processorID_) {
        sendBuffer = MPI_IN_PLACE;
        receiveBuffer = buffer;
    } else {
        sendBuffer = buffer;
        receiveBuffer = nullptr;
    }
    MPI_Reduce(sendBuffer, receiveBuffer, count, datatype, op, root,
               communicator_);
}

MPI_Comm &MPIContainer::getComm() { return communicator_; }
#endif

}  // namespace Base
