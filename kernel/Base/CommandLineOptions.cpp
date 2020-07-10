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

#include "CommandLineOptions.h"
//#include "MpiContainer.h"
#include "Logger.h"
#include <cstring>

#ifdef HPGEM_USE_MPI
#include <mpi.h>
#endif

#ifdef HPGEM_USE_ANY_PETSC
#include <petscsys.h>
#endif

#ifdef HPGEM_USE_SLEPC
#include <slepcsys.h>
#endif
namespace hpgem {
// special flags
auto& isDone = Base::register_argument<bool>(
    '\0', "",
    "Signals the end of the arguments passed to hpGEM, the rest will be passed "
    "to linked libraries",
    false, false);
auto& warnCrash = Base::register_argument<bool>(
    '\0', "crashOnWarn",
    "Tell hpGEM to crash when the logger issues warnings \n(useful in "
    "combination with a debugger to quickly find the source of warnings)",
    false, false);

std::map<std::string, Base::Detail::CommandLineOptionBase*>&
    Base::Detail::getCLOMapping_long() {
    static std::map<std::string, CommandLineOptionBase*> mapping;
    return mapping;
}
std::map<char, Base::Detail::CommandLineOptionBase*>&
    Base::Detail::getCLOMapping_short() {
    static std::map<char, CommandLineOptionBase*> mapping;
    return mapping;
}
std::vector<Base::Detail::CommandLineOptionBase*>& Base::Detail::getCLOList() {
    static std::vector<Base::Detail::CommandLineOptionBase*> list;
    return list;
}

static auto& printHelp = Base::register_argument<bool>(
    '?', "help", "Prints this help message", false, false);

static bool hasParsed = false;
bool Base::parse_isDone() { return hasParsed; }
void Base::parse_options(int argc, char** argv) {
    logger.assert_debug(!hasParsed, "Arguments have already been parsed");
#ifdef HPGEM_USE_MPI

    int init;
    MPI_Initialized(&init);

    if (init == 0) {
        // somebody might have been kind enough to actually
        // initialize MPI for us... so yeah. Let's prevent
        // initializing it twice, as this will end in chaos and mayhem
        MPI_Init(&argc, &argv);
        ;

        // We should call MPI::Finalize() when we quit. As we have no
        // clue when people are actually going to call this, we just
        // register an on-exit handler.
        std::atexit([]() { MPI_Finalize(); });
    }

#endif

    Base::Detail::CLOParser parser(argc, argv);
    hasParsed = true;
    int count = parser.go();

    if (warnCrash.getValue()) {
        loggerOutput->onWarn = loggerOutput->onFatal;
    }

    // move the name of the program to the new 'beginning' of the arguments
    argv[count] = argv[0];

    argc -= count;
    argv += count;

#if defined(HPGEM_USE_ANY_PETSC)
    // block so variables can be created without affecting the rest of the
    // function
    {
        PetscBool initialised;
        PetscInitialized(&initialised);
        if (initialised == PETSC_FALSE) {
#ifdef HPGEM_USE_MPI
            // if we know mpi exists make sure PETSc based communication does
            // not happen on COMM_WORLD communicating on COMM_WORLD is a bad
            // idea if you are a library and are not sure who else might use MPI
            //(PETSc CLAIMS this is not needed, but also does not provide this
            // safeguard itself)
            MPI_Comm_dup(MPI_COMM_WORLD, &PETSC_COMM_WORLD);
#endif
#ifdef HPGEM_USE_SLEPC
            SlepcInitialize(&argc, &argv, PETSC_NULL, "PETSc help\n");
#else
            PetscInitialize(&argc, &argv, PETSC_NULL, "PETSc help\n");
#endif
            // please do not catch signals PETSc, you are confusing users
            PetscPopSignalHandler();
            std::atexit([]() {
#ifdef HPGEM_USE_SLEPC
                SlepcFinalize();
#else
                PetscFinalize();
#endif
            });
        }
    }
#endif
}
int Base::Detail::CLOParser::go() {
    auto& lmapping = getCLOMapping_long();
    auto& smapping = getCLOMapping_short();

    char** firstArg = this->pData;
    bool done = false;

    ++(*this);
    try {
        while (remaining() && !done) {
            std::string tag = *(*this);

            Base::Detail::CommandLineOptionBase* pBase;

            logger.assert_always(tag.size() >= 2 && tag[0] == '-',
                                 "Unknown argument % found", tag);
            if (tag[1] == '-') {
                std::string realtag = tag.substr(2);
                auto it = lmapping.find(realtag);
                logger.assert_always(it != lmapping.end(),
                                     "Argument % not found.", tag);
                logger.assert_always(
                    !(it->second->hasArgument() && remaining() <= 1),
                    "Argument % has an argument, but was not provided one.",
                    tag);

                pBase = it->second;

                pBase->parse(*this);
            } else {
                for (std::string::size_type i = 1; i < tag.size(); i++) {
                    auto it = smapping.find(tag[i]);
                    logger.assert_always(it != smapping.end(),
                                         "Argument % not found.", tag);
                    logger.assert_always(
                        !(it->second->hasArgument() &&
                          (i < (tag.size() - 1) || remaining() <= 1)),
                        "Argument % has an argument, but was not provided one.",
                        tag);

                    pBase = it->second;

                    pBase->parse(*this);
                }
            }
            if (isDone.isUsed() && isDone.getValue()) {
                *(this->pData) = *firstArg;
                done = true;
            } else {
                ++(*this);
            }
        }
    } catch (const std::string& explain) {
        std::cerr << explain << std::endl;
        for (std::size_t i = 0; i < count; i++) {
            std::cerr << pData[i] << ' ';
        }
        std::cerr << std::endl;
        for (std::size_t i = 0; i < currCount; i++) {
            for (std::size_t j = 0; j < std::strlen(pData[i]); j++) {
                std::cerr << ' ';
            }
            std::cerr << ' ';
        }
        for (std::size_t j = 0; j < std::strlen(pData[currCount]); j++) {
            std::cerr << '^';
        }
        std::cerr << std::endl;
        std::exit(1);
    }

    bool kill = false;
    for (Base::Detail::CommandLineOptionBase* base :
         Base::Detail::getCLOList()) {
        // do not use logger, because also print help
        if (base->isRequired() && !base->isUsed()) {
            std::cerr << "Argument \'" << base->getLongTag()
                      << "\' required.\n";
            kill = true;
        }
    }
    if (printHelp.getValue() || kill) {
        std::cout << "Usage: " << pData[0] << " [args...] with args: \n\n";
        for (Base::Detail::CommandLineOptionBase* base :
             Base::Detail::getCLOList()) {
            std::cout << '\t';
            if (base->getTag() != '\0')
                std::cout << "'-" << base->getTag() << "', ";
            std::cout << "'--" << base->getLongTag() << "' ";
            if (base->hasArgument()) std::cout << "[args...]";
            if (base->isRequired()) std::cout << " - REQUIRED";
            std::cout << "\n\t  " << base->getDescription() << "\n\n";
        }
        kill = true;

        std::cout << "\n\n---------------------------\nFeatures:\n";
#ifdef HPGEM_USE_MPI
        std::cout << "\tMPI\n";
#endif
#ifdef HPGEM_USE_METIS
        std::cout << "\tMETIS\n";
#endif
#ifdef HPGEM_USE_PETSC
        std::cout << "\tPETSC\n";
#endif
#ifdef HPGEM_USE_COMPLEX_PETSC
        std::cout << "\tComplex PETSC\n";
#endif
#ifdef HPGEM_USE_SLEPC
        std::cout << "\tSLEPC\n";
#endif
#ifdef HPGEM_USE_QHULL
        std::cout << "\tQHull\n";
#endif
    }
    logger.assert_always(!kill, "Help printed, terminating program...");
    return currCount;
}

}  // namespace hpgem
