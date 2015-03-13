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

#include "Logger.h"

#ifdef HPGEM_STACKTRACE_ENABLE
//To create stacktraces, we need the
//  backtrace(3)
//function calls.. (really, we don't want to do this by hand)
#include <execinfo.h>

#ifdef HPGEM_STACKTRACE_DEMANGLE
//However, we would end up with mangled function names...
//So, instead we'll use the abi::__cxa_demangle function call..
#include <cxxabi.h>
//Oh, and we really want to get the function names as well..
#include <dlfcn.h>
#endif

#endif
#include <cstdlib>
#include <iostream>
#include <csignal>
/*
 *  We need these to actually exists. These are used as tags in the template metaprogramming for
 *  the Logger class.
 */
LL<Log::FATAL> FATAL;
LL<Log::ERROR> ERROR;
LL<Log::WARN> WARN;
LL<Log::INFO> INFO;
LL<Log::DEFAULT> DEFAULT;
LL<Log::VERBOSE> VERBOSE;
LL<Log::DEBUG> DEBUG;

/* Actual definition of the default logger. */
Logger<HPGEM_LOGLEVEL> logger("hpGEM Kernel");

/* Default implementation for logging warnings / messages */
static void printMessage(std::string module, std::string msg)
{
    std::cout << "Module " << module << ":\n" << msg << std::endl;
}

/* Default implementation for logging errors / fatals */
static void printError(std::string module, std::string msg)
{
    std::cerr << "Module " << module << ":\n" << msg << std::endl;
#ifdef HPGEM_STACKTRACE_ENABLE
    std::cerr << "\n-----------------[Stack Trace]-----------------\n";
    
    void* stackBuffer[64]; //This should be enough for all purposes..
    //First, we retrieve the addresses of the entire stack...
    int nStackFrames = backtrace(stackBuffer, 64);
#ifndef HPGEM_STACKTRACE_DEMANGLE
    //We don't have the demangling infra, so just use backtrace_symbols.
    char** functionNames = backtrace_symbols(stackBuffer, nStackFrames);
    for( int i = 0; i < nStackFrames; i++ )
    {   
        std::cerr << '\t' << functionNames[i] << '\n';
    }
    std::cerr << "Exiting.\n" << std::endl;

    //DO NOT USE DELETE HERE. THIS SHOULD BE free()'d!
    // -- dducks
    free(functionNames);
#else
    //We request the symbol information ourselves, in order to be able to demangle it.
    //And request the function names using dladdr.
    Dl_info infoStruct;
    for (int i = 4; i < nStackFrames; i++)
    {
        if (dladdr(stackBuffer[i], &infoStruct))
        { // We succesfully loaded the address...
            int demangleStatus;
            char* fnDemangled = abi::__cxa_demangle(infoStruct.dli_sname, NULL, NULL, &demangleStatus);
            if (infoStruct.dli_sname == nullptr)
                continue;
            
            //We even succesfully demangled the symbol...
            if (demangleStatus == 0)
            {
                std::cerr << fnDemangled << " +" << (void*) ((char*) stackBuffer[i] - (char*) infoStruct.dli_saddr) << "\t(" << infoStruct.dli_fname << ")\n";
                free(fnDemangled);
            }
            else
            { //Well, we tried. Lets output at least our raw symbol name.
                std::cerr << infoStruct.dli_sname << " +" << (void*) ((char*) stackBuffer[i] - (char*) infoStruct.dli_saddr) << "\t(" << infoStruct.dli_fname << ")\n";
            }
        }
        else
        { //Name lookup failed.
            std::cerr << stackBuffer[i] << ": ?????" << std::endl;
        }
    }
#endif
#endif
    //send a signal first, in case a debugger can catch it
    std::raise(SIGTERM);
    std::exit(2);
}

// Default output methods.
LoggerOutput loggerOutputDefaultImpl = {printError, //onFatal
        printError, //onError
        printMessage, //onWarn
        printMessage, //onInfo
        printMessage, //onVerbose
        printMessage //onDebug
        };

//And we assign them.
LoggerOutput* loggerOutput = &loggerOutputDefaultImpl;
