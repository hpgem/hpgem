#==========[ assert_debug Config ]============
set(hpGEM_DISABLE_ASSERTS OFF CACHE BOOL "Disables asserts - use at your own risk!")
mark_as_advanced( FORCE hpGEM_DISABLE_ASSERTS )
if (hpGEM_DISABLE_ASSERTS)
    add_definitions( -DHPGEM_NO_ASSERTS )
endif()

set(hpGEM_ENFORCE_ASSERTS OFF CACHE BOOL "Enables asserts even in release mode - may be slow!")
mark_as_advanced( FORCE hpGEM_ENFORCE_ASSERTS )
if (hpGEM_ENFORCE_ASSERTS)
    if(hpGEM_DISABLE_ASSERTS)
        message(FATAL_ERROR, "Cannot disable the assertions when they are enforced!!!")
    endif()
    add_definitions( -DHPGEM_FORCE_ASSERTS )
endif()

#=======[ The default loglevel ]========
set(hpGEM_LOGLEVEL "DEFAULT" CACHE STRING "Verbosity of hpGEM. Default is recommended.")
set_property( CACHE hpGEM_LOGLEVEL PROPERTY STRINGS FATAL ERROR WARN INFO DEFAULT VERBOSE DEBUG )
mark_as_advanced( FORCE hpGEM_LOGLEVEL )
add_definitions( -DHPGEM_LOGLEVEL=Log::${hpGEM_LOGLEVEL} )

#===========[ Stacktraces ]=============
#we need this to check if we can
#even make stacktraces
include(CheckIncludeFiles)
include(CheckIncludeFileCXX)
include(CheckLibraryExists)

#Now, we need execinfo() to do it a bit
#in a sane way (nobody likes assembly)
CHECK_INCLUDE_FILES(execinfo.h  HAS_EXECINFO_H)

#Now, we would also like to have some
#name demangling going on...
#dynlinker
CHECK_INCLUDE_FILES(dlfcn.h     HAS_DLFCN_H)
#CXX demangler api
CHECK_INCLUDE_FILE_CXX(cxxabi.h HAS_CXXABI_H)
CHECK_LIBRARY_EXISTS(dl dladdr "" HAS_DLADDR_LIB)

#So, now see what we've end up with.
if (${HAS_EXECINFO_H})
    #the bare minimum is there...
    if (CMAKE_BUILD_TYPE)
        if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
            set( hpGEM_STACKTRACE_ENABLE ON CACHE BOOL "Enable stacktraces for log messages")
        else()
            set( hpGEM_STACKTRACE_ENABLE OFF CHACHE BOOL "Enable stacktraces for log messages")
        endif()
    endif()
    #now, check whether or not we can do demangling
    set( USE_DEMANGLE ON )
    if (NOT HAS_DLFCN_H)
        set( USE_DEMANGLE OFF )
    endif()
    if (NOT HAS_DLADDR_LIB)
        set( USE_DEMANGLE OFF )
    endif()
    if (NOT HAS_CXXABI_H)
        set( USE_DEMANGLE OFF )
    endif()
    set( hpGEM_STACKTRACE_DEMANGLE ${USE_DEMANGLE} CACHE BOOL "Enable name demangling for stacktraces")
else()
    set( hpGEM_STACKTRACE_ENABLE   OFF CACHE INTERNAL "execinfo.h not found - no stacktraces possible")
    set( hpGEM_STACKTRACE_DEMANGLE OFF CACHE INTERNAL "demangling is not going to happen, no bt!")
endif()
mark_as_advanced( FORCE hpGEM_STACKTRACE_ENABLE hpGEM_STACKTRACE_DEMANGLE )
#Well, that should work now...

if (hpGEM_STACKTRACE_ENABLE)
    add_definitions(-DHPGEM_STACKTRACE_ENABLE)
endif()
if (hpGEM_STACKTRACE_DEMANGLE)
    add_definitions(-DHPGEM_STACKTRACE_DEMANGLE)
    link_libraries(dl)
endif()