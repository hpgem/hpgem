# Online developers guide

The developer notes contain general information for developers that is useful no matter what part of the code you are editing. As such it does
not provide arguments for or against specific design choices. If you are
interested in an introduction about how to use hpGEM, you may find the readme more useful.

## Coding conventions

### Header

Please add the following text to each top each source file:

>
>/*
>This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at 
>http://hpgem.org/about-the-code/team 
>
>This code is distributed using BSD 3-Clause License. A copy of which can found below.
>
>
>Copyright (c) 2014, Univesity of Twenete
>All rights reserved.
>
>Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
>
>1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
>
>2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
>
>3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
>
>THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A >PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT >LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT >(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
>*/

### Formatting

The code is formatted using [clang-format](https://clang.llvm.org/docs/ClangFormat.html) using the Google formatting conventions. The appropriate ``.clang-format`` is in the main directory. 


### General Recommendations

- Any violation to the guide is allowed if it enhances readability.

### Naming Conventions

- File names must be of the form ``FileName``.
- Directory names must be of the form ``DirectoryName``.
- Names representing types must be in mixed case starting with upper case, i.e. ``TypeName``.
- Variable names must be in mixed case starting with lower case, i.e. ``variableName``.
- Named constants (including enumeration values) must be all uppercase using underscore to separate words.
- Names representing methods or functions must be verbs and written in mixed case starting with lower case, i.e. ``functionName``.
- Names representing Classes must be of the form ``ClassName``.
- Names representing template types should be a single uppercase letter.
- Private class variables should have ``_`` suffix.
- Never use identifiers which begin with one or two underscores (`_' or `__'). In most cases those variables are used in imported libraries.
- Plural form should be used on names representing a collection of objects.
- A variable with a large scope should have a long name.
- Choose variable names that suggest the usage.
- The prefix ``numberOf`` should be used for variables representing a number of objects.
- The suffix ``Id`` should be used for variables representing an entity number.
- The prefix ``is`` or ``do`` or ``use`` or conjugations thereof should be used for boolean variables and methods.
- Negated boolean variable names must be avoided.

### Files

- C++ header files should have the extension ``.h``.
- Source files should have the extension ``.cpp``.
- Application files should have the extension ``.cpp``.
- Implementation files of templated classes should have the extension ``_Impl.h`` and be included at the end of the header file.
- Large source files can be split in multiple files by using the same name, then an ``_`` and a description of its usage.
- A class should be declared in a header file and defined in a source file where the name of the files match the name of the class.
- All definitions should reside in source files.
- Use spaces for indenting, the use of tab should be avoided.
- Header files must contain an include guard. The name convention resembles the location of the file inside the source tree. (all capitals i.e. BASE_ELEMENT_HPP)
- Include statements must be located at the top of a file only; except for template class implementations.

### Statements

- The parts of a class must be sorted public, protected and private.
- Within the parts, methods should go before variables.
- Use the ``explicit`` keyword for single argument constructors.
- Type conversions must always be done explicitly, never rely on implicit type conversion. The only exception is upcasting, which may be done implicitly or explicitly.
- For every ``new``, there must be a ``delete`` in the same class. If the ``delete`` is in another class, it must be documented carefully with both the ``new`` and ``delete``.
- Variables must never have dual meaning.
- Multiple variables should not be declared on the same line; except for the primitive data types e.g. ``int``.
- Use of global variables should be minimized.
- Class variables should never be declared public.
- C++ pointers and references should have their reference symbol next to the name rather than to the type, i.e. ``Element *element``.
- Implicit test for 0 should not be used other than for boolean variables.
- Variables should be declared in the smallest scope possible.
- Only loop control statements must be included in the ``for()`` construction.
- Range based for-loops are preferred over other types of for-loops. For example, a for-loop over vector vec can look like ``for ( double d: vec) { }``.
- Loop variables should be initialized immediately before the loop.
- Do-while loops can be avoided.
- The use of ``break`` and ``continue`` in loops should be avoided, go for readability.
- The form ``while (true)`` should be used for infinite loops.
- In cases of an if-statement, the nominal case should be put in the if-part and the exception in the else-part.
- Floating point constants should always be written with decimal point and at least one decimal. (i.e. ``x = 0.0``).
- Use ``nullptr`` instead of ``NULL`` or ``0``.
- Try to avoid smart pointers. If a smart pointer is necessary, use unique_ptr or shared_ptr.
- Use the type ``std::size_t`` for non-negative integer variables.
- Encapsulate global variables and constants, enumerated types, and type definitions in a class or namespace.
- The ``auto``-keyword for type declaration is only allowed if it is instantly clear what the type is.
- When using a template, use typename T and NOT class T.
- Avoid using ``typename`` as a return specifier.
- Try to avoid ``typedef`` for type definitions. Use the ``using``-keyword if either the type is likely to change, or the typename is very long. This should be done in the class declaration and not globally.
- Avoid returning by final function argument (for example ``void myFun (int a, int b, Matrix mat)`` for computing mat).



## About assertions

The intended use of assertions in general and ``logger.assert`` for hpGEM is to check blatant assumptions on the arguments provided to the
functions. For example, if you are trying to request entry 12 of an array of size 4, this is clearly wrong and you should change your function.

Assertions should be used freely, but checking everything, all the time, is expensive. Therefore assertions are compiled out of the code
completely in Release configuration. However, this also means that if you accidentally do something that has side effects in an assertion, the
behaviour of your code may change between configurations. For example, NEVER DO ``logger.assert(x=4, "/// WRONG!!");``

Preparing data for the assertions is also expensive. If you want to assert something that depends on data that is not yet computed and not
needed elsewhere in the function, do NOT compute it before the assertion. Rather write a dedicated function checkThisCondition that return
bool and call it inside the assertion. (make sure to think of a better name for the function!)


## About optimisation

Always assume compiler developers are better at optimising code than you are. They are hired specifically to optimise code. Moreover they
don't have to worry about the readability of the code they produce and they can use tricks that are not available in non-compiled c++.
Compiler developers have to make assumptions about the code they are going to process. These tend towards processing code that is easy
to read and write. Because of this, if you break these assumptions, they will instead assume you know what you are doing and are after a
specific effect. THIS MAY MEAN THAT WRITING CLEVER CODE BECAUSE YOU THINK IT IS FASTER PRODUCES WORST
PERFORMANCE THAN A NAIVE IMPLEMENTATION!

That being said, there are things that can be done to make an application faster, you just have to be careful that they actually make the
application faster.

Optimising code only makes sense if you are building in Release mode, before doing anything else, set ``CMAKE_BUILD_TYPE`` to ``Release``
and make sure ``hpGEM_ENFORCE_ASSERTIONS`` is ``OFF``. All effort spend in other configurations is almost certainly wasted and probably
even counterproductive. Make sure you have some timing information to compare before and after, to see if you made any improvements.

Before changing anything, you want to know where the bottlenecks are. To do this, use a profiler, such as ``gprof``(linux) or ``Instruments``(mac),
on a representative application. With this information in hand either make sure the most expensive methods get called less often or try to
improve the efficiency of this method. Note that saving information will most likely speed up your application, but comes at the cost of a larger
memory footprint. This may mean large applications no longer run. This is NOT an ideal situation.

Between steps and after you are done, make sure to compare with your old timings. If things went worse, discard your changes and try
again.

## About debugging

Debugging code only makes sense if you are building in Debug mode, before doing anything else, set ``CMAKE_BUILD_TYPE`` to ``Debug``. If
this changes the behaviour of your code, first toggle ``hpGEM_DISABLE_ASSERTIONS`` to see if someone made a mistake in an assertion. If
this did not alter the behaviour you have to hope a combination of ``CMAKE_BUILD_TYPE=RelWithDebInfo`` and
``hpGEM_ENFORCE_ASSERTIONS`` provides enough useful information

If you feel it is necessary to see additional debug output, you can set the advanced CMake option ``hpGEM_LOGLEVEL`` to ``DEBUG``. There is
quite a lot of additional debug info, so you probably want to pipe your output into grep or a similar text search tool.

If PETSc is reporting errors this usually means something went wrong in the interaction between hpGEM and PETSc. In this case it may be
useful to add the PETSc option -on_error_abort. This will make your code crash on the first PETSc error.

If you are working on a mac, entering ``env DYLD_INSERT_LIBRARIES=/usr/lib/libgmalloc.dylib`` in lldb before typing ``process launch`` will
enable lldb to trace memory corruption (for linux, valgrind can do this)

If something unexpected happens in the kernel of hpGEM this likely has one of the four following causes:

### Case 1: Nonsense arguments were passed to a function and therefore it returned nonsense

This happens for example when some function tries to request entry 12 of an array of size 4. This situation should be debugged by adding
assertions until the function successfully detects the erroneous arguments. In the example this means adding the line ``logger.assert(entry <
4, "This array only has 4 entries")``. Your code will now crash when 12 is passed to the array, making it more easy to find the places where
stuff actually went wrong. Other developers may have the same issue, please commit the extra assertions.

### Case 2: Some fields in a class contain nonsense and therefore it produces nonsense

This happens for example when someone constructs a face with a weird element. This means you are detecting one of the other two causes
in an unfortunate place. Check the other functions that interact with the erroneous field for the actual source of the error.

### Case 3: All arguments and the state of the class are valid, but it still breaks

This happens because of a programming error. Add a unit-test that detects the current problem and change things in the kernel until all unit
tests pass. Once they do please commit the bugfix.

### Case 4: Everything is fine, but there is an assertion or a unit test that fails

Unit-tests and assertions are great for detecting something is wrong, but unfortunately it may also be the unit-test or the assertion itself that is
wrong. Please be extra careful to make sure it is actually the unit-test or the assertion that is wrong and when in doubt ask. When you
commit the fix, clearly state in you commit message that you changed assertions or unit-tests because they were broken.

## FAQ

Q: Your 'performance optimisations' made the test suite much slower


A: The test suite tries to test parts of hpGEM in isolation. This means it sometimes has to take some non-standard actions to create or
manipulate data. It often also takes an action once for many different shapes/basisfunctions/etc. while a typical application takes this action
many times for the same shapes/basisfunctions/etc.. This means performance considerations for the test suite are different from performance
considerations for performance critical applications and improving speed for the latter may have adverse effect on the former.

Q: This bit of code seems like it is wasteful to me
A: Did you have a look at a profiler already? Did you read the section about performance optimisations already?

