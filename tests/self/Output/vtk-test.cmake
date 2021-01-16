# This file forms part of hpGEM. This package has been developed over a number of
# years by various people at the University of Twente and a full list of
# contributors can be found at http://hpgem.org/about-the-code/team
#
# This code is distributed using BSD 3-Clause License. A copy of which can found
# below.
#
# Copyright (c) 2021, University of Twente
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors
# may be used to endorse or promote products derived from this software without
# specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Testing helper script to generate output and compare the result against some
# known good files. This is combined so that it can be used as single command to
# add_test().
#
# This script expects two inputs variables
# generator_cmd: A command to generate the output, e.g. $<TARGET_FILE:[application]>
# reference_path: A path with reference files. Each files is compared against the output.

# Path to the executable generating the output
if (NOT generator_cmd)
    message(FATAL_ERROR "No generator command")
endif ()
# Path with references
if (NOT reference_path)
    message(FATAL_ERROR "No path with references")
endif ()

# Run command to get
execute_process(COMMAND ${generator_cmd} RESULT_VARIABLE CMD_RESULT)
if (CMD_RESULT)
    message(FATAL_ERROR "Error running ${generator_cmd}")
endif ()

# Compare with output
file(GLOB reference_files "${reference_path}/*")
if (NOT reference_files)
    message(FATAL_ERROR "No reference files")
endif ()
# Check correctness of each file
foreach (reference_file ${reference_files})
    get_filename_component(test_file ${reference_file} NAME)
    execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files ${reference_file} ${test_file}
            RESULT_VARIABLE test_not_successful OUTPUT_QUIET ERROR_QUIET)
    if (test_not_successful)
        message(SEND_ERROR "${test_file} differs from ${reference_file}")
    endif ()
endforeach ()