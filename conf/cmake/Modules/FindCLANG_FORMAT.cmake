# Find clang-format
#
# CLANG_FORMAT_EXECUTABLE   - Path to clang-format executable
# CLANG_FORMAT_FOUND        - True if the clang-format executable was found.
# CLANG_FORMAT_VERSION      - The version of clang-format found
#
# Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# from http://github.com/votca/votca/blob/master/CMakeModules/FindCLANG_FORMAT.cmake

find_program(CLANG_FORMAT_EXECUTABLE
             NAMES clang-format-11
                   clang-format-10
	               clang-format
                   clang-format-9
                   clang-format-8
                   clang-format-7

             DOC "clang-format executable")
mark_as_advanced(CLANG_FORMAT_EXECUTABLE)

message("Clang format executable found at ${CLANG_FORMAT_EXECUTABLE}")
# Extract version from command "clang-format -version"
if(CLANG_FORMAT_EXECUTABLE)
  execute_process(COMMAND ${CLANG_FORMAT_EXECUTABLE} -version
                  OUTPUT_VARIABLE clang_format_version
                  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
  message("Found clang-format reports version ${clang_format_version}")
  if(clang_format_version MATCHES "^.*clang-format version .*")
    # clang_format_version sample: "clang-format version 3.9.1-4ubuntu3~16.04.1
    # (tags/RELEASE_391/rc2)"
    string(REGEX
           REPLACE "^.*clang-format version ([.0-9]+).*"
                   "\\1"
                   CLANG_FORMAT_VERSION
                   "${clang_format_version}")
    # CLANG_FORMAT_VERSION sample: "3.9.1"
  else()
    set(CLANG_FORMAT_VERSION 0.0)
  endif()
else()
  set(CLANG_FORMAT_VERSION 0.0)
endif()

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set CLANG_FORMAT_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(CLANG_FORMAT REQUIRED_VARS CLANG_FORMAT_EXECUTABLE VERSION_VAR CLANG_FORMAT_VERSION)
