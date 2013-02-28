# CMake generated Testfile for 
# Source directory: /Users/antRthorn/Dropbox/hpgem/Phoenix/tests/unit/LinearAlgebra
# Build directory: /Users/antRthorn/Dropbox/hpgem/Phoenix/tests/unit/LinearAlgebra
# 
# This file includes the relevent testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
ADD_TEST(BuildTest "/Applications/CMake 2.8-10.app/Contents/bin/cmake" "--build" "/Users/antRthorn/Dropbox/hpgem/Phoenix" "--target" "MatrixTest.out")
ADD_TEST(MatrixTest "MatrixTest.out")
SET_TESTS_PROPERTIES(MatrixTest PROPERTIES  DEPENDS "BuildTest")
ADD_TEST(VectorTest "VectorTest.out")
