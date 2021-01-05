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

#ifndef HPGEM_KERNEL_VTKSPECIFICTIMEWRITER_H
#define HPGEM_KERNEL_VTKSPECIFICTIMEWRITER_H

#include <functional>
#include <fstream>
#include "Base/MeshManipulator.h"
namespace hpgem {
namespace Output {

///\class VTKSpecificTimeWriter
///\brief writes steady state data or a single time level
///
/// this produces multiple files, you do not have to append anything, just load
/// the .pvtu into paraview VTK makes the assumption that all data is 3D data,
/// this class will provide conversions where necessary, but not from 4D to 3D
///
/// Implementation note: The output can in general be discontinuous at the
/// element boundaries. To handle such a function requires a bit of special
/// care. The VTK fileformat associates data values with cells (elements) or
/// points. If we would store a single function value f(x) for a node at x, then
/// the result is necessarily continuous. After all, independently from which
/// element one approaches x, the value will be f(x).
/// To represent discontinuous data we will therefore need to duplicate the node
/// x to x_E1, x_E2, ... x_EN for the N attached elements. This way we can store
/// multiple values f(x_E1), f(x_E2), ..., f(x_EN) and represent data that is
/// discontinuous at the element boundaries.
///
/// \tparam DIM The dimension of the mesh
// class is final because the destructor would be the only virtual function
template <std::size_t DIM>
class VTKSpecificTimeWriter final {
   public:
    ///\brief write front matter and open file stream
    ///\param baseName name of the file WITHOUT extentions
    ///\param mesh the mesh containing the data you want to output
    /// if you want to write from multiple meshes, simply have paraview load
    /// both output files
    VTKSpecificTimeWriter(const std::string& baseName,
                          const Base::MeshManipulator<DIM>* mesh,
                          std::size_t timelevel = 0);

    ///\brief write end matter and close the file stream
    ~VTKSpecificTimeWriter();

    ///\brief write a scalar field
    void write(
        std::function<double(
            Base::Element*, const Geometry::PointReference<DIM>&, std::size_t)>,
        const std::string& name);

    ///\brief write a vector field
    void write(
        std::function<LinearAlgebra::SmallVector<DIM>(
            Base::Element*, const Geometry::PointReference<DIM>&, std::size_t)>,
        const std::string& name);

    ///\brief write an order 2 tensor field
    void write(
        std::function<LinearAlgebra::SmallMatrix<DIM, DIM>(
            Base::Element*, const Geometry::PointReference<DIM>&, std::size_t)>,
        const std::string& name);

    ///\brief do not copy the writer to prevent havoc when destructing all the
    /// copies
    VTKSpecificTimeWriter(const VTKSpecificTimeWriter& orig) = delete;
    VTKSpecificTimeWriter operator=(const VTKSpecificTimeWriter& orig) = delete;

   private:
    /**
     * Write the start of the master pvtu file, this will contain the
     * indirection to the data in the local vtu files.
     * @param baseName The base name for the file (actual file will have .pvtu
     * extension)
     */
    void writeMasterFileHeader(const std::string& baseName);
    /**
     * Write the start of the processor local vtu file.
     * @param baseName The base name for the file (actual files will have .[proc
     * num].vtu extension)
     */
    void writeLocalFileHeader(const std::string& baseName);
    /**
     * Write binary data for a <DataArray> to the localFile stream.
     *
     * The data will be encoded in base64 and will be given the appropriate
     * length header.
     * @param data A pointer to the data.
     * @param byteLength The length in bytes to write
     */
    void writeLocalFileBinaryData(void* data, std::uint32_t byteLength);

    std::ofstream localFile_;
    std::ofstream masterFile_;
    /// Number of points in the local file
    std::uint32_t totalPoints_;
    /// Number of elements in the local file
    std::uint32_t totalElements_;
    const Base::MeshManipulator<DIM>* mesh_;
    std::size_t timelevel_;
};
}  // namespace Output
}  // namespace hpgem
#include "VTKSpecificTimeWriter_Impl.h"

#endif  // HPGEM_KERNEL_VTKSPECIFICTIMEWRITER_H
