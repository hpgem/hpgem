/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2017, University of Twente
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef HPGEM_ELEMENTSHAPE_H
#define HPGEM_ELEMENTSHAPE_H

#include<cstddef>
#include<vector>
#include "Logger.h"
#include "customVector.h"

namespace Preprocessor {

    template<std::size_t dimension>
    class ElementShape;

    namespace Detail {

        template<std::size_t entityDimension, std::size_t dimension>
        struct EntityData: public EntityData<entityDimension - 1, dimension> {
            template<typename... superArgs>
            EntityData(std::vector<stackVector<std::size_t>> adjacentNodes,
                       std::vector<const ElementShape<entityDimension>*> entityShapes, superArgs... args) :
                       EntityData<entityDimension - 1, dimension>(args...), adjacentShapes(), entityShapes(entityShapes) {
                adjacentShapes[0] = adjacentNodes;
            }

            EntityData() = default;

            ~EntityData() = default;

            EntityData(const EntityData&) = default;

            EntityData(EntityData&&) noexcept = default;

            EntityData& operator=(const EntityData&) = default;

            EntityData& operator=(EntityData&&) noexcept = default;

            std::array<std::vector<stackVector<std::size_t>>, dimension> adjacentShapes;
            std::vector<const ElementShape<entityDimension>*> entityShapes;
            
            template<std::size_t SUB_DIM>
            std::vector<const ElementShape<SUB_DIM>*>& getEntityShapes()
            {
                return EntityData<SUB_DIM, dimension>::entityShapes;
            }

            template<std::size_t SUB_DIM>
            const std::vector<const ElementShape<SUB_DIM>*>& getEntityShapes() const
            {
                return EntityData<SUB_DIM, dimension>::entityShapes;
            }

            template<std::size_t SUB_DIM>
            std::array<std::vector<stackVector<std::size_t>>, dimension>& getAdjacentShapes()
            {
                return this->EntityData<SUB_DIM, dimension>::adjacentShapes;
            }

            template<std::size_t SUB_DIM>
            const std::array<std::vector<stackVector<std::size_t>>, dimension>& getAdjacentShapes() const
            {
                return EntityData<SUB_DIM, dimension>::adjacentShapes;
            }
        };

        template<std::size_t dimension>
        struct EntityData<0, dimension> {
            EntityData(std::vector<stackVector<std::size_t>> adjacentNodes,
                       std::vector<const ElementShape<0>*> entityShapes):
                       adjacentShapes(), entityShapes(entityShapes) {
                adjacentShapes[0] = adjacentNodes;
            }

            EntityData() = default;

            ~EntityData() = default;

            EntityData(const EntityData&) = default;

            EntityData(EntityData&&) noexcept = default;

            EntityData& operator=(const EntityData&) = default;

            EntityData& operator=(EntityData&&) noexcept = default;

            std::array<std::vector<stackVector<std::size_t>>, dimension> adjacentShapes;
            std::vector<const ElementShape<0>*> entityShapes;

            template<std::size_t SUB_DIM>
            std::vector<const ElementShape<SUB_DIM>*>& getEntityShapes()
            {
                return EntityData<SUB_DIM, dimension>::entityShapes;
            }

            template<std::size_t SUB_DIM>
            std::array<std::vector<stackVector<std::size_t>>, dimension>& getAdjacentShapes()
            {
                return this->EntityData<SUB_DIM, dimension>::adjacentShapes;
            }

            template<std::size_t SUB_DIM>
            const std::vector<const ElementShape<SUB_DIM>*>& getEntityShapes() const
            {
                return EntityData<SUB_DIM, dimension>::entityShapes;
            }

            template<std::size_t SUB_DIM>
            const std::array<std::vector<stackVector<std::size_t>>, dimension>& getAdjacentShapes() const
            {
                return this->EntityData<SUB_DIM, dimension>::adjacentShapes;
            }
        };
    }

    template<std::size_t dimension>
    class ElementShape {
        template<std::size_t d>
        struct tag{};
    public:
        ElementShape() = default;

        ~ElementShape() = default;

        ElementShape(const ElementShape&) = default;

        ElementShape(ElementShape&&) noexcept = default;

        ElementShape& operator=(const ElementShape&) = default;

        ElementShape& operator=(ElementShape&&) noexcept = default;

        template<typename... Args>
        ElementShape(Args... args) : shapeData(args...) {
            completeSubShapes(tag<dimension - 1>{}, tag<dimension - 1>{});
            logger.assert(checkShape(), "Input generated an inconsistent shape");
        }


        std::size_t getNumberOfNodes() const {
            return getNumberOfEntities<0>();
        }

        std::size_t getNumberOfEdges() const {
            return getNumberOfEntities<1>();
        }

        std::size_t getNumberOfFaces() const {
            return getNumberOfEntities<-1>();
        }

        //if (entityDimension >= 0) ElementShape<entityDimension>
        //else if(entityDimension + dimension >= 0) ElementShape<entityDimension + dimension> //codim case
        //else ElementShape<0> //make the compiler not crash while we invoke the logger
        template<int entityDimension>
        using ShapeType = const ElementShape<
        (entityDimension < 0 ? (entityDimension + dimension < 0 ? 0 : entityDimension + dimension) : entityDimension)>;

        ShapeType<-1> getFaceShape(std::size_t faceIndex) const {
            return getBoundaryShape<-1>(faceIndex);
        }

        template<int entityDimension>
        stackVector<std::size_t> getNodesOfEntity(std::size_t entityIndex) const {
            return getAdjacentEntities<entityDimension, 0>(entityIndex);
        }

        /**
         * @brief returns the number of entities of the specified dimension. Negative dimensions are treated as codimensions
         */
        template<int entityDimension>
        std::enable_if_t<(entityDimension < 0), std::size_t> getNumberOfEntities() const;
        template<int entityDimension>
        std::enable_if_t<(entityDimension >= dimension), std::size_t> getNumberOfEntities() const;
        template<int entityDimension>
        std::enable_if_t<(entityDimension >= 0 && entityDimension < dimension), std::size_t> getNumberOfEntities() const;

        template<int entityDimension>
        std::enable_if_t<(entityDimension < 0), const ShapeType<entityDimension>*> getBoundaryShape(std::size_t entityIndex) const;
        template<int entityDimension>
        std::enable_if_t<(entityDimension >= dimension), const ShapeType<entityDimension>*> getBoundaryShape(std::size_t entityIndex) const;
        template<int entityDimension>
        std::enable_if_t<(entityDimension >= 0 && entityDimension < dimension), const ShapeType<entityDimension>*> getBoundaryShape(std::size_t entityIndex) const;

        /**
         * @brief returns the indices of adjacent entities in the index space of this shape. If the target dimension is
         * negative, it will be taken with respect to the full shape
         */
        template<int entityDimension, int targetDimension>
        std::enable_if_t<(entityDimension < 0), stackVector<std::size_t>> getAdjacentEntities(std::size_t entityIndex) const;
        template<int entityDimension, int targetDimension>
        std::enable_if_t<(targetDimension < 0 && entityDimension >= 0), stackVector<std::size_t>> getAdjacentEntities(std::size_t entityIndex) const;

        template<int entityDimension, int targetDimension>
        std::enable_if_t<(entityDimension >= 0 && targetDimension >= 0 && (entityDimension >= dimension || targetDimension >= dimension)), stackVector<std::size_t>>
        getAdjacentEntities(std::size_t entityIndex) const;

        template<int entityDimension, int targetDimension>
        std::enable_if_t<(entityDimension >= 0 && entityDimension < dimension && targetDimension >= 0 && targetDimension < dimension), stackVector<std::size_t>>
        getAdjacentEntities(std::size_t entityIndex) const;

        bool checkShape() const {
            return checkBoundaryShape(tag<dimension - 1>{});
        }

    private:

        template<std::size_t d, std::size_t shapeDimension>
        bool checkSingleShape(const ElementShape<shapeDimension>* boundingShape, stackVector<std::size_t> boundaryNodes, std::size_t currentIndex, tag<d>) const;

        template<std::size_t shapeDimension>
        bool checkSingleShape(const ElementShape<shapeDimension>* boundingShape, stackVector<std::size_t> boundaryNodes, std::size_t currentIndex, tag<0>) const;

        template<std::size_t d>
        bool checkBoundaryShape(tag<d>) const;

        bool checkBoundaryShape(tag<0>) const;

        template<std::size_t entityDimension, std::size_t targetDimension>
        void completeSubShapes(tag<entityDimension>, tag<targetDimension>);

        template<std::size_t entityDimension>
        void completeSubShapes(tag<entityDimension>, tag<0>);

        void completeSubShapes(tag<0>, tag<0>) {}

        Detail::EntityData<dimension - 1, dimension> shapeData;

    };

    template<>
    class ElementShape<0> {
    public:
        ElementShape() = default;

        ~ElementShape() = default;

        ElementShape(const ElementShape&) = default;

        ElementShape(ElementShape&&) noexcept = default;

        ElementShape& operator=(const ElementShape&) = default;

        ElementShape& operator=(ElementShape&&) noexcept = default;


        std::size_t getNumberOfNodes() const {
            return getNumberOfEntities<0>();
        }

        std::size_t getNumberOfEdges() const {
            return getNumberOfEntities<1>();
        }

        std::size_t getNumberOfFaces() const {
            return getNumberOfEntities<-1>();
        }

        //if (entityDimension >= 0) ElementShape<entityDimension>
        //else if(entityDimension + dimension >= 0) ElementShape<entityDimension + dimension> //codim case
        //else ElementShape<0> //make the compiler not crash while we invoke the logger
        template<int entityDimension>
        using ShapeType = ElementShape<0>;

        ShapeType<-1> getFaceShape(std::size_t faceIndex) const {
            return getBoundaryShape<-1>(faceIndex);
        }

        template<int entityDimension>
        stackVector<std::size_t> getNodesOfEntity(std::size_t entityIndex) const {
            return getAdjacentEntities<entityDimension, 0>(entityIndex);
        }

        /**
         * @brief returns the number of entities of the specified dimension. Negative dimensions are treated as codimensions
         */
        template<int entityDimension>
        std::size_t getNumberOfEntities() const {
            if(entityDimension == 0) return 1; else return 0;
        }

        template<int entityDimension>
        ShapeType<entityDimension> getBoundaryShape(std::size_t entityIndex) const {
            logger.assert(entityDimension == 0, "A point is not bounded by shapes of other dimensions");
            logger.assert(entityIndex == 0, "A point consists of only 1 shape, but you asked for shape %", entityIndex);
            return *this;
        }

        /**
         * @brief returns the indices of adjacent entities in the index space of this shape. If the target dimension is
         * negative, it will be taken with respect to the full shape
         */
        template<int entityDimension, int targetDimension>
        stackVector<std::size_t> getAdjacentEntities(std::size_t entityIndex) const {
            logger.assert(entityDimension == 0, "A point is not bounded by shapes of other dimensions");
            logger.assert(entityIndex == 0, "A point consists of only 1 shape, but you asked for shape %", entityIndex);
            if(targetDimension == 0) return {0}; else return {};
        };

        //points are hardcoded to be correct
        bool checkShape() const{
            return true;
        }
    };

    //note when debugging new shapes: invoking the logger during static initialisation can be a bit messy
    //so you may need to use a debugger to see the error message
    const ElementShape<0> point{};
    const ElementShape<1> line(std::vector<stackVector<std::size_t>>{{0}, {1}}, std::vector<const ElementShape<0>*>{&point, &point});
    const ElementShape<2> triangle(std::vector<stackVector<std::size_t>>{{0, 1}, {0, 2}, {1, 2}}, std::vector<const ElementShape<1>*>{&line, &line, &line},
                                   std::vector<stackVector<std::size_t>>{{0}, {1}, {2}}, std::vector<const ElementShape<0>*>{&point, &point, &point});
    const ElementShape<2> square(std::vector<stackVector<std::size_t>>{{0, 1}, {0, 2}, {1, 3}, {2, 3}}, std::vector<const ElementShape<1>*>{&line, &line, &line, &line},
                                 std::vector<stackVector<std::size_t>>{{0}, {1}, {2}, {3}}, std::vector<const ElementShape<0>*>{&point, &point, &point, &point});
    const ElementShape<3> tetrahedron(std::vector<stackVector<std::size_t>>{{0, 3, 2}, {0, 1, 3}, {0, 2, 1}, {1, 2, 3}}, std::vector<const ElementShape<2>*>{&triangle, &triangle, &triangle, &triangle},
                                      std::vector<stackVector<std::size_t>>{{0, 1}, {0, 2}, {0, 3}, {2, 3}, {1, 3}, {1, 2}}, std::vector<const ElementShape<1>*>{&line, &line, &line, &line, &line, &line},
                                      std::vector<stackVector<std::size_t>>{{0}, {1}, {2}, {3}}, std::vector<const ElementShape<0>*>{&point, &point, &point, &point});
    const ElementShape<3> cube(std::vector<stackVector<std::size_t>>{{0, 1, 2, 3}, {0, 1, 4, 5}, {0, 2, 4, 6}, {1, 3, 5, 7}, {2, 3, 6, 7}, {4, 5, 6, 7}}, std::vector<const ElementShape<2>*>{&square, &square, &square, &square, &square, &square},
                               std::vector<stackVector<std::size_t>>{{0, 1}, {2, 3}, {4, 5}, {6, 7}, {0, 2}, {1, 3}, {4, 6}, {5, 7}, {0, 4}, {1, 5}, {2, 6}, {3, 7}}, std::vector<const ElementShape<1>*>{&line, &line, &line, &line, &line, &line, &line, &line, &line, &line, &line, &line},
                               std::vector<stackVector<std::size_t>>{{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}}, std::vector<const ElementShape<0>*>{&point, &point, &point, &point, &point, &point, &point, &point});
    const ElementShape<3> triangularPrism(std::vector<stackVector<std::size_t>>{{0, 2, 1}, {3, 4, 5}, {2, 0, 5, 3}, {0, 1, 3, 4}, {1, 2, 4, 5}}, std::vector<const ElementShape<2>*>{&triangle, &triangle, &square, &square, &square},
                                          std::vector<stackVector<std::size_t>>{{0, 1}, {0, 2}, {1, 2}, {3, 4}, {3, 5}, {4, 5}, {0, 3}, {1, 4}, {2, 5}}, std::vector<const ElementShape<1>*>{&line, &line, &line, &line, &line, &line, &line, &line, &line},
                                          std::vector<stackVector<std::size_t>>{{0}, {1}, {2}, {3}, {4}, {5}}, std::vector<const ElementShape<0>*>{&point, &point, &point, &point, &point, &point});
    const ElementShape<3> pyramid(std::vector<stackVector<std::size_t>>{{3, 4, 1, 2}, {3, 1, 0}, {2, 4, 0}, {1, 2, 0}, {4, 3, 0}}, std::vector<const ElementShape<2>*>{&square, &triangle, &triangle, &triangle, &triangle},
                                  std::vector<stackVector<std::size_t>>{{0, 1}, {0, 2}, {0, 3}, {0, 4}, {1, 2}, {2, 4}, {4, 3}, {3, 1}}, std::vector<const ElementShape<1>*>{&line, &line, &line, &line, &line, &line, &line, &line},
                                  std::vector<stackVector<std::size_t>>{{0}, {1}, {2}, {3}, {4}}, std::vector<const ElementShape<0>*>{&point, &point, &point, &point, &point});

    template<std::size_t dimension>
    const std::vector<const ElementShape<dimension>*> defaultShapes;

    template<>
    const std::vector<const ElementShape<0>*> defaultShapes<0> = {&point};

    template<>
    const std::vector<const ElementShape<1>*> defaultShapes<1> = {&line};

    template<>
    const std::vector<const ElementShape<2>*> defaultShapes<2> = {&triangle, &square};

    template<>
    const std::vector<const ElementShape<3>*> defaultShapes<3> = {&tetrahedron, &cube, &triangularPrism, &pyramid};

}


#include<elementShape_impl.h>

#endif //HPGEM_ELEMENTSHAPE_H
