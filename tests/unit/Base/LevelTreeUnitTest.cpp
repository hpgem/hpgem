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

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>
#include <typeinfo>
#include "Base/LevelTree.h"
#include "Base/TreeIterator.h"
#include "Base/TreeEntry.h"

#define CATCH_CONFIG_MAIN
#include "../catch.hpp"

using namespace hpgem;
//========================================================= Here is the main
// function
TEST_CASE("LevelTreeUnitTest", "[LevelTreeUnitTest]") {
    Base::LevelTree<std::size_t> testBalancedBinary, testSingleRoot,
        testMultiRoot;

    INFO("initial tree nonempty");
    CHECK(testBalancedBinary.empty());
    INFO("initial tree nonempty");
    CHECK(testSingleRoot.empty());
    INFO("initial tree nonempty");
    CHECK(testMultiRoot.empty());
    INFO("initial tree nonempty");
    CHECK(testBalancedBinary.size() == 0);
    INFO("initial tree nonempty");
    CHECK(testSingleRoot.size() == 0);
    INFO("initial tree nonempty");
    CHECK(testMultiRoot.size() == 0);
    INFO("initial tree nonempty");
    CHECK(testBalancedBinary.maxLevel() == 0);
    INFO("initial tree nonempty");
    CHECK(testSingleRoot.maxLevel() == 0);
    INFO("initial tree nonempty");
    CHECK(testMultiRoot.maxLevel() == 0);

    for (std::size_t i : testBalancedBinary) {
        logger(ERROR, "should be able to iterator over an empty tree");
    }

    testBalancedBinary.addRootEntry(0);
    testSingleRoot.addRootEntry(0);
    for (std::size_t i = 0; i < 5; ++i) {
        testMultiRoot.addRootEntry(i);
    }
    INFO("tree should have entries now");
    CHECK(!testBalancedBinary.empty());
    INFO("tree should have entries now");
    CHECK(!testSingleRoot.empty());
    INFO("tree should have entries now");
    CHECK(!testMultiRoot.empty());
    INFO("tree should have entries now");
    CHECK(testBalancedBinary.size() == 1);
    INFO("tree should have entries now");
    CHECK(testSingleRoot.size() == 1);
    INFO("tree should have entries now");
    CHECK(testMultiRoot.size() == 5);
    INFO("tree should only have root entries now");
    CHECK(testBalancedBinary.maxLevel() == 1);
    INFO("tree should only have root entries now");
    CHECK(testSingleRoot.maxLevel() == 1);
    INFO("tree should only have root entries now");
    CHECK(testMultiRoot.maxLevel() == 1);

    testBalancedBinary.setAllLevelTraversal();
    testSingleRoot.setTraversalMethod(Base::TreeTraversalMethod::ALLLEVEL);
    Base::TreeIterator<std::size_t> iterator = testBalancedBinary.begin();
    for (std::size_t i = 1; i < 31; i += 2) {
        testBalancedBinary.addChild(iterator, i);
        testBalancedBinary.addChild(iterator, i + 1);
        ++iterator;
    }
    iterator = testSingleRoot.begin();
    testSingleRoot.addChildren(iterator, std::vector<std::size_t>{1, 2, 3});
    iterator++;
    testSingleRoot.addChildren(iterator, std::vector<std::size_t>{14, 15}, 3);
    iterator++;
    testSingleRoot.addChild(iterator, 4);
    testSingleRoot.addChildren(iterator, std::vector<std::size_t>{5, 6});
    iterator++;
    testSingleRoot.addChildren(iterator, std::vector<std::size_t>{7, 8, 9}, 2);
    testSingleRoot.addChildren(iterator, std::vector<std::size_t>{10, 11, 12});
    testSingleRoot.addChild(iterator, 13, 2);
    iterator++;
    iterator++;
    testSingleRoot.addChildren(iterator, std::vector<std::size_t>(1, 16));
    iterator++;
    iterator++;
    testSingleRoot.addChildren(iterator, std::vector<std::size_t>{17, 18});
    while (*iterator < 15) {
        ++iterator;
    }
    for (std::size_t i = 19; i < 31; ++i) {
        testSingleRoot.addChild(iterator, i);
    }
    iterator = testMultiRoot.beginAllLevel();
    testMultiRoot.addChildren(iterator, std::vector<std::size_t>{5, 6, 7});
    iterator++;
    testMultiRoot.addChild(iterator, 8);
    iterator++;
    testMultiRoot.addChildren(iterator, std::vector<std::size_t>{9, 10});
    iterator++;
    iterator++;
    testMultiRoot.addChild(iterator, 11);
    iterator++;
    iterator++;
    iterator++;
    testMultiRoot.addChild(iterator, 12, 3);
    iterator++;
    for (std::size_t i = 13; i < 21; ++i) {
        testMultiRoot.addChild(iterator, i);
    }
    iterator++;
    for (std::size_t i = 21; i < 25; ++i) {
        testMultiRoot.addChild(iterator, i);
    }
    iterator++;
    for (std::size_t i = 25; i < 31; ++i) {
        testMultiRoot.addChild(iterator, i);
    }
    INFO("tree should have entries now");
    CHECK(!testBalancedBinary.empty());
    INFO("tree should have entries now");
    CHECK(!testSingleRoot.empty());
    INFO("tree should have entries now");
    CHECK(!testMultiRoot.empty());
    INFO("tree should have entries now");
    CHECK(testBalancedBinary.size() == 1);
    INFO("tree should have entries now");
    CHECK(testSingleRoot.size() == 1);
    INFO("tree should have entries now");
    CHECK(testMultiRoot.size() == 5);
    INFO("tree should have multiple levels now");
    CHECK(testBalancedBinary.maxLevel() == 5);
    INFO("tree should have multiple levels now");
    CHECK(testSingleRoot.maxLevel() == 5);
    INFO("tree should have multiple levels now");
    CHECK(testMultiRoot.maxLevel() == 4);

    std::vector<std::size_t> binaryPostOrder = {
        15, 16, 7,  17, 18, 8, 3,  19, 20, 9,  21, 22, 10, 4, 1, 23,
        24, 11, 25, 26, 12, 5, 27, 28, 13, 29, 30, 14, 6,  2, 0};
    std::vector<std::size_t> binaryPreOrder = {
        0, 1, 3,  7,  15, 16, 8,  17, 18, 4,  9,  19, 20, 10, 21, 22,
        2, 5, 11, 23, 24, 12, 25, 26, 6,  13, 27, 28, 14, 29, 30};
    std::vector<std::size_t> singlePostOrder = {
        14, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 15, 1, 1,
        16, 4,  5,  17, 18, 6,  2,  7,  8,  9,  10, 11, 12, 13, 3, 0};
    std::vector<std::size_t> singlePreOrder = {
        0,  1, 1, 14, 15, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
        30, 2, 4, 16, 5,  6,  17, 18, 3,  7,  8,  9,  10, 11, 12, 13};
    std::vector<std::size_t> singleAllLevel = {
        0,  1,  2,  3,  1,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14,
        15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};
    std::vector<std::size_t> multiPostOrder = {
        5,  6,  12, 7,  7, 0,  13, 14, 15, 16, 17, 18, 19, 20, 8,  1,
        21, 22, 23, 24, 9, 25, 26, 27, 28, 29, 30, 10, 2,  3,  11, 4};
    std::vector<std::size_t> multiPreOrder = {
        0, 5, 6,  7,  7,  12, 1,  8,  13, 14, 15, 16, 17, 18, 19, 20,
        2, 9, 21, 22, 23, 24, 10, 25, 26, 27, 28, 29, 30, 3,  4,  11};
    std::vector<std::size_t> multiAllLevel = {
        0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 7,  13, 14, 15,
        16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 12};

    std::size_t i = 0;
    for (std::size_t test : testBalancedBinary) {
        INFO("ordering of the entries got messed up");
        CHECK(test == i);
        i++;
    }
    INFO("amount of entries got messed up");
    CHECK(i == 31);
    i = 7;
    testBalancedBinary.setSingleLevelTraversal(3);
    for (std::size_t test : testBalancedBinary) {
        INFO("entries on level 3 got messed up");
        CHECK(test == i);
        i++;
    }
    INFO("entries on level 3 got messed up");
    CHECK(i == 15);
    i = 14;
    testSingleRoot.setSingleLevelTraversal(3);
    for (std::size_t test : testSingleRoot) {
        INFO("entries on level 3 got messed up");
        CHECK(test == i);
        i++;
    }
    INFO("entries on level 3 got messed up");
    CHECK(i == 19);
    i = 12;
    testMultiRoot.setSingleLevelTraversal(3);
    for (std::size_t test : testMultiRoot) {
        INFO("entries on level 3 got messed up");
        CHECK(test == i);
        i++;
    }
    INFO("entries on level 3 got messed up");
    CHECK(i == 13);
    // range-for also works with the other traversal methods, but we have to do
    // a loop over two ranges to test for sanity
    auto testIterator = binaryPreOrder.begin();
    for (iterator = testBalancedBinary.beginPreOrder();
         iterator != testBalancedBinary.end(); ++iterator, ++testIterator) {
        INFO("pre-ordering of the entries got messed up");
        CHECK(*iterator == *testIterator);
    }
    testSingleRoot.setPreOrderTraversal();
    testIterator = singlePreOrder.begin();
    for (iterator = testSingleRoot.begin();
         testIterator != singlePreOrder.end(); ++iterator, ++testIterator) {
        INFO("pre-ordering of the entries got messed up");
        CHECK(*iterator == *testIterator);
    }
    testMultiRoot.setTraversalMethod(Base::TreeTraversalMethod::PREORDER);
    testIterator = multiPreOrder.begin();
    for (iterator = testMultiRoot.begin(); iterator != testMultiRoot.end();
         ++iterator, ++testIterator) {
        INFO("pre-ordering of the entries got messed up");
        CHECK(*iterator == *testIterator);
    }
    testBalancedBinary.setTraversalMethod(Base::TreeTraversalMethod::POSTORDER);
    testIterator = binaryPostOrder.begin();
    for (iterator = testBalancedBinary.begin();
         testIterator != binaryPostOrder.end(); ++iterator, ++testIterator) {
        INFO("post-ordering of the entries got messed up");
        CHECK(*iterator == *testIterator);
    }
    testIterator = singlePostOrder.begin();
    for (iterator = testSingleRoot.beginPostOrder();
         iterator != testSingleRoot.end(); ++iterator, ++testIterator) {
        INFO("post-ordering of the entries got messed up");
        CHECK(*iterator == *testIterator);
    }
    testMultiRoot.setPostOrderTraversal();
    testIterator = multiPostOrder.begin();
    for (iterator = testMultiRoot.begin(); testIterator != multiPostOrder.end();
         ++iterator, ++testIterator) {
        INFO("post-ordering of the entries got messed up");
        CHECK(*iterator == *testIterator);
    }
    testIterator = singleAllLevel.begin();
    for (iterator = testSingleRoot.beginAllLevel();
         iterator != testSingleRoot.end(); ++iterator, ++testIterator) {
        INFO("level-ordering of the entries got messed up");
        CHECK(*iterator == *testIterator);
    }
    testMultiRoot.setAllLevelTraversal();
    testIterator = multiAllLevel.begin();
    for (iterator = testMultiRoot.begin(); testIterator != multiAllLevel.end();
         ++iterator, ++testIterator) {
        INFO("level-ordering of the entries got messed up");
        CHECK(*iterator == *testIterator);
    }
    // and again for reverse iteration
    auto testReverseIterator = binaryPreOrder.rbegin();
    testBalancedBinary.setPreOrderTraversal();
    auto reverseIterator = testBalancedBinary.rbegin();
    for (; reverseIterator != testBalancedBinary.rend();
         ++reverseIterator, ++testReverseIterator) {
        INFO("pre-ordering of the entries got messed up");
        CHECK(*reverseIterator == *testReverseIterator);
    }
    testSingleRoot.setPreOrderTraversal();
    testReverseIterator = singlePreOrder.rbegin();
    for (reverseIterator = testSingleRoot.rbegin();
         reverseIterator != testSingleRoot.rend();
         ++reverseIterator, ++testReverseIterator) {
        INFO("pre-ordering of the entries got messed up");
        CHECK(*reverseIterator == *testReverseIterator);
    }
    testMultiRoot.setTraversalMethod(Base::TreeTraversalMethod::PREORDER);
    testReverseIterator = multiPreOrder.rbegin();
    for (reverseIterator = testMultiRoot.rbegin();
         reverseIterator != testMultiRoot.rend();
         ++reverseIterator, ++testReverseIterator) {
        INFO("pre-ordering of the entries got messed up");
        CHECK(*reverseIterator == *testReverseIterator);
    }
    testBalancedBinary.setTraversalMethod(Base::TreeTraversalMethod::POSTORDER);
    testReverseIterator = binaryPostOrder.rbegin();
    for (reverseIterator = testBalancedBinary.rbegin();
         reverseIterator != testBalancedBinary.rend();
         ++reverseIterator, ++testReverseIterator) {
        INFO("post-ordering of the entries got messed up");
        CHECK(*reverseIterator == *testReverseIterator);
    }
    testReverseIterator = singlePostOrder.rbegin();
    testSingleRoot.setPostOrderTraversal();
    for (reverseIterator = testSingleRoot.rbegin();
         reverseIterator != testSingleRoot.rend();
         ++reverseIterator, ++testReverseIterator) {
        INFO("post-ordering of the entries got messed up");
        CHECK(*reverseIterator == *testReverseIterator);
    }
    testMultiRoot.setPostOrderTraversal();
    testReverseIterator = multiPostOrder.rbegin();
    for (reverseIterator = testMultiRoot.rbegin();
         reverseIterator != testMultiRoot.rend();
         ++reverseIterator, ++testReverseIterator) {
        INFO("post-ordering of the entries got messed up");
        CHECK(*reverseIterator == *testReverseIterator);
    }
    i = 30;
    testBalancedBinary.setAllLevelTraversal();
    for (reverseIterator = testBalancedBinary.rbegin();
         reverseIterator != testBalancedBinary.rend(); ++reverseIterator, --i) {
        INFO("level-ordering of the entries got messed up");
        CHECK(*reverseIterator == i);
    }

    // test some special cases
    testSingleRoot.fillToLevel(4);
    std::vector<std::size_t> level3 = {14, 15, 16, 5,  17, 18, 7,
                                       8,  9,  10, 11, 12, 13};
    std::vector<std::size_t> level4 = {14, 19, 20, 21, 22, 23, 24, 25,
                                       26, 27, 28, 29, 30, 16, 5,  17,
                                       18, 7,  8,  9,  10, 11, 12, 13};

    INFO("tree should still have entries");
    CHECK(!testSingleRoot.empty());
    INFO("tree should still have entries");
    CHECK(testSingleRoot.size() == 1);
    INFO("tree should still have multiple levels");
    CHECK(testSingleRoot.maxLevel() == 5);
    testSingleRoot.setSingleLevelTraversal(3);
    testIterator = level3.begin();
    for (iterator = testSingleRoot.begin(); iterator != testSingleRoot.end();
         ++iterator, ++testIterator) {
        INFO("level-ordering of the entries got messed up");
        CHECK(*iterator == *testIterator);
    }
    std::cout << std::endl;
    // leaf traversal
    testSingleRoot.setSingleLevelTraversal(4);
    testIterator = level4.begin();
    for (iterator = testSingleRoot.begin(); iterator != testSingleRoot.end();
         ++iterator, ++testIterator) {
        INFO("level-ordering of the entries got messed up");
        CHECK(*iterator == *testIterator);
    }
    testMultiRoot.setAllLevelTraversal();
    iterator = testMultiRoot.begin();
    while (*iterator != 8) {
        ++iterator;
    }
    testMultiRoot.eraseChilds(iterator);
    ++iterator;
    testMultiRoot.eraseChilds(iterator);
    ++iterator;
    testMultiRoot.eraseChilds(iterator);
    --iterator;
    --iterator;
    --iterator;
    INFO("tree should still have entries");
    CHECK(!testMultiRoot.empty());
    INFO("tree should still have entries");
    CHECK(testMultiRoot.size() == 5);
    INFO("tree should still have multiple levels");
    CHECK(testMultiRoot.maxLevel() == 4);
    testMultiRoot.eraseChilds(iterator);
    INFO("tree should still have entries");
    CHECK(!testMultiRoot.empty());
    INFO("tree should still have entries");
    CHECK(testMultiRoot.size() == 5);
    INFO("tree should have less levels");
    CHECK(testMultiRoot.maxLevel() == 2);
}
