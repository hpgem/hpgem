//--------------------------------------
// Copyright 2013 M.T. Julianto
//
#ifndef LevelTree_hpp
#define LevelTree_hpp
//------------------------------------------------------------------------------
// System includes and names imported from them:
#include <list>
#include <vector>
#include <iostream>
//------------------------------------------------------------------------------
// Package includes:
#include "TreeBase.hpp"
#include "TreeIterator.hpp"
#include "TreeEntry.hpp"

//------------------------------------------------------------------------------


namespace Base 
{

    template <class V>
    class LevelTree
    {

    public:
        typedef V valueT;
        typedef TreeEntry<V> treeEntryT;
        //  this is our (rather special) iterator 
        typedef Base::TreeIterator<V,TreeEntry<V> > iterator;
        typedef unsigned int DimT;

        LevelTree();

        ~LevelTree();

        bool empty() const;

        //! Number of entries in the LevelTree
        int size() const;

        int maxLevel() const;

      //protected:
        void setActiveLevel(const DimT level);

        void resetActiveLevel();
        
        int getActiveLevel() const;

      public:
        //! Getting the beginning of traversal
        iterator begin();
        
        //! Getting the end of traversal
        iterator end();

        //! Add new entry
        iterator addEntry(const valueT& newEl, const bool preserveLinks=false);
        
        //! Add new dummy TreeEntry
        iterator addDummyEntry(iterator it);

        //! Add children of an entry.
        iterator addChildren(iterator parentEl, const std::vector<valueT>& subEntries);

        //! Set the current entries as part of the coarsest mesh.
        void setAsTheCoarsestEntries();

        //! Getting the beginning of tree-level traversal
        iterator beginLevel(const int level);

        //! Getting the beginning of leaves traversal
        iterator beginLeaf();

        //! Getting the beginning of pre-order traversal
        iterator beginPreOrder();

        //! Getting the beginning of post-order traversal
        iterator beginPostOrder();

        //! Erase all descendants of an entry
        void eraseChildsOf(iterator fci);

        //! Describe the LevelTree
        friend std::ostream& operator<<(std::ostream& os, const LevelTree<V>& e)
        {
            os<< "LevelTree: ";
            os<< "noRootEntries_= "<< e.noRootEntries_<< " ";
            os<< "entries_.size()= "<< e.entries_.size() << " ";
            os<< "minLevel_= "<< e.minLevel_<< " ";
            os<< "maxLevel_= "<< e.maxLevel_<< " ";
            os<< "activeLevel_= "<< e.activeLevel_<< " ";
            os<< "coarsestEntriesSet_= "<< e.coarsestEntriesSet_<< " ";
            return os;
        }


    private:
        //! Add new TreeEntry
        iterator addTreeEntry(treeEntryT* const newEnt, const bool preserveLinks=false);

        //! Add tree-children of an entry.
        iterator addTreeChildren(iterator parentEl, const std::vector<treeEntryT>& subEntries);

        //! Erase a leaf entry
        iterator eraseLastLeaf(iterator fci);

    private:
        int                     noRootEntries_;
        std::list<treeEntryT*>  entries_;
        int                     minLevel_;
        int                     maxLevel_;
        int                     activeLevel_;
        bool                    coarsestEntriesSet_;
    };

} // close namespace Base

// merge the implementation file here.
#include "Base/LevelTree_Impl.hpp"

#endif //LevelTree_hpp
