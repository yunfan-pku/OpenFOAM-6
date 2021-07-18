/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::ISATbinaryTree

Description

    Data storage of the data according to a binary
    tree structure.

            0 (root node)
         /     \
        0       0
      /   \   /   \
     L     R L     0
    / \
   L   R

    L: leafLeft_
    R: leafRight_

\*---------------------------------------------------------------------------*/

#ifndef ISATbinaryTree_H
#define ISATbinaryTree_H


#include "dictionary.H"
//#include "ISATleaf.H"
#include "ISATNode.H"
#include "List2D.H"
//#include "binaryNode.H"
//#include "chemPointISAT.H"

namespace Foam
{


    //template<class CompType, class ThermoType>
    //class TDACChemistryModel;

    class ISATleaf;

    class ISATNode;


    /*
        class Heap
        {
        private:
            List<ISATleaf*> heap_;
            label size_;
            label maxsize_;
        public:
            Heap(label n) :heap_(n), size_(0), maxsize_(n) {}
            void reset(label n) { heap_.resize(n);size_ = 0;maxsize_ = n; }
            void insert(ISATleaf* pin) { if (size_ >= maxsize_)FatalErrorInFunction << "size_ >=maxsize_" << exit(FatalError);heap_[size_++] = pin;shiftUp(size_ - 1); }
            void shiftUp(int start) {
                int j = start, i = (j - 1) / 2; ISATleaf* temp = heap_[j];
                while (j > 0) {
                    if (heap_[i]->numRetrieve() >= temp->numRetrieve()) break;
                    else { heap_[j] = heap_[i];  j = i; i = (j - 1) / 2; }
                }
                heap_[j] = temp;
            }

            void ShiftDown(int start) {
                int i = start, j = 2 * i + 1; ISATleaf* temp = heap_[i];
                while (j < size_) {
                    if (j + 1 < size_ && heap_[j]->numRetrieve() < heap_[j + 1]->numRetrieve()) j++;
                    if (temp->numRetrieve() >= heap_[j]->numRetrieve())  break;
                    else { heap_[i] = heap_[j]; i = j; j = 2 * j + 1; }
                }
                heap_[i] = temp;
            }
            ISATleaf* pop() { if (size_ <= 0)FatalErrorInFunction << "size_<=0" << exit(FatalError);ISATleaf* temp = heap_[0];heap_[0] = heap_[--size_];ShiftDown(0);return temp; }
        };*/

        //template<class CompType, class ThermoType>

    class ISATbinaryTree
    {

    public:
        //typedef binaryNode<CompType, ThermoType> bn;
        //typedef chemPointISAT<CompType, ThermoType> chP;

    private:
        label n_in_;
        label n_out_;

        //- Reference to the chemistryModel
       // TDACChemistryModel<CompType, ThermoType>& chemistry_;

        //- Root node of the binary tree
        ISATNode* root_;

        //- Maximum number of elements in the binary tree
        label maxNLeafs_;

        //- Size of the BST (= number of chemPoint stored)
        label size_;

        List2D<ISATleaf> timeTagList_;

        //- Secondary retrieve search variables
        //label n2ndSearch_;
        //label max2ndSearch_;

        //- Insert new node at the position of phi0
        //  phi0 should be already attached to another node or the pointer to it
        //  will be lost

        void insertNode
        (
            ISATleaf*& phi0,
            ISATNode*& newNode
        );

        ISATleaf* Sibling(ISATleaf* x);
        ISATNode* nodeSibling(ISATleaf* x);
        //- Replace the binaryNode u with v
        void transplant(ISATNode* u, ISATNode* v);
        //- Perform a search in the subtree starting from the subtree node y
        //  This search continue to use the hyperplan to walk the tree
        //  If covering EOA is found return true and x points to the chemPoint
        /*bool inSubTree
        (
            const scalarField& phiq,
            bn* y,
            chP* x
        );
        */
        /*







        chP* chemPSibling(bn* y);


        bn* nodeSibling(bn* y);

        void deleteAllNode(bn* subTreeRoot);

        dictionary coeffsDict_;
        */
    public:
        //- Constructors

            //- Construct from dictionary and chemistryOnLineLibrary
        ISATbinaryTree
        (
            //TDACChemistryModel<CompType, ThermoType>& chemistry,
            //dictionary coeffsDict
            label in_n,
            label out_n,
            label maxNLeafs,
            label timeTagList
        ) :n_in_(in_n), n_out_(out_n), root_(nullptr), maxNLeafs_(maxNLeafs), size_(0), timeTagList_(timeTagList) {};



        inline ISATNode* root()
        {
            return root_;
        }

        inline const ISATNode* root()const
        {
            return root_;
        }


        inline  label n_in()const
        {
            return n_in_;
        }
        inline  label n_out()const
        {
            return n_out_;
        }

        inline  List2D<ISATleaf>& timeTagList()
        {
            return timeTagList_;
        }
        void deleteSubTree(ISATNode* subTreeRoot);

        inline void deleteSubTree()
        {
            deleteSubTree(root_);
        }


        inline void deleteSubTree_Node()
        {
            deleteSubTree_Node(root_);
        }

        void deleteSubTree_Node(ISATNode* subTreeRoot);
        
        inline label size() const
        {
            return size_;
        }
        //- Computes iteratively the depth of the subTree
        label depth(ISATNode* subTreeRoot)const;

        inline label depth() const
        {
            return depth(root_);
        }
        void deleteLeaf(ISATleaf* pleaf);


        //- Member functions



        inline label maxNLeafs()
        {
            return maxNLeafs_;
        }
        inline void setmaxNLeafs(label n)
        {
            maxNLeafs_ = n;
        }


        // Insert a new leaf starting from the parent node of phi0
        // Parameters: phi0 the leaf to replace by a node
        // phiq the new composition to store
        // Rphiq the mapping of the new composition point
        // A the mapping gradient matrix
        // B the matrix used to initialize the EOA
        // nCols the size of the matrix
        // Returns: void
        // Description :
        //1) Create a new leaf with the data to initialize the EOA and to
        // retrieve the mapping by linear interpolation (the EOA is
        // initialize in the chemPoint constructor)
        //2) Get the parent node of phi0 and connect a new node in place of the
        // leaf of phi0. This new node is constructed with phi0 on the left
        // and phiq on the right (the hyperplane is computed inside the
        // binaryNode constructor)
        ISATleaf* insertNewLeaf
        (
            const scalarList& value,
            const scalarList& data
        );



        // Search the ISATbinaryTree until the nearest leaf of a specified
        // leaf is found.
        void binaryTreeSearch
        (
            const scalarList& value,
            ISATNode* node,
            ISATleaf*& nearest
        );

        void eval(const scalarList& value, scalarList& ret);
        /*void print
        (
            Ostream& OFout,
            ISATNode* node = nullptr,
            int a = 0
        );
        */

        //- Removes every entries of the tree and delete the associated objects
        void clear();

        ~ISATbinaryTree();

        // Perform a secondary binary tree search starting from a failed
        // chemPoint x, with a depth-first search algorithm
        // If another candidate is found return true and x points to the chemP
        //bool secondaryBTSearch(const scalarField& phiq, chP*& x);

        //- Delete a leaf from the binary tree and reshape the binary tree for
        //  the following binary tree search
        //  Return the index in the nodeList of the removed node
        //  (-1 when no node)
        //void deleteLeaf(chP*& phi0);

        //- Cheap balance function
        //  This function just roughly separate the space in two parts
        //  with a hyperplane which separate the two extreme chemPoint in the
        //  direction of the maximum the variance
        //  Then, it repopulate the tree with this hyperplane stored at the root
        //  and by inserting the chemPoint in increasing order of value in that
        //  direction



        ISATNode* balance_build(scalarRectangularMatrix& scaleIn, ISATNode* pNode, int start, int end, List<ISATleaf*>& arrays);
        void balance(scalarRectangularMatrix& scaleIn);
        //void balance();
        void getMiddle(List<ISATleaf*>& arrays, int ic, int left, int right, int posleft, int posright);
        //inline void balance(scalarRectangularMatrix & scaleIn);
/*
        inline void deleteAllNode()
        {
            deleteAllNode(root_);
        }

        chP* treeMin(bn* subTreeRoot);

        inline chP* treeMin()
        {
            return treeMin(root_);
        }

        chP* treeSuccessor(chP* x);





        void resetNumRetrieve();
        */

        //- ListFull
        bool isFull();


        template <typename T>
        friend class ISATmanager;
    };
    /*
        Ostream& operator<<(Ostream& OFout, const ISATbinaryTree& tree)
        {
            return  OFout << (*tree.root());

        }
        */


        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//#include "ISATbinaryTree.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
