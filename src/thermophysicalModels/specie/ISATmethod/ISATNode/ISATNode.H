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
    Foam::ISATNode

Description
    Node of the binary tree

SourceFile
    ISATNode.C

\*---------------------------------------------------------------------------*/

#ifndef ISATNode_H
#define ISATNode_H

#include "ISATleaf.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                               Class ISATNode Declaration
    \*---------------------------------------------------------------------------*/


    class ISATNode
    {

    public:

        //- Element on the left
        ISATleaf* leafLeft_;

        //- Element on the right
        ISATleaf* leafRight_;

        //- Node which follows on the left
        ISATNode* nodeLeft_;

        //- Node which follows on the right
        ISATNode* nodeRight_;

        //- Parent node
        ISATNode* parent_;

        //- Number of equations in addition to the species eqs.
        //label nAdditionalEqns_;

        scalarList v_;

        scalar a_;

        //- Compute vector v:
        //  Let E be the ellipsoid which covers the region of accuracy of
        //  the left leaf (previously defined). E is described by
        //  E={phi| ||L^T.(phi-phi0)|| <= 1}, (see chemPoint for more details).
        //  let E' be the transformation of E in a space where E' is a hypersphere
        //  centered at the origin, in this space y=L^T.(phi-phi0) and then
        //  E'={y| ||y||<=1}
        //  let u be the unit vector joining the center of E' and the newly added
        //  composition point in the transformed space
        //  (y2=L^T.(phiq-phi0)),u = y2/||y2|
        //  Then the hyperplane separating the two points is defined as the
        //  perpendicular bisector of the segment linking 0 to y2
        //      H' = {y| u^T.(y-yh) = 0},
        //  where yh = y2/2.
        //  In the original composition space, the hyperplane H is defined by
        //      H = {y| v^T(phi-phih) = 0},
        //  where phih = phi0 + L^-T.yh = (phi0 + phiq) / 2 and v is
        //            L.L^T (phiq-phi0)
        //      v = -------------------- .
        //          ||L.L^T (phiq-phi0)||
        //
        //  Note :  v is not normalised in this implementation since it is used
        //              on both side of an equality to know if one should go on the
        //              left or the right in the binary tree
        //  Parameter:
        //      elementLeft : chemPoint of the left element
        //      elementRight: chemPoint of the right element
        //       v : empty scalar field to store v
        //  Returnq: void (v is stored in the empty scalarField)
        void calcV
        (
            ISATleaf*& elementLeft,
            ISATleaf*& elementRight,
            scalarList& v
        );

        //- Compute a the product v^T.phih, with phih = (phi0 + phiq)/2.
        //  When travelling in the binary tree,
        //  to know in which part of the composition space the query point 'phi'
        //  belongs to, v^T.phi is performed. If the result is "> a" then it belongs
        //  to the right part (where phiq is), otherwise it belongs to the left
        // part (where phi0 is).
        scalar calcA
        (
            ISATleaf* elementLeft,
            ISATleaf* elementRight
        );

        // Constructors

            //- Construct null
        ISATNode();

        //- Construct from components
        ISATNode
        (
            ISATleaf* elementLeft,
            ISATleaf* elementRight,
            ISATNode* parent = nullptr
        );

        ISATNode
        (
            ISATleaf* elementLeft,
            ISATleaf* elementRight,
            ISATNode* parent,
            int dir
        );


        // Member functions

            //- Access

        inline ISATleaf*& leafLeft()
        {
            return leafLeft_;
        }

        inline ISATleaf*& leafRight()
        {
            return leafRight_;
        }

        inline ISATNode*& nodeLeft()
        {
            return nodeLeft_;
        }

        inline ISATNode*& nodeRight()
        {
            return nodeRight_;
        }

        inline const  ISATleaf* leafLeft() const
        {
            return leafLeft_;
        }

        inline const ISATleaf* leafRight() const
        {
            return leafRight_;
        }

        inline const ISATNode* nodeLeft() const
        {
            return nodeLeft_;
        }

        inline const  ISATNode* nodeRight() const
        {
            return nodeRight_;
        }

        inline ISATNode*& parent()
        {
            return parent_;
        }

        //- Topology

        inline const scalarList& v() const
        {
            return v_;
        }

        inline scalarList& v()
        {
            return v_;
        }

        inline const scalar& a() const
        {
            return a_;
        }

        inline scalar& a()
        {
            return a_;
        }
    };
    /*
    Foam::Ostream& operator<<(Foam::Ostream& OFout, const Foam::ISATNode& node)
    {
        if (nullptr != node.leafLeft())// phi0 is on the right
        {
            OFout<<(*node.leafLeft());
        }
        if (nullptr != node.nodeLeft())// phi0 is on the right
        {
            OFout<<(*node.nodeLeft());
        }
        if (nullptr != node.leafRight())// phi0 is on the right
        {
            OFout<<*node.leafRight();
        }
        if (nullptr != node.nodeRight())// phi0 is on the right
        {
            OFout<<*node.nodeRight();
        }
        return OFout;
    }
    */

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//    #include "ISATNode.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
