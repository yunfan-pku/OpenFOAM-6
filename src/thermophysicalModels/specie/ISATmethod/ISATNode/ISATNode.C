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

\*---------------------------------------------------------------------------*/

#include "ISATNode.H"
#include "IOstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ISATNode::ISATNode()
    : leafLeft_(nullptr),
    leafRight_(nullptr),
    nodeLeft_(nullptr),
    nodeRight_(nullptr),
    parent_(nullptr)
{
}

Foam::ISATNode::ISATNode(
    ISATleaf* elementLeft,
    ISATleaf* elementRight,
    ISATNode* parent)
    : leafLeft_(elementLeft),
    leafRight_(elementRight),
    nodeLeft_(nullptr),
    nodeRight_(nullptr),
    parent_(parent),
    v_(0)
{


    calcV(elementLeft, elementRight, v_);
    scalarList v(v_.size());
    a_ = calcA(elementLeft, elementRight);
    forAll(v, i)
    {
        v[i] = 0;
    }
    v[0] = 1;
    scalar s = 0;
    forAll(v, i)
    {
        s += v[0] * v_[0];
    }

    if (s < 0)
    {
        calcV(elementRight, elementLeft, v_);
        a_ = calcA(elementRight, elementLeft);
        leafLeft_ = elementRight;
        leafRight_ = elementLeft;

    }


}


Foam::ISATNode::ISATNode(
    ISATleaf* elementLeft,
    ISATleaf* elementRight,
    ISATNode* parent,
    int dir)
    : leafLeft_(nullptr),
    leafRight_(nullptr),
    nodeLeft_(nullptr),
    nodeRight_(nullptr),
    parent_(parent),
    v_(0)
{
    v_.resize(elementRight->value().size());

    for (int i = 0;i < v_.size();i++)
        v_[i] = 0;
    v_[dir] = 1.0;
    a_ = calcA(elementLeft, elementRight);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ISATNode::calcV(
    ISATleaf*& elementLeft,
    ISATleaf*& elementRight,
    scalarList& v) {
    scalarList d = elementRight->value() - elementLeft->value();

    v.resize(d.size());
    scalar sum = 0;
    forAll(elementRight->value(), i)
    {
        sum += d[i] * d[i];
    }
    v = d / sqrt(sum);
}

Foam::scalar Foam::ISATNode::calcA(
    ISATleaf* elementLeft,
    ISATleaf* elementRight)
{
    scalarList temp = (elementLeft->value() + elementRight->value()) / 2;
    scalar sum = 0;
    forAll(temp, i)
    {
        sum += temp[i] * v_[i];
    }
    return sum;
}



// ************************************************************************* //
