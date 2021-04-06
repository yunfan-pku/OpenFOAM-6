/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "multispecie.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    template <class Thermo>
    Foam::multispecie<Thermo>::multispecie(
        const word& name,
        PtrList<Thermo>& speciesData,
        const speciesTable& specieNames,
        const dictionary& thermoDict)
        : N_(specieNames.size()),
        Y_(specieNames.size()),
        X_(specieNames.size()),
        speciesData_(speciesData),
        species_(specieNames) //,
                              //changed(1)
    {
    }

    template <class Thermo>
    Foam::multispecie<Thermo>::multispecie(const dictionary& thermoDict, PtrList<Thermo>& speciesData)
        : Y_(0),
        X_(0),
        speciesData_(0),
        species_(speciesData) //,
                              //changed(1)

    {
        wordList s(thermoDict.lookup("species"));
        //species_.transfer(s);
        N_ = species_.size();
        speciesData_.setSize(N_);
        Y_.setSize(N_);
        X_.setSize(N_);
        forAll(species_, i)
        {
            speciesData_.set(
                i,
                new Thermo(thermoDict.subDict(species_[i])));
        }
    }
    template <class Thermo>
    inline void Foam::multispecie<Thermo>::setY(const scalarList& Y)
    {
        //changed=1;
        scalar sum = 0;
        forAll(speciesData_, i)
        {
            Y_[i] = Y[i];
            sum += Y[i] / speciesData_[i].W();
        }
        forAll(speciesData_, i)
        {
            X_[i] = Y[i] / speciesData_[i].W() / sum;
        }
    }

    template <class Thermo>
    inline void Foam::multispecie<Thermo>::setX(const scalarList& X)
    {
        //changed=1;
        scalar sum = 0, sn = 0;

        forAll(speciesData_, i)
        {
            X_[i] = X[i];
            sum += X[i] * speciesData_[i].W();
            sn += X[i];
        }
        forAll(speciesData_, i)
        {
            Y_[i] = X[i] * speciesData_[i].W() / sum;
            X_[i] /= sn;
        }
    }

    // * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
    /*
    template<class Thermo>
    inline multispecie<Thermo>::multispecie(
        const word &name,
        const scalar Y,
        const scalar molWeight)
        : name_(name),
          Y_(Y),
          molWeight_(molWeight)
    {
    }
    */
    /*
template<class Thermo>
    inline multispecie<Thermo>::multispecie(
        const scalar Y,
        const scalar molWeight)
        : Y_(Y),
          molWeight_(molWeight)
    {
    }

    forAll(species_, i)
    {
        speciesData_.set(
            i,
            new Foam::thermotable<ThermoType>(thermoDict.subDict(species_[i])));
    }
    */

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
    /*
    inline multispecie::multispecie(const word &name, const multispecie &st)
        : name_(name),
          Y_(st.Y_),
          molWeight_(st.molWeight_)
    {
    }
*/
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
    template <class Thermo>
    inline const Thermo& Foam::multispecie<Thermo>::operator[](const label i) const
    {
        return speciesData_[i];
    }
    template <class Thermo>
    inline const PtrList<Thermo>& Foam::multispecie<Thermo>::speciesData() const
    {
        return speciesData_;
    }

    template <class Thermo>
    inline const speciesTable& Foam::multispecie<Thermo>::species() const
    {
        return species_;
    }

    template <class Thermo>
    inline const scalarList& Foam::multispecie<Thermo>::Y() const
    {
        return Y_;
    }

    template <class Thermo>
    inline const scalarList& Foam::multispecie<Thermo>::X() const
    {
        return X_;
    }
    /*
    inline const word &multispecie::name() const
    {
        return name_;
    }

    inline scalar multispecie::W() const
    {
        return molWeight_;
    }

    inline scalar multispecie::Y() const
    {
        return Y_;
    }

    inline scalar multispecie::R() const
    {
        return RR / molWeight_;
    }
*/
// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //
/*
inline void multispecie::operator=(const multispecie &st)
{
    // name_ = st.name_;
    Y_ = st.Y_;
    molWeight_ = st.molWeight_;
}

inline void multispecie::operator+=(const multispecie &st)
{
    const scalar sumY = Y_ + st.Y_;
    if (mag(sumY) > small)
    {
        molWeight_ = sumY / (Y_ / molWeight_ + st.Y_ / st.molWeight_);
    }

    Y_ = sumY;
}

inline void multispecie::operator*=(const scalar s)
{
    Y_ *= s;
}
*/
// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //
/*
inline multispecie operator+(const multispecie &st1, const multispecie &st2)
{
    const scalar sumY = max(st1.Y_ + st2.Y_, small);

    if (mag(sumY) > small)
    {
        return multispecie(
            sumY,
            sumY / (st1.Y_ / st1.molWeight_ + st2.Y_ / st2.molWeight_));
    }
    else
    {
        return st1;
    }
}

inline multispecie operator*(const scalar s, const multispecie &st)
{
    return multispecie(
        s * st.Y_,
        st.molWeight_);
}

inline multispecie operator==(const multispecie &st1, const multispecie &st2)
{
    scalar diffY = st2.Y_ - st1.Y_;
    if (mag(diffY) < small)
    {
        diffY = small;
    }

    const scalar diffRW =
        st2.Y_ / st2.molWeight_ - st1.Y_ / st1.molWeight_;

    scalar molWeight = great;
    if (mag(diffRW) > small)
    {
        molWeight = diffY / diffRW;
    }

    return multispecie(diffY, molWeight);
}
*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
