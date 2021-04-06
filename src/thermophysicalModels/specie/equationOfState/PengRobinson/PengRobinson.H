/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
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
    Foam::PengRobinson

Description
    PengRobinson mixture equation of state.

SourceFiles
    PengRobinsonI.H
    PengRobinson.C

\*---------------------------------------------------------------------------*/

#ifndef PengRobinson_H
#define PengRobinson_H

#include "autoPtr.H"
//#include "scalarList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // Forward declaration of friend functions and operators

    template <class Specie>
    class PengRobinson;

    template <class Specie>
    inline PengRobinson<Specie> operator+(
        const PengRobinson<Specie>&,
        const PengRobinson<Specie>&);

    template <class Specie>
    inline PengRobinson<Specie> operator*(
        const scalar,
        const PengRobinson<Specie>&);

    template <class Specie>
    inline PengRobinson<Specie> operator==(
        const PengRobinson<Specie>&,
        const PengRobinson<Specie>&);

    template <class Specie>
    Ostream& operator<<(
        Ostream&,
        const PengRobinson<Specie>&);

    /*---------------------------------------------------------------------------*\
                           Class PengRobinson Declaration
\*---------------------------------------------------------------------------*/

    template <class Specie>
    class PengRobinson
        : public Specie
    {
    private:
    public:
        //- Critical Temperature [K]
        scalar Tc_;
        //- Critical volume [m^3/kmol]
        scalar Vc_;
        //- Critical compression factor [-]
        scalar Zc_;
        //- Critical Pressure [Pa]
        scalar Pc_;
        //- Acentric factor [-]
        scalar omega_;

        scalarList Hig_phase_;
        scalarList Hig2_phase_;

    public:
        label phaseFlag; //gas 1,liq 0

    public:
        // Constructors
        typedef PengRobinson<Specie> EOS;
        //- Construct from components
        inline PengRobinson(
            const Specie& sp,
            const scalar& Tc,
            const scalar& Vc,
            const scalar& Pc,
            const scalar& omega);

        //- Construct from dictionary
        PengRobinson(const dictionary& dict);

        //- Construct as named copy
        inline PengRobinson(const word& name, const PengRobinson&);

        //- Construct and return a clone
        inline autoPtr<PengRobinson> clone() const;

        // Selector from dictionary
        inline static autoPtr<PengRobinson> New(
            const dictionary& dict);

        // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "PengRobinson<" + word(Specie::typeName_()) + '>';
        }

        // Fundamental properties

        //- Is the equation of state is incompressible i.e. rho != f(p)
        static const bool incompressible = false;

        //- Is the equation of state is isochoric i.e. rho = const
        static const bool isochoric = false;

        //- Return density [kg/m^3]
        inline scalar rho(scalar p, scalar T) const;

        //- Return enthalpy departure [J/kg]
        inline scalar H(const scalar p, const scalar T) const;

        //- Return Cp departure [J/(kg K]
        inline scalar Cp(scalar p, scalar T) const;

        //- Return entropy [J/(kg K)]
        inline scalar S(const scalar p, const scalar T) const;

        //- Return compressibility rho/p [s^2/m^2]
        inline scalar psi(scalar p, scalar T) const;

        //- Return compression factor []
        inline scalar Z(scalar p, scalar T) const;

        //- Return (Cp - Cv) [J/(kg K]
        inline scalar CpMCv(scalar p, scalar T) const;

        inline scalar A(scalar p, scalar T) const;
        inline scalar dAdT(scalar p, scalar T) const;
        inline scalar dAdP(scalar p, scalar T) const;
        inline scalar B(scalar p, scalar T) const;
        inline scalar dBdT(scalar p, scalar T) const;
        inline scalar dBdP(scalar p, scalar T) const;
        inline scalar a(scalar p, scalar T) const;
        inline scalar b(scalar p, scalar T) const;
        inline scalar dadt(scalar p, scalar T) const;
        inline scalar dadp(scalar p, scalar T) const;
        inline scalar dbdT(scalar p, scalar T) const;
        inline scalar dbdP(scalar p, scalar T) const;
        inline scalar d2adt2(scalar p, scalar T) const;
        inline scalar EOSalpha(scalar P, scalar T) const;

        //Temp
        inline scalar Hig_phase(scalar T) const;
        inline scalar Hig2_phase(scalar T) const;
        //solve cubic equatoin
        static inline scalar cubicSolve(scalar a0, scalar a1, scalar a2, label flag); //flag 0 smallest>0 ,flag 1 largest
        // IO

        //- Write to Ostream
        void write(Ostream& os) const;

        // Member operators

        inline void operator+=(const PengRobinson&);
        inline void operator*=(const scalar);

        // Friend operators

        friend PengRobinson operator+<Specie>(
            const PengRobinson&,
            const PengRobinson&);

        friend PengRobinson operator*<Specie>(
            const scalar s,
            const PengRobinson&);

        friend PengRobinson operator==<Specie>(
            const PengRobinson&,
            const PengRobinson&);

        // Ostream Operator

        friend Ostream& operator<<<Specie>(
            Ostream&,
            const PengRobinson&);
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PengRobinsonI.H"

#ifdef NoRepository
#include "PengRobinson.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
