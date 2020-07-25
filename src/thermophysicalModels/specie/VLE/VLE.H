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
    Foam::VLE

Description
    VLE mixture equation of state.

SourceFiles
    VLEI.H
    VLE.C

\*---------------------------------------------------------------------------*/

#ifndef VLE_H
#define VLE_H

#include "autoPtr.H"
//#include "scalarList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // Forward declaration of friend functions and operators

    template <class ThermoMixture>
    class VLE;

    template <class ThermoMixture>
    inline VLE<ThermoMixture> operator==(
        const VLE<ThermoMixture> &,
        const VLE<ThermoMixture> &);

    template <class ThermoMixture>
    inline VLE<ThermoMixture> operator*(
        const scalar s,
        const VLE<ThermoMixture> &);

    template <class ThermoMixture>
    Ostream &operator<<(
        Ostream &,
        const VLE<ThermoMixture> &);

    /*---------------------------------------------------------------------------*\
                           Class VLE Declaration
\*---------------------------------------------------------------------------*/

    template <class ThermoMixture>
    class VLE
        : public ThermoMixture
    {

    public:
        typedef typename ThermoMixture::SingleThermoType SingleThermoType;
        // Constructors

        //- Construct from components
        inline VLE(
            const ThermoMixture &sp,
            const scalar &Tc,
            const scalar &Vc,
            const scalar &Pc,
            const scalar &omega);

        //- Construct from dictionary
        VLE(const dictionary &dict, PtrList<ThermoMixture> &speciesData);

        //- Construct as named copy
        inline VLE(const word &name, const VLE &);

        inline VLE(const word &name, PtrList<SingleThermoType> &, const speciesTable &, const dictionary &);
        //- Construct and return a clone
        inline autoPtr<VLE> clone() const;

        // Selector from dictionary
        inline static autoPtr<VLE> New(
            const dictionary &dict);

        // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "VLE" + word(ThermoMixture::typeName());
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
        inline scalar B(scalar p, scalar T) const;
        inline scalar Ha(scalar p, scalar T) const;
        inline scalar Hc() const;
        inline scalar Hs(scalar p, scalar T) const;
        inline scalar dGdT(scalar p, scalar T) const;
        inline scalar dCpdT(scalar p, scalar T) const;
        inline scalar limit(scalar T) const;
        //solve cubic equatoin
        static inline scalar cubicSolve(scalar a0, scalar a1, scalar a2, label flag); //flag 0 smallest>0 ,flag 1 largest

        struct solution
        {
            autoPtr<scalarList> X_liq, X_gas, equalconstant;
            scalar vaporfra;
            solution() : X_liq(), X_gas(), equalconstant(), vaporfra(0){};
            solution(const solution &ref) : X_liq(ref.X_liq), X_gas(ref.X_gas), equalconstant(ref.equalconstant), vaporfra(ref.vaporfra){};
        };
        inline scalar vaporfra(scalar p, scalar T) const;
        inline scalar alpha(scalar p, scalar T, solution &sol) const;
        inline label solve_vaporf(scalar p, scalar T, solution &sol) const;
        inline label noramlize(scalarList &X_in) const;
        inline autoPtr<scalarList> wilsoneq(scalar p, scalar T) const;
        inline autoPtr<scalarList> myinitial(scalar p, scalar T) const;
        autoPtr<solution> TPn_flash(scalar p, scalar T) const;
        void  TPn_validation(scalar p, scalar T) const;
        

         //mutable solution temp_sol;
         //mutable scalar temp_p,temp_T;
        /* friend VLE operator==<ThermoMixture>(
            const VLE &,
            const VLE &);

        friend VLE operator*<ThermoMixture>(
            const scalar &,
            const VLE &);
            */
        // IO

        //- Write to Ostream
        void write(Ostream &os) const;

        // Ostream Operator

        friend Ostream &operator<<<ThermoMixture>(
            Ostream &,
            const VLE &);
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "VLEI.H"

#ifdef NoRepository
#include "VLE.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
