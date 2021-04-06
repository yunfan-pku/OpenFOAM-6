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

Class
    Foam::species::ISATmultithermo

Description
    Basic thermodynamics type based on the use of fitting functions for
    cp, h, s obtained from the template argument type ISATmultithermo.  All other
    properties are derived from these primitive functions.

SourceFiles
    ISATmultithermoI.H
    ISATmultithermo.C

\*---------------------------------------------------------------------------*/

#ifndef ISATmultithermo_H
#define ISATmultithermo_H

#include "thermodynamicConstants.H"
#include "ISATmanager.H"
#include <tuple>
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    //TODO finish this!!!!
/*
    struct ISATfunction
    {
        ISATmanager* p;
        typedef double (ISATmanager::* VLEfunc)(double, double);
        VLEfunc pfunc;
        ISATfunction(ISATmanager p_in, VLEfunc pfunc_in) :p(p_in), pfunc(pfunc_in) {}
        void  value(const scalarList& in, scalarList& out)
        {
            out.resize(1);
            out[0] = (*p).*pfunc(in);//+2000*::Foam::sin(in[1]/10);
        }
        void  derive(const scalarList& in, scalarRectangularMatrix& out)
        {
            scalarList temp = in;
            scalarList tempv;
            value(in, tempv);
            scalarList tempv2;
            scalarList dd;
            scalar dx = 1e-5;
            for (int i = 0;i < in.size();i++)
            {
                temp[i] += dx;
                value(temp, tempv2);
                dd = (tempv2 - tempv) / dx;
                for (int j = 0;j < tempv2.size();j++)
                {
                    out[j][i] = dd[j];
                }
                temp[i] -= dx;
            }
        }
    };
    */


    namespace species
    {

        // Forward declaration of friend functions and operators

        template<class ThermoMixture, template<class> class Type> class ISATmultithermo;

        template<class ThermoMixture, template<class> class Type>
        inline ISATmultithermo<ThermoMixture, Type> operator+
            (
                const ISATmultithermo<ThermoMixture, Type>&,
                const ISATmultithermo<ThermoMixture, Type>&
                );

        template<class ThermoMixture, template<class> class Type>
        inline ISATmultithermo<ThermoMixture, Type> operator*
            (
                const scalar,
                const ISATmultithermo<ThermoMixture, Type>&
                );

        template<class ThermoMixture, template<class> class Type>
        inline ISATmultithermo<ThermoMixture, Type> operator==
            (
                const ISATmultithermo<ThermoMixture, Type>&,
                const ISATmultithermo<ThermoMixture, Type>&
                );


        template<class ThermoMixture, template<class> class Type>
        Ostream& operator<<
            (
                Ostream&,
                const ISATmultithermo<ThermoMixture, Type>&
                );

        template<class Type>
        struct Hnflashfunction
        {
            Type* p;
            Hnflashfunction(Type* pin) :p(pin) {}
            void value(const scalarList& in, scalarList& out, scalar T0) const   //ToDo ADd t0 as sep agru
            {
                scalar H = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                scalarList temp(p->X_.size());
                out.resize(1);
                double sum = 0;
                int n_sp = p->X_.size();
                for (int i = 0;i < n_sp - 1;i++)
                {
                    temp[i] = p->X_[i];
                    p->X_[i] = in[i];
                    sum += in[i];
                }
                temp[n_sp - 1] = p->X_[n_sp - 1];
                p->X_[n_sp - 1] = 1 - sum;
                out[0] = p->T(H, P, T0, &Type::Hs, &Type::Cp, &Type::limit);
                //Info << out << endl;
                for (int i = 0;i < n_sp;i++)
                {
                    p->X_[i] = temp[i];
                }

            }

            void  derive(const scalarList& in, scalarRectangularMatrix& out, scalar T0)
            {
                scalarList temp = in;
                scalarList tempv;
                value(in, tempv, T0);
                scalarList tempv2;
                scalarList dd;
                double sum = 0;
                for (int i = 0;i < in.size() - 2;i++)
                    sum += in[i];
                scalar dx = 1e-6;
                for (int i = 0;i < in.size();i++)
                {
                    temp[i] += dx;
                    value(temp, tempv2, T0);
                    dd = (tempv2 - tempv) / dx;
                    for (int j = 0;j < tempv2.size();j++)
                    {
                        out[i][j] = dd[j];
                    }
                    temp[i] -= dx;
                }
            }
        };

        template<class Type>
        struct thermofunction
        {
            Type* p;
            typedef scalar(Type::* Tfunc) (const scalar p, const scalar T, autoPtr<typename Type::solution>& sol) const;
            Tfunc pfunc;
            thermofunction(Type* pin, Tfunc fin) :p(pin), pfunc(fin) {}
            void value(const scalarList& in, scalarList& out) const
            {
                scalarList temp(p->X_.size());
                out.resize(1);
                double sum = 0;
                int n_sp = p->X_.size();
                for (int i = 0;i < n_sp - 1;i++)
                {
                    temp[i] = p->X_[i];
                    p->X_[i] = in[i];
                    sum += in[i];
                }
                temp[n_sp - 1] = p->X_[n_sp - 1];
                p->X_[n_sp - 1] = 1 - sum;
                autoPtr<typename Type::solution> sol(p->TPn_flash(in[in.size() - 2], in[in.size() - 1]));
                out[0] = (p->*pfunc)(in[in.size() - 2], in[in.size() - 1], sol);
                //Info << out << endl;
                for (int i = 0;i < n_sp;i++)
                {
                    p->X_[i] = temp[i];
                }

            }

            void  derive(const scalarList& in, scalarRectangularMatrix& out)
            {
                scalarList temp = in;
                scalarList tempv;
                value(in, tempv);
                scalarList tempv2;
                scalarList dd;
                double sum = 0;
                for (int i = 0;i < in.size() - 2;i++)
                    sum += in[i];
                scalar dx = 1e-6;
                for (int i = 0;i < in.size();i++)
                {
                    temp[i] += dx;
                    value(temp, tempv2);
                    dd = (tempv2 - tempv) / dx;
                    for (int j = 0;j < tempv2.size();j++)
                    {
                        out[i][j] = dd[j];
                    }
                    temp[i] -= dx;
                }
            }
        };
        template<class Type>
        struct thermofunction_psi
        {
            Type* p;
            typedef scalar(Type::* Tfunc) (const scalar p, const scalar T, autoPtr<typename Type::solution>& sol) const;
            thermofunction_psi(Type* pin) :p(pin) {}
            void value(const scalarList& in, scalarList& out, scalar T0) const
            {
                scalar H = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                scalar T = p->T(H, P, T0, &Type::Hs, &Type::Cp, &Type::limit);
                //scalarList temp(p->X_.size());
                out.resize(1);
                /*
                for (int i = 0;i < p->X_.size();i++)
                {
                    temp[i] = in[i];
                }
                */
                p->setX(in);
                autoPtr<typename Type::solution> sol(p->TPn_flash(P, T));

                out[0] = (p->Type::thermoVLE::psi)(P, T, sol);
            }

            void  derive(const scalarList& in, scalarRectangularMatrix& out, scalar T0)
            {
                //static double maxx = 0;
                scalar H = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                scalar T = p->T(H, P, T0, &Type::Hs, &Type::Cp, &Type::limit);
                /*
                scalarList temp(p->X_.size());
                int n_sp = p->X_.size();
                for (int i = 0;i < n_sp;i++)
                {
                    temp[i] = in[i];
                }
                */
                p->setX(in);
                autoPtr<typename Type::solution> sol(p->TPn_flash(P, T));
                autoPtr<scalarList> grad(p->drhoPdXHP_HsP(P, T, sol()));
                for (int i = 0;i < grad->size();i++)
                {
                    out[i][0] = grad()[i];
                }
                /*
                for (int i = 0;i < p->X_.size();i++)
                {
                    out[i][0] = p->drhoPdXi_HsP(P, T, i, sol());
                }
                out[p->X_.size()][0] = p->drhoPdH_HsP(P, T, sol());
                out[p->X_.size() + 1][0] = p->drhoPdP_HsP(P, T, sol());
                */
                //Info << grad() << endl;
                //Info << out << endl;
                /*
                for (int i = 0;i < 4;i++)
                {
                    if (fabs((out[i][0] - grad()[i]) / grad()[i]) > maxx)
                        maxx = fabs((out[i][0] - grad()[i]) / grad()[i]);
                    //Info << (out[i][0] - grad()[i]) / grad()[i] << ",";
                }
                if (maxx > 1e-06)
                {
                    Info << grad() << endl;
                    Info << out << endl;
                    Info << maxx << endl;
                }*/
                //Info << maxx << endl;
                //FatalErrorInFunction << "vaporf 1" << exit(FatalError);
            }
        };
        template<class Type>
        struct thermofunction_Tpsi
        {
            Type* p;
            thermofunction_Tpsi(Type* pin) :p(pin) {}
            void value(const scalarList& in, scalarList& out, scalar T0) const
            {
                scalar H = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                scalar T = p->T(H, P, T0, &Type::Hs, &Type::Cp, &Type::limit);
                //scalarList temp(p->X_.size());
                out.resize(2);
                /*
                for (int i = 0;i < p->X_.size();i++)
                {
                    temp[i] = in[i];
                }
                */
                p->setX(in);
                autoPtr<typename Type::solution> sol(p->TPn_flash(P, T));
                out[0] = T;
                out[1] = (p->Type::thermoVLE::psi)(P, T, sol);
            }

            void  derive(const scalarList& in, scalarRectangularMatrix& out, scalar T0)
            {
                scalar H = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                scalar T = p->T(H, P, T0, &Type::Hs, &Type::Cp, &Type::limit);

                p->setX(in);
                autoPtr<typename Type::solution> sol(p->TPn_flash(P, T));
                autoPtr<scalarRectangularMatrix> grad(p->dTrhoPdXHP_HsP(P, T, sol()));
                out = grad();
            }
        };

        template<class Type>
        struct thermofunction_Tpsivf
        {
            Type* p;
            thermofunction_Tpsivf(Type* pin) :p(pin) {}
            void value(const scalarList& in, scalarList& out, scalar T0) const
            {
                scalar H = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                scalar T = p->T(H, P, T0, &Type::Hs, &Type::Cp, &Type::limit);
                //scalarList temp(p->X_.size());
                out.resize(3);
                /*
                for (int i = 0;i < p->X_.size();i++)
                {
                    temp[i] = in[i];
                }
                */
                p->setX(in);
                autoPtr<typename Type::solution> sol(p->TPn_flash(P, T));
                out[0] = T;
                out[1] = (p->Type::thermoVLE::psi)(P, T, sol);
                out[2] = sol().vaporfra;
            }

            void  derive(const scalarList& in, scalarRectangularMatrix& out, scalar T0)
            {
                scalar H = in[in.size() - 2];
                scalar P = in[in.size() - 1];
                scalar T = p->T(H, P, T0, &Type::Hs, &Type::Cp, &Type::limit);

                p->setX(in);
                autoPtr<typename Type::solution> sol(p->TPn_flash(P, T));
                autoPtr<scalarRectangularMatrix> grad(p->dTrhoPvfdXHP_HsP(P, T, sol()));
                out = grad();
            }
        };



        template<class Type>
        struct thermofunction2
        {
            Type* p;
            typedef scalar(Type::* Tfunc) (const scalar p, const scalar T) const;
            Tfunc pfunc;
            thermofunction2(Type* pin, Tfunc fin) :p(pin), pfunc(fin) {}
            void value(const scalarList& in, scalarList& out) const
            {
                scalarList temp(p->X_.size());
                out.resize(1);
                double sum = 0;
                int n_sp = p->X_.size();
                for (int i = 0;i < n_sp - 1;i++)
                {
                    temp[i] = p->X_[i];
                    p->X_[i] = in[i];
                    sum += in[i];
                }
                temp[n_sp - 1] = p->X_[n_sp - 1];
                p->X_[n_sp - 1] = 1 - sum;
                out[0] = (p->*pfunc)(in[in.size() - 2], in[in.size() - 1]);
                //Info << out << endl;
                for (int i = 0;i < n_sp;i++)
                {
                    p->X_[i] = temp[i];
                }

            }

            void  derive(const scalarList& in, scalarRectangularMatrix& out)
            {
                scalarList temp = in;
                scalarList tempv;
                value(in, tempv);
                scalarList tempv2;
                scalarList dd;
                double sum = 0;
                for (int i = 0;i < in.size() - 2;i++)
                    sum += in[i];
                scalar dx = 1e-6;
                for (int i = 0;i < in.size();i++)
                {
                    temp[i] += dx;
                    value(temp, tempv2);
                    dd = (tempv2 - tempv) / dx;
                    for (int j = 0;j < tempv2.size();j++)
                    {
                        out[i][j] = dd[j];
                    }
                    temp[i] -= dx;
                }
            }
        };

        template<class Type>
        struct psi_vffunction
        {
            Type* p;
            typedef scalar(Type::* Tfunc) (const scalar p, const scalar T, autoPtr<typename Type::solution>& sol) const;
            Tfunc ppsi, pvf;
            psi_vffunction(Type* pin, Tfunc ppsi_in, Tfunc pvf_in) :p(pin), ppsi(ppsi_in), pvf(pvf_in) {}
            void value(const scalarList& in, scalarList& out) const
            {
                scalarList temp(p->X_.size());
                out.resize(2);
                double sum = 0;
                int n_sp = p->X_.size();
                for (int i = 0;i < n_sp - 1;i++)
                {
                    temp[i] = p->X_[i];
                    p->X_[i] = in[i];
                    sum += in[i];
                }
                temp[n_sp - 1] = p->X_[n_sp - 1];
                p->X_[n_sp - 1] = 1 - sum;
                autoPtr<typename Type::solution> sol(p->TPn_flash(in[in.size() - 2], in[in.size() - 1]));
                out[0] = (p->*ppsi)(in[in.size() - 2], in[in.size() - 1], sol);
                out[1] = (p->*pvf)(in[in.size() - 2], in[in.size() - 1], sol);
                //Info << out << endl;
                for (int i = 0;i < n_sp;i++)
                {
                    p->X_[i] = temp[i];
                }

            }

            void  derive(const scalarList& in, scalarRectangularMatrix& out)
            {
                scalarList temp = in;
                scalarList tempv;
                value(in, tempv);
                scalarList tempv2;
                scalarList dd;
                double sum = 0;
                for (int i = 0;i < in.size() - 2;i++)
                    sum += in[i];
                scalar dx = 1e-6;
                for (int i = 0;i < in.size();i++)
                {
                    temp[i] += dx;
                    value(temp, tempv2);
                    dd = (tempv2 - tempv) / dx;
                    for (int j = 0;j < tempv2.size();j++)
                    {
                        out[i][j] = dd[j];
                    }
                    temp[i] -= dx;
                }
            }
        };


        /*---------------------------------------------------------------------------*\
                                   Class ISATmultithermo Declaration
        \*---------------------------------------------------------------------------*/

        template<class ThermoMixture, template<class> class Type>
        class ISATmultithermo
            :
            public ThermoMixture,
            public Type<ISATmultithermo<ThermoMixture, Type>>
        {
            // Private data
            
                //- Convergence tolerance of energy -> temperature inversion functions
            static const scalar tol_;

            //- Max number of iterations in energy->temperature inversion functions
            static const int maxIter_;
            dictionary ISATDict_;
            //Hnflashfunction<ISATmultithermo<ThermoMixture, Type>> Hnfun;
            //thermofunction<ThermoMixture> psifun;
            //thermofunction_psi<ISATmultithermo<ThermoMixture, Type>> psifun_HP;
            //thermofunction_Tpsi<ISATmultithermo<ThermoMixture, Type>> Tpsifun_HP;
            thermofunction_Tpsivf<ISATmultithermo<ThermoMixture, Type>> Tpsivffun_HP;
            thermofunction2<ThermoMixture> mufun;
            thermofunction2<ThermoMixture> alphahfun;
            //thermofunction<ThermoMixture> vaporfrafun;
            //psi_vffunction<ThermoMixture> psivffun;
            
            
            //mutable ISATmanager<Hnflashfunction<ISATmultithermo<ThermoMixture, Type>>> Ttree;
            //mutable ISATmanager< thermofunction<ThermoMixture>> psitree;
            //mutable ISATmanager< thermofunction_psi<ISATmultithermo<ThermoMixture, Type>>> psitree_HP;
            //mutable ISATmanager< thermofunction_Tpsi<ISATmultithermo<ThermoMixture, Type>>> Tpsitree_HP;
            mutable ISATmanager< thermofunction_Tpsivf<ISATmultithermo<ThermoMixture, Type>>> Tpsivftree_HP;
            mutable ISATmanager< thermofunction2<ThermoMixture>> mutree;
            mutable ISATmanager< thermofunction2<ThermoMixture>> alphahtree;
            //mutable ISATmanager< thermofunction<ThermoMixture>> vaporfratree;
            //mutable ISATmanager< psi_vffunction<ThermoMixture>> psivftree;

            //ISATmanage<> Ttree;
            //ISATmanager<> psitree;
            //ISATmanager<> mutree;
            //ISATmanager<> alphahtree;





            // Private Member Functions

                //- Return the temperature corresponding to the value of the
                //  thermodynamic property f, given the function f = F(p, T)
                //  and dF(p, T)/dT
            inline scalar T
            (
                scalar f,
                scalar p,
                scalar T0,
                scalar(ISATmultithermo::* F)(const scalar, const scalar) const,
                scalar(ISATmultithermo::* dFdT)(const scalar, const scalar) const,
                scalar(ISATmultithermo::* limit)(const scalar) const
            ) const;
            /*
                    inline scalar T_IQI
                    (
                        scalar f,
                        scalar p,
                        scalar T0,
                        scalar (ISATmultithermo::*F)(const scalar, const scalar) const,
                        scalar (ISATmultithermo::*dFdT)(const scalar, const scalar) const,
                        scalar (ISATmultithermo::*limit)(const scalar) const
                    ) const;
            */

        public:

            //- The thermodynamics of the individual species'
            typedef ThermoMixture thermoVLE;
            typedef ISATmultithermo<ThermoMixture, Type> thermoType;
            typedef typename ThermoMixture::SingleThermoType SingleThermoType;

            // Constructors

                //- Construct from components
            inline ISATmultithermo(const ThermoMixture& sp);

            //- Construct from dictionary
            ISATmultithermo(const dictionary& dict, PtrList<ThermoMixture>& speciesData);

            //- Construct as named copy
            inline ISATmultithermo(const word& name, const ISATmultithermo&);


            inline ISATmultithermo(const word& name, PtrList<SingleThermoType>&, const speciesTable&, const dictionary&);


            // Member Functions

                //- Return the instantiated type name
            static word typeName()
            {
                return ThermoMixture::typeName();
            }

            //- Name of Enthalpy/Internal energy
            static inline word heName();

            inline void newTimeStep();
            // Fundamental properties
            // (These functions must be provided in derived types)

                // Heat capacity at constant pressure [J/(kg K)]
                // inline scalar Cp(const scalar p, const scalar T) const;

                // Sensible enthalpy [J/kg]
                // inline scalar Hs(const scalar p, const scalar T) const;

                // Chemical enthalpy [J/kg]
                // inline scalar Hc() const;

                // Absolute Enthalpy [J/kg]
                // inline scalar Ha(const scalar p, const scalar T) const;

                // Entropy [J/(kg K)]
                // inline scalar S(const scalar p, const scalar T) const;
            inline scalar psi(const scalar p, const scalar T) const;
            inline scalar psi_a(const scalar p, const scalar T, const scalar T0) const;
            inline scalar psi_HP(const scalar p, const scalar T, const scalar T0) const;
            inline std::tuple<scalar, scalar> Tpsi_HP(const scalar p, const scalar T, const scalar T0) const;
            inline std::tuple<scalar, scalar, scalar> Tpsivf_HP(const scalar p, const scalar T, const scalar T0) const;

            inline scalar mu(const scalar p, const scalar T) const;
            inline scalar alphah(const scalar p, const scalar T) const;
            inline scalar vaporfra(const scalar p, const scalar T) const;
            inline void psi_vf(const scalar p, const scalar T, scalar& psi_out, scalar& vf_out) const;
            //inline scalar psi(const scalar p, const scalar T) const;

            // Mass specific derived properties

                //- Heat capacity at constant volume [J/(kg K)]
            inline scalar Cv(const scalar p, const scalar T) const;

            //- Heat capacity at constant pressure/volume [J/(kg K)]
            inline scalar Cpv(const scalar p, const scalar T) const;

            //- Gamma = Cp/Cv []
            inline scalar gamma(const scalar p, const scalar T) const;

            //- Ratio of heat capacity at constant pressure to that at
            //  constant pressure/volume []
            inline scalar CpByCpv(const scalar p, const scalar T) const;

            //- Enthalpy/Internal energy [J/kg]
            inline scalar HE(const scalar p, const scalar T) const;

            //- Sensible internal energy [J/kg]
            inline scalar Es(const scalar p, const scalar T) const;

            //- Absolute internal energy [J/kg]
            inline scalar Ea(const scalar p, const scalar T) const;

            //- Gibbs free energy [J/kg]
            inline scalar G(const scalar p, const scalar T) const;

            //- Helmholtz free energy [J/kg]
            inline scalar A(const scalar p, const scalar T) const;


            // Mole specific derived properties

                //- Heat capacity at constant pressure [J/(kmol K)]
            inline scalar cp(const scalar p, const scalar T) const;

            //- Absolute Enthalpy [J/kmol]
            inline scalar ha(const scalar p, const scalar T) const;

            //- Sensible enthalpy [J/kmol]
            inline scalar hs(const scalar p, const scalar T) const;

            //- Chemical enthalpy [J/kmol]
            inline scalar hc() const;

            //- Entropy [J/(kmol K)]
            inline scalar s(const scalar p, const scalar T) const;

            //- Enthalpy/Internal energy [J/kmol]
            inline scalar he(const scalar p, const scalar T) const;

            //- Heat capacity at constant volume [J/(kmol K)]
            inline scalar cv(const scalar p, const scalar T) const;

            //- Sensible internal energy [J/kmol]
            inline scalar es(const scalar p, const scalar T) const;

            //- Absolute internal energy [J/kmol]
            inline scalar ea(const scalar p, const scalar T) const;

            //- Gibbs free energy [J/kmol]
            inline scalar g(const scalar p, const scalar T) const;

            //- Helmholtz free energy [J/kmol]
            inline scalar a(const scalar p, const scalar T) const;


            // Equilibrium reaction thermodynamics

                //- Equilibrium constant [] i.t.o fugacities
                //  = PIi(fi/Pstd)^nui
            inline scalar K(const scalar p, const scalar T) const;

            //- Equilibrium constant [] i.t.o. partial pressures
            //  = PIi(pi/Pstd)^nui
            //  For low pressures (where the gas mixture is near perfect) Kp = K
            inline scalar Kp(const scalar p, const scalar T) const;

            //- Equilibrium constant i.t.o. molar concentration
            //  = PIi(ci/cstd)^nui
            //  For low pressures (where the gas mixture is near perfect)
            //  Kc = Kp(pstd/(RR*T))^nu
            inline scalar Kc(const scalar p, const scalar T) const;

            //- Equilibrium constant [] i.t.o. mole-fractions
            //  For low pressures (where the gas mixture is near perfect)
            //  Kx = Kp(pstd/p)^nui
            inline scalar Kx
            (
                const scalar p,
                const scalar T
            ) const;

            //- Equilibrium constant [] i.t.o. number of moles
            //  For low pressures (where the gas mixture is near perfect)
            //  Kn = Kp(n*pstd/p)^nui where n = number of moles in mixture
            inline scalar Kn
            (
                const scalar p,
                const scalar T,
                const scalar n
            ) const;


            // Energy->temperature  inversion functions

                //- Temperature from enthalpy or internal energy
                //  given an initial temperature T0
            inline scalar THE
            (
                const scalar H,
                const scalar p,
                const scalar T0
            ) const;

            //- Temperature from sensible enthalpy given an initial T0
            inline scalar THs
            (
                const scalar Hs,
                const scalar p,
                const scalar T0
            ) const;

            //- Temperature from absolute enthalpy
            //  given an initial temperature T0
            inline scalar THa
            (
                const scalar H,
                const scalar p,
                const scalar T0
            ) const;

            //- Temperature from sensible internal energy
            //  given an initial temperature T0
            inline scalar TEs
            (
                const scalar E,
                const scalar p,
                const scalar T0
            ) const;

            //- Temperature from absolute internal energy
            //  given an initial temperature T0
            inline scalar TEa
            (
                const scalar E,
                const scalar p,
                const scalar T0
            ) const;


            // Derivative term used for Jacobian

                //- Derivative of B (acooding to Niemeyer et al.) w.r.t. temperature
            inline scalar dKcdTbyKc(const scalar p, const scalar T) const;

            //- Derivative of cp w.r.t. temperature
            inline scalar dcpdT(const scalar p, const scalar T) const;


            // I-O

                //- Write to Ostream
            void write(Ostream& os) const;


            // Member operators

            inline void operator+=(const ISATmultithermo&);
            inline void operator*=(const scalar);


            // Friend operators

            friend ISATmultithermo operator+ <ThermoMixture, Type>
                (
                    const ISATmultithermo&,
                    const ISATmultithermo&
                    );

            friend ISATmultithermo operator* <ThermoMixture, Type>
                (
                    const scalar s,
                    const ISATmultithermo&
                    );

            friend ISATmultithermo operator== <ThermoMixture, Type>
                (
                    const ISATmultithermo&,
                    const ISATmultithermo&
                    );



            // Ostream Operator

            friend Ostream& operator<< <ThermoMixture, Type>
                (
                    Ostream&,
                    const ISATmultithermo&
                    );

            template<class TypeName>
            friend struct Hnflashfunction;

            template<class TypeName>
            friend struct thermofunction_psi;
            template<class TypeName>
            friend struct thermofunction_Tpsi;
            template<class TypeName>
            friend struct thermofunction_Tpsivf;
        };




        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace species
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ISATmultithermoI.H"

#ifdef NoRepository
#include "ISATmultithermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
