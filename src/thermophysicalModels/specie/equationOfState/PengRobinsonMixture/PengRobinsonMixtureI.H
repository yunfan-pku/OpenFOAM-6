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

\*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* PengRobinson EOS                                                          */
/* Thermal solver                                                            */
/* Developed by   Ping Yi,  03/01/2019,   UMN                                */
/* mix means each phase mixture; ave means two-phase mixture                 */
/*---------------------------------------------------------------------------*/

#include "PengRobinsonMixture.H"
#include "PengRobinson.H"
#include "specie.H"
#include "mathematicalConstants.H"
#include "OFstream.H"
#include "IFstream.H"
#include "IOmanip.H"
#include "string.H"
#include "tmp.H"

using namespace Foam;
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class ThermoMixture>
inline Foam::PengRobinsonMixture<ThermoMixture>::PengRobinsonMixture(
    const word& name,
    PtrList<SingleThermoType>& speciesData,
    const speciesTable& specieNames,
    const dictionary& thermoDict)
    : ThermoMixture(name, speciesData, specieNames, thermoDict),
    bico(this->N_),
    phaseFlag(1)
{
    forAll(this->X_, i)
    {
        bico[i].resize(this->N_);
    }

    speciesTable ts(thermoDict.subDict("BinaryInteractionParameter").lookup("specielist"));
    scalarList bilist(thermoDict.subDict("BinaryInteractionParameter").lookup("matrix"));
    int indexi, indexj, len = ts.size();
    forAll(this->X_, i)
    {
        indexi = ts[this->species()[i]];
        forAll(this->X_, j)
        {
            indexj = ts[this->species()[j]];
            bico[i][j] = bilist[indexi * len + indexj];
            //Info<<bico[i][j]<<" ";
        }
        //Info<<endl;
    }
}

template <class ThermoMixture>
Foam::PengRobinsonMixture<ThermoMixture>::PengRobinsonMixture(
    const dictionary& dict, PtrList<SingleThermoType>& speciesData)
    : ThermoMixture(dict, speciesData),
    bico(this->N_),
    phaseFlag(1)
{
    forAll(this->X_, i)
    {
        bico[i].resize(this->N_);
    }

    speciesTable ts(dict.subDict("BinaryInteractionParameter").lookup("specielist"));
    scalarList bilist(dict.subDict("BinaryInteractionParameter").lookup("matrix"));
    int indexi, indexj, len = ts.size();
    forAll(this->X_, i)
    {
        indexi = ts[this->species()[i]];
        forAll(this->X_, j)
        {
            indexj = ts[this->species()[j]];
            bico[i][j] = bilist[indexi * len + indexj];
            //Info<<bico[i][j]<<" ";
        }
        //Info<<endl;
    }
}

template <class ThermoMixture>
inline Foam::autoPtr<Foam::PengRobinsonMixture<ThermoMixture>>
Foam::PengRobinsonMixture<ThermoMixture>::clone() const
{
    return autoPtr<PengRobinsonMixture<ThermoMixture>>(
        new PengRobinsonMixture<ThermoMixture>(*this));
}

template <class ThermoMixture>
inline Foam::autoPtr<Foam::PengRobinsonMixture<ThermoMixture>>
Foam::PengRobinsonMixture<ThermoMixture>::New(
    const dictionary& dict)
{
    return autoPtr<PengRobinsonMixture<ThermoMixture>>(
        new PengRobinsonMixture<ThermoMixture>(dict));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::W // kg/mol
(const scalarList* Xinp) const
{
    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
        Xp = Xinp;
    const scalarList& X = *Xp;

    const PengRobinsonMixture<ThermoMixture>& This = *this;
    scalar Wmix = 0.0;
    forAll(X, i)
    {
        Wmix += X[i] * This[i].W(); //PY: Y() is mole number of component in mixture
    }
    return Wmix * 1.0e-03;
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::rho // PY:phase mixture density
(
    scalar p,
    scalar T,
    const label pf,
    const scalarList* Xinp) const
{                                                                //RR: 10e+03*m3 Pa K-1 mol-1
    return W(Xinp) * p / (Z(p, T, pf, Xinp) * RR * 1.0e-03 * T); // kg/m3
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::drhodT // PY:phase mixture density
(
    scalar p,
    scalar T,
    const label pf,
    const scalarList* Xinp) const
{
    scalar z = Z(p, T, pf, Xinp);
    //RR: 10e+03*m3 Pa K-1 mol-1
    return -W(Xinp) * p / (RR * 1.0e-03) * (dZdT(p, T, pf, Xinp) * T + z) / (z * z * T * T); // kg/m3
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::drhodP // PY:phase mixture density
(
    scalar p,
    scalar T,
    const label pf,
    const scalarList* Xinp) const
{
    scalar z = Z(p, T, pf, Xinp);
    //RR: 10e+03*m3 Pa K-1 mol-1
    return W(Xinp) / (RR * 1.0e-03 * T) * (z - p * dZdP(p, T, pf, Xinp)) / (z * z); // kg/m3
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::V //m3/mol
(
    scalar p,
    scalar T,
    const label pf,
    const scalarList* Xinp) const
{
    return W(Xinp) / rho(p, T, pf, Xinp);
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::A(
    scalar p,
    scalar T,
    const scalarList* Xinp) const
{
    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
        Xp = Xinp;
    const scalarList& X = *Xp;

    const PengRobinsonMixture<ThermoMixture>& This = *this;
    scalar Amix = 0.0;
    forAll(X, i)
    {
        forAll(X, j)
        {
            Amix += X[i] * X[j] * sqrt(This[i].EOS::A(p, T) * This[j].EOS::A(p, T)) * (1.0 - bico[i][j]);
        }
    }
    return Amix;
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::dAdT(
    scalar p,
    scalar T,
    const scalarList* Xinp) const
{
    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
        Xp = Xinp;
    const scalarList& X = *Xp;

    const PengRobinsonMixture<ThermoMixture>& This = *this;
    scalar Amix = 0.0;
    forAll(X, i)
    {
        forAll(X, j)
        {
            Amix += X[i] * X[j] * (sqrt(This[i].EOS::A(p, T) / This[j].EOS::A(p, T)) * This[j].EOS::dAdT(p, T) + sqrt(This[j].EOS::A(p, T) / This[i].EOS::A(p, T)) * This[i].EOS::dAdT(p, T)) * 0.5 * (1.0 - bico[i][j]);
        }
    }
    return Amix;
}
template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::dAdP(
    scalar p,
    scalar T,
    const scalarList* Xinp) const
{
    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
        Xp = Xinp;
    const scalarList& X = *Xp;

    const PengRobinsonMixture<ThermoMixture>& This = *this;
    scalar Amix = 0.0;
    forAll(X, i)
    {
        forAll(X, j)
        {
            Amix += X[i] * X[j] * (sqrt(This[i].EOS::A(p, T) / This[j].EOS::A(p, T)) * This[j].EOS::dAdP(p, T) + sqrt(This[j].EOS::A(p, T) / This[i].EOS::A(p, T)) * This[i].EOS::dAdP(p, T)) * 0.5 * (1.0 - bico[i][j]);
        }
    }
    return Amix;
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::dAdxi(
    scalar p,
    scalar T,
    label di,
    const scalarList* Xinp) const
{
    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
        Xp = Xinp;
    const scalarList& X = *Xp;

    const PengRobinsonMixture<ThermoMixture>& This = *this;
    scalar dAmixdxi = 0.0;
    forAll(X, i)
    {
        dAmixdxi += 2.0 * X[i] * sqrt(This[i].EOS::A(p, T) * This[di].EOS::A(p, T)) * (1.0 - bico[i][di]);
    }
    return dAmixdxi;
}




template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::B(
    scalar p,
    scalar T,
    const scalarList* Xinp) const
{
    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
        Xp = Xinp;
    const scalarList& X = *Xp;

    const PengRobinsonMixture<ThermoMixture>& This = *this;
    scalar Bmix = 0.0;
    forAll(X, i)
    {
        Bmix += X[i] * This[i].B(p, T); //0.31845932267269694,0.22670292099211661,0.32089979844514827
    }
    return Bmix;
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::dBdT(
    scalar p,
    scalar T,
    const scalarList* Xinp) const
{
    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
        Xp = Xinp;
    const scalarList& X = *Xp;

    const PengRobinsonMixture<ThermoMixture>& This = *this;
    scalar Bmix = 0.0;
    forAll(X, i)
    {
        Bmix += X[i] * This[i].dBdT(p, T); //0.31845932267269694,0.22670292099211661,0.32089979844514827
    }
    return Bmix;
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::dBdP(
    scalar p,
    scalar T,
    const scalarList* Xinp) const
{
    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
        Xp = Xinp;
    const scalarList& X = *Xp;

    const PengRobinsonMixture<ThermoMixture>& This = *this;
    scalar Bmix = 0.0;
    forAll(X, i)
    {
        Bmix += X[i] * This[i].dBdP(p, T); //0.31845932267269694,0.22670292099211661,0.32089979844514827
    }
    return Bmix;
}
template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::dBdxi(
    scalar p,
    scalar T,
    label di,
    const scalarList* Xinp) const
{
    //const scalarList* Xp = &this->X_;
    //if (Xinp != NULL)
    //    Xp = Xinp;
    //const scalarList& X = *Xp;

    const PengRobinsonMixture<ThermoMixture>& This = *this;
    scalar Bmix = 0.0;
    Bmix += This[di].B(p, T); //0.31845932267269694,0.22670292099211661,0.32089979844514827
    return Bmix;
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::a(
    scalar p,
    scalar T,
    const scalarList* Xinp) const
{
    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
        Xp = Xinp;
    const scalarList& X = *Xp;

    const PengRobinsonMixture<ThermoMixture>& This = *this;
    scalar amix = 0.0;
    forAll(X, i)
    {
        forAll(X, j)
        {
            amix += X[i] * X[j] * sqrt(This[i].EOS::a(p, T) * This[j].EOS::a(p, T)) * (1.0 - bico[i][j]);
        }
    }
    return amix;
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::dadt(
    scalar p,
    scalar T,
    const scalarList* Xinp) const
{
    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
        Xp = Xinp;
    const scalarList& X = *Xp;

    const PengRobinsonMixture<ThermoMixture>& This = *this;
    scalar dadtmix = 0.0, ai, aj;
    forAll(X, i)
    {
        forAll(X, j)
        {
            ai = This[i].EOS::a(p, T);
            aj = This[j].EOS::a(p, T);

            dadtmix += X[i] * X[j] * 0.5 * (sqrt(ai / aj) * This[j].dadt(p, T) + sqrt(aj / ai) * This[i].dadt(p, T)) * (1.0 - bico[i][j]);
        }
    }
    return dadtmix;
}


template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::dadp(
    scalar p,
    scalar T,
    const scalarList* Xinp) const
{
    return 0.0;
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::dadxi(
    scalar p,
    scalar T,
    label di,
    const scalarList* Xinp) const
{
    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
        Xp = Xinp;
    const scalarList& X = *Xp;

    const PengRobinsonMixture<ThermoMixture>& This = *this;
    scalar dadtmix = 0.0;//, ai, adi;
    //adi = This[di].EOS::a(p, T);
    forAll(X, i)
    {
        //ai = This[i].EOS::a(p, T);
        dadtmix += 2 * X[i] * sqrt(This[i].EOS::a(p, T) * This[di].EOS::a(p, T)) * (1.0 - bico[i][di]);
    }
    return dadtmix;
}
template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::d2adtdxi(
    scalar p,
    scalar T,
    label di,
    const scalarList* Xinp) const
{
    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
        Xp = Xinp;
    const scalarList& X = *Xp;

    const PengRobinsonMixture<ThermoMixture>& This = *this;
    scalar dadtmix = 0.0, ai, adi;
    adi = This[di].EOS::a(p, T);
    forAll(X, i)
    {
        ai = This[i].EOS::a(p, T);
        dadtmix += X[i] * (sqrt(ai / adi) * This[di].dadt(p, T) + sqrt(adi / ai) * This[i].dadt(p, T)) * (1.0 - bico[i][di]);
    }
    return dadtmix;
}




template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::d2adt2(
    scalar p,
    scalar T,
    const scalarList* Xinp) const
{
    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
        Xp = Xinp;
    const scalarList& X = *Xp;

    const PengRobinsonMixture<ThermoMixture>& This = *this;
    scalar d2adt2mix = 0.0, ai, aj, dai, daj;
    forAll(X, i)
    {
        forAll(X, j)
        {
            ai = This[i].EOS::a(p, T);
            aj = This[j].EOS::a(p, T);
            dai = This[i].dadt(p, T);
            daj = This[j].dadt(p, T);
            d2adt2mix += X[i] * X[j] * 0.25 / (ai * aj * sqrt(ai * aj)) * (-sqr(aj * dai - ai * daj) + 2 * ai * aj * (ai * This[j].d2adt2(p, T) + aj * This[i].d2adt2(p, T))) * (1.0 - bico[i][j]);
        }
    }
    return d2adt2mix;
}
template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::b(
    scalar p,
    scalar T,
    const scalarList* Xinp) const
{
    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
        Xp = Xinp;
    const scalarList& X = *Xp;

    const PengRobinsonMixture<ThermoMixture>& This = *this;
    scalar bmix = 0.0;
    forAll(X, i)
    {
        bmix += X[i] * This[i].b(p, T);
    }
    return bmix;
}



template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::cubicSolve(
    scalar a2,
    scalar a1,
    scalar a0,
    label flag)
{

    scalar Q = (3.0 * a1 - a2 * a2) / 9.0;
    scalar Rl = (9.0 * a2 * a1 - 27.0 * a0 - 2.0 * a2 * a2 * a2) / 54.0;
    scalar Q3 = Q * Q * Q;
    scalar D = Q3 + Rl * Rl;
    scalar root = -1;
    if (D <= 0)
    {
        scalar th = ::acos(Rl / sqrt(-Q3));
        scalar qm = 2 * sqrt(-Q);
        scalar r1 = qm * cos(th / 3) - a2 / 3;
        scalar r2 = qm * cos((th + 2 * constant::mathematical::pi) / 3) - a2 / 3;
        scalar r3 = qm * cos((th + 4 * constant::mathematical::pi) / 3) - a2 / 3;

        // three root;
        if (flag == 1) //gas
        {
            root = max(r1, max(r2, r3));
        }
        else //liquid
        {
            root = max(r1, max(r2, r3));
            root = min(r1 > 1e-10 ? r1 : root, root);
            root = min(r2 > 1e-10 ? r2 : root, root);
            root = min(r3 > 1e-10 ? r3 : root, root);
        }
    }
    else
    {
        // One root is real: in supercritical zone
        scalar D05 = sqrt(D);
        scalar S = sign(Rl + D05) * pow(mag(Rl + D05), 1.0 / 3.0);
        scalar Tl = sign(Rl - D05) * pow(mag(Rl - D05), 1.0 / 3.0);

        root = S + Tl - a2 / 3;
    }

    if (root <= 0.0)
        FatalErrorInFunction << "Z is negative! Z = " << root << exit(FatalError);
    return root;
}


template <class ThermoMixture>
inline std::tuple<int, double, double, double> Foam::PengRobinsonMixture<ThermoMixture>::cubicSolve(
    scalar a2,
    scalar a1,
    scalar a0)
{

    scalar Q = (3.0 * a1 - a2 * a2) / 9.0;
    scalar Rl = (9.0 * a2 * a1 - 27.0 * a0 - 2.0 * a2 * a2 * a2) / 54.0;
    scalar Q3 = Q * Q * Q;
    scalar D = Q3 + Rl * Rl;
    scalar root = -1;
    if (D <= 0)
    {
        scalar th = ::acos(Rl / sqrt(-Q3));
        scalar qm = 2 * sqrt(-Q);
        scalar r1 = qm * cos(th / 3) - a2 / 3;
        scalar r2 =
            qm * cos((th + 2 * constant::mathematical::pi) / 3) - a2 / 3;
        scalar r3 =
            qm * cos((th + 4 * constant::mathematical::pi) / 3) - a2 / 3;


        return std::make_tuple(3, r1, r2, r3);
    }
    else
    {
        // One root is real: in supercritical zone
        scalar D05 = sqrt(D);
        scalar S = sign(Rl + D05) * pow(mag(Rl + D05), 1.0 / 3.0);
        scalar Tl = sign(Rl - D05) * pow(mag(Rl - D05), 1.0 / 3.0);
        root = S + Tl - a2 / 3;
        return std::make_tuple(1, root, root, root);
    }
    return std::make_tuple(1, root, root, root);
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::Z(
    scalar p,
    scalar T,
    const label pf,
    const scalarList* Xinp) const
{
    scalar A = this->A(p, T, Xinp), B = this->B(p, T, Xinp);
    scalar a2 = B - 1.0;
    scalar a1 = A - 2.0 * B - 3.0 * B * B;
    scalar a0 = -A * B + B * B + B * B * B;
    scalar r1, r2, r3;
    int nr;
    std::tie(nr, r1, r2, r3) = cubicSolve(a2, a1, a0);
    if (nr == 1)
        return r1;
    scalar zmax = max(r1, max(r2, r3));
    scalar zmin = min(r1, min(r2, r3));
    if (zmin < B)
    {
        if (zmax < B)
            FatalErrorInFunction << "zmax<B,zmax= " << zmax << exit(FatalError);
        return zmax;
    }
    if (pf == 1)
        return zmax;
    if (pf == 0)
        return zmin;
    FatalErrorInFunction << "pf  = " << pf << exit(FatalError);
    return zmin;
    /*
    if (pf != 2)
        return cubicSolve(a2, a1, a0, pf);
    else
        return cubicSolve(a2, a1, a0, phaseFlag);
        */
}


template <class ThermoMixture>
inline std::tuple<int, double, double, double> Foam::PengRobinsonMixture<ThermoMixture>::Z_full(
    scalar p,
    scalar T,
    const scalarList* Xinp) const
{
    scalar A = this->A(p, T, Xinp), B = this->B(p, T, Xinp);
    scalar a2 = B - 1.0;
    scalar a1 = A - 2.0 * B - 3.0 * B * B;
    scalar a0 = -A * B + B * B + B * B * B;
    return cubicSolve(a2, a1, a0);
}


template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::Z_gibbs(
    scalar p,
    scalar T,
    const scalarList* Xinp) const
{
    int nr;
    scalar  r1, r2, r3;
    std::tie(nr, r1, r2, r3) = Z_full(p, T, Xinp);
    scalar z_ret;
    scalar sqrt2 = sqrt(2.0);
    if (nr == 1)
    {
        z_ret = r1;
    }
    else
    {
        scalar z_v = max(r1, max(r2, r3));
        scalar z_l = min(r1, min(r2, r3));
        scalar A = this->A(p, T, Xinp), B = this->B(p, T, Xinp);
        if (z_l < B)
        {
            z_ret = z_v;
            if (z_v < B)
                FatalErrorInFunction << "z_v < B,z_v= " << z_v << exit(FatalError);
        }
        else
        {
            scalar dg = log((z_l - B) / (z_v - B)) + A / (B * 2 * sqrt2) * log((z_l + B * (1 + sqrt2)) * (z_v + B * (1 - sqrt2)) / (z_l + B * (1 - sqrt2)) / (z_v + B * (1 + sqrt2))) - (z_l - z_v);

            if (dg < 0.0)
                z_ret = z_v;
            else
                z_ret = z_l;
        }
    }
    return z_ret;
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::vol_gibbs(
    scalar p,
    scalar T,
    const scalarList* Xinp) const
{

    return Z_gibbs(p, T, Xinp) * RR * 1.0e-03 * T / p;
}



template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::dZdT(
    scalar p,
    scalar T,
    const label pf,
    const scalarList* Xinp) const
{
    scalar A = this->A(p, T, Xinp), B = this->B(p, T, Xinp);
    scalar dAdT = this->dAdT(p, T, Xinp), dBdT = this->dBdT(p, T, Xinp);
    scalar Z0 = this->Z(p, T, pf, Xinp);
    scalar a2 = B - 1.0;
    scalar da2dT = dBdT;
    scalar a1 = A - 2.0 * B - 3.0 * B * B;
    scalar da1dT = dAdT - 2.0 * dBdT - 6.0 * B * dBdT;
    //scalar a0 = -A * B + B * B + B * B * B;
    scalar da0dT = -(dAdT * B + A * dBdT) + 2 * B * dBdT + 3 * B * B * dBdT;

    return -(da2dT * Z0 * Z0 + da1dT * Z0 + da0dT) / (3 * Z0 * Z0 + 2 * Z0 * a2 + a1);
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::dZdP(
    scalar p,
    scalar T,
    const label pf,
    const scalarList* Xinp) const
{
    scalar A = this->A(p, T, Xinp), B = this->B(p, T, Xinp);
    scalar dAdP = this->dAdP(p, T, Xinp), dBdP = this->dBdP(p, T, Xinp);
    scalar Z0 = this->Z(p, T, pf, Xinp);
    scalar a2 = B - 1.0;
    scalar da2dP = dBdP;
    scalar a1 = A - 2.0 * B - 3.0 * B * B;
    scalar da1dP = dAdP - 2.0 * dBdP - 6.0 * B * dBdP;
    //scalar a0 = -A * B + B * B + B * B * B;
    scalar da0dP = -(dAdP * B + A * dBdP) + 2 * B * dBdP + 3 * B * B * dBdP;

    return -(da2dP * Z0 * Z0 + da1dP * Z0 + da0dP) / (3 * Z0 * Z0 + 2 * Z0 * a2 + a1);
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::dZdxi(
    scalar p,
    scalar T,
    label di,
    const label pf,
    const scalarList* Xinp) const
{
    scalar A = this->A(p, T, Xinp), B = this->B(p, T, Xinp);
    scalar dAdxi = this->dAdxi(p, T, di, Xinp), dBdxi = this->dBdxi(p, T, di, Xinp);
    scalar Z0 = this->Z(p, T, pf, Xinp);
    scalar a2 = B - 1.0;
    scalar da2dxi = dBdxi;
    scalar a1 = A - 2.0 * B - 3.0 * B * B;
    scalar da1dxi = dAdxi - 2.0 * dBdxi - 6.0 * B * dBdxi;
    //scalar a0 = -A * B + B * B + B * B * B;
    scalar da0dxi = -(dAdxi * B + A * dBdxi) + 2 * B * dBdxi + 3 * B * B * dBdxi;

    return -(da2dxi * Z0 * Z0 + da1dxi * Z0 + da0dxi) / (3 * Z0 * Z0 + 2 * Z0 * a2 + a1);
}

template <class ThermoMixture>
inline Foam::autoPtr<scalarList> Foam::PengRobinsonMixture<ThermoMixture>::A2(
    scalar p, //inp
    scalar T, //inp
    const scalarList* Xinp) const
{
    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
    {
        Xp = Xinp;
    }
    const scalarList& X = *Xp;

    scalarList& Amix2 = *(new scalarList(this->N_));

    forAll(X, i)
    {
        Amix2[i] = 0.0;
        forAll(X, j)
        {
            Amix2[i] += X[j] * sqrt((*this)[i].EOS::A(p, T) * (*this)[j].EOS::A(p, T)) * (1.0 - bico[j][i]);
        }
    }
    return autoPtr<scalarList>(&Amix2);
}//01,0.13700000000000001,02,0.10299999999999999,12,0.17000000000000001,21,0.031199999999999999

template <class ThermoMixture>
inline Foam::autoPtr<scalarList> Foam::PengRobinsonMixture<ThermoMixture>::dA2dT(
    scalar p, //inp
    scalar T, //inp
    const scalarList* Xinp) const
{
    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
    {
        Xp = Xinp;
    }
    const scalarList& X = *Xp;

    scalarList& dAmix2dT = *(new scalarList(this->N_));

    forAll(X, i)
    {
        dAmix2dT[i] = 0.0;
        forAll(X, j)
        {
            dAmix2dT[i] += X[j] * 0.5 * (sqrt((*this)[i].EOS::A(p, T) / (*this)[j].EOS::A(p, T)) * (*this)[j].EOS::dAdT(p, T) + sqrt((*this)[j].EOS::A(p, T) / (*this)[i].EOS::A(p, T)) * (*this)[i].EOS::dAdT(p, T)) * (1.0 - bico[j][i]);
        }
    }
    return autoPtr<scalarList>(&dAmix2dT);
}
template <class ThermoMixture>
inline Foam::autoPtr<scalarList> Foam::PengRobinsonMixture<ThermoMixture>::dA2dP(
    scalar p, //inp
    scalar T, //inp
    const scalarList* Xinp) const
{
    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
    {
        Xp = Xinp;
    }
    const scalarList& X = *Xp;

    scalarList& dAmix2dT = *(new scalarList(this->N_));

    forAll(X, i)
    {
        dAmix2dT[i] = 0.0;
        forAll(X, j)
        {
            dAmix2dT[i] += X[j] * 0.5 * (sqrt((*this)[i].EOS::A(p, T) / (*this)[j].EOS::A(p, T)) * (*this)[j].EOS::dAdP(p, T) + sqrt((*this)[j].EOS::A(p, T) / (*this)[i].EOS::A(p, T)) * (*this)[i].EOS::dAdP(p, T)) * (1.0 - bico[j][i]);
        }
    }
    return autoPtr<scalarList>(&dAmix2dT);
}

template <class ThermoMixture>
inline Foam::autoPtr<scalarList> Foam::PengRobinsonMixture<ThermoMixture>::dA2dxi(
    scalar p, //inp
    scalar T, //inp
    label di,
    const scalarList* Xinp) const
{
    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
    {
        Xp = Xinp;
    }
    const scalarList& X = *Xp;

    scalarList& Amix2 = *(new scalarList(this->N_));

    forAll(X, i)
    {
        //Amix2[i] = 0.0;

        Amix2[i] = sqrt((*this)[i].EOS::A(p, T) * (*this)[di].EOS::A(p, T)) * (1.0 - bico[di][i]);

    }
    return autoPtr<scalarList>(&Amix2);
}


// PY: fugacity coefficient of species in each phase
template <class ThermoMixture>
inline Foam::autoPtr<scalarList> Foam::PengRobinsonMixture<ThermoMixture>::fugacityCoefficient(
    scalar p, //input
    scalar T,
    const label pf,
    const scalarList* Xinp) const
{
    label pflag = this->phaseFlag;
    if (pf != 2)
        pflag = pf;

    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
    {
        Xp = Xinp;
    }
    const scalarList& X = *Xp;

    scalarList& fugcoef = *(new scalarList(this->N_));
    scalar Amixz = this->A(p, T, Xinp);
    scalar Bmixz = this->B(p, T, Xinp);
    scalar Zmix = 0.0;
    autoPtr<scalarList> pAmix2(A2(p, T, Xinp));
    const scalarList& Amix2 = pAmix2();


    const scalar small8 = 1.e-08; //PY: need to be noticed

    //this->ABZmix(p, T, comp_inp, Amixz, Bmixz, Amix2, Bz);

    Zmix = this->Z(p, T, pf, Xinp);


    scalar Zmix0 = this->Z(p, T, 1, Xinp);

    scalar comptrifile = 0.7;
    scalarList Xs(this->X_.size());

    if (pflag == 0)
    {
        while (Zmix < 0.0)
        {

            //FatalErrorInFunction << "Zmix < 0.0" << exit(FatalError);
            Xs[0] = comptrifile;
            for (int i = 1; i < this->X_.size(); i++)
                Xs[i] = (1.0 - comptrifile) / (Xs.size() - 1);
            Zmix = this->Z(p, T, pflag, &Xs);
            comptrifile = comptrifile + 0.01;
            if (comptrifile >= 1.0)
            {
                Zmix = Zmix0;
            }
        }
    }

    scalar cc0 = 2.0 * sqrt(2.0);
    scalar cc1 = 1.0 + sqrt(2.0);
    scalar cc2 = 1.0 - sqrt(2.0);

    scalar zfactor00 = Zmix;
    scalar ZMB = zfactor00 - Bmixz;
    if (ZMB < 0.0)
    {
        //FatalErrorInFunction << "ZMB < 0.0" << exit(FatalError);
        zfactor00 = Bmixz + small8;
        ZMB = small8;
    }
    forAll(X, spid)
    {
        scalar b = (*this)[spid].B(p, T);
        fugcoef[spid] = exp(b / Bmixz * (zfactor00 - 1.0) - log(ZMB) - Amixz / (cc0 * Bmixz) * (2.0 * Amix2[spid] / Amixz - b / Bmixz) * log((zfactor00 + cc1 * Bmixz) / (zfactor00 + cc2 * Bmixz)));
    }
    return autoPtr<scalarList>(&fugcoef);
}

template <class ThermoMixture>
inline Foam::autoPtr<scalarList> Foam::PengRobinsonMixture<ThermoMixture>::ln_fugacityCoefficient(
    scalar p, //input
    scalar T,
    scalar z,
    const scalarList* Xinp) const
{

    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
    {
        Xp = Xinp;
    }
    const scalarList& X = *Xp;

    scalarList& fugcoef = *(new scalarList(this->N_));
    scalar Amixz = this->A(p, T, Xinp);
    scalar Bmixz = this->B(p, T, Xinp);
    scalar Zmix = z;
    autoPtr<scalarList> pAmix2(A2(p, T, Xinp));
    const scalarList& Amix2 = pAmix2();




    scalarList Xs(this->X_.size());


    scalar cc0 = 2.0 * sqrt(2.0);
    scalar cc1 = 1.0 + sqrt(2.0);
    scalar cc2 = 1.0 - sqrt(2.0);

    scalar zfactor00 = Zmix;
    scalar ZMB = zfactor00 - Bmixz;
    if (ZMB < 0.0)
    {
        FatalErrorInFunction << "Z-B < 0.0" << exit(FatalError);
        const scalar small8 = 1.e-08; //PY: need to be noticed
        zfactor00 = Bmixz + small8;
        ZMB = small8;
    }
    forAll(X, spid)
    {
        scalar b = (*this)[spid].B(p, T);
        fugcoef[spid] = b / Bmixz * (zfactor00 - 1.0) - log(ZMB) - Amixz / (cc0 * Bmixz) * (2.0 * Amix2[spid] / Amixz - b / Bmixz) * log((zfactor00 + cc1 * Bmixz) / (zfactor00 + cc2 * Bmixz));
    }
    return autoPtr<scalarList>(&fugcoef);
}

template <class ThermoMixture>
inline Foam::autoPtr<scalarList> Foam::PengRobinsonMixture<ThermoMixture>::ddT_Ln_fugacityCoefficient(
    scalar p, //input
    scalar T,
    const label pf,
    const scalarList* Xinp) const
{
    label pflag = this->phaseFlag;
    if (pf != 2)
        pflag = pf;

    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
    {
        Xp = Xinp;
    }
    const scalarList& X = *Xp;

    scalarList& ddT_Ln_fugcoef = *(new scalarList(this->N_));
    scalar Amixz = this->A(p, T, Xinp);
    scalar dAmixzdT = this->dAdT(p, T, Xinp);
    scalar Bmixz = this->B(p, T, Xinp);
    scalar dBmixzdT = this->dBdT(p, T, Xinp);
    scalar Zmix = 0.0;
    autoPtr<scalarList> pAmix2(A2(p, T, Xinp));
    const scalarList& Amix2 = pAmix2();

    autoPtr<scalarList> pdAmix2dT(dA2dT(p, T, Xinp));
    const scalarList& dAmix2dT = pdAmix2dT();
    const scalar small8 = 1.e-08; //PY: need to be noticed


    Zmix = this->Z(p, T, pf, Xinp);
    scalar dZmixdT = this->dZdT(p, T, pf, Xinp);

    scalar Zmix0 = this->Z(p, T, 1, Xinp);

    scalar comptrifile = 0.7;
    scalarList Xs(this->X_.size());

    if (pflag == 0)
    {
        while (Zmix < 0.0)
        {   //TODO fix Zmix<0
            FatalErrorInFunction << "Zmix < 0.0" << exit(FatalError);
            Xs[0] = comptrifile;
            for (int i = 1; i < this->X_.size(); i++)
                Xs[i] = (1.0 - comptrifile) / (Xs.size() - 1);
            Zmix = this->Z(p, T, pflag, &Xs);
            comptrifile = comptrifile + 0.01;
            if (comptrifile >= 1.0)
            {
                Zmix = Zmix0;
            }
        }
    }

    scalar cc0 = 2.0 * sqrt(2.0);
    scalar cc1 = 1.0 + sqrt(2.0);
    scalar cc2 = 1.0 - sqrt(2.0);

    scalar zfactor00 = Zmix;
    scalar ZMB = zfactor00 - Bmixz;

    if (ZMB <= 0.0)
    {
        FatalErrorInFunction << "ZMB <= 0.0" << exit(FatalError);
        zfactor00 = Bmixz + small8;
        ZMB = small8;
    }
    forAll(X, spid)
    {
        scalar Bi = (*this)[spid].B(p, T);
        scalar dBidT = (*this)[spid].dBdT(p, T);
        //ddT_Ln_fugcoef[spid] = Bi / Bmixz * (Zmix - 1);
        ddT_Ln_fugcoef[spid] = dBidT / Bmixz * (Zmix - 1) - Bi / (Bmixz * Bmixz) * dBmixzdT * (Zmix - 1) + Bi / Bmixz * dZmixdT;
        ddT_Ln_fugcoef[spid] -= (dZmixdT - dBmixzdT) / (Zmix - Bmixz);
        ddT_Ln_fugcoef[spid] -= (dAmix2dT[spid] / (sqrt(2.0) * Bmixz) - Amix2[spid] * dBmixzdT / (sqrt(2.0) * Bmixz * Bmixz) - (dAmixzdT * Bi + Amixz * dBidT) / (cc0 * Bmixz * Bmixz) + (Amixz * Bi * dBmixzdT) / (sqrt(2.0) * Bmixz * Bmixz * Bmixz)) * log((zfactor00 + cc1 * Bmixz) / (zfactor00 + cc2 * Bmixz));

        ddT_Ln_fugcoef[spid] -= (Amix2[spid] / (sqrt(2.0) * Bmixz) - Amixz * Bi / (cc0 * Bmixz * Bmixz)) * ((dZmixdT + cc1 * dBmixzdT) / (Zmix + cc1 * Bmixz) - (dZmixdT + cc2 * dBmixzdT) / (Zmix + cc2 * Bmixz));
    }
    return autoPtr<scalarList>(&ddT_Ln_fugcoef);
}

template <class ThermoMixture>
inline Foam::autoPtr<scalarList> Foam::PengRobinsonMixture<ThermoMixture>::ddP_Ln_fugacityCoefficient(
    scalar p, //input
    scalar T,
    const label pf,
    const scalarList* Xinp) const
{
    label pflag = this->phaseFlag;
    if (pf != 2)
        pflag = pf;

    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
    {
        Xp = Xinp;
    }
    const scalarList& X = *Xp;

    scalarList& ddT_Ln_fugcoef = *(new scalarList(this->N_));
    scalar Amixz = this->A(p, T, Xinp);
    scalar dAmixzdT = this->dAdP(p, T, Xinp);
    scalar Bmixz = this->B(p, T, Xinp);
    scalar dBmixzdT = this->dBdP(p, T, Xinp);
    scalar Zmix = 0.0;
    autoPtr<scalarList> pAmix2(A2(p, T, Xinp));
    const scalarList& Amix2 = pAmix2();

    autoPtr<scalarList> pdAmix2dT(dA2dP(p, T, Xinp));
    const scalarList& dAmix2dT = pdAmix2dT();
    const scalar small8 = 1.e-08; //PY: need to be noticed


    Zmix = this->Z(p, T, pf, Xinp);
    scalar dZmixdT = this->dZdP(p, T, pf, Xinp);

    scalar Zmix0 = this->Z(p, T, 1, Xinp);

    scalar comptrifile = 0.7;
    scalarList Xs(this->X_.size());

    if (pflag == 0)
    {
        while (Zmix < 0.0)
        {
            FatalErrorInFunction << "Zmix < 0.0" << exit(FatalError);
            Xs[0] = comptrifile;
            for (int i = 1; i < this->X_.size(); i++)
                Xs[i] = (1.0 - comptrifile) / (Xs.size() - 1);
            Zmix = this->Z(p, T, pflag, &Xs);
            comptrifile = comptrifile + 0.01;
            if (comptrifile >= 1.0)
            {
                Zmix = Zmix0;
            }
        }
    }

    scalar cc0 = 2.0 * sqrt(2.0);
    scalar cc1 = 1.0 + sqrt(2.0);
    scalar cc2 = 1.0 - sqrt(2.0);

    scalar zfactor00 = Zmix;
    scalar ZMB = zfactor00 - Bmixz;

    if (ZMB <= 0.0)
    {
        FatalErrorInFunction << "ZMB <= 0.0" << exit(FatalError);
        zfactor00 = Bmixz + small8;
        ZMB = small8;
    }
    forAll(X, spid)
    {
        scalar Bi = (*this)[spid].B(p, T);
        scalar dBidT = (*this)[spid].dBdP(p, T);
        //ddT_Ln_fugcoef[spid] = Bi / Bmixz * (Zmix - 1);
        ddT_Ln_fugcoef[spid] = dBidT / Bmixz * (Zmix - 1) - Bi / (Bmixz * Bmixz) * dBmixzdT * (Zmix - 1) + Bi / Bmixz * dZmixdT;
        ddT_Ln_fugcoef[spid] -= (dZmixdT - dBmixzdT) / (Zmix - Bmixz);
        ddT_Ln_fugcoef[spid] -= (dAmix2dT[spid] / (sqrt(2.0) * Bmixz) - Amix2[spid] * dBmixzdT / (sqrt(2.0) * Bmixz * Bmixz) - (dAmixzdT * Bi + Amixz * dBidT) / (cc0 * Bmixz * Bmixz) + (Amixz * Bi * dBmixzdT) / (sqrt(2.0) * Bmixz * Bmixz * Bmixz)) * log((zfactor00 + cc1 * Bmixz) / (zfactor00 + cc2 * Bmixz));

        ddT_Ln_fugcoef[spid] -= (Amix2[spid] / (sqrt(2.0) * Bmixz) - Amixz * Bi / (cc0 * Bmixz * Bmixz)) * ((dZmixdT + cc1 * dBmixzdT) / (Zmix + cc1 * Bmixz) - (dZmixdT + cc2 * dBmixzdT) / (Zmix + cc2 * Bmixz));
    }
    return autoPtr<scalarList>(&ddT_Ln_fugcoef);
}

template <class ThermoMixture>
inline Foam::autoPtr<scalarList> Foam::PengRobinsonMixture<ThermoMixture>::ddxi_Ln_fugacityCoefficient(
    scalar p, //input
    scalar T,
    label di,
    const label pf,
    const scalarList* Xinp) const
{
    label pflag = this->phaseFlag;
    if (pf != 2)
        pflag = pf;

    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
    {
        Xp = Xinp;
    }
    const scalarList& X = *Xp;

    scalarList& ddxi_Ln_fugcoef = *(new scalarList(this->N_));
    scalar Amixz = this->A(p, T, Xinp);
    scalar dAmixzdxi = this->dAdxi(p, T, di, Xinp);
    scalar Bmixz = this->B(p, T, Xinp);
    scalar dBmixzdxi = this->dBdxi(p, T, di, Xinp);
    scalar Zmix = 0.0;
    autoPtr<scalarList> pAmix2(A2(p, T, Xinp));
    const scalarList& Amix2 = pAmix2();

    autoPtr<scalarList> pdAmix2dxi(dA2dxi(p, T, di, Xinp));
    const scalarList& dAmix2dxi = pdAmix2dxi();
    const scalar small8 = 1.e-08; //PY: need to be noticed


    Zmix = this->Z(p, T, pf, Xinp);
    scalar dZmixdxi = this->dZdxi(p, T, di, pf, Xinp);


    scalar Zmix0 = this->Z(p, T, 1, Xinp);

    scalar comptrifile = 0.7;
    scalarList Xs(this->X_.size());

    if (pflag == 0)
    {
        while (Zmix < 0.0)
        {
            FatalErrorInFunction << "Zmix < 0.0" << exit(FatalError);
            Xs[0] = comptrifile;
            for (int i = 1; i < this->X_.size(); i++)
                Xs[i] = (1.0 - comptrifile) / (Xs.size() - 1);
            Zmix = this->Z(p, T, pflag, &Xs);
            comptrifile = comptrifile + 0.01;
            if (comptrifile >= 1.0)
            {
                Zmix = Zmix0;
            }
        }
    }

    scalar cc0 = 2.0 * sqrt(2.0);
    scalar cc1 = 1.0 + sqrt(2.0);
    scalar cc2 = 1.0 - sqrt(2.0);

    scalar zfactor00 = Zmix;
    scalar ZMB = zfactor00 - Bmixz;

    if (ZMB <= 0.0)
    {
        FatalErrorInFunction << "ZMB <= 0.0" << exit(FatalError);
        zfactor00 = Bmixz + small8;
        ZMB = small8;
    }
    forAll(X, spid)
    {
        scalar Bi = (*this)[spid].B(p, T);
        scalar dBidxi = 0;
        //ddT_Ln_fugcoef[spid] = Bi / Bmixz * (Zmix - 1);
        ddxi_Ln_fugcoef[spid] = dBidxi / Bmixz * (Zmix - 1) - Bi / (Bmixz * Bmixz) * dBmixzdxi * (Zmix - 1) + Bi / Bmixz * dZmixdxi;
        ddxi_Ln_fugcoef[spid] -= (dZmixdxi - dBmixzdxi) / (Zmix - Bmixz);
        ddxi_Ln_fugcoef[spid] -= (dAmix2dxi[spid] / (sqrt(2.0) * Bmixz) - Amix2[spid] * dBmixzdxi / (sqrt(2.0) * Bmixz * Bmixz) - (dAmixzdxi * Bi + Amixz * dBidxi) / (cc0 * Bmixz * Bmixz) + (Amixz * Bi * dBmixzdxi) / (sqrt(2.0) * Bmixz * Bmixz * Bmixz)) * log((zfactor00 + cc1 * Bmixz) / (zfactor00 + cc2 * Bmixz));

        ddxi_Ln_fugcoef[spid] -= (Amix2[spid] / (sqrt(2.0) * Bmixz) - Amixz * Bi / (cc0 * Bmixz * Bmixz)) * ((dZmixdxi + cc1 * dBmixzdxi) / (Zmix + cc1 * Bmixz) - (dZmixdxi + cc2 * dBmixzdxi) / (Zmix + cc2 * Bmixz));
    }
    return autoPtr<scalarList>(&ddxi_Ln_fugcoef);
}



template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::H(scalar p, scalar T) const
{
    return 0;
}
/*
template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::Cp(scalar p, scalar T) const
{
    return 0;
}
*/
template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::S(
    scalar p,
    scalar T) const
{
    return 0;
}
template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::psi(
    scalar p,
    scalar T) const
{
    const scalar Z = this->Z(p, T);

    return 1.0 / (Z * this->R() * T);
}
template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::Hs // internal energy with unit of J/kg
(
    scalar p,
    scalar T,
    const label pf,
    const scalarList* Xinp) const
{

    return Ha(p, T, pf, Xinp) - this->Hideal(1, Tstd, Xinp); //+11990129.829480369 +3991.2*1000;
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::Ha // internal energy with unit of J/kg
(
    scalar p,
    scalar T,
    const label pf,
    const scalarList* Xinp) const
{
    scalar bmix = this->b(p, T, Xinp);
    scalar cc1 = 1.0 + sqrt(2.0);
    scalar cc2 = 1.0 - sqrt(2.0);

    scalar VV = this->W(Xinp) / this->rho(p, T, pf, Xinp); //m3*mol-1, phase specific volume

    scalar xxx = ((VV + cc1 * bmix) / (VV + cc2 * bmix));
    scalar dH;
    //scalar a_b=this->a(p, T, Xinp);
    //scalar dadt_b=this->dadt(p, T, Xinp);
    //scalar z_b=Z(p, T, pf, Xinp);
    //Info<<"!!!dH1="<<Z(p, T, pf, Xinp) <<"!!!\n";
    if (xxx > 0)
    {
        dH = (Z(p, T, pf, Xinp) - 1) * RR * 1e-3 * T + ((T * this->dadt(p, T, Xinp) - this->a(p, T, Xinp)) / (2.0 * sqrt(2.0) * bmix)) * log(xxx);
        //dH = (T * this->dadt(p, T, Xinp) - this->a(p, T, Xinp))/bmix;
        //dH = 1;
        //dH = 0 * (Z(p, T, pf, Xinp) - 1) * RR * 1e-3 * T + this->dadt(p, T, Xinp);
    }
    else
    {
        dH = (Z(p, T, pf, Xinp) - 1) * RR * 1e-3 * T;
    }

    dH = dH / this->W(Xinp); // J/kg
    // Based on Eq.12, Internal energy is estimated based on ideal(close to absolute 0 T) + departure
    //Info<<"!!!dH="<<dH<<"!!!Hi="<<this->Hideal(p, T, Xinp)<<"!!!\n";
    return dH + this->Hideal(p, T, Xinp);
    //+11990129.829480369 +3991.2*1000;
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::dHadT // internal energy with unit of J/kg
(
    scalar p,
    scalar T,
    const label pf,
    const scalarList* Xinp) const
{
    scalar bmix = this->b(p, T, Xinp);
    //scalar dbmixdT = 0;
    scalar cc1 = 1.0 + sqrt(2.0);
    scalar cc2 = 1.0 - sqrt(2.0);

    scalar rho_ = this->rho(p, T, pf, Xinp);
    scalar VV = this->W(Xinp) / rho_; //m3*mol-1, phase specific volume
    scalar dVVdT = -this->W(Xinp) / (rho_ * rho_) * this->drhodT(p, T, pf, Xinp);
    scalar xxx = ((VV + cc1 * bmix) / (VV + cc2 * bmix));
    scalar ddT_dH;
    scalar z = Z(p, T, pf, Xinp);
    //scalar a_b=this->a(p, T, Xinp);
    //scalar dadt_b=this->dadt(p, T, Xinp);
    //scalar z_b=Z(p, T, pf, Xinp);
    //Info<<"!!!dH1="<<Z(p, T, pf, Xinp) <<"!!!\n";
    if (xxx > 0)
    {
        ddT_dH = dZdT(p, T, pf, Xinp) * RR * 1e-3 * T + (z - 1) * RR * 1e-3;
        ddT_dH += ((T * this->d2adt2(p, T, Xinp)) / (2.0 * sqrt(2.0) * bmix)) * log(xxx);
        //ddT_dH = this->d2adt2(p, T, Xinp);
        ddT_dH += ((T * this->dadt(p, T, Xinp) - this->a(p, T, Xinp)) / (2.0 * sqrt(2.0) * bmix)) * (dVVdT / (VV + cc1 * bmix) - dVVdT / (VV + cc2 * bmix));
    }
    else
    {
        ddT_dH = dZdT(p, T, pf, Xinp) * RR * 1e-3 * T + (z - 1) * RR * 1e-3;
    }
    //Info << "haha=" << setw(20) << setprecision(20) << this->dHadT(p, T, &flag, &(sol().X_gas())) << endl;
    ddT_dH = ddT_dH / this->W(Xinp); // J/kg
    // Based on Eq.12, Internal energy is estimated based on ideal(close to absolute 0 T) + departure
    //Info<<"!!!dH="<<dH<<"!!!Hi="<<this->Hideal(p, T, Xinp)<<"!!!\n";
    return ddT_dH + this->dHidealdT(p, T, Xinp);
}
template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::dHsdT // internal energy with unit of J/kg
(
    scalar p,
    scalar T,
    const label pf,
    const scalarList* Xinp) const
{
    return dHadT(p, T, pf, Xinp); //+11990129.829480369 +3991.2*1000;
}
template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::dHadP // internal energy with unit of J/kg
(
    scalar p,
    scalar T,
    const label pf,
    const scalarList* Xinp) const
{
    scalar bmix = this->b(p, T, Xinp);
    //scalar dbmixdT = 0;
    scalar cc1 = 1.0 + sqrt(2.0);
    scalar cc2 = 1.0 - sqrt(2.0);

    scalar rho_ = this->rho(p, T, pf, Xinp);
    scalar VV = this->W(Xinp) / rho_; //m3*mol-1, phase specific volume
    scalar dVVdT = -this->W(Xinp) / (rho_ * rho_) * this->drhodP(p, T, pf, Xinp);
    scalar xxx = ((VV + cc1 * bmix) / (VV + cc2 * bmix));
    scalar ddT_dH;
    //scalar z = Z(p, T, pf, Xinp);
    //scalar a_b=this->a(p, T, Xinp);
    //scalar dadt_b=this->dadt(p, T, Xinp);
    //scalar z_b=Z(p, T, pf, Xinp);
    //Info<<"!!!dH1="<<Z(p, T, pf, Xinp) <<"!!!\n";
    if (xxx > 0)
    {
        ddT_dH = dZdP(p, T, pf, Xinp) * RR * 1e-3 * T;
        ddT_dH += ((T * this->dadt(p, T, Xinp) - this->a(p, T, Xinp)) / (2.0 * sqrt(2.0) * bmix)) * (dVVdT / (VV + cc1 * bmix) - dVVdT / (VV + cc2 * bmix));
    }
    else
    {
        ddT_dH = dZdP(p, T, pf, Xinp) * RR * 1e-3 * T;
    }
    //Info << "haha=" << setw(20) << setprecision(20) << this->dHadT(p, T, &flag, &(sol().X_gas())) << endl;
    ddT_dH = ddT_dH / this->W(Xinp); // J/kg
    // Based on Eq.12, Internal energy is estimated based on ideal(close to absolute 0 T) + departure
    //Info<<"!!!dH="<<dH<<"!!!Hi="<<this->Hideal(p, T, Xinp)<<"!!!\n";
    return ddT_dH;
}
template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::dHsdP // internal energy with unit of J/kg
(
    scalar p,
    scalar T,
    const label pf,
    const scalarList* Xinp) const
{
    return dHadP(p, T, pf, Xinp);  //+11990129.829480369 +3991.2*1000;
}
template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::dHadxi // internal energy with unit of J/kg
(
    scalar p,
    scalar T,
    label di,
    const label pf,
    const scalarList* Xinp) const
{
    scalar bmix = this->b(p, T, Xinp);
    const PengRobinsonMixture<ThermoMixture>& This = *this;
    scalar dbmixdxi = This[di].b(p, T);

    scalar cc1 = 1.0 + sqrt(2.0);
    scalar cc2 = 1.0 - sqrt(2.0);

    scalar VV = this->W(Xinp) / this->rho(p, T, pf, Xinp); //m3*mol-1, phase specific volume
    scalar dVVdxi = dZdxi(p, T, di, pf, Xinp) * RR * 1.0e-03 * T / p;
    scalar xxx = ((VV + cc1 * bmix) / (VV + cc2 * bmix));
    scalar ddxi_dH, dH;
    //scalar a_b=this->a(p, T, Xinp);
    //scalar dadt_b=this->dadt(p, T, Xinp);
    //scalar z_b=Z(p, T, pf, Xinp);
    //Info<<"!!!dH1="<<Z(p, T, pf, Xinp) <<"!!!\n";
    if (xxx > 0)
    {
        dH = (Z(p, T, pf, Xinp) - 1) * RR * 1e-3 * T + ((T * this->dadt(p, T, Xinp) - this->a(p, T, Xinp)) / (2.0 * sqrt(2.0) * bmix)) * log(xxx);
        //dH = 1;
        ddxi_dH = dZdxi(p, T, di, pf, Xinp) * RR * 1e-3 * T;
        //ddxi_dH = 0;
        ddxi_dH += ((T * this->d2adtdxi(p, T, di, Xinp) - this->dadxi(p, T, di, Xinp)) / (2.0 * sqrt(2.0) * bmix)) * log(xxx);
        ddxi_dH -= ((T * this->dadt(p, T, Xinp) - this->a(p, T, Xinp)) / (2.0 * sqrt(2.0) * bmix * bmix) * dbmixdxi) * log(xxx);
        ddxi_dH += ((T * this->dadt(p, T, Xinp) - this->a(p, T, Xinp)) / (2.0 * sqrt(2.0) * bmix)) * ((dVVdxi + cc1 * dbmixdxi) / (VV + cc1 * bmix) - (dVVdxi + cc2 * dbmixdxi) / (VV + cc2 * bmix));
        //dH = (T * this->dadt(p, T, Xinp) - this->a(p, T, Xinp))/bmix;
        //ddxi_dH = (T * this->d2adtdxi(p, T, di, Xinp) - this->dadxi(p, T, di, Xinp))/bmix+(T * this->dadt(p, T, Xinp) - this->a(p, T, Xinp))/bmix/bmix*dbmixdxi;
    }
    else
    {
        ddxi_dH = dZdxi(p, T, di, pf, Xinp) * RR * 1e-3 * T;
        dH = (Z(p, T, pf, Xinp) - 1) * RR * 1e-3 * T;
    }
    scalar w = this->W(Xinp);
    ddxi_dH = ddxi_dH / w - dH / (w * w) * This[di].W() * 1e-3;
    //dH = dH / this->W(Xinp); // J/kg
    // Based on Eq.12, Internal energy is estimated based on ideal(close to absolute 0 T) + departure
    //Info<<"!!!dH="<<dH<<"!!!Hi="<<this->Hideal(p, T, Xinp)<<"!!!\n";
    return  ddxi_dH + this->dHidealdxi(p, T, di, Xinp);
}
template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::dHsdxi // internal energy with unit of J/kg
(
    scalar p,
    scalar T,
    label di,
    const label pf,
    const scalarList* Xinp) const
{
    return  dHadxi(p, T, di, pf, Xinp) - this->dHidealdxi(1, Tstd, di, Xinp); //+11990129.829480369 +3991.2*1000;
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::Hc // internal energy with unit of J/kg
(
    const scalarList* Xinp) const
{
    return this->Hideal(1, Tstd, Xinp); //+11990129.829480369 +3991.2*1000;
}

//temp E

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::Hideal(
    scalar p,
    scalar T,
    const scalarList* Xinp //note, it is phase compositions
) const
{
    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
    {
        Xp = Xinp;
    }
    const scalarList& X = *Xp;

    scalar H_ideal = 0, sum = 0;
    forAll(X, i)
    {
        sum += X[i] * (*this)[i].W();
        H_ideal += X[i] * (*this)[i].Ha(p, T) * (*this)[i].W();
    }
    return H_ideal / sum;
}


template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::Hideal_Mole(
    scalar p,
    scalar T,
    const scalarList* Xinp //note, it is phase compositions
) const
{
    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
    {
        Xp = Xinp;
    }
    const scalarList& X = *Xp;

    scalar H_ideal = 0, sum = 0;
    forAll(X, i)
    {
        H_ideal += X[i] * (*this)[i].Ha(p, T) * (*this)[i].W();
    }
    return H_ideal * 1e-3;
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::Gideal(
    scalar p,
    scalar T,
    const scalarList* Xinp //note, it is phase compositions
) const
{
    /*
    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
    {
        Xp = Xinp;
    }
    const scalarList& X = *Xp;
    scalar g1 = ((*this)[0].Ha(p, T) - (*this)[0].S(p, T) * T);
    scalar g2 = ((*this)[1].Ha(p, T) - (*this)[1].S(p, T) * T);
    Info << setprecision(20) << "H=" << Hideal(p, T, Xinp) - T * Sideal(p, T, Xinp) << ",Hmy=" << (X[0] * g1 * (*this)[0].W() + X[1] * g2 * (*this)[1].W()) / this->W(Xinp) * 1.0e-03 << ",g1=" << g1 << ",g2=" << g2 << endl;
    */
    return Hideal(p, T, Xinp) - T * Sideal(p, T, Xinp);
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::Gideal_Mole(
    scalar p,
    scalar T,
    const scalarList* Xinp //note, it is phase compositions
) const
{

    return Hideal_Mole(p, T, Xinp) - T * Sideal_Mole(p, T, Xinp);
}


template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::muideal_Mole(     //!! for test [J/mol]
    scalar p,
    scalar T,
    label di, //note, it is phase compositions
    const scalarList* Xinp
) const
{
    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
    {
        Xp = Xinp;
    }
    const scalarList& X = *Xp;
    return 1.0e-03 * (*this)[di].W() * ((*this)[di].Ha(p, T) - T * (*this)[di].S(p, T)) + T * log(X[di]) * RR * 1e-3;
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::Aideal(
    scalar p,
    scalar T,
    const scalarList* Xinp //note, it is phase compositions
) const
{
    return Hideal(p, T, Xinp) - T * Sideal(p, T, Xinp) - RR * 1e-3 * T / this->W(Xinp);
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::Aideal_Mole(
    scalar p,
    scalar T,
    const scalarList* Xinp //note, it is phase compositions
) const
{
    return Hideal_Mole(p, T, Xinp) - T * Sideal_Mole(p, T, Xinp) - RR * 1e-3 * T;
}
/*
template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::G // Gibbs energy with unit of J/kg
(
    scalar p,
    scalar T,
    const label pf,
    const scalarList* Xinp) const
{
    scalar bmix = this->b(p, T, Xinp);
    scalar cc1 = 1.0 + sqrt(2.0);
    scalar cc2 = 1.0 - sqrt(2.0);
    scalar Zmix = this->Z(p, T, pf, Xinp);
    scalar Bmix = this->B(p, T, Xinp);
    scalar VV = this->W(Xinp) / this->rho(p, T, pf, Xinp); //m3*mol-1, phase specific volume

    scalar xxx = ((VV + cc1 * bmix) / (VV + cc2 * bmix));
    scalar dG;
    //scalar a_b=this->a(p, T, Xinp);
    //scalar dadt_b=this->dadt(p, T, Xinp);
    //scalar z_b=Z(p, T, pf, Xinp);
    //Info<<"!!!dH1="<<Z(p, T, pf, Xinp) <<"!!!\n";
    if (xxx > 0)
    {
        dG = (Z(p, T, pf, Xinp) - 1) * RR * 1e-3 * T - log(Zmix - Bmix) * RR * 1e-3 * T + ((-this->a(p, T, Xinp)) / (2.0 * sqrt(2.0) * bmix)) * log(xxx);
        //dH = (T * this->dadt(p, T, Xinp) - this->a(p, T, Xinp))/bmix;
        //dH = 1;
        //dH = 0 * (Z(p, T, pf, Xinp) - 1) * RR * 1e-3 * T + this->dadt(p, T, Xinp);
    }
    else
    {
        dG = (Z(p, T, pf, Xinp) - 1) * RR * 1e-3 * T - log(Zmix - Bmix) * RR * 1e-3 * T;
    }

    dG = dG / this->W(Xinp); // J/kg
    // Based on Eq.12, Internal energy is estimated based on ideal(close to absolute 0 T) + departure
    //Info<<"!!!dH="<<dH<<"!!!Hi="<<this->Hideal(p, T, Xinp)<<"!!!\n";
    return dG + this->Gideal(p, T, Xinp);
    //+11990129.829480369 +3991.2*1000;
}
*/


template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::G_departure_Mole // Gibbs energy with unit of J/kg
(
    scalar p,
    scalar T,
    const scalarList* Xinp) const
{
    scalar bmix = this->b(p, T, Xinp);
    scalar cc1 = 1.0 + sqrt(2.0);
    scalar cc2 = 1.0 - sqrt(2.0);
    scalar Zmix = this->Z_gibbs(p, T, Xinp);
    scalar Bmix = this->B(p, T, Xinp);
    scalar VV = (Zmix * RR * 1.0e-03 * T) / p; //m3*mol-1, phase specific volume

    scalar xxx = ((VV + cc1 * bmix) / (VV + cc2 * bmix));
    scalar dG;
    //scalar a_b=this->a(p, T, Xinp);
    //scalar dadt_b=this->dadt(p, T, Xinp);
    //scalar z_b=Z(p, T, pf, Xinp);
    //Info<<"!!!dH1="<<Z(p, T, pf, Xinp) <<"!!!\n";
    if (xxx > 0)
    {
        dG = (Zmix - 1) * RR * 1e-3 * T - log(Zmix - Bmix) * RR * 1e-3 * T + ((-this->a(p, T, Xinp)) / (2.0 * sqrt(2.0) * bmix)) * log(xxx);
        //dH = (T * this->dadt(p, T, Xinp) - this->a(p, T, Xinp))/bmix;
        //dH = 1;
        //dH = 0 * (Z(p, T, pf, Xinp) - 1) * RR * 1e-3 * T + this->dadt(p, T, Xinp);
    }
    else
    {
        dG = (Zmix - 1) * RR * 1e-3 * T - log(Zmix - Bmix) * RR * 1e-3 * T;
    }
    // Based on Eq.12, Internal energy is estimated based on ideal(close to absolute 0 T) + departure
    //Info<<"!!!dH="<<dH<<"!!!Hi="<<this->Hideal(p, T, Xinp)<<"!!!\n";
    return dG;
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::G_Mole // Gibbs energy with unit of J/kg
(
    scalar p,
    scalar T,
    const scalarList* Xinp) const
{
    return G_departure_Mole(p, T, Xinp) + this->Gideal_Mole(p, T, Xinp);
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::G // Gibbs energy with unit of J/kg
(
    scalar p,
    scalar T,
    const label pf,
    const scalarList* Xinp) const
{
    scalar bmix = this->b(p, T, Xinp);
    scalar cc1 = 1.0 + sqrt(2.0);
    scalar cc2 = 1.0 - sqrt(2.0);
    scalar Zmix = this->Z(p, T, pf, Xinp);
    scalar Bmix = this->B(p, T, Xinp);
    scalar VV = this->W(Xinp) / this->rho(p, T, pf, Xinp); //m3*mol-1, phase specific volume

    scalar xxx = ((VV + cc1 * bmix) / (VV + cc2 * bmix));
    scalar dG;
    //scalar a_b=this->a(p, T, Xinp);
    //scalar dadt_b=this->dadt(p, T, Xinp);
    //scalar z_b=Z(p, T, pf, Xinp);
    //Info<<"!!!dH1="<<Z(p, T, pf, Xinp) <<"!!!\n";
    if (xxx > 0)
    {
        dG = (Z(p, T, pf, Xinp) - 1) * RR * 1e-3 * T - log(Zmix - Bmix) * RR * 1e-3 * T + ((-this->a(p, T, Xinp)) / (2.0 * sqrt(2.0) * bmix)) * log(xxx);
        //dH = (T * this->dadt(p, T, Xinp) - this->a(p, T, Xinp))/bmix;
        //dH = 1;
        //dH = 0 * (Z(p, T, pf, Xinp) - 1) * RR * 1e-3 * T + this->dadt(p, T, Xinp);
    }
    else
    {
        dG = (Z(p, T, pf, Xinp) - 1) * RR * 1e-3 * T - log(Zmix - Bmix) * RR * 1e-3 * T;
    }

    dG = dG / this->W(Xinp); // J/kg
    // Based on Eq.12, Internal energy is estimated based on ideal(close to absolute 0 T) + departure
    //Info<<"!!!dH="<<dH<<"!!!Hi="<<this->Hideal(p, T, Xinp)<<"!!!\n";
    return dG + this->Gideal(p, T, Xinp);
    //+11990129.829480369 +3991.2*1000;
}


template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::G_TPD // Gibbs energy with unit of J/kg
(
    scalar p,
    scalar T,
    const scalarList* Xinp) const
{
    scalar bmix = this->b(p, T, Xinp);
    scalar cc1 = 1.0 + sqrt(2.0);
    scalar cc2 = 1.0 - sqrt(2.0);
    scalar Zmix = this->Z_gibbs(p, T, Xinp);
    scalar Bmix = this->B(p, T, Xinp);
    scalar VV = (Zmix * RR * 1.0e-03 * T) / p; //m3*mol-1, phase specific volume

    scalar xxx = ((VV + cc1 * bmix) / (VV + cc2 * bmix));
    scalar dG;
    //Info << (*Xinp)[0] << ",Zmix =" << Zmix << endl;
    if (xxx > 0)
    {
        dG = (Zmix - 1) * RR * 1e-3 * T - log(Zmix - Bmix) * RR * 1e-3 * T - this->a(p, T, Xinp) / (2.0 * sqrt(2.0) * bmix) * log(xxx);
        //dG = log(Zmix - Bmix) * RR * 1e-3 * T - this->a(p, T, Xinp) / (2.0 * sqrt(2.0) * bmix) * log(xxx);
    }
    else
    {
        dG = (Zmix - 1) * RR * 1e-3 * T - log(Zmix - Bmix) * RR * 1e-3 * T;
        //dG = log(Zmix - Bmix) * RR * 1e-3 * T;
    }

    dG = dG / this->W(Xinp); // J/kg
    return dG + this->Gideal(p, T, Xinp);
    //return  this->Gideal(p, T, Xinp);
}


template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::A_TPD // Gibbs energy with unit of J/kg
(
    scalar p,
    scalar T,
    const scalarList* Xinp) const
{
    scalar bmix = this->b(p, T, Xinp);
    scalar cc1 = 1.0 + sqrt(2.0);
    scalar cc2 = 1.0 - sqrt(2.0);
    scalar Zmix = this->Z_gibbs(p, T, Xinp);
    scalar Bmix = this->B(p, T, Xinp);
    scalar VV = (Zmix * RR * 1.0e-03 * T) / p; //m3*mol-1, phase specific volume

    scalar xxx = ((VV + cc1 * bmix) / (VV + cc2 * bmix));
    scalar dG;
    Info << (*Xinp)[0] << ",Zmix =" << Zmix << endl;
    if (xxx > 0)
    {
        dG = -log(Zmix - Bmix) * RR * 1e-3 * T - this->a(p, T, Xinp) / (2.0 * sqrt(2.0) * bmix) * log(xxx);
        //dG = log(Zmix - Bmix) * RR * 1e-3 * T - this->a(p, T, Xinp) / (2.0 * sqrt(2.0) * bmix) * log(xxx);
    }
    else
    {
        dG = -log(Zmix - Bmix) * RR * 1e-3 * T;
        //dG = log(Zmix - Bmix) * RR * 1e-3 * T;
    }

    dG = dG / this->W(Xinp); // J/kg
    //return dG + this->Gideal(p, T, Xinp);
    return  dG + this->Aideal(p, T, Xinp);
}



template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::Sideal(
    scalar p,
    scalar T,
    const scalarList* Xinp //note, it is phase compositions
) const
{
    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
    {
        Xp = Xinp;
    }
    const scalarList& X = *Xp;

    scalar S_ideal = 0, sum = 0;
    forAll(X, i)
    {
        sum += X[i] * (*this)[i].W();
        S_ideal += X[i] * (*this)[i].S(p, T) * (*this)[i].W();
        if (X[i] > 0.0)
            S_ideal -= X[i] * log(X[i]) * RR;
    }
    //Info << "S_ideal=" << S_ideal / sum << endl;
    return S_ideal / sum;
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::Sideal_Mole(
    scalar p,
    scalar T,
    const scalarList* Xinp //note, it is phase compositions
) const
{
    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
    {
        Xp = Xinp;
    }
    const scalarList& X = *Xp;

    scalar S_ideal = 0;
    forAll(X, i)
    {
        S_ideal += X[i] * (*this)[i].S(p, T) * (*this)[i].W();
        if (X[i] > 0.0)
            S_ideal -= X[i] * log(X[i]) * RR;
    }
    //Info << "S_ideal=" << S_ideal / sum << endl;
    return S_ideal * 1e-3;
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::dHidealdT(
    scalar p,
    scalar T,
    const scalarList* Xinp //note, it is phase compositions
) const
{
    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
    {
        Xp = Xinp;
    }
    const scalarList& X = *Xp;

    scalar dH_idealdT = 0, sum = 0;
    forAll(X, i)
    {
        sum += X[i] * (*this)[i].W();
        dH_idealdT += X[i] * (*this)[i].Cp(p, T) * (*this)[i].W();
    }
    return dH_idealdT / sum;
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::dHidealdxi(
    scalar p,
    scalar T,
    label di,
    const scalarList* Xinp //note, it is phase compositions
) const
{
    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
    {
        Xp = Xinp;
    }
    const scalarList& X = *Xp;
    scalar dH_idealdxi = 0;
    scalar H_ideal = 0, sum = 0;
    forAll(X, i)
    {
        sum += X[i] * (*this)[i].W();
        H_ideal += X[i] * (*this)[i].Ha(p, T) * (*this)[i].W();
    }

    dH_idealdxi = (*this)[di].Ha(p, T) * (*this)[di].W() / sum - H_ideal / (sum * sum) * (*this)[di].W();

    return dH_idealdxi;
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::dHidealdP(
    scalar p,
    scalar T,
    const scalarList* Xinp //note, it is phase compositions
) const
{
    return 0.0;
}

template<class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::CpMCv
(
    scalar p,
    scalar T,
    const label pf,
    const scalarList* Xinp) const
{
    scalar amix = this->a(p, T, Xinp);
    scalar bmix = this->b(p, T, Xinp);
    scalar dadt = this->dadt(p, T, Xinp);

    scalar mw_phase = this->W(Xinp);
    scalar rho_phase = this->rho(p, T, pf, Xinp);
    scalar VV = this->W(Xinp) / this->rho(p, T, pf, Xinp);
    scalar den = sqr(VV) + 2.0 * bmix * VV - sqr(bmix);
    scalar dpdt = RR * 1.0e-03 / (VV - bmix) - dadt / den;
    scalar dpdv = -RR * 1.0e-03 * T / sqr(VV - bmix) + 2.0 * amix * (VV + bmix) / sqr(den);
    scalar ct2 = -dpdv * mw_phase / sqr(rho_phase);
    return  T * sqr(dpdt) / (ct2 * sqr(rho_phase));
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::Cp // PY:phase mixture heat capability
(
    scalar p,
    scalar T,
    const label pf,
    const scalarList* Xinp) const
{
    scalar cp_phase;
    /*
    if (cvcpflag == 0) //ideal
    {
        scalar cpig_phase = 0.0, cvig_phase = 0.0;
        this->cvcpig_phase(T, comp_phase, cpig_phase, cvig_phase);
        cv_phase = cvig_phase;
    }
    else if (cvcpflag == 1) //real 1
    {
        */
    scalar h_phase = this->Ha(p, T, pf, Xinp);
    scalar deltat = 0.1;
    scalar T_plus = T + deltat;
    scalar h_phase_plus = this->Ha(p, T_plus, pf, Xinp);

    scalar delta_sie = h_phase_plus - h_phase;
    cp_phase = delta_sie / deltat;
    return cp_phase;
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::S(
    scalar p,
    scalar T,
    const label pf,
    const scalarList* Xinp) const //species entropy
{
    const scalar cc1 = 1.0 + sqrt(2.0);
    const scalar cc2 = 1.0 - sqrt(2.0);
    const scalar cc0 = 2.0 * sqrt(2.0);

    scalar bmix = this->b(p, T, Xinp);
    scalar VV = this->W(Xinp) / this->rho(p, T, pf, Xinp); // m3*mol-1,phase specific volume
    scalar xxx = (VV + cc1 * bmix) / (VV + cc2 * bmix);
    scalar dS_phase;
    if (xxx > 0.0)
    {
        dS_phase = RR * log((p * VV - bmix * p) / (RR * 1.0e-03 * T)) + this->dadt(p, T, Xinp) / (cc0 * bmix) * log(xxx);
    }
    else
    {
        dS_phase = 0.0;
    }

    dS_phase = dS_phase / this->W(Xinp);

    // Remark here, the ideal gas properties can be known, but how about the ideal liquid properites?????


    return dS_phase + this->Sideal(p, T, Xinp);
}
/*
template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::Sideal(
    scalar p,
    scalar T,
    const scalarList* Xinp //note, it is phase compositions
) const
{
    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
    {
        Xp = Xinp;
    }
    const scalarList& X = *Xp;

    scalar S_ideal = 0;
    forAll(X, i)
    {
        S_ideal += X[i] * (*this)[i].S(p, T);
    }
    return S_ideal;
}
*/
template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::Hig_phase(
    scalar T,
    const scalarList* Xinp) const
{
    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
    {
        Xp = Xinp;
    }
    const scalarList& X = *Xp;
    scalar H_ideal = 0;
    forAll(X, i)
    {
        H_ideal += X[i] * (*this)[i].Hig_phase(T);
    }
    return H_ideal;
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::Hig2_phase(
    scalar T,
    const scalarList* Xinp) const
{
    const scalarList* Xp = &this->X_;
    if (Xinp != NULL)
    {
        Xp = Xinp;
    }
    const scalarList& X = *Xp;
    scalar H_ideal = 0;
    forAll(X, i)
    {
        H_ideal += X[i] * (*this)[i].Hig2_phase(T);
    }

    return H_ideal * 4.18680;
}
// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// ************************************************************************* //
