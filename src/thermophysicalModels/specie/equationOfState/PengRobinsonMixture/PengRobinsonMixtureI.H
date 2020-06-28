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
#include "string.H"
#include "tmp.H"
using namespace Foam;
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class ThermoMixture>
inline Foam::PengRobinsonMixture<ThermoMixture>::PengRobinsonMixture(
    const word &name,
    PtrList<SingleThermoType> &speciesData,
    const speciesTable &specieNames,
    const dictionary &thermoDict)
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
    const dictionary &dict,PtrList<SingleThermoType> &speciesData)
    : ThermoMixture(dict,speciesData),
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
    const dictionary &dict)
{
    return autoPtr<PengRobinsonMixture<ThermoMixture>>(
        new PengRobinsonMixture<ThermoMixture>(dict));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::W // kg/mol
    (const scalarList *Xinp) const
{
    const scalarList *Xp = &this->X_;
    if (Xinp != NULL)
        Xp = Xinp;
    const scalarList &X = *Xp;

    const PengRobinsonMixture<ThermoMixture> &This = *this;
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
        const label *pf,
        const scalarList *Xinp) const
{                                                                //RR: 10e+03*m3 Pa K-1 mol-1
    return W(Xinp) * p / (Z(p, T, pf, Xinp) * RR * 1.0e-03 * T); // kg/m3
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::V //m3/mol
    (
        scalar p,
        scalar T,
        const label *pf,
        const scalarList *Xinp) const
{
    return W(Xinp) / rho(p, T, pf, Xinp);
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::A(
    scalar p,
    scalar T,
    const scalarList *Xinp) const
{
    const scalarList *Xp = &this->X_;
    if (Xinp != NULL)
        Xp = Xinp;
    const scalarList &X = *Xp;

    const PengRobinsonMixture<ThermoMixture> &This = *this;
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
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::B(
    scalar p,
    scalar T,
    const scalarList *Xinp) const
{
    const scalarList *Xp = &this->X_;
    if (Xinp != NULL)
        Xp = Xinp;
    const scalarList &X = *Xp;

    const PengRobinsonMixture<ThermoMixture> &This = *this;
    scalar Bmix = 0.0;
    forAll(X, i)
    {
        Bmix += X[i] * This[i].B(p, T);
    }
    return Bmix;
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::a(
    scalar p,
    scalar T,
    const scalarList *Xinp) const
{
    const scalarList *Xp = &this->X_;
    if (Xinp != NULL)
        Xp = Xinp;
    const scalarList &X = *Xp;

    const PengRobinsonMixture<ThermoMixture> &This = *this;
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
    const scalarList *Xinp) const
{
    const scalarList *Xp = &this->X_;
    if (Xinp != NULL)
        Xp = Xinp;
    const scalarList &X = *Xp;

    const PengRobinsonMixture<ThermoMixture> &This = *this;
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
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::d2adt2(
    scalar p,
    scalar T,
    const scalarList *Xinp) const
{
    const scalarList *Xp = &this->X_;
    if (Xinp != NULL)
        Xp = Xinp;
    const scalarList &X = *Xp;

    const PengRobinsonMixture<ThermoMixture> &This = *this;
    scalar d2adt2mix = 0.0, ai, aj, dai, daj;
    forAll(X, i)
    {
        forAll(X, j)
        {
            ai = This[i].EOS::a(p, T);
            aj = This[j].EOS::a(p, T);
            dai = This[i].dadt(p, T);
            daj = This[j].dadt(p, T);
            d2adt2mix += X[i] * X[j] * 0.25 / (ai * aj * sqrt(ai * aj)) * (-sqr(aj * dai + ai * daj) + 2 * ai * aj * (2 * dai * daj + ai * This[j].d2adt2(p, T) + aj * This[i].d2adt2(p, T))) * (1.0 - bico[i][j]);
        }
    }
    return d2adt2mix;
}
template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::b(
    scalar p,
    scalar T,
    const scalarList *Xinp) const
{
    const scalarList *Xp = &this->X_;
    if (Xinp != NULL)
        Xp = Xinp;
    const scalarList &X = *Xp;

    const PengRobinsonMixture<ThermoMixture> &This = *this;
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
        scalar r2 =
            qm * cos((th + 2 * constant::mathematical::pi) / 3) - a2 / 3;
        scalar r3 =
            qm * cos((th + 4 * constant::mathematical::pi) / 3) - a2 / 3;

        // three root;
        if (flag == 1) //gas
        {
            root = max(r1, max(r2, r3));
        }
        else //liquid
        {
            root = max(r1, max(r2, r3));
            root = min(r1 > 1e-3 ? r1 : root, root);
            root = min(r2 > 1e-3 ? r2 : root, root);
            root = min(r3 > 1e-3 ? r3 : root, root);
        }
    }
    else
    {
        // One root is real: in supercritical zone
        scalar D05 = sqrt(D);
        scalar S = sign(Rl + D05) * pow(mag(Rl + D05), 1.0 / 3.0);
        scalar Tl = sign(Rl - D05) * pow(mag(Rl - D05), 1.0 / 3.0);
        /*
        scalar RlD = Rl + D05;
        scalar S = 0.0;
        if (RlD > 0.0)
        {
            S = sign(Rl + D05) pow(mag(Rl + D05), 1.0 / 3.0);
        }
        else
        {
            S = -pow(-RlD, 1.0 / 3.0);
        }
        scalar Tl = 0;
        if (D05 > Rl)
        {
            Tl = -pow(mag(Rl - D05), 1.0 / 3.0);
        }
        else
        {
            Tl = pow(Rl - D05, 1.0 / 3.0);
        }
        */
        root = S + Tl - a2 / 3;
    }

    if (root <= 0.0)
        FatalErrorInFunction << "Z is negative! Z = " << root << exit(FatalError);
    return root;
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::Z(
    scalar p,
    scalar T,
    const label *pf,
    const scalarList *Xinp) const
{
    scalar A = this->A(p, T, Xinp), B = this->B(p, T, Xinp);
    scalar a2 = B - 1.0;
    scalar a1 = A - 2.0 * B - 3.0 * B * B;
    scalar a0 = -A * B + B * B + B * B * B;
    if (pf != NULL)
        return cubicSolve(a2, a1, a0, *pf);
    else
        return cubicSolve(a2, a1, a0, phaseFlag);
}

template <class ThermoMixture>
inline Foam::autoPtr<scalarList> Foam::PengRobinsonMixture<ThermoMixture>::A2(
    scalar p, //inp
    scalar T, //inp
    const scalarList *Xinp) const
{
    const scalarList *Xp = &this->X_;
    if (Xinp != NULL)
    {
        Xp = Xinp;
    }
    const scalarList &X = *Xp;

    scalarList &Amix2 = *(new scalarList(this->N_));

    forAll(X, i)
    {
        Amix2[i] = 0.0;
        forAll(X, j)
        {
            Amix2[i] += X[j] * sqrt((*this)[i].EOS::A(p, T) * (*this)[j].EOS::A(p, T)) * (1.0 - bico[i][j]);
        }
    }
    return autoPtr<scalarList>(&Amix2);
}

// PY: fugacity coefficient of species in each phase
template <class ThermoMixture>
inline Foam::autoPtr<scalarList> Foam::PengRobinsonMixture<ThermoMixture>::fugacityCoefficient(
    scalar p, //input
    scalar T,
    const label *pf,
    const scalarList *Xinp) const
{
    label pflag = this->phaseFlag;
    if (pf != NULL)
        pflag = *pf;

    const scalarList *Xp = &this->X_;
    if (Xinp != NULL)
    {
        Xp = Xinp;
    }
    const scalarList &X = *Xp;

    scalarList &fugcoef = *(new scalarList(this->N_));
    scalar Amixz = this->A(p, T, Xinp);
    scalar Bmixz = this->B(p, T, Xinp);
    scalar Zmix = 0.0;
    autoPtr<scalarList> pAmix2(A2(p, T, Xinp));
    const scalarList &Amix2 = pAmix2();
    const scalar small8 = 1.e-08; //PY: need to be noticed

    //this->ABZmix(p, T, comp_inp, Amixz, Bmixz, Amix2, Bz);

    Zmix = this->Z(p, T, pf, Xinp);

    label temp = 1;
    scalar Zmix0 = this->Z(p, T, &temp, Xinp);

    scalar comptrifile = 0.7;
    scalarList Xs(this->X_.size());

    if (pflag == 0)
    {
        while (Zmix < 0.0)
        {
            Xs[0] = comptrifile;
            for (int i = 1; i < this->X_.size(); i++)
                Xs[i] = (1.0 - comptrifile) / (Xs.size() - 1);
            Zmix = this->Z(p, T, &pflag, &Xs);
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

    scalar zfactor00 = Zmix;//0.0020455756600868957
    scalar ZMB = zfactor00 - Bmixz;//0.0002558026230427609

    if (ZMB < 0.0)
    {
        zfactor00 = Bmixz + small8;
        ZMB = small8;
    }
    forAll(X, spid)
    {
        scalar b = (*this)[spid].B(p, T);//0.0019591145739235173,0.0013946427843255767
        fugcoef[spid] = exp(b / Bmixz * (zfactor00 - 1.0) - log(ZMB) - Amixz / (cc0 * Bmixz) * (2.0 * Amix2[spid] / Amixz - b / Bmixz) * log((zfactor00 + cc1 * Bmixz) / (zfactor00 + cc2 * Bmixz)));//1.4644880820052337,
    }
    return autoPtr<scalarList>(&fugcoef);
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
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::Ha // internal energy with unit of J/kg
    (
        scalar p,
        scalar T,
        const label *pf,
        const scalarList *Xinp) const
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
    }
    else
    {
        dH = (Z(p, T, pf, Xinp) - 1) * RR * 1e-3 * T;
    }

    dH = dH / this->W(Xinp); // J/kg
    // Based on Eq.12, Internal energy is estimated based on ideal(close to absolute 0 T) + departure
    //Info<<"!!!dH="<<dH<<"!!!Hi="<<this->Hideal(p, T, Xinp)<<"!!!\n";
    return dH + this->Hideal(p, T, Xinp); //+11990129.829480369 +3991.2*1000;
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::Hc // internal energy with unit of J/kg
    (
        const scalarList *Xinp) const
{
    return this->Hideal(1, Tstd, Xinp); //+11990129.829480369 +3991.2*1000;
}

//temp E

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::Hideal(
    scalar p,
    scalar T,
    const scalarList *Xinp //note, it is phase compositions
) const
{
    const scalarList *Xp = &this->X_;
    if (Xinp != NULL)
    {
        Xp = Xinp;
    }
    const scalarList &X = *Xp;

    scalar H_ideal = 0,sum=0;
    forAll(X, i)
    {
        sum+= X[i]*(*this)[i].W();
        H_ideal += X[i] * (*this)[i].Ha(p, T) *(*this)[i].W();
    }
    return H_ideal/sum;
}

template<class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::CpMCv
(
        scalar p,
        scalar T,
        const label *pf,
        const scalarList *Xinp) const
{
    scalar amix = this->a(p, T, Xinp);
    scalar bmix = this->b(p, T, Xinp);
    scalar dadt = this->dadt(p, T, Xinp);

    scalar mw_phase = this->W(Xinp);
    scalar rho_phase = this->rho(p, T, pf, Xinp);
    scalar VV = this->W(Xinp) / this->rho(p, T, pf, Xinp);
    scalar den  = sqr(VV) + 2.0*bmix*VV-sqr(bmix);
    scalar dpdt = RR*1.0e-03/(VV-bmix) - dadt/den;
    scalar dpdv = -RR*1.0e-03*T/sqr(VV-bmix) + 2.0*amix*(VV+bmix)/sqr(den); 
    scalar ct2  = -dpdv*mw_phase/sqr(rho_phase);
    return  T*sqr(dpdt)/(ct2*sqr(rho_phase));
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::Cp // PY:phase mixture heat capability
    (
        scalar p,
        scalar T,
        const label *pf,
        const scalarList *Xinp) const
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
    const label *pf,
    const scalarList *Xinp) const //species entropy
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


    return dS_phase + this->Sideal(p,T, Xinp);
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::Sideal(
    scalar p,
    scalar T,
    const scalarList *Xinp //note, it is phase compositions
) const
{
    const scalarList *Xp = &this->X_;
    if (Xinp != NULL)
    {
        Xp = Xinp;
    }
    const scalarList &X = *Xp;

    scalar S_ideal = 0;
    forAll(X, i)
    {
        S_ideal += X[i] * (*this)[i].S(p, T);
    }
    return S_ideal;
}

template <class ThermoMixture>
inline Foam::scalar Foam::PengRobinsonMixture<ThermoMixture>::Hig_phase(
    scalar T,
    const scalarList *Xinp) const
{
    const scalarList *Xp = &this->X_;
    if (Xinp != NULL)
    {
        Xp = Xinp;
    }
    const scalarList &X = *Xp;
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
    const scalarList *Xinp) const
{
    const scalarList *Xp = &this->X_;
    if (Xinp != NULL)
    {
        Xp = Xinp;
    }
    const scalarList &X = *Xp;
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
