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

#include "ISATmultithermo.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class ThermoMixture, template <class> class Type>
inline Foam::species::ISATmultithermo<ThermoMixture, Type>::ISATmultithermo(
    const ThermoMixture& sp)
    : ThermoMixture(sp)
{
}
template <class ThermoMixture, template <class> class Type>
inline Foam::species::ISATmultithermo<ThermoMixture, Type>::ISATmultithermo(
    const word& name,
    PtrList<SingleThermoType>& speciesData,
    const speciesTable& specieNames,
    const dictionary& thermoDict)
    : ThermoMixture(name, speciesData, specieNames, thermoDict),
    Hnfun(this), psifun(this, &ThermoMixture::psi), mufun(this, &ThermoMixture::mu),
    alphahfun(this, &ThermoMixture::alphah), vaporfrafun(this, &ThermoMixture::vaporfra),
    Ttree(speciesData.size() + 1, 1, Hnfun), psitree(speciesData.size() + 1, 1, psifun),//,Ttree(specieNames.size()+1,1),psimualphahtree(specieNames.size()+1,3)
    mutree(speciesData.size() + 1, 1, mufun), alphahtree(speciesData.size() + 1, 1, alphahfun),
    vaporfratree(speciesData.size() + 1, 1, vaporfrafun)
{
    psitree.epsilon() = 1e-8;
    psitree.relepsilon() = 5e-2;
    Ttree.relepsilon() = 2e-3;
    mutree.relepsilon() = 2e-3;
    alphahtree.relepsilon() = 2e-3;
    vaporfratree.epsilon() = 1e-3;
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar Foam::species::ISATmultithermo<ThermoMixture, Type>::T(
    scalar f,
    scalar p,
    scalar T0,
    scalar(ISATmultithermo<ThermoMixture, Type>::* F)(const scalar, const scalar) const,
    scalar(ISATmultithermo<ThermoMixture, Type>::* dFdT)(const scalar, const scalar)
    const,
    scalar(ISATmultithermo<ThermoMixture, Type>::* limit)(const scalar) const) const
{

    if (T0 < 0)
    {
        FatalErrorInFunction
            << "Negative initial temperature T0: " << T0
            << abort(FatalError);
    }
    //static int call_n=0; 
    //static int inter_n=0;
    scalar T1 = T0;
    scalar F1 = (this->*F)(p, T1), F2 = F1;
    scalar T2 = (this->*limit)(T0 + 1);
    // if(F1>f)
    //     T2=(this->*limit)(T0-30);
     //scalar F1=(this->*F)(p, T1);
     //scalar F2=(this->*F)(p, T2);
     //scalar Tnew = T0;
    scalar Ttol = T0 * tol_;
    scalar dF = 0, dfdt = 0;
    int iter = 0;
    scalar relax = 1;
    do
    {
        if (fabs(F2) > fabs(F1))
            relax *= 0.5;
        F1 = F2;
        F2 = (this->*F)(p, T2);

        dF = (F2 - f);
        dfdt = (F2 - F1) / (T2 - T1);//(this->*dFdT)(p, Test);
        T1 = T2;
        if (mag(relax * dF / dfdt) > 30)
            T2 = (this->*limit)(T1 - sign(relax * dF / dfdt) * 30);
        else
            T2 = (this->*limit)(T1 - relax * dF / dfdt);
        //Tnew =
         //   (this->*limit)(Test - ((this->*F)(p, Test) - f) / (this->*dFdT)(p, Test));

        if (T2 == 0)
        {
            FatalErrorInFunction << setprecision(50)
                << "Maximum number of iterations exceeded: " << maxIter_
                << " f: " << f
                << " p: " << p
                << " T0: " << T0
                << "  Tg: " << T2
                << " X: " << this->X_
                << abort(FatalError);
        }
        if (iter++ > maxIter_)
        {
            FatalErrorInFunction << setprecision(50)
                << "Maximum number of iterations exceeded: " << maxIter_
                << " f: " << f
                << " p: " << p
                << " T0: " << T0
                << "  Tg: " << T2
                << " X: " << this->X_
                << abort(FatalError);
        }

    } while (mag(T2 - T1) > Ttol && mag(f - F2) > 0.0001);
    //call_n++;
    //inter_n+=iter;
    //if(call_n%10000==0)
    //Info<<"-------------------call N="<<call_n<<",average iter="<<inter_n*1.0/call_n<<endl;
    return T2;


    /*IQI wrong
       if (T0 < 0)
        {
            FatalErrorInFunction
                << "Negative initial temperature T0: " << T0
                << abort(FatalError);
        }

        static int call_n=0;
        static int inter_n=0;
        scalar T1 = T0;
        scalar T2 = T0+10;
        scalar T3 = T0+20;

        scalar F1=0, F2=0,F3=0;
        F2=(this->*F)(p, T1)-f;
        if(F2==0)
        return T1;
        F3=(this->*F)(p, T2)-f;
        if(F3==0)
        return T2;
        F1=F2;

        scalar Ttol = T0 * tol_;
        int iter = 0;
        scalar r,q,s;

        do
        {

            F1=F2;
            F2=F3;
            F3=(this->*F)(p, T3)-f;
            if(F3==0)
            return T3;

            q=F1/F2;
            r=F3/F2;
            s=F3/F1;

            T1 = T2;
            T2 = T3;
            T3=  T3 - (r*(r-q)*(T3-T2)+(1-r)*s*(T3-T1))/((q-1)*(r-1)*(s-1));

            if (iter++ > maxIter_)
            {
                FatalErrorInFunction
                    << "Maximum number of iterations exceeded: " << maxIter_
                    << " f: " << f
                    << " p: " << p
                    << " T0: " << T0
                    << "  Tg: " << T2
                    <<" X: "<<this->X_
                    << abort(FatalError);
            }

        } while (mag(T3 - T2) > Ttol);
        call_n++;
        inter_n+=iter;
        if(call_n%10000==0)
        Info<<"-------------------call N="<<call_n<<",average iter="<<inter_n*1.0/call_n<<endl;

        return T3;
        */
}
/*
template <class ThermoMixture, template <class> class Type>
inline Foam::scalar Foam::species::ISATmultithermo<ThermoMixture, Type>::T(
    scalar f,
    scalar p,
    scalar T0,
    scalar (ISATmultithermo<ThermoMixture, Type>::*F)(const scalar, const scalar) const,
    scalar (ISATmultithermo<ThermoMixture, Type>::*dFdT)(const scalar, const scalar)
        const,
    scalar (ISATmultithermo<ThermoMixture, Type>::*limit)(const scalar) const) const
{
    if (T0 < 0)
    {
        FatalErrorInFunction
            << "Negative initial temperature T0: " << T0
            << abort(FatalError);
    }
    static int call_n=0;
    static int inter_n=0;
    scalar T1 = T0-10;
    scalar T2 = T0;
    scalar T3 = T0+10;
    scalar F1=0, F2=0,F3=0;
    F2=(this->*F)(p, T1);
    F3=(this->*F)(p, T2);
    F1=F2;
    scalar Ttol = T0 * tol_;
    int iter = 0;
    scalar r,q,s;
    do
    {
        F1=F2;
        F2=F3;
        F3=(this->*F)(p, T3);

        q=F1/F2;
        r=F3/F2;
        s=F3/F1;

        T1 = T2;
        T2 = T3;
        T3=  T3 - (r*(r-q)*(T3-T2)+(1-r)*s*(T3-T1))/((q-1)*(r-1)*(s-1));
        if (iter++ > maxIter_)
        {
            FatalErrorInFunction
                << "Maximum number of iterations exceeded: " << maxIter_
                << " f: " << f
                << " p: " << p
                << " T0: " << T0
                << "  Tg: " << T2
                <<" X: "<<this->X_
                << abort(FatalError);
        }

    } while (mag(T3 - T2) > Ttol && iter++<30);
    call_n++;
    inter_n+=iter;
    if(call_n%10000==0)
    Info<<"-------------------call N="<<call_n<<",average iter="<<inter_n*1.0/call_n<<endl;
    return T3;
}
*/
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class ThermoMixture, template <class> class Type>
inline Foam::species::ISATmultithermo<ThermoMixture, Type>::ISATmultithermo(
    const word& name,
    const ISATmultithermo& st)
    : ThermoMixture(name, st)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class ThermoMixture, template <class> class Type>
inline Foam::word
Foam::species::ISATmultithermo<ThermoMixture, Type>::heName()
{
    return Type<ISATmultithermo<ThermoMixture, Type>>::energyName();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::Cv(const scalar p, const scalar T) const
{
    return this->Cp(p, T) - this->CpMCv(p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::Cpv(const scalar p, const scalar T) const
{
    return Type<ISATmultithermo<ThermoMixture, Type>>::Cpv(*this, p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::gamma(const scalar p, const scalar T) const
{
    const scalar Cp = this->Cp(p, T);
    return Cp / (Cp - this->CpMCv(p, T));
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::CpByCpv(
    const scalar p,
    const scalar T) const
{
    return Type<ISATmultithermo<ThermoMixture, Type>>::CpByCpv(*this, p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::HE(const scalar p, const scalar T) const
{
    return Type<ISATmultithermo<ThermoMixture, Type>>::HE(*this, p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::Es(const scalar p, const scalar T) const
{
    return this->Hs(p, T) - p / this->rho(p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::Ea(const scalar p, const scalar T) const
{
    return this->Ha(p, T) - p / this->rho(p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::G(const scalar p, const scalar T) const
{
    return this->Ha(p, T) - T * this->S(p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::A(const scalar p, const scalar T) const
{
    return this->Ea(p, T) - T * this->S(p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::cp(const scalar p, const scalar T) const
{
    return this->Cp(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::ha(const scalar p, const scalar T) const
{
    return this->Ha(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::hs(const scalar p, const scalar T) const
{
    return this->Hs(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::hc() const
{
    return this->Hc() * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::s(const scalar p, const scalar T) const
{
    return this->S(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::he(const scalar p, const scalar T) const
{
    return this->HE(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::cv(const scalar p, const scalar T) const
{
    return this->Cv(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::es(const scalar p, const scalar T) const
{
    return this->Es(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::ea(const scalar p, const scalar T) const
{
    return this->Ea(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::g(const scalar p, const scalar T) const
{
    return this->G(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::a(const scalar p, const scalar T) const
{
    return this->A(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::K(const scalar p, const scalar T) const
{
    scalar arg = -1; //this->Y()*this->G(Pstd, T)/(RR*T);

    if (arg < 600)
    {
        return exp(arg);
    }
    else
    {
        return rootVGreat;
    }
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::Kp(const scalar p, const scalar T) const
{
    return K(p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::Kc(const scalar p, const scalar T) const
{
    const scalar nm = 1; //this->Y()/this->W();

    if (equal(nm, small))
    {
        return Kp(p, T);
    }
    else
    {
        return Kp(p, T) * pow(Pstd / (RR * T), nm);
    }
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar Foam::species::ISATmultithermo<ThermoMixture, Type>::Kx(
    const scalar p,
    const scalar T) const
{
    const scalar nm = this->Y() / this->W();

    if (equal(nm, small))
    {
        return Kp(p, T);
    }
    else
    {
        return Kp(p, T) * pow(Pstd / p, nm);
    }
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar Foam::species::ISATmultithermo<ThermoMixture, Type>::Kn(
    const scalar p,
    const scalar T,
    const scalar n) const
{
    const scalar nm = this->Y() / this->W();

    if (equal(nm, small))
    {
        return Kp(p, T);
    }
    else
    {
        return Kp(p, T) * pow(n * Pstd / p, nm);
    }
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar Foam::species::ISATmultithermo<ThermoMixture, Type>::THE(
    const scalar he,
    const scalar p,
    const scalar T0) const
{
    return Type<ISATmultithermo<ThermoMixture, Type>>::THE(*this, he, p, T0);
}
//This one !!!
template <class ThermoMixture, template <class> class Type>
inline Foam::scalar Foam::species::ISATmultithermo<ThermoMixture, Type>::THs(
    const scalar hs,
    const scalar p,
    const scalar T0) const
{
    double sum = 0;
    scalarList in(this->X().size() + 2), out(1);
    for (int i = 0;i < this->X().size() - 1;i++)
    {
        in[i] = this->X_[i];
        sum += this->X_[i];
    }
    in[this->X().size() - 1] = hs;
    in[this->X().size()] = p;
    in[this->X().size() + 1] = T0;

    if (1 - sum < 1e-5)
    {
        Ttree.notCall++;
        return T(
            hs,
            p,
            T0,
            &ISATmultithermo<ThermoMixture, Type>::Hs,
            &ISATmultithermo<ThermoMixture, Type>::Cp,
            &ISATmultithermo<ThermoMixture, Type>::limit);
    }

    //Hnfun.value(in, out);
    Ttree.call(in, out);

    return out[0];
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar Foam::species::ISATmultithermo<ThermoMixture, Type>::THa(
    const scalar ha,
    const scalar p,
    const scalar T0) const
{
    return T(
        ha,
        p,
        T0,
        &ISATmultithermo<ThermoMixture, Type>::Ha,
        &ISATmultithermo<ThermoMixture, Type>::Cp,
        &ISATmultithermo<ThermoMixture, Type>::limit);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar Foam::species::ISATmultithermo<ThermoMixture, Type>::TEs(
    const scalar es,
    const scalar p,
    const scalar T0) const
{
    return T(
        es,
        p,
        T0,
        &ISATmultithermo<ThermoMixture, Type>::Es,
        &ISATmultithermo<ThermoMixture, Type>::Cv,
        &ISATmultithermo<ThermoMixture, Type>::limit);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar Foam::species::ISATmultithermo<ThermoMixture, Type>::TEa(
    const scalar ea,
    const scalar p,
    const scalar T0) const
{
    return T(
        ea,
        p,
        T0,
        &ISATmultithermo<ThermoMixture, Type>::Ea,
        &ISATmultithermo<ThermoMixture, Type>::Cv,
        &ISATmultithermo<ThermoMixture, Type>::limit);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::dKcdTbyKc(
    const scalar p,
    const scalar T) const
{
    const scalar nm = 1; //this->Y()/this->W();

    if (equal(nm, small))
    {
        return -this->dGdT(Pstd, T); //*this->Y()/RR;
    }
    else
    {
        return -(nm / T + this->dGdT(Pstd, T)); //*this->Y()/RR);
    }
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::dcpdT(const scalar p, const scalar T) const
{
    return this->dCpdT(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::psi(const scalar p, const scalar T) const
{

    double sum = 0;
    scalarList in(this->X().size() + 1), out(1);
    for (int i = 0;i < this->X().size() - 1;i++)
    {
        in[i] = this->X_[i];
        sum += this->X_[i];
    }
    in[this->X().size() - 1] = p;
    in[this->X().size()] = T;

    if (1 - sum < 1e-5)
    {
        psitree.notCall++;
        return this->ThermoMixture::psi(p, T);
    }
    // Info<<"psihaha"<<endl;

    psitree.call(in, out);
    return out[0];

}
template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::mu(const scalar p, const scalar T) const
{

    double sum = 0;
    scalarList in(this->X().size() + 1), out(1);
    for (int i = 0;i < this->X().size() - 1;i++)
    {
        in[i] = this->X_[i];
        sum += this->X_[i];
    }
    in[this->X().size() - 1] = p;
    in[this->X().size()] = T;

    if (1 - sum < 1e-5)
    {
        mutree.notCall++;
        return this->ThermoMixture::mu(p, T);
    }
    // Info<<"psihaha"<<endl;

    mutree.call(in, out);
    return out[0];
}
template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::alphah(const scalar p, const scalar T) const
{

    double sum = 0;
    scalarList in(this->X().size() + 1), out(1);
    for (int i = 0;i < this->X().size() - 1;i++)
    {
        in[i] = this->X_[i];
        sum += this->X_[i];
    }
    in[this->X().size() - 1] = p;
    in[this->X().size()] = T;

    if (1 - sum < 1e-5)
    {
        alphahtree.notCall++;
        return this->ThermoMixture::alphah(p, T);
    }
    // Info<<"psihaha"<<endl;

    alphahtree.call(in, out);
    return out[0];

}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::ISATmultithermo<ThermoMixture, Type>::vaporfra(const scalar p, const scalar T) const
{
    /*
    double sum = 0;
    scalarList in(this->X().size() + 1), out(1);
    for (int i = 0;i < this->X().size() - 1;i++)
    {
        in[i] = this->X_[i];
        sum += this->X_[i];
    }
    in[this->X().size() - 1] = p;
    in[this->X().size()] = T;

    if (1 - sum < 1e-5)*/
    return this->ThermoMixture::vaporfra(p, T);
    // Info<<"psihaha"<<endl;
    /*
    vaporfratree.call(in, out);
    return out[0];
    */
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template <class ThermoMixture, template <class> class Type>
inline void Foam::species::ISATmultithermo<ThermoMixture, Type>::operator+=(
    const ISATmultithermo<ThermoMixture, Type>& st)
{
    ThermoMixture::operator+=(st);
}

template <class ThermoMixture, template <class> class Type>
inline void Foam::species::ISATmultithermo<ThermoMixture, Type>::operator*=(const scalar s)
{
    ThermoMixture::operator*=(s);
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template <class ThermoMixture, template <class> class Type>
inline Foam::species::ISATmultithermo<ThermoMixture, Type> Foam::species::operator+(
    const ISATmultithermo<ThermoMixture, Type>& st1,
    const ISATmultithermo<ThermoMixture, Type>& st2)
{
    return ISATmultithermo<ThermoMixture, Type>(
        static_cast<const ThermoMixture&>(st1) + static_cast<const ThermoMixture&>(st2));
}

template <class ThermoMixture, template <class> class Type>
inline Foam::species::ISATmultithermo<ThermoMixture, Type> Foam::species::operator*(
    const scalar s,
    const ISATmultithermo<ThermoMixture, Type>& st)
{
    return ISATmultithermo<ThermoMixture, Type>(
        s * static_cast<const ThermoMixture&>(st));
}

template <class ThermoMixture, template <class> class Type>
inline Foam::species::ISATmultithermo<ThermoMixture, Type> Foam::species::operator==(
    const ISATmultithermo<ThermoMixture, Type>& st1,
    const ISATmultithermo<ThermoMixture, Type>& st2)
{
    return ISATmultithermo<ThermoMixture, Type>(
        static_cast<const ThermoMixture&>(st1) == static_cast<const ThermoMixture&>(st2));
}

// ************************************************************************* //
