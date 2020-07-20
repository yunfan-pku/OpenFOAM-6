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

#include "multithermo.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class ThermoMixture, template <class> class Type>
inline Foam::species::multithermo<ThermoMixture, Type>::multithermo(
    const ThermoMixture &sp)
    : ThermoMixture(sp)
{
}
template <class ThermoMixture, template <class> class Type>
inline Foam::species::multithermo<ThermoMixture, Type>::multithermo(
    const word &name,
    PtrList<SingleThermoType> &speciesData,
    const speciesTable &specieNames,
    const dictionary &thermoDict)
    : ThermoMixture(name, speciesData, specieNames, thermoDict)
{
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar Foam::species::multithermo<ThermoMixture, Type>::T(
    scalar f,
    scalar p,
    scalar T0,
    scalar (multithermo<ThermoMixture, Type>::*F)(const scalar, const scalar) const,
    scalar (multithermo<ThermoMixture, Type>::*dFdT)(const scalar, const scalar)
        const,
    scalar (multithermo<ThermoMixture, Type>::*limit)(const scalar) const) const
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
    scalar F1=(this->*F)(p, T1), F2=F1;
    scalar T2 =(this->*limit)(T0+1);
   // if(F1>f)
   //     T2=(this->*limit)(T0-30);
    //scalar F1=(this->*F)(p, T1);
    //scalar F2=(this->*F)(p, T2);
    //scalar Tnew = T0;
    scalar Ttol = T0 * tol_;
    scalar dF=0,dfdt=0;
    int iter = 0;
    scalar relax=1;
    do
    {
        if(fabs(F2)>fabs(F1))
            relax*=0.5;
        F1=F2;
        F2=(this->*F)(p, T2);

        dF= (F2 - f);
        dfdt=(F2 - F1)/(T2-T1);//(this->*dFdT)(p, Test);
        T1 = T2;
        if (mag(relax*dF / dfdt)>30)
        T2=  (this->*limit)(T1 - sign(relax*dF / dfdt)*30);
        else 
        T2=  (this->*limit)(T1 - relax*dF / dfdt);
        //Tnew =
         //   (this->*limit)(Test - ((this->*F)(p, Test) - f) / (this->*dFdT)(p, Test));
        
        if(T2==0)
        {
             FatalErrorInFunction<<setprecision(50)
                << "Maximum number of iterations exceeded: " << maxIter_
                << " f: " << f
                << " p: " << p
                << " T0: " << T0
                << "  Tg: " << T2
                <<" X: "<<this->X_
                << abort(FatalError);
        }
        if (iter++ > maxIter_)
        {
            FatalErrorInFunction<<setprecision(50)
                << "Maximum number of iterations exceeded: " << maxIter_
                << " f: " << f
                << " p: " << p
                << " T0: " << T0
                << "  Tg: " << T2
                <<" X: "<<this->X_
                << abort(FatalError);
        }

    } while (mag(T2 - T1) > Ttol && mag(f-F2) >0.0001 );
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
inline Foam::scalar Foam::species::multithermo<ThermoMixture, Type>::T(
    scalar f,
    scalar p,
    scalar T0,
    scalar (multithermo<ThermoMixture, Type>::*F)(const scalar, const scalar) const,
    scalar (multithermo<ThermoMixture, Type>::*dFdT)(const scalar, const scalar)
        const,
    scalar (multithermo<ThermoMixture, Type>::*limit)(const scalar) const) const
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
inline Foam::species::multithermo<ThermoMixture, Type>::multithermo(
    const word &name,
    const multithermo &st)
    : ThermoMixture(name, st)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class ThermoMixture, template <class> class Type>
inline Foam::word
Foam::species::multithermo<ThermoMixture, Type>::heName()
{
    return Type<multithermo<ThermoMixture, Type>>::energyName();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::multithermo<ThermoMixture, Type>::Cv(const scalar p, const scalar T) const
{
    return this->Cp(p, T) - this->CpMCv(p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::multithermo<ThermoMixture, Type>::Cpv(const scalar p, const scalar T) const
{
    return Type<multithermo<ThermoMixture, Type>>::Cpv(*this, p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::multithermo<ThermoMixture, Type>::gamma(const scalar p, const scalar T) const
{
    const scalar Cp = this->Cp(p, T);
    return Cp / (Cp - this->CpMCv(p, T));
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::multithermo<ThermoMixture, Type>::CpByCpv(
    const scalar p,
    const scalar T) const
{
    return Type<multithermo<ThermoMixture, Type>>::CpByCpv(*this, p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::multithermo<ThermoMixture, Type>::HE(const scalar p, const scalar T) const
{
    return Type<multithermo<ThermoMixture, Type>>::HE(*this, p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::multithermo<ThermoMixture, Type>::Es(const scalar p, const scalar T) const
{
    return this->Hs(p, T) - p / this->rho(p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::multithermo<ThermoMixture, Type>::Ea(const scalar p, const scalar T) const
{
    return this->Ha(p, T) - p / this->rho(p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::multithermo<ThermoMixture, Type>::G(const scalar p, const scalar T) const
{
    return this->Ha(p, T) - T * this->S(p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::multithermo<ThermoMixture, Type>::A(const scalar p, const scalar T) const
{
    return this->Ea(p, T) - T * this->S(p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::multithermo<ThermoMixture, Type>::cp(const scalar p, const scalar T) const
{
    return this->Cp(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::multithermo<ThermoMixture, Type>::ha(const scalar p, const scalar T) const
{
    return this->Ha(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::multithermo<ThermoMixture, Type>::hs(const scalar p, const scalar T) const
{
    return this->Hs(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::multithermo<ThermoMixture, Type>::hc() const
{
    return this->Hc() * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::multithermo<ThermoMixture, Type>::s(const scalar p, const scalar T) const
{
    return this->S(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::multithermo<ThermoMixture, Type>::he(const scalar p, const scalar T) const
{
    return this->HE(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::multithermo<ThermoMixture, Type>::cv(const scalar p, const scalar T) const
{
    return this->Cv(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::multithermo<ThermoMixture, Type>::es(const scalar p, const scalar T) const
{
    return this->Es(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::multithermo<ThermoMixture, Type>::ea(const scalar p, const scalar T) const
{
    return this->Ea(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::multithermo<ThermoMixture, Type>::g(const scalar p, const scalar T) const
{
    return this->G(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::multithermo<ThermoMixture, Type>::a(const scalar p, const scalar T) const
{
    return this->A(p, T) * this->W();
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::multithermo<ThermoMixture, Type>::K(const scalar p, const scalar T) const
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
Foam::species::multithermo<ThermoMixture, Type>::Kp(const scalar p, const scalar T) const
{
    return K(p, T);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::multithermo<ThermoMixture, Type>::Kc(const scalar p, const scalar T) const
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
inline Foam::scalar Foam::species::multithermo<ThermoMixture, Type>::Kx(
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
inline Foam::scalar Foam::species::multithermo<ThermoMixture, Type>::Kn(
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
inline Foam::scalar Foam::species::multithermo<ThermoMixture, Type>::THE(
    const scalar he,
    const scalar p,
    const scalar T0) const
{
    return Type<multithermo<ThermoMixture, Type>>::THE(*this, he, p, T0);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar Foam::species::multithermo<ThermoMixture, Type>::THs(
    const scalar hs,
    const scalar p,
    const scalar T0) const
{
    return T(
        hs,
        p,
        T0,
        &multithermo<ThermoMixture, Type>::Hs,
        &multithermo<ThermoMixture, Type>::Cp,
        &multithermo<ThermoMixture, Type>::limit);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar Foam::species::multithermo<ThermoMixture, Type>::THa(
    const scalar ha,
    const scalar p,
    const scalar T0) const
{
    return T(
        ha,
        p,
        T0,
        &multithermo<ThermoMixture, Type>::Ha,
        &multithermo<ThermoMixture, Type>::Cp,
        &multithermo<ThermoMixture, Type>::limit);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar Foam::species::multithermo<ThermoMixture, Type>::TEs(
    const scalar es,
    const scalar p,
    const scalar T0) const
{
    return T(
        es,
        p,
        T0,
        &multithermo<ThermoMixture, Type>::Es,
        &multithermo<ThermoMixture, Type>::Cv,
        &multithermo<ThermoMixture, Type>::limit);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar Foam::species::multithermo<ThermoMixture, Type>::TEa(
    const scalar ea,
    const scalar p,
    const scalar T0) const
{
    return T(
        ea,
        p,
        T0,
        &multithermo<ThermoMixture, Type>::Ea,
        &multithermo<ThermoMixture, Type>::Cv,
        &multithermo<ThermoMixture, Type>::limit);
}

template <class ThermoMixture, template <class> class Type>
inline Foam::scalar
Foam::species::multithermo<ThermoMixture, Type>::dKcdTbyKc(
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
Foam::species::multithermo<ThermoMixture, Type>::dcpdT(const scalar p, const scalar T) const
{
    return this->dCpdT(p, T) * this->W();
    ;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template <class ThermoMixture, template <class> class Type>
inline void Foam::species::multithermo<ThermoMixture, Type>::operator+=(
    const multithermo<ThermoMixture, Type> &st)
{
    ThermoMixture::operator+=(st);
}

template <class ThermoMixture, template <class> class Type>
inline void Foam::species::multithermo<ThermoMixture, Type>::operator*=(const scalar s)
{
    ThermoMixture::operator*=(s);
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template <class ThermoMixture, template <class> class Type>
inline Foam::species::multithermo<ThermoMixture, Type> Foam::species::operator+(
    const multithermo<ThermoMixture, Type> &st1,
    const multithermo<ThermoMixture, Type> &st2)
{
    return multithermo<ThermoMixture, Type>(
        static_cast<const ThermoMixture &>(st1) + static_cast<const ThermoMixture &>(st2));
}

template <class ThermoMixture, template <class> class Type>
inline Foam::species::multithermo<ThermoMixture, Type> Foam::species::operator*(
    const scalar s,
    const multithermo<ThermoMixture, Type> &st)
{
    return multithermo<ThermoMixture, Type>(
        s * static_cast<const ThermoMixture &>(st));
}

template <class ThermoMixture, template <class> class Type>
inline Foam::species::multithermo<ThermoMixture, Type> Foam::species::operator==(
    const multithermo<ThermoMixture, Type> &st1,
    const multithermo<ThermoMixture, Type> &st2)
{
    return multithermo<ThermoMixture, Type>(
        static_cast<const ThermoMixture &>(st1) == static_cast<const ThermoMixture &>(st2));
}

// ************************************************************************* //
