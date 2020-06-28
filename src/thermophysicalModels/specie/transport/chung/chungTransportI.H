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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class Thermo>
inline Foam::chungTransport<Thermo>::chungTransport(
    const Thermo &t)
    : Thermo(t)
{
}

template <class Thermo>
inline Foam::chungTransport<Thermo>::chungTransport(
    const word &name,
    const chungTransport &ct)
    : Thermo(name, ct)
{
}

template <class Thermo>
inline Foam::autoPtr<Foam::chungTransport<Thermo>>
Foam::chungTransport<Thermo>::clone() const
{
    return autoPtr<chungTransport<Thermo>>(
        new chungTransport<Thermo>(*this));
}

template <class Thermo>
inline Foam::autoPtr<Foam::chungTransport<Thermo>>
Foam::chungTransport<Thermo>::New(
    const dictionary &dict)
{
    return autoPtr<chungTransport<Thermo>>(
        new chungTransport<Thermo>(dict));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*
template <class Thermo>
inline Foam::scalar Foam::chungTransport<Thermo>::mu(
    scalar p,
    scalar T,
    const scalarList *Xinp) const
{
    const scalarList *Xp = &this->X_;
    if (Xinp != NULL)
    {
        Xp = Xinp;
    }
    const scalarList &X = *Xp;

    scalarListList kij(this->N_);
    scalarListList mwij(this->N_);
    scalarListList w(this->N_);
    scalarListList sigma(this->N_);
    scalarListList e_k0(this->N_);
    forAll(X, i)
    {
        kij[i].resize(this->N_);
        mwij[i].resize(this->N_);
        w[i].resize(this->N_);
        sigma[i].resize(this->N_);
        e_k0[i].resize(this->N_);
    }

    scalarList k(this->N_, Zero); //need to be sure
    scalarList e_k(this->N_, Zero);
    scalarList mw_sp3(this->N_, Zero);
    scalarList Vc_sp3(this->N_, Zero);

    forAll(X, spid)
    {
        e_k[spid] = (*this)[spid].Tc_/ 1.2593;
        mw_sp3[spid] = (*this)[spid].W();            //g/mol
        Vc_sp3[spid] = (*this)[spid].Vc_ * 1.0e+03; //Vc is m3/kmol, Vc_sp3 is cm3/mol
    }
    forAll(X, spid)
    {
        forAll(X, spjd)
        {
            w[spid][spjd] = 0.50 * ((*this)[spid].omega_ + (*this)[spjd].omega_;
            sigma[spid][spjd] = sqrt((0.8090 * pow(Vc_sp3[spid], 1.0 / 3.0)) * (0.8090 * pow(Vc_sp3[spjd], 1.0 / 3.0)));
            //sigma(i,j)= bip_kij(i,j)*sqrt((0.809d0*vc(i)**(1.d0/3.d0))*(0.809d0*vc(j)**(1.d0/3.d0)))
            e_k0[spid][spjd] = sqrt(e_k[spid] * e_k[spjd]);
            //e_k0(i,j)= bip_kij(i,j)*sqrt(e_k(i)*e_k(j))
            //kij[spid][spjd] = sqrt(k[spid] * k[spjd]);
            kij[spid][spjd] = sqrt((*this)[spid].kappa_ * (*this)[spid].kappa_;
            mwij[spid][spjd] = 2.0 * mw_sp3[spid] * mw_sp3[spjd] / (mw_sp3[spid] + mw_sp3[spjd]);
        }
    }

    scalar sigma_mix3 = 0.0;
    scalar sigma_mix2 = 0.0;

    //scalar sigma_mix_23 = sigma_mix2/sigma_mix3;
    scalar e_k_mix = 0.0;
    scalar w_mix = 0.0;
    scalar k_mix = 0.0;
    scalar mu4 = 0.0;

    scalar mw_binary = 0.0;
    forAll(X, spid)
    {
        forAll(X, spjd)
        {
            sigma_mix3 += X[spid] * X[spjd] * pow3(sigma[spid][spjd]);
            sigma_mix2 += X[spid] * X[spjd] * sqr(sigma[spid][spjd]);
            e_k_mix += X[spid] * X[spjd] * e_k0[spid][spjd] * pow3(sigma[spid][spjd]);
        }
    }

    forAll(X, spid)
    {
        forAll(X, spjd)
        {

            w_mix += X[spid] * X[spjd] * w[spid][spjd] * pow3(sigma[spid][spjd]);
            mu4 += X[spid] * X[spjd] * sqr(single_specie[spid].mu()) * sqr(single_specie[spjd].mu()) / (pow3(sigma[spid][spjd]) * e_k0[spid][spjd]);
            k_mix += X[spid] * X[spjd] * kij[spid][spjd];
            mw_binary += (X[spid] * X[spjd] * e_k0[spid][spjd] * sqr(sigma[spid][spjd]) * sqrt(mwij[spid][spjd])); // (e_k_mix * sigma_mix2);
        }
    }
    e_k_mix /= sigma_mix3;
    w_mix /= sigma_mix3;
    mu4 *= sigma_mix3 * e_k_mix;
    mw_binary /= (e_k_mix * sigma_mix2);
    mw_binary = sqr(mw_binary);

    scalar volu_mix = sigma_mix3 / pow3(0.8090); //Eq.(9-5.43);Eq.(9-4.8) bug??
    //scalar volu_mix= pow3(sigma_mix3/0.8090);
    scalar tempc_mix = 1.25930 * e_k_mix;
    scalar mur_mix4 = pow4(131.30) * mu4 / sqr(tempc_mix * volu_mix); //volu_mix in cm3/mol, mu4 in debyes
    scalar t_mix = T / e_k_mix;
    scalar Fcm = 1.0 - 0.2756 * w_mix + 0.0590350 * mur_mix4 + k_mix;
    scalar omegav = 1.161450 * (pow(t_mix, -0.148740) + 0.524870 * exp(-0.773200 * t_mix) + 2.161780 * exp(-2.437870 * t_mix)) + (-6.435e-4) * pow(t_mix, 0.148740) * sin(18.0323 * pow(t_mix, -0.7683) - 7.27371); //???GT...

    //1. mixture phases viscosity 40.785 ? 26.69 there is error in 9-5.24
    return 1e-07 * 26.69 * Fcm * sqrt(mw_binary * T) / (sigma_mix2 * omegav);
}
*/



// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
