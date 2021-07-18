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

template <class ThermoMixture>
inline Foam::chungTransportMixture<ThermoMixture>::chungTransportMixture(
    const ThermoMixture &t)
    : ThermoMixture(t)
{
}

template <class ThermoMixture>
inline Foam::chungTransportMixture<ThermoMixture>::chungTransportMixture(
    const word &name,
    const chungTransportMixture &ct)
    : ThermoMixture(name, ct)
{
}

template <class ThermoMixture>
inline Foam::autoPtr<Foam::chungTransportMixture<ThermoMixture>>
Foam::chungTransportMixture<ThermoMixture>::clone() const
{
    return autoPtr<chungTransportMixture<ThermoMixture>>(
        new chungTransportMixture<ThermoMixture>(*this));
}

template <class ThermoMixture>
inline Foam::autoPtr<Foam::chungTransportMixture<ThermoMixture>>
Foam::chungTransportMixture<ThermoMixture>::New(
    const dictionary &dict)
{
    return autoPtr<chungTransportMixture<ThermoMixture>>(
        new chungTransportMixture<ThermoMixture>(dict));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class ThermoMixture>
inline Foam::scalar Foam::chungTransportMixture<ThermoMixture>::mu(
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
        e_k[spid] = (*this)[spid].Tc_ / 1.2593;
        mw_sp3[spid] = (*this)[spid].W();           //g/mol
        Vc_sp3[spid] = (*this)[spid].Vc_ * 1.0e+03; //Vc is m3/kmol, Vc_sp3 is cm3/mol
    }
    forAll(X, spid)
    {
        forAll(X, spjd)
        {
            w[spid][spjd] = 0.50 * ((*this)[spid].omega_ + (*this)[spjd].omega_);
            sigma[spid][spjd] = sqrt((0.8090 * pow(Vc_sp3[spid], 1.0 / 3.0)) * (0.8090 * pow(Vc_sp3[spjd], 1.0 / 3.0)));
            //sigma(i,j)= bip_kij(i,j)*sqrt((0.809d0*vc(i)**(1.d0/3.d0))*(0.809d0*vc(j)**(1.d0/3.d0)))
            e_k0[spid][spjd] = sqrt(e_k[spid] * e_k[spjd]);
            //e_k0(i,j)= bip_kij(i,j)*sqrt(e_k(i)*e_k(j))
            //kij[spid][spjd] = sqrt(k[spid] * k[spjd]);
            kij[spid][spjd] = sqrt((*this)[spid].kappa_ * (*this)[spid].kappa_);
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
            mu4 += X[spid] * X[spjd] * sqr((*this)[spid].mu_) * sqr((*this)[spjd].mu_) / (pow3(sigma[spid][spjd]) * e_k0[spid][spjd]);
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

template <class Specie>
inline Foam::scalar Foam::chungTransportMixture<Specie>::kappa(
    scalar p,
    scalar T,
    const label pf,
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
        e_k[spid] = (*this)[spid].Tc_ / 1.2593;
        mw_sp3[spid] = (*this)[spid].W();           //g/mol
        Vc_sp3[spid] = (*this)[spid].Vc_ * 1.0e+03; //Vc is m3/kmol, Vc_sp3 is cm3/mol
    }
    forAll(X, spid)
    {
        forAll(X, spjd)
        {
            w[spid][spjd] = 0.50 * ((*this)[spid].omega_ + (*this)[spjd].omega_);
            sigma[spid][spjd] = sqrt((0.8090 * pow(Vc_sp3[spid], 1.0 / 3.0)) * (0.8090 * pow(Vc_sp3[spjd], 1.0 / 3.0)));
            //sigma(i,j)= bip_kij(i,j)*sqrt((0.809d0*vc(i)**(1.d0/3.d0))*(0.809d0*vc(j)**(1.d0/3.d0)))
            e_k0[spid][spjd] = sqrt(e_k[spid] * e_k[spjd]);
            //e_k0(i,j)= bip_kij(i,j)*sqrt(e_k(i)*e_k(j))
            //kij[spid][spjd] = sqrt(k[spid] * k[spjd]);
            kij[spid][spjd] = sqrt((*this)[spid].kappa_ * (*this)[spid].kappa_);
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
            mu4 += X[spid] * X[spjd] * sqr((*this)[spid].mu_) * sqr((*this)[spjd].mu_) / (pow3(sigma[spid][spjd]) * e_k0[spid][spjd]);
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
    //scalar mur_mix4 = pow4(131.30) * mu4 / sqr(tempc_mix * volu_mix); //volu_mix in cm3/mol, mu4 in debyes
    //scalar t_mix = T / e_k_mix;
    //scalar omegav = 1.161450 * (pow(t_mix, -0.148740) + 0.524870 * exp(-0.773200 * t_mix) + 2.161780 * exp(-2.437870 * t_mix)) + (-6.435e-4) * pow(t_mix, 0.148740) * sin(18.0323 * pow(t_mix, -0.7683) - 7.27371); //???GT...

    //1. mixture phases viscosity 40.785 ? 26.69 there is error in 9-5.24
    //viscos = 1e-07 * 26.69 * Fcm * sqrt(mw_binary * T) / (sigma_mix2 * omegav);

    scalar mw_phase = this->W(Xinp);
    scalar Cv_phase = this->Cp(p, T, pf, Xinp) - this->CpMCv(p, T, pf, Xinp); //J/kgK
    scalar alpha = Cv_phase * mw_phase / (RR * 1.0e-03) - 1.50;
    scalar tempr = T / tempc_mix;

    //compute specific volume, m3/mol
    //scalar amix_phase = this->a(p,T,Xinp), bmix_phase = this->b(p,T,Xinp), dadt_phase = this->dadt(p,T,Xinp), dadt2 = this->d2adt2(p,T,Xinp);
    scalar sv_phase = this->V(p, T, pf, Xinp); //m3/mol

    //2, heat conductivity,W/m.K
    scalar beta = 0.78620 - 0.71090 * w_mix + 1.31680 * sqr(w_mix);
    scalar mw_binary0 = mw_binary / 1.0e+3;
    scalar zeta = 2.00 + 10.50 * sqr(tempr);
    scalar psi = 1.0 + alpha * ((0.2150 + 0.282880 * alpha - 1.0610 * beta + 0.266650 * zeta) / (0.63660 + beta * zeta + 1.0610 * alpha * beta));

    scalar qq = 3.5860 * 1e-03 * sqrt(tempc_mix / mw_binary0) / pow(volu_mix, 2.0 / 3.0);
    scalar yy = volu_mix / (6.0e+06 * sv_phase);
    scalar G1 = (1.0 - 0.50 * yy) / pow3(1.0 - yy);

    // Table 10-3
    scalar a[7] = {2.41660, -5.0924e-01, 6.61070, 1.4543e+01, 7.9274e-01, -5.86340, 9.1098e+01};
    scalar bb[7] = {7.4824e-01, -1.50940, 5.62070, -8.91390, 8.2019e-01, 1.2801e+01, 1.2811e+02};
    scalar c[7] = {-9.1858e-01, -4.9991e+01, 6.4760e+01, -5.63790, -6.9369e-01, 9.58930, -5.4217e+01};
    scalar d[7] = {1.2712e+02, 6.9983e+01, 2.7039e+01, 7.4344e+01, 6.31730, 6.5529e+01, 5.2381e+02};

    scalar B[8] = {0.0};
    for (label ii = 0; ii < 7; ii++)
    {
        B[ii] = a[ii] + bb[ii] * w_mix + c[ii] * mu4 + d[ii] * k_mix;
    }

    scalar G2 = ((B[1] / yy) * (1.0 - exp(-B[4] * yy)) + B[2] * G1 * exp(B[5] * yy) + B[3] * G1) / (B[1] * B[4] + B[2] + B[3]);
    //Info<< "G2* * * *"<<G2<<"* * "<< yy <<"* * "<< volu_mix <<"**"<<sv_phase << endl;
    //Info<< "comp_phase* * * *"<<comp_phase[0] <<"* * "<< comp_phase[1] << endl;

    return 31.20 * mu(p, T, Xinp) * psi * (1.0 / G2 + B[6] * yy) / mw_binary0 + qq * B[7] * sqr(yy) * sqrt(tempr) * G2;
    //printf("kappa_phase_thcond_phase %5.8f\n\n",thcond_phase);
}

template <class ThermoMixture>
inline Foam::scalar Foam::chungTransportMixture<ThermoMixture>::alphah(
    const scalar p,
    const scalar T) const
{
    return kappa(p, T,1) / this->Cp(p, T,1);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //
template <class Specie>
Foam::chungTransportMixture<Specie> Foam::operator==(
    const chungTransportMixture<Specie> &a,
    const chungTransportMixture<Specie> &b)
{
    return chungTransportMixture<Specie>(a);
}

template <class Specie>
Foam::chungTransportMixture<Specie> Foam::operator*(
    const scalar s,
    const chungTransportMixture<Specie> &b)
{
    return chungTransportMixture<Specie>(b);
}

template <class Specie>
inline void Foam::chungTransportMixture<Specie>::operator+=(const chungTransportMixture<Specie> &pg)
{
    return;
}

// ************************************************************************* //
