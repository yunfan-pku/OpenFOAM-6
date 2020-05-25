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

#include "PengRobinsonS.H"
#include "specie.H"
#include "mathematicalConstants.H"
#include "OFstream.H"
#include "IFstream.H"
#include "string.H"
#define PRMIX

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class Specie>
inline Foam::PengRobinsonS<Specie>::PengRobinsonS(
    const word &name,
    const PengRobinsonS &pg)
    : Specie(name, pg),
      Tc_(pg.Tc_),
      Vc_(pg.Vc_),
      Zc_(pg.Zc_),
      Pc_(pg.Pc_),
      omega_(pg.omega_)
{
}

template <class Specie>
inline Foam::autoPtr<Foam::PengRobinsonS<Specie>>
Foam::PengRobinsonS<Specie>::clone() const
{
    return autoPtr<PengRobinsonS<Specie>>(
        new PengRobinsonS<Specie>(*this));
}

template <class Specie>
inline Foam::autoPtr<Foam::PengRobinsonS<Specie>>
Foam::PengRobinsonS<Specie>::New(
    const dictionary &dict)
{
    return autoPtr<PengRobinsonS<Specie>>(
        new PengRobinsonS<Specie>(dict));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class Specie>
inline Foam::scalar Foam::PengRobinsonS<Specie>::Tc() const
{
    return Tc_;
}
template <class Specie>
inline Foam::scalar Foam::PengRobinsonS<Specie>::Pc() const
{
    return Pc_;
}
template <class Specie>
inline Foam::scalar Foam::PengRobinsonS<Specie>::Vc() const
{
    return Vc_;
}
template <class Specie>
inline Foam::scalar Foam::PengRobinsonS<Specie>::omega() const
{
    return omega_;
}
template <class Specie>
inline Foam::scalar Foam::PengRobinsonS<Specie>::Zc() const
{
    return Zc_;
}

template <class Specie>
inline void Foam::PengRobinsonS<Specie>::AzBz(
    scalar Pr,
    scalar Tr,
    scalar &Az,
    scalar &Bz) const
{
    scalar m_omega;
    scalar alpha;
    if (omega_ > 0.49)
    {
        m_omega = omega_ * (omega_ * (omega_ * (0.016666) - 0.164423) + 1.48503) + 0.379642;
    }
    else
    {
        m_omega = omega_ * (omega_ * (-0.26992) + 1.54226) + 0.37464;
    }
    alpha = sqr(1.0 + m_omega * (1.0 - sqrt(Tr)));
    Az = 0.45724 * alpha * Pr / sqr(Tr);
    Bz = 0.07780 * Pr / Tr;
    return;
}

template <class Specie>
inline void Foam::PengRobinsonS<Specie>::azbz(
    scalar Tr,
    scalar &m_omega,
    scalar &az,
    scalar &bz,
    scalar &alpha,
    scalar &comp4
    ) const
{
    if (omega_ > 0.49)
    {
        m_omega = omega_ * (omega_ * (omega_ * (0.016666) - 0.164423) + 1.48503) + 0.379642;
    }
    else
    {
        m_omega = omega_ * (omega_ * (-0.26992) + 1.54226) + 0.37464;
    }
    alpha = sqr(1.0 + m_omega * (1.0 - sqrt(Tr)));
    az = 0.45724 * alpha * sqr(RR * 1.0e-03 * Tc_) / Pc_;
    bz = 0.07780 * RR * 1.0e-03 * Tc_ / Pc_;
    comp4 = m_omega/ sqrt(Tc_);
    return;
}

// ======================================================================//
// Ideal phase mixture enthalpy						 //
// Input: temperature, phase compositions				 //
// Output: ideal phase mixture H			 		 //
//									 //
// Source: F.A. Aly, L.L. Lee, 						 //
// Self-consistent equations for calculating the ideal gas heat capacity&//
// Fluid Phase Equilibria 6 (3-4) (1981) 169–179.			 //
// Coefficients are in NASA Technical Momorandum 4513: 			 //
// Calculating Thermodynamic and Transport Properties			 //
// Reference temperature is 0 K, tref = 0.d0				 //
//									 //
// History :   Creation.	     P. YI			2019     //
// ======================================================================//
template<class Specie>
inline Foam::scalar Foam::PengRobinsonS<Specie>::Hig_phase 
(
    scalar T
) const
{
    return RR*1.0e-03*T*(T*(T*(T*(T*Hig_phase_[4]/5.0+Hig_phase_[3]/4.0)+Hig_phase_[2]/3.0)+Hig_phase_[1]/2.0)+Hig_phase_[0]);
}


// ======================================================================//
// Ideal phase mixture enthalpy						 //
// Input: temperature, phase compositions				 //
// Output: ideal phase mixture H			 		 //
//									 //
// Source: F.A. Aly, L.L. Lee, 						 //
// Self-consistent equations for calculating the ideal gas heat capacity&//
// Fluid Phase Equilibria 6 (3-4) (1981) 169–179.			 //
// Coefficients are in above F.A. Aly paper:				 //
// Calculating Thermodynamic and Transport Properties			 //
// Reference temperature is 0 K, tref = 0.d0				 //
//									 //
// History :   Creation.	     P. YI			2019     //
// ======================================================================//
template<class Specie>
inline Foam::scalar Foam::PengRobinsonS<Specie>::Hig2_phase 
(
    scalar T
) const
{
	return   Hig2_phase_[0]*T +Hig2_phase_[1]*Hig2_phase_[2]*(2.0/(exp(2.0*Hig2_phase_[2]/T)-1.0))+ 2.0*Hig2_phase_[3]*Hig2_phase_[4]*(1.0/(exp(2.0*Hig2_phase_[4]/T)+1.0));
}


// ======================================================================//
// Ideal specific heat capacity						 //
// Input: temperature, phase compositions				 //
// Output: ideal phase mixture H			 		 //
//									 //
// Source: F.A. Aly, L.L. Lee, 						 //
// Self-consistent equations for calculating the ideal gas heat capacity&//
// Fluid Phase Equilibria 6 (3-4) (1981) 169–179.			 //
// Coefficients are in NASA Technical Momorandum 4513: 			 //
// Calculating Thermodynamic and Transport Properties			 //
// Reference temperature is 0 K, tref = 0.d0				 //
//									 //
// History :   Creation.	     P. YI			2019     //
// ======================================================================//
template<class Specie>
inline Foam::scalar Foam::PengRobinsonS<Specie>::cp_phase
(
    scalar T
) const
{   
    return   RR*1.0e-03*(T*(T*(T*(T*Hig_phase_[4]+Hig_phase_[3])+Hig_phase_[2]) +Hig_phase_[1])+Hig_phase_[0]);
}

template<class Specie>
inline Foam::scalar Foam::PengRobinsonS<Specie>::mu() const
{
    return mu_;
}

template<class Specie>
inline Foam::scalar Foam::PengRobinsonS<Specie>::kappa() const
{
    return kappa_;
}

template <class Specie>
inline Foam::scalar Foam::PengRobinsonS<Specie>::rho(
    scalar p,
    scalar T) const
{
    const scalar Z = this->Z(p, T);
    return p / (Z * this->R() * T);
}

template <class Specie>
inline Foam::scalar Foam::PengRobinsonS<Specie>::H(scalar p, scalar T) const
{
    const scalar Pr = p / Pc_;
    const scalar Tr = T / Tc_;
    const scalar B = 0.07780 * Pr / Tr;
    const scalar kappa = 0.37464 + 1.54226 * omega_ - 0.26992 * sqr(omega_);
    const scalar alpha = sqr(1 + kappa * (1 - sqrt(Tr)));

    const scalar Z = this->Z(p, T);

    return this->R() * Tc_ * (Tr * (Z - 1) - 2.078 * (1 + kappa) * sqrt(alpha) * log((Z + 2.414 * B) / (Z - 0.414 * B)));
}

template <class Specie>
inline Foam::scalar Foam::PengRobinsonS<Specie>::Cp(scalar p, scalar T) const
{
    const scalar Tr = T / Tc_;
    const scalar a = 0.45724 * sqr(RR * Tc_) / Pc_;
    const scalar b = 0.07780 * RR * Tc_ / Pc_;
    const scalar kappa = 0.37464 + 1.54226 * omega_ - 0.26992 * sqr(omega_);
    const scalar alpha = sqr(1 + kappa * (1 - sqrt(Tr)));

    const scalar A = a * alpha * p / sqr(RR * T);
    const scalar B = b * p / (RR * T);

    const scalar Z = this->Z(p, T);

    const scalar ap = kappa * a * (kappa / Tc_ - (1 + kappa) / sqrt(T * Tc_));
    const scalar app = kappa * a * (1 + kappa) / (2 * sqrt(pow3(T) * Tc_));

    const scalar M = (sqr(Z) + 2 * B * Z - sqr(B)) / (Z - B);
    const scalar N = ap * B / (b * RR);

    const scalar root2 = sqrt(2.0);

    return (
               app * (T / (2 * root2 * b)) * log((Z + (root2 + 1) * B) / (Z - (root2 - 1) * B)) + RR * sqr(M - N) / (sqr(M) - 2 * A * (Z + B)) - RR) /
           this->W();
}

template <class Specie>
inline Foam::scalar Foam::PengRobinsonS<Specie>::S(
    scalar p,
    scalar T) const
{
    const scalar Pr = p / Pc_;
    const scalar Tr = T / Tc_;
    const scalar B = 0.07780 * Pr / Tr;
    const scalar kappa = 0.37464 + 1.54226 * omega_ - 0.26992 * sqr(omega_);

    const scalar Z = this->Z(p, T);

    return this->R() * (-log(p / Pstd) + (log(Z - B) - 2.078 * kappa * ((1 + kappa) / sqrt(Tr) - kappa) * log((Z + 2.414 * B) / (Z - 0.414 * B))));
}

template <class Specie>
inline Foam::scalar Foam::PengRobinsonS<Specie>::psi(
    scalar p,
    scalar T) const
{
    const scalar Z = this->Z(p, T);

    return 1.0 / (Z * this->R() * T);
}

template <class Specie>
inline Foam::scalar Foam::PengRobinsonS<Specie>::Z(
    scalar p,
    scalar T) const
{
    const scalar Tr = T / Tc_;
    const scalar a = 0.45724 * sqr(RR * Tc_) / Pc_;
    const scalar b = 0.07780 * RR * Tc_ / Pc_;
    const scalar kappa = 0.37464 + 1.54226 * omega_ - 0.26992 * sqr(omega_);
    const scalar alpha = sqr(1 + kappa * (1 - sqrt(Tr)));

    const scalar A = a * alpha * p / sqr(RR * T);
    const scalar B = b * p / (RR * T);

    const scalar a2 = B - 1;
    const scalar a1 = A - 2 * B - 3 * sqr(B);
    const scalar a0 = -A * B + sqr(B) + pow3(B);

    const scalar Q = (3 * a1 - a2 * a2) / 9.0;
    const scalar Rl = (9 * a2 * a1 - 27 * a0 - 2 * a2 * a2 * a2) / 54;

    const scalar Q3 = Q * Q * Q;
    const scalar D = Q3 + Rl * Rl;

    scalar root = -1;

    if (D <= 0)
    {
        const scalar th = ::acos(Rl / sqrt(-Q3));
        const scalar qm = 2 * sqrt(-Q);
        const scalar r1 = qm * cos(th / 3) - a2 / 3;
        const scalar r2 =
            qm * cos((th + 2 * constant::mathematical::pi) / 3) - a2 / 3;
        const scalar r3 =
            qm * cos((th + 4 * constant::mathematical::pi) / 3) - a2 / 3;

        root = max(r1, max(r2, r3));
    }
    else
    {
        // One root is real
        const scalar D05 = sqrt(D);
        const scalar S = pow(Rl + D05, 1.0 / 3.0);
        scalar Tl = 0;
        if (D05 > Rl)
        {
            Tl = -pow(mag(Rl - D05), 1.0 / 3.0);
        }
        else
        {
            Tl = pow(Rl - D05, 1.0 / 3.0);
        }

        root = S + Tl - a2 / 3;
    }

    return root;
}

double phasestate(double vaporfra)
{
    if (vaporfra > 1.0e-06 && vaporfra < 0.999999)
    {
        return 2.0;
    }
    else if (vaporfra <= 1.0e-06)
    {
        return 1.0;
    }
    else if (vaporfra >= 0.999999)
    {
        return 0.0;
    }
    return 0;
}

// ======================================================================//
// Evaluate the phase volume fraction, alpha				 //
// Input: pressure, temperature, vapor fraction				 //
// Output: phase compositions 				 		 //
//									 //
// History :   Creation.	     P. YI			2019     //
// ======================================================================//
template <class Specie>
inline Foam::scalar Foam::PengRobinsonS<Specie>::Evaluate_alpha(
    scalar press,     //inp
    scalar temp,      //inp
    scalar vaporfra,  //inp
    scalar *comp_liq, //inp
    scalar *comp_gas, //inp
    scalar *comp_inp) const
{
    scalar sv_gas = this->volmmix_phase(0, press, temp, comp_gas); //gas phase specific volume
    scalar sv_liq = this->volmmix_phase(1, press, temp, comp_liq); //liquid phase specific volume
                                                                   /*
    scalar mw_mixture = this->mwmix(comp_inp);
    scalar mw_gas = this->mwmix(comp_gas);
    scalar msgas = mw_gas/sv_gas;
    scalar ygas = vaporfra*mw_gas/mw_mixture;
    scalar mw_liq = this->mwmix(comp_liq);
    scalar msliq = mw_liq/sv_liq;
    scalar alphagas = ygas*msliq/(msgas-ygas*(msgas-msliq));
*/
    scalar alphagas = vaporfra * sv_gas / (vaporfra * sv_gas + (1.0 - vaporfra) * sv_liq);
    if (alphagas >= 0.9999999)
    {
        alphagas = 0.9999999;
    }
    else if (alphagas <= 1.0e-07)
    {
        alphagas = 1.0e-07;
    }

    return alphagas;
}

// ======================================================================//
// Ideal phase mixture internal energy => U = H - PV			 //
// Input: temperature, phase compositions				 //
// Output: ideal phase mixture sie			 		 //
//									 //
// History :   Creation.	     P. YI			2019     //
// ======================================================================//
template <class Specie>
inline Foam::scalar Foam::PengRobinsonS<Specie>::sieig_phase(
    scalar T,
    scalar *comp_phase //note, it is phase compositions
    ) const
{
    scalar mw_phase = this->mwmix(comp_phase);
    int sieflag = 2;
    scalar sie_igkg = 0.0;
    if (sieflag == 1)
    {
        sie_igkg = this->Hig_phase(T, comp_phase);
    }
    else if (sieflag == 2)
    {
        sie_igkg = this->Hig2_phase(T, comp_phase);
    }

    sie_igkg = sie_igkg - RR * 1.0e-03 * T;

    return sie_igkg / mw_phase;
}

//can be used for phases or mixtures, depending on the Azin and Bzin, three roots
template <class Specie>
inline void Foam::PengRobinsonS<Specie>::Root_Z(
    int phasei, //input
    scalar Azin,
    scalar Bzin,
    scalar *rootz) const
{
    const scalar A = Azin;
    const scalar B = Bzin;

    const scalar a2 = B - 1.0;
    const scalar a1 = A - 2.0 * B - 3.0 * sqr(B);
    const scalar a0 = -A * B + sqr(B) + pow3(B);

    const scalar Q = (3.0 * a1 - a2 * a2) / 9.0;
    const scalar Rl = (9.0 * a2 * a1 - 27.0 * a0 - 2.0 * a2 * a2 * a2) / 54.0;

    const scalar Q3 = Q * Q * Q;
    const scalar D = Q3 + Rl * Rl;
    //printf("Zmix_root00= %+5.10f,%+5.8f,%+5.8f\n", D,Q3,Rl*Rl);
    scalar root1 = -1, root2 = -1, root3 = -1;
    if (D <= 0)
    {
        const scalar th = ::acos(Rl / sqrt(-Q3));
        const scalar qm = 2 * sqrt(-Q);
        root1 = qm * cos(th / 3) - a2 / 3;
        root2 = qm * cos((th + 2 * constant::mathematical::pi) / 3) - a2 / 3;
        root3 = qm * cos((th + 4 * constant::mathematical::pi) / 3) - a2 / 3;
        rootz[0] = root1;
        rootz[1] = root2;
        rootz[2] = root3;
    }
    else
    {
        // One root is real
        const scalar D05 = sqrt(D);
        const scalar RlD = Rl + D05;
        scalar S = 0.0;
        if (RlD > 0.0)
        {
            S = pow(Rl + D05, 1.0 / 3.0);
        }
        else
        {
            S = -pow(-RlD, 1.0 / 3.0);
        }
        scalar Tl = 0;
        //printf("Zmix_root200= %+5.8f,%+5.8f,%+5.8f\n", S,Rl,D05);
        if (D05 > Rl)
        {
            Tl = -pow(mag(Rl - D05), 1.0 / 3.0);
        }
        else
        {
            Tl = pow(Rl - D05, 1.0 / 3.0);
        }
        scalar root1 = S + Tl - a2 / 3;
        rootz[0] = root1;
        rootz[1] = root1;
        rootz[2] = root1;
    }
}

template <class Specie>
inline Foam::scalar Foam::PengRobinsonS<Specie>::volmmix_phase //m3/mol
    (
        int phasei, //input
        scalar p,
        scalar T,
        scalar *comp_inp) const
{
    scalar volspecific = 0.0;
    scalar mw_phase = this->mwmix(comp_inp);                             //kg/mol
    const scalar rho_phase = this->rhomix_phase(phasei, p, T, comp_inp); //kg/m3
    volspecific = mw_phase / rho_phase;

    /*
    scalar root1 = 0.0, root2 = 0.0, root3 = 0.0;
    scalar amixvol = 0.0, bmixvol = 0.0;
    scalar rootzvol[3] = {0.0};
    scalar dadtvol = 0.0, dadt2vol = 0.0;

    this->abmix(T, comp_inp, amixvol, bmixvol, dadtvol, dadt2vol);  

    if (phasei == 0) 
    {
	this->Root_Z(0,amixvol,bmixvol,rootzvol); //three roots
	root1 = rootzvol[0];
	root2 = rootzvol[1];
	root3 = rootzvol[2];
        volspecific = max(root1, max(root2, root3));
    }
    else
    {
	this->Root_Z(1,amixvol,bmixvol,rootzvol); //three roots
	root1 = rootzvol[0];
	root2 = rootzvol[1];
	root3 = rootzvol[2];
        volspecific = min(root1, min(root2, root3));
    }
*/
    return volspecific;
}

template <class Specie>
inline void Foam::PengRobinsonS<Specie>::ab_spd(
    scalar p,       //input
    scalar T,       //input
    scalar Tcsp,    //input
    scalar Pcsp,    //input
    scalar omegasp, //input

    scalar &alphaz,
    scalar &az,
    scalar &bz,
    scalar &dadtz) const
{

    scalar Trsp = T / Tcsp;
    scalar m_omega_z = 0.0;

    if (omegasp > 0.49)
    {
        m_omega_z = 0.379642 + 1.48503 * omegasp - 0.164423 * sqr(omegasp) + 0.016666 * pow3(omegasp);
    }
    else
    {
        m_omega_z = 0.37464 + 1.54226 * omegasp - 0.26992 * sqr(omegasp);
    }

    alphaz = sqr(1.0 + omegasp * (1.0 - sqrt(Trsp)));
    az = 0.45724 * alphaz * sqr(RR * 1.0e-03 * Tcsp) / Pcsp;
    bz = 0.07780 * RR * 1.0e-03 * Tcsp / Pcsp;
    dadtz = (2.0 * 0.45725 * sqr(RR * 1.0e-03 * Tcsp) / Pcsp) * sqrt(alphaz) * (-m_omega_z * 1.0 / (2.0 * sqrt(T * Tcsp)));
}

template <class Specie>
inline Foam::scalar Foam::PengRobinsonS<Specie>::siemix_phase // internal energy with unit of J/kg
    (
        int phasei, //input
        scalar p,
        scalar T,
        scalar *comp_phase) const
{
    scalar amix = 0.0;
    scalar bmix = 0.0;
    scalar dadt = 0.0;
    scalar dadt2 = 0.0;
    scalar mw_phase = this->mwmix(comp_phase);                       //kg*mol-1
    scalar rho_phase = this->rhomix_phase(phasei, p, T, comp_phase); //kg*m-3

    const scalar cc1 = 1.0 + sqrt(2.0);
    const scalar cc2 = 1.0 - sqrt(2.0);

    scalar VV = mw_phase / rho_phase; //m3*mol-1, phase specific volume

    this->abmix(T, comp_phase, amix, bmix, dadt, dadt2);
    scalar xxx = ((VV + cc1 * bmix) / (VV + cc2 * bmix));
    scalar dsie_phase = 0.0;
    if (xxx > 0)
    {
        dsie_phase = ((T * dadt - amix) / (2.0 * sqrt(2.0) * bmix)) * log((VV + cc1 * bmix) / (VV + cc2 * bmix));
    }
    else
    {
        dsie_phase = 0.0;
    }

    dsie_phase = dsie_phase / mw_phase; // J/kg
    // Based on Eq.12, Internal energy is estimated based on ideal(close to absolute 0 T) + departure
    scalar sieReal_phase = 0.0;
    scalar sieig = this->sieig_phase(T, comp_phase);
    sieReal_phase = dsie_phase + sieig;

    return sieReal_phase;
}

template <class Specie>
inline Foam::scalar Foam::PengRobinsonS<Specie>::cvmix_phase // PY:phase mixture heat capability
    (
        int phasei, //input
        scalar p,
        scalar T,
        scalar *comp_phase) const
{
    scalar cv_phase;
    int cvcpflag = 1;

    if (cvcpflag == 0) //ideal
    {
        scalar cpig_phase = 0.0, cvig_phase = 0.0;
        this->cvcpig_phase(T, comp_phase, cpig_phase, cvig_phase);
        cv_phase = cvig_phase;
    }
    else if (cvcpflag == 1) //real 1
    {
        scalar sie_phase = this->siemix_phase(phasei, p, T, comp_phase);
        scalar deltat = 10.0;
        scalar T_plus = T + deltat;
        scalar sie_phase_plus = this->siemix_phase(phasei, p, T_plus, comp_phase);

        scalar delta_sie = sie_phase_plus - sie_phase;
        cv_phase = delta_sie / deltat;
    }
    return cv_phase;
}

template <class Specie>
inline Foam::scalar Foam::PengRobinsonS<Specie>::cpmix_phase // PY:phase mixture heat capability
    (
        int phasei, //input
        scalar p,
        scalar T,
        scalar *comp_phase) const
{
    scalar cv_phase = this->cvmix_phase(phasei, p, T, comp_phase);

    scalar amix = 0.0;
    scalar bmix = 0.0;
    scalar dadt = 0.0;
    scalar dadt2 = 0.0;

    this->abmix(T, comp_phase, amix, bmix, dadt, dadt2);
    scalar mw_phase = this->mwmix(comp_phase);
    scalar rho_phase = this->rhomix_phase(phasei, p, T, comp_phase);
    scalar VV = mw_phase / rho_phase;
    scalar den = sqr(VV) + 2.0 * bmix * VV - sqr(bmix);
    scalar dpdt = RR * 1.0e-03 / (VV - bmix) - dadt / den;
    scalar dpdv = -RR * 1.0e-03 * T / sqr(VV - bmix) + 2.0 * amix * (VV + bmix) / sqr(den);
    scalar ct2 = -dpdv * mw_phase / sqr(rho_phase);
    scalar cp_phase = cv_phase + T * sqr(dpdt) / (ct2 * sqr(rho_phase));
    //printf("VV,cp_phase= %+5.8f %+5.8f %+5.8f %+5.8f\n",VV,rho_phase,cv_phase,cp_phase);
    //scalar cp_phase = cv_phase;

    if (cp_phase <= 100)
    {
        cp_phase = 100;
    }

    return cp_phase;
}

template <class Specie>
inline void Foam::PengRobinsonS<Specie>::soundspeedmix // m/s
    (
        int phasei, //input
        scalar p,
        scalar T,
        scalar *comp_inp,
        scalar &ct_phase, //output
        scalar &cs_phase) const
{
    scalar cv_phase = this->cvmix_phase(phasei, p, T, comp_inp);
    scalar mw_phase = this->mwmix(comp_inp);
    scalar rho_phase = this->rhomix_phase(phasei, p, T, comp_inp);
    scalar VV = mw_phase / rho_phase;
    scalar amix = 0.0, bmix = 0.0;
    scalar dadt = 0.0, dadt2 = 0.0;
    this->abmix(T, comp_inp, amix, bmix, dadt, dadt2);
    scalar den = sqr(VV) + 2.0 * bmix * VV - sqr(bmix);
    scalar dpdt = RR * 1.0e-03 / (VV - bmix) - dadt / den;
    scalar dpdv = -RR * 1.0e-03 * T / sqr(VV - bmix) + 2.0 * amix * (VV + bmix) / sqr(den);
    scalar ct2 = -dpdv * mw_phase / sqr(rho_phase);

    // PY: Deal with the negative sound speed
    if (ct2 < 0.0)
    {
        ct2 = 1.0e-04;
    }
    scalar cs2 = ct2 + T * sqr(dpdt) / (cv_phase * sqr(rho_phase));

    cs_phase = sqrt(cs2);
    ct_phase = sqrt(ct2);
}

// ======================================================================//
// Correction for mass diffusion in HIGH pressure                        //
// based on Ping's critical point	                                 //
// Harstad and Bellan Ref. Ind. Eng. Chem. Res. 2004, 43, 645-654        //
// Balenovic et al., J. Chem. Phys., 1970                                //
// Input: feed, p, T            					 //
// Output: Correction for Dij under high pressure conditions		 //
// Range: 27.6 - 132 MPa                                                 //
//									 //
// History :   Creation.	     P. YI		2019.07		 //
// ======================================================================//

template <class Specie>
inline Foam::scalar Foam::PengRobinsonS<Specie>::Dij_correct //for gas
    (
        scalar p,
        scalar T,
        scalar *comp_inp) const
{
    // Tcij and Vcmix are interpolated;
    scalar Tcij = this->TcBi(comp_inp);
    scalar Pcij = this->PcBi(comp_inp);
    scalar Vcij = 1.0e+06 * this->volmmix_phase(0, Pcij, Tcij, comp_inp);
    scalar Vmol = 1.0e+06 * this->volmmix_phase(0, p, T, comp_inp);
    // Reduced density
    scalar rho_reduced = Vcij / Vmol;
    // Fit parameters
    scalar c_adjust = 0.58;
    scalar b_adjust = 3.0 / 2;
    // Scaling function
    scalar theta_D = c_adjust * pow(rho_reduced, b_adjust);
    return theta_D;
}

template <class Specie>
inline Foam::scalar Foam::PengRobinsonS<Specie>::Dij_highP //for dense gas
    (
        scalar p,
        scalar T,
        scalar *comp_inp) const
{
    scalar DijL = this->Dij_lowP(p, T, comp_inp);
    scalar WD = this->Dij_correct(p, T, comp_inp) + 1.0;
    return DijL / WD;
}

template <class Specie>
inline Foam::scalar Foam::PengRobinsonS<Specie>::Dij_highP_pesudo //for dense gas
    (
        scalar p,
        scalar T,
        scalar *comp_inp) const
{
    scalar DijL = this->Dij_lowP_pesudo(p, T, comp_inp);
    scalar WD = this->Dij_correct_pesudo(p, T, comp_inp) + 1.0;

    return DijL / WD;
}

// ======================================================================//
// No phase change: From sie to Temp					 //
// Input: feed, p, T, and sie, ro					 //
// Output: relaxed p and T 				 		 //
//									 //
// History :   Creation.	     P. YI		2019.06		 //
// ======================================================================//

template <class Specie>
inline void Foam::PengRobinsonS<Specie>::from_sie_ro_to_PT(
    scalar sie_inp,
    scalar ro_inp,
    scalar *comp_inp, //input:feed
    scalar &press,
    scalar &temp) const
{
    scalar amix = 0.0, bmix = 0.0, dadt = 0.0, dadt2 = 0.0;

    scalar mw_mixture = this->mwmix(comp_inp);
    scalar sv_mixture = mw_mixture / ro_inp;
    int iter = 0;
    int niter_max = 100;
    scalar delta_sie = 1.0e-08;
    while (iter < niter_max)
    {
        this->abmix(temp, comp_inp, amix, bmix, dadt, dadt2);
        scalar del = sqr(sv_mixture) + 2.0 * bmix * sv_mixture - sqr(bmix);
        const scalar cc0 = 2.0 * sqrt(2.0);
        const scalar cc1 = 1.0 + sqrt(2.0);
        const scalar cc2 = 1.0 - sqrt(2.0);

        scalar sie_d = ((temp * dadt - amix) / (cc0 * bmix)) * log((sv_mixture + cc1 * bmix) / (sv_mixture + cc2 * bmix));

        scalar sieig = this->sieig_phase(temp, comp_inp);
        //call sie_ig(temp0,comp,spdfre(i,1:nspfre),spdfue(i,1:nbffuel),sieig,cv_ig,h_ig,cp_ig)
        // change the units
        scalar sie_test = sie_d / mw_mixture + sieig;
        scalar epsilon = sqr((sie_inp - sie_test) / sie_inp);

        if (epsilon < delta_sie)
        {
            press = RR * 1.0e-03 * temp / (sv_mixture - bmix) - amix / del;
            break;
        }
        else
        {
            // update T,P, using NEWton Rapson
            scalar press = RR * 1.0e-03 * temp / (sv_mixture - bmix) - amix / del;
            // compute dudt_v, dudt_v=cv_real
            // need another subrotine to determin Liq or gas;
            scalar cv_real = this->cvmix_phase(0, press, temp, comp_inp);
            // Rewton Rapson
            temp = temp - (sie_inp - sie_test) / (-cv_real);
        }
        iter++;
    }
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template <class Specie>
inline void Foam::PengRobinsonS<Specie>::operator+=(
    const PengRobinsonS<Specie> &pg)
{
    scalar Y1 = this->Y();  //YP: Y() is mole number of component in mixture
    Specie::operator+=(pg); //YP: is it mean that sum all Y?  pg means all species Y, so after this, Y() means sumY

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();                      // Y1 means mole fraction
        const scalar Y2 = pg.Y() / this->Y(); // pg.Y() means all species corresponding value; Y2 means others mole fraction

        Tc_ = Y1 * Tc_ + Y2 * pg.Tc_;
        Vc_ = Y1 * Vc_ + Y2 * pg.Vc_;
        Zc_ = Y1 * Zc_ + Y2 * pg.Zc_;
        Pc_ = RR * Zc_ * Tc_ / Vc_;
        omega_ = Y1 * omega_ + Y2 * pg.omega_;
    }
}

template <class Specie>
inline void Foam::PengRobinsonS<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template <class Specie>
Foam::PengRobinsonS<Specie> Foam::operator+(
    const PengRobinsonS<Specie> &pg1,
    const PengRobinsonS<Specie> &pg2)
{
    Specie sp(
        static_cast<const Specie &>(pg1) + static_cast<const Specie &>(pg2));

    if (mag(sp.Y()) < small)
    {
        return PengRobinsonS<Specie>(
            sp,
            pg1.Tc_,
            pg1.Vc_,
            pg1.Zc_,
            pg1.Pc_,
            pg1.omega_

        );
    }
    else
    {
        const scalar Y1 = pg1.Y() / sp.Y();
        const scalar Y2 = pg2.Y() / sp.Y();

        const scalar Tc = Y1 * pg1.Tc_ + Y2 * pg2.Tc_;
        const scalar Vc = Y1 * pg1.Vc_ + Y2 * pg2.Vc_;
        const scalar Zc = Y1 * pg1.Zc_ + Y2 * pg2.Zc_;

        return PengRobinsonS<Specie>(
            sp,
            Tc,
            Vc,
            Zc,
            RR * Zc * Tc / Vc,
            Y1 * pg1.omega_ + Y2 * pg2.omega_);
    }
}

template <class Specie>
Foam::PengRobinsonS<Specie> Foam::operator*(
    const scalar s,
    const PengRobinsonS<Specie> &pg)
{
    return PengRobinsonS<Specie>(
        s * static_cast<const Specie &>(pg),
        pg.Tc_,
        pg.Vc_,
        pg.Zc_,
        pg.Pc_,
        pg.omega_);
}

template <class Specie>
Foam::PengRobinsonS<Specie> Foam::operator==(
    const PengRobinsonS<Specie> &pg1,
    const PengRobinsonS<Specie> &pg2)
{
    Specie sp(
        static_cast<const Specie &>(pg1) == static_cast<const Specie &>(pg2));

    const scalar Y1 = pg1.Y() / sp.Y();
    const scalar Y2 = pg2.Y() / sp.Y();

    const scalar Tc = Y2 * pg2.Tc_ - Y1 * pg1.Tc_;
    const scalar Vc = Y2 * pg2.Vc_ - Y1 * pg1.Vc_;
    const scalar Zc = Y2 * pg2.Zc_ - Y1 * pg1.Zc_;

    return PengRobinsonS<Specie>(
        sp,
        Tc,
        Vc,
        Zc,
        RR * Zc * Tc / Vc,
        Y2 * pg2.omega_ - Y1 * pg1.omega_);
}

// ************************************************************************* //
