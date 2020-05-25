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

#include "PengRobinsonM.H"
#include "specie.H"
#include "mathematicalConstants.H"
#include "OFstream.H"
#include "IFstream.H"
#include "string.H"
#define PRMIX

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class Specie>
inline Foam::autoPtr<Foam::PengRobinsonM<Specie>>
Foam::PengRobinsonM<Specie>::clone() const
{
    return autoPtr<PengRobinsonM<Specie>>(
        new PengRobinsonM<Specie>(*this));
}

template <class Specie>
inline Foam::autoPtr<Foam::PengRobinsonM<Specie>>
Foam::PengRobinsonM<Specie>::New(
    const dictionary &dict)
{
    return autoPtr<PengRobinsonM<Specie>>(
        new PengRobinsonM<Specie>(dict));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class Specie>
inline Foam::scalar Foam::PengRobinsonM<Specie>::rho(
    scalar p,
    scalar T) const
{
    const scalar Z = this->Z(p, T);
    return p / (Z * this->R() * T);
}

template <class Specie>
inline Foam::scalar Foam::PengRobinsonM<Specie>::psi(
    scalar p,
    scalar T) const
{
    const scalar Z = this->Z(p, T);

    return 1.0 / (Z * this->R() * T);
}

template <class Specie>
inline void Foam::PengRobinsonM<Specie>::ABZmix(
    scalar p,             //inp
    scalar T,             //inp
    scalarList &comp_inp, //inp
    scalar &Amixz,        //outp
    scalar &Bmixz,        //outp
    scalarList &Amix2,    //outp
    scalarList &Bz        //outp
    ) const
{
    scalarList Pr_sp(comp_inp.size(), Zero);
    scalarList Tr_sp(comp_inp.size(), Zero);
    scalarList Az(comp_inp.size(), Zero);
    scalarListList Aijz(comp_inp.size());
    forAll(Aijz, i)
    {
        Aijz[i].resize(comp_inp.size());
    }

    forAll(comp_inp, spid)
    {
        Pr_sp[spid] = p / single_specie[spid].Pc();
        Tr_sp[spid] = T / single_specie[spid].Tc();
    }

    forAll(comp_inp, spid)
    {
        single_specie[spid].AzBz(Pr_sp[spid], Tr_sp[spid], Az[spid], Bz[spid]);
    }
    forAll(comp_inp, i)
    {
        forAll(comp_inp, j)
        {
            Aijz[i][j] = sqrt(Az[i] * Az[j]) * (1.0 - bico[i][j]);
        }
    }
    Amixz = 0.0;
    Bmixz = 0.0;
    forAll(comp_inp, i)
    {
        Bmixz += comp_inp[i] * Bz[i];
        Amix2[i] = 0.0;
        forAll(comp_inp, j)
        {
            Amixz += comp_inp[i] * comp_inp[j] * Aijz[i][j];
            Amix2[i] += comp_inp[j] * Aijz[j][i];
        }
    }
}

//can be used for phases or mixtures, depending on the Azin and Bzin, only one root
template <class Specie>
inline Foam::scalar Foam::PengRobinsonM<Specie>::Zmix(
    int phasei, //input
    scalar Azin,
    scalar Bzin) const
{
    scalar a2 = Bzin - 1.0;
    scalar a1 = Azin - 2.0 * Bzin - 3.0 * sqr(Bzin);
    scalar a0 = -Azin * Bzin + sqr(Bzin) + pow3(Bzin);

    scalar Q = (3.0 * a1 - a2 * a2) / 9.0;
    scalar Rl = (9.0 * a2 * a1 - 27.0 * a0 - 2.0 * a2 * a2 * a2) / 54.0;

    scalar Q3 = Q * Q * Q;
    scalar D = Q3 + Rl * Rl;
    //printf("Zmix_root00= %+5.2i,%+5.10f,%+5.8f,%+5.8f\n", phasei,D,Azin,Bzin);
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
        if (phasei == 0) //gas
        {
            root = max(r1, max(r2, r3));
        }
        else //liquid
        {
            if (0) //Neg_rootFlag == 0)//!!!!!!!!!!!!!!!!
            {
                scalar root1 = min(r1, min(r2, r3));
                if (root1 <= 0.0)
                {
                    root1 = 1.0e-03;
                }
                root = root1;
            }
            else
            {
                scalar root1 = min(r2, r3);
                if (root1 <= 0.0)
                {
                    root1 = max(r2, r3);
                }
                scalar root2 = min(r1, root1);
                if (root2 <= 0.0)
                {
                    root2 = max(r1, root1);
                }
                root = root2;
            }
            if (root < 1.0e-03)
            {
                root = max(r1, max(r2, r3));
            }
        }
        //printf("Zmix_root11= %+5.2i,%+5.8f,%+5.8f,%+5.8f,%+5.8f\n", phasei,r1,r2,r3,root);
    }
    else
    {
        // One root is real: in supercritical zone
        scalar D05 = sqrt(D);
        scalar RlD = Rl + D05;
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
        if (D05 > Rl)
        {
            Tl = -pow(mag(Rl - D05), 1.0 / 3.0);
        }
        else
        {
            Tl = pow(Rl - D05, 1.0 / 3.0);
        }
        root = S + Tl - a2 / 3;
        //printf("Zmix_root22= %+5.2i,%+5.8f\n", phasei,root);
    }

    if (root <= 0.0)
    {
        root = 1.0;
        printf("Zmix_root4444444= %+5.2i,%+5.8f\n", phasei, root);
        FatalErrorInFunction
            << exit(FatalError);
    }
    return root;
}

// PY: fugacity coefficient of species in each phase
template <class Specie>
inline void Foam::PengRobinsonM<Specie>::fugcoef_phase(
    int phasei,           //input
    scalar p,             //input
    scalar T,             //input
    scalarList &comp_inp, //input
    scalarList &fugcoef   //output
    ) const
{
    scalar Amixz = 0.0;
    scalar Bmixz = 0.0;
    scalar Zmix = 0.0;
    scalarList Amix2(comp_inp.size(), Zero);
    scalarList Bz(comp_inp.size(), Zero);
    scalarList comp1(comp_inp.size(), Zero);
    const scalar small8 = 1.e-08; //PY: need to be noticed
    //printf("fugcoef_phase_0000= %+5.8f,%+5.8f\n", comp_inp[0],comp_inp[1]);

    this->ABZmix(p, T, comp_inp, Amixz, Bmixz, Amix2, Bz);

    Zmix = this->Zmix(phasei, Amixz, Bmixz);

    scalar Zmix0 = this->Zmix(0, Amixz, Bmixz);
    scalar comptrifile = 0.7;

    if (phasei == 1)
    {
        while (Zmix < 0.0)
        {
            comp1[0] = comptrifile;
            for (int i = 1; i < comp_inp.size(); i++)
                comp1[i] = (1.0 - comptrifile) / (comp_inp.size() - 1);
            this->ABZmix(p, T, comp1, Amixz, Bmixz, Amix2, Bz);
            Zmix = this->Zmix(phasei, Amixz, Bmixz);
            comptrifile = comptrifile + 0.01;
            if (comptrifile >= 1.0)
            {
                Zmix = Zmix0;
            }
        }
    }
    //printf("fugcoef_phase_1111= %+5.8f,%+5.8f,%+5.8f\n", Amixz,Bmixz,Zmix);
    /*
    do j=1,3
      if(zfactor(j)>Bmixz) then
         zfactor_0= zfactor(j)

 	 call fugacitycoef(zfactor_0,Amixz,Bmixz,Amix2,Bz,fugcoef)
         
  	 call gibbs(comp_inp,fugcoef,gibbsenergy)
  	 gibbsenergy0(j)=gibbsenergy
      else 
         gibbsenergy0(j)=1.d10
      end if  
   end do
   !
   g_min = 1.d20
   do i=1,3
      if (gibbsenergy0(i) < g_min ) then
         g_min = gibbsenergy0(i)
         i_min = i
      end if 
   end do
   zfactor_overall=zfactor(i_min)
*/

    scalar cc0 = 2.0 * sqrt(2.0);
    scalar cc1 = 1.0 + sqrt(2.0);
    scalar cc2 = 1.0 - sqrt(2.0);

    scalar zfactor00 = Zmix;
    scalar zdiff = zfactor00 - Bmixz;

    if (zdiff < 0.0)
    {
        zfactor00 = Bmixz + small8;
        zdiff = small8;
    }

    forAll(comp_inp, spid)
    {
        fugcoef[spid] = Bz[spid] / Bmixz * (zfactor00 - 1.0) - log(zdiff) - Amixz / (cc0 * Bmixz) * (2.0 * Amix2[spid] / Amixz - Bz[spid] / Bmixz) * log((zfactor00 + cc1 * Bmixz) / (zfactor00 + cc2 * Bmixz));

        fugcoef[spid] = exp(fugcoef[spid]);
    }
    //printf("fugcoef_phase_2222= %+5.15f,%+5.15f\n", fugcoef[0], fugcoef[1] );
}

template <class Specie>
inline void Foam::PengRobinsonM<Specie>::wilsoneq(
    scalar p,                 //inp
    scalar T,                 //inp
    scalarList &equl_constant //outp
    ) const
{
    scalarList dummy(equl_constant.size(), Zero);

    forAll(equl_constant, spid)
    {
        //dummy[spid] = (7.0/3.0)*(1.0 - Tc_sp[spid]/T)* (1.0 + omega_sp[spid]);
        //equl_constant[spid] = Pc_sp[spid]*pow(10.0,dummy[spid])/p ;

        dummy[spid] = 5.373 * (1.0 + single_specie[spid].omega()) * (1.0 - single_specie[spid].Tc() / T);
        equl_constant[spid] = (single_specie[spid].Pc() / p) * exp(dummy[spid]);
    }
}

template <class Specie>
inline Foam::scalar Foam::PengRobinsonM<Specie>::OAD217_function( //????? useless
    scalar XX,
    scalar YY) const
{
    int test0 = 1;
    int OAD217 = 0;

    scalar dumdif = (XX - YY);
    scalar dumsum = (XX + YY) / 2.0;
    scalar dumdif_abs = dumdif;
    if (dumdif <= 0.0)
    {
        dumdif_abs = -dumdif;
    }
    scalar dumsum_abs = dumsum;
    if (dumsum <= 0.0)
    {
        dumsum_abs = -dumsum;
    }

    if (dumdif_abs < (test0 * dumsum_abs) || dumdif == 0.0)
    {
        OAD217 = 0;
    }
    else if (dumdif > 0.0)
    {
        OAD217 = 1;
    }
    else if (dumdif < 0.0)
    {
        OAD217 = -1;
    }
    return OAD217;
}

template <class Specie>
inline void Foam::PengRobinsonM<Specie>::OAD013_function( //noramlize
    scalarList &comp_in,  // input
    scalarList &comp_out, // output
    int &ierr             //outpur
    ) const
{
    ierr = 0;
    scalar sum = 0.0;
    forAll(comp_in, spid)
    {
        sum += comp_in[spid];
    }

    int OAD217_val = this->OAD217_function(sum, 0.0);
    if (OAD217_val == 0)
    {
        ierr = 1;
        return;
    }

    forAll(comp_in, spid)
    {
        comp_out[spid] = comp_in[spid] / sum;
    }
}

template <class Specie>
inline void Foam::PengRobinsonM<Specie>::OAD012_function( //solve vaporf
    scalar p,                  //input
    scalar T,                  //input
    scalarList &comp_inp,      //input
    scalarList &equalconstant, //input
    scalarList &comp_liq,      //output
    scalarList &comp_gas,      //output
    scalar &vaporf,            //output
    int &ierr                  //output
    ) const
{
    int miter = 1;

    ierr = 0;
    scalar crita1 = -1.234567e+10;
    scalar crita2 = 1.234567e+10;
    scalar tolerance = 1.0e-14;
    scalar coeff1 = 0.1;
    scalar coeff2 = 0.5; //if the output vaporf=0.5, try to enlarge crita2 and give more time for tpn convergence;
    int nitmax = 500;    //100
    scalar dumm = 0.0;
    scalarList comp_liq0(comp_inp.size(), Zero);
    scalarList comp_gas0(comp_inp.size(), Zero);

    vaporf = 0.9;
    ierr = 0;

    scalar ad = -crita1, ag = crita1, tetad = -crita1, tetag = crita1;
    //which can be replaced by input comp_gas or comp_liq
    forAll(comp_inp, spid)
    {
        comp_liq[spid] = equalconstant[spid] - 1.0;
        //printf("OAD012function0******= %+5.8f %+5.8f \n", comp_liq[spid],equalconstant[spid]);
        if (comp_inp[spid] > 0.0)
        {
            if (comp_liq[spid] > 0.0)
            {
                ag = max(ag, -1.0 / comp_liq[spid]);
                dumm = (comp_inp[spid] - 1.0) / comp_liq[spid] + comp_inp[spid];
                tetag = max(tetag, dumm);
            }
            else if (comp_liq[spid] < 0.0)
            {
                ad = min(ad, -1.0 / comp_liq[spid]);
                dumm = (comp_inp[spid] - 1.0) / comp_liq[spid];
                tetad = min(tetad, dumm);
            }
            // Remark !!!!! what we will do for comp_liq[spid] = 0.0
        }
    }
    //printf("OAD012function1******= %+5.8f %+5.8f %+5.8f %+5.8f %+5.8f\n", ag,ad,dumm,tetag,tetad);
    if (ad > crita2 && ag < -crita2) //Not very sure, no physics, please double check
    {
        vaporf = 0.5;
        forAll(comp_inp, spid)
        {
            comp_liq[spid] = comp_inp[spid];
            comp_gas[spid] = comp_inp[spid];
        }
        ierr = 10; //change to 2
    }
    else if (ad == -crita1) //Not very sure, no physics, please double check
    {
        vaporf = 0.9999999;
        forAll(comp_inp, spid)
        {
            comp_gas[spid] = comp_inp[spid];
            comp_liq0[spid] = comp_gas[spid] / equalconstant[spid];
        }
        this->OAD013_function(comp_liq0, comp_liq, ierr);
        ierr = 4; //change to 0
    }
    else if (ag == crita1) //Not very sure, no physics, please double check
    {
        vaporf = 1.0e-07;
        forAll(comp_inp, spid)
        {
            comp_liq[spid] = comp_inp[spid];
            comp_gas0[spid] = comp_liq[spid] * equalconstant[spid];
        }
        this->OAD013_function(comp_gas0, comp_gas, ierr);
        ierr = 3; //change to 1
    }
    else
    {
        // the initialization indicates the state is far from two-phase regime;
        if (tetag < 0.0)
        {
            scalar f0 = -1.0;
            forAll(comp_inp, spid)
            {
                f0 += comp_inp[spid] * equalconstant[spid];
            }
            if (f0 > 0.0)
            {
                tetag = 0.0;
            }
            else if (f0 < 0.0)
            {
                tetad = 0.0;
            }
            else
            {
                vaporf = 1.0e-09;
                miter = 2;
            }
        }

        /* 
        if (tetad > 1.0) //Not very sure, no physics, please double check
        {
            scalar f1 = 1.0;
            for (label spid = 0; spid < num_sp; spid++)
            {
                f1 -= comp_inp[spid] / equalconstant[spid];
            }
            if (f1 < 0.0) 
            {
		tetad = 1.0;
            }
	    else if (f1 > 0.0)
            {
		tetag = 1.0;
	    }
            else
            {
		vaporf = 0.9999999;
		miter = 2;
	    }
        }
*/
        if (miter == 1)
        {
            vaporf = (tetag + tetad) / 2.0;
            scalar dx = tetad - tetag;
            //printf("OAD012function_00000= %+5.8f %+5.8f\n", vaporf,dx);
            //solve Rachford-rice equation
            int niter = 0;
            scalar fonc = 0.0, deri = 0.0, test = 0.0;
            while (niter < nitmax)
            {
                fonc = 0.0, deri = 0.0;
                forAll(comp_inp, spid)
                {
                    dumm = comp_liq[spid] / (1.0 + comp_liq[spid] * vaporf);
                    fonc = fonc + dumm * comp_inp[spid];
                    deri = deri - sqr(dumm) * comp_inp[spid];
                }
                dumm = (vaporf - ag) * (ad - vaporf);
                deri = dumm * deri - (2.0 * vaporf - ag - ad) * fonc;
                fonc = fonc * dumm;
                test = ((vaporf - tetad) * deri - fonc) * ((vaporf - tetag) * deri - fonc);

                if (fonc > 0.0)
                {
                    tetag = vaporf;
                }
                else
                {
                    tetad = vaporf;
                }

                if (test >= 0.0)
                {
                    dx = coeff1 * (tetad - tetag); //if TPn flash converged not very well, reduce dx
                    vaporf = tetag + dx;
                }
                else
                {
                    dx = coeff2 * fonc / deri;
                    vaporf = vaporf - dx;
                }

                scalar absdx = dx;
                if (dx <= 0.0)
                {
                    absdx = -dx;
                }

                if (absdx < tolerance)
                {
                    niter = nitmax;
                    //printf("OAD012function_iteration= %+5.3i %+5.20f %+5.6f\n", niter,absdx,vaporf);
                    break;
                }
                niter++;
            }
            if (niter < nitmax)
            {
                ierr = 2;
            }
        }

        forAll(comp_inp, spid)
        {
            comp_liq[spid] = comp_inp[spid] / (1.0 + comp_liq[spid] * vaporf);
            comp_gas[spid] = equalconstant[spid] * comp_liq[spid];
        }

        if (ierr == 0)
        {
            this->OAD013_function(comp_liq0, comp_liq, ierr);
            this->OAD013_function(comp_gas0, comp_gas, ierr);
        }
    }
    //printf("OAD012function_equalconstant= %+5.8f %+5.8f %+5.8f\n", equalconstant[0],equalconstant[1],equalconstant[2]);
    //FatalErrorInFunction
    //<< exit(FatalError);
}

template <class Specie>
inline void Foam::PengRobinsonM<Specie>::gibbs_phase(
    int phasei, //inp
    scalar p,
    scalar T,
    scalarList &comp_inp,
    scalar &gibbsenergy //outp
    ) const
{
    scalarList fugcoef(comp_inp.size(), Zero);
    gibbsenergy = 0.0;
    this->fugcoef_phase(phasei, p, T, comp_inp, fugcoef);

    forAll(comp_inp, spid)
    {
        gibbsenergy += comp_inp[spid] * log(comp_inp[spid] * fugcoef[spid]);
    }
    //printf("gibbs_phase= %+5.2i,%+5.8f\n", phasei,gibbsenergy);
}

template <class Specie>
inline void Foam::PengRobinsonM<Specie>::TPn_flash // PY:with a course initialization
    (
        scalar p,                 //input
        scalar T,                 //input
        scalarList &comp,         //input
        scalarList &comp_liq,     //output
        scalarList &comp_gas,     //output
        scalar &vaporfra,         //output
        scalarList &equalconstant //output
        ) const
{
    const int maxiter = 1000;
    const int igas = 0, iliquid = 1;
    int ierr = 0;
    int iaccelation = 1;
    scalar tol_fug = 1.0e-13; //Convergence criteria
    scalar tol_fug2 = 1.0;
    scalar coeff_tpn = 1.0;
    scalar epsilon1 = 1.0e+10;

    forAll(comp, i)
    {
        equalconstant[i] = 0;
    }
    scalar gibbsenergy = 0.0;
    scalarList fugcoef_liq(comp.size(), Zero), fugcoef_gas(comp.size(), Zero); //, fugcoef_mix[num_sp] = {0.0};
    scalarList comp_liq0(comp.size(), Zero), comp_gas0(comp.size(), Zero);     //, fugcoef_mix[num_sp] = {0.0};

    forAll(comp, i)
    {
        comp_liq[i] = comp[i]; //use this in PT_x flash, cancle this in hn flash
        comp_gas[i] = comp[i]; //use this in PT_x flash, cancle this in hn flash
        comp_liq0[i] = comp_liq[i];
        comp_gas0[i] = comp_gas[i];
    }
    //scalar tpdvalue = 0;
    //this->TPD(p, T, comp, comp_liq, comp_gas, tpdvalue);
    this->fugcoef_phase(iliquid, p, T, comp_liq, fugcoef_liq);
    this->fugcoef_phase(igas, p, T, comp_gas, fugcoef_gas);
    //this->fugcoef_mixture(p, T, comp, fugcoef_mix, gibbsenergy);
    this->wilsoneq(p, T, equalconstant);
    //printf("TPflash_OAD012_function0 = %+5.8f %+5.8f %+5.8f %+5.8f\n", p,T,equalconstant[0],equalconstant[1]);

    /*for (label spid = 0; spid < num_sp; spid++)//0725
    {
	if (fugcoef_liq[spid] < 1.0e-20)
	{
	     fugcoef_liq[spid] = 1.0e-20;
	}
	if (fugcoef_gas[spid] < 1.0e-20)
	{
	     fugcoef_gas[spid] = 1.0e-20;
	}	
	if (equalconstant[spid] < 1.0e-20)
	{
	     equalconstant[spid] = 1.0e-20;
	}		
	if (equalconstant[spid] > 1.0e+20)
	{
	     equalconstant[spid] = 1.0e+20;
	}			
    }*/

    scalarList re_iter1(comp.size());
    forAll(re_iter1, i)
        re_iter1[i] = 1;
    scalar gibbs_iter = gibbsenergy;
    scalar vaporfra_iter1 = 0.0;
    scalar vaporfra_iter0 = 0.0;
    //scalar vaporf0 = vaporfra;
    for (label iter = 0; iter < maxiter; iter++)
    {
        /*
	scalar equalKmax = 0.0;
        for (label spid = 0; spid < num_sp; spid++)
        {
	    equalKmax = max(equalKmax, equalconstant[spid]);
	}
	if (equalKmax > 1.0e+10)
	{
	    this->wilsoneq(p, T, equalconstant);
	    this->OAD012_function(p, T, comp, equalconstant, comp_liq, comp_gas, vaporfra, ierr);
	    scalar detvap = vaporfra - vaporf0;
	    if (detvap < 0.0)
	    {
		detvap = -detvap;
	    }
	    if (detvap > 0.1)
	    {
		vaporfra = vaporf0 + 0.1*(vaporfra);
	    }

	    for (label spid = 0; spid < num_sp; spid++)
	    {
		comp_liq[spid]=comp[spid]*(1.0/equalconstant[spid]);
		//comp_gas[spid]=equalconstant[spid]*comp[spid]; 
		comp_gas[spid]=1.0-comp_liq[spid];           
	    }
	    break;
	}
        */
        vaporfra_iter0 = vaporfra;
        this->OAD012_function(p, T, comp, equalconstant, comp_liq, comp_gas, vaporfra, ierr);
        //printf("TPflash_OAD012_function1 = %+5.2i %+5.8f %+5.8f %+5.8f\n", ierr,vaporfra,equalconstant[0],equalconstant[1]);

        // only for Hn flash at high temperature and pressure conditions
        if (ierr == 4 || ierr == 3)
        {
            forAll(comp, spid)
            {
                comp_liq[spid] = comp_liq0[spid];
                comp_gas[spid] = comp_gas0[spid];
            }
            if (vaporfra > 1.0)
            {
                vaporfra = 0.9999999;
            }
            else if (vaporfra < 0.0)
            {
                vaporfra = 1.0e-07;
            }
            break;
        }

        if (vaporfra > 0.99999999) //!!!!&& flag_tp > 0
        {
            vaporfra = 0.99999999;
            forAll(comp, spid)
            {
                comp_liq[spid] = comp_liq0[spid];
                comp_gas[spid] = comp_gas0[spid];
            }
            break;
        }
        if (vaporfra < 1.0e-08) //!!!!&& flag_tp > 0
        {
            vaporfra = 1.0e-08;
            forAll(comp, spid)
            {
                comp_liq[spid] = comp_liq0[spid];
                comp_gas[spid] = comp_gas0[spid];
            }
            break;
        }

        vaporfra_iter1 = vaporfra;
        if ((ierr == 10) && iter > 2)
        {
            //this->TPD(p, T, comp, comp_liq, comp_gas, tpdvalue);
            this->fugcoef_phase(iliquid, p, T, comp_liq, fugcoef_liq);
            this->fugcoef_phase(igas, p, T, comp_gas, fugcoef_gas);
            //this->fugcoef_mixture(p, T, comp, fugcoef_mix, gibbsenergy);
            forAll(comp, spid)
            {
                equalconstant[spid] = fugcoef_liq[spid] / fugcoef_gas[spid];
            }
            this->OAD012_function(p, T, comp, equalconstant, comp_liq, comp_gas, vaporfra, ierr);
            // Here, for tp diagram, please reconsider it
            forAll(comp, spid)
            {
                comp_liq[spid] = comp_liq0[spid];
                comp_gas[spid] = comp_gas0[spid];
            }

            if (vaporfra > 1.0)
            {
                vaporfra = 0.9999999;
            }
            if (vaporfra < 0.0)
            {
                vaporfra = 1.0e-07;
            }
            break;
        }

        this->fugcoef_phase(iliquid, p, T, comp_liq, fugcoef_liq);
        this->fugcoef_phase(igas, p, T, comp_gas, fugcoef_gas);
        //printf("TPflash_OAD012_function3 = %+5.8f %+5.8f %+5.8f %+5.8f %+5.8f %+5.8f\n",fugcoef_liq[0], fugcoef_gas[0],equalconstant[0],equalconstant[1],comp_liq[0],comp_gas[0]);
        scalar errr = 0.0;
        scalarList re_iter0(comp.size(), Zero);
        forAll(comp, spid)
        {
            re_iter0[spid] = re_iter1[spid];
            fugcoef_gas[spid] = comp_gas[spid] * fugcoef_gas[spid] * p;
            fugcoef_liq[spid] = comp_liq[spid] * fugcoef_liq[spid] * p;
            errr = errr + sqr(fugcoef_liq[spid] / fugcoef_gas[spid] - 1.0);
            equalconstant[spid] = equalconstant[spid] * ((fugcoef_liq[spid] / fugcoef_gas[spid]) * coeff_tpn);
            re_iter1[spid] = fugcoef_liq[spid] / fugcoef_gas[spid];
        }
        epsilon1 = errr;
        /*for (label spid = 0; spid < num_sp; spid++)//0725
        {
	     if (fugcoef_liq[spid] < 1.0e-20)
	     {
	        fugcoef_liq[spid] = 1.0e-20;
	     }
	     if (fugcoef_gas[spid] < 1.0e-20)
	     {
	        fugcoef_gas[spid] = 1.0e-20;
	     }	
	     if (equalconstant[spid] < 1.0e-20)
	     {
	        equalconstant[spid] = 1.0e-20;
	     }		
	     if (equalconstant[spid] > 1.0e+20)
	     {
	        equalconstant[spid] = 1.0e+20;
	     }			
        }*/
        //printf("TPflash_epsilon1****** = %+5.10f %+5.10f %+5.10f %+5.10f\n",equalconstant[0],equalconstant[1], comp_liq[0], comp_gas[0]);
        scalar epsilon2 = 1.0e+10;
        if (errr >= 1.0e+05 && iter > 10)
        {
            this->fugcoef_phase(iliquid, p, T, comp_liq, fugcoef_liq);
            this->fugcoef_phase(igas, p, T, comp_gas, fugcoef_gas);
            //this->fugcoef_mixture(p, T, comp, fugcoef_mix, gibbsenergy);
            forAll(comp, spid)
            {
                equalconstant[spid] = fugcoef_liq[spid] / fugcoef_gas[spid];
            }

            this->OAD012_function(p, T, comp, equalconstant, comp_liq, comp_gas, vaporfra, ierr);
            scalar errr2 = 0.0;
            forAll(comp, spid)
            {
                re_iter0[spid] = re_iter1[spid];
                fugcoef_gas[spid] = comp_gas[spid] * fugcoef_gas[spid] * p;
                fugcoef_liq[spid] = comp_liq[spid] * fugcoef_liq[spid] * p;
                errr2 = errr2 + sqr(fugcoef_liq[spid] / fugcoef_gas[spid] - 1.0);
                equalconstant[spid] = equalconstant[spid] * (fugcoef_liq[spid] / fugcoef_gas[spid]);
                re_iter1[spid] = fugcoef_liq[spid] / fugcoef_gas[spid];
            }
            epsilon2 = errr2;
            //printf("TPflash_epsilon44444444444 = %+5.15f\n",epsilon2);
        }

        if (epsilon1 < tol_fug || epsilon2 < tol_fug2)
        {
            //printf("TPflash_epsilon****** = %+5.2i %+5.8f %+5.6f %+5.15f\n",ierr, vaporfra, T, epsilon1);
            break;
        }

        //ASSM accelerating
        if (iaccelation == 2) //1)
        {
            if (iter > 2)
            {
                scalar sum_iter0 = 0.0, sum_iter1 = 0.0;
                scalarList lamda_coef(comp.size(), Zero);
                forAll(comp, spid)
                {
                    sum_iter0 += sqr(re_iter0[spid] - 1.0);
                    sum_iter1 += sqr(re_iter1[spid] - 1.0);
                }

                scalar ratioiter = sum_iter1 / sum_iter0;
                scalar detiter = vaporfra_iter1 - vaporfra_iter0;
                if (ratioiter > 0.8 && detiter < 0.1 && sum_iter1 > 1.0e-5 && sum_iter1 < 1.0e-3 && vaporfra_iter1 > 1.0e-07 && vaporfra_iter1 < 0.9999999)
                {
                    forAll(comp, spid)
                    {
                        lamda_coef[spid] = (re_iter0[spid] - 1.0) / (re_iter0[spid] - re_iter1[spid]);
                        equalconstant[spid] = equalconstant[spid] * pow(re_iter1[spid], lamda_coef[spid]);
                    }
                }
            }
        }
        scalar gibbsenergyl = 0.0, gibbsenergyg = 0.0;
        this->gibbs_phase(iliquid, p, T, comp_liq, gibbsenergyl);
        this->gibbs_phase(igas, p, T, comp_gas, gibbsenergyg);

        scalar gibbs_iter0 = gibbs_iter;
        gibbs_iter = gibbsenergyl + gibbsenergyg;

        if ((gibbs_iter > gibbs_iter0) && iter > 1)
        {
            iaccelation = 0;
        }
        else
        {
            iaccelation = 0; //1;
        }
    }
    //printf("TPflashEnd= %+5.8f %+5.8f\n", comp_liq[0], comp_gas[0]);
    //FatalErrorInFunction
    //<< exit(FatalError);
}

template <class Specie>
inline Foam::scalar Foam::PengRobinsonM<Specie>::mwmix // kg/mol
    (
        scalarList &comp_inp) const
{

    scalar mwmix = 0.0;
    forAll(comp_inp, spid)
    {
        mwmix += comp_inp[spid] * single_specie[spid].W(); //PY: Y() is mole number of component in mixture
    }
    return mwmix * 1.0e-03;
}

template <class Specie>
inline Foam::scalar Foam::PengRobinsonM<Specie>::rhomix_phase // PY:phase mixture density
    (
        int phasei, //input
        scalar p,
        scalar T,
        scalarList &comp_inp) const
{
    scalar rhomixPhase = 0.0;
    scalar Amixz = 0.0;
    scalar Bmixz = 0.0;
    scalarList Amix2(comp_inp.size(), Zero);
    scalarList Bz(comp_inp.size(), Zero);
    this->ABZmix(p, T, comp_inp, Amixz, Bmixz, Amix2, Bz);
    scalar Zmix = this->Zmix(phasei, Amixz, Bmixz);
    Info<<Amixz<<" "<<Bmixz<<endl;
    scalar mw_phase = this->mwmix(comp_inp); //kg/mol
    //RR: 10e+03*m3 Pa K-1 mol-1
    rhomixPhase = mw_phase * p / (Zmix * RR * 1.0e-03 * T); // kg/m3
    //printf("rhomixPhase1= %+5.8f %+5.8f %+5.8f %+5.8f %+5.8f\n", p, T, Zmix, mw_phase, rhomixPhase);

    return rhomixPhase;
}

template <class Specie>
inline Foam::scalar Foam::PengRobinsonM<Specie>::volmmix_phase //m3/mol
    (
        int phasei, //input
        scalar p,
        scalar T,
        scalarList &comp_inp) const
{
    scalar volspecific = 0.0;
    //!!!!!!
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

// ======================================================================//
// Evaluate the phase volume fraction, alpha				 //
// Input: pressure, temperature, vapor fraction				 //
// Output: phase compositions 				 		 //
//									 //
// History :   Creation.	     P. YI			2019     //
// ======================================================================//
template <class Specie>
inline Foam::scalar Foam::PengRobinsonM<Specie>::Evaluate_alpha(
    scalar press,         //inp
    scalar temp,          //inp
    scalar vaporfra,      //inp
    scalarList &comp_liq, //inp
    scalarList &comp_gas, //inp
    scalarList &comp_inp) const
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

// PY: need to use dynamc array to replace static array
template <class Specie>
inline void Foam::PengRobinsonM<Specie>::abmix(
    scalar T,
    scalarList &comp_inp,
    scalar &amix,
    scalar &bmix,
    scalar &dadt,
    scalar &dadt2) const
{
    scalarList m_omega_sp(comp_inp.size(), Zero);
    scalarList az(comp_inp.size(), Zero);
    scalarList bz(comp_inp.size(), Zero);
    scalarListList aijz(comp_inp.size());
    forAll(comp_inp, i)
    {
        aijz[i].resize(comp_inp.size());
    }
    scalarList alpha(comp_inp.size(), Zero);

    scalarList Pr_sp(comp_inp.size(), Zero);
    scalarList Tr_sp(comp_inp.size(), Zero);
    scalarList Az(comp_inp.size(), Zero);
    scalarListList Aijz(comp_inp.size());
    forAll(comp_inp, i)
    {
        Aijz[i].resize(comp_inp.size());
    }

    forAll(comp_inp, spid)
    {
        Tr_sp[spid] = T / single_specie[spid].Tc();
        ;
    }

    scalarList comp4(comp_inp.size(), Zero);
    forAll(comp_inp, spid)
    {
        single_specie[spid].azbz(Tr_sp[spid], m_omega_sp[spid], az[spid], bz[spid], alpha[spid], comp4[spid]);
    }

    scalarListList comp5(comp_inp.size());
    scalarListList comp6(comp_inp.size());
    scalarListList comp7(comp_inp.size());
    forAll(comp_inp, i)
    {
        comp5[i].resize(comp_inp.size());
        comp6[i].resize(comp_inp.size());
        comp7[i].resize(comp_inp.size());
    }
    forAll(comp_inp, spid)
    {
        forAll(comp_inp, spjd)
        {
            aijz[spid][spjd] = sqrt(az[spid] * az[spjd]) * (1.0 - bico[spid][spjd]);
            comp5[spid][spjd] = comp4[spid] * sqrt(alpha[spjd]);
            comp6[spid][spjd] = sqrt(alpha[spid]) * comp4[spjd];
            comp7[spid][spjd] = comp5[spid][spjd] + comp6[spid][spjd];
        }
    }

    // return these values
    amix = 0.0;
    bmix = 0.0;
    scalarList comp1(comp_inp.size(), Zero);
    scalarListList comp2(comp_inp.size());
    forAll(comp_inp, i)
    {
        comp2[i].resize(comp_inp.size());
    }

    forAll(comp_inp, spid)
    {
        forAll(comp_inp, spjd)
        {
            comp1[spid] = comp_inp[spid] * aijz[spjd][spid];
            comp2[spid][spjd] = 0.457240 * sqr(RR * 1.0e-03) * (bico[spid][spjd] - 1.0);
        }
        amix += comp1[spid] * comp_inp[spid]; //a*alpha
        bmix += comp_inp[spid] * bz[spid];    //b
    }

    scalarListList comp3(comp_inp.size());
    scalarListList dadt1(comp_inp.size());
    scalarListList dadt12(comp_inp.size());
    scalarListList dadt22(comp_inp.size());
    forAll(comp_inp, i)
    {
        comp3[i].resize(comp_inp.size());
        dadt1[i].resize(comp_inp.size());
        dadt12[i].resize(comp_inp.size());
        dadt22[i].resize(comp_inp.size());
    }
    forAll(comp_inp, spid)
    {
        forAll(comp_inp, spjd)
        {
            comp2[spid][spjd] = comp2[spid][spjd] * single_specie[spid].Tc() * single_specie[spjd].Tc();
            comp3[spid][spjd] = comp2[spid][spjd] * 1.0 / sqrt(single_specie[spid].Pc() * single_specie[spjd].Pc());
            dadt1[spid][spjd] = comp3[spid][spjd] * (1.0 / (2.0 * sqrt(T)));
            dadt12[spid][spjd] = dadt1[spid][spjd] * comp7[spid][spjd];
            dadt22[spid][spjd] = -0.457240 * sqr(RR * 1.0e-03) * (1.0 - bico[spid][spjd]) /
                                 (2.0 * sqrt(T)) * (m_omega_sp[spjd] * sqrt(alpha[spid]) * single_specie[spid].Tc() * sqrt(single_specie[spjd].Tc()) / sqrt(single_specie[spid].Pc() * single_specie[spjd].Pc()) * (-m_omega_sp[spjd] / (2.0 * sqrt(single_specie[spjd].Tc() * T * alpha[spid])) - 1.0 / (2.0 * T)) + m_omega_sp[spid] * sqrt(alpha[spjd]) * single_specie[spjd].Tc() * sqrt(single_specie[spid].Tc()) / sqrt(single_specie[spid].Pc() * single_specie[spjd].Pc()) * (-m_omega_sp[spjd] / (2.0 * sqrt(single_specie[spjd].Tc() * T * alpha[spjd])) - 1.0 / (2.0 * T)));
        }
    }

    dadt = 0.0;
    dadt2 = 0.0;
    scalarList comp8(comp_inp.size(), Zero);
    scalarList comp88(comp_inp.size(), Zero);
    forAll(comp_inp, spid)
    {
        forAll(comp_inp, spjd)
        {
            comp8[spid] = comp_inp[spid] * dadt12[spjd][spid];
            comp88[spid] = comp_inp[spid] * dadt22[spjd][spid];
        }
        dadt = dadt + comp8[spid] * comp_inp[spid];
        dadt2 = dadt2 + comp88[spid] * comp_inp[spid];
    }
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
template <class Specie>
inline Foam::scalar Foam::PengRobinsonM<Specie>::Hig_phase(
    scalar T,
    scalarList &comp_phase) const
{
    scalar H_ideal = 0;
    forAll(comp_phase, i)
    {
        H_ideal = H_ideal + comp_phase[i] * single_specie[i].Hig_phase(T);
    }

    return H_ideal;
    //return H_ideal/this->mwmix(comp_phase);
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
template <class Specie>
inline Foam::scalar Foam::PengRobinsonM<Specie>::Hig2_phase(
    scalar T,
    scalarList &comp_phase) const
{
    scalar H_ideal = 0;
    forAll(comp_phase, i)
    {
        H_ideal = H_ideal + comp_phase[i] * single_specie[i].Hig2_phase(T);
    }

    return H_ideal * 4.18680;;
}

// ======================================================================//
// Ideal phase mixture internal energy => U = H - PV			 //
// Input: temperature, phase compositions				 //
// Output: ideal phase mixture sie			 		 //
//									 //
// History :   Creation.	     P. YI			2019     //
// ======================================================================//
template <class Specie>
inline Foam::scalar Foam::PengRobinsonM<Specie>::sieig_phase(
    scalar T,
    scalarList &comp_phase //note, it is phase compositions
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

template <class Specie>
inline Foam::scalar Foam::PengRobinsonM<Specie>::siemix_phase // internal energy with unit of J/kg
    (
        int phasei, //input
        scalar p,
        scalar T,
        scalarList &comp_phase) const
{
    scalar amix = 0.0;
    scalar bmix = 0.0;
    scalar dadt = 0.0;
    scalar dadt2 = 0.0;
    scalar mw_phase = this->mwmix(comp_phase);                       //kg*mol-1
    scalar rho_phase = this->rhomix_phase(phasei, p, T, comp_phase); //kg*m-3

    scalar cc1 = 1.0 + sqrt(2.0);
    scalar cc2 = 1.0 - sqrt(2.0);

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
template <class Specie>
inline void Foam::PengRobinsonM<Specie>::cvcpig_phase(
    scalar T,
    scalarList &comp_phase,
    scalar &cpig,
    scalar &cvig) const
{
    scalar cp_ideal = 0.0;
    forAll(comp_phase, i)
    {
        cp_ideal += comp_phase[i] * single_specie[i].cp_phase(T);
    }
    scalar cv_ideal = cp_ideal - RR * 1.0 - 03;
    cpig = cp_ideal / this->mwmix(comp_phase);
    cvig = cv_ideal / this->mwmix(comp_phase);
}

template <class Specie>
inline Foam::scalar Foam::PengRobinsonM<Specie>::cvmix_phase // PY:phase mixture heat capability
    (
        int phasei, //input
        scalar p,
        scalar T,
        scalarList &comp_phase) const
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
inline Foam::scalar Foam::PengRobinsonM<Specie>::cpmix_phase // PY:phase mixture heat capability
    (
        int phasei, //input
        scalar p,
        scalar T,
        scalarList &comp_phase) const
{
    //!!!
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
inline void Foam::PengRobinsonM<Specie>::soundspeedmix // m/s
    (
        int phasei, //input
        scalar p,
        scalar T,
        scalarList &comp_inp,
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
// Phase mixture thermal conductivity, W/m.K: 				 //
// Input: temperature, pressure, phase compositions			 //
// Output: ideal phase mixture H			 		 //
//									 //
// Source: Chung,et al, 1988, 						 //
// Generalized Multiparameter Correlation for Nonpolar and Polar &	 //
// Fluid Transport Properties.			 			 //
// Book, THE PROPERTIES OF GASES & LIQUIDS, &  			 	 //
// Robert C Reid, John M. Prausnitz, Bruce E. Poling, P 522:		 //
//									 //
// History :   Creation.	     P. YI			2019     //
// ======================================================================//
template <class Specie>
inline void Foam::PengRobinsonM<Specie>::kappa_phase(
    int phasei,
    scalar p,
    scalar T,
    scalarList &comp_phase,
    scalar &thcond_phase,
    scalar &viscos) const
{

    scalarListList kij(comp_phase.size());
    scalarListList mwij(comp_phase.size());
    scalarListList w(comp_phase.size());
    scalarListList sigma(comp_phase.size());
    scalarListList e_k0(comp_phase.size());
    forAll(comp_phase, i)
    {
        kij[i].resize(comp_phase.size());
        mwij[i].resize(comp_phase.size());
        w[i].resize(comp_phase.size());
        sigma[i].resize(comp_phase.size());
        e_k0[i].resize(comp_phase.size());
    }

    scalarList k(comp_phase.size(), Zero); //need to be sure
    scalarList e_k(comp_phase.size(), Zero);
    scalarList mw_sp3(comp_phase.size(), Zero);
    scalarList Vc_sp3(comp_phase.size(), Zero);

    forAll(comp_phase, spid)
    {
        e_k[spid] = single_specie[spid].Tc() / 1.2593;
        mw_sp3[spid] = single_specie[spid].W();            //g/mol
        Vc_sp3[spid] = single_specie[spid].Vc() * 1.0e+03; //Vc is m3/kmol, Vc_sp3 is cm3/mol
    }
    forAll(comp_phase, spid)
    {
        forAll(comp_phase, spjd)
        {
            w[spid][spjd] = 0.50 * (single_specie[spid].omega() +single_specie[spjd].omega());
            sigma[spid][spjd] = sqrt((0.8090 * pow(Vc_sp3[spid], 1.0 / 3.0)) * (0.8090 * pow(Vc_sp3[spjd], 1.0 / 3.0)));
            //sigma(i,j)= bip_kij(i,j)*sqrt((0.809d0*vc(i)**(1.d0/3.d0))*(0.809d0*vc(j)**(1.d0/3.d0)))
            e_k0[spid][spjd] = sqrt(e_k[spid] * e_k[spjd]);
            //e_k0(i,j)= bip_kij(i,j)*sqrt(e_k(i)*e_k(j))
            //kij[spid][spjd] = sqrt(k[spid] * k[spjd]);
            kij[spid][spjd] = sqrt(single_specie[spid].kappa() * single_specie[spid].kappa());
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
    forAll(comp_phase, spid)
    {
        forAll(comp_phase, spjd)
        {
            sigma_mix3 += comp_phase[spid] * comp_phase[spjd] * pow3(sigma[spid][spjd]);
            sigma_mix2 += comp_phase[spid] * comp_phase[spjd] * sqr(sigma[spid][spjd]);
            e_k_mix +=  comp_phase[spid] * comp_phase[spjd] * e_k0[spid][spjd] * pow3(sigma[spid][spjd]);
        }
    }

        forAll(comp_phase, spid)
    {
        forAll(comp_phase, spjd)
        {
            
            w_mix   +=  comp_phase[spid] * comp_phase[spjd] * w[spid][spjd] * pow3(sigma[spid][spjd]);
            mu4     +=  comp_phase[spid] * comp_phase[spjd] * sqr(single_specie[spid].mu())*sqr(single_specie[spjd].mu() ) / (pow3(sigma[spid][spjd])*e_k0[spid][spjd] );
            k_mix   +=  comp_phase[spid] * comp_phase[spjd] * kij[spid][spjd];
            mw_binary += (comp_phase[spid] * comp_phase[spjd] * e_k0[spid][spjd] * sqr(sigma[spid][spjd]) * sqrt(mwij[spid][spjd])) ;// (e_k_mix * sigma_mix2);
        }
    }
    e_k_mix   /=sigma_mix3;
    w_mix     /=sigma_mix3;
    mu4       *=sigma_mix3* e_k_mix;
    mw_binary /=(e_k_mix * sigma_mix2);
    mw_binary = sqr(mw_binary);


    scalar volu_mix = sigma_mix3 / pow3(0.8090); //Eq.(9-5.43);Eq.(9-4.8) bug??
    //scalar volu_mix= pow3(sigma_mix3/0.8090);
    scalar tempc_mix = 1.25930 * e_k_mix;
    scalar mur_mix4 = pow4(131.30) * mu4 / sqr(tempc_mix * volu_mix); //volu_mix in cm3/mol, mu4 in debyes
    scalar t_mix = T / e_k_mix;
    scalar Fcm = 1.0 - 0.2756 * w_mix + 0.0590350 * mur_mix4 + k_mix;
    scalar omegav = 1.161450 * (pow(t_mix, -0.148740) + 0.524870 * exp(-0.773200 * t_mix) + 2.161780 * exp(-2.437870 * t_mix))+(-6.435e-4)*pow(t_mix, 0.148740)*sin(18.0323*pow(t_mix, -0.7683)-7.27371);//???GT...

    //1. mixture phases viscosity 40.785 ? 26.69 there is error in 9-5.24
    viscos = 1e-07 * 26.69 * Fcm * sqrt(mw_binary * T) / (sigma_mix2 * omegav);

    scalar mw_phase = this->mwmix(comp_phase);
    scalar Cv_phase = this->cvmix_phase(phasei, p, T, comp_phase); //J/kgK
    scalar alpha = Cv_phase * mw_phase / (RR * 1.0e-03) - 1.50;
    scalar tempr = T / tempc_mix;

    //compute specific volume, m3/mol
    scalar amix_phase = 0.0, bmix_phase = 0.0, dadt_phase = 0.0, dadt2 = 0.0;
    this->abmix(T, comp_phase, amix_phase, bmix_phase, dadt_phase, dadt2);
    scalar sv_phase = this->volmmix_phase(phasei, p, T, comp_phase); //m3/mol

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

    thcond_phase = 31.20 * viscos * psi * (1.0 / G2 + B[6] * yy) / mw_binary0 + qq * B[7] * sqr(yy) * sqrt(tempr) * G2;
    //printf("kappa_phase_thcond_phase %5.8f\n\n",thcond_phase);
}
