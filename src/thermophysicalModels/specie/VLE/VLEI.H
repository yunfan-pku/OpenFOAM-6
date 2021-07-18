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

#include "VLE.H"
#include "specie.H"
#include "mathematicalConstants.H"
#include "OFstream.H"
#include "IFstream.H"
#include "IOmanip.H"
#include "string.H"
#include "labelList.H"
#include "scalarMatrices.H"
#include <random>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class ThermoMixture>
inline Foam::autoPtr<Foam::VLE<ThermoMixture>>
Foam::VLE<ThermoMixture>::clone() const
{
    return autoPtr<VLE<ThermoMixture>>(
        new VLE<ThermoMixture>(*this));
}

template <class ThermoMixture>
inline Foam::autoPtr<Foam::VLE<ThermoMixture>>
Foam::VLE<ThermoMixture>::New(
    const dictionary& dict)
{
    return autoPtr<VLE<ThermoMixture>>(
        new VLE<ThermoMixture>(dict));
}

template <class ThermoMixture>
inline Foam::VLE<ThermoMixture>::VLE(
    const word& name,
    PtrList<SingleThermoType>& speciesData,
    const speciesTable& specieNames,
    const dictionary& thermoDict)
    : ThermoMixture(name, speciesData, specieNames, thermoDict)
{
}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class ThermoMixture>
inline autoPtr<scalarList> Foam::VLE<ThermoMixture>::wilsoneq(
    scalar p,
    scalar T) const
{
    autoPtr<scalarList> pequl_constant(new scalarList(this->N_));
    forAll(this->X_, i)
    {
        pequl_constant()[i] = ((*this)[i].Pc_ / p) * exp(5.373 * (1.0 + (*this)[i].omega_) * (1.0 - (*this)[i].Tc_ / T));
    }
    return pequl_constant;
}

template <class ThermoMixture>
inline autoPtr<scalarList> Foam::VLE<ThermoMixture>::myinitial(
    scalar p,
    scalar T) const
{
    std::default_random_engine e;
    std::uniform_real_distribution<double> u(0, 1);
    autoPtr<scalarList> pequl_constant(new scalarList(this->N_));
    scalarList x(this->N_), y(this->N_);
    double lowb = -1 + 2 * this->X_[0], upb = 1;
    double lb, ub;
    double sx = 0, sy = 0;
    for (int i = 0; i < this->N_ - 1; i++)
    {
        lb = lowb > 0 ? lowb : 0;
        ub = upb > 2 * this->X_[i] ? 2 * this->X_[i] : upb;
        x[i] = u(e) * (ub - lb) + lb;
        y[i] = 2 * this->X_[i] - x[i];
        //if(x[i]<0||y[i]<0||x[i]>1||y[i]>1)
        //cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        lowb += (2 * this->X_[i + 1] - this->X_[i]);
        upb -= x[i];
        sx += x[i];
        sy += y[i];
        if (x[i] != 0)
            pequl_constant()[i] = y[i] / x[i];
        else
        {
            pequl_constant()[i] = 100000000;
            /* code */
        }

        //cout<<z[i]<<","<<x[i]<<","<<y[i]<<","<<2*z[i]-x[i]-y[i]<<endl;
    }

    x[this->N_ - 1] = 1 - sx;
    y[this->N_ - 1] = 1 - sy;
    if (x[this->N_ - 1] != 0)
        pequl_constant()[this->N_ - 1] = y[this->N_ - 1] / x[this->N_ - 1];
    else
    {
        pequl_constant()[this->N_ - 1] = 100000000;
        /* code */
    }
    //cout<<z[n-1]<<","<<x[n-1]<<","<<y[n-1]<<","<<2*z[n-1]-x[n-1]-y[n-1]<<endl;
    return pequl_constant;
}

template <class ThermoMixture>
inline label Foam::VLE<ThermoMixture>::noramlize(scalarList& X_in) const //noramlize, error return 1
{
    scalar sum = 0.0;
    forAll(X_in, i)
    {
        sum += X_in[i];
    }
    if (sum == 0)
    {
        return 1;
    }
    forAll(X_in, i)
    {
        X_in[i] /= sum;
    }
    return 0;
}

template <class ThermoMixture>
inline label Foam::VLE<ThermoMixture>::solve_vaporf( //solve vaporf return error message,requir solution.equalconstant
    scalar p,
    scalar T,
    solution& sol) const
{
    int miter = 1;

    label ierr = 0;
    const scalar ninfty = -great;
    const scalar pinfty = great;
    const scalar tolerance = 1.0e-14;
    const scalar coeff1 = 0.1;
    const scalar coeff2 = 0.5; //if the output vaporf=0.5, try to enlarge crita2 and give more time for tpn convergence;
    const label nitmax = 500;  //100

    scalar temp = 0.0;

    scalarList X_liq0(this->N_);
    scalarList X_gas0(this->N_);

    scalarList KM1(this->N_);

    sol.vaporfra = 0.9;
    ierr = 0;

    scalar upperv = pinfty, lowerv = ninfty, upperl = pinfty, lowerl = ninfty;

    sol.X_liq.reset(new scalarList(this->N_));
    sol.X_gas.reset(new scalarList(this->N_));

    //which can be replaced by input comp_gas or comp_liq
    forAll(this->X_, i)
    {
        KM1[i] = sol.equalconstant()[i] - 1.0;
        //printf("OAD012function0******= %+5.8f %+5.8f \n", comp_liq[spid],equalconstant[spid]);
        if (this->X_[i] > 0.0)
        {
            if (KM1[i] > 0.0)
            {
                lowerv = max(lowerv, -1.0 / KM1[i]);
                lowerl = max(lowerl, (this->X_[i] - 1.0) / KM1[i] + this->X_[i]);
            }
            else if (KM1[i] < 0.0)
            {
                upperv = min(upperv, -1.0 / KM1[i]);
                upperl = min(upperl, (this->X_[i] - 1.0) / KM1[i]);
            }
            // Remark !!!!! what we will do for comp_liq[spid] = 0.0
        }
    }
    if (upperv > pinfty && lowerv < ninfty) //Not very sure, no physics, please double check
    {
        FatalErrorInFunction << "vaporf 1" << exit(FatalError);
    }
    else if (upperv == pinfty) //Not very sure, no physics, please double check
    {

        sol.vaporfra = 0.9999999;
        forAll(this->X_, i)
        {
            sol.X_gas()[i] = this->X_[i];
            sol.X_liq()[i] = sol.X_gas()[i] / sol.equalconstant()[i];
        }
        this->noramlize(sol.X_liq());
        ierr = 4; //change to 0

        //FatalErrorInFunction << "vaporf 2" << exit(FatalError);
    }
    else if (lowerv == ninfty) //Not very sure, no physics, please double check
    {

        sol.vaporfra = 1.0e-07;
        forAll(this->X_, i)
        {
            sol.X_liq()[i] = this->X_[i];
            sol.X_gas()[i] = sol.X_liq()[i] * sol.equalconstant()[i];
        }
        this->noramlize(sol.X_gas());
        ierr = 3; //change to 1
        //FatalErrorInFunction << "vaporf 3" << exit(FatalError);
    }
    else
    {
        // the initialization indicates the state is far from two-phase regime;
        if (lowerl < 0.0)
        {
            scalar f0 = -1.0;
            forAll(this->X_, i)
            {
                f0 += this->X_[i] * sol.equalconstant()[i];
            }
            if (f0 > 0.0) //??
            {
                lowerl = 0.0;
            }
            else if (f0 < 0.0) //??
            {
                upperl = 0.0;
            }
            else //??
            {
                sol.vaporfra = 1.0e-09;
                miter = 2;
            }
        }

        if (miter == 1)
        {
            sol.vaporfra = (upperl + lowerl) / 2.0;
            scalar dx = upperl - lowerl;

            //solve Rachford-rice equation
            int niter = 0;
            scalar fonc = 0.0, deri = 0.0, test = 0.0;
            while (niter < nitmax)
            {
                fonc = 0.0, deri = 0.0;
                forAll(this->X_, i)
                {
                    temp = KM1[i] / (1.0 + KM1[i] * sol.vaporfra);
                    fonc = fonc + temp * this->X_[i];
                    deri = deri - temp * temp * this->X_[i];
                }
                temp = (sol.vaporfra - lowerv) * (upperv - sol.vaporfra);
                deri = temp * deri - (2.0 * sol.vaporfra - lowerv - upperv) * fonc;
                fonc = fonc * temp; //guess: prevent from denominator equal to zero
                test = ((sol.vaporfra - upperl) * deri - fonc) * ((sol.vaporfra - lowerl) * deri - fonc);

                if (fonc > 0.0)
                {
                    lowerl = sol.vaporfra;
                }
                else
                {
                    upperl = sol.vaporfra;
                }

                if (test >= 0.0)
                {
                    dx = coeff1 * (upperl - lowerl); //if TPn flash converged not very well, reduce dx
                    sol.vaporfra = lowerl + dx;
                }
                else
                {
                    dx = coeff2 * fonc / deri;
                    sol.vaporfra = sol.vaporfra - dx;
                }

                scalar absdx = dx > 0 ? dx : -dx;

                if (absdx < tolerance)
                {
                    // niter = nitmax;
                    break;
                }
                niter++;
            }
            // if (niter < nitmax)
            // {
            //     ierr = 2;
            // }
        }

        forAll(this->X_, i)
        {
            sol.X_liq()[i] = this->X_[i] / (1.0 + KM1[i] * sol.vaporfra);
            sol.X_gas()[i] = sol.equalconstant()[i] * sol.X_liq()[i];
        }

        if (ierr == 0)
        {
            ierr = this->noramlize(sol.X_gas());
            ierr = this->noramlize(sol.X_liq());
        }
    }
    return ierr;
}

template <class ThermoMixture>
inline autoPtr<typename Foam::VLE<ThermoMixture>::solution> Foam::VLE<ThermoMixture>::TPn_flash // PY:with a course initialization
(
    scalar p, //input
    scalar T  //input
) const
{
    //Info<<"VLE!\n";
    //autoPtr<solution> pS(&temp_sol);
    autoPtr<solution> pS(new solution());
    //if(this->changed==0&&temp_p==p&&temp_T==T)
    //return pS;
    solution& Sol = *pS;
    const int maxiter = 1000;
    int ierr = 0;
    //int iaccelation = 1;
    scalar tol_fug = 1.0e-13; //Convergence criteria
    scalar tol_fug2 = 1.0;
    scalar coeff_tpn = 1;
    scalar epsilon1 = 1.0e+10;

    //scalar gibbsenergy = 0.0;
    autoPtr<scalarList> fugcoef_liq, fugcoef_gas;              //, fugcoef_mix[num_sp] = {0.0};
    scalarList X_liq0(this->N_, Zero), X_gas0(this->N_, Zero); //, fugcoef_mix[num_sp] = {0.0};
    scalarList temp(this->N_, Zero);
    int newN = 0;
    labelList mapi(this->N_);
    forAll(this->X_, i)
    {
        if (this->X_[i] > 1e-10)
        {
            mapi[newN++] = i;
        }
    }
    if (newN == 0)
    {
        FatalErrorInFunction << "newN=0" << exit(FatalError);
    }

    Sol.X_liq.reset(new scalarList(this->N_));
    Sol.X_gas.reset(new scalarList(this->N_));
    forAll(this->X_, i)
    {
        Sol.X_liq()[i] = this->X_[i]; //use this in PT_x flash, cancle this in hn flash
        Sol.X_gas()[i] = this->X_[i]; //use this in PT_x flash, cancle this in hn flash
        X_liq0[i] = Sol.X_liq()[i];
        X_gas0[i] = Sol.X_gas()[i];
    }

    fugcoef_liq.reset(this->fugacityCoefficient(p, T, 0).ptr());
    fugcoef_gas.reset(this->fugacityCoefficient(p, T, 1).ptr());

    //Info <<"fugcoef_liq="<<fugcoef_liq()<<endl;
    //Info <<"fugcoef_gas="<<fugcoef_gas()<<endl;
    Sol.equalconstant.reset(this->wilsoneq(p, T).ptr());
    //Sol.equalconstant.reset(this->myinitial(p, T).ptr());
    //Info << "equalconstant= " << Sol.equalconstant() << endl;

    Sol.equalconstant()[0] = 300.98332688687367;
    Sol.equalconstant()[1] = 0.012591909335566309;
    Sol.equalconstant()[2] = 3302.6955988917484;

    /*
        Sol.equalconstant()[0] = 135.508;
        Sol.equalconstant()[1] = 0.0560828;
        Sol.equalconstant()[2] = 648.055;
    */

    //Sol.equalconstant()[0] = 4000.006;
    //Sol.equalconstant()[1] = 0.240333;
    //Sol.equalconstant()[2] = 1;

    if (newN == 1)
    {
        //Info << "gas:" << fugcoef_gas()[0] << " liq:" << fugcoef_liq()[0] << endl;
        if (fugcoef_gas()[mapi[0]] > fugcoef_liq()[mapi[0]])
            Sol.vaporfra = 0;
        else
            Sol.vaporfra = 1;
        Sol.equalconstant()[mapi[0]] = 1;
        return pS;
    }
    scalarList re_iter1(this->N_);
    forAll(re_iter1, i)
        re_iter1[i] = 1;
    //scalar gibbs_iter = gibbsenergy;
    //scalar vaporfra_iter1 = 0.0;
    //scalar vaporfra_iter0 = 0.0;

    for (label iter = 0; iter < maxiter; iter++)
    {

        //vaporfra_iter0 = Sol.vaporfra;
        ierr = this->solve_vaporf(p, T, Sol);
        //Info << "vaporfra= " << Sol.vaporfra << endl;
        // only for Hn flash at high temperature and pressure conditions

        if (ierr == 4 || ierr == 3)
        {
            //WarningInFunction << "VLE "
            //                  << "p=" << p << " T=" << T << " X=" << this->X_ << endl;
            forAll(this->X_, i)
            {
                Sol.X_liq()[i] = X_liq0[i];
                Sol.X_gas()[i] = X_gas0[i];
            }
            if (Sol.vaporfra > 1.0)
            {
                Sol.vaporfra = 0.9999999;
            }
            else if (Sol.vaporfra < 0.0)
            {
                Sol.vaporfra = 1.0e-07;
            }
            break;
        }
        //Info << "vaporfra4= " << pS().vaporfra << endl;

        /* if (Sol.vaporfra > 0.99999999) //!!!!&& flag_tp > 0
        {
            Sol.vaporfra = 0.99999999;
            forAll(this->X_, i)
            {
                Sol.X_liq()[i] = X_liq0[i];
                Sol.X_gas()[i] = X_gas0[i];
            }
            break;
        }*/
        /*if (Sol.vaporfra < 1.0e-08) //!!!!&& flag_tp > 0
        {
            Sol.vaporfra = 1.0e-08;
            forAll(this->X_, i)
            {
                Sol.X_liq()[i] = X_liq0[i];
                Sol.X_gas()[i] = X_gas0[i];
            }
            break;
        }*/
        //Info << "vaporfra3= " << pS().vaporfra << endl;
        //vaporfra_iter1 = Sol.vaporfra;
        if ((ierr == 10) && iter > 2)
        {
            FatalErrorInFunction << "VLE!! "
                << "p=" << p << " T=" << T << " X=" << this->X_ << exit(FatalError);
            //this->TPD(p, T, comp, comp_liq, comp_gas, tpdvalue);

            fugcoef_liq.reset(this->fugacityCoefficient(p, T, 0, &(Sol.X_liq())).ptr());
            fugcoef_gas.reset(this->fugacityCoefficient(p, T, 1, &(Sol.X_gas())).ptr());

            forAll(this->X_, i)
            {
                Sol.equalconstant()[i] = fugcoef_liq()[i] / fugcoef_gas()[i];
            }
            ierr = this->solve_vaporf(p, T, Sol);

            // Here, for tp diagram, please reconsider it
            forAll(this->X_, i)
            {
                Sol.X_liq()[i] = X_liq0[i];
                Sol.X_gas()[i] = X_gas0[i];
            }

            if (Sol.vaporfra > 1.0)
            {
                Sol.vaporfra = 0.9999999;
            }
            if (Sol.vaporfra < 0.0)
            {
                Sol.vaporfra = 1.0e-07;
            }
            break;
        }

        fugcoef_liq.reset(this->fugacityCoefficient(p, T, 0, &(Sol.X_liq())).ptr()); //0.79210807054894983,0.10280737903377266
        fugcoef_gas.reset(this->fugacityCoefficient(p, T, 1, &(Sol.X_gas())).ptr());    //0.64085593096880322,0.2817704554420456
        /*
        Info <<"--------------------------"<<endl;
        Info <<"fugcoef_liq="<<fugcoef_liq()<<endl;
        Info <<"fugcoef_gas="<<fugcoef_gas()<<endl;
        Info <<"X_liq="<<Sol.X_liq()<<endl;
        Info <<"X_gas="<<Sol.X_gas()<<endl;
        Info <<"--------------------------"<<endl;
        */
        //Info << "vaporfra2= " << pS().vaporfra << endl;
        scalar errr = 0.0;
        scalarList re_iter0(this->N_);
        //scalarList save_k(this->N_);
        //scalar rate=0.1,temp_vf;
        //rate=1;

        /*
        forAll(this->X_, i)
        {
            re_iter0[i] = re_iter1[i];
            fugcoef_gas()[i] *= Sol.X_gas()[i] * p; //14702443.333238741,16374.974213553087
            fugcoef_liq()[i] *= Sol.X_liq()[i] * p; //12945283.239314515,684406.75748597819
            errr = errr + sqr(fugcoef_liq()[i] / fugcoef_gas()[i] - 1.0);
            //save_k[i]=Sol.equalconstant()[i];
            Sol.equalconstant()[i] *= (fugcoef_liq()[i] / fugcoef_gas()[i]) * coeff_tpn; //1.2360158223884925,0.36486216722930281

            re_iter1[i] = fugcoef_liq()[i] / fugcoef_gas()[i];
        }
*/
        for (int i = 0;i < newN; i++)
        {
            re_iter0[mapi[i]] = re_iter1[mapi[i]];
            fugcoef_gas()[mapi[i]] *= Sol.X_gas()[mapi[i]] * p;
            fugcoef_liq()[mapi[i]] *= Sol.X_liq()[mapi[i]] * p;
            errr = errr + sqr(fugcoef_liq()[mapi[i]] / fugcoef_gas()[mapi[i]] - 1.0);
            //save_k[i]=Sol.equalconstant()[i];
            Sol.equalconstant()[mapi[i]] *= (fugcoef_liq()[mapi[i]] / fugcoef_gas()[mapi[i]]) * coeff_tpn;
            //Sol.equalconstant()[mapi[i]] *=pow( (fugcoef_liq()[mapi[i]] / fugcoef_gas()[mapi[i]]) ,0.1);
            re_iter1[mapi[i]] = fugcoef_liq()[mapi[i]] / fugcoef_gas()[mapi[i]];
        }



        /*
        ierr=this->solve_vaporf(p, T, Sol);
        temp_vf=Sol.vaporfra;
        while (temp_vf<0||temp_vf>1||ierr == 4 || ierr == 3)
        {
            rate/=2;
            forAll(this->X_, i)
            {
                Sol.equalconstant()[i] =save_k[i]*max( min(1+rate,(fugcoef_liq()[i] / fugcoef_gas()[i]) * coeff_tpn),1/(1+rate));
            }
            ierr=this->solve_vaporf(p, T, Sol);
            temp_vf=Sol.vaporfra;
        }
        */
        epsilon1 = errr;

        scalar epsilon2 = 1.0e+10;
        if (errr >= 1.0e+05 && iter > 10)
        {
            fugcoef_liq.reset(this->fugacityCoefficient(p, T, 0, &(Sol.X_liq())).ptr());
            fugcoef_gas.reset(this->fugacityCoefficient(p, T, 1, &(Sol.X_gas())).ptr());
            /*
                        forAll(this->X_, i)
                        {
                            Sol.equalconstant()[i] = fugcoef_liq()[i] / fugcoef_gas()[i];
                        }
                        */
            for (int i = 0;i < newN; i++)
            {
                Sol.equalconstant()[mapi[i]] = fugcoef_liq()[mapi[i]] / fugcoef_gas()[mapi[i]];
            }

            ierr = this->solve_vaporf(p, T, Sol);
            scalar errr2 = 0.0;

            /*forAll(this->X_, i)
            {
                re_iter0[i] = re_iter1[i];
                fugcoef_gas()[i] *= Sol.X_gas()[i] * p;
                fugcoef_liq()[i] *= Sol.X_liq()[i] * p;
                errr2 = errr2 + sqr(fugcoef_liq()[i] / fugcoef_gas()[i] - 1.0);
                Sol.equalconstant()[i] *= (fugcoef_liq()[i] / fugcoef_gas()[i]);
                re_iter1[i] = fugcoef_liq()[i] / fugcoef_gas()[i];
            }*/
            for (int i = 0;i < newN; i++)
            {
                re_iter0[mapi[i]] = re_iter1[mapi[i]];
                fugcoef_gas()[mapi[i]] *= Sol.X_gas()[mapi[i]] * p;
                fugcoef_liq()[mapi[i]] *= Sol.X_liq()[mapi[i]] * p;
                errr2 = errr2 + sqr(fugcoef_liq()[mapi[i]] / fugcoef_gas()[mapi[i]] - 1.0);
                Sol.equalconstant()[mapi[i]] *= (fugcoef_liq()[mapi[i]] / fugcoef_gas()[mapi[i]]);
                re_iter1[mapi[i]] = fugcoef_liq()[mapi[i]] / fugcoef_gas()[mapi[i]];
            }
            epsilon2 = errr2;
        }

        if (epsilon1 < tol_fug || epsilon2 < tol_fug2)
        {
            if (Sol.vaporfra > 0.99999999)
            {
                Sol.vaporfra = 0.99999999;
                forAll(this->X_, i)
                {
                    Sol.X_liq()[i] = X_liq0[i];
                    Sol.X_gas()[i] = X_gas0[i];
                }
            }
            else if (Sol.vaporfra < 1.0e-010) //!!!!&& flag_tp > 0
            {
                Sol.vaporfra = 1.0e-010;
                forAll(this->X_, i)
                {
                    Sol.X_liq()[i] = X_liq0[i];
                    Sol.X_gas()[i] = X_gas0[i];
                }
            }
            break;
        }
        /*
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
        */
        //Info << "vaporfra1= " << pS().vaporfra << endl;
    }
    //printf("TPflashEnd= %+5.8f %+5.8f\n", comp_liq[0], comp_gas[0]);
    //FatalErrorInFunction
    //<< exit(FatalError);
    //Info << "vaporfra= " << pS().vaporfra << endl;
    //this->changed=0;
    //temp_p=p;
    //temp_T=T;
    return pS;
}

template <class ThermoMixture>
inline autoPtr<typename Foam::VLE<ThermoMixture>::solution> Foam::VLE<ThermoMixture>::TPn_flash_New // PY:with a course initialization
(
    scalar p, //input
    scalar T  //input
) const
{

    autoPtr<solution> pS(new solution());
    //if(this->changed==0&&temp_p==p&&temp_T==T)
    //return pS;
    solution& Sol = *pS;
    const int maxiter = 1000;
    int ierr = 0;
    //int iaccelation = 1;
    scalar tol_fug = 1.0e-13; //Convergence criteria
    scalar tol_fug2 = 1.0;
    scalar coeff_tpn = 1;
    scalar epsilon1 = 1.0e+10;

    //scalar gibbsenergy = 0.0;
    autoPtr<scalarList> fugcoef_liq, fugcoef_gas;              //, fugcoef_mix[num_sp] = {0.0};
    scalarList X_liq0(this->N_, Zero), X_gas0(this->N_, Zero); //, fugcoef_mix[num_sp] = {0.0};
    scalarList temp(this->N_, Zero);
    int newN = 0;
    labelList mapi(this->N_);
    forAll(this->X_, i)
    {
        if (this->X_[i] > 1e-10)
        {
            mapi[newN++] = i;
        }
    }
    if (newN == 0)
    {
        FatalErrorInFunction << "newN=0" << exit(FatalError);
    }

    Sol.X_liq.reset(new scalarList(this->N_));
    Sol.X_gas.reset(new scalarList(this->N_));
    forAll(this->X_, i)
    {
        Sol.X_liq()[i] = this->X_[i]; //use this in PT_x flash, cancle this in hn flash
        Sol.X_gas()[i] = this->X_[i]; //use this in PT_x flash, cancle this in hn flash
        X_liq0[i] = Sol.X_liq()[i];
        X_gas0[i] = Sol.X_gas()[i];
    }

    bool stable, isVapor;
    Foam::autoPtr<Foam::scalarList> K_init;
    std::tie(stable, K_init, isVapor) = solveTPD_BFGS(p, T);

    if (stable != true)
        Sol.equalconstant.reset(K_init.ptr());
    else
    {
        if (isVapor == true)
        {
            Sol.vaporfra = 1.0;
        }
        else
        {
            Sol.vaporfra = 0.0;
        }
        Sol.equalconstant.reset(new scalarList(this->N_));
        return pS;
    }

    fugcoef_liq.reset(this->fugacityCoefficient(p, T, 0).ptr());
    fugcoef_gas.reset(this->fugacityCoefficient(p, T, 1).ptr());

    if (newN == 1)
    {
        //Info << "gas:" << fugcoef_gas()[0] << " liq:" << fugcoef_liq()[0] << endl;
        if (fugcoef_gas()[mapi[0]] > fugcoef_liq()[mapi[0]])
            Sol.vaporfra = 0;
        else
            Sol.vaporfra = 1;
        Sol.equalconstant()[mapi[0]] = 1;
        return pS;
    }
    scalarList re_iter1(this->N_);
    forAll(re_iter1, i)
        re_iter1[i] = 1;
    //scalar gibbs_iter = gibbsenergy;
    //scalar vaporfra_iter1 = 0.0;
    //scalar vaporfra_iter0 = 0.0;

    for (label iter = 0; iter < maxiter; iter++)
    {

        //vaporfra_iter0 = Sol.vaporfra;
        ierr = this->solve_vaporf(p, T, Sol);
        //Info << "vaporfra= " << Sol.vaporfra << endl;
        // only for Hn flash at high temperature and pressure conditions


        if (ierr == 4 || ierr == 3)
        {
            FatalErrorInFunction << "VLE " << "p=" << p << " T=" << T << " X=" << this->X_ << endl;
            forAll(this->X_, i)
            {
                Sol.X_liq()[i] = X_liq0[i];
                Sol.X_gas()[i] = X_gas0[i];
            }
            if (Sol.vaporfra > 1.0)
            {
                Sol.vaporfra = 0.9999999;
            }
            else if (Sol.vaporfra < 0.0)
            {
                Sol.vaporfra = 1.0e-07;
            }
            break;
        }

        if ((ierr == 10) && iter > 2)
        {
            FatalErrorInFunction << "VLE!! "
                << "p=" << p << " T=" << T << " X=" << this->X_ << exit(FatalError);
            //this->TPD(p, T, comp, comp_liq, comp_gas, tpdvalue);

            fugcoef_liq.reset(this->fugacityCoefficient(p, T, 0, &(Sol.X_liq())).ptr());
            fugcoef_gas.reset(this->fugacityCoefficient(p, T, 1, &(Sol.X_gas())).ptr());

            forAll(this->X_, i)
            {
                Sol.equalconstant()[i] = fugcoef_liq()[i] / fugcoef_gas()[i];
            }
            ierr = this->solve_vaporf(p, T, Sol);

            // Here, for tp diagram, please reconsider it
            forAll(this->X_, i)
            {
                Sol.X_liq()[i] = X_liq0[i];
                Sol.X_gas()[i] = X_gas0[i];
            }

            if (Sol.vaporfra > 1.0)
            {
                Sol.vaporfra = 0.9999999;
            }
            if (Sol.vaporfra < 0.0)
            {
                Sol.vaporfra = 1.0e-07;
            }
            break;
        }

        fugcoef_liq.reset(this->fugacityCoefficient(p, T, 0, &(Sol.X_liq())).ptr());
        fugcoef_gas.reset(this->fugacityCoefficient(p, T, 1, &(Sol.X_gas())).ptr());

        scalar errr = 0.0;
        scalarList re_iter0(this->N_);

        for (int i = 0;i < newN; i++)
        {
            re_iter0[mapi[i]] = re_iter1[mapi[i]];
            fugcoef_gas()[mapi[i]] *= Sol.X_gas()[mapi[i]] * p;
            fugcoef_liq()[mapi[i]] *= Sol.X_liq()[mapi[i]] * p;
            errr = errr + sqr(fugcoef_liq()[mapi[i]] / fugcoef_gas()[mapi[i]] - 1.0);
            //save_k[i]=Sol.equalconstant()[i];
            temp[i] = Sol.equalconstant()[mapi[i]];
            //Sol.equalconstant()[mapi[i]] *= (fugcoef_liq()[mapi[i]] / fugcoef_gas()[mapi[i]]) * coeff_tpn;
            //Sol.equalconstant()[mapi[i]] *=pow( (fugcoef_liq()[mapi[i]] / fugcoef_gas()[mapi[i]]) ,0.1);
            re_iter1[mapi[i]] = fugcoef_liq()[mapi[i]] / fugcoef_gas()[mapi[i]];
        }
        double powi = 2.0;
        for (ierr = this->solve_vaporf(p, T, Sol);Sol.vaporfra > 1.0;ierr = this->solve_vaporf(p, T, Sol))
        {
            for (int i = 0;i < newN; i++)
            {
                Sol.equalconstant()[mapi[i]] = temp[i] / (fugcoef_liq()[mapi[i]] / fugcoef_gas()[mapi[i]]);
                temp[i] /= (fugcoef_liq()[mapi[i]] / fugcoef_gas()[mapi[i]]);
            }
        }
        do
        {
            ;
            powi /= 2.0;
            for (int i = 0;i < newN; i++)
            {
                Sol.equalconstant()[mapi[i]] = temp[i] * pow((fugcoef_liq()[mapi[i]] / fugcoef_gas()[mapi[i]]), powi);
            }
            ierr = this->solve_vaporf(p, T, Sol);
        } while (Sol.vaporfra > 1.0);



        epsilon1 = errr;

        scalar epsilon2 = 1.0e+10;
        if (errr >= 1.0e+05 && iter > 10)
        {
            fugcoef_liq.reset(this->fugacityCoefficient(p, T, 0, &(Sol.X_liq())).ptr());
            fugcoef_gas.reset(this->fugacityCoefficient(p, T, 1, &(Sol.X_gas())).ptr());
            /*
                        forAll(this->X_, i)
                        {
                            Sol.equalconstant()[i] = fugcoef_liq()[i] / fugcoef_gas()[i];
                        }
                        */
            for (int i = 0;i < newN; i++)
            {
                Sol.equalconstant()[mapi[i]] = fugcoef_liq()[mapi[i]] / fugcoef_gas()[mapi[i]];
            }

            ierr = this->solve_vaporf(p, T, Sol);
            scalar errr2 = 0.0;

            /*forAll(this->X_, i)
            {
                re_iter0[i] = re_iter1[i];
                fugcoef_gas()[i] *= Sol.X_gas()[i] * p;
                fugcoef_liq()[i] *= Sol.X_liq()[i] * p;
                errr2 = errr2 + sqr(fugcoef_liq()[i] / fugcoef_gas()[i] - 1.0);
                Sol.equalconstant()[i] *= (fugcoef_liq()[i] / fugcoef_gas()[i]);
                re_iter1[i] = fugcoef_liq()[i] / fugcoef_gas()[i];
            }*/
            for (int i = 0;i < newN; i++)
            {
                re_iter0[mapi[i]] = re_iter1[mapi[i]];
                fugcoef_gas()[mapi[i]] *= Sol.X_gas()[mapi[i]] * p;
                fugcoef_liq()[mapi[i]] *= Sol.X_liq()[mapi[i]] * p;
                errr2 = errr2 + sqr(fugcoef_liq()[mapi[i]] / fugcoef_gas()[mapi[i]] - 1.0);
                Sol.equalconstant()[mapi[i]] *= (fugcoef_liq()[mapi[i]] / fugcoef_gas()[mapi[i]]);
                re_iter1[mapi[i]] = fugcoef_liq()[mapi[i]] / fugcoef_gas()[mapi[i]];
            }
            epsilon2 = errr2;
        }

        if (epsilon1 < tol_fug || epsilon2 < tol_fug2)
        {
            if (Sol.vaporfra > 0.99999999)
            {
                Sol.vaporfra = 0.99999999;
                forAll(this->X_, i)
                {
                    Sol.X_liq()[i] = X_liq0[i];
                    Sol.X_gas()[i] = X_gas0[i];
                }
            }
            else if (Sol.vaporfra < 1.0e-010) //!!!!&& flag_tp > 0
            {
                Sol.vaporfra = 1.0e-010;
                forAll(this->X_, i)
                {
                    Sol.X_liq()[i] = X_liq0[i];
                    Sol.X_gas()[i] = X_gas0[i];
                }
            }
            break;
        }
    }

    return pS;
}


template <class ThermoMixture>
inline autoPtr<typename Foam::VLE<ThermoMixture>::solution> Foam::VLE<ThermoMixture>::TPn_flash_New2 // PY:with a course initialization
(
    scalar p, //input
    scalar T,  //input
    bool stable,
    bool isVapor,
    Foam::scalarList& K_init
) const
{

    autoPtr<solution> pS(new solution());
    //if(this->changed==0&&temp_p==p&&temp_T==T)
    //return pS;
    solution& Sol = *pS;
    const int maxiter = 1000;
    int ierr = 0;
    //int iaccelation = 1;
    scalar tol_fug = 1.0e-13; //Convergence criteria
    scalar tol_fug2 = 1.0;
    scalar coeff_tpn = 1;
    scalar epsilon1 = 1.0e+10;

    //scalar gibbsenergy = 0.0;
    autoPtr<scalarList> fugcoef_liq, fugcoef_gas;              //, fugcoef_mix[num_sp] = {0.0};
    scalarList X_liq0(this->N_, Zero), X_gas0(this->N_, Zero); //, fugcoef_mix[num_sp] = {0.0};
    scalarList temp(this->N_, Zero);
    int newN = 0;
    labelList mapi(this->N_);
    forAll(this->X_, i)
    {
        if (this->X_[i] > 1e-10)
        {
            mapi[newN++] = i;
        }
    }
    if (newN == 0)
    {
        FatalErrorInFunction << "newN=0" << exit(FatalError);
    }

    Sol.X_liq.reset(new scalarList(this->N_));
    Sol.X_gas.reset(new scalarList(this->N_));
    forAll(this->X_, i)
    {
        Sol.X_liq()[i] = this->X_[i]; //use this in PT_x flash, cancle this in hn flash
        Sol.X_gas()[i] = this->X_[i]; //use this in PT_x flash, cancle this in hn flash
        X_liq0[i] = Sol.X_liq()[i];
        X_gas0[i] = Sol.X_gas()[i];
    }



    if (stable != true)
    {
        Sol.equalconstant.reset(new scalarList(this->N_));
        forAll(this->X_, i)
        {
            Sol.equalconstant()[i] = K_init[i];
        }
    }
    else
    {
        if (isVapor == true)
        {
            Sol.vaporfra = 1.0;
        }
        else
        {
            Sol.vaporfra = 0.0;
        }
        Sol.equalconstant.reset(new scalarList(this->N_));
        return pS;
    }

    fugcoef_liq.reset(this->fugacityCoefficient(p, T, 0).ptr());
    fugcoef_gas.reset(this->fugacityCoefficient(p, T, 1).ptr());

    if (newN == 1)
    {
        //Info << "gas:" << fugcoef_gas()[0] << " liq:" << fugcoef_liq()[0] << endl;
        if (fugcoef_gas()[mapi[0]] > fugcoef_liq()[mapi[0]])
            Sol.vaporfra = 0;
        else
            Sol.vaporfra = 1;
        Sol.equalconstant()[mapi[0]] = 1;
        return pS;
    }
    scalarList re_iter1(this->N_);
    forAll(re_iter1, i)
        re_iter1[i] = 1;
    //scalar gibbs_iter = gibbsenergy;
    //scalar vaporfra_iter1 = 0.0;
    //scalar vaporfra_iter0 = 0.0;

    for (label iter = 0; iter < maxiter; iter++)
    {

        //vaporfra_iter0 = Sol.vaporfra;
        ierr = this->solve_vaporf(p, T, Sol);
        //Info << "vaporfra= " << Sol.vaporfra << endl;
        // only for Hn flash at high temperature and pressure conditions


        if (ierr == 4 || ierr == 3)
        {
            FatalErrorInFunction << "VLE " << "p=" << p << " T=" << T << " X=" << this->X_ << endl;
            forAll(this->X_, i)
            {
                Sol.X_liq()[i] = X_liq0[i];
                Sol.X_gas()[i] = X_gas0[i];
            }
            if (Sol.vaporfra > 1.0)
            {
                Sol.vaporfra = 0.9999999;
            }
            else if (Sol.vaporfra < 0.0)
            {
                Sol.vaporfra = 1.0e-07;
            }
            break;
        }

        if ((ierr == 10) && iter > 2)
        {
            FatalErrorInFunction << "VLE!! "
                << "p=" << p << " T=" << T << " X=" << this->X_ << exit(FatalError);
            //this->TPD(p, T, comp, comp_liq, comp_gas, tpdvalue);

            fugcoef_liq.reset(this->fugacityCoefficient(p, T, 0, &(Sol.X_liq())).ptr());
            fugcoef_gas.reset(this->fugacityCoefficient(p, T, 1, &(Sol.X_gas())).ptr());

            forAll(this->X_, i)
            {
                Sol.equalconstant()[i] = fugcoef_liq()[i] / fugcoef_gas()[i];
            }
            ierr = this->solve_vaporf(p, T, Sol);

            // Here, for tp diagram, please reconsider it
            forAll(this->X_, i)
            {
                Sol.X_liq()[i] = X_liq0[i];
                Sol.X_gas()[i] = X_gas0[i];
            }

            if (Sol.vaporfra > 1.0)
            {
                Sol.vaporfra = 0.9999999;
            }
            if (Sol.vaporfra < 0.0)
            {
                Sol.vaporfra = 1.0e-07;
            }
            break;
        }

        fugcoef_liq.reset(this->fugacityCoefficient(p, T, 0, &(Sol.X_liq())).ptr());
        fugcoef_gas.reset(this->fugacityCoefficient(p, T, 1, &(Sol.X_gas())).ptr());

        scalar errr = 0.0;
        scalarList re_iter0(this->N_);

        for (int i = 0;i < newN; i++)
        {
            re_iter0[mapi[i]] = re_iter1[mapi[i]];
            fugcoef_gas()[mapi[i]] *= Sol.X_gas()[mapi[i]] * p;
            fugcoef_liq()[mapi[i]] *= Sol.X_liq()[mapi[i]] * p;
            errr = errr + sqr(fugcoef_liq()[mapi[i]] / fugcoef_gas()[mapi[i]] - 1.0);
            //save_k[i]=Sol.equalconstant()[i];
            temp[i] = Sol.equalconstant()[mapi[i]];
            //Sol.equalconstant()[mapi[i]] *= (fugcoef_liq()[mapi[i]] / fugcoef_gas()[mapi[i]]) * coeff_tpn;
            //Sol.equalconstant()[mapi[i]] *=pow( (fugcoef_liq()[mapi[i]] / fugcoef_gas()[mapi[i]]) ,0.1);
            re_iter1[mapi[i]] = fugcoef_liq()[mapi[i]] / fugcoef_gas()[mapi[i]];
        }
        double powi = 2.0;
        for (ierr = this->solve_vaporf(p, T, Sol);Sol.vaporfra > 1.0;ierr = this->solve_vaporf(p, T, Sol))
        {
            for (int i = 0;i < newN; i++)
            {
                Sol.equalconstant()[mapi[i]] = temp[i] / (fugcoef_liq()[mapi[i]] / fugcoef_gas()[mapi[i]]);
                temp[i] /= (fugcoef_liq()[mapi[i]] / fugcoef_gas()[mapi[i]]);
            }
        }
        do
        {
            ;
            powi /= 2.0;
            for (int i = 0;i < newN; i++)
            {
                Sol.equalconstant()[mapi[i]] = temp[i] * pow((fugcoef_liq()[mapi[i]] / fugcoef_gas()[mapi[i]]), powi);
            }
            ierr = this->solve_vaporf(p, T, Sol);
        } while (Sol.vaporfra > 1.0);



        epsilon1 = errr;

        scalar epsilon2 = 1.0e+10;
        if (errr >= 1.0e+05 && iter > 10)
        {
            fugcoef_liq.reset(this->fugacityCoefficient(p, T, 0, &(Sol.X_liq())).ptr());
            fugcoef_gas.reset(this->fugacityCoefficient(p, T, 1, &(Sol.X_gas())).ptr());
            /*
                        forAll(this->X_, i)
                        {
                            Sol.equalconstant()[i] = fugcoef_liq()[i] / fugcoef_gas()[i];
                        }
                        */
            for (int i = 0;i < newN; i++)
            {
                Sol.equalconstant()[mapi[i]] = fugcoef_liq()[mapi[i]] / fugcoef_gas()[mapi[i]];
            }

            ierr = this->solve_vaporf(p, T, Sol);
            scalar errr2 = 0.0;

            /*forAll(this->X_, i)
            {
                re_iter0[i] = re_iter1[i];
                fugcoef_gas()[i] *= Sol.X_gas()[i] * p;
                fugcoef_liq()[i] *= Sol.X_liq()[i] * p;
                errr2 = errr2 + sqr(fugcoef_liq()[i] / fugcoef_gas()[i] - 1.0);
                Sol.equalconstant()[i] *= (fugcoef_liq()[i] / fugcoef_gas()[i]);
                re_iter1[i] = fugcoef_liq()[i] / fugcoef_gas()[i];
            }*/
            for (int i = 0;i < newN; i++)
            {
                re_iter0[mapi[i]] = re_iter1[mapi[i]];
                fugcoef_gas()[mapi[i]] *= Sol.X_gas()[mapi[i]] * p;
                fugcoef_liq()[mapi[i]] *= Sol.X_liq()[mapi[i]] * p;
                errr2 = errr2 + sqr(fugcoef_liq()[mapi[i]] / fugcoef_gas()[mapi[i]] - 1.0);
                Sol.equalconstant()[mapi[i]] *= (fugcoef_liq()[mapi[i]] / fugcoef_gas()[mapi[i]]);
                re_iter1[mapi[i]] = fugcoef_liq()[mapi[i]] / fugcoef_gas()[mapi[i]];
            }
            epsilon2 = errr2;
        }

        if (epsilon1 < tol_fug || epsilon2 < tol_fug2)
        {
            if (Sol.vaporfra > 0.99999999)
            {
                Sol.vaporfra = 0.99999999;
                forAll(this->X_, i)
                {
                    Sol.X_liq()[i] = X_liq0[i];
                    Sol.X_gas()[i] = X_gas0[i];
                }
            }
            else if (Sol.vaporfra < 1.0e-010) //!!!!&& flag_tp > 0
            {
                Sol.vaporfra = 1.0e-010;
                forAll(this->X_, i)
                {
                    Sol.X_liq()[i] = X_liq0[i];
                    Sol.X_gas()[i] = X_gas0[i];
                }
            }
            break;
        }
    }

    return pS;
}



template <class ThermoMixture>
inline void Foam::VLE<ThermoMixture>::TPn_validation(
    scalar p, //input
    scalar T  //input
) const
{
    autoPtr<solution> sol(TPn_flash(p, T));

    label pfg = 1, pfl = 0;
    scalarList newX(this->N_);
    scalarList fg(this->N_), fl(this->N_);
    autoPtr<scalarList> fg_ref = this->fugacityCoefficient(p, T, &pfg, &(sol().X_gas()));
    autoPtr<scalarList> fl_ref = this->fugacityCoefficient(p, T, &pfl, &(sol().X_liq()));
    Info << this->N_ << endl;
    forAll(newX, i)
    {
        newX[i] = sol().X_gas()[i] * sol().vaporfra + (1 - sol().vaporfra) * sol().X_liq()[i];
        fg[i] = sol().X_gas()[i] * fg_ref()[i];
        fl[i] = sol().X_liq()[i] * fl_ref()[i];
    }
    Info << "c\n";
    Info << "fugacity_g=" << fg << endl;

    Info << "fugacity_l=" << fl << endl;
    Info << "vaporfra=" << sol().vaporfra << endl;
    Info << "z=" << this->X_ << endl;

    Info << "z_new=" << newX << endl;
}


template <class ThermoMixture>
inline autoPtr<typename Foam::VLE<ThermoMixture>::solution> Foam::VLE<ThermoMixture>::TPn_flash_New_TPD // PY:with a course initialization
(
    scalar p, //input
    scalar T  //input
) const
{
    bool stable, isVapor;
    scalar gmin = DBL_MAX;
    const scalar delta = 1e-8;
    autoPtr<List<scalarList>> p_K_init;
    autoPtr<solution> sol_min;
    autoPtr<scalarList> p_TPD_star;
    tie(stable, p_K_init, isVapor, p_TPD_star) = solveTPD_BFGS_v2(p, T);
    labelList order;
    sortedOrder(p_TPD_star(), order);
    for (int i = 0;i < p_TPD_star->size();i++)
    {
        //Info << p_TPD_star()[order[i]] << " trying2:" << i << endl;
        if (i == 0 || (p_TPD_star()[order[i]] < delta && fabs(p_TPD_star()[order[i]] - p_TPD_star()[order[i - 1]]) > delta))
        {
            //Info << "tryed2:" << i << endl;
            //Info << p_TPD_star()[order[i]] << "," << p_TPD_star()[order[i]] - p_TPD_star()[order[i == 0 ? i : i - 1]] << endl;
            autoPtr<solution> sol = TPn_flash_New2(p, T, stable, isVapor, p_K_init()[order[i]]);
            scalar ggas = this->G_Mole(p, T, &(sol().X_gas()));
            scalar gliq = this->G_Mole(p, T, &(sol().X_liq()));
            scalar g = sol().vaporfra * ggas + (1 - sol().vaporfra) * gliq;
            //Info << sol().vaporfra << " g:" << g << endl;
            //Info << ggas << "," << gliq << "," << sol().vaporfra << "," << sol().vaporfra * ggas + (1 - sol().vaporfra) * gliq << endl;
            if (g < gmin)
            {
                gmin = g;
                sol_min.reset(sol);
            }
        }
    }

    //Info << "trying" << endl;
    //sol_min();
    //Info << "tryed" << endl;
    //autoPtr<solution> sol(TPn_flash(p, T));
    return sol_min;
}

template <class ThermoMixture>
inline tuple<bool, Foam::autoPtr<scalarList>, bool>  Foam::VLE<ThermoMixture>::solveTPD_BFGS
(
    scalar p, //input
    scalar T  //input
) const
{
    const scalar be_negative = -1.E-8;
    const scalar be_pure_substance = 1. - 1.E-8;
    const scalar be_very_large = 1.E10;

    scalar pow_v[4] = { 1.0,1.0 / 3,-1.0,-1.0 / 3.0 };
    scalar pow_l[4] = { -1.0,-1.0 / 3,1.0,1.0 / 3.0 };

    bool stable = true;
    bool  isVapor;
    scalar TPDmin = 10.0;
    autoPtr<scalarList> p_K_init_min(new scalarList(this->X_.size()));
    scalarList& K_init_min = p_K_init_min();
    int n = this->X_.size();
    int ntrial = 4 + n;
    scalar tol = 1.E-9;
    int iter_max = 1000;
    scalarList y_sp[ntrial];
    scalar z = this->Z_gibbs(p, T, &this->X_);
    scalarList TPD_star(ntrial), TPD(ntrial), mask_stationary(ntrial);
    autoPtr<scalarList> lnphi(this->ln_fugacityCoefficient(p, T, z, &this->X_));
    scalarList d(this->X_.size()), K(this->X_.size());
    for (int i = 0;i < this->X_.size();i++)
    {
        d[i] = log(this->X_[i]) + lnphi()[i];
        K[i] = (*this)[i].Pc_ / p * exp(5.373 * (1.0 + (*this)[i].omega_) * (1.0 - (*this)[i].Tc_ / T));
    }
    scalar B = this->B(p, T, &this->X_);
    scalarList Y_init[ntrial];
    for (int i = 0;i < ntrial;i++)
    {
        y_sp[i].resize(this->X_.size());
        Y_init[i].resize(this->X_.size());
    }
    if (z < 3.5 * B)
    {
        isVapor = false;
        for (int i = 0;i < this->X_.size();i++)
        {
            for (int j = 0;j < 4;j++)
            {
                Y_init[j][i] = this->X_[i] * pow(K[i], pow_v[j]);
            }
        }
    }
    else
    {
        isVapor = true;
        for (int i = 0;i < this->X_.size();i++)
        {
            for (int j = 0;j < ntrial;j++)
            {
                Y_init[j][i] = this->X_[i] * pow(K[i], pow_l[j]);
            }
        }
    }

    for (int i = 0;i < this->X_.size();i++)
    {
        for (int j = 0;j < this->X_.size();j++)
        {
            Y_init[4 + i][j] = 0.0001 / (n - 1);
        }
        Y_init[4 + i][i] = 0.9999;
    }
    scalar sum = 0, z_y, beta;
    scalarList Y(this->X_.size()), y_trial(this->X_.size());
    scalarList alpha_old(this->X_.size()), F_old(this->X_.size());
    scalarList alpha_new(this->X_.size()), F_new(this->X_.size());
    scalarList g(this->X_.size()), s(this->X_.size());
    scalarList dalpha(this->X_.size());
    scalar F_times_s, F_times_g, s_times_g, g_times_g, eps, r;

    for (int itrial = 0;itrial < ntrial;itrial++)
    {
        sum = 0;
        for (int i = 0;i < this->X_.size();i++)
        {
            Y[i] = Y_init[itrial][i];
            sum += Y[i];
        }
        for (int i = 0;i < this->X_.size();i++)
        {
            y_trial[i] = Y[i] / sum;
        }
        z_y = this->Z_gibbs(p, T, &y_trial);
        autoPtr<scalarList> lnphi_y(this->ln_fugacityCoefficient(p, T, z_y, &y_trial));
        for (int i = 0;i < this->X_.size();i++)
        {
            alpha_old[i] = 2 * sqrt(Y[i]);
            F_old[i] = alpha_old[i] / 2 * (lnphi_y()[i] + log(Y[i]) - d[i]);
            alpha_new[i] = 2 * exp(-0.5 * (lnphi_y()[i] - d[i]));
        }
        for (int iter = 0;iter < iter_max;iter++)
        {
            sum = 0;
            for (int i = 0;i < this->X_.size();i++)
            {
                Y[i] = 0.25 * alpha_new[i] * alpha_new[i];
                sum += Y[i];
            }
            for (int i = 0;i < this->X_.size();i++)
            {
                y_trial[i] = Y[i] / sum;
            }
            z_y = this->Z_gibbs(p, T, &y_trial);

            TPD_star[itrial] = 1.;
            beta = 0;

            lnphi_y.reset(this->ln_fugacityCoefficient(p, T, z_y, &y_trial));
            for (int i = 0;i < this->X_.size();i++)
            {
                F_new[i] = alpha_new[i] / 2 * (lnphi_y()[i] + log(Y[i]) - d[i]);
                TPD_star[itrial] += +Y[i] * (lnphi_y()[i] + log(Y[i]) - d[i] - 1);
                beta += (Y[i] - this->X_[i]) * (lnphi_y()[i] + log(Y[i]) - d[i]);

            }
            for (int i = 0;i < this->X_.size();i++)
            {
                g[i] = F_new[i] - F_old[i];
                F_old[i] = F_new[i];
                s[i] = alpha_new[i] - alpha_old[i];
                alpha_old[i] = alpha_new[i];
            }
            F_times_s = 0;
            F_times_g = 0;
            s_times_g = 0;
            g_times_g = 0;
            for (int i = 0;i < this->X_.size();i++)
            {
                F_times_s += F_old[i] * s[i];
                F_times_g += F_old[i] * g[i];
                s_times_g += s[i] * g[i];
                g_times_g += g[i] * g[i];
            }
            eps = 0;
            for (int i = 0;i < this->X_.size();i++)
            {
                dalpha[i] = (g[i] * F_times_s - s[i] * (F_times_s - F_times_g + g_times_g * F_times_s / s_times_g)) / s_times_g - F_old[i];
                alpha_new[i] = alpha_old[i] + dalpha[i];
                eps += dalpha[i] * dalpha[i];
            }
            eps = sqrt(eps);
            r = 2 * TPD_star[itrial] / beta;

            if (eps < tol)
            {

                sum = 0;
                for (int i = 0;i < this->X_.size();i++)
                {
                    sum += Y[i];
                }
                for (int i = 0;i < this->X_.size();i++)
                {
                    y_sp[itrial][i] = Y[i] / sum;
                }
                TPD[itrial] = -log(sum);
                TPD_star[itrial] = 1. - sum;
                mask_stationary[itrial] = true;
                if (TPD_star[itrial] < be_negative)
                {

                    stable = false;
                    if (TPDmin > TPD_star[itrial])
                    {
                        TPDmin = TPD_star[itrial];
                        scalar z2 = this->Z_gibbs(p, T, &y_sp[itrial]);
                        if (z < z2)
                            for (int i = 0;i < this->X_.size();i++)
                            {
                                //K_init_min[i] = Y_init[itrial][i] / this->X_[i];
                                K_init_min[i] = y_sp[itrial][i] / this->X_[i];
                            }
                        else
                            for (int i = 0;i < this->X_.size();i++)
                            {
                                //K_init_min[i] = Y_init[itrial][i] / this->X_[i];
                                K_init_min[i] = this->X_[i] / y_sp[itrial][i];
                            }
                    }
                }
                break;

            }
            if (fabs(r - 1) < 0.1)
            {
                TPD[itrial] = 10.;
                TPD_star[itrial] = 10.;
                mask_stationary[itrial] = false;
                for (int i = 0;i < this->X_.size();i++)
                {
                    y_sp[itrial][i] = -1;
                }
                break;
            }
        }
    }
    return make_tuple(stable, p_K_init_min, isVapor);

}


template <class ThermoMixture>
inline tuple<bool, Foam::autoPtr<List<scalarList>>, bool, autoPtr<scalarList>>  Foam::VLE<ThermoMixture>::solveTPD_BFGS_v2
(
    scalar p, //input
    scalar T  //input
) const
{
    const scalar be_negative = -1.E-8;
    const scalar be_pure_substance = 1. - 1.E-8;
    const scalar be_very_large = 1.E10;

    scalar pow_v[4] = { 1.0,1.0 / 3,-1.0,-1.0 / 3.0 };
    scalar pow_l[4] = { -1.0,-1.0 / 3,1.0,1.0 / 3.0 };

    bool stable = true;
    bool  isVapor;
    scalar TPDmin = 10.0;
    int n = this->X_.size();
    int ntrial = 4 + n;
    autoPtr<List<scalarList>> p_K_init(new List<scalarList>(ntrial));
    List<scalarList>& K_init = p_K_init();
    for (int i = 0;i < ntrial;i++)
    {
        K_init[i].resize(this->X_.size());
    }
    scalar tol = 1.E-9;
    int iter_max = 1000;
    scalarList y_sp[ntrial];
    scalar z = this->Z_gibbs(p, T, &this->X_);
    autoPtr<scalarList> p_TPD_star(new scalarList(ntrial));
    scalarList& TPD_star = p_TPD_star();
    scalarList TPD(ntrial), mask_stationary(ntrial);
    autoPtr<scalarList> lnphi(this->ln_fugacityCoefficient(p, T, z, &this->X_));
    scalarList d(this->X_.size()), K(this->X_.size());
    for (int i = 0;i < this->X_.size();i++)
    {
        d[i] = log(this->X_[i]) + lnphi()[i];
        K[i] = (*this)[i].Pc_ / p * exp(5.373 * (1.0 + (*this)[i].omega_) * (1.0 - (*this)[i].Tc_ / T));
    }
    scalar B = this->B(p, T, &this->X_);
    scalarList Y_init[ntrial];
    for (int i = 0;i < ntrial;i++)
    {
        y_sp[i].resize(this->X_.size());
        Y_init[i].resize(this->X_.size());
    }
    if (z < 3.5 * B)
    {
        isVapor = false;
        for (int i = 0;i < this->X_.size();i++)
        {
            for (int j = 0;j < 4;j++)
            {
                Y_init[j][i] = this->X_[i] * pow(K[i], pow_v[j]);
            }
        }
    }
    else
    {
        isVapor = true;
        for (int i = 0;i < this->X_.size();i++)
        {
            for (int j = 0;j < ntrial;j++)
            {
                Y_init[j][i] = this->X_[i] * pow(K[i], pow_l[j]);
            }
        }
    }

    for (int i = 0;i < this->X_.size();i++)
    {
        for (int j = 0;j < this->X_.size();j++)
        {
            Y_init[4 + i][j] = 0.001 / (n - 1);
        }
        Y_init[4 + i][i] = 0.990;
    }
    scalar sum = 0, z_y, beta;
    scalarList Y(this->X_.size()), y_trial(this->X_.size());
    scalarList alpha_old(this->X_.size()), F_old(this->X_.size());
    scalarList alpha_new(this->X_.size()), F_new(this->X_.size());
    scalarList g(this->X_.size()), s(this->X_.size());
    scalarList dalpha(this->X_.size());
    scalar F_times_s, F_times_g, s_times_g, g_times_g, eps, r;

    for (int itrial = 0;itrial < ntrial;itrial++)
    {
        sum = 0;
        for (int i = 0;i < this->X_.size();i++)
        {
            Y[i] = Y_init[itrial][i];
            sum += Y[i];
        }
        for (int i = 0;i < this->X_.size();i++)
        {
            y_trial[i] = Y[i] / sum;
        }
        z_y = this->Z_gibbs(p, T, &y_trial);
        autoPtr<scalarList> lnphi_y(this->ln_fugacityCoefficient(p, T, z_y, &y_trial));
        for (int i = 0;i < this->X_.size();i++)
        {
            alpha_old[i] = 2 * sqrt(Y[i]);
            F_old[i] = alpha_old[i] / 2 * (lnphi_y()[i] + log(Y[i]) - d[i]);
            alpha_new[i] = 2 * exp(-0.5 * (lnphi_y()[i] - d[i]));
        }
        for (int iter = 0;iter < iter_max;iter++)
        {
            sum = 0;
            for (int i = 0;i < this->X_.size();i++)
            {
                Y[i] = 0.25 * alpha_new[i] * alpha_new[i];
                sum += Y[i];
            }
            for (int i = 0;i < this->X_.size();i++)
            {
                y_trial[i] = Y[i] / sum;
            }
            z_y = this->Z_gibbs(p, T, &y_trial);

            TPD_star[itrial] = 1.;
            beta = 0;

            lnphi_y.reset(this->ln_fugacityCoefficient(p, T, z_y, &y_trial));
            for (int i = 0;i < this->X_.size();i++)
            {
                F_new[i] = alpha_new[i] / 2 * (lnphi_y()[i] + log(Y[i]) - d[i]);
                TPD_star[itrial] += +Y[i] * (lnphi_y()[i] + log(Y[i]) - d[i] - 1);
                beta += (Y[i] - this->X_[i]) * (lnphi_y()[i] + log(Y[i]) - d[i]);

            }
            for (int i = 0;i < this->X_.size();i++)
            {
                g[i] = F_new[i] - F_old[i];
                F_old[i] = F_new[i];
                s[i] = alpha_new[i] - alpha_old[i];
                alpha_old[i] = alpha_new[i];
            }
            F_times_s = 0;
            F_times_g = 0;
            s_times_g = 0;
            g_times_g = 0;
            for (int i = 0;i < this->X_.size();i++)
            {
                F_times_s += F_old[i] * s[i];
                F_times_g += F_old[i] * g[i];
                s_times_g += s[i] * g[i];
                g_times_g += g[i] * g[i];
            }
            eps = 0;
            for (int i = 0;i < this->X_.size();i++)
            {
                dalpha[i] = (g[i] * F_times_s - s[i] * (F_times_s - F_times_g + g_times_g * F_times_s / s_times_g)) / s_times_g - F_old[i];
                alpha_new[i] = alpha_old[i] + dalpha[i];
                eps += dalpha[i] * dalpha[i];
            }
            eps = sqrt(eps);
            r = 2 * TPD_star[itrial] / beta;

            if (eps < tol)
            {

                sum = 0;
                for (int i = 0;i < this->X_.size();i++)
                {
                    sum += Y[i];
                }
                for (int i = 0;i < this->X_.size();i++)
                {
                    y_sp[itrial][i] = Y[i] / sum;
                }
                TPD[itrial] = -log(sum);
                TPD_star[itrial] = 1. - sum;
                mask_stationary[itrial] = true;
                if (TPD_star[itrial] < be_negative)
                {

                    stable = false;

                    scalar z2 = this->Z_gibbs(p, T, &y_sp[itrial]);
                    if (z < z2)
                        for (int i = 0;i < this->X_.size();i++)
                        {
                            //K_init_min[i] = Y_init[itrial][i] / this->X_[i];
                            K_init[itrial][i] = y_sp[itrial][i] / this->X_[i];
                        }
                    else
                        for (int i = 0;i < this->X_.size();i++)
                        {
                            //K_init_min[i] = Y_init[itrial][i] / this->X_[i];
                            K_init[itrial][i] = this->X_[i] / y_sp[itrial][i];
                        }

                }
                break;

            }
            if (fabs(r - 1) < 0.1 || TPD_star[itrial] > 1e40)
            {
                TPD[itrial] = 10.;
                TPD_star[itrial] = 10.;
                mask_stationary[itrial] = false;
                for (int i = 0;i < this->X_.size();i++)
                {
                    y_sp[itrial][i] = -1;
                }
                break;
            }
        }
    }
    return make_tuple(stable, p_K_init, isVapor, p_TPD_star);

}


template <class ThermoMixture>
inline autoPtr<scalarList> Foam::VLE<ThermoMixture>::dvidT(scalar p, scalar T, solution& sol) const
{
    scalarSquareMatrix A(this->X_.size(), Zero);
    scalarList& b = *(new scalarList(this->X_.size()));
    scalar vf = sol.vaporfra;
    autoPtr<scalarList> ddT_Ln_fugcoef_gas, ddT_Ln_fugcoef_liq;
    List<autoPtr<scalarList>> ddxi_Ln_fugcoef_gas(this->X_.size()), ddxi_Ln_fugcoef_liq(this->X_.size());
    ddT_Ln_fugcoef_gas.reset((this->ddT_Ln_fugacityCoefficient(p, T, 1, &sol.X_gas())).ptr());
    for (int i = 0;i < this->X_.size();i++)
        ddxi_Ln_fugcoef_gas[i].reset((this->ddxi_Ln_fugacityCoefficient(p, T, i, 1, &sol.X_gas())).ptr());
    ddT_Ln_fugcoef_liq.reset((this->ddT_Ln_fugacityCoefficient(p, T, 0, &sol.X_liq())).ptr());
    for (int i = 0;i < this->X_.size();i++)
        ddxi_Ln_fugcoef_liq[i].reset((this->ddxi_Ln_fugacityCoefficient(p, T, i, 0, &sol.X_liq())).ptr());
    for (int i = 0;i < this->X_.size();i++)
    {
        for (int j = 0;j < this->X_.size();j++)
        {
            A(i, j) = 1;
            if (i == j)
                A(i, j) -= this->X_[i] / (sol.X_gas()[i] * sol.X_liq()[i]);
            for (int k = 0;k < this->X_.size();k++)
            {
                if (j == k)
                    A(i, j) += (1 - vf) * ddxi_Ln_fugcoef_gas[k]()[i] * (sol.X_gas()[k] - 1);
                else
                    A(i, j) += (1 - vf) * ddxi_Ln_fugcoef_gas[k]()[i] * sol.X_gas()[k];
                if (j == k)
                    A(i, j) += vf * ddxi_Ln_fugcoef_liq[k]()[i] * (sol.X_liq()[k] - 1);
                else
                    A(i, j) += vf * ddxi_Ln_fugcoef_liq[k]()[i] * sol.X_liq()[k];
            }
        }
    }


    for (int i = 0;i < this->X_.size();i++)
        b[i] = (1 - vf) * vf * (ddT_Ln_fugcoef_gas()[i] - ddT_Ln_fugcoef_liq()[i]);

    LUsolve(A, b);
    return autoPtr<scalarList>(&b);
}
template <class ThermoMixture>
inline autoPtr<scalarList> Foam::VLE<ThermoMixture>::dXdT_gas(scalar p, scalar T, solution& sol) const
{
    scalarList& dXdT = *(new scalarList(this->X_.size()));
    autoPtr<scalarList> dvidt;
    dvidt.reset((this->dvidT(p, T, sol)).ptr());
    scalar vf = sol.vaporfra;
    for (int i = 0;i < this->X_.size();i++)
    {
        dXdT[i] = 0;
        for (int j = 0;j < this->X_.size();j++)
        {
            if (i == j)
                dXdT[i] -= 1.0 / vf * (sol.X_gas()[i] - 1) * dvidt()[j];
            else
                dXdT[i] -= 1.0 / vf * sol.X_gas()[i] * dvidt()[j];

        }
    }
    return autoPtr<scalarList>(&dXdT);
}

template <class ThermoMixture>
inline autoPtr<scalarList> Foam::VLE<ThermoMixture>::dXdP_gas(scalar p, scalar T, solution& sol) const
{
    scalarList& dXdT = *(new scalarList(this->X_.size()));
    autoPtr<scalarList> dvidt;
    dvidt.reset((this->dvidP(p, T, sol)).ptr());
    scalar vf = sol.vaporfra;
    for (int i = 0;i < this->X_.size();i++)
    {
        dXdT[i] = 0;
        for (int j = 0;j < this->X_.size();j++)
        {
            if (i == j)
                dXdT[i] -= 1.0 / vf * (sol.X_gas()[i] - 1) * dvidt()[j];
            else
                dXdT[i] -= 1.0 / vf * sol.X_gas()[i] * dvidt()[j];

        }
    }
    return autoPtr<scalarList>(&dXdT);
}
template <class ThermoMixture>
inline autoPtr<scalarList> Foam::VLE<ThermoMixture>::dXdXi_gas(scalar p, scalar T, label di, solution& sol) const
{
    scalarList& dXdT = *(new scalarList(this->X_.size()));
    autoPtr<scalarList> dvidt;
    dvidt.reset((this->dvidXi(p, T, di, sol)).ptr());
    scalar vf = sol.vaporfra;
    for (int i = 0;i < this->X_.size();i++)
    {
        dXdT[i] = 0;
        for (int j = 0;j < this->X_.size();j++)
        {
            dXdT[i] -= dvidt()[j];

        }
        dXdT[i] *= sol.X_gas()[i];
        dXdT[i] += dvidt()[i];
        dXdT[i] /= vf;
    }
    //Info << di << "#########=" << setprecision(20) << (dvidt()[0] - sol.X_gas()[0] * dvidt()[0] - sol.X_gas()[0] * dvidt()[1]) / vf << endl;
    //Info << di << "#########=" << setprecision(15) << -((sol.X_gas()[0] - 1) * dvidt()[0] + sol.X_gas()[0] * dvidt()[1]) / vf << endl;

    return autoPtr<scalarList>(&dXdT);
}

template <class ThermoMixture>
inline autoPtr<scalarList> Foam::VLE<ThermoMixture>::dXdT_liq(scalar p, scalar T, solution& sol) const
{
    scalarList& dXdT = *(new scalarList(this->X_.size()));
    autoPtr<scalarList> dvidt;
    dvidt.reset((this->dvidT(p, T, sol)).ptr());
    scalar vf = sol.vaporfra;
    for (int i = 0;i < this->X_.size();i++)
    {
        dXdT[i] = 0;
        for (int j = 0;j < this->X_.size();j++)
        {
            if (i == j)
                dXdT[i] += 1.0 / (1 - vf) * (sol.X_liq()[i] - 1) * dvidt()[j];
            else
                dXdT[i] += 1.0 / (1 - vf) * sol.X_liq()[i] * dvidt()[j];

        }
    }
    return autoPtr<scalarList>(&dXdT);
}
template <class ThermoMixture>
inline autoPtr<scalarList> Foam::VLE<ThermoMixture>::dXdP_liq(scalar p, scalar T, solution& sol) const
{
    scalarList& dXdT = *(new scalarList(this->X_.size()));
    autoPtr<scalarList> dvidt;
    dvidt.reset((this->dvidP(p, T, sol)).ptr());
    scalar vf = sol.vaporfra;
    for (int i = 0;i < this->X_.size();i++)
    {
        dXdT[i] = 0;
        for (int j = 0;j < this->X_.size();j++)
        {
            if (i == j)
                dXdT[i] += 1.0 / (1 - vf) * (sol.X_liq()[i] - 1) * dvidt()[j];
            else
                dXdT[i] += 1.0 / (1 - vf) * sol.X_liq()[i] * dvidt()[j];

        }
    }
    return autoPtr<scalarList>(&dXdT);
}

template <class ThermoMixture>
inline autoPtr<scalarList> Foam::VLE<ThermoMixture>::dXdXi_liq(scalar p, scalar T, label di, solution& sol) const
{
    scalarList& dXdT = *(new scalarList(this->X_.size()));
    autoPtr<scalarList> dvidxi;
    dvidxi.reset((this->dvidXi(p, T, di, sol)).ptr());
    scalar vf = sol.vaporfra;
    scalar dNvfdNi = 0;
    for (int i = 0;i < this->X_.size();i++)
    {
        dNvfdNi += dvidxi()[i];
    }
    for (int i = 0;i < this->X_.size();i++)
    {
        if (i == di)
        {
            dXdT[i] = 1 / (1 - vf) * (1 - dvidxi()[i]) - sol.X_liq()[i] / (1 - vf) * (1 - dNvfdNi);
        }
        else
        {
            dXdT[i] = 1 / (1 - vf) * (-dvidxi()[i]) - sol.X_liq()[i] / (1 - vf) * (1 - dNvfdNi);
        }

    }
    return autoPtr<scalarList>(&dXdT);
}
/*
template <class ThermoMixture>
inline autoPtr<scalarList> Foam::VLE<ThermoMixture>::dXdXi_liq(scalar p, scalar T, label di, solution& sol) const
{
    scalarList& dXdT = *(new scalarList(this->X_.size()));
    autoPtr<scalarList> dvidt;
    dvidt.reset((this->dvidXi(p, T, di, sol)).ptr());
    scalar vf = sol.vaporfra;
    for (int i = 0;i < this->X_.size();i++)
    {
        dXdT[i] = 0;
        for (int j = 0;j < this->X_.size();j++)
        {
            if (i == j)
                dXdT[i] += 1.0 / (1 - vf) * (sol.X_liq()[i] - 1) * dvidt()[j];
            else
                dXdT[i] += 1.0 / (1 - vf) * sol.X_liq()[i] * dvidt()[j];

        }
    }
    return autoPtr<scalarList>(&dXdT);
}
*/

template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::dWdT_gas(scalar p, scalar T, solution& sol) const
{
    scalar dWdT = 0;
    const VLE<ThermoMixture>& This = *this;
    autoPtr<scalarList> dXdT;
    dXdT.reset((this->dXdT_gas(p, T, sol)).ptr());
    for (int i = 0;i < this->X_.size();i++)
    {
        dWdT += dXdT()[i] * This[i].W() * 1.0e-03;
    }
    return dWdT;
}

template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::dWdP_gas(scalar p, scalar T, solution& sol) const
{
    scalar dWdT = 0;
    const VLE<ThermoMixture>& This = *this;
    autoPtr<scalarList> dXdT;
    dXdT.reset((this->dXdP_gas(p, T, sol)).ptr());
    for (int i = 0;i < this->X_.size();i++)
    {
        dWdT += dXdT()[i] * This[i].W() * 1.0e-03;
    }
    return dWdT;
}
template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::dWdXi_gas(scalar p, scalar T, label di, solution& sol) const
{
    scalar dWdT = 0;
    const VLE<ThermoMixture>& This = *this;
    autoPtr<scalarList> dXdT;
    dXdT.reset((this->dXdXi_gas(p, T, di, sol)).ptr());
    for (int i = 0;i < this->X_.size();i++)
    {
        dWdT += dXdT()[i] * This[i].W() * 1.0e-03;
    }
    return dWdT;
}

template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::dWdXi(label di) const
{
    const VLE<ThermoMixture>& This = *this;
    scalar dWdT = This[di].W() * 1.0e-03;
    for (int i = 0;i < this->X_.size();i++)
    {
        dWdT -= this->X_[i] * This[i].W() * 1.0e-03;
    }
    return dWdT;
}

template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::dWdT_liq(scalar p, scalar T, solution& sol) const
{
    scalar dWdT = 0;
    const PengRobinsonMixture<ThermoMixture>& This = *this;
    autoPtr<scalarList> dXdT;
    dXdT.reset((this->dXdT_liq(p, T, sol)).ptr());
    scalar vf = sol.vaporfra;
    for (int i = 0;i < this->X_.size();i++)
    {
        dWdT += dXdT()[i] * This[i].W() * 1.0e-03;
    }
    return dWdT;
}
template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::dWdXi_liq(scalar p, scalar T, label di, solution& sol) const
{
    scalar dWdT = 0;
    const PengRobinsonMixture<ThermoMixture>& This = *this;
    autoPtr<scalarList> dXdT;
    dXdT.reset((this->dXdXi_liq(p, T, di, sol)).ptr());
    scalar vf = sol.vaporfra;
    for (int i = 0;i < this->X_.size();i++)
    {
        dWdT += dXdT()[i] * This[i].W() * 1.0e-03;
    }
    return dWdT;
}

template <class ThermoMixture>
inline autoPtr<scalarList> Foam::VLE<ThermoMixture>::dvidP(scalar p, scalar T, solution& sol) const
{
    scalarSquareMatrix A(this->X_.size(), Zero);
    scalarList& b = *(new scalarList(this->X_.size()));
    scalar vf = sol.vaporfra;
    autoPtr<scalarList> ddT_Ln_fugcoef_gas, ddT_Ln_fugcoef_liq;
    List<autoPtr<scalarList>> ddxi_Ln_fugcoef_gas(this->X_.size()), ddxi_Ln_fugcoef_liq(this->X_.size());

    ddT_Ln_fugcoef_gas.reset((this->ddP_Ln_fugacityCoefficient(p, T, 1, &sol.X_gas())).ptr());

    for (int i = 0;i < this->X_.size();i++)
        ddxi_Ln_fugcoef_gas[i].reset((this->ddxi_Ln_fugacityCoefficient(p, T, i, 1, &sol.X_gas())).ptr());

    ddT_Ln_fugcoef_liq.reset((this->ddP_Ln_fugacityCoefficient(p, T, 0, &sol.X_liq())).ptr());

    for (int i = 0;i < this->X_.size();i++)
        ddxi_Ln_fugcoef_liq[i].reset((this->ddxi_Ln_fugacityCoefficient(p, T, i, 0, &sol.X_liq())).ptr());
    for (int i = 0;i < this->X_.size();i++)
    {
        for (int j = 0;j < this->X_.size();j++)
        {
            A(i, j) = 1;
            if (i == j)
                A(i, j) -= this->X_[i] / (sol.X_gas()[i] * sol.X_liq()[i]);
            for (int k = 0;k < this->X_.size();k++)
            {
                if (j == k)
                    A(i, j) += (1 - vf) * ddxi_Ln_fugcoef_gas[k]()[i] * (sol.X_gas()[k] - 1);
                else
                    A(i, j) += (1 - vf) * ddxi_Ln_fugcoef_gas[k]()[i] * sol.X_gas()[k];
                if (j == k)
                    A(i, j) += vf * ddxi_Ln_fugcoef_liq[k]()[i] * (sol.X_liq()[k] - 1);
                else
                    A(i, j) += vf * ddxi_Ln_fugcoef_liq[k]()[i] * sol.X_liq()[k];
            }
        }
    }
    for (int i = 0;i < this->X_.size();i++)
        b[i] = (1 - vf) * vf * (ddT_Ln_fugcoef_gas()[i] - ddT_Ln_fugcoef_liq()[i]);
    LUsolve(A, b);
    return autoPtr<scalarList>(&b);
}
template <class ThermoMixture>
inline autoPtr<scalarList> Foam::VLE<ThermoMixture>::dvidXi(scalar p, scalar T, label di, solution& sol) const
{
    scalarSquareMatrix A(this->X_.size(), Zero);
    scalarList& b = *(new scalarList(this->X_.size()));

    scalar vf = sol.vaporfra;
    autoPtr<scalarList> ddT_Ln_fugcoef_gas, ddT_Ln_fugcoef_liq;
    List<autoPtr<scalarList>> ddxi_Ln_fugcoef_gas(this->X_.size()), ddxi_Ln_fugcoef_liq(this->X_.size());

    ddT_Ln_fugcoef_gas.reset((this->ddT_Ln_fugacityCoefficient(p, T, 1, &sol.X_gas())).ptr());

    for (int i = 0;i < this->X_.size();i++)
        ddxi_Ln_fugcoef_gas[i].reset((this->ddxi_Ln_fugacityCoefficient(p, T, i, 1, &sol.X_gas())).ptr());

    ddT_Ln_fugcoef_liq.reset((this->ddT_Ln_fugacityCoefficient(p, T, 0, &sol.X_liq())).ptr());

    for (int i = 0;i < this->X_.size();i++)
        ddxi_Ln_fugcoef_liq[i].reset((this->ddxi_Ln_fugacityCoefficient(p, T, i, 0, &sol.X_liq())).ptr());
    for (int i = 0;i < this->X_.size();i++)
    {
        for (int j = 0;j < this->X_.size();j++)
        {
            A(i, j) = -1.0 / (vf * (1 - vf));
            if (i == j)
                A(i, j) += 1 / (sol.X_liq()[i] * (1 - vf)) + 1 / (sol.X_gas()[i] * vf);
            for (int k = 0;k < this->X_.size();k++)
            {
                if (j == k)
                    A(i, j) += ddxi_Ln_fugcoef_gas[k]()[i] * (1 - sol.X_gas()[k]) / vf;
                else
                    A(i, j) += ddxi_Ln_fugcoef_gas[k]()[i] * (0 - sol.X_gas()[k]) / vf;
                if (j == k)
                    A(i, j) += ddxi_Ln_fugcoef_liq[k]()[i] * (1 - sol.X_liq()[k]) / (1 - vf);
                else
                    A(i, j) += ddxi_Ln_fugcoef_liq[k]()[i] * (0 - sol.X_liq()[k]) / (1 - vf);
            }
        }
    }

    for (int k = 0;k < this->X_.size();k++)
    {
        b[k] = 0;
        for (int j = 0;j < this->X_.size();j++)
        {
            if (j == di)
            {
                for (int s = 0;s < this->X_.size();s++)
                {
                    if (s == j)
                        b[k] += ddxi_Ln_fugcoef_liq[s]()[k] * (1 - sol.X_liq()[s]) / (1 - vf);
                    else
                        b[k] += ddxi_Ln_fugcoef_liq[s]()[k] * (0 - sol.X_liq()[s]) / (1 - vf);
                }

                if (k == j)
                {
                    b[k] += 1 / sol.X_liq()[k] * (1 - sol.X_liq()[k]) / (1 - vf);
                }
                else
                {
                    b[k] += 1 / sol.X_liq()[k] * (0 - sol.X_liq()[k]) / (1 - vf);
                }

            }

        }
    }


    LUsolve(A, b);
    return autoPtr<scalarList>(&b);
}
template <class ThermoMixture>
inline Foam::scalar Foam::VLE<ThermoMixture>::alpha(
    scalar p,
    scalar T,
    solution& sol) const
{
    scalar sv_gas = this->V(p, T, 1, &sol.X_gas()); //gas phase specific volume
    scalar sv_liq = this->V(p, T, 0, &sol.X_liq()); //liquid phase specific volume

    scalar alphagas = sol.vaporfra * sv_gas / (sol.vaporfra * sv_gas + (1.0 - sol.vaporfra) * sv_liq);
    if (alphagas > 1.0)//0.9999999)
    {
        FatalErrorInFunction << "alphagas=" << alphagas << exit(FatalError);
        alphagas = 0.9999999;
    }
    else if (alphagas < 0.0)// 1.0e-07)
    {
        FatalErrorInFunction << "alphagas=" << alphagas << exit(FatalError);
        alphagas = 1.0e-07;
    }

    return alphagas;
}
//-mixing rule
template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::rho(scalar p, scalar T) const
{
    autoPtr<solution> sol(TPn_flash(p, T));
    scalar rhoMixture = p * this->W() / (Z(p, T) * RR * 1.0e-03 * T);
    return rhoMixture;
}
/*
template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::drhodT(scalar p, scalar T) const
{
    autoPtr<solution> sol(TPn_flash(p, T));
    label pfg = 1, pfl = 0;
    scalar alpha = this->alpha(p, T, sol());

    scalar ZGas = this->ThermoMixture::Z(p, T, &pfg, &(sol().X_gas()));
    scalar ZLiq = this->ThermoMixture::Z(p, T, &pfl, &(sol().X_liq()));
    scalar ZMixture = ZGas * alpha + ZLiq * (1.0 - alpha);
    scalar rhoMixture = p * this->W() / (ZMixture * RR * 1.0e-03 * T);
    return rhoMixture;
}*/

template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::rho(scalar p, scalar T, solution& sol) const
{
    scalar alpha = this->alpha(p, T, sol);

    scalar ZGas = this->ThermoMixture::Z(p, T, 1, &(sol.X_gas()));
    scalar ZLiq = this->ThermoMixture::Z(p, T, 0, &(sol.X_liq()));
    scalar ZMixture = ZGas * alpha + ZLiq * (1.0 - alpha);
    scalar rhoMixture = p * this->W() / (ZMixture * RR * 1.0e-03 * T);
    return rhoMixture;
}

template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::drhodT(scalar p, scalar T, solution& sol) const
{
    scalar ZMixture = Z(p, T, sol);
    scalar dZMixturedT = dZdT(p, T, sol);
    scalar rhoMixture = -p * this->W() / (ZMixture * ZMixture * RR * 1.0e-03 * T * T) * (dZMixturedT * T + ZMixture);
    return rhoMixture;
}
template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::drhodP(scalar p, scalar T, solution& sol) const
{
    scalar ZMixture = Z(p, T, sol);
    scalar dZMixturedT = dZdP(p, T, sol);
    scalar rhoMixture = this->W() / (ZMixture * RR * 1.0e-03 * T) - p * this->W() / (ZMixture * ZMixture * RR * 1.0e-03 * T) * dZMixturedT;
    return rhoMixture;
}

template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::drhodXi(scalar p, scalar T, label di, solution& sol) const
{
    scalar ZMixture = Z(p, T, sol);
    scalar dZMixturedT = dZdXi(p, T, di, sol);
    scalar rhoMixture = p * this->dWdXi(di) / (ZMixture * RR * 1.0e-03 * T) - p * this->W() / (ZMixture * ZMixture * RR * 1.0e-03 * T) * dZMixturedT;

    return rhoMixture;
}

template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::Z(scalar p, scalar T) const
{
    autoPtr<solution> sol(TPn_flash(p, T));
    return Z(p, T, sol());
}
template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::Z(scalar p, scalar T, solution& sol) const
{
    scalar vf = sol.vaporfra;
    scalar ZGas = this->ThermoMixture::Z(p, T, 1, &(sol.X_gas()));
    scalar ZLiq = this->ThermoMixture::Z(p, T, 0, &(sol.X_liq()));
    return (vf * ZGas * ZGas + (1 - vf) * ZLiq * ZLiq) / (vf * ZGas + (1 - vf) * ZLiq);
}

template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::dZdT(scalar p, scalar T, solution& sol) const
{

    autoPtr<scalarList> dXdt_gas(this->dXdT_gas(p, T, sol));
    autoPtr<scalarList> dXdt_liq(this->dXdT_liq(p, T, sol));
    scalar dvfdt = this->dvfdT(p, T, sol);
    scalar vf = sol.vaporfra;
    scalar ZGas = this->ThermoMixture::Z(p, T, 1, &(sol.X_gas()));
    scalar ZLiq = this->ThermoMixture::Z(p, T, 0, &(sol.X_liq()));
    scalar dZgasdT = this->ThermoMixture::dZdT(p, T, 1, &(sol.X_gas()));
    scalar dZliqdT = this->ThermoMixture::dZdT(p, T, 0, &(sol.X_liq()));
    for (int i = 0;i < this->X_.size();i++)
    {
        dZgasdT += this->ThermoMixture::dZdxi(p, T, i, 1, &(sol.X_gas())) * dXdt_gas()[i];
        dZliqdT += this->ThermoMixture::dZdxi(p, T, i, 0, &(sol.X_liq())) * dXdt_liq()[i];
    }
    scalar temp = (vf * ZGas + (1 - vf) * ZLiq);
    return (dvfdt * ZGas * ZGas + vf * 2 * ZGas * dZgasdT + 2 * (1 - vf) * ZLiq * dZliqdT - dvfdt * ZLiq * ZLiq) / temp - (vf * ZGas * ZGas + (1 - vf) * ZLiq * ZLiq) / (temp * temp) * (dvfdt * ZGas + vf * dZgasdT + (1 - vf) * dZliqdT - dvfdt * ZLiq);
}

template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::dZdP(scalar p, scalar T, solution& sol) const
{
    autoPtr<scalarList> dXdt_gas(this->dXdP_gas(p, T, sol));
    autoPtr<scalarList> dXdt_liq(this->dXdP_liq(p, T, sol));
    scalar dvfdt = this->dvfdP(p, T, sol);
    scalar vf = sol.vaporfra;
    scalar ZGas = this->ThermoMixture::Z(p, T, 1, &(sol.X_gas()));
    scalar ZLiq = this->ThermoMixture::Z(p, T, 0, &(sol.X_liq()));
    scalar dZgasdT = this->ThermoMixture::dZdP(p, T, 1, &(sol.X_gas()));
    scalar dZliqdT = this->ThermoMixture::dZdP(p, T, 0, &(sol.X_liq()));
    for (int i = 0;i < this->X_.size();i++)
    {
        dZgasdT += this->ThermoMixture::dZdxi(p, T, i, 1, &(sol.X_gas())) * dXdt_gas()[i];
        dZliqdT += this->ThermoMixture::dZdxi(p, T, i, 0, &(sol.X_liq())) * dXdt_liq()[i];
    }
    scalar temp = (vf * ZGas + (1 - vf) * ZLiq);
    return (dvfdt * ZGas * ZGas + vf * 2 * ZGas * dZgasdT + 2 * (1 - vf) * ZLiq * dZliqdT - dvfdt * ZLiq * ZLiq) / temp - (vf * ZGas * ZGas + (1 - vf) * ZLiq * ZLiq) / (temp * temp) * (dvfdt * ZGas + vf * dZgasdT + (1 - vf) * dZliqdT - dvfdt * ZLiq);
}

template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::dZdXi(scalar p, scalar T, label di, solution& sol) const
{

    autoPtr<scalarList> dXdt_gas(this->dXdXi_gas(p, T, di, sol));
    autoPtr<scalarList> dXdt_liq(this->dXdXi_liq(p, T, di, sol));
    scalar dvfdt = this->dvfdXi(p, T, di, sol);
    scalar vf = sol.vaporfra;
    scalar ZGas = this->ThermoMixture::Z(p, T, 1, &(sol.X_gas()));
    scalar ZLiq = this->ThermoMixture::Z(p, T, 0, &(sol.X_liq()));
    scalar dZgasdT = 0;
    scalar dZliqdT = 0;
    for (int i = 0;i < this->X_.size();i++)
    {
        dZgasdT += this->ThermoMixture::dZdxi(p, T, i, 1, &(sol.X_gas())) * dXdt_gas()[i];
        dZliqdT += this->ThermoMixture::dZdxi(p, T, i, 0, &(sol.X_liq())) * dXdt_liq()[i];
    }
    scalar temp = (vf * ZGas + (1 - vf) * ZLiq);
    //Info << di << "#########=" << setprecision(20) << dXdt_gas()[1] << endl;
    //Info << di << "#########=" << setprecision(20) << dZgasdT << endl;
    return (dvfdt * ZGas * ZGas + vf * 2 * ZGas * dZgasdT + 2 * (1 - vf) * ZLiq * dZliqdT - dvfdt * ZLiq * ZLiq) / temp - (vf * ZGas * ZGas + (1 - vf) * ZLiq * ZLiq) / (temp * temp) * (dvfdt * ZGas + vf * dZgasdT + (1 - vf) * dZliqdT - dvfdt * ZLiq);
}

template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::dvfdT(scalar p, scalar T, solution& sol) const
{
    autoPtr<scalarList> dvidT(this->dvidT(p, T, sol));
    scalar dvfdt = 0;
    for (int i = 0;i < this->X_.size();i++)
    {
        dvfdt += dvidT()[i];
    }
    return dvfdt;
}

template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::dvfdP(scalar p, scalar T, solution& sol) const
{
    autoPtr<scalarList> dvidT(this->dvidP(p, T, sol));

    scalar dvfdt = 0;
    for (int i = 0;i < this->X_.size();i++)
    {
        dvfdt += dvidT()[i];
    }
    return dvfdt;
}
template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::dvfdXi(scalar p, scalar T, label di, solution& sol) const
{
    autoPtr<scalarList> dvidT(this->dvidXi(p, T, di, sol));
    scalar dvfdt = -sol.vaporfra;
    for (int i = 0;i < this->X_.size();i++)
    {
        dvfdt += dvidT()[i];
    }

    return dvfdt;
}



template <class ThermoMixture>
inline Foam::scalar Foam::VLE<ThermoMixture>::psi(
    scalar p,
    scalar T) const
{
    return rho(p, T) / p;
}

template <class ThermoMixture>
inline Foam::scalar Foam::VLE<ThermoMixture>::psi(
    scalar p,
    scalar T,
    autoPtr<solution>& sol) const
{
    return rho(p, T, sol()) / p;
}
template <class ThermoMixture>
inline Foam::scalar Foam::VLE<ThermoMixture>::vaporfra(
    scalar p,
    scalar T) const
{
    autoPtr<solution> sol(TPn_flash(p, T));
    return sol().vaporfra;
}

template <class ThermoMixture>
inline Foam::scalar Foam::VLE<ThermoMixture>::vaporfra(
    scalar p,
    scalar T,
    autoPtr<solution>& sol) const
{
    return sol().vaporfra;
}

template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::Ha(scalar p, scalar T) const
{
    autoPtr<solution> sol(TPn_flash(p, T));

    scalar mw_gas = this->W(&(sol().X_gas())); //kg/mol
    scalar mw_mixture = this->W();
    scalar ygas = sol().vaporfra * mw_gas / mw_mixture;
    //Info << sol().vaporfra << endl;
    //Info << this->ThermoMixture::Ha(p, T, 1, &(sol().X_gas())) << endl;
    //Info << this->ThermoMixture::Ha(p, T, 0, &(sol().X_liq())) << endl;
    //Info << "mw_gas / mw_mixture=" << mw_gas / mw_mixture << endl;
    //Info << ygas << endl;
    return this->ThermoMixture::Ha(p, T, 1, &(sol().X_gas())) * ygas + this->ThermoMixture::Ha(p, T, 0, &(sol().X_liq())) * (1.0 - ygas);
}

template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::Hc() const
{

    return this->ThermoMixture::Hc(&(this->X_));
}
template <class ThermoMixture>
inline Foam::scalar Foam::VLE<ThermoMixture>::Hs(
    const scalar p,
    const scalar T) const
{
    return Ha(p, T) - Hc();
}
template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::Cp(scalar p, scalar T) const
{
    autoPtr<solution> sol(TPn_flash(p, T));
    return Cp(p, T, sol());
}
template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::Cp(scalar p, scalar T, solution& sol) const
{
    scalar mw_gas = this->W(&(sol.X_gas())); //kg/mol
    scalar dmw_gasdT = this->dWdT_gas(p, T, sol);
    scalar mw_mixture = this->W();
    scalar ygas = sol.vaporfra * mw_gas / mw_mixture;
    autoPtr<scalarList> dvidt(this->dvidT(p, T, sol));
    scalar dvfdT = 0;
    for (int i = 0;i < this->X_.size();i++)
    {
        dvfdT += dvidt()[i];
    }
    scalar dygasdT = (dvfdT * mw_gas + sol.vaporfra * dmw_gasdT) / mw_mixture;
    scalar dHadT_gas = this->dHadT(p, T, 1, &(sol.X_gas()));
    autoPtr<scalarList> dXdt_gas(this->dXdT_gas(p, T, sol));
    for (int i = 0;i < this->X_.size();i++)
    {
        dHadT_gas += this->dHadxi(p, T, i, 1, &(sol.X_gas())) * dXdt_gas()[i];
    }
    scalar dHadT_liq = this->dHadT(p, T, 0, &(sol.X_liq()));
    autoPtr<scalarList> dXdt_liq(this->dXdT_liq(p, T, sol));
    for (int i = 0;i < this->X_.size();i++)
    {
        dHadT_liq += this->dHadxi(p, T, i, 0, &(sol.X_liq())) * dXdt_liq()[i];
    }
    return this->ThermoMixture::Ha(p, T, 1, &(sol.X_gas())) * dygasdT + dHadT_gas * ygas
        - this->ThermoMixture::Ha(p, T, 0, &(sol.X_liq())) * dygasdT + dHadT_liq * (1.0 - ygas);
}

template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::Cp_Hs(scalar p, scalar T, solution& sol) const
{

    scalar mw_gas = this->W(&(sol.X_gas())); //kg/mol
    scalar dmw_gasdT = this->dWdT_gas(p, T, sol);
    scalar mw_mixture = this->W();
    scalar ygas = sol.vaporfra * mw_gas / mw_mixture;
    autoPtr<scalarList> dvidt(this->dvidT(p, T, sol));
    scalar dvfdT = 0;
    for (int i = 0;i < this->X_.size();i++)
    {
        dvfdT += dvidt()[i];
    }
    scalar dygasdT = (dvfdT * mw_gas + sol.vaporfra * dmw_gasdT) / mw_mixture;

    scalar dHadT_gas = this->dHsdT(p, T, 1, &(sol.X_gas()));
    autoPtr<scalarList> dXdt_gas(this->dXdT_gas(p, T, sol));
    for (int i = 0;i < this->X_.size();i++)
    {
        dHadT_gas += this->dHsdxi(p, T, i, 1, &(sol.X_gas())) * dXdt_gas()[i];
    }
    scalar dHadT_liq = this->dHsdT(p, T, 0, &(sol.X_liq()));
    autoPtr<scalarList> dXdt_liq(this->dXdT_liq(p, T, sol));
    for (int i = 0;i < this->X_.size();i++)
    {
        dHadT_liq += this->dHsdxi(p, T, i, 0, &(sol.X_liq())) * dXdt_liq()[i];
    }
    return this->ThermoMixture::Hs(p, T, 1, &(sol.X_gas())) * dygasdT + dHadT_gas * ygas
        - this->ThermoMixture::Hs(p, T, 0, &(sol.X_liq())) * dygasdT + dHadT_liq * (1.0 - ygas);
}


template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::dHadXi(scalar p, scalar T, label di, solution& sol) const
{

    scalar mw_gas = this->W(&(sol.X_gas())); //kg/mol
    scalar dmw_gasdT = this->dWdXi_gas(p, T, di, sol);
    scalar mw_mixture = this->W();
    scalar ygas = sol.vaporfra * mw_gas / mw_mixture;
    autoPtr<scalarList> dvidt(this->dvidXi(p, T, di, sol));
    scalar dvfdT = 0;
    for (int i = 0;i < this->X_.size();i++)
    {
        dvfdT += dvidt()[i];
    }
    const VLE<ThermoMixture>& This = *this;
    scalar dygasdT = (dvfdT * mw_gas + sol.vaporfra * dmw_gasdT) / mw_mixture - ygas / mw_mixture * This[di].W() * 1.0e-03;

    scalar dHadT_gas = 0.0;
    autoPtr<scalarList> dXdt_gas(this->dXdXi_gas(p, T, di, sol));
    for (int i = 0;i < this->X_.size();i++)
    {
        dHadT_gas += this->dHadxi(p, T, i, 1, &(sol.X_gas())) * dXdt_gas()[i];
    }

    scalar dHadT_liq = 0.0;
    autoPtr<scalarList> dXdt_liq(this->dXdXi_liq(p, T, di, sol));
    for (int i = 0;i < this->X_.size();i++)
    {
        dHadT_liq += this->dHadxi(p, T, i, 0, &(sol.X_liq())) * dXdt_liq()[i];
    }

    return this->ThermoMixture::Ha(p, T, 1, &(sol.X_gas())) * dygasdT + dHadT_gas * ygas
        - this->ThermoMixture::Ha(p, T, 0, &(sol.X_liq())) * dygasdT + dHadT_liq * (1.0 - ygas);

}
template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::dHsdXi(scalar p, scalar T, label di, solution& sol) const
{

    scalar mw_gas = this->W(&(sol.X_gas())); //kg/mol
    scalar dmw_gasdT = this->dWdXi_gas(p, T, di, sol);
    scalar mw_mixture = this->W();
    scalar ygas = sol.vaporfra * mw_gas / mw_mixture;
    autoPtr<scalarList> dvidt(this->dvidXi(p, T, di, sol));
    scalar dvfdT = 0;
    for (int i = 0;i < this->X_.size();i++)
    {
        dvfdT += dvidt()[i];
    }
    const VLE<ThermoMixture>& This = *this;
    scalar dygasdT = (dvfdT * mw_gas + sol.vaporfra * dmw_gasdT) / mw_mixture - ygas / mw_mixture * This[di].W() * 1.0e-03;

    scalar dHadT_gas = 0.0;
    autoPtr<scalarList> dXdt_gas(this->dXdXi_gas(p, T, di, sol));
    for (int i = 0;i < this->X_.size();i++)
    {
        dHadT_gas += this->dHsdxi(p, T, i, 1, &(sol.X_gas())) * dXdt_gas()[i];
    }

    scalar dHadT_liq = 0.0;
    autoPtr<scalarList> dXdt_liq(this->dXdXi_liq(p, T, di, sol));
    for (int i = 0;i < this->X_.size();i++)
    {
        dHadT_liq += this->dHsdxi(p, T, i, 0, &(sol.X_liq())) * dXdt_liq()[i];
    }
    //Info << ",dygasd=" << this->dHsdxi(p, T, 0, 1, &(sol.X_gas())) << endl;
    return this->ThermoMixture::Hs(p, T, 1, &(sol.X_gas())) * dygasdT + dHadT_gas * ygas
        - this->ThermoMixture::Hs(p, T, 0, &(sol.X_liq())) * dygasdT + dHadT_liq * (1.0 - ygas);

}


template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::dHadP(scalar p, scalar T, solution& sol) const
{

    scalar mw_gas = this->W(&(sol.X_gas())); //kg/mol
    scalar dmw_gasdT = this->dWdP_gas(p, T, sol);
    scalar mw_mixture = this->W();
    scalar ygas = sol.vaporfra * mw_gas / mw_mixture;
    autoPtr<scalarList> dvidt(this->dvidP(p, T, sol));
    scalar dvfdT = 0;
    for (int i = 0;i < this->X_.size();i++)
    {
        dvfdT += dvidt()[i];
    }
    scalar dygasdT = (dvfdT * mw_gas + sol.vaporfra * dmw_gasdT) / mw_mixture;

    scalar dHadT_gas = this->ThermoMixture::dHadP(p, T, 1, &(sol.X_gas()));
    autoPtr<scalarList> dXdt_gas(this->dXdP_gas(p, T, sol));
    for (int i = 0;i < this->X_.size();i++)
    {
        dHadT_gas += this->dHadxi(p, T, i, 1, &(sol.X_gas())) * dXdt_gas()[i];
    }

    scalar dHadT_liq = this->ThermoMixture::dHadP(p, T, 0, &(sol.X_liq()));
    autoPtr<scalarList> dXdt_liq(this->dXdP_liq(p, T, sol));
    for (int i = 0;i < this->X_.size();i++)
    {
        dHadT_liq += this->dHadxi(p, T, i, 0, &(sol.X_liq())) * dXdt_liq()[i];
    }

    return this->ThermoMixture::Ha(p, T, 1, &(sol.X_gas())) * dygasdT + dHadT_gas * ygas
        - this->ThermoMixture::Ha(p, T, 0, &(sol.X_liq())) * dygasdT + dHadT_liq * (1.0 - ygas);

}
template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::dHsdP(scalar p, scalar T, solution& sol) const
{

    scalar mw_gas = this->W(&(sol.X_gas())); //kg/mol
    scalar dmw_gasdT = this->dWdP_gas(p, T, sol);
    scalar mw_mixture = this->W();
    scalar ygas = sol.vaporfra * mw_gas / mw_mixture;
    autoPtr<scalarList> dvidt(this->dvidP(p, T, sol));
    scalar dvfdT = 0;
    for (int i = 0;i < this->X_.size();i++)
    {
        dvfdT += dvidt()[i];
    }
    scalar dygasdT = (dvfdT * mw_gas + sol.vaporfra * dmw_gasdT) / mw_mixture;

    scalar dHadT_gas = this->ThermoMixture::dHsdP(p, T, 1, &(sol.X_gas()));
    autoPtr<scalarList> dXdt_gas(this->dXdP_gas(p, T, sol));
    for (int i = 0;i < this->X_.size();i++)
    {
        dHadT_gas += this->dHsdxi(p, T, i, 1, &(sol.X_gas())) * dXdt_gas()[i];
    }

    scalar dHadT_liq = this->ThermoMixture::dHsdP(p, T, 0, &(sol.X_liq()));
    autoPtr<scalarList> dXdt_liq(this->dXdP_liq(p, T, sol));
    for (int i = 0;i < this->X_.size();i++)
    {
        dHadT_liq += this->dHsdxi(p, T, i, 0, &(sol.X_liq())) * dXdt_liq()[i];
    }

    return this->ThermoMixture::Hs(p, T, 1, &(sol.X_gas())) * dygasdT + dHadT_gas * ygas
        - this->ThermoMixture::Hs(p, T, 0, &(sol.X_liq())) * dygasdT + dHadT_liq * (1.0 - ygas);

}
template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::dTdP_HP(scalar p, scalar T, solution& sol) const
{

    return -dHadP(p, T, sol) / Cp(p, T, sol);

}
template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::dTdH_HP(scalar p, scalar T, solution& sol) const
{
    return 1 / Cp(p, T, sol);
}
template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::dTdXi_HP(scalar p, scalar T, label di, solution& sol) const
{
    return -dHadXi(p, T, di, sol) / Cp(p, T, sol);
}


template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::dTdP_HsP(scalar p, scalar T, solution& sol) const
{

    return -dHsdP(p, T, sol) / Cp_Hs(p, T, sol);

}
template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::dTdH_HsP(scalar p, scalar T, solution& sol) const
{
    return 1 / Cp_Hs(p, T, sol);
}
template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::dTdXi_HsP(scalar p, scalar T, label di, solution& sol) const
{

    return -dHsdXi(p, T, di, sol) / Cp_Hs(p, T, sol);
}


template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::drhodH_HP(scalar p, scalar T, solution& sol) const
{
    //return drhodT(p,T,sol)*dTdH(p,T,sol)+drhodP(p,T,sol)*dPdH(p,T,sol)+drhodXi(p,T,di,sol)*dPdXi(p,T,di,sol);
    return drhodT(p, T, sol) * dTdH_HP(p, T, sol);
}
template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::drhodP_HP(scalar p, scalar T, solution& sol) const
{

    return drhodT(p, T, sol) * dTdP_HP(p, T, sol) + drhodP(p, T, sol);
}

template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::drhodXi_HP(scalar p, scalar T, label di, solution& sol) const
{
    return drhodT(p, T, sol) * dTdXi_HP(p, T, di, sol) + drhodXi(p, T, di, sol);
}


template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::drhodH_HsP(scalar p, scalar T, solution& sol) const
{

    return drhodT(p, T, sol) * dTdH_HsP(p, T, sol);
}

template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::drhodP_HsP(scalar p, scalar T, solution& sol) const
{

    return drhodT(p, T, sol) * dTdP_HsP(p, T, sol) + drhodP(p, T, sol);
}

template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::drhodXi_HsP(scalar p, scalar T, label di, solution& sol) const
{

    return drhodT(p, T, sol) * dTdXi_HsP(p, T, di, sol) + drhodXi(p, T, di, sol);
}



template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::drhoPdH_HP(scalar p, scalar T, solution& sol) const
{

    return  drhodH_HP(p, T, sol) / p;
}
template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::drhoPdP_HP(scalar p, scalar T, solution& sol) const
{
    scalar r = rho(p, T, sol);
    return drhodP_HP(p, T, sol) / p - r / (p * p);
}

template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::drhoPdXi_HP(scalar p, scalar T, label di, solution& sol) const
{
    return  drhodXi_HP(p, T, di, sol) / p;
}


template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::drhoPdH_HsP(scalar p, scalar T, solution& sol) const
{
    return  drhodH_HsP(p, T, sol) / p;
}
template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::drhoPdP_HsP(scalar p, scalar T, solution& sol) const
{
    scalar r = rho(p, T, sol);
    return drhodP_HsP(p, T, sol) / p - r / (p * p);
}

template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::drhoPdXi_HsP(scalar p, scalar T, label di, solution& sol) const
{

    return  drhodXi_HsP(p, T, di, sol) / p;
}


template <class ThermoMixture>
inline autoPtr<scalarList> Foam::VLE<ThermoMixture>::drhoPdXHP_HsP(scalar p, scalar T, solution& sol) const
{

    int nsp = this->X_.size();
    scalar temp;
    scalarList& grad = *(new scalarList(nsp + 2));
    const VLE<ThermoMixture>& This = *this;
    scalar vf = sol.vaporfra;
    scalarList& X_gas = sol.X_gas();
    scalarList& X_liq = sol.X_liq();
    scalar rho_value = rho(p, T, sol);





    autoPtr<scalarList> dvidP_value(this->dvidP(p, T, sol));
    autoPtr<scalarList> dvidT_value(this->dvidT(p, T, sol));
    List<autoPtr<scalarList>> dvidXi_value(nsp);
    for (int i = 0;i < nsp;i++)
    {
        dvidXi_value[i] = this->dvidXi(p, T, i, sol);
    }
    scalar dvfdP_value = 0, dvfdT_value = 0;
    scalarList dvfdXi_value(nsp);
    scalarList dHsdxi_gas_value(nsp), dHsdxi_liq_value(nsp);
    for (int i = 0;i < nsp;i++)
    {
        dvfdT_value += dvidT_value()[i];
        dvfdP_value += dvidP_value()[i];
        dvfdXi_value[i] = -vf;
        for (int j = 0;j < nsp;j++)
        {
            dvfdXi_value[i] += dvidXi_value[i]()[j];
        }
        dHsdxi_gas_value[i] = this->dHsdxi(p, T, i, 1, &X_gas);
        dHsdxi_liq_value[i] = this->dHsdxi(p, T, i, 0, &X_liq);
    }

    scalar dHsdT_gas = this->ThermoMixture::dHsdT(p, T, 1, &X_gas);
    scalar dHsdT_liq = this->ThermoMixture::dHsdT(p, T, 0, &X_liq);
    scalar dHsdP_gas = this->ThermoMixture::dHsdP(p, T, 1, &X_gas);
    scalar dHsdP_liq = this->ThermoMixture::dHsdP(p, T, 0, &X_liq);
    scalarList dHsdXi_gas(nsp);
    scalarList dHsdXi_liq(nsp);

    //autoPtr<scalarList> dXdT_gas_value(this->dXdT_gas(p, T, sol));
    scalarList dXdT_gas_value(nsp);
    //autoPtr<scalarList> dXdT_liq_value(this->dXdT_liq(p, T, sol));
    scalarList dXdT_liq_value(nsp);
    //autoPtr<scalarList> dXdP_gas_value(this->dXdP_gas(p, T, sol));
    scalarList dXdP_gas_value(nsp);
    //autoPtr<scalarList> dXdP_liq_value(this->dXdP_liq(p, T, sol));
    scalarList dXdP_liq_value(nsp);
    for (int i = 0;i < nsp;i++)
    {
        dXdT_gas_value[i] = (dvidT_value()[i] - X_gas[i] * dvfdT_value) / vf;
        dXdT_liq_value[i] = -(dvidT_value()[i] - X_liq[i] * dvfdT_value) / (1 - vf);
        dXdP_gas_value[i] = (dvidP_value()[i] - X_gas[i] * dvfdP_value) / vf;
        dXdP_liq_value[i] = -(dvidP_value()[i] - X_liq[i] * dvfdP_value) / (1 - vf);
    }
    //List<autoPtr<scalarList>> dXdXi_gas_value(nsp);
    //List<autoPtr<scalarList>> dXdXi_liq_value(nsp);
    List<scalarList> dXdXi_gas_value(nsp);
    List<scalarList> dXdXi_liq_value(nsp);
    /*
    for (int i = 0;i < nsp;i++)
    {
        dXdXi_gas_value[i] = this->dXdXi_gas(p, T, i, sol);
        dXdXi_liq_value[i] = this->dXdXi_liq(p, T, i, sol);
    }*/
    for (int i = 0;i < nsp;i++)
    {
        dXdXi_gas_value[i].resize(nsp);
        dXdXi_liq_value[i].resize(nsp);
        for (int j = 0;j < nsp;j++)
        {
            dXdXi_gas_value[i][j] = (dvidXi_value[i]()[j] - X_gas[j] * (dvfdXi_value[i] + vf)) / vf;
            if (i == j)
            {
                dXdXi_liq_value[i][j] = 1 / (1 - vf) * (1 - dvidXi_value[i]()[j]) - X_liq[j] / (1 - vf) * (1 - dvfdXi_value[i] - vf);
            }
            else
            {
                dXdXi_liq_value[i][j] = 1 / (1 - vf) * (-dvidXi_value[i]()[j]) - X_liq[j] / (1 - vf) * (1 - dvfdXi_value[i] - vf);
            }
        }

    }
    //Info << "!!!!!!=" << setprecision(20) << (dvidXi_value[0]()[0] - X_gas[0] * (dvfdXi_value[0] + vf)) / vf << endl;
    //Info << "!!!!!!=" << setprecision(15) << -((X_gas[0] - 1) * dvidXi_value[0]()[0] + X_gas[0] * dvidXi_value[0]()[1]) / vf << endl;
    for (int i = 0;i < nsp;i++)
    {
        dHsdT_gas += dHsdxi_gas_value[i] * dXdT_gas_value[i];
        dHsdT_liq += dHsdxi_liq_value[i] * dXdT_liq_value[i];
        dHsdP_gas += dHsdxi_gas_value[i] * dXdP_gas_value[i];
        dHsdP_liq += dHsdxi_liq_value[i] * dXdP_liq_value[i];
        dHsdXi_gas[i] = 0;
        dHsdXi_liq[i] = 0;
        for (int j = 0;j < nsp;j++)
        {
            dHsdXi_gas[i] += dHsdxi_gas_value[j] * dXdXi_gas_value[i][j];
            dHsdXi_liq[i] += dHsdxi_liq_value[j] * dXdXi_liq_value[i][j];
        }
    }

    scalar mw_gas = this->W(&X_gas); //kg/mol
    scalar mw_mixture = this->W();
    scalar ygas = vf * mw_gas / mw_mixture;
    //scalar dmw_gasdT = this->dWdT_gas(p, T, sol);
    scalar dmw_gasdT = 0;
    //scalar dmw_gasdP = this->dWdP_gas(p, T, sol);
    scalar dmw_gasdP = 0;
    for (int i = 0;i < nsp;i++)
    {
        dmw_gasdT += dXdT_gas_value[i] * This[i].W() * 1.0e-03;
        dmw_gasdP += dXdP_gas_value[i] * This[i].W() * 1.0e-03;
    }
    scalarList dmw_gasdXi(nsp);
    for (int i = 0;i < nsp;i++)
    {
        dmw_gasdXi[i] = 0;
        for (int j = 0;j < nsp;j++)
        {
            dmw_gasdXi[i] += dXdXi_gas_value[i][j] * This[j].W() * 1.0e-03;

        }
        //dmw_gasdXi[i] = this->dWdXi_gas(p, T, i, sol);
    }
    scalar dygasdT = (dvfdT_value * mw_gas + vf * dmw_gasdT) / mw_mixture;
    scalar dygasdP = (dvfdP_value * mw_gas + vf * dmw_gasdP) / mw_mixture;
    scalarList dygasdXi(nsp);
    for (int i = 0;i < nsp;i++)
    {
        dygasdXi[i] = ((dvfdXi_value[i] + vf) * mw_gas + vf * dmw_gasdXi[i]) / mw_mixture - ygas / mw_mixture * This[i].W() * 1.0e-03;
    }

    scalar Hs_gas = this->ThermoMixture::Hs(p, T, 1, &X_gas);
    scalar Hs_liq = this->ThermoMixture::Hs(p, T, 0, &X_liq);
    scalar Cp_Hs_value = Hs_gas * dygasdT + dHsdT_gas * ygas
        - Hs_liq * dygasdT + dHsdT_liq * (1.0 - ygas);
    scalar dHsdP_value = Hs_gas * dygasdP + dHsdP_gas * ygas
        - Hs_liq * dygasdP + dHsdP_liq * (1.0 - ygas);
    scalarList dHsdXi_value(nsp);
    for (int i = 0;i < nsp;i++)
    {
        dHsdXi_value[i] = Hs_gas * dygasdXi[i] + dHsdXi_gas[i] * ygas
            - Hs_liq * dygasdXi[i] + dHsdXi_liq[i] * (1.0 - ygas);
    }

    //Info << "aa=" << dHsdxi_gas_value[0] << "," << dHsdXi(p, T, 0, sol) << endl;



    scalar ZMixture = Z(p, T, sol);
    scalar ZGas = this->ThermoMixture::Z(p, T, 1, &X_gas);
    scalar ZLiq = this->ThermoMixture::Z(p, T, 0, &X_liq);
    scalar dZgasdT = this->ThermoMixture::dZdT(p, T, 1, &X_gas);
    scalar dZliqdT = this->ThermoMixture::dZdT(p, T, 0, &X_liq);
    scalar dZgasdP = this->ThermoMixture::dZdP(p, T, 1, &X_gas);
    scalar dZliqdP = this->ThermoMixture::dZdP(p, T, 0, &X_liq);
    scalarList dZgasdXi(nsp), dZliqdXi(nsp);
    scalarList dZdxi_gas_value(nsp), dZdxi_liq_value(nsp);
    for (int i = 0;i < nsp;i++)
    {
        dZdxi_gas_value[i] = this->ThermoMixture::dZdxi(p, T, i, 1, &X_gas);
        dZdxi_liq_value[i] = this->ThermoMixture::dZdxi(p, T, i, 0, &X_liq);
    }
    for (int i = 0;i < nsp;i++)
    {
        dZgasdT += dZdxi_gas_value[i] * dXdT_gas_value[i];
        dZliqdT += dZdxi_liq_value[i] * dXdT_liq_value[i];
        dZgasdP += dZdxi_gas_value[i] * dXdP_gas_value[i];
        dZliqdP += dZdxi_liq_value[i] * dXdP_liq_value[i];
        dZgasdXi[i] = 0;
        dZliqdXi[i] = 0;
        for (int j = 0;j < nsp;j++)
        {
            dZgasdXi[i] += dZdxi_gas_value[j] * dXdXi_gas_value[i][j];
            dZliqdXi[i] += dZdxi_liq_value[j] * dXdXi_liq_value[i][j];
        }
    }
    temp = (vf * ZGas + (1 - vf) * ZLiq);
    scalar dZMixturedT = (dvfdT_value * ZGas * ZGas + vf * 2 * ZGas * dZgasdT + 2 * (1 - vf) * ZLiq * dZliqdT - dvfdT_value * ZLiq * ZLiq) / temp - (vf * ZGas * ZGas + (1 - vf) * ZLiq * ZLiq) / (temp * temp) * (dvfdT_value * ZGas + vf * dZgasdT + (1 - vf) * dZliqdT - dvfdT_value * ZLiq);

    scalar dZMixturedP = (dvfdP_value * ZGas * ZGas + vf * 2 * ZGas * dZgasdP + 2 * (1 - vf) * ZLiq * dZliqdP - dvfdP_value * ZLiq * ZLiq) / temp - (vf * ZGas * ZGas + (1 - vf) * ZLiq * ZLiq) / (temp * temp) * (dvfdP_value * ZGas + vf * dZgasdP + (1 - vf) * dZliqdP - dvfdP_value * ZLiq);

    scalarList dZMixturedXi(nsp);

    for (int i = 0;i < nsp;i++)
    {
        dZMixturedXi[i] = (dvfdXi_value[i] * ZGas * ZGas + vf * 2 * ZGas * dZgasdXi[i] + 2 * (1 - vf) * ZLiq * dZliqdXi[i] - dvfdXi_value[i] * ZLiq * ZLiq) / temp - (vf * ZGas * ZGas + (1 - vf) * ZLiq * ZLiq) / (temp * temp) * (dvfdXi_value[i] * ZGas + vf * dZgasdXi[i] + (1 - vf) * dZliqdXi[i] - dvfdXi_value[i] * ZLiq);
    }

    //Info << "!!!!!!=" << setprecision(20) << dXdXi_gas_value[0][1] << endl;


    //scalar dZMixturedT = dZdT(p, T, sol);
    //scalar dZMixturedP = dZdP(p, T, sol);
    //scalarList dZMixturedXi(nsp);
    /*for (int i = 0;i < nsp;i++)
    {
        dZMixturedXi[i] = dZdXi(p, T, i, sol);
    }
    */
    scalar drhodT_value = -p * mw_mixture / (ZMixture * ZMixture * RR * 1.0e-03 * T * T) * (dZMixturedT * T + ZMixture);
    scalar drhodP_value = mw_mixture / (ZMixture * RR * 1.0e-03 * T) - p * mw_mixture / (ZMixture * ZMixture * RR * 1.0e-03 * T) * dZMixturedP;
    scalarList drhodXi_value(nsp);
    for (int i = 0;i < nsp;i++)
    {
        drhodXi_value[i] = p * (This[i].W() * 1.0e-03 - mw_mixture) / (ZMixture * RR * 1.0e-03 * T) - p * mw_mixture / (ZMixture * ZMixture * RR * 1.0e-03 * T) * dZMixturedXi[i];
        //dHsdXi_value[i] = dHsdXi(p, T, i, sol);
    }//this->dWdXi(i)

    //scalar drhodT_value = drhodT(p, T, sol);
    //scalar drhodP_value = drhodP(p, T, sol);
    // scalar Cp_Hs_value = Cp_Hs(p, T, sol);
    //scalar dHsdP_value = dHsdP(p, T, sol);
    //scalarList drhodXi_value(nsp);
    //scalarList dHsdXi_value(nsp);
    /*
    for (int i = 0;i < nsp;i++)
    {
        drhodXi_value[i] = drhodXi(p, T, i, sol);
        dHsdXi_value[i] = dHsdXi(p, T, i, sol);
    }*/


    for (int i = 0;i < nsp;i++)
    {
        grad[i] = (-drhodT_value * dHsdXi_value[i] / Cp_Hs_value + drhodXi_value[i]) / p;
    }
    grad[nsp] = (drhodT_value / Cp_Hs_value) / p;
    grad[nsp + 1] = (-drhodT_value * dHsdP_value / Cp_Hs_value + drhodP_value) / p - rho_value / (p * p);
    return autoPtr<scalarList>(&grad);
}


template <class ThermoMixture>
inline autoPtr<scalarRectangularMatrix> Foam::VLE<ThermoMixture>::dTrhoPdXHP_HsP(scalar p, scalar T, solution& sol) const
{

    int nsp = this->X_.size();
    scalar temp;
    scalarRectangularMatrix& grad = *(new scalarRectangularMatrix(nsp + 2, 2));
    //scalarList& grad = *(new scalarList(nsp + 2));
    const VLE<ThermoMixture>& This = *this;
    scalar vf = sol.vaporfra;
    scalarList& X_gas = sol.X_gas();
    scalarList& X_liq = sol.X_liq();
    scalar rho_value = rho(p, T, sol);



    autoPtr<scalarList> dvidP_value(this->dvidP(p, T, sol));
    autoPtr<scalarList> dvidT_value(this->dvidT(p, T, sol));
    List<autoPtr<scalarList>> dvidXi_value(nsp);
    for (int i = 0;i < nsp;i++)
    {
        dvidXi_value[i] = this->dvidXi(p, T, i, sol);
    }
    scalar dvfdP_value = 0, dvfdT_value = 0;
    scalarList dvfdXi_value(nsp);
    scalarList dHsdxi_gas_value(nsp), dHsdxi_liq_value(nsp);
    for (int i = 0;i < nsp;i++)
    {
        dvfdT_value += dvidT_value()[i];
        dvfdP_value += dvidP_value()[i];
        dvfdXi_value[i] = -vf;
        for (int j = 0;j < nsp;j++)
        {
            dvfdXi_value[i] += dvidXi_value[i]()[j];
        }
        dHsdxi_gas_value[i] = this->dHsdxi(p, T, i, 1, &X_gas);
        dHsdxi_liq_value[i] = this->dHsdxi(p, T, i, 0, &X_liq);
    }

    scalar dHsdT_gas = this->ThermoMixture::dHsdT(p, T, 1, &X_gas);
    scalar dHsdT_liq = this->ThermoMixture::dHsdT(p, T, 0, &X_liq);
    scalar dHsdP_gas = this->ThermoMixture::dHsdP(p, T, 1, &X_gas);
    scalar dHsdP_liq = this->ThermoMixture::dHsdP(p, T, 0, &X_liq);
    scalarList dHsdXi_gas(nsp);
    scalarList dHsdXi_liq(nsp);

    //autoPtr<scalarList> dXdT_gas_value(this->dXdT_gas(p, T, sol));
    scalarList dXdT_gas_value(nsp);
    //autoPtr<scalarList> dXdT_liq_value(this->dXdT_liq(p, T, sol));
    scalarList dXdT_liq_value(nsp);
    //autoPtr<scalarList> dXdP_gas_value(this->dXdP_gas(p, T, sol));
    scalarList dXdP_gas_value(nsp);
    //autoPtr<scalarList> dXdP_liq_value(this->dXdP_liq(p, T, sol));
    scalarList dXdP_liq_value(nsp);
    for (int i = 0;i < nsp;i++)
    {
        dXdT_gas_value[i] = (dvidT_value()[i] - X_gas[i] * dvfdT_value) / vf;
        dXdT_liq_value[i] = -(dvidT_value()[i] - X_liq[i] * dvfdT_value) / (1 - vf);
        dXdP_gas_value[i] = (dvidP_value()[i] - X_gas[i] * dvfdP_value) / vf;
        dXdP_liq_value[i] = -(dvidP_value()[i] - X_liq[i] * dvfdP_value) / (1 - vf);
    }
    //List<autoPtr<scalarList>> dXdXi_gas_value(nsp);
    //List<autoPtr<scalarList>> dXdXi_liq_value(nsp);
    List<scalarList> dXdXi_gas_value(nsp);
    List<scalarList> dXdXi_liq_value(nsp);
    /*
    for (int i = 0;i < nsp;i++)
    {
        dXdXi_gas_value[i] = this->dXdXi_gas(p, T, i, sol);
        dXdXi_liq_value[i] = this->dXdXi_liq(p, T, i, sol);
    }*/
    for (int i = 0;i < nsp;i++)
    {
        dXdXi_gas_value[i].resize(nsp);
        dXdXi_liq_value[i].resize(nsp);
        for (int j = 0;j < nsp;j++)
        {
            dXdXi_gas_value[i][j] = (dvidXi_value[i]()[j] - X_gas[j] * (dvfdXi_value[i] + vf)) / vf;
            if (i == j)
            {
                dXdXi_liq_value[i][j] = 1 / (1 - vf) * (1 - dvidXi_value[i]()[j]) - X_liq[j] / (1 - vf) * (1 - dvfdXi_value[i] - vf);
            }
            else
            {
                dXdXi_liq_value[i][j] = 1 / (1 - vf) * (-dvidXi_value[i]()[j]) - X_liq[j] / (1 - vf) * (1 - dvfdXi_value[i] - vf);
            }
        }

    }
    //Info << "!!!!!!=" << setprecision(20) << (dvidXi_value[0]()[0] - X_gas[0] * (dvfdXi_value[0] + vf)) / vf << endl;
    //Info << "!!!!!!=" << setprecision(15) << -((X_gas[0] - 1) * dvidXi_value[0]()[0] + X_gas[0] * dvidXi_value[0]()[1]) / vf << endl;
    for (int i = 0;i < nsp;i++)
    {
        dHsdT_gas += dHsdxi_gas_value[i] * dXdT_gas_value[i];
        dHsdT_liq += dHsdxi_liq_value[i] * dXdT_liq_value[i];
        dHsdP_gas += dHsdxi_gas_value[i] * dXdP_gas_value[i];
        dHsdP_liq += dHsdxi_liq_value[i] * dXdP_liq_value[i];
        dHsdXi_gas[i] = 0;
        dHsdXi_liq[i] = 0;
        for (int j = 0;j < nsp;j++)
        {
            dHsdXi_gas[i] += dHsdxi_gas_value[j] * dXdXi_gas_value[i][j];
            dHsdXi_liq[i] += dHsdxi_liq_value[j] * dXdXi_liq_value[i][j];
        }
    }

    scalar mw_gas = this->W(&X_gas); //kg/mol
    scalar mw_mixture = this->W();
    scalar ygas = vf * mw_gas / mw_mixture;
    //scalar dmw_gasdT = this->dWdT_gas(p, T, sol);
    scalar dmw_gasdT = 0;
    //scalar dmw_gasdP = this->dWdP_gas(p, T, sol);
    scalar dmw_gasdP = 0;
    for (int i = 0;i < nsp;i++)
    {
        dmw_gasdT += dXdT_gas_value[i] * This[i].W() * 1.0e-03;
        dmw_gasdP += dXdP_gas_value[i] * This[i].W() * 1.0e-03;
    }
    scalarList dmw_gasdXi(nsp);
    for (int i = 0;i < nsp;i++)
    {
        dmw_gasdXi[i] = 0;
        for (int j = 0;j < nsp;j++)
        {
            dmw_gasdXi[i] += dXdXi_gas_value[i][j] * This[j].W() * 1.0e-03;

        }
        //dmw_gasdXi[i] = this->dWdXi_gas(p, T, i, sol);
    }
    scalar dygasdT = (dvfdT_value * mw_gas + vf * dmw_gasdT) / mw_mixture;
    scalar dygasdP = (dvfdP_value * mw_gas + vf * dmw_gasdP) / mw_mixture;
    scalarList dygasdXi(nsp);
    for (int i = 0;i < nsp;i++)
    {
        dygasdXi[i] = ((dvfdXi_value[i] + vf) * mw_gas + vf * dmw_gasdXi[i]) / mw_mixture - ygas / mw_mixture * This[i].W() * 1.0e-03;
    }

    scalar Hs_gas = this->ThermoMixture::Hs(p, T, 1, &X_gas);
    scalar Hs_liq = this->ThermoMixture::Hs(p, T, 0, &X_liq);
    scalar Cp_Hs_value = Hs_gas * dygasdT + dHsdT_gas * ygas
        - Hs_liq * dygasdT + dHsdT_liq * (1.0 - ygas);
    scalar dHsdP_value = Hs_gas * dygasdP + dHsdP_gas * ygas
        - Hs_liq * dygasdP + dHsdP_liq * (1.0 - ygas);
    scalarList dHsdXi_value(nsp);
    for (int i = 0;i < nsp;i++)
    {
        dHsdXi_value[i] = Hs_gas * dygasdXi[i] + dHsdXi_gas[i] * ygas
            - Hs_liq * dygasdXi[i] + dHsdXi_liq[i] * (1.0 - ygas);
    }

    //Info << "aa=" << dHsdxi_gas_value[0] << "," << dHsdXi(p, T, 0, sol) << endl;



    scalar ZMixture = Z(p, T, sol);
    scalar ZGas = this->ThermoMixture::Z(p, T, 1, &X_gas);
    scalar ZLiq = this->ThermoMixture::Z(p, T, 0, &X_liq);
    scalar dZgasdT = this->ThermoMixture::dZdT(p, T, 1, &X_gas);
    scalar dZliqdT = this->ThermoMixture::dZdT(p, T, 0, &X_liq);
    scalar dZgasdP = this->ThermoMixture::dZdP(p, T, 1, &X_gas);
    scalar dZliqdP = this->ThermoMixture::dZdP(p, T, 0, &X_liq);
    scalarList dZgasdXi(nsp), dZliqdXi(nsp);
    scalarList dZdxi_gas_value(nsp), dZdxi_liq_value(nsp);
    for (int i = 0;i < nsp;i++)
    {
        dZdxi_gas_value[i] = this->ThermoMixture::dZdxi(p, T, i, 1, &X_gas);
        dZdxi_liq_value[i] = this->ThermoMixture::dZdxi(p, T, i, 0, &X_liq);
    }
    for (int i = 0;i < nsp;i++)
    {
        dZgasdT += dZdxi_gas_value[i] * dXdT_gas_value[i];
        dZliqdT += dZdxi_liq_value[i] * dXdT_liq_value[i];
        dZgasdP += dZdxi_gas_value[i] * dXdP_gas_value[i];
        dZliqdP += dZdxi_liq_value[i] * dXdP_liq_value[i];
        dZgasdXi[i] = 0;
        dZliqdXi[i] = 0;
        for (int j = 0;j < nsp;j++)
        {
            dZgasdXi[i] += dZdxi_gas_value[j] * dXdXi_gas_value[i][j];
            dZliqdXi[i] += dZdxi_liq_value[j] * dXdXi_liq_value[i][j];
        }
    }
    temp = (vf * ZGas + (1 - vf) * ZLiq);
    scalar dZMixturedT = (dvfdT_value * ZGas * ZGas + vf * 2 * ZGas * dZgasdT + 2 * (1 - vf) * ZLiq * dZliqdT - dvfdT_value * ZLiq * ZLiq) / temp - (vf * ZGas * ZGas + (1 - vf) * ZLiq * ZLiq) / (temp * temp) * (dvfdT_value * ZGas + vf * dZgasdT + (1 - vf) * dZliqdT - dvfdT_value * ZLiq);

    scalar dZMixturedP = (dvfdP_value * ZGas * ZGas + vf * 2 * ZGas * dZgasdP + 2 * (1 - vf) * ZLiq * dZliqdP - dvfdP_value * ZLiq * ZLiq) / temp - (vf * ZGas * ZGas + (1 - vf) * ZLiq * ZLiq) / (temp * temp) * (dvfdP_value * ZGas + vf * dZgasdP + (1 - vf) * dZliqdP - dvfdP_value * ZLiq);

    scalarList dZMixturedXi(nsp);

    for (int i = 0;i < nsp;i++)
    {
        dZMixturedXi[i] = (dvfdXi_value[i] * ZGas * ZGas + vf * 2 * ZGas * dZgasdXi[i] + 2 * (1 - vf) * ZLiq * dZliqdXi[i] - dvfdXi_value[i] * ZLiq * ZLiq) / temp - (vf * ZGas * ZGas + (1 - vf) * ZLiq * ZLiq) / (temp * temp) * (dvfdXi_value[i] * ZGas + vf * dZgasdXi[i] + (1 - vf) * dZliqdXi[i] - dvfdXi_value[i] * ZLiq);
    }

    //Info << "!!!!!!=" << setprecision(20) << dXdXi_gas_value[0][1] << endl;


    //scalar dZMixturedT = dZdT(p, T, sol);
    //scalar dZMixturedP = dZdP(p, T, sol);
    //scalarList dZMixturedXi(nsp);
    /*for (int i = 0;i < nsp;i++)
    {
        dZMixturedXi[i] = dZdXi(p, T, i, sol);
    }
    */
    scalar drhodT_value = -p * mw_mixture / (ZMixture * ZMixture * RR * 1.0e-03 * T * T) * (dZMixturedT * T + ZMixture);
    scalar drhodP_value = mw_mixture / (ZMixture * RR * 1.0e-03 * T) - p * mw_mixture / (ZMixture * ZMixture * RR * 1.0e-03 * T) * dZMixturedP;
    scalarList drhodXi_value(nsp);
    for (int i = 0;i < nsp;i++)
    {
        drhodXi_value[i] = p * (This[i].W() * 1.0e-03 - mw_mixture) / (ZMixture * RR * 1.0e-03 * T) - p * mw_mixture / (ZMixture * ZMixture * RR * 1.0e-03 * T) * dZMixturedXi[i];
        //dHsdXi_value[i] = dHsdXi(p, T, i, sol);
    }//this->dWdXi(i)

    //scalar drhodT_value = drhodT(p, T, sol);
    //scalar drhodP_value = drhodP(p, T, sol);
    // scalar Cp_Hs_value = Cp_Hs(p, T, sol);
    //scalar dHsdP_value = dHsdP(p, T, sol);
    //scalarList drhodXi_value(nsp);
    //scalarList dHsdXi_value(nsp);
    /*
    for (int i = 0;i < nsp;i++)
    {
        drhodXi_value[i] = drhodXi(p, T, i, sol);
        dHsdXi_value[i] = dHsdXi(p, T, i, sol);
    }*/


    for (int i = 0;i < nsp;i++)
    {
        grad[i][0] = -dHsdXi_value[i] / Cp_Hs_value;
        grad[i][1] = (-drhodT_value * dHsdXi_value[i] / Cp_Hs_value + drhodXi_value[i]) / p;
    }
    grad[nsp][0] = 1 / Cp_Hs_value;
    grad[nsp][1] = (drhodT_value / Cp_Hs_value) / p;
    grad[nsp + 1][0] = -dHsdP_value / Cp_Hs_value;
    grad[nsp + 1][1] = (-drhodT_value * dHsdP_value / Cp_Hs_value + drhodP_value) / p - rho_value / (p * p);
    return autoPtr<scalarRectangularMatrix>(&grad);
}

template <class ThermoMixture>
inline autoPtr<scalarRectangularMatrix> Foam::VLE<ThermoMixture>::dTrhoPvfdXHP_HsP(scalar p, scalar T, solution& sol) const
{

    int nsp = this->X_.size();
    scalar temp;
    scalarRectangularMatrix& grad = *(new scalarRectangularMatrix(nsp + 2, 3));
    //scalarList& grad = *(new scalarList(nsp + 2));
    const VLE<ThermoMixture>& This = *this;
    scalar vf = sol.vaporfra;
    scalarList& X_gas = sol.X_gas();
    scalarList& X_liq = sol.X_liq();
    scalar rho_value = rho(p, T, sol);



    autoPtr<scalarList> dvidP_value(this->dvidP(p, T, sol));
    autoPtr<scalarList> dvidT_value(this->dvidT(p, T, sol));
    List<autoPtr<scalarList>> dvidXi_value(nsp);
    for (int i = 0;i < nsp;i++)
    {
        dvidXi_value[i] = this->dvidXi(p, T, i, sol);
    }
    scalar dvfdP_value = 0, dvfdT_value = 0;
    scalarList dvfdXi_value(nsp);
    scalarList dHsdxi_gas_value(nsp), dHsdxi_liq_value(nsp);
    for (int i = 0;i < nsp;i++)
    {
        dvfdT_value += dvidT_value()[i];
        dvfdP_value += dvidP_value()[i];
        dvfdXi_value[i] = -vf;
        for (int j = 0;j < nsp;j++)
        {
            dvfdXi_value[i] += dvidXi_value[i]()[j];
        }
        dHsdxi_gas_value[i] = this->dHsdxi(p, T, i, 1, &X_gas);
        dHsdxi_liq_value[i] = this->dHsdxi(p, T, i, 0, &X_liq);
    }

    scalar dHsdT_gas = this->ThermoMixture::dHsdT(p, T, 1, &X_gas);
    scalar dHsdT_liq = this->ThermoMixture::dHsdT(p, T, 0, &X_liq);
    scalar dHsdP_gas = this->ThermoMixture::dHsdP(p, T, 1, &X_gas);
    scalar dHsdP_liq = this->ThermoMixture::dHsdP(p, T, 0, &X_liq);
    scalarList dHsdXi_gas(nsp);
    scalarList dHsdXi_liq(nsp);

    //autoPtr<scalarList> dXdT_gas_value(this->dXdT_gas(p, T, sol));
    scalarList dXdT_gas_value(nsp);
    //autoPtr<scalarList> dXdT_liq_value(this->dXdT_liq(p, T, sol));
    scalarList dXdT_liq_value(nsp);
    //autoPtr<scalarList> dXdP_gas_value(this->dXdP_gas(p, T, sol));
    scalarList dXdP_gas_value(nsp);
    //autoPtr<scalarList> dXdP_liq_value(this->dXdP_liq(p, T, sol));
    scalarList dXdP_liq_value(nsp);
    for (int i = 0;i < nsp;i++)
    {
        dXdT_gas_value[i] = (dvidT_value()[i] - X_gas[i] * dvfdT_value) / vf;
        dXdT_liq_value[i] = -(dvidT_value()[i] - X_liq[i] * dvfdT_value) / (1 - vf);
        dXdP_gas_value[i] = (dvidP_value()[i] - X_gas[i] * dvfdP_value) / vf;
        dXdP_liq_value[i] = -(dvidP_value()[i] - X_liq[i] * dvfdP_value) / (1 - vf);
    }
    //List<autoPtr<scalarList>> dXdXi_gas_value(nsp);
    //List<autoPtr<scalarList>> dXdXi_liq_value(nsp);
    List<scalarList> dXdXi_gas_value(nsp);
    List<scalarList> dXdXi_liq_value(nsp);
    /*
    for (int i = 0;i < nsp;i++)
    {
        dXdXi_gas_value[i] = this->dXdXi_gas(p, T, i, sol);
        dXdXi_liq_value[i] = this->dXdXi_liq(p, T, i, sol);
    }*/
    for (int i = 0;i < nsp;i++)
    {
        dXdXi_gas_value[i].resize(nsp);
        dXdXi_liq_value[i].resize(nsp);
        for (int j = 0;j < nsp;j++)
        {
            dXdXi_gas_value[i][j] = (dvidXi_value[i]()[j] - X_gas[j] * (dvfdXi_value[i] + vf)) / vf;
            if (i == j)
            {
                dXdXi_liq_value[i][j] = 1 / (1 - vf) * (1 - dvidXi_value[i]()[j]) - X_liq[j] / (1 - vf) * (1 - dvfdXi_value[i] - vf);
            }
            else
            {
                dXdXi_liq_value[i][j] = 1 / (1 - vf) * (-dvidXi_value[i]()[j]) - X_liq[j] / (1 - vf) * (1 - dvfdXi_value[i] - vf);
            }
        }

    }
    //Info << "!!!!!!=" << setprecision(20) << (dvidXi_value[0]()[0] - X_gas[0] * (dvfdXi_value[0] + vf)) / vf << endl;
    //Info << "!!!!!!=" << setprecision(15) << -((X_gas[0] - 1) * dvidXi_value[0]()[0] + X_gas[0] * dvidXi_value[0]()[1]) / vf << endl;
    for (int i = 0;i < nsp;i++)
    {
        dHsdT_gas += dHsdxi_gas_value[i] * dXdT_gas_value[i];
        dHsdT_liq += dHsdxi_liq_value[i] * dXdT_liq_value[i];
        dHsdP_gas += dHsdxi_gas_value[i] * dXdP_gas_value[i];
        dHsdP_liq += dHsdxi_liq_value[i] * dXdP_liq_value[i];
        dHsdXi_gas[i] = 0;
        dHsdXi_liq[i] = 0;
        for (int j = 0;j < nsp;j++)
        {
            dHsdXi_gas[i] += dHsdxi_gas_value[j] * dXdXi_gas_value[i][j];
            dHsdXi_liq[i] += dHsdxi_liq_value[j] * dXdXi_liq_value[i][j];
        }
    }

    scalar mw_gas = this->W(&X_gas); //kg/mol
    scalar mw_mixture = this->W();
    scalar ygas = vf * mw_gas / mw_mixture;
    //scalar dmw_gasdT = this->dWdT_gas(p, T, sol);
    scalar dmw_gasdT = 0;
    //scalar dmw_gasdP = this->dWdP_gas(p, T, sol);
    scalar dmw_gasdP = 0;
    for (int i = 0;i < nsp;i++)
    {
        dmw_gasdT += dXdT_gas_value[i] * This[i].W() * 1.0e-03;
        dmw_gasdP += dXdP_gas_value[i] * This[i].W() * 1.0e-03;
    }
    scalarList dmw_gasdXi(nsp);
    for (int i = 0;i < nsp;i++)
    {
        dmw_gasdXi[i] = 0;
        for (int j = 0;j < nsp;j++)
        {
            dmw_gasdXi[i] += dXdXi_gas_value[i][j] * This[j].W() * 1.0e-03;

        }
        //dmw_gasdXi[i] = this->dWdXi_gas(p, T, i, sol);
    }
    scalar dygasdT = (dvfdT_value * mw_gas + vf * dmw_gasdT) / mw_mixture;
    scalar dygasdP = (dvfdP_value * mw_gas + vf * dmw_gasdP) / mw_mixture;
    scalarList dygasdXi(nsp);
    for (int i = 0;i < nsp;i++)
    {
        dygasdXi[i] = ((dvfdXi_value[i] + vf) * mw_gas + vf * dmw_gasdXi[i]) / mw_mixture - ygas / mw_mixture * This[i].W() * 1.0e-03;
    }

    scalar Hs_gas = this->ThermoMixture::Hs(p, T, 1, &X_gas);
    scalar Hs_liq = this->ThermoMixture::Hs(p, T, 0, &X_liq);
    scalar Cp_Hs_value = Hs_gas * dygasdT + dHsdT_gas * ygas
        - Hs_liq * dygasdT + dHsdT_liq * (1.0 - ygas);
    scalar dHsdP_value = Hs_gas * dygasdP + dHsdP_gas * ygas
        - Hs_liq * dygasdP + dHsdP_liq * (1.0 - ygas);
    scalarList dHsdXi_value(nsp);
    for (int i = 0;i < nsp;i++)
    {
        dHsdXi_value[i] = Hs_gas * dygasdXi[i] + dHsdXi_gas[i] * ygas
            - Hs_liq * dygasdXi[i] + dHsdXi_liq[i] * (1.0 - ygas);
    }

    //Info << "aa=" << dHsdxi_gas_value[0] << "," << dHsdXi(p, T, 0, sol) << endl;



    scalar ZMixture = Z(p, T, sol);
    scalar ZGas = this->ThermoMixture::Z(p, T, 1, &X_gas);
    scalar ZLiq = this->ThermoMixture::Z(p, T, 0, &X_liq);
    scalar dZgasdT = this->ThermoMixture::dZdT(p, T, 1, &X_gas);
    scalar dZliqdT = this->ThermoMixture::dZdT(p, T, 0, &X_liq);
    scalar dZgasdP = this->ThermoMixture::dZdP(p, T, 1, &X_gas);
    scalar dZliqdP = this->ThermoMixture::dZdP(p, T, 0, &X_liq);
    scalarList dZgasdXi(nsp), dZliqdXi(nsp);
    scalarList dZdxi_gas_value(nsp), dZdxi_liq_value(nsp);
    for (int i = 0;i < nsp;i++)
    {
        dZdxi_gas_value[i] = this->ThermoMixture::dZdxi(p, T, i, 1, &X_gas);
        dZdxi_liq_value[i] = this->ThermoMixture::dZdxi(p, T, i, 0, &X_liq);
    }
    for (int i = 0;i < nsp;i++)
    {
        dZgasdT += dZdxi_gas_value[i] * dXdT_gas_value[i];
        dZliqdT += dZdxi_liq_value[i] * dXdT_liq_value[i];
        dZgasdP += dZdxi_gas_value[i] * dXdP_gas_value[i];
        dZliqdP += dZdxi_liq_value[i] * dXdP_liq_value[i];
        dZgasdXi[i] = 0;
        dZliqdXi[i] = 0;
        for (int j = 0;j < nsp;j++)
        {
            dZgasdXi[i] += dZdxi_gas_value[j] * dXdXi_gas_value[i][j];
            dZliqdXi[i] += dZdxi_liq_value[j] * dXdXi_liq_value[i][j];
        }
    }
    temp = (vf * ZGas + (1 - vf) * ZLiq);
    scalar dZMixturedT = (dvfdT_value * ZGas * ZGas + vf * 2 * ZGas * dZgasdT + 2 * (1 - vf) * ZLiq * dZliqdT - dvfdT_value * ZLiq * ZLiq) / temp - (vf * ZGas * ZGas + (1 - vf) * ZLiq * ZLiq) / (temp * temp) * (dvfdT_value * ZGas + vf * dZgasdT + (1 - vf) * dZliqdT - dvfdT_value * ZLiq);

    scalar dZMixturedP = (dvfdP_value * ZGas * ZGas + vf * 2 * ZGas * dZgasdP + 2 * (1 - vf) * ZLiq * dZliqdP - dvfdP_value * ZLiq * ZLiq) / temp - (vf * ZGas * ZGas + (1 - vf) * ZLiq * ZLiq) / (temp * temp) * (dvfdP_value * ZGas + vf * dZgasdP + (1 - vf) * dZliqdP - dvfdP_value * ZLiq);

    scalarList dZMixturedXi(nsp);

    for (int i = 0;i < nsp;i++)
    {
        dZMixturedXi[i] = (dvfdXi_value[i] * ZGas * ZGas + vf * 2 * ZGas * dZgasdXi[i] + 2 * (1 - vf) * ZLiq * dZliqdXi[i] - dvfdXi_value[i] * ZLiq * ZLiq) / temp - (vf * ZGas * ZGas + (1 - vf) * ZLiq * ZLiq) / (temp * temp) * (dvfdXi_value[i] * ZGas + vf * dZgasdXi[i] + (1 - vf) * dZliqdXi[i] - dvfdXi_value[i] * ZLiq);
    }

    //Info << "!!!!!!=" << setprecision(20) << dXdXi_gas_value[0][1] << endl;


    //scalar dZMixturedT = dZdT(p, T, sol);
    //scalar dZMixturedP = dZdP(p, T, sol);
    //scalarList dZMixturedXi(nsp);
    /*for (int i = 0;i < nsp;i++)
    {
        dZMixturedXi[i] = dZdXi(p, T, i, sol);
    }
    */
    scalar drhodT_value = -p * mw_mixture / (ZMixture * ZMixture * RR * 1.0e-03 * T * T) * (dZMixturedT * T + ZMixture);
    scalar drhodP_value = mw_mixture / (ZMixture * RR * 1.0e-03 * T) - p * mw_mixture / (ZMixture * ZMixture * RR * 1.0e-03 * T) * dZMixturedP;
    scalarList drhodXi_value(nsp);
    for (int i = 0;i < nsp;i++)
    {
        drhodXi_value[i] = p * (This[i].W() * 1.0e-03 - mw_mixture) / (ZMixture * RR * 1.0e-03 * T) - p * mw_mixture / (ZMixture * ZMixture * RR * 1.0e-03 * T) * dZMixturedXi[i];
        //dHsdXi_value[i] = dHsdXi(p, T, i, sol);
    }//this->dWdXi(i)

    //scalar drhodT_value = drhodT(p, T, sol);
    //scalar drhodP_value = drhodP(p, T, sol);
    // scalar Cp_Hs_value = Cp_Hs(p, T, sol);
    //scalar dHsdP_value = dHsdP(p, T, sol);
    //scalarList drhodXi_value(nsp);
    //scalarList dHsdXi_value(nsp);
    /*
    for (int i = 0;i < nsp;i++)
    {
        drhodXi_value[i] = drhodXi(p, T, i, sol);
        dHsdXi_value[i] = dHsdXi(p, T, i, sol);
    }*/


    for (int i = 0;i < nsp;i++)
    {
        grad[i][0] = -dHsdXi_value[i] / Cp_Hs_value;
        grad[i][1] = (-drhodT_value * dHsdXi_value[i] / Cp_Hs_value + drhodXi_value[i]) / p;
        grad[i][2] = -dvfdT_value * dHsdXi_value[i] / Cp_Hs_value + dvfdXi_value[i];

    }
    grad[nsp][0] = 1 / Cp_Hs_value;
    grad[nsp][1] = (drhodT_value / Cp_Hs_value) / p;
    grad[nsp][2] = dvfdT_value / Cp_Hs_value;
    grad[nsp + 1][0] = -dHsdP_value / Cp_Hs_value;
    grad[nsp + 1][1] = (-drhodT_value * dHsdP_value / Cp_Hs_value + drhodP_value) / p - rho_value / (p * p);
    grad[nsp + 1][2] = -dvfdT_value * dHsdP_value / Cp_Hs_value + dvfdP_value;
    return autoPtr<scalarRectangularMatrix>(&grad);
}












template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::dvfdH_HP(scalar p, scalar T, solution& sol) const
{
    return dvfdT(p, T, sol) * dTdH_HP(p, T, sol);
}

template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::dvfdP_HP(scalar p, scalar T, solution& sol) const
{
    return dvfdT(p, T, sol) * dTdP_HP(p, T, sol) + dvfdP(p, T, sol);
}

template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::dvfdXi_HP(scalar p, scalar T, label di, solution& sol) const
{
    return dvfdT(p, T, sol) * dTdXi_HP(p, T, di, sol) + dvfdXi(p, T, di, sol);
}


template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::dvfdH_HsP(scalar p, scalar T, solution& sol) const
{
    return dvfdT(p, T, sol) * dTdH_HsP(p, T, sol);
}

template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::dvfdP_HsP(scalar p, scalar T, solution& sol) const
{
    return dvfdT(p, T, sol) * dTdP_HsP(p, T, sol) + dvfdP(p, T, sol);
}

template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::dvfdXi_HsP(scalar p, scalar T, label di, solution& sol) const
{
    return dvfdT(p, T, sol) * dTdXi_HsP(p, T, di, sol) + dvfdXi(p, T, di, sol);
}

template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::Cp_num(scalar p, scalar T) const
{
    scalar dt = 0.01;
    return (Hs(p, T + dt) - Hs(p, T)) / dt;
}
template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::CpMCv(scalar p, scalar T) const
{
    autoPtr<solution> sol(TPn_flash(p, T));
    //label pfg = 1, pfl = 0;
    scalar mw_gas = this->W(&(sol().X_gas())); //kg/mol
    scalar mw_mixture = this->W();
    scalar ygas = sol().vaporfra * mw_gas / mw_mixture;
    return this->ThermoMixture::CpMCv(p, T, 1, &(sol().X_gas())) * ygas + this->ThermoMixture::CpMCv(p, T, 0, &(sol().X_liq())) * (1.0 - ygas);
}

template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::S(scalar p, scalar T) const
{
    autoPtr<solution> sol(TPn_flash(p, T));
    //label pfg = 1, pfl = 0;
    scalar mw_gas = this->W(&(sol().X_gas())); //kg/mol
    scalar mw_mixture = this->W();
    scalar ygas = sol().vaporfra * mw_gas / mw_mixture;
    return this->ThermoMixture::S(p, T, 1, &(sol().X_gas())) * ygas + this->ThermoMixture::S(p, T, 0, &(sol().X_liq())) * (1.0 - ygas);
}

template <class ThermoMixture>
inline Foam::scalar Foam::VLE<ThermoMixture>::limit(
    const scalar T) const
{
    scalar T_t = T;
    for (int i = 0; i < this->N_; i++)
    {
        T_t = (*this)[i].limit(T_t);
    }
    return T_t;
}

template <class ThermoMixture>
inline Foam::scalar Foam::VLE<ThermoMixture>::dGdT //todo
(
    const scalar p,
    const scalar T) const
{
    return 0;
}

template <class ThermoMixture>
inline Foam::scalar Foam::VLE<ThermoMixture>::dCpdT //todo
(
    const scalar p,
    const scalar T) const
{
    return 0;
}

//Temp E

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template <class ThermoMixture>
Foam::VLE<ThermoMixture> Foam::operator==(
    const VLE<ThermoMixture>& a,
    const VLE<ThermoMixture>& b)
{
    return VLE<ThermoMixture>(a);
}

template <class ThermoMixture>
Foam::VLE<ThermoMixture> Foam::operator*(
    const scalar s,
    const VLE<ThermoMixture>& b)
{
    return VLE<ThermoMixture>(b);
}

// ************************************************************************* //
