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
#include "string.H"
#include "labelList.H"
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
    const dictionary &dict)
{
    return autoPtr<VLE<ThermoMixture>>(
        new VLE<ThermoMixture>(dict));
}

template <class ThermoMixture>
inline Foam::VLE<ThermoMixture>::VLE(
    const word &name,
    PtrList<SingleThermoType> &speciesData,
    const speciesTable &specieNames,
    const dictionary &thermoDict)
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
inline label Foam::VLE<ThermoMixture>::noramlize(scalarList &X_in) const //noramlize, error return 1
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
    solution &sol) const
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
    solution &Sol = *pS;
    const int maxiter = 1000;
    const int igas = 1, iliquid = 0;
    int ierr = 0;
    //int iaccelation = 1;
    scalar tol_fug = 1.0e-13; //Convergence criteria
    scalar tol_fug2 = 1.0;
    scalar coeff_tpn = 1;
    scalar epsilon1 = 1.0e+10;

    //scalar gibbsenergy = 0.0;
    autoPtr<scalarList> fugcoef_liq, fugcoef_gas;              //, fugcoef_mix[num_sp] = {0.0};
    scalarList X_liq0(this->N_, Zero), X_gas0(this->N_, Zero); //, fugcoef_mix[num_sp] = {0.0};

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

    fugcoef_liq.reset(this->fugacityCoefficient(p, T, &iliquid).ptr()); //0.30637832246552882,0.057000001844270946,0.74015085224153332
    fugcoef_gas.reset(this->fugacityCoefficient(p, T, &igas).ptr());    //0.30637832246552882,0.057000001844270946,0.74015085224153332

    //Info <<"fugcoef_liq="<<fugcoef_liq()<<endl;
    //Info <<"fugcoef_gas="<<fugcoef_gas()<<endl;
    Sol.equalconstant.reset(this->wilsoneq(p, T).ptr()); //1.4037894924124943,0.008729616583312786
    //Sol.equalconstant.reset(this->myinitial(p, T).ptr());
    //Info << "equalconstant= " << Sol.equalconstant() << endl;
    /*
    Sol.equalconstant()[0] = 300.98332688687367;
    Sol.equalconstant()[1] = 0.012591909335566309;
    Sol.equalconstant()[2] = 3302.6955988917484;
    */
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
        ierr = this->solve_vaporf(p, T, Sol); //v=-0.73386572585804888,
        //Info << "vaporfra= " << Sol.vaporfra << endl;
        // only for Hn flash at high temperature and pressure conditions

        if (ierr == 4 || ierr == 3)
        {
            WarningInFunction << "VLE "
                              << "p=" << p << " T=" << T << " X=" << this->X_ << endl;
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

            fugcoef_liq.reset(this->fugacityCoefficient(p, T, &iliquid, &(Sol.X_liq())).ptr());
            fugcoef_gas.reset(this->fugacityCoefficient(p, T, &igas, &(Sol.X_gas())).ptr());

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

        fugcoef_liq.reset(this->fugacityCoefficient(p, T, &iliquid, &(Sol.X_liq())).ptr()); //0.79210807054894983,0.10280737903377266
        fugcoef_gas.reset(this->fugacityCoefficient(p, T, &igas, &(Sol.X_gas())).ptr());    //0.64085593096880322,0.2817704554420456
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
            fugcoef_liq.reset(this->fugacityCoefficient(p, T, &iliquid, &(Sol.X_liq())).ptr());
            fugcoef_gas.reset(this->fugacityCoefficient(p, T, &igas, &(Sol.X_gas())).ptr());

            forAll(this->X_, i)
            {
                Sol.equalconstant()[i] = fugcoef_liq()[i] / fugcoef_gas()[i];
            }

            ierr = this->solve_vaporf(p, T, Sol);
            scalar errr2 = 0.0;
            forAll(this->X_, i)
            {
                re_iter0[i] = re_iter1[i];
                fugcoef_gas()[i] *= Sol.X_gas()[i] * p;
                fugcoef_liq()[i] *= Sol.X_liq()[i] * p;
                errr2 = errr2 + sqr(fugcoef_liq()[i] / fugcoef_gas()[i] - 1.0);
                Sol.equalconstant()[i] *= (fugcoef_liq()[i] / fugcoef_gas()[i]);
                re_iter1[i] = fugcoef_liq()[i] / fugcoef_gas()[i];
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
inline Foam::scalar Foam::VLE<ThermoMixture>::alpha(
    scalar p,
    scalar T,
    solution &sol) const
{
    label pfg = 1, pfl = 0;
    scalar sv_gas = this->V(p, T, &pfg, &sol.X_gas()); //gas phase specific volume
    scalar sv_liq = this->V(p, T, &pfl, &sol.X_liq()); //liquid phase specific volume

    scalar alphagas = sol.vaporfra * sv_gas / (sol.vaporfra * sv_gas + (1.0 - sol.vaporfra) * sv_liq);
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
//-mixing rule
template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::rho(scalar p, scalar T) const
{
    autoPtr<solution> sol(TPn_flash(p, T));
    label pfg = 1, pfl = 0;
    scalar alpha = this->alpha(p, T, sol());

    scalar ZGas = this->ThermoMixture::Z(p, T, &pfg, &(sol().X_gas()));
    scalar ZLiq = this->ThermoMixture::Z(p, T, &pfl, &(sol().X_liq()));
    scalar ZMixture = ZGas * alpha + ZLiq * (1.0 - alpha);
    scalar rhoMixture = p * this->W() / (ZMixture * RR * 1.0e-03 * T);
    return rhoMixture;
}

template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::Z(scalar p, scalar T) const
{
    autoPtr<solution> sol(TPn_flash(p, T));
    label pfg = 1, pfl = 0;
    scalar alpha = this->alpha(p, T, sol());

    scalar ZGas = this->ThermoMixture::Z(p, T, &pfg, &(sol().X_gas()));
    scalar ZLiq = this->ThermoMixture::Z(p, T, &pfl, &(sol().X_liq()));
    return ZGas * alpha + ZLiq * (1.0 - alpha);
}

template <class ThermoMixture>
inline Foam::scalar Foam::VLE<ThermoMixture>::psi(
    scalar p,
    scalar T) const
{
    return rho(p, T) / p;
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
inline scalar Foam::VLE<ThermoMixture>::Ha(scalar p, scalar T) const
{
    autoPtr<solution> sol(TPn_flash(p, T));
    label pfg = 1, pfl = 0;
    scalar mw_gas = this->W(&(sol().X_gas())); //kg/mol
    scalar mw_mixture = this->W();
    scalar ygas = sol().vaporfra * mw_gas / mw_mixture;
    //Info<<"p="<<p<<",T="<<T<<",H="<<this->ThermoMixture::Ha(p, T, &pfg,&this->X_)- Hc()<<endl;
    //scalar hy1=this->ThermoMixture::Ha(p, T, &pfg,&this->X_);
    //Info<<"----------------------"<<endl;
    //Info<<"zg="<<this->ThermoMixture::Z(p, T, &pfg, &(sol().X_gas()))<<endl;
    //Info<<"fugg="<<this->fugacityCoefficient(p, T, &pfg, &(sol().X_gas()))()<<endl;
    //Info<<"zl="<<this->ThermoMixture::Z(p, T, &pfl, &(sol().X_liq()))<<endl;
    //Info<<"fugl="<<this->fugacityCoefficient(p, T, &pfl, &(sol().X_liq()))()<<endl;
    //Info<<"vaporfra="<< sol().vaporfra<<endl;
    //Info<<"----------------------"<<endl;
    //scalar hg=this->ThermoMixture::Ha(p, T, &pfg, &(sol().X_gas()));
    //scalar hl=this->ThermoMixture::Ha(p, T, &pfl, &(sol().X_liq()));
    //Info<<"--------------------\n";
    //scalar hl1=this->ThermoMixture::Ha(p, 341.064, &pfl, &(sol().X_liq()));
    //scalar hl2=this->ThermoMixture::Ha(p, 341.066, &pfl, &(sol().X_liq()));
    //Info<<"--------------------\n";
    //Info<<"-----------h1="<<hl1<<",h2="<<hl2<<endl;
    //Info<<"vaporfra="<<sol().vaporfra<<",ygas="<<ygas<<",hg="<<hg<<",hl="<<hl<<endl;
    //Info<<"X="<<sol().X_liq()<<endl;
    //Info<<"Ha="<<this->ThermoMixture::Ha(p, T, &pfg, &(sol().X_gas())) * ygas + this->ThermoMixture::Ha(p, T, &pfl, &(sol().X_liq())) * (1.0 - ygas)<<endl;
    //scalar ht=this->ThermoMixture::Ha(p, T, &pfg, &(sol().X_gas())) * ygas + this->ThermoMixture::Ha(p, T, &pfl, &(sol().X_liq())) * (1.0 - ygas);
    return this->ThermoMixture::Ha(p, T, &pfg, &(sol().X_gas())) * ygas + this->ThermoMixture::Ha(p, T, &pfl, &(sol().X_liq())) * (1.0 - ygas);
    //return this->ThermoMixture::Ha(p, T, &pfg,&this->X_);
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
    /*
    autoPtr<solution> sol(TPn_flash(p, T));
    label pfg = 1, pfl = 0;
    //scalar alpha = this->alpha(p, T, sol());
    scalar mw_gas = this->W(&(sol().X_gas())); //kg/mol
    scalar mw_mixture = this->W();
    scalar ygas = sol().vaporfra * mw_gas / mw_mixture;
    //Info<<this->ThermoMixture::Cp(p, T, &pfg, &(sol().X_gas()));
    return this->ThermoMixture::Cp(p, T, &pfg, &(sol().X_gas())) * ygas + this->ThermoMixture::Cp(p, T, &pfl, &(sol().X_liq())) * (1.0 - ygas);
    */
    scalar dt = 0.01;
    return (Hs(p, T + dt) - Hs(p, T)) / dt;
}

template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::CpMCv(scalar p, scalar T) const
{
    autoPtr<solution> sol(TPn_flash(p, T));
    label pfg = 1, pfl = 0;
    scalar mw_gas = this->W(&(sol().X_gas())); //kg/mol
    scalar mw_mixture = this->W();
    scalar ygas = sol().vaporfra * mw_gas / mw_mixture;
    return this->ThermoMixture::CpMCv(p, T, &pfg, &(sol().X_gas())) * ygas + this->ThermoMixture::CpMCv(p, T, &pfl, &(sol().X_liq())) * (1.0 - ygas);
}

template <class ThermoMixture>
inline scalar Foam::VLE<ThermoMixture>::S(scalar p, scalar T) const
{
    autoPtr<solution> sol(TPn_flash(p, T));
    label pfg = 1, pfl = 0;
    scalar mw_gas = this->W(&(sol().X_gas())); //kg/mol
    scalar mw_mixture = this->W();
    scalar ygas = sol().vaporfra * mw_gas / mw_mixture;
    return this->ThermoMixture::S(p, T, &pfg, &(sol().X_gas())) * ygas + this->ThermoMixture::S(p, T, &pfl, &(sol().X_liq())) * (1.0 - ygas);
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
    const VLE<ThermoMixture> &a,
    const VLE<ThermoMixture> &b)
{
    return VLE<ThermoMixture>(a);
}

template <class ThermoMixture>
Foam::VLE<ThermoMixture> Foam::operator*(
    const scalar s,
    const VLE<ThermoMixture> &b)
{
    return VLE<ThermoMixture>(b);
}

// ************************************************************************* //
