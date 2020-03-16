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

#include "PengRobinsonMix.H"
#include "specie.H"
#include "mathematicalConstants.H"
#include "OFstream.H"
#include "IFstream.H"
#include "string.H"
#define PRMIX


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::PengRobinsonMix<Specie>::PengRobinsonMix
(
    const Specie& sp,
    const scalar& Tc,
    const scalar& Vc,
    const scalar& Zc,
    const scalar& Pc,
    const scalar& omega,
    const scalar& mwsp,

    const scalar& Tcc12,
    const scalar& Pcc12,
    const scalar& Vcc12,
    const scalar& omegac12,
    const scalar& mwspc12,

    const scalar& Tcc7,
    const scalar& Pcc7,
    const scalar& Vcc7,
    const scalar& omegac7,
    const scalar& mwspc7,

    const scalar& Tcco2,
    const scalar& Pcco2,
    const scalar& Vcco2,
    const scalar& omegaco2,
    const scalar& mwspco2,

    const scalar& Tcch4,
    const scalar& Pcch4,
    const scalar& Vcch4,
    const scalar& omegach4,
    const scalar& mwspch4,

    const scalar& Tco2,
    const scalar& Pco2,
    const scalar& Vco2,
    const scalar& omegao2,
    const scalar& mwspo2,

    const scalar& Tch2o,
    const scalar& Pch2o,
    const scalar& Vch2o,
    const scalar& omegah2o,
    const scalar& mwsph2o,

    const scalar& TcH2,
    const scalar& PcH2,
    const scalar& VcH2,
    const scalar& omegaH2,
    const scalar& mwspH2,

    const scalar& TcH2S,
    const scalar& PcH2S,
    const scalar& VcH2S,
    const scalar& omegaH2S,
    const scalar& mwspH2S,
    
    const scalar& system,
    const scalar& bipc12n2,   
    const scalar& bipco2ch4,   
    const scalar& bipco2h2o,   
    const scalar& bipch4o2,  
    const scalar& bipc7n2 
)
:
    Specie(sp), //YP: initial list
    Tc_(Tc),
    Vc_(Vc),
    Zc_(Zc),
    Pc_(Pc),
    omega_(omega),
    mwsp_(mwsp),

    Tcc12_(Tcc12),
    Pcc12_(Pcc12),
    Vcc12_(Vcc12),
    omegac12_(omegac12),
    mwspc12_(mwspc12),

    Tcc7_(Tcc7),
    Pcc7_(Pcc7),
    Vcc7_(Vcc7),
    omegac7_(omegac7),
    mwspc7_(mwspc7),

    Tcco2_(Tcco2),
    Pcco2_(Pcco2),
    Vcco2_(Vcco2),
   omegaco2_(omegaco2),
    mwspco2_(mwspco2),

    Tcch4_(Tcch4),
    Pcch4_(Pcch4),
    Vcch4_(Vcch4),
    omegach4_(omegach4),
    mwspch4_(mwspch4),

    Tco2_(Tco2),
    Pco2_(Pco2),
    Vco2_(Vco2),
    omegao2_(omegao2),
    mwspo2_(mwspo2),

    Tch2o_(TcH2),
    Pch2o_(PcH2),
    Vch2o_(VcH2),
    omegah2o_(omegaH2),
    mwsph2o_(mwspH2),

    TcH2_(TcH2),
    PcH2_(PcH2),
    VcH2_(VcH2),
    omegaH2_(omegaH2),
    mwspH2_(mwspH2),

    TcH2S_(TcH2S),
    PcH2S_(PcH2S),
    VcH2S_(VcH2S),
    omegaH2S_(omegaH2S),
    mwspH2S_(mwspH2S),
    
    system_(system),
    bipc12n2_(bipc12n2),
    bipco2ch4_(bipco2ch4),
    bipco2h2o_(bipco2h2o),
    bipch4o2_(bipch4o2),
    bipc7n2_(bipc7n2) 

{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::PengRobinsonMix<Specie>::PengRobinsonMix 
(
    const word& name,
    const PengRobinsonMix& pg 
)
:
    Specie(name, pg),
    Tc_(pg.Tc_),
    Vc_(pg.Vc_),
    Zc_(pg.Zc_),
    Pc_(pg.Pc_),
    omega_(pg.omega_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::PengRobinsonMix <Specie>>
Foam::PengRobinsonMix<Specie>::clone() const
{
    return autoPtr<PengRobinsonMix<Specie>>
    (
        new PengRobinsonMix<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::PengRobinsonMix<Specie>>
Foam::PengRobinsonMix<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<PengRobinsonMix<Specie>>
    (
        new PengRobinsonMix<Specie>(dict)
    );
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::rho
(
    scalar p,
    scalar T
) const
{
    const scalar Z = this->Z(p, T);
    return p/(Z*this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::H(scalar p, scalar T) const
{
    const scalar Pr = p/Pc_;
    const scalar Tr = T/Tc_;
    const scalar B = 0.07780*Pr/Tr;
    const scalar kappa = 0.37464 + 1.54226*omega_ - 0.26992*sqr(omega_);
    const scalar alpha = sqr(1 + kappa*(1 - sqrt(Tr)));

    const scalar Z = this->Z(p, T);

    return
        this->R()
       *Tc_
       *(
           Tr*(Z - 1)
         - 2.078*(1 + kappa)*sqrt(alpha)
          *log((Z + 2.414*B)/(Z - 0.414*B))
        );
}


template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::Cp(scalar p, scalar T) const
{
    const scalar Tr = T/Tc_;
    const scalar a = 0.45724*sqr(RR*Tc_)/Pc_;
    const scalar b = 0.07780*RR*Tc_/Pc_;
    const scalar kappa = 0.37464 + 1.54226*omega_ - 0.26992*sqr(omega_);
    const scalar alpha = sqr(1 + kappa*(1 - sqrt(Tr)));

    const scalar A = a*alpha*p/sqr(RR*T);
    const scalar B = b*p/(RR*T);

    const scalar Z = this->Z(p, T);

    const scalar ap = kappa*a*(kappa/Tc_ - (1 + kappa)/sqrt(T*Tc_));
    const scalar app = kappa*a*(1 + kappa)/(2*sqrt(pow3(T)*Tc_));

    const scalar M = (sqr(Z) + 2*B*Z - sqr(B))/(Z - B);
    const scalar N = ap*B/(b*RR);

    const scalar root2 = sqrt(2.0);

    return
    (
        app*(T/(2*root2*b))*log((Z + (root2 + 1)*B)/(Z - (root2 - 1)*B))
      + RR*sqr(M - N)/(sqr(M) - 2*A*(Z + B))
      - RR
    )/this->W();
}


template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::S
(
    scalar p,
    scalar T
) const
{
    const scalar Pr = p/Pc_;
    const scalar Tr = T/Tc_;
    const scalar B = 0.07780*Pr/Tr;
    const scalar kappa = 0.37464 + 1.54226*omega_ - 0.26992*sqr(omega_);

    const scalar Z = this->Z(p, T);

    return
        this->R()
       *(
          - log(p/Pstd)
          + (
                log(Z - B)
              - 2.078*kappa*((1 + kappa)/sqrt(Tr) - kappa)
               *log((Z + 2.414*B)/(Z - 0.414*B))
            )
        );
}


template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::psi
(
    scalar p,
    scalar T
) const
{
    const scalar Z = this->Z(p, T);

    return 1.0/(Z*this->R()*T);
}


template<class Specie> 
inline Foam::scalar Foam::PengRobinsonMix<Specie>::Z
(
    scalar p,
    scalar T
) const
{
    const scalar Tr = T/Tc_;
    const scalar a = 0.45724*sqr(RR*Tc_)/Pc_;
    const scalar b = 0.07780*RR*Tc_/Pc_;
    const scalar kappa = 0.37464 + 1.54226*omega_ - 0.26992*sqr(omega_);
    const scalar alpha = sqr(1 + kappa*(1 - sqrt(Tr)));

    const scalar A = a*alpha*p/sqr(RR*T);
    const scalar B = b*p/(RR*T);

    const scalar a2 = B - 1;
    const scalar a1 = A - 2*B - 3*sqr(B);
    const scalar a0 = -A*B + sqr(B) + pow3(B);

    const scalar Q = (3*a1 - a2*a2)/9.0;
    const scalar Rl = (9*a2*a1 - 27*a0 - 2*a2*a2*a2)/54;

    const scalar Q3 = Q*Q*Q;
    const scalar D = Q3 + Rl*Rl;

    scalar root = -1;

    if (D <= 0)
    {
        const scalar th = ::acos(Rl/sqrt(-Q3));
        const scalar qm = 2*sqrt(-Q);
        const scalar r1 = qm*cos(th/3) - a2/3;
        const scalar r2 =
            qm*cos((th + 2*constant::mathematical::pi)/3) - a2/3;
        const scalar r3 =
            qm*cos((th + 4*constant::mathematical::pi)/3) - a2/3;

        root = max(r1, max(r2, r3));
    }
    else
    {
        // One root is real
        const scalar D05 = sqrt(D);
        const scalar S = pow(Rl + D05, 1.0/3.0);
        scalar Tl = 0;
        if (D05 > Rl)
        {
            Tl = -pow(mag(Rl - D05), 1.0/3.0);
        }
        else
        {
            Tl = pow(Rl - D05, 1.0/3.0);
        }

        root = S + Tl - a2/3;
    }

    return root;
}


// ======================================================================//
// Ideal properties
// History :   Creation.	     P. YI			2019     //
// ======================================================================//
template<class Specie>
inline void Foam::PengRobinsonMix<Specie>:: initialization //PY: initialize Yi, species mass fraction
(
    scalar p,   //inp
    scalar T    //inp
) const
{
    scalar comp[num_sp] = {0.0};
    if (TPnFlag == 1)
    {
	scalar comp_gas[num_sp] = {0.0};
	scalar comp_liq[num_sp] = {0.0};
	this->initPTXdiagram(p, T, comp, comp_gas, comp_liq);
    }
    else if(TabFlag == 1)
    {
	this->genTabulation(p, T, comp);    
    }	
    else if(TabFlag == 2) //For binary diffusion analyzation
    {
	this->genDiffusion(p, T, comp);    
    }	    
    else if(HnFlag == 1)
    {
	this->Hn_flash(p, T, comp);  
    }	
    else if(TPDFlag == 1)
    {
	scalar tpdtest;
	scalar comp_gas[num_sp] = {0.0};
	scalar comp_liq[num_sp] = {0.0};
	this->initTPD(p, T, comp, comp_gas, comp_liq, tpdtest); 
    }	
    Info<< "Initialization_Thermo_done" << "###### @_@ ###### ^_^ ###### $_$ ######" << endl;
}


double phasestate(double vaporfra)
{
    if (vaporfra > 1.0e-06 && vaporfra < 0.999999)
	{
		return  2.0;
	}
	else if (vaporfra <= 1.0e-06)
	{
		return  1.0;
	}
	else if (vaporfra >= 0.999999)
	{
		return 0.0;
	}
    return  0;
}

// ======================================================================//
// Initilization for TP and PX diagrams					 //
// Input: pressure, temperature, feed					 //
// Output: phase compositions						 //
//									 //
// History :   Creation.	     P. YI		2019.04		 //
// ======================================================================//
template<class Specie>
inline void Foam::PengRobinsonMix<Specie>:: initPTXdiagram
(
    scalar press, //inp
    scalar temp,  //inp
    scalar *comp, //inp
    scalar *comp_liq, //outp
    scalar *comp_gas  //outp
) const
{
    scalar vaporfra = 0.9;
    scalar equalconstant[num_sp] = {0.0};

    scalar tpdtest = 0.0;
    scalar tempphase_pt[10] = {0.0};
    scalar presphase = 0.0;
    scalar tpdlast = 1.0;
    scalar dettpd = 0.0;
    int ndettpd = 0;

    ofstream output;
    if (flag_px == 1)
    {
	output.open("PXdiagram.txt");
    }
    if (flag_px == 2)
    {
	output.open("TXdiagram.txt");
    }
    if (flag_tp == 1)
    {
	output.open("PTdiagram.txt");
    }

    //for (label ncomp = 0; ncomp < num_comp; ncomp++)
    //{
	if (num_sp == 2)//binary:which is used to calcualte the TX diagram
	{        
	    comp[0] = compinp1; 
            comp[1] = compinp2;
	}
	else if (num_sp == 3)
	{
    	    comp[0] = compinp1;
	    comp[1] = compinp2;
            comp[2] = compinp3;
	}
	else if (num_sp == 4)
	{
    	    comp[0] = compinp1;
	    comp[1] = compinp2;
            comp[2] = compinp3;
	    comp[3] = compinp4;
	}
	else if (num_sp == 5)
	{
    	    comp[0] = compinp1;
	    comp[1] = compinp2;
            comp[2] = compinp3;
	    comp[3] = compinp4;
	    comp[4] = compinp5;
	}  

	for (label npres = 0; npres < num_pres; npres++)
	{
	    tpdlast = 1.0;
	    press = press0 + (npres - 0)*det_pres;
	    ndettpd = 0;

	    for (label ntemp = 0; ntemp < num_temp; ntemp++)
	    {
		temp = temp0 + (ntemp - 0)*det_temp;
		Info<<"PY:Screen_TPnFlash"<<", z="<<comp[0]<<", p(bar)="<<press*1.e-05<<", temp="<<temp<<" ^_^"<< endl;
		this->TPn_flash(press, temp, comp, comp_liq, comp_gas, vaporfra, equalconstant);
		//alphagas = this->Evaluate_alpha(press, temp, vaporfra, comp_liq, comp_gas,comp);

		scalar detCompPhase = 0.0;
		if (flag_px == 1)
		{
		    if (system_ == 3)
		    {
		    output<<comp_liq[0]<<' '<<' '<<' '<<comp_gas[0]<<' '<<' '<<' '<<press*1.0e-06; 
		    output<< "\n";
		    }
		    else
		    {
		    output<<comp_liq[1]<<' '<<' '<<' '<<comp_gas[1]<<' '<<' '<<' '<<press*1.0e-06; 
		    output<< "\n";
		    }
		}

		if (flag_px == 2) //TX diagram
		{
		    if (system_ == 1 || system_ == 5)
		    {
		    	detCompPhase = comp_liq[1] - comp_gas[1];
		        if (detCompPhase > 1.0e-04 || detCompPhase < -1.0e-04)
		        {
			    output<<comp_liq[1]<<' '<<' '<<' '<<comp_gas[1]<<' '<<' '<<' '<<temp; 
			    output<< "\n";
		        }
		    }
		    else if (system_ != 1 || system_ != 5)
		    {
		    	detCompPhase = comp_liq[0] - comp_gas[0];
		        if (detCompPhase > 1.0e-04 || detCompPhase < -1.0e-04)
		        {
			    output<<comp_liq[1]<<' '<<' '<<' '<<comp_gas[1]<<' '<<' '<<' '<<temp; 
			    output<< "\n";
		        }
		    }
		    //printf("Screen_observe2_detCompPhase %5.10f \n",detCompPhase); 
		}

		if (flag_tp == 1)
		{
		    if (vaporfra > 1.0e-06 && vaporfra < 0.999999)
		    {
			tpdtest = 2.0;
		    }
		    else if (vaporfra <= 1.0e-06)
		    {
			tpdtest = 1.0;
		    }
		    else if (vaporfra >= 0.999999)
		    {
			tpdtest = 0.0;
		    }
		    dettpd = tpdtest - tpdlast; //if tpd change first time, abs(dettpd) > 0;
		    scalar maxtpd = max(tpdtest,tpdlast);
            /*
		    if (dettpd != 0.0 && maxtpd == 2.0)
		    {
			ndettpd = ndettpd + 1;
			int ndettpd0 = ndettpd-1;
			for (label nwrite = ndettpd0; nwrite < ndettpd; nwrite++)
			{
			    tempphase_pt[nwrite] = temp;
			}
		    }
            */
           if (dettpd != 0.0)
           {
               double templ=temp-det_temp;
               double tempr=temp;
               double tempm=(templ+tempr)/2;
               double tpdtest_t=0;
               while(tempr-templ>1e-3)
               {
                   tempm=(templ+tempr)/2;
                   this->TPn_flash(press, tempm, comp, comp_liq, comp_gas, vaporfra, equalconstant);
                   tpdtest_t  = phasestate(vaporfra);
                   if(tpdtest_t ==tpdlast)
                   {
                       templ=tempm;
                   }
                   if (tpdtest_t == tpdtest)
                   {
                       tempr=tempm;
                   }
               }                                                              
               tempphase_pt[ndettpd++] = tempm;
           }
		    presphase = press;

		    //if (dettpd != 0.0 && ndettpd == 2)
            if (ntemp == num_temp-1)
		    {
            for(int i=0;i<ndettpd;i++)
			output<<tempphase_pt[i]<<' '<<' '<<' ';
            output<<presphase*1.0e-05; 
			output<< "\n";
		    }
		    tpdlast = tpdtest;
		}

		if (flag_tp == 2)
		{
		    detCompPhase = comp_liq[1] - comp_gas[1];
		    if ((detCompPhase > 1.0e-04 || detCompPhase < -1.0e-04) && vaporfra <= 1.0e-06)
		    {
			tpdtest = 1.0;
		    }

		    else if ((detCompPhase > 1.0e-04 || detCompPhase < -1.0e-04) && vaporfra >= 0.999999)
		    {
			tpdtest = 0.0;
		    }
		    else 
		    {
			tpdtest = 2.0;
		    }

		    dettpd = tpdtest - tpdlast; //if tpd change first time, abs(dettpd) > 0;
		    scalar maxtpd = max(tpdtest,tpdlast);
		    if (dettpd != 0.0 && maxtpd == 2.0)
		    {
			ndettpd = ndettpd + 1;
			int ndettpd0 = ndettpd-1;
			for (label nwrite = ndettpd0; nwrite < ndettpd; nwrite++)
			{
			    tempphase_pt[nwrite] = temp;
			}
		    }

		    presphase = press;
		    if (dettpd != 0.0 && ndettpd == 2)
		    {
			output<<tempphase_pt[0]<<' '<<' '<<' '<<tempphase_pt[1]<<' '<<' '<<' '<<presphase*1.0e-05; 
			output<< "\n";
		    }
		    tpdlast = tpdtest;
		}
	        Info<<"PY:Screen_TPnFlash"<<", x="<<comp_liq[0]<<", y="<<comp_gas[0]<<", vap="<<vaporfra<<"\n"<< endl;
	    }
	}
    //}
    output.close();
}


// ======================================================================//
// Generate thermo table for tabulation1 				 //
// Input: pressure, temperature, feed					 //
// Output: thermo and transport properties 				 //
//									 //
// History :   Creation.	     P. YI		2019.05		 //
// ======================================================================//
template<class Specie>
inline void Foam::PengRobinsonMix<Specie>::genTabulation
(
    scalar press, //inp
    scalar temp,  //inp
    scalar *comp //inp
) const
{
    string ThermoTable;    
    scalar vaporfra = 0.5;
    scalar equalconstant[num_sp] = {0.0};

    scalar comp_liq[num_sp] = {0.0};
    scalar comp_gas[num_sp] = {0.0};
    scalar comp_liq0[num_sp] = {0.0};
    scalar comp_gas0[num_sp] = {0.0};

    ofstream output;
    output.open("ThermoTable.txt");
    if (num_sp == 2)
    {        
	comp[0] = compinp1;
        comp[1] = compinp2;
    }
    else if (num_sp == 3)
    {
    	comp[0] = compinp1;
	comp[1] = compinp2;
        comp[2] = compinp3;
    }
    else if (num_sp == 4)
    {
    	comp[0] = compinp1;
	comp[1] = compinp2;
        comp[2] = compinp3;
	comp[3] = compinp4;
    }
    else if (num_sp == 5)
    {
    	comp[0] = compinp1;
	comp[1] = compinp2;
        comp[2] = compinp3;
	comp[3] = compinp4;
	comp[4] = compinp5;
    }  
    for (label isp = 0; isp < num_sp; isp++)
    {
	comp_liq0[isp] = comp[isp];
	comp_gas0[isp] = comp[isp];
    }

    scalar x_fuel  = 0.0;
    for (label ncomp = 0; ncomp < num_comp; ncomp++)
    {
	x_fuel = x_fuel0 + (ncomp - 0)*det_com1;
	if (num_sp == 2) //C7n2 OR c12n2
	{
	    comp[0] = x_fuel;
	    comp[1] = 1.0-comp[0];
	}
	else if (num_sp == 3) //CO2ch4H2O 
	{
	    comp[1] = x_fuel;           //CH4
	    comp[0] = x_fuel;           //CO2
	    comp[2] = 1.0 - 2.0*x_fuel; //H2O
	}	
	else if (num_sp == 4) //CO2ch4o2H2O
	{
	    comp[1] = x_fuel;//CH4
	    comp[0] = (1.0-comp[1])*1.0/6.0;//CO2
	    comp[2] = (1.0-comp[1])*1.0/6.0;//O2
	    comp[3] = (1.0-comp[1])*4.0/6.0;//H2O
	}
	else if (num_sp == 5) //CO2ch4n2o2H2O
	{
	    comp[1] = x_fuel;//CH4
	    comp[0] = (1.0-comp[1])*5.0/10.1;//CO2=O2
	    comp[2] = (1.0-comp[1])*1.0/10.1;
	    comp[3] = (1.0-comp[1])*4.0/10.1;
	    comp[4] = (1.0-comp[1])*0.1/10.1;//H2O=O2
	}
	
	scalar mw_mixture = 0.0;  //kg/mol
	scalar rhoMixture = 0.0;  //kg/m3
	scalar sieMixture = 0.0;  //J/kg
	scalar CpMixture  = 0.0;  //J/kgK
	scalar CvMixture  = 0.0;  //J/kgK
	scalar CsMixture  = 0.0;  //cm/s	
	scalar ZMixture   = 0.0;  //[-]
	//scalar ZMixture2  = 0.0;  //[-]
	//scalar GammaMixture = 0.0; //[-]	
		
	scalar kappaMixture     = 0.0; 
	scalar muMixture        = 0.0; 
	scalar rho_gas          = 0.0,          rho_liq = 0.0; 
	scalar Cp_gas           = 0.0,          Cp_liq = 0.0; 
	scalar Cv_gas           = 0.0,          Cv_liq = 0.0; 
	scalar alphagas         = 0.0,          ygas = 0.0;
	scalar kappa_gas        = 0.0,          kappa_liq = 0.0;
	scalar mu_gas           = 0.0,          mu_liq = 0.0;
	scalar Ct_gas           = 0.0,          Ct_liq = 0.0;
	scalar Cs_gas           = 0.0,          Cs_liq = 0.0;
	scalar VspecificGas     = 0.0,          VspecificLiq = 0.0;
	scalar ZGas             = 0.0,          ZLiq = 0.0;
        scalar Dij_binary_high  = 0.0;	
	
	for (label npres = 0; npres < num_pres; npres++)
	{
	    press = press0 + (npres - 0)*det_pres;
	    
	    for (label ntemp = 0; ntemp < num_temp; ntemp++)
	    {
		temp = temp0 + (ntemp - 0)*det_temp;
		Info<<"PY:Screen_Tabulation"<<", z="<<comp[0]<<", p(bar)="<<press*1.e-05<<", temp="<<temp<<" ^_^XYY"<< endl;
		
		if (PCFlag == 0)
		{
		    for (label isp = 0; isp < num_sp; isp++)
		    {
			comp_liq[isp] = comp_liq0[isp];
			comp_gas[isp] = comp_gas0[isp];
		    }

		    this->TPn_flash(press, temp, comp, comp_liq, comp_gas, vaporfra, equalconstant);
//Info<<"PY:Screen_1"<< ", temp="<<temp<<", comp_liq="<<comp_liq[0]<<", comp_gas[0]="<<comp_gas[0]<<"\n"<< endl;		    
		    for (label isp = 0; isp < num_sp; isp++)
		    {
			comp_liq0[isp] = comp_liq[isp];
			comp_gas0[isp] = comp_gas[isp];
		    }

	            if (vaporfra > 1.0)
	            {
	                vaporfra = 0.9999999;
	            }
	            else if (vaporfra < 0.0)
	            {
	                vaporfra = 1.0e-07;
	            }
	            alphagas		= this->Evaluate_alpha(press, temp, vaporfra, comp_liq, comp_gas,comp);
		    scalar mw_gas	= this->mwmix(comp_gas);//kg/mol
		    //mw_mixture	        = this->mwmix(comp); 		    	    		    
		    
		    scalar mw_liq	= this->mwmix(comp_liq);//kg/mol	
		    mw_mixture	        = mw_gas*vaporfra + mw_liq*(1.0-vaporfra);
		    ygas		= vaporfra*mw_gas/mw_mixture;
		    
                    //Thermal properties
		    VspecificGas        = this->volmmix_phase(0, press, temp, comp_gas);//m3/mol
		    VspecificLiq        = this->volmmix_phase(1, press, temp, comp_liq);	
		    ZGas                = press*VspecificGas/(RR*1.0e-03*temp);
		    ZLiq                = press*VspecificLiq/(RR*1.0e-03*temp);			    
		    ZMixture	        = ZGas*alphagas + ZLiq*(1.0 - alphagas);
		    rhoMixture		= press*mw_mixture/(ZMixture*RR*1.0e-03*temp);	
		            
		    rho_gas		= this->rhomix_phase(0, press, temp, comp_gas);  //kg/m3
		    rho_liq		= this->rhomix_phase(1, press, temp, comp_liq); 
		    //rhoMixture	        = rho_gas*alphagas + rho_liq*(1.0 - alphagas);	        

		    scalar sieReal_gas	= this->siemix_phase(0, press, temp, comp_gas); //J/kg
		    scalar sieReal_liq	= this->siemix_phase(1, press, temp, comp_liq);
		    sieMixture		= sieReal_gas*ygas + sieReal_liq*(1.0 - ygas);
		    
		    Cp_gas		= this->cpmix_phase(0, press, temp, comp_gas); //J/kgK
		    Cp_liq		= this->cpmix_phase(1, press, temp, comp_liq); 
		    CpMixture		= Cp_gas*ygas + Cp_liq*(1.0 - ygas);

		    Cv_gas		= this->cvmix_phase(0, press, temp, comp_gas); //J/kgK
		    Cv_liq		= this->cvmix_phase(1, press, temp, comp_liq); 
		    CvMixture		= Cv_gas*ygas + Cv_liq*(1.0 - ygas);			     		    
		    //GammaMixture        = CpMixture/CvMixture;

                                 this->soundspeedmix(0, press, temp, comp_gas, Ct_gas, Cs_gas);
                                 this->soundspeedmix(1, press, temp, comp_liq, Ct_liq, Cs_liq);  
 
                    CsMixture           = sqrt(1.0/(rhoMixture*(alphagas/(rho_gas*sqr(Cs_gas)) + 
                                          (1.0-alphagas)/(rho_liq*sqr(Cs_liq))))); 
                    //CtMixture           = CsMixture * sqr(GammaMixture);	            

                    //Transport	Properties    
		                          this->kappa_phase(0, press, temp, comp_gas, kappa_gas, mu_gas); //W/mK 		    
		                          this->kappa_phase(1, press, temp, comp_liq, kappa_liq, mu_liq); 
		    kappaMixture	= kappa_gas*ygas + kappa_liq*(1.0 - ygas); //W/mK 

		    muMixture	        = 1.0e+06*(mu_gas*ygas + mu_liq*(1.0 - ygas)); //up 

		    //Dij_binary_high	= this->Dij_highP(press, temp, comp);//cm2/s

		    //ZMixture2	        = press*mw_mixture/(rhoMixture*RR*1.0e-03*temp);
		}
		else
		{
		    rhoMixture		= this->rhomix_phase(0, press, temp, comp); 
		    sieMixture		= this->siemix_phase(0, press, temp, comp); 
		    CpMixture		= this->cpmix_phase(0, press, temp, comp); 
		    CvMixture		= this->cvmix_phase(0, press, temp, comp); 		    
		    this->kappa_phase(0, press, temp, comp, kappaMixture, muMixture); 
		    
		    //Dij_binary_high	= this->Dij_highP_pesudo(press, temp, comp);
		}

	        //output<<temp<<' '<<' '<<press*1.0e-05<<' '<<' '<<//comp[0]<<' '<<' '<<
	        
	                //rhoMixture;<<' '<<' '<<CpMixture<<' '<<' '<<
	                
	                //kappaMixture;//<<' '<<' '<<muMixture<<' '<<' '<<
	                
	                //ZMixture;<<' '<<' '<<Dij_binary_high*1.0e+06; 

	        output<<temp<<' '<<' '<<press<<' '<<' '<<comp[0]<<' '<<' '<<
	        
	                rhoMixture<<' '<<' '<<sieMixture<<' '<<' '<<
	                
	                CpMixture<<' '<<' '<<CvMixture<<' '<<' '<<CsMixture<<' '<<' '<<
	                
	                kappaMixture<<' '<<' '<<muMixture<<' '<<' '<<
	                
	                ZMixture<<' '<<' '<<Dij_binary_high<<' '<<' '<<
	                
	                comp_liq[0]<<' '<<' '<<comp_gas[0]<<' '<<' '<<alphagas; 
                    
	        output<< "\n";
	        //Info<<"PY:Screen_Tabulation"<<", kappa_gas="<<kappa_gas<<", ygas="<<ygas<<", kappa_liq="<<kappa_liq<<", kappaMixture="<<kappaMixture<<"\n"<< endl;
	    }
	}
    }
    output.close();
}



// ======================================================================//
// Generate binary mass diffusion          				 //
// Input: pressure, temperature, feed					 //
// Output: mass diffusion                				 //
//									 //
// History :   Creation.	     P. YI		2019.07		 //
// ======================================================================//
template<class Specie>
inline void Foam::PengRobinsonMix<Specie>::genDiffusion
(
    scalar press, //inp
    scalar temp,  //inp
    scalar *comp  //inp
) const
{
    string DijTable;    
    scalar vaporfra = 0.5;
    scalar equalconstant[num_sp] = {0.0};

    scalar comp_liq[num_sp] = {0.0};
    scalar comp_gas[num_sp] = {0.0};
    scalar comp_liq0[num_sp] = {0.0};
    scalar comp_gas0[num_sp] = {0.0};

    ofstream output;
    output.open("DijTable.txt");
    for (label isp = 0; isp < num_sp; isp++)
    {
	comp_liq0[isp] = comp[isp];
	comp_gas0[isp] = comp[isp];
    }

    scalar x_fuel0 = 0.5;//0.999;
    scalar x_fuel  = 0.0;
    for (label ncomp = 0; ncomp < num_comp; ncomp++)
    {
	x_fuel = x_fuel0 + (ncomp - 0)*det_com1;
	comp[0] = x_fuel;
        comp[1] = 1.0-comp[0];
	
	scalar Dij_binary_low = 0.0;
        scalar Dij_binary_high = 0.0;
        scalar rho_reduced = 0.0;
        scalar Vcij = 0.0, Vmol = 0.0;
        scalar wdd = 0.0;        
	for (label npres = 0; npres < num_pres; npres ++)
	{
	    press = press0 + (npres - 0)*det_pres;
	    temp = temp0;
	    printf("Screen_tabulation_z,p,T %5.3f %5.3f %5.3f \n",comp[0],press*1.e-05,temp); 
	    for (label ntemp = 0; ntemp < num_temp; ntemp++)
	    {
		temp = temp0 + (ntemp - 0)*det_temp;
		for (label isp = 0; isp < num_sp; isp++)
		{
	            comp_liq[isp] = comp_liq0[isp];
		    comp_gas[isp] = comp_gas0[isp];
		}
		this->TPn_flash(press, temp, comp, comp_liq, comp_gas, vaporfra, equalconstant);
		for (label isp = 0; isp < num_sp; isp++)
		{
		    comp_liq0[isp] = comp_liq[isp];
		    comp_gas0[isp] = comp_gas[isp];
		}

		if (PCFlag == 0)
		{		    
		    //Ping's C.P.
                    scalar Tcij = this->TcBi(comp);
                    scalar Pcij = this->PcBi(comp);
                    Vcij = 1.0e+06*this->volmmix_phase(0, Pcij, Tcij, comp);
                    Vmol = 1.0e+06*this->volmmix_phase(0, press, temp, comp);
                    // Reduced density
                    rho_reduced = Vcij/Vmol;
		    
		    Dij_binary_low	= this->Dij_lowP(press, temp, comp);
		    Dij_binary_high	= this->Dij_highP(press, temp, comp);
		    wdd                 = Dij_binary_high/Dij_binary_low;
		    printf("Diffusion_low_high_Ping= %5.6f %5.6f %5.6f %5.6f\n\n",
                        rho_reduced, wdd, Vcij, Vmol);  
		}
		else
		{
                    scalar Vc_sp3[num_sp] = {0.0};
                    for (label spid = 0; spid < num_sp; spid++)
                    {
	                Vc_sp3[spid] = Vc_sp[spid]*1.0e+03;//cm3/mol
                    }
                    scalar sumVij = pow(Vc_sp3[0],1.0/3) + pow(Vc_sp3[1],1.0/3);
                    Vcij = pow(sumVij,3.0)/8.0;
                    Vmol = RR*1.0e+03*temp/press;//cm3/mol
                    rho_reduced = Vcij/Vmol;

		    //Pesudo C.P.
		    Dij_binary_low	= this->Dij_lowP_pesudo(press, temp, comp);
		    Dij_binary_high	= this->Dij_highP_pesudo(press, temp, comp);
		    wdd                 = Dij_binary_high/Dij_binary_low;		    
		    printf("Diffusion_low_high_Bellan= %5.6f %5.6f %5.6f %5.6f\n\n",
                        rho_reduced, wdd, Vcij, Vmol);  
		}
	        output<<temp<<' '<<' '<<press*1.0e-05<<' '<<' '<<rho_reduced<<' '<<' '<<Dij_binary_low<<' '<<' '<<Dij_binary_high*1.0e+03;
	        //output<<press*1.0e-05<<' '<<' '<<rho_reduced<<' '<<' '<<wdd<<' '<<' '<<Dij_binary_low*1.0e+03; 
	        //output<<press*1.0e-05<<' '<<' '<<wdd; 	
	                 
	        output<< "\n";
	    }
	}
    }
    output.close();
}



// ======================================================================//
// Hn Flash				 				 //
// Input: pressure, enthalpy and feed					 //
// Output: Phase compostions and temperature 				 //
//									 //
// History :   Creation.	     P. YI		2019.06		 //
// ======================================================================//
template<class Specie>
inline void Foam::PengRobinsonMix<Specie>::Hn_flash // PY:with a course initialization
(
    scalar press,	//input
    scalar temp,	//input
    scalar *comp	//input
) const
{
    string HTxTable;    
    scalar vaporfra = 0.5;
    scalar equalconstant[num_sp] = {0.0};

    scalar comp_liq[num_sp] = {0.0};
    scalar comp_gas[num_sp] = {0.0};

    ofstream output;
    output.open("HTxTable.txt");

    if (num_sp == 2)
    {        
	comp[0] = compinp1;
	comp[1] = compinp2;
    }
    else if (num_sp == 3)
    {
	comp[0] = compinp1;
	comp[1] = compinp2;
	comp[2] = compinp3;
    }
    else if (num_sp == 4)
    {
	comp[0] = compinp1;
	comp[1] = compinp2;
	comp[2] = compinp3;
	comp[3] = compinp4;
    }
    else if (num_sp == 5)
    {
	comp[0] = compinp1;
	comp[1] = compinp2;
	comp[2] = compinp3;
	comp[3] = compinp4;
	comp[4] = compinp5;
    }  

    scalar comp_liq0[num_sp] = {0.0}, comp_gas0[num_sp] = {0.0};

    scalar templiq = 300.0, tempgas = 1100;

    for (label npres = 0; npres < num_pres; npres++)
    {
	press = press0 + (npres - 0)*det_pres;
	//this->TPn_flash(press, templiq, comp, comp_liq, comp_gas, vaporfra, equalconstant);

	vaporfra = 1.0e-06;
        scalar h_target = 0.0;
	scalar x_fuel0 = 0.0;
        scalar x_fuel = 0.0;
	scalar comp_out[num_sp]= {0.0};
	scalar det_x = -0.01, num_fuel = 100;
	//scalar det_x = -0.02, num_fuel = 50;

	scalar vaporfra0 = vaporfra;
        //temp = templiq;
	for (label nfuel = 0; nfuel < num_fuel; nfuel++)
	{
           temp = templiq;

	   if (num_sp == 2) //C7n2 OR c12n2
	   {
		x_fuel0 = 0.999999;
		x_fuel = x_fuel0 + (nfuel - 0)*det_x;

		scalar h_target_ch = 0.0, h_target_n2 = 0.0;
	   	if (system_ == 1 || system_ == 5 || system_ == 3) //C7n2 OR c12n2
	   	{
		    comp_out[0] = x_fuel;//CH
		    comp_out[1] = 1.0-comp_out[0];
		    h_target_ch = this->Hspd_mix(press, templiq, comp_out, 0); //J/mol: ch enthalpy
		    h_target_n2 = this->Hspd_mix(press, tempgas, comp_out, 1); 
	   	    h_target    = h_target_ch * comp_out[0] + h_target_n2 * (1.0-comp_out[0]);
		}
	   	if (system_ == 2) //CO2CH4
	   	{
		    comp_out[1] = x_fuel;//CH4
		    comp_out[0] = 1.0-comp_out[1];
		    h_target_ch = this->Hspd_mix(press, templiq, comp_out, 1); //J/mol: ch enthalpy
		    h_target_n2 = this->Hspd_mix(press, tempgas, comp_out, 0); 
	   	    h_target    = h_target_ch * comp_out[1] + h_target_n2 * (1.0-comp_out[1]);
		}
//printf("detH11111************ = %5.6f %5.6f %5.6f %5.6f\n",h_target,h_target_ch,h_target_n2,comp_out[0] ); 
	   }
	   else if (num_sp == 4) //CO2ch4n2o2 or CO2CH4O2H2O
	   {
		x_fuel0 = 0.99;
		x_fuel = x_fuel0 + (nfuel - 0)*det_x;
		
		comp_out[1] = x_fuel;//CH4
		comp_out[0] = (1.0-comp_out[1])*1.0/6.0;//CO2=O2
		comp_out[2] = (1.0-comp_out[1])*1.0/6.0;
		comp_out[3] = (1.0-comp_out[1])*4.0/6.0;
		scalar h_target_ch4 = this->Hspd_mix(press, templiq, comp_out, 0); 
		scalar h_target_res = this->Hspd_mix(press, tempgas, comp_out, 1); 
		       h_target    = h_target_ch4 + h_target_res;
	   }
	   else if (num_sp == 5) //CO2ch4n2o2H2O
	   {
		x_fuel0 = 0.9999;
		x_fuel = x_fuel0 + (nfuel - 0)*det_x;
		
		comp_out[1] = x_fuel;//CH4
		//comp_out[0] = (1.0-comp_out[1])*1.0/6.1;//CO2=O2
		//comp_out[2] = (1.0-comp_out[1])*1.0/6.1;
		//comp_out[3] = (1.0-comp_out[1])*4.0/6.1;
		//comp_out[4] = (1.0-comp_out[1])*0.1/6.1;//H2O=O2

		comp_out[0] = (1.0-comp_out[1])*5.0/10.1;//CO2=O2
		comp_out[2] = (1.0-comp_out[1])*1.0/10.1;
		comp_out[3] = (1.0-comp_out[1])*4.0/10.1;
		comp_out[4] = (1.0-comp_out[1])*0.1/10.1;//H2O=O2

		scalar h_target_co2 = this->Hspd_mix(press, templiq, comp_out, 0); 
		scalar h_target_res = this->Hspd_mix(press, tempgas, comp_out, 1); 
			//h_target    = h_target_co2 * comp_out[0] + h_target_res * (1.0-comp_out[0]);
			h_target    = h_target_co2 + h_target_res;
	   }
           scalar mw_mixture = 0.0;
	   scalar cp_m = 0.0, cp_gas = 0.0, cp_liq = 0.0;
	   scalar mw_gas = 0.0;
	   scalar y_gas = 0.0;
	   scalar delta_h  = 1.0e+10;
	   scalar deta_T   = 1000.0;
	   scalar C        = 1000.0;
	   scalar alpha_dT = 0.1;

	   int tpn_flag = 0;
	   label nx = 0, nx_max = 1000;//50000; 
	   for (label spid = 0; spid < num_sp; spid++)
	   {
		comp_liq0[spid] = comp_out[spid];
		comp_gas0[spid] = comp_out[spid];
	   }

	   // ***note the position 
	   vaporfra = vaporfra0;
	   for (label spid = 0; spid < num_sp; spid++)
	   {
	       comp_liq[spid] = comp_liq0[spid];
	       comp_gas[spid] = comp_gas0[spid];
	   }

	   while (nx < nx_max && (delta_h >= 1.0e-02 || delta_h <= -1.0e-02))
	   {
		if (tpn_flag == 0)
		{
	   	    // ***note the position 
		    this->TPn_flash(press, temp, comp_out, comp_liq, comp_gas, vaporfra, equalconstant);
		}

		/*
		if (vaporfra >= 1.0 && nx > 0)
		{
		    tpn_flag = 1;
		}
		if (vaporfra <= -1.0 && nx > 0)
		{
		    tpn_flag = 2;
		}

		if (tpn_flag == 1)
		{
		    for (label spid = 0; spid < num_sp; spid++)
		    {
			comp_liq[spid] = comp_liq0[spid];
			comp_gas[spid] = comp_gas0[spid];
		    }
		    vaporfra = vaporfra0;
		    this->TPn_simple(press, temp, comp_out, comp_liq, comp_gas, vaporfra, equalconstant);
		}

		if (tpn_flag == 2)
		{
		    //break;
		    for (label spid = 0; spid < num_sp; spid++)
		    {
			comp_liq[spid] = comp_liq0[spid];
			comp_gas[spid] = comp_gas0[spid];
		    }
		    vaporfra = vaporfra0;
		    this->TPn_simple(press, temp, comp_out, comp_liq, comp_gas, vaporfra, equalconstant);
		}
		*/

		mw_gas = this->mwmix(comp_gas);
		//mw_liq = this->mwmix(comp_liq);
		//mw_mixture = mw_gas*vaporfra + mw_liq*vaporfra;
                mw_mixture = this->mwmix(comp_out);
		y_gas = vaporfra*mw_gas/mw_mixture;
	        if (y_gas > 1.0)
		{
		   y_gas = 1.0;
		}
		//y_gas = mw_gas/mw_mixture;
		cp_gas = this->cpmix_phase(0, press, temp, comp_gas); 
		cp_liq = this->cpmix_phase(1, press, temp, comp_liq); 
		cp_m = cp_gas*y_gas + cp_liq*(1.0-y_gas); 

		scalar hReal_gas = this->Hmix(0, press, temp, comp_gas);
		scalar hReal_liq = this->Hmix(1, press, temp, comp_liq);
		//scalar hReal_mix = hReal_gas*y_gas + hReal_liq*(1.0-y_gas);
		scalar hReal_mix = hReal_gas*vaporfra + hReal_liq*(1.0-vaporfra);
		delta_h = hReal_mix - h_target;

		deta_T  = - delta_h / cp_m; 
		scalar ddeta_T = deta_T * alpha_dT;
		if (ddeta_T <= 0.0)
		{
		    ddeta_T = - ddeta_T;
		}

		C       = 1.0 / ( 1.0 + ddeta_T );
		temp = temp + deta_T*C;
//printf("detH11111************ = %5.6f %5.6f %5.6f %5.6f %5.8f\n",h_target,hReal_mix,hReal_gas,hReal_liq,vaporfra); 
//printf("detH22222 = %5.3i %5.6f %5.6f %5.6f %5.6f %5.6f %5.6f %5.6f\n", nx,h_target,delta_h, temp, cp_m,vaporfra,y_gas,hReal_liq); 
//printf("detH33333 = %5.3i %5.6f %5.6f %5.6f %5.6f %5.6f %5.6f %5.6f\n", nx, cp_liq, y_gas, comp_liq[0],comp_liq[1],comp_liq[2],comp_liq[3],comp_liq[4]); 
		nx++;
	    }
	    //scalar x_amb = 1.0 - x_fuel;
	    output<<x_fuel<<' '<<' '<<temp<<' '<<' '<<vaporfra<<' '<<' '<<delta_h; 
	    output<< "\n";
 	    
	    vaporfra0 = vaporfra;
	    for (label spid = 0; spid < num_sp; spid++)
	    {
		comp_liq0[spid] = comp_liq[spid];
		comp_gas0[spid] = comp_gas[spid];
		//comp_liq0[spid] = comp_out[spid];
		//comp_gas0[spid] = comp_out[spid];
	    }
printf("HTxTable_*******comp= %5.6f %5.6f %5.6f %5.8f %5.6f %5.8f %5.8f\n\n",x_fuel,temp,cp_m,vaporfra,mw_gas,comp_gas[0],y_gas); 
	}
    }
    output.close();
}




// ======================================================================//
// Get mixture critical point						 //
// Input: feed								 //
// Output: Tcmix, Pcmix, Vcmix, Zcmix					 //
//									 //
// History :   Creation.	     P. YI		2019.06		 //
// ======================================================================//
template<class Specie>
inline void Foam::PengRobinsonMix<Specie>:: initMix_critical
(
    scalar *comp //inp
) const
{
    ofstream output;
    string Mix_critical;    
    output.open("Mix_critical.txt");
    if (num_sp == 2)
    {        
	comp[0] = compinp1;
        comp[1] = compinp2;
    }
    else if (num_sp == 3)
    {
    	comp[0] = compinp1;
	comp[1] = compinp2;
        comp[2] = compinp3;
    }
    else if (num_sp == 4)
    {
    	comp[0] = compinp1;
	comp[1] = compinp2;
        comp[2] = compinp3;
	comp[3] = compinp4;
    }
    else if (num_sp == 5)
    {
    	comp[0] = compinp1;
	comp[1] = compinp2;
        comp[2] = compinp3;
	comp[3] = compinp4;
	comp[4] = compinp5;
    }  
    scalar x_fuel0 = 0.999;
    scalar x_fuel  = 0.0;

    for (label ncomp = 0; ncomp < num_comp; ncomp++)
    {
	x_fuel = x_fuel0 + (ncomp - 0)*det_com1;
	if (num_sp == 2) //C7n2 OR c12n2
	{
	    comp[0] = x_fuel;
	    comp[1] = 1.0-comp[0];
	}
	else if (num_sp == 4) //CO2ch4n2o2
	{
	    comp[1] = x_fuel;//CH4
	    comp[0] = (1.0-comp[1])*1.0/6.0;//CO2=O2
	    comp[2] = (1.0-comp[1])*1.0/6.0;
	    comp[3] = (1.0-comp[1])*4.0/6.0;
	}
	else if (num_sp == 5) //CO2ch4n2o2H2O
	{
	    comp[1] = x_fuel;//CH4
	    //comp[0] = (1.0-comp[1])*1.0/6.1;//CO2=O2
	    //comp[2] = (1.0-comp[1])*1.0/6.1;
	    //comp[3] = (1.0-comp[1])*4.0/6.1;
	    //comp[4] = (1.0-comp[1])*0.1/6.1;//H2O=O2
	    comp[0] = (1.0-comp[1])*5.0/10.1;//CO2=O2
	    comp[2] = (1.0-comp[1])*1.0/10.1;
	    comp[3] = (1.0-comp[1])*4.0/10.1;
	    comp[4] = (1.0-comp[1])*0.1/10.1;//H2O=O2
	}

	scalar Tcmix = 0.0, Pcmix = 0.0, Vcmix = 0.0, Zcmix = 0.0;
	this->Mixcritic_Bellan(Tcmix, Pcmix, Vcmix, Zcmix);

	output<<comp[0]<<' '<<' '<<Pcmix*1.0e-05<<' '<<' '<<Tcmix<<' '<<' '<<Vcmix<<' '<<' '<<Zcmix; 
	output<< "\n";
printf("Screen_Bellan_*******P,T,V,Zcmix= %5.6f %5.6f %5.6f %5.6f %5.6f\n\n",comp[0],Pcmix*1.0e-06,Tcmix,Vcmix,Zcmix); 
    }
    output.close();
}





// ======================================================================//
// Evaluate the phase volume fraction, alpha				 //
// Input: pressure, temperature, vapor fraction				 //
// Output: phase compositions 				 		 //
//									 //
// History :   Creation.	     P. YI			2019     //
// ======================================================================//
template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::Evaluate_alpha
(
    scalar press, 	//inp
    scalar temp,  	//inp
    scalar vaporfra,	//inp
    scalar *comp_liq, 	//inp
    scalar *comp_gas,  	//inp
    scalar *comp_inp
) const
{
    scalar sv_gas = this->volmmix_phase(0, press, temp, comp_gas);  //gas phase specific volume
    scalar sv_liq = this->volmmix_phase(1, press, temp, comp_liq);  //liquid phase specific volume
/*
    scalar mw_mixture = this->mwmix(comp_inp);
    scalar mw_gas = this->mwmix(comp_gas);
    scalar msgas = mw_gas/sv_gas;
    scalar ygas = vaporfra*mw_gas/mw_mixture;
    scalar mw_liq = this->mwmix(comp_liq);
    scalar msliq = mw_liq/sv_liq;
    scalar alphagas = ygas*msliq/(msgas-ygas*(msgas-msliq));
*/
    scalar alphagas = vaporfra*sv_gas/(vaporfra*sv_gas+(1.0-vaporfra)*sv_liq);
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
// Inilization for TPD							 //
// Input: pressure, temperature, vapor fraction				 //
// Output: TPD		 				 		 //
// Althogrim need to be improved, which is not been used, 2019.05   	 //
//									 //
// History :   Creation.	     P. YI			2019     //
// ======================================================================//
template<class Specie>
inline void Foam::PengRobinsonMix<Specie>:: initTPD
(
    scalar press,	//inp
    scalar temp,	//inp
    scalar *comp,	//inp
    scalar *comp_liq,	//outp
    scalar *comp_gas,	//outp
    scalar &tpdtest	//outp
) const
{
    string TPD;
     
    ofstream output;
    output.open("TPD.txt");
    //comp[0] = compinp1;
    //for (label ncomp = 0; ncomp < num_comp; ncomp++)
    //{
	if (num_sp == 2)
	{        
	    comp[0] = compinp1;
            comp[1] = compinp2;
	}
	else if (num_sp == 3)
	{
	    comp[0] = compinp1;
	    comp[1] = compinp2;
            comp[2] = compinp3;
	}
	else if (num_sp == 4)
	{
	    comp[0] = compinp1;
	    comp[1] = compinp2;
            comp[2] = compinp3;
	    comp[3] = compinp4;
	}
	else if (num_sp == 5)
	{
	    comp[0] = compinp1;
	    comp[1] = compinp2;
            comp[2] = compinp3;
	    comp[3] = compinp4;
	    comp[4] = compinp5;
	}
	for (label spid = 0; spid < num_sp; spid++)
	{
	    comp_liq[spid] = comp[spid];
	    comp_gas[spid] = comp[spid];
	} 
	for (label npres = 0; npres < num_pres; npres++)
    	{
	    press = press0 + (npres - 0)*det_pres;
	    for (label ntemp = 0; ntemp < num_temp; ntemp++)
    	    {
	        //printf("Screen_observe_000 %5.8f  %5.8f\n",comp[0],comp[1]); 
	    	temp = temp0 + (ntemp - 0)*det_temp;
		this->TPD(press, temp, comp, comp_liq, comp_gas, tpdtest);
		//output<<comp[0]<<' '<<' '<<' '<<press*1.0e-05<<' '<<' '<<' '<<temp<<' '<<' '<<' '<<tpdtest;
		output<<press*1.0e-05<<' '<<' '<<' '<<temp<<' '<<' '<<' '<<tpdtest;  
		output<< "\n";
printf("Screen_observe_TPD %5.3i %5.3i %5.5f %5.2f %5.2f %5.2f\n\n",npres,ntemp,comp[1],press,temp,tpdtest); 
	    }
	}
    //}
    output.close();
}




// ======================================================================//
// Ideal phase mixture internal energy => U = H - PV			 //
// Input: temperature, phase compositions				 //
// Output: ideal phase mixture sie			 		 //
//									 //
// History :   Creation.	     P. YI			2019     //
// ======================================================================//
template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::sieig_phase 
(
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

    sie_igkg = sie_igkg - RR*1.0e-03*T;  

    return sie_igkg/mw_phase;
}



// ======================================================================//
// Ideal phase mixture entropy						 //
// Input: temperature, phase compositions				 //
// Output: ideal phase mixture S			 		 //
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
inline Foam::scalar Foam::PengRobinsonMix<Specie>::Sig_phase //PY: phase mixture entropy
(
    scalar T, 
    scalar *comp_phase 
) const
{
    //co2_, ch4_, o2_, N2_, h2o_; 
    //coefficients
				//co2		ch4		o2		n2		h2o
    scalar a1s1000[5] = { 2.35677352e+00,  5.14911468e+00,  3.78246636e+09,  3.53100528e+00,  7.25575005e+01 };
    scalar a2s1000[5] = { 8.98459677e-03, -1.36622009e-02, -2.99673416e-03, -1.23660988e-04, -6.62445402e-01 };
    scalar a3s1000[5] = {-7.12356269e-06,  4.91453921e-05,  9.84730200e-06, -5.02999433e-07,  2.56198746e-03 };
    scalar a4s1000[5] = { 2.45919022e-09, -4.84246767e-08, -9.68129608e-09,  2.43530612e-09, -4.36591923e-06 };
    scalar a5s1000[5] = {-1.43699548e-13,  1.66603441e-11,  3.24372836e-12, -1.40881235e-12,  2.78178981e-09 };

    scalar s_ig[num_sp]  = { 0.0 };
    scalar S_ideal = 0.0;
    int num_sp0 = 0, num_spend = num_sp;
    if (system_ == 1.0 || system_ == 5.0)
    {
	num_sp0 = 3, num_spend = 4;
    }
    else if (system_ == 6.0)
    {
	num_sp0 = 1, num_spend = 4;
    }
    else if (system_ == 7.0)
    {
	num_sp0 = 0, num_spend = 3;
    }
    else if (system_ == 8.0)
    {
	num_sp0 = 0, num_spend = 4;
    }
    else if (system_ == 10.0 )
    {
	num_sp0 = 0, num_spend = 4;
	a1s1000[3] = 7.25575005e+01;
	a2s1000[3] = -6.62445402e-01;
	a3s1000[3] = 2.56198746e-03;
	a4s1000[3] = -4.36591923e-06;
	a5s1000[3] = 2.78178981e-09;
    }
    else if (system_ == 9.0)
    {
	num_sp0 = 0, num_spend = 5;
    }

    for (label spid = num_sp0; spid < num_spend; spid++)
    {
        s_ig[spid] = RR*1.0e-03*(a1s1000[spid]*log(T) + a2s1000[spid]*T + a3s1000[spid]*sqr(T)/2.0 + a4s1000[spid]*pow3(T)/3.0 + a5s1000[spid]*pow4(T)/4.0);
	S_ideal = S_ideal + s_ig[spid]*comp_phase[spid];
    }    
				//C12H26	C7H16 //Update the coeffecients for C7
    scalar a1sfue1000[2] = { 2.13264480e+01,  1.85354704e+01 };
    scalar a2sfue1000[2] = {-3.86394002e-02, -6.62445402e-01 };
    scalar a3sfue1000[2] = { 3.99476113e-04,  2.56198746e-03 };
    scalar a4sfue1000[2] = {-5.06681097e-07, -4.36591923e-06 };
    scalar a5sfue1000[2] = { 2.00697878e-10,  2.78178981e-09 };

    scalar sfue_ig[num_sp]  = { 0.0 };
    if (system_ == 1.0)
    {
	for (label spid = 0; spid < 1; spid++)
	{
	    sfue_ig[spid] = RR*1.0e-03*(a1sfue1000[spid]*log(T) + a2sfue1000[spid]*T + a3sfue1000[spid]*sqr(T)/2.0 + a4sfue1000[spid]*pow3(T)/3.0 + a5sfue1000[spid]*pow4(T)/4.0); 
	}    
        s_ig[3] = RR*1.0e-03*(a1s1000[3]*log(T) + a2s1000[3]*T + a3s1000[3]*sqr(T)/2.0 + a4s1000[3]*pow3(T)/3.0 + a5s1000[3]*pow4(T)/4.0);
	S_ideal = sfue_ig[0]*comp_phase[0] + s_ig[3]*comp_phase[1]; //C12 + N2
    }

    if (system_ == 5.0)
    {
	for (label spid = 1; spid < 2; spid++)
	{
	    sfue_ig[spid] = RR*1.0e-03*(a1sfue1000[spid]*log(T) + a2sfue1000[spid]*T + a3sfue1000[spid]*sqr(T)/2.0 + a4sfue1000[spid]*pow3(T)/3.0 + a5sfue1000[spid]*pow4(T)/4.0); 
	}    
        s_ig[3] = RR*1.0e-03*(a1s1000[3]*log(T) + a2s1000[3]*T + a3s1000[3]*sqr(T)/2.0 + a4s1000[3]*pow3(T)/3.0 + a5s1000[3]*pow4(T)/4.0);

	S_ideal = sfue_ig[1]*comp_phase[0] + s_ig[3]*comp_phase[1]; //C12 + N2
    }

    return S_ideal/this->mwmix(comp_phase);
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
inline Foam::scalar Foam::PengRobinsonMix<Specie>::Hig_phase 
(
    scalar T, 
    scalar *comp_phase 
) const
{
    //co2_, ch4_, o2_, N2_, h2o_; 
    //coefficients
				//co2		ch4		o2		n2		h2o
    scalar a1s1000[5] = { 2.35677352e+00,  5.14911468e+00,  3.78246636e+00,  3.53100528e+00,  7.25575005e+01 };
    scalar a2s1000[5] = { 8.98459677e-03, -1.36622009e-02, -2.99673416e-03, -1.23660988e-04, -6.62445402e-01 };
    scalar a3s1000[5] = {-7.12356269e-06,  4.91453921e-05,  9.84730200e-06, -5.02999433e-07,  2.56198746e-03 };
    scalar a4s1000[5] = { 2.45919022e-09, -4.84246767e-08, -9.68129608e-09,  2.43530612e-09, -4.36591923e-06 };
    scalar a5s1000[5] = {-1.43699548e-13,  1.66603441e-11,  3.24372836e-12, -1.40881235e-12,  2.78178981e-09 };

    scalar h_ig[num_sp]  = { 0.0 };
    scalar H_ideal = 0.0;
    int num_sp0 = 0, num_spend = num_sp;
    if (system_ == 1.0 || system_ == 5.0)
    {
	num_sp0 = 3, num_spend = 4;
    }
    else if (system_ == 6.0)
    {
	num_sp0 = 1, num_spend = 4;
    }
    else if (system_ == 7.0)
    {
	num_sp0 = 0, num_spend = 3;
    }
    else if (system_ == 8.0)
    {
	num_sp0 = 0, num_spend = 4;
    }
    else if (system_ == 10.0 )
    {
	num_sp0 = 0, num_spend = 4;
	a1s1000[3] = 7.25575005e+01;
	a2s1000[3] = -6.62445402e-01;
	a3s1000[3] = 2.56198746e-03;
	a4s1000[3] = -4.36591923e-06;
	a5s1000[3] = 2.78178981e-09;
    }
    else if (system_ == 9.0)
    {
	num_sp0 = 0, num_spend = 5;
    }

    for (label spid = num_sp0; spid < num_spend; spid++)
    {
        h_ig[spid] = RR*1.0e-03*T* (a1s1000[spid] + a2s1000[spid]*T/2.0 + a3s1000[spid]*sqr(T)/3.0 + a4s1000[spid]*pow3(T)/4.0 + a5s1000[spid]*pow4(T)/5.0);
	H_ideal = H_ideal + h_ig[spid]*comp_phase[spid];
    }    
				//C12H26	C7H16 
    scalar a1sfue1000[2] = { 2.13264480e+01,  1.11532484e+01 };
    scalar a2sfue1000[2] = {-3.86394002e-02, -9.49415433e-03 };
    scalar a3sfue1000[2] = { 3.99476113e-04,  1.95571181e-04 };
    scalar a4sfue1000[2] = {-5.06681097e-07, -2.49752520e-07 };
    scalar a5sfue1000[2] = { 2.00697878e-10,  9.84873213e-11 };

    scalar hfuel_ig[num_sp]  = { 0.0 };
    if (system_ == 1.0 )
    {
	for (label spid = 0; spid < 1; spid++)
	{
	    hfuel_ig[spid] = RR*1.0e-03*T*(a1sfue1000[spid] + a2sfue1000[spid]*T/2.0 + a3sfue1000[spid]*sqr(T)/3.0 + a4sfue1000[spid]*pow3(T)/4.0 + a5sfue1000[spid]*pow4(T)/5.0);
	}   
        h_ig[3] = RR*1.0e-03*T* (a1s1000[3] + a2s1000[3]*T/2.0 + a3s1000[3]*sqr(T)/3.0 + a4s1000[3]*pow3(T)/4.0 + a5s1000[3]*pow4(T)/5.0);
 
	H_ideal = hfuel_ig[0]*comp_phase[0] + h_ig[3]*comp_phase[1]; //C12 + N2
    }
    if (system_ == 5.0)
    {
	for (label spid = 1; spid < 2; spid++)
	{
	    hfuel_ig[spid] = RR*1.0e-03*T*(a1sfue1000[spid] + a2sfue1000[spid]*T/2.0 + a3sfue1000[spid]*sqr(T)/3.0 + a4sfue1000[spid]*pow3(T)/4.0 + a5sfue1000[spid]*pow4(T)/5.0);
	}   
        h_ig[3] = RR*1.0e-03*T* (a1s1000[3] + a2s1000[3]*T/2.0 + a3s1000[3]*sqr(T)/3.0 + a4s1000[3]*pow3(T)/4.0 + a5s1000[3]*pow4(T)/5.0);
 
	H_ideal = hfuel_ig[1]*comp_phase[0] + h_ig[3]*comp_phase[1]; //C12 + N2
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
template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::Hig2_phase 
(
    scalar T, 
    scalar *comp_phase 
) const
{
    //co2_, ch4_, o2_, N2_, h2o_; 
    //coefficients
		    //co2	ch4	  o2	   n2	   h2o
    scalar BB[5] = { 7.54056, 8.003180, 6.96302, 6.95808, 7.97183 };
    scalar CC[5] = { 7.51625, 19.26330, 2.40013, 2.03952, 6.27078 };
    scalar DD[5] = { 1442.70, 2148.790, 2522.05, 1681.60, 2572.63 };
    scalar EE[5] = { 5.38023, 10.44230, 2.21752, 0.506863, 2.0501 };
    scalar FF[5] = { 647.502, 1017.070, 1154.15, 6535.68, 1156.72 };

    scalar h_ig[num_sp]  = { 0.0 };
    scalar H_ideal = 0.0;
    int num_sp0 = 0, num_spend = num_sp;
    if (system_ == 1.0 || system_ == 5.0)
    {
	num_sp0 = 3, num_spend = 4;
    }
    else if (system_ == 6.0)
    {
	num_sp0 = 1, num_spend = 4;
    }
    else if (system_ == 7.0)
    {
	num_sp0 = 0, num_spend = 3;
    }
    else if (system_ == 8.0)
    {
	num_sp0 = 0, num_spend = 4;
    }
    else if (system_ == 10.0 )
    {
	num_sp0 = 0, num_spend = 4;
	BB[3] = 7.97183;
	CC[3] = 6.27078;
	DD[3] = 2572.63;
	EE[3] = 2.0501;
	FF[3] = 1156.72;
    }
    else if (system_ == 9.0)
    {
	num_sp0 = 0, num_spend = 5;
    }

    for (label spid = num_sp0; spid < num_spend; spid++)
    {
        h_ig[spid] = BB[spid]*T +CC[spid]*DD[spid]*(2.0/(exp(2.0*DD[spid]/T)-1.0))+ 2.0*EE[spid]*FF[spid]*(1.0/(exp(2.0*FF[spid]/T)+1.0));
	H_ideal = H_ideal + h_ig[spid]*comp_phase[spid];
    } 
    H_ideal = H_ideal * 4.18680;
			   //C12H26	C7H16   C20H42
    scalar BBf[3] = { 51.02700, 30.4029, 82.7586 };
    scalar CCf[3] = { 143.8250, 90.6941, 226.791 };
    scalar DDf[3] = { 1664.760, 1669.32, 1587.50 };
    scalar EEf[3] = { 109.4510, 63.2028, 175.085 };
    scalar FFf[3] = { 789.2970, 764.611, 769.249 };

    scalar hfuel_ig[num_sp]  = { 0.0 };
    if (system_ == 1.0)
    {
	for (label spid = 0; spid < 1; spid++)
	{
	    hfuel_ig[spid] = BBf[spid]*T +CCf[spid]*DDf[spid]*(2.0/(exp(2.0*DDf[spid]/T)-1.0))+ 2.0*EEf[spid]*FFf[spid]*(1.0/(exp(2.0*FFf[spid]/T)+1.0));
	}    
	h_ig[3] = BB[3]*T +CC[3]*DD[3]*(2.0/(exp(2.0*DD[3]/T)-1.0))+ 2.0*EE[3]*FF[3]*(1.0/(exp(2.0*FF[3]/T)+1.0));
	H_ideal = (hfuel_ig[0]*comp_phase[0] + h_ig[3]*comp_phase[1])*4.18680;
    }

    if (system_ == 5.0)
    {
	for (label spid = 1; spid < 2; spid++)
	{
	    hfuel_ig[spid] = BBf[spid]*T +CCf[spid]*DDf[spid]*(2.0/(exp(2.0*DDf[spid]/T)-1.0))+ 2.0*EEf[spid]*FFf[spid]*(1.0/(exp(2.0*FFf[spid]/T)+1.0));
	}    
	h_ig[3] = BB[3]*T +CC[3]*DD[3]*(2.0/(exp(2.0*DD[3]/T)-1.0))+ 2.0*EE[3]*FF[3]*(1.0/(exp(2.0*FF[3]/T)+1.0));
	H_ideal = (hfuel_ig[1]*comp_phase[0] + h_ig[3]*comp_phase[1])*4.18680;
    }

    return H_ideal;
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
inline void Foam::PengRobinsonMix<Specie>::cvcpig_phase 
(
    scalar T, 
    scalar *comp_phase,
    scalar& cpig,  
    scalar& cvig
) const
{
    //co2_, ch4_, o2_, N2_, h2o_; 
    //coefficients
				//co2		ch4		o2		n2		h2o
    scalar a1s1000[5] = { 2.35677352e+00,  5.14911468e+00,  3.78246636e+00,  3.53100528e+00,  7.25575005e+01 };
    scalar a2s1000[5] = { 8.98459677e-03, -1.36622009e-02, -2.99673416e-03, -1.23660988e-04, -6.62445402e-01 };
    scalar a3s1000[5] = {-7.12356269e-06,  4.91453921e-05,  9.84730200e-06, -5.02999433e-07,  2.56198746e-03 };
    scalar a4s1000[5] = { 2.45919022e-09, -4.84246767e-08, -9.68129608e-09,  2.43530612e-09, -4.36591923e-06 };
    scalar a5s1000[5] = {-1.43699548e-13,  1.66603441e-11,  3.24372836e-12, -1.40881235e-12,  2.78178981e-09 };

    scalar cp_ig[num_sp]  = { 0.0 };
    scalar cp_ideal = 0.0;
    int num_sp0 = 0, num_spend = num_sp;

    if (system_ == 1.0 || system_ == 5.0)
    {
	num_sp0 = 3, num_spend = 4;
    }
    else if (system_ == 6.0)
    {
	num_sp0 = 1, num_spend = 4;
    }
    else if (system_ == 7.0)
    {
	num_sp0 = 0, num_spend = 3;
    }
    else if (system_ == 8.0)
    {
	num_sp0 = 0, num_spend = 4;
    }
    else if (system_ == 10.0 )
    {
	num_sp0 = 0, num_spend = 4;
	a1s1000[3] = 7.25575005e+01;
	a2s1000[3] = -6.62445402e-01;
	a3s1000[3] = 2.56198746e-03;
	a4s1000[3] = -4.36591923e-06;
	a5s1000[3] = 2.78178981e-09;
    }
    else if (system_ == 9.0)
    {
	num_sp0 = 0, num_spend = 5;
    }

    for (label spid = num_sp0; spid < num_spend; spid++)
    {
        cp_ig[spid] = RR*1.0e-03*(a1s1000[spid] + a2s1000[spid]*T + a3s1000[spid]*sqr(T) + a4s1000[spid]*pow3(T) + a5s1000[spid]*pow4(T));
	cp_ideal = cp_ideal + cp_ig[spid]*comp_phase[spid];
    }    
				//C12H26	C7H16 //Update the coeffecients for C7
    scalar a1sfue1000[2] = { 2.13264480e+01,  1.85354704e+01 };
    scalar a2sfue1000[2] = {-3.86394002e-02, -6.62445402e-01 };
    scalar a3sfue1000[2] = { 3.99476113e-04,  2.56198746e-03 };
    scalar a4sfue1000[2] = {-5.06681097e-07, -4.36591923e-06 };
    scalar a5sfue1000[2] = { 2.00697878e-10,  2.78178981e-09 };

    scalar cpfuel_ig[num_sp]  = { 0.0 };
    if (system_ == 1.0)
    {
	for (label spid = 0; spid < 1; spid++)
	{
	    cpfuel_ig[spid] = RR*1.0e-03*(a1sfue1000[spid] + a2sfue1000[spid]*T + a3sfue1000[spid]*sqr(T) + a4sfue1000[spid]*pow3(T) + a5sfue1000[spid]*pow4(T));
	    cp_ig[3] = RR*1.0e-03*(a1s1000[3] + a2s1000[3]*T + a3s1000[3]*sqr(T) + a4s1000[3]*pow3(T) + a5s1000[3]*pow4(T));
	}    
	cp_ideal = cpfuel_ig[0]*comp_phase[0] + cp_ig[3]*comp_phase[1]; //C12 + N2
    }
    if (system_ == 5.0)
    {
	for (label spid = 1; spid < 2; spid++)
	{
	    cpfuel_ig[spid] = RR*1.0e-03*(a1sfue1000[spid] + a2sfue1000[spid]*T + a3sfue1000[spid]*sqr(T) + a4sfue1000[spid]*pow3(T) + a5sfue1000[spid]*pow4(T));
	    cp_ig[3] = RR*1.0e-03*(a1s1000[3] + a2s1000[3]*T + a3s1000[3]*sqr(T) + a4s1000[3]*pow3(T) + a5s1000[3]*pow4(T));
	}    
	cp_ideal = cpfuel_ig[1]*comp_phase[0] + cp_ig[3]*comp_phase[1]; //C12 + N2
    }

    scalar cv_ideal = cp_ideal - RR*1.0-03;
    //printf("cvcpig_phase %6.8f %6.8f %6.8f %6.8f %6.8f\n",cp_ideal,cpfuel_ig[0],cp_ig[3],comp_phase[0],comp_phase[1]); 
    cpig = cp_ideal/this->mwmix(comp_phase);
    cvig = cv_ideal/this->mwmix(comp_phase);
}






// PY: Refer to liquidMixtureProperties.C and icoPolynomial.C to locate critical properties
template<class Specie>
inline void Foam::PengRobinsonMix<Specie>::calccomp //PY:calculate the compositions
(
    scalar *comp //output, YP: how to add the species name
) const
{
    scalar mole_tol = 0.0;
    scalar mol_num[num_sp] = {0.0}; 
    comp[num_sp] = {0.0};
    
    for (label spid = 0; spid < num_sp; spid++)
    {
        mol_num[spid] = this->initialization();
	mole_tol += mol_num[spid]; 
    }    
    
    for (label spid = 0; spid < num_sp; spid++)
    {
        comp[spid] = mol_num[spid]/mole_tol;
    } 
}



template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::mwmix // kg/mol
(
    scalar *comp_inp
) const
{

    scalar mwmix = 0.0;
    for (label spid = 0; spid < num_sp; spid++)
    {
        mwmix += comp_inp[spid]*mw_sp[spid];  //PY: Y() is mole number of component in mixture
    }    
    return mwmix*1.0e-03;
}



template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::yTox_phase
(
    scalar yGasInp,	//Yi, input
    scalar *comp_gas,	//yi, input
    scalar *comp_liq	//xi, input
) const
{
    scalar mwGasInp = this->mwmix(comp_gas);
    scalar mwLiqInp = this->mwmix(comp_liq);
    scalar gasratio = yGasInp/mwGasInp;
    scalar liqratio = (1.0-yGasInp)/mwLiqInp; 
    scalar xGas = gasratio/(gasratio+liqratio); 

    return xGas;
}



template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::S_phase(int phasei,scalar p, scalar T, scalar *comp) const //species entropy
{
    scalar Pr_sp[num_sp] = {0.0};
    scalar Tr_sp[num_sp] = {0.0};

    for (label spid = 0; spid < num_sp; spid++)
    {
	Pr_sp[spid] = p/Pc_sp[spid];
	Tr_sp[spid] = T/Tc_sp[spid];
    }   
    
    scalar amix = 0.0, bmix = 0.0, dadt = 0.0, dadt2 = 0.0;
    this->abmix(T, comp, amix, bmix, dadt, dadt2); 

    scalar mw_phase = this->mwmix(comp);//kg*mol-1
    scalar rho_phase = this->rhomix_phase(phasei, p, T, comp);//kg*m-3

    const scalar cc1 = 1.0 + sqrt(2.0);
    const scalar cc2 = 1.0 - sqrt(2.0);
    const scalar cc0 = 2.0*sqrt(2.0);

    scalar VV = mw_phase/rho_phase; // m3*mol-1,phase specific volume
    scalar xxx = (VV+cc1*bmix)/(VV+cc2*bmix);
    scalar dS_phase = 0.0;
    if (xxx > 0.0)
    {
	dS_phase = RR*log((p*VV-bmix*p)/(RR*1.0e-03*T))+dadt/(cc0*bmix)*log((VV+cc1*bmix)/(VV+cc2*bmix));
    }
    else
    {
	dS_phase = 0.0;
    }

    dS_phase = dS_phase/mw_phase;

    // Remark here, the ideal gas properties can be known, but how about the ideal liquid properites?????
    scalar Sig = this->sieig_phase(T, comp);
    S_phase = dS_phase + Sig;

    return S_phase;
}


template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::H_phase(int phasei, scalar p, scalar T, scalar *comp) const //species enthalpy
{
    const scalar Pr = p/Pc_;
    const scalar Tr = T/Tc_;
    const scalar B = 0.07780*Pr/Tr;
    const scalar kappa = 0.37464 + 1.54226*omega_ - 0.26992*sqr(omega_);
    const scalar alpha = sqr(1 + kappa*(1 - sqrt(Tr)));
  
    scalar Amixz = 0.0, Bmixz = 0.0, Amix2 = 0.0;
    scalar Bz[num_sp] = {0.0};
    this->ABZmix(p, T, comp, Amixz, Bmixz, Amix2, Bz);
    const scalar Zmix = this->Zmix(phasei,Amixz, Bmixz);

    return
        this->R()
       *Tc_
       *(
           Tr*(Zmix - 1.0)
         - 2.078*(1 + kappa)*sqrt(alpha)
          *log((Zmix + 2.414*B)/(Zmix - 0.414*B))
        );
}


template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::psi_phase(int phasei, scalar p, scalar T, scalar *comp) const
{
      
    scalar Amixz = 0.0, Bmixz = 0.0, Amix2 = 0.0;
    scalar Bz[num_sp] = {0.0};
    this->ABZmix(p, T, comp, Amixz, Bmixz, Amix2, Bz);
    const scalar Zmix = this->Zmix(phasei,Amixz, Bmixz);

    return 1.0/(Zmix*this->R()*T); //there is a bug
}


template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::CpMCv
(
    scalar p,
    scalar T
) const
{
    const scalar Tr = T/Tc_;
    const scalar a = 0.45724*sqr(RR*Tc_)/Pc_;
    const scalar b = 0.07780*RR*Tc_/Pc_;
    const scalar kappa = 0.37464 + 1.54226*omega_ - 0.26992*sqr(omega_);
    const scalar alpha = sqr(1 + kappa*(1 - sqrt(Tr)));

    const scalar A = alpha*a*p/sqr(RR*T);
    const scalar B = b*p/(RR*T);

    const scalar Z = this->Z(p, T);

    const scalar ap = kappa*a*(kappa/Tc_ - (1 + kappa)/sqrt(T*Tc_));
    const scalar M = (sqr(Z) + 2*B*Z - sqr(B))/(Z - B);
    const scalar N = ap*B/(b*RR);

    return this->R()*sqr(M - N)/(sqr(M) - 2*A*(Z + B));
}



template<class Specie>
inline void Foam::PengRobinsonMix<Specie>::ABZmix 
(
    scalar p,		//inp
    scalar T,		//inp
    scalar *comp_inp,	//inp
    scalar &Amixz,	//outp
    scalar &Bmixz,	//outp
    scalar *Amix2,	//outp
    scalar *Bz		//outp
) const
{
    scalar bipc12n2 = bipc12n2_;   
    scalar bipco2ch4 = bipco2ch4_;   
    scalar bipco2h2o = bipco2h2o_;   
    scalar bipch4o2 = bipch4o2_;  
    scalar bipc7n2 = bipc7n2_;  

    scalar Pr_sp[num_sp] = {0.0};
    scalar Tr_sp[num_sp] = {0.0};
    scalar bip_kij[num_sp][num_sp] = { {0.0} };
    scalar m_omega_sp[num_sp] = {0.0};
    scalar Az[num_sp] = {0.0};
    scalar Aijz[num_sp][num_sp] = { {0.0} };    
    scalar alpha[num_sp] = {0.0};

    for (label spid = 0; spid < num_sp; spid++)
    {
	Pr_sp[spid] = p/Pc_sp[spid];
	Tr_sp[spid] = T/Tc_sp[spid];
	//Info<< "ABZmix *****" <<spid <<","<<Tr_sp[spid] << ", " << T << endl;
    } 

    if (system_ < 6.0)
    {
	for (label spid = 0; spid < num_sp; spid++)
	{
	    for (label spjd = 0; spjd < num_sp; spjd++)
            {	  
		if (system_ == 1)
		{
		    bip_kij[spid][spjd] = bipc12n2;
		}
		if (system_ == 2)
		{
		    bip_kij[spid][spjd] = bipco2ch4;
		}
		if (system_ == 3)
		{
		    bip_kij[spid][spjd] = bipco2h2o;
		    //bip_kij[0][1] = -0.0814;
		    //bip_kij[1][0] = 0.2186;
		}
		if (system_ == 4)
		{
		    bip_kij[spid][spjd] = bipch4o2;
		}
		if (system_ == 5)
		{
		    bip_kij[spid][spjd] = bipc7n2;
		}

	        if (spid == spjd) 
	        {
		    bip_kij[spid][spjd] = 0.0;
		}
	    }	  	
	}    
    }  
    if (system_ == 6.0)
    {
	string bipkij3_CH4O2N2;
	ifstream input; //infile;
	input.open("bipkij3_CH4O2N2.txt");
	for (label spid = 0; spid < num_sp; spid++)
	{
	    for (label spjd = 0; spjd < num_sp; spjd++)
            {	  
		input >> bip_kij[spid][spjd];
	    }
	}
	input.close();
    }
    if (system_ == 7.0)
    {
	//string bipkij3;
	//ifstream input; //infile;
	//input.open("bipkij3.txt");//co2+ch4+o2
	string bipkij3_CO2CH4H2O;
	ifstream input; //infile;	
	input.open("bipkij3_CO2CH4H2O.txt");//co2+ch4+H2O
	for (label spid = 0; spid < num_sp; spid++)
	{
	    for (label spjd = 0; spjd < num_sp; spjd++)
            {	  
		input >> bip_kij[spid][spjd];
	    }
	}
	input.close();
    }  
    if (system_ == 8.0)
    {
	string bipkij4;
	ifstream input; //infile;
	input.open("bipkij4.txt");
	for (label spid = 0; spid < num_sp; spid++)
	{
	    for (label spjd = 0; spjd < num_sp; spjd++)
            {	  
		input >> bip_kij[spid][spjd];
	    }
	}
	input.close();
    }  

    if (system_ == 9.0)
    {
	string bipkij5;
	ifstream input; //infile;
	input.open("bipkij5.txt");
	for (label spid = 0; spid < num_sp; spid++)
	{
	    for (label spjd = 0; spjd < num_sp; spjd++)
            {	  
		input >> bip_kij[spid][spjd];
	    }
	}
	input.close();
    } 
    if (system_ == 10.0)
    {
	string bipkij4_h2o;
	ifstream input; //infile;
	input.open("bipkij4_h2o.txt");
	for (label spid = 0; spid < num_sp; spid++)
	{
	    for (label spjd = 0; spjd < num_sp; spjd++)
            {	  
		input >> bip_kij[spid][spjd];
	    }
	}
	input.close();
    }  
    if (system_ == 11.0)
    {
	string bipkij_ch;
	ifstream input; 
	input.open("bipkij_ch.txt"); //Ref: Li, H., Thermodynamic properties of CO2 mixtures and their applications in advanced power cycles with CO2 capture processes. 2008
	for (label spid = 0; spid < num_sp; spid++)
	{
	    for (label spjd = 0; spjd < num_sp; spjd++)
            {	  
		input >> bip_kij[spid][spjd];
	    }
	}
	input.close();
    } 

    if (system_ == 12.0)//CH4+n2
    {
	for (label spid = 0; spid < num_sp; spid++)
	{
	    for (label spjd = 0; spjd < num_sp; spjd++)
            {	  
		bip_kij[spid][spjd] = 0.0305;
	        if (spid == spjd) 
	        {
		    bip_kij[spid][spjd] = 0.0;
		}
	    }
	}
    }
    if (system_ == 13.0)//H2+n2
    {
	for (label spid = 0; spid < num_sp; spid++)
	{
	    for (label spjd = 0; spjd < num_sp; spjd++)
            {	  
		bip_kij[spid][spjd] = 0.103;
	        if (spid == spjd) 
	        {
		    bip_kij[spid][spjd] = 0.0;
		}
	    }
	}
    }
    if (system_ == 14.0)//CH4-H2S
    {
	for (label spid = 0; spid < num_sp; spid++)
	{
	    for (label spjd = 0; spjd < num_sp; spjd++)
            {	  
		bip_kij[spid][spjd] = 0.08;
	        if (spid == spjd) 
	        {
		    bip_kij[spid][spjd] = 0.0;
		}
	    }
	}
    }
    	//Info<< "bipkij3_CO2CH4H2O 222*****" << endl;  
    for (label spid = 0; spid < num_sp; spid++)
    {
         if (omega_sp[spid] > 0.49) 
	 {
             m_omega_sp[spid] = 0.379642 + 1.48503*omega_sp[spid] - 0.164423*sqr(omega_sp[spid]) + 0.016666*pow3(omega_sp[spid]);
	 }
         else
	 {
             m_omega_sp[spid] = 0.37464 + 1.54226*omega_sp[spid] - 0.26992*sqr(omega_sp[spid]);
         }
         alpha[spid] = sqr(1.0 + m_omega_sp[spid]*(1.0 - sqrt(Tr_sp[spid])));
         Az[spid] = 0.45724*alpha[spid]*Pr_sp[spid]/sqr(Tr_sp[spid]);
         Bz[spid] = 0.07780*Pr_sp[spid]/Tr_sp[spid];
    }
    	//Info<< "bipkij3_CO2CH4H2O 333*****" << endl;  
    for (label spid = 0; spid < num_sp; spid++)
    {
        for (label spjd = 0; spjd < num_sp; spjd++)
        {
            Aijz[spid][spjd] = sqrt(Az[spid]*Az[spjd])*(1.0 - bip_kij[spid][spjd]);
        }
    }   
     //Info<< "Aijz *****" <<Az[0] <<",,"<<Az[1] << ", " << Aijz[0][0] <<", " << Aijz[0][1] <<", " << Aijz[1][0] <<", " << Aijz[1][1] << endl;        
    // return these values
    Amixz = 0.0; 
    Bmixz = 0.0; 

    for (label spid = 0; spid < num_sp; spid++)
    {
        Bmixz += comp_inp[spid] * Bz[spid]; 
	Amix2[spid] = 0.0;       
        for (label spjd = 0; spjd < num_sp; spjd++)
        {            
            Amixz += comp_inp[spid]*comp_inp[spjd]*Aijz[spid][spjd];
            Amix2[spid] += comp_inp[spjd]*Aijz[spjd][spid];
        }   
    }
//printf("Amixz_Bmixz= %+5.8f,%+5.8f, %+5.8f,%+5.8f\n", Amixz,Bmixz,Amix2[0],Amix2[1]); 
//FatalErrorInFunction
//<< exit(FatalError);    
}



// PY: need to use dynamc array to replace static array
template<class Specie>
inline void Foam::PengRobinsonMix<Specie>::abmix 
(
    scalar T,
    scalar *comp_inp,
    scalar &amix,
    scalar &bmix,
    scalar &dadt,
    scalar &dadt2
) const
{
    scalar bipc12n2 = bipc12n2_;   
    scalar bipco2ch4 = bipco2ch4_;   
    scalar bipco2h2o = bipco2h2o_;   
    scalar bipch4o2 = bipch4o2_;  
    scalar bipc7n2 = bipc7n2_;  

    scalar bip_kij[num_sp][num_sp] = { {0.0} };
    scalar m_omega_sp[num_sp] = {0.0};
    scalar az[num_sp] = {0.0};
    scalar bz[num_sp] = {0.0};
    scalar aijz[num_sp][num_sp] = { {0.0} };    
    scalar alpha[num_sp] = {0.0};
    scalar Tr_sp[num_sp] = {0.0};

    for (label spid = 0; spid < num_sp; spid++)
    {
        Tr_sp[spid] = T/Tc_sp[spid];
    }    
    
    if (system_ < 6.0)
    {
	for (label spid = 0; spid < num_sp; spid++)
	{
	    for (label spjd = 0; spjd < num_sp; spjd++)
            {	  
		if (system_ == 1.0)
		{
		    bip_kij[spid][spjd] = bipc12n2;
		}
		if (system_ == 2.0)
		{
		    bip_kij[spid][spjd] = bipco2ch4;
		}
		if (system_ == 3.0)
		{
		    bip_kij[spid][spjd] = bipco2h2o;
		    //bip_kij[0][1] = -0.0814;
		    //bip_kij[1][0] = 0.2186;
		}
		if (system_ == 4.0)
		{
		    bip_kij[spid][spjd] = bipch4o2;
		}
		if (system_ == 5.0)
		{
		    bip_kij[spid][spjd] = bipc7n2;
		}

	        if (spid == spjd) 
	        {
		    bip_kij[spid][spjd] = 0.0;
		}
	    }	  	
	}    
    }  
    if (system_ == 6.0)
    {
	string bipkij3_CH4O2N2;
	ifstream input; //infile;
	input.open("bipkij3_CH4O2N2.txt");
	for (label spid = 0; spid < num_sp; spid++)
	{
	    for (label spjd = 0; spjd < num_sp; spjd++)
            {	  
		input >> bip_kij[spid][spjd];
	    }
	}
	input.close();
    }  
    if (system_ == 7.0)
    {
	//string bipkij3;
	//ifstream input; //infile;
	//input.open("bipkij3.txt");//co2+ch4+o2
	string bipkij3_CO2CH4H2O;
	ifstream input; //infile;	
	input.open("bipkij3_CO2CH4H2O.txt");//co2+ch4+H2O

	for (label spid = 0; spid < num_sp; spid++)
	{
	    for (label spjd = 0; spjd < num_sp; spjd++)
            {	  
		input >> bip_kij[spid][spjd];
	    }
	}
	input.close();
    }  
    if (system_ == 8.0)
    {
	string bipkij4;
	ifstream input; //infile;
	input.open("bipkij4.txt");
	for (label spid = 0; spid < num_sp; spid++)
	{
	    for (label spjd = 0; spjd < num_sp; spjd++)
            {	  
		input >> bip_kij[spid][spjd];
	    }
	}
	input.close();
    }    
    if (system_ == 9.0)
    {
	string bipkij5;
	ifstream input; //infile;
	input.open("bipkij5.txt");
	for (label spid = 0; spid < num_sp; spid++)
	{
	    for (label spjd = 0; spjd < num_sp; spjd++)
            {	  
		input >> bip_kij[spid][spjd];
	    }
	}
	input.close();
    }  
    if (system_ == 10.0)
    {
	string bipkij4_h2o;
	ifstream input; //infile;
	input.open("bipkij4_h2o.txt");
	for (label spid = 0; spid < num_sp; spid++)
	{
	    for (label spjd = 0; spjd < num_sp; spjd++)
            {	  
		input >> bip_kij[spid][spjd];
	    }
	}
	input.close();
    }  

    if (system_ == 11.0)
    {
	string bipkij_ch;
	ifstream input; //infile;
	input.open("bipkij_ch.txt");
	for (label spid = 0; spid < num_sp; spid++)
	{
	    for (label spjd = 0; spjd < num_sp; spjd++)
            {	  
		input >> bip_kij[spid][spjd];
	    }
	}
	input.close();
    }   
    if (system_ == 12.0)//CH4+n2
    {
	for (label spid = 0; spid < num_sp; spid++)
	{
	    for (label spjd = 0; spjd < num_sp; spjd++)
            {	  
		bip_kij[spid][spjd] = 0.0305;
	        if (spid == spjd) 
	        {
		    bip_kij[spid][spjd] = 0.0;
		}
	    }
	}
    }
    if (system_ == 13.0)//H2+n2
    {
	for (label spid = 0; spid < num_sp; spid++)
	{
	    for (label spjd = 0; spjd < num_sp; spjd++)
            {	  
		bip_kij[spid][spjd] = 0.103;
	        if (spid == spjd) 
	        {
		    bip_kij[spid][spjd] = 0.0;
		}
	    }
	}
    }
    if (system_ == 14.0)//CH4+H2S
    {
	for (label spid = 0; spid < num_sp; spid++)
	{
	    for (label spjd = 0; spjd < num_sp; spjd++)
            {	  
		bip_kij[spid][spjd] = 0.08;
	        if (spid == spjd) 
	        {
		    bip_kij[spid][spjd] = 0.0;
		}
	    }
	}
    }
       
    scalar comp4[num_sp] = {0.0};
    for (label spid = 0; spid < num_sp; spid++)
    {
        if (omega_sp[spid] > 0.49) 
	{
            m_omega_sp[spid] = 0.379642 + 1.48503*omega_sp[spid] - 0.164423*sqr(omega_sp[spid]) + 0.016666*pow3(omega_sp[spid]);
	}
        else
	{
            m_omega_sp[spid] = 0.37464 + 1.54226*omega_sp[spid] - 0.26992*sqr(omega_sp[spid]);
        }
        alpha[spid] = sqr(1.0 + m_omega_sp[spid] * (1.0 - sqrt(Tr_sp[spid])));
        az[spid] = 0.45724 * alpha[spid] * sqr(RR*1.0e-03*Tc_sp[spid])/Pc_sp[spid];
        bz[spid] = 0.07780 * RR * 1.0e-03* Tc_sp[spid]/Pc_sp[spid];
        comp4[spid] = m_omega_sp[spid]/sqrt(Tc_sp[spid]);
    }
    
    scalar comp5[num_sp][num_sp] = { {0.0} };
    scalar comp6[num_sp][num_sp] = { {0.0} };
    scalar comp7[num_sp][num_sp] = { {0.0} };
    for (label spid = 0; spid < num_sp; spid++)
    {
        for (label spjd = 0; spjd < num_sp; spjd++)
        {
            aijz[spid][spjd] = sqrt(az[spid]*az[spjd]) * (1.0 - bip_kij[spid][spjd]);
            comp5[spid][spjd] = comp4[spid] * sqrt(alpha[spjd]);
            comp6[spid][spjd] = sqrt(alpha[spid]) * comp4[spjd];
            comp7[spid][spjd] = comp5[spid][spjd] + comp6[spid][spjd];
        }
    }   

    // return these values
    amix = 0.0; 
    bmix = 0.0; 
    scalar comp1[num_sp] = {0.0};
    scalar comp2[num_sp][num_sp] = { {0.0} };
    for (label spid = 0; spid < num_sp; spid++)
    {
        for (label spjd = 0; spjd < num_sp; spjd++)
        {            
            comp1[spid] = comp_inp[spid] * aijz[spjd][spid];
            comp2[spid][spjd] = 0.457240 * sqr(RR*1.0e-03) * (bip_kij[spid][spjd] - 1.0);
        }    
        amix += comp1[spid] * comp_inp[spid];//a*alpha
        bmix += comp_inp[spid] * bz[spid];//b
    }   

    scalar comp3[num_sp][num_sp] = { {0.0} };
    scalar dadt1[num_sp][num_sp] = { {0.0} };
    scalar dadt12[num_sp][num_sp] = { {0.0} };
    scalar dadt22[num_sp][num_sp] = { {0.0} };
    for (label spid = 0; spid < num_sp; spid++)
    {
        for (label spjd = 0; spjd < num_sp; spjd++)
        {            
            comp2[spid][spjd]  = comp2[spid][spjd] * Tc_sp[spid] * Tc_sp[spjd];
            comp3[spid][spjd]  = comp2[spid][spjd] * 1.0/sqrt(Pc_sp[spid]*Pc_sp[spjd]);
            dadt1[spid][spjd]  = comp3[spid][spjd] * (1.0/(2.0*sqrt(T)));
	    dadt12[spid][spjd] = dadt1[spid][spjd] * comp7[spid][spjd];
	    dadt22[spid][spjd] = -0.457240*sqr(RR*1.0e-03) * (1.0-bip_kij[spid][spjd])/
               (2.0*sqrt(T)) * (m_omega_sp[spjd]*sqrt(alpha[spid])*Tc_sp[spid]*sqrt(Tc_sp[spjd])
               / sqrt(Pc_sp[spid]*Pc_sp[spjd]) * (-m_omega_sp[spjd]/(2.0*sqrt(Tc_sp[spjd]*T*alpha[spid]))
               - 1.0/(2.0*T)) + m_omega_sp[spid] * sqrt(alpha[spjd])*Tc_sp[spjd]*sqrt(Tc_sp[spid])
               / sqrt(Pc_sp[spid]*Pc_sp[spjd])*(-m_omega_sp[spjd]/(2.0*sqrt(Tc_sp[spjd]*T*alpha[spjd]))
               -1.0/(2.0*T)));
        }    
    }   

    dadt = 0.0;
    dadt2 = 0.0;
    scalar comp8[num_sp] = {0.0};
    scalar comp88[num_sp] = {0.0};
    for (label spid = 0; spid < num_sp; spid++)
    {
        for (label spjd = 0; spjd < num_sp; spjd++)
        {            
            comp8[spid] = comp_inp[spid]*dadt12[spjd][spid];
            comp88[spid] = comp_inp[spid]*dadt22[spjd][spid];
        }    
        dadt  = dadt + comp8[spid]*comp_inp[spid];
        dadt2 = dadt2 + comp88[spid]*comp_inp[spid];
    }   
}



//can be used for phases or mixtures, depending on the Azin and Bzin, three roots
template<class Specie>
inline void Foam::PengRobinsonMix<Specie>::Root_Z
(
    int phasei, //input
    scalar Azin,
    scalar Bzin,
    scalar *rootz
) const
{
    const scalar A = Azin;
    const scalar B = Bzin;

    const scalar a2 = B - 1.0;
    const scalar a1 = A - 2.0*B - 3.0*sqr(B);
    const scalar a0 = -A*B + sqr(B) + pow3(B);

    const scalar Q = (3.0*a1 - a2*a2)/9.0;
    const scalar Rl = (9.0*a2*a1 - 27.0*a0 - 2.0*a2*a2*a2)/54.0;

    const scalar Q3 = Q*Q*Q;
    const scalar D = Q3 + Rl*Rl;
//printf("Zmix_root00= %+5.10f,%+5.8f,%+5.8f\n", D,Q3,Rl*Rl); 
    scalar root1 = -1, root2 = -1, root3 = -1;
    if (D <= 0)
    {
        const scalar th = ::acos(Rl/sqrt(-Q3));
        const scalar qm = 2*sqrt(-Q);
        root1 = qm*cos(th/3) - a2/3;
        root2 = qm*cos((th + 2*constant::mathematical::pi)/3) - a2/3;
        root3 = qm*cos((th + 4*constant::mathematical::pi)/3) - a2/3;
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
	    S = pow(Rl + D05, 1.0/3.0);
	}
	else  
	{
	    S = -pow(-RlD, 1.0/3.0);
	}
        scalar Tl = 0;
//printf("Zmix_root200= %+5.8f,%+5.8f,%+5.8f\n", S,Rl,D05); 
        if (D05 > Rl)
        {
            Tl = -pow(mag(Rl - D05), 1.0/3.0);
        }
        else
        {
            Tl = pow(Rl - D05, 1.0/3.0);
        }
        scalar root1 = S + Tl - a2/3;
        rootz[0] = root1;
	rootz[1] = root1;
	rootz[2] = root1;
    } 
}




//can be used for phases or mixtures, depending on the Azin and Bzin, only one root
template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::Zmix 
(
    int phasei, //input
    scalar Azin,
    scalar Bzin
) const
{
    const scalar A = Azin;
    const scalar B = Bzin;

    const scalar a2 = B - 1.0;
    const scalar a1 = A - 2.0*B - 3.0*sqr(B);
    const scalar a0 = -A*B + sqr(B) + pow3(B);

    const scalar Q = (3.0*a1 - a2*a2)/9.0;
    const scalar Rl = (9.0*a2*a1 - 27.0*a0 - 2.0*a2*a2*a2)/54.0;

    const scalar Q3 = Q*Q*Q;
    const scalar D = Q3 + Rl*Rl;
    //printf("Zmix_root00= %+5.2i,%+5.10f,%+5.8f,%+5.8f\n", phasei,D,Azin,Bzin); 
    scalar root = -1;
    if (D <= 0)
    {
        const scalar th = ::acos(Rl/sqrt(-Q3));
        const scalar qm = 2*sqrt(-Q);
        const scalar r1 = qm*cos(th/3) - a2/3;
        const scalar r2 =
            qm*cos((th + 2*constant::mathematical::pi)/3) - a2/3;
        const scalar r3 =
            qm*cos((th + 4*constant::mathematical::pi)/3) - a2/3;

        // three root;
    	if (phasei == 0) 
    	{
            root = max(r1, max(r2, r3));
	}
	else
	{
	    if (Neg_rootFlag == 0)
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
        const scalar D05 = sqrt(D);
	const scalar RlD = Rl + D05;
        scalar S = 0.0;
	if (RlD > 0.0) 
	{
	    S = pow(Rl + D05, 1.0/3.0);
	}
	else  
	{
	    S = -pow(-RlD, 1.0/3.0);
	}
        scalar Tl = 0;
        if (D05 > Rl)
        {
            Tl = -pow(mag(Rl - D05), 1.0/3.0);
        }
        else
        {
            Tl = pow(Rl - D05, 1.0/3.0);
        }
        root = S + Tl - a2/3;
        //printf("Zmix_root22= %+5.2i,%+5.8f\n", phasei,root); 
    } 

    if (root <= 0.0) 
    {
	root = 1.0;
printf("Zmix_root4444444= %+5.2i,%+5.8f\n", phasei,root); 
FatalErrorInFunction
<< exit(FatalError);
    }
    return root;
}



template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::rhomix_phase // PY:phase mixture density
(
    int phasei, //input
    scalar p,	
    scalar T,
    scalar *comp_inp
) const
{
    scalar rhomixPhase = 0.0;
    scalar Amixz = 0.0;
    scalar Bmixz = 0.0;
    scalar Amix2[num_sp] = {0.0};
    scalar Bz[num_sp] = {0.0}; 
    this->ABZmix(p, T, comp_inp, Amixz, Bmixz, Amix2, Bz);
    const scalar Zmix = this->Zmix(phasei,Amixz, Bmixz);
    scalar mw_phase = this->mwmix(comp_inp);//kg/mol
    //RR: 10e+03*m3 Pa K-1 mol-1 
    rhomixPhase = mw_phase*p/(Zmix*RR*1.0e-03*T);// kg/m3
    //printf("rhomixPhase1= %+5.8f %+5.8f %+5.8f %+5.8f %+5.8f\n", p, T, Zmix, mw_phase, rhomixPhase); 

    return rhomixPhase;
}



template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::volmmix_phase //m3/mol
(
    int phasei, //input
    scalar p,
    scalar T,
    scalar *comp_inp
) const
{
    scalar volspecific = 0.0;
    scalar mw_phase = this->mwmix(comp_inp);//kg/mol
    const scalar rho_phase = this->rhomix_phase(phasei, p, T, comp_inp);  //kg/m3
    volspecific = mw_phase/rho_phase;  

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



// PY: fugacity coefficient of species in each phase
template<class Specie>
inline void Foam::PengRobinsonMix<Specie>::fugcoef_phase 
(
    int phasei, 	//input
    scalar p,		//input
    scalar T,		//input
    scalar *comp_inp, 	//input
    scalar *fugcoef 	//output
) const
{
    scalar Amixz = 0.0;
    scalar Bmixz = 0.0;
    scalar Zmix = 0.0;
    scalar Amix2[num_sp] = {0.0};
    scalar Bz[num_sp] = {0.0};
    scalar comp1[num_sp] = {0.0};
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
	    if (num_sp == 2)
	    {
		comp1[0] = comptrifile, comp1[1] = (1.0-comp1[0]); //binary
	    }
	    else if (num_sp == 3)    
	    {
		comp1[0] = comptrifile, comp1[1] = (1.0-comp1[0])/2.0, comp1[2] = (1.0-comp1[0])/2.0; //ternary
	    }
	    else if (num_sp == 4)    
	    {
		comp1[0] = comptrifile, comp1[1] = (1.0-comp1[0])/3.0, comp1[2] = (1.0-comp1[0])/3.0, comp1[3] = (1.0-comp1[0])/3.0;
	    }
	    else if (num_sp == 5)    
	    {
		comp1[0] = comptrifile, comp1[1] = (1.0-comp1[0])/4.0, comp1[2] = (1.0-comp1[0])/4.0, comp1[3] = (1.0-comp1[0])/4.0, comp1[4] = (1.0-comp1[0])/4.0;
	    }
	
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

    const scalar cc0 = 2.0*sqrt(2.0);
    const scalar cc1 = 1.0 + sqrt(2.0);
    const scalar cc2 = 1.0 - sqrt(2.0);

    scalar zfactor00 = Zmix;
    scalar zdiff = zfactor00 - Bmixz;

    if (zdiff < 0.0) 
    {
        zfactor00 = Bmixz + small8;
        zdiff = small8;
    }

    for (label spid = 0; spid < num_sp; spid++)
    {
        fugcoef[spid] = Bz[spid]/Bmixz*(zfactor00 - 1.0) - log(zdiff) 
	        - Amixz/(cc0*Bmixz) * (2.0*Amix2[spid]/Amixz - Bz[spid]/Bmixz) 
	        *log((zfactor00 + cc1*Bmixz)/(zfactor00 + cc2*Bmixz));

	fugcoef[spid] = exp(fugcoef[spid]);
    }    
//printf("fugcoef_phase_2222= %+5.15f,%+5.15f\n", fugcoef[0], fugcoef[1] );      
}


template<class Specie>
inline void Foam::PengRobinsonMix<Specie>::gibbs_phase 
(
    int phasei, //inp
    scalar p,
    scalar T,
    scalar *comp_inp,
    scalar &gibbsenergy  //outp
) const
{
    scalar fugcoef[num_sp] = {0.0};
    gibbsenergy = 0.0;
    this->fugcoef_phase(phasei, p, T, comp_inp, fugcoef);

    for (label spid = 0; spid < num_sp; spid++)
    {
        gibbsenergy += comp_inp[spid]*log(comp_inp[spid]*fugcoef[spid]);
    }  
//printf("gibbs_phase= %+5.2i,%+5.8f\n", phasei,gibbsenergy);  
}



// PY: fugacity coefficient of species in mixture
template<class Specie>
inline void Foam::PengRobinsonMix<Specie>::fugcoef_mixture 
(
    scalar p,     //inp
    scalar T,     //inp
    scalar *comp_inp, //inp 
    scalar *fugcoef_mix, //outp
    scalar& gibbsenergy  //outp
) const
{
    int igas = 0;
    int iliquid = 1;
    
    //scalar Zmixture  = 0.0;
    scalar Amixz = 0.0;
    scalar Bmixz = 0.0;
    scalar Amix2[num_sp] = {0.0};
    scalar Bz[num_sp] = {0.0};
    scalar Zmin = 0.0;
    scalar Zmax = 0.0;
    scalar gibbsenergyl = 0.0;
    scalar gibbsenergyg = 0.0;
    scalar fugcoef_mixl[num_sp] = {0.0};
    scalar fugcoef_mixg[num_sp] = {0.0};
    scalar comp1[num_sp] = {0.0};
    gibbsenergy = 0.0;   
    // Remark!!!!!!!!!!!!!! 
    // Reconsider the middle root for Z_ave, which is very important
    // For feed, we could select the Z_ave to make the system gibbs energy minimum
    this->ABZmix(p, T, comp_inp, Amixz, Bmixz, Amix2, Bz);
    Zmax = this->Zmix(igas, Amixz, Bmixz);
    Zmin = this->Zmix(iliquid, Amixz, Bmixz);

    scalar Zmin0 = Zmax;
    scalar comptrifile = 0.7;
    while (Zmin < 0.0)
    {
	if (num_sp == 2)
	    {
		comp1[0] = comptrifile, comp1[1] = (1.0-comp1[0]); //binary
	    }
	    else if (num_sp == 3)    
	    {
		comp1[0] = comptrifile, comp1[1] = (1.0-comp1[0])/2.0, comp1[2] = (1.0-comp1[0])/2.0; //ternary
	    }
	    else if (num_sp == 4)    
	    {
		comp1[0] = comptrifile, comp1[1] = (1.0-comp1[0])/3.0, comp1[2] = (1.0-comp1[0])/3.0, comp1[3] = (1.0-comp1[0])/3.0;
	    }
	    else if (num_sp == 5)    
	    {
		comp1[0] = comptrifile, comp1[1] = (1.0-comp1[0])/4.0, comp1[2] = (1.0-comp1[0])/4.0, comp1[3] = (1.0-comp1[0])/4.0, comp1[4] = (1.0-comp1[0])/4.0;
	    }
	this->ABZmix(p, T, comp1, Amixz, Bmixz, Amix2, Bz);      
	Zmin = this->Zmix(iliquid, Amixz, Bmixz);
	comptrifile = comptrifile + 0.01;
	if (comptrifile >= 1.0)
	{
	    Zmin = Zmin0;
	}
    }

    if (Zmin > Bmixz) 
    {
        this->fugcoef_phase(iliquid, p, T, comp_inp, fugcoef_mixl);
        for (label spid = 0; spid < num_sp; spid++)
        {
            if (fugcoef_mixl[spid] < 0.0)  
            {
                FatalErrorInFunction
                << "trying to insert a node with a wrong pointer to a chemPoint"
                << exit(FatalError);
            }
        } 
        this->gibbs_phase(iliquid, p, T, comp_inp, gibbsenergyl);
    } 
    else 
    {
        gibbsenergyl = 1.0e+10;
    } 

    if (Zmax > Bmixz) 
    {
        this->fugcoef_phase(igas, p, T, comp_inp, fugcoef_mixg);
        for (label spid = 0; spid < num_sp; spid++)
        {
            if (fugcoef_mixg[spid] < 0.0)  
            {
                FatalErrorInFunction
                << "trying to insert a node with a wrong pointer to a chemPoint"
                << exit(FatalError);
            }
        }      
        this->gibbs_phase(igas, p, T, comp_inp, gibbsenergyg);
    } 
    else 
    {
        gibbsenergyg = 1.0e+10;    
    }    

    if (gibbsenergyg < 0.0) 
    {
	gibbsenergyg = -gibbsenergyg;
    }
    if (gibbsenergyl < 0.0) 
    {
	gibbsenergyl = -gibbsenergyl;
    }

    if (gibbsenergyl >= gibbsenergyg) 
    {
        gibbsenergy = gibbsenergyg;
        for (label spid = 0; spid < num_sp; spid++)
        {
            fugcoef_mix[spid] =  fugcoef_mixl[spid];  
            //Zmixture = Zmin;   
        } 
    } 
    else 
    {
        gibbsenergy = gibbsenergyl;
        for (label spid = 0; spid < num_sp; spid++)
        {
            fugcoef_mix[spid] =  fugcoef_mixg[spid];     
            //Zmixture = Zmax;   
        } 
    }
//printf("fugcoef5555= %+5.8f,%+5.8f,%+5.8f\n", gibbsenergy,gibbsenergyl,gibbsenergyg);    
//printf("fugcoef6666= %+5.8f,%+5.8f\n", fugcoef_mix[0],fugcoef_mix[1]);    
//printf("fugcoef7777= %+5.8f,%+5.8f\n", fugcoef_mixl[0],fugcoef_mixl[1]);  
//printf("fugcoef8888= %+5.8f,%+5.8f\n", fugcoef_mixg[0],fugcoef_mixg[1]);      
}





/*
template<class Specie>
inline void Foam::PengRobinsonMix<Specie>::gibbs_mix 
(
    scalar p,
    scalar T,
    scalar *comp,
    scalar gibbsenergy 
) const
{
    phasei = 0;
    scalar fugcoef[num_sp] = 0.0;
    gibbsenergy = 0.0;

    this->fugcoef_mixture(p, T, comp, fugcoef, gibbsenergy);

    for (label spid = 0; spid < num_sp; spid++)
    {
        gibbsenergy += comp[spid]*log(comp[spid]*fugcoef[spid]);
    }    
}
*/


template<class Specie>
inline void Foam::PengRobinsonMix<Specie>::wilsoneq 
(
    scalar p, //inp
    scalar T, //inp
    scalarList *equl_constant  //outp
) const
{
    scalar dummy[num_sp] = {0.0};

    for (label spid = 0; spid < num_sp; spid++)
    {
        //dummy[spid] = (7.0/3.0)*(1.0 - Tc_sp[spid]/T)* (1.0 + omega_sp[spid]);
        //equl_constant[spid] = Pc_sp[spid]*pow(10.0,dummy[spid])/p ;

	dummy[spid] = 5.373*(1.0 + omega_sp[spid])*(1.0 - Tc_sp[spid]/T);
        equl_constant[spid] = (Pc_sp[spid]/p)*exp(dummy[spid]);
    }    
}


// PY: each species density
// Which is different with before, before we use the ro_phase * Yi_phase
// PY: I prefer to the present method
template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::rho_spd 
(
    int phasei, //input
    label spid, //input
    scalar p,	//input
    scalar T	//input
) const
{
    scalar Tr_sp[num_sp] = {0.0};

    scalar a[num_sp] = {0.0};
    scalar b[num_sp] = {0.0};
    scalar alpha[num_sp] = {0.0};
    scalar dadt_sp[num_sp] = {0.0};

    scalar A_spdz[num_sp] = {0.0};
    scalar B_spdz[num_sp] = {0.0};
    scalar Z_spd[num_sp] = {0.0};
    scalar AAx = {0.0};
    scalar BBx = {0.0};
    //for (label spid = 0; spid < num_sp; spid++)
    //{
	this->ab_spd(p, T, Tc_sp[spid], Pc_sp[spid], omega_sp[spid], alpha[spid], a[spid], b[spid], dadt_sp[spid]);
	Tr_sp[spid] = T/Tc_sp[spid];
 
        A_spdz[spid] = a[spid]*alpha[spid]*p/sqr(RR*1.0e-03*T);
        B_spdz[spid] = b[spid]*p/(RR*1.0e-03*T);
        AAx = A_spdz[spid];
        BBx = B_spdz[spid];
        Z_spd[spid] = this->Zmix(phasei, AAx, BBx);    
    //}  
    scalar rhosp = (mw_sp[spid]*p)/(Z_spd[spid]*RR*T);// kg/m3
    //printf("rho_sp= %+5.8f,%+5.8f,%+5.8f\n", mw_sp[spid],Z_spd[spid],rhosp);    
    return rhosp; 
}


template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::vol_spd // PY:phase mixture specific volume
(
    int phasei, //input
    label spid, //input
    scalar p, 	//input
    scalar T	//input
) const
{
    //scalar rho_sp[num_sp] = {0.0};
    //this->rho_spd(phasei, p, T, rho_sp);
    scalar rhosp = this->rho_spd(phasei, spid, p, T);
    scalar VVsp = mw_sp[spid]*1.0e-03/rhosp;//m3/mol
    return VVsp;
}




template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::Hmix // internal energy with unit of J/kg
(
    int phasei, //input
    scalar p,
    scalar T,
    scalar *comp_inp
) const
{

    //co2_, ch4_, o2_, N2_, h2o_; 
    //coefficients
		           	//co2		ch4	   	o2	  	 n2	   	h2o
    scalar a1s1000[5] = { 2.35677352e+00,  5.14911468e+00,  3.78246636e+00,  3.53100528e+00,  7.25575005e+01 };
    scalar a2s1000[5] = { 8.98459677e-03, -1.36622009e-02, -2.99673416e-03, -1.23660988e-04, -6.62445402e-01 };
    scalar a3s1000[5] = {-7.12356269e-06,  4.91453921e-05,  9.84730200e-06, -5.02999433e-07,  2.56198746e-03 };
    scalar a4s1000[5] = { 2.45919022e-09, -4.84246767e-08, -9.68129608e-09,  2.43530612e-09, -4.36591923e-06 };
    scalar a5s1000[5] = {-1.43699548e-13,  1.66603441e-11,  3.24372836e-12, -1.40881235e-12,  2.78178981e-09 };

    scalar h_ig[num_sp]  = { 0.0 };
    scalar H_ideal_fuel = 0.0;
    scalar H_ideal_amb = 0.0;

    int num_sp0 = 0, num_spend = num_sp;

    if (system_ == 1.0 || system_ == 5.0)
    {
	num_sp0 = 3, num_spend = 4;
    }
    else if (system_ == 6.0)
    {
	num_sp0 = 1, num_spend = 4;
    }
    else if (system_ == 7.0)
    {
	num_sp0 = 0, num_spend = 3;
    }
    else if (system_ == 8.0)
    {
	num_sp0 = 0, num_spend = 4;
    }
    else if (system_ == 9.0)
    {
	num_sp0 = 0, num_spend = 5;
    }
    else if (system_ == 10.0 )
    {
	num_sp0 = 0, num_spend = 4;
	a1s1000[3] = 7.25575005e+01;
	a2s1000[3] = -6.62445402e-01;
	a3s1000[3] = 2.56198746e-03;
	a4s1000[3] = -4.36591923e-06;
	a5s1000[3] = 2.78178981e-09;
    }

    for (label spid = num_sp0; spid < num_spend; spid++)
    {
        h_ig[spid] = RR*1.0e-03*T* (a1s1000[spid] + a2s1000[spid]*T/2.0 + a3s1000[spid]*sqr(T)/3.0 + a4s1000[spid]*pow3(T)/4.0 + a5s1000[spid]*pow4(T)/5.0);
	H_ideal_amb = H_ideal_amb + h_ig[spid]*comp_inp[spid];
        //printf("dH_mix0 = %+5.8f %+5.8f  %+5.8f \n", H_ideal_amb, h_ig[spid],comp_inp[spid]);  
    }    

    scalar hfuel_ig[num_sp]  = { 0.0 };
    if (system_ == 1.0 || system_ == 5.0 )
    {
					//C12H26	C7H16 
	scalar a1sfue1000[2] = { 2.13264480e+01,  1.11532484e+01 };
    	scalar a2sfue1000[2] = {-3.86394002e-02, -9.49415433e-03 };
    	scalar a3sfue1000[2] = { 3.99476113e-04,  1.95571181e-04 };
    	scalar a4sfue1000[2] = {-5.06681097e-07, -2.49752520e-07 };
    	scalar a5sfue1000[2] = { 2.00697878e-10,  9.84873213e-11 };

	if (system_ == 1.0)
	{
	    num_sp0 = 0, num_spend = 1;
	}
	else if (system_ == 5.0)
	{
	    num_sp0 = 1, num_spend = 2;
	}

	for (label spid = num_sp0; spid < num_spend; spid++)
	{
	    hfuel_ig[spid] = RR*1.0e-03*T*(a1sfue1000[spid] + a2sfue1000[spid]*T/2.0 + a3sfue1000[spid]*sqr(T)/3.0 + a4sfue1000[spid]*pow3(T)/4.0 + a5sfue1000[spid]*pow4(T)/5.0);
	    H_ideal_fuel = H_ideal_fuel + (hfuel_ig[spid]*comp_inp[0]);
	}    

 	H_ideal_fuel = H_ideal_fuel;
	H_ideal_amb = h_ig[3]*comp_inp[1];
    }
    scalar H_ideal = H_ideal_fuel + H_ideal_amb;
    scalar mw_phase = this->mwmix(comp_inp);//kg*mol-1

    scalar amix  = 0.0;
    scalar bmix  = 0.0;
    scalar dadt  = 0.0;
    scalar dadt2 = 0.0;

    scalar rho_phase = this->rhomix_phase(phasei, p, T, comp_inp);//kg*m-3
    scalar VV = mw_phase/rho_phase; //m3*mol-1, phase specific volume
    this->abmix(T, comp_inp, amix, bmix, dadt, dadt2);
    scalar CC1 =  2.0*sqrt(2.0), CC2 = 1.0+sqrt(2.0), CC3 = 1.0-sqrt(2.0);

    scalar xxx = (VV+bmix*CC3)/(VV+bmix*CC2);
    scalar KE = 0.0;
    if (xxx > 0)
    {
        KE = log((VV+bmix*CC3)/(VV+bmix*CC2))/(CC1*bmix);
    }
    else
    {
        KE = 0.0;
    }

    scalar detE = (amix-T*dadt)*KE;
    scalar dH_phase = p*VV-RR*1.0e-03*T+detE;  
    //scalar B_mix = bmix*p/(T*RR*1.0e-03);
    //scalar dH_phase = p*VV-RR*1.0e-03*T+(T*dadt-amix)/(CC1*bmix)*log((VV+CC2*bmix)/(VV+CC3*bmix));
    //scalar dH_phase = p*VV-RR*1.0e-03*T+(T*dadt-amix)/(CC1*bmix)*log((VV+CC2*B_mix)/(VV+CC3*B_mix));  
/*
	scalar xxx = p*VV;
	scalar yyy = (T*dadt-amix)/(CC1*bmix);
	scalar zzz = log((VV+CC2*bmix)/(VV+CC3*bmix)); 
    	printf("dH_mix0 = %+5.8f %+5.8f %+5.8f\n", xxx,yyy,zzz);  
    	printf("dH_mix00 = %+5.8f %+5.8f %+5.8f\n", amix,bmix,dadt);  
*/
    scalar HReal_phase = H_ideal + dH_phase;
    //printf("dH_mix0 = %+5.8f %+5.8f\n", comp_inp[0], comp_inp[1]);  
    //printf("dH_mix2 = %+5.8f %+5.8f  %+5.8f \n", HReal_phase, H_ideal, dH_phase);  

    //return HReal_phase/mw_phase;
    return HReal_phase;
}





template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::Hmix2 // internal energy with unit of J/kg
(
    int phasei, //input
    scalar p,
    scalar T,
    scalar *comp_inp
) const
{
    //co2_, ch4_, o2_, N2_, h2o_; 
    //coefficients
		    //co2	ch4	  o2	   n2	   h2o
    scalar BB[5] = { 7.54056, 8.003180, 6.96302, 6.95808, 7.97183 };
    scalar CC[5] = { 7.51625, 19.26330, 2.40013, 2.03952, 6.27078 };
    scalar DD[5] = { 1442.70, 2148.790, 2522.05, 1681.60, 2572.63 };
    scalar EE[5] = { 5.38023, 10.44230, 2.21752, 0.506863, 2.0501 };
    scalar FF[5] = { 647.502, 1017.070, 1154.15, 6535.68, 1156.72 };

    scalar h_ig[num_sp]  = { 0.0 };
    scalar H_ideal_fuel = 0.0;
    scalar H_ideal_amb = 0.0;
    int num_sp0 = 0, num_spend = num_sp;

    if (system_ == 1.0 || system_ == 5.0)
    {
	num_sp0 = 3, num_spend = 4;
    }
    else if (system_ == 6.0)
    {
	num_sp0 = 1, num_spend = 4;
    }
    else if (system_ == 7.0)
    {
	num_sp0 = 0, num_spend = 3;
    }
    else if (system_ == 8.0)
    {
	num_sp0 = 0, num_spend = 4;
    }
    else if (system_ == 9.0)
    {
	num_sp0 = 0, num_spend = 5;
    }
    else if (system_ == 10.0 )
    {
	num_sp0 = 0, num_spend = 4;
	BB[3] = 7.97183;
	CC[3] = 6.27078;
	DD[3] = 2572.63;
	EE[3] = 2.0501;
	FF[3] = 1156.72;
    }

    for (label spid = num_sp0; spid < num_spend; spid++)
    {
        h_ig[spid] = BB[spid]*T +CC[spid]*DD[spid]*(2.0/(exp(2.0*DD[spid]/T)-1.0))+ 2.0*EE[spid]*FF[spid]*(1.0/(exp(2.0*FF[spid]/T)+1.0));
	H_ideal_amb = H_ideal_amb + h_ig[spid]*comp_inp[1];
    }    
    H_ideal_amb = H_ideal_amb *4.18680;

    scalar hfuel_ig[num_sp]  = { 0.0 };
    if (system_ == 1.0 || system_ == 5.0 )
    {
			   //C12H26	C7H16   C20H42
	scalar BBf[3] = { 51.02700, 30.4029, 82.7586 };
	scalar CCf[3] = { 143.8250, 90.6941, 226.791 };
	scalar DDf[3] = { 1664.760, 1669.32, 1587.50 };
	scalar EEf[3] = { 109.4510, 63.2028, 175.085 };
	scalar FFf[3] = { 789.2970, 764.611, 769.249 };
	if (system_ == 1.0)
	{
	    num_sp0 = 0, num_spend = 1;
	}
	else if (system_ == 5.0)
	{
	    num_sp0 = 1, num_spend = 2;
	}

	for (label spid = num_sp0; spid < num_spend; spid++)
	{
	    hfuel_ig[spid] = BBf[spid]*T +CCf[spid]*DDf[spid]*(2.0/(exp(2.0*DDf[spid]/T)-1.0))+ 2.0*EEf[spid]*FFf[spid]*(1.0/(exp(2.0*FFf[spid]/T)+1.0));
	    H_ideal_fuel = H_ideal_fuel + (hfuel_ig[spid]*comp_inp[0]);
	}    
	H_ideal_fuel = H_ideal_fuel *4.18680;
	H_ideal_amb = h_ig[3]*comp_inp[1]*4.18680;
    }
    scalar H_ideal = H_ideal_fuel + H_ideal_amb;
    scalar mw_phase = this->mwmix(comp_inp);//kg*mol-1

    scalar amix  = 0.0;
    scalar bmix  = 0.0;
    scalar dadt  = 0.0;
    scalar dadt2 = 0.0;

    scalar rho_phase = this->rhomix_phase(phasei, p, T, comp_inp);//kg*m-3
    scalar VV = mw_phase/rho_phase; //m3*mol-1, phase specific volume
    this->abmix(T, comp_inp, amix, bmix, dadt, dadt2);
    scalar CC1 =  2.0*sqrt(2.0), CC2 = 1.0+sqrt(2.0), CC3 = 1.0-sqrt(2.0);

    scalar xxx = (VV+CC2*bmix)/(VV+CC3*bmix);
    scalar dH_phase = 0.0;
    if (xxx > 0)
    {
        dH_phase = p*VV-RR*1.0e-03*T+(T*dadt-amix)/(CC1*bmix)*log((VV+CC2*bmix)/(VV+CC3*bmix)); 
    }
    else
    {
        dH_phase = 0.0;
    }

    scalar HReal_phase = H_ideal + dH_phase;
    return HReal_phase/mw_phase;
}




template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::Hspd_mix // internal energy with unit of J/mol
(
    scalar p,
    scalar T,
    scalar *comp_inp,
    int spflag //Compositions in Mixture;
) const
{
    //co2_, ch4_, o2_, N2_, h2o_; 
    //coefficients
		           	//co2		ch4	   	o2	  	 n2	   	h2o
    scalar a1s1000[5] = { 2.35677352e+00,  5.14911468e+00,  3.78246636e+00,  3.53100528e+00,  7.25575005e+01 };
    scalar a2s1000[5] = { 8.98459677e-03, -1.36622009e-02, -2.99673416e-03, -1.23660988e-04, -6.62445402e-01 };
    scalar a3s1000[5] = {-7.12356269e-06,  4.91453921e-05,  9.84730200e-06, -5.02999433e-07,  2.56198746e-03 };
    scalar a4s1000[5] = { 2.45919022e-09, -4.84246767e-08, -9.68129608e-09,  2.43530612e-09, -4.36591923e-06 };
    scalar a5s1000[5] = {-1.43699548e-13,  1.66603441e-11,  3.24372836e-12, -1.40881235e-12,  2.78178981e-09 };

    scalar h_ig[num_sp] = { 0.0 };
    scalar H_ideal_sp[num_sp] = { 0.0 };

    int num_sp0 = 0, num_spend = num_sp;

    if (system_ == 1.0 || system_ == 5.0)
    {
	num_sp0 = 3, num_spend = 4;
    }
    else if (system_ == 6.0)
    {
	num_sp0 = 1, num_spend = 4;
    }
    else if (system_ == 7.0)
    {
	num_sp0 = 0, num_spend = 3;
    }
    else if (system_ == 8.0)
    {
	num_sp0 = 0, num_spend = 4;
    }
    else if (system_ == 9.0)
    {
	num_sp0 = 0, num_spend = 5;
    }
    else if (system_ == 10.0 )
    {
	num_sp0 = 0, num_spend = 4;
	a1s1000[3] = 7.25575005e+01;
	a2s1000[3] = -6.62445402e-01;
	a3s1000[3] = 2.56198746e-03;
	a4s1000[3] = -4.36591923e-06;
	a5s1000[3] = 2.78178981e-09;
    }

    for (label spid = num_sp0; spid < num_spend; spid++)
    {
        h_ig[spid] = RR*1.0e-03*T* (a1s1000[spid] + a2s1000[spid]*T/2.0 + a3s1000[spid]*sqr(T)/3.0 + a4s1000[spid]*pow3(T)/4.0 + a5s1000[spid]*pow4(T)/5.0);
	H_ideal_sp[spid] = h_ig[spid];
    }    

    scalar a_sp[num_sp] = {0.0};
    scalar b_sp[num_sp] = {0.0};
    scalar dadt_sp[num_sp] = {0.0};
    scalar alpha_sp[num_sp] = {0.0};
    scalar dH_sp[num_sp] = {0.0};
    scalar VV_sp[num_sp] = {0.0};
    scalar H_real_sp[num_sp] = {0.0};
    scalar H_real_presudo[num_sp] = {0.0};
    scalar VVsp = 0.0;
    scalar CC1 =  2.0*sqrt(2.0), CC2 = 1.0+sqrt(2.0), CC3 = 1.0-sqrt(2.0);  

// ************************************************************************* //
    scalar hfuel_ig[num_sp]  = { 0.0 };
    scalar H_realsp = 0.0;
    if ( system_ == 1.0 || system_ == 5.0 )
    {
					//C12H26	C7H16 
	scalar a1sfue1000[2] = { 2.13264480e+01,  1.11532484e+01 };
    	scalar a2sfue1000[2] = {-3.86394002e-02, -9.49415433e-03 };
    	scalar a3sfue1000[2] = { 3.99476113e-04,  1.95571181e-04 };
    	scalar a4sfue1000[2] = {-5.06681097e-07, -2.49752520e-07 };
    	scalar a5sfue1000[2] = { 2.00697878e-10,  9.84873213e-11 };

	if (system_ == 1.0)
	{
	    num_sp0 = 0, num_spend = 1;
	}
	else if (system_ == 5.0)
	{
	    num_sp0 = 1, num_spend = 2;
	}

	for (label spid = num_sp0; spid < num_spend; spid++)
	{
	    hfuel_ig[spid] = RR*1.0e-03*T*(a1sfue1000[spid] + a2sfue1000[spid]*T/2.0 + a3sfue1000[spid]*sqr(T)/3.0 + a4sfue1000[spid]*pow3(T)/4.0 + a5sfue1000[spid]*pow4(T)/5.0);
	}    
	H_ideal_sp[0] = hfuel_ig[num_sp0];
	H_ideal_sp[1] = h_ig[3];

    	//printf("Hspd_mix00000\n");  
        for (label spid = 0; spid < num_sp; spid++)
        { 
	    if (spid == 0)
	    {
	    	VVsp = this->vol_spd(1, spid, p, T);//m3/mol
	    }
	    else
	    {
	    	VVsp = this->vol_spd(0, spid, p, T);
	    }
	    VV_sp[spid] = VVsp;//m3/mol
	    this->ab_spd(p, T, Tc_sp[spid], Pc_sp[spid], omega_sp[spid], alpha_sp[spid],a_sp[spid], b_sp[spid], dadt_sp[spid]);
	    //unit: RR*1.0e-03 = J/Kmol; RR*1.0e-03*T = J/mol
	    //unit: VV_sp = m3/mol; p*VV_sp = J/mol
	    dH_sp[spid] = p*VV_sp[spid] - RR*1.0e-03*T + (T*dadt_sp[spid]-a_sp[spid])/(CC1*b_sp[spid])*log((VV_sp[spid]+CC2*b_sp[spid])/(VV_sp[spid]+CC3*b_sp[spid]));
	    H_real_sp[spid] = (H_ideal_sp[spid] + 0.1*dH_sp[spid]);//J/mol
    	    //printf("Hspd_mix11 = %+5.8f %+5.8f %+5.8f %+5.8f\n", H_real_sp[spid],H_ideal_sp[spid],dH_sp[spid],VV_sp[spid]);  
        }

        if (spflag == 0)
        {
	    H_realsp = H_real_sp[0]; 
        }
        else if (spflag == 1)
        {
	    H_realsp = H_real_sp[1];
        }
    }
// ************************************************************************* //

    if ( system_ == 2.0 || system_ == 3.0)
    {
        H_realsp = 0.0;
        for (label spid = 0; spid < num_sp; spid++)
        { 
	    VVsp = this->vol_spd(1, spid, p, T);//m3/mol
	    VV_sp[spid] = VVsp;
	    //this->vol_spd(0, p, T, VVsp);//m3/mol
	    this->ab_spd(p, T, Tc_sp[spid], Pc_sp[spid], omega_sp[spid], alpha_sp[spid],a_sp[spid], b_sp[spid], dadt_sp[spid]);
	    dH_sp[spid] = p*VV_sp[spid] - RR*1.0e-03*T + (T*dadt_sp[spid]-a_sp[spid])/(CC1*b_sp[spid])*log((VV_sp[spid]+CC2*b_sp[spid])/(VV_sp[spid]+CC3*b_sp[spid]));
	    H_real_presudo[spid] = (H_ideal_sp[spid] + dH_sp[spid]);//J/mol
//printf("Hspd_mix22 = %+5.8f %+5.8f\n", H_ideal_sp[spid],dH_sp[spid]);       
        }
        if (spflag == 0)
        {
	    H_realsp = H_real_presudo[0]; 
    	}
    	else if (spflag == 1)
    	{
	    H_realsp = H_real_presudo[1]; 
    	}
    }
// ************************************************************************* //

    if ( system_ == 7.0 || system_ == 8.0 || system_ == 9.0 || system_ == 10.0)
    {
        H_realsp = 0.0;
    	scalar H_real_res = 0.0;
    	for (label spid = 0; spid < num_sp; spid++)
    	{ 
	    VVsp = this->vol_spd(1, spid, p, T);//m3/mol
	    VV_sp[spid] = VVsp;
	    this->ab_spd(p, T, Tc_sp[spid], Pc_sp[spid], omega_sp[spid], alpha_sp[spid],a_sp[spid], b_sp[spid], dadt_sp[spid]);
	    scalar xxx = ((VV_sp[spid]+CC2*b_sp[spid])/(VV_sp[spid]+CC3*b_sp[spid]));
	    if (xxx > 0)
	    {
		dH_sp[spid] = p*VV_sp[spid] - RR*1.0e-03*T + (T*dadt_sp[spid]-a_sp[spid])/(CC1*b_sp[spid])*log((VV_sp[spid]+CC2*b_sp[spid])/(VV_sp[spid]+CC3*b_sp[spid]));
	    }
	    else
	    {
		dH_sp[spid] = 0.0;
	    }

	    H_real_presudo[spid] = (H_ideal_sp[spid] + dH_sp[spid]);//J/mol
    	}
    	for (label spid = 2; spid < num_sp; spid++)
    	{ 
	    H_real_res = H_real_res + H_real_presudo[spid]*comp_inp[spid];
    	} 
	H_real_res = H_real_res + H_real_presudo[0]*comp_inp[0];

    	if (spflag == 0)
   	{
	    H_realsp = H_real_presudo[1]*comp_inp[1]; 
    	}
    	else if (spflag == 1)
    	{
	    H_realsp = H_real_res;
    	}
    }

//printf("Hspd_mix22 = %+5.8f \n", H_realsp);       
    return H_realsp;
}



template<class Specie>
inline void Foam::PengRobinsonMix<Specie>::ab_spd 
(
    scalar p,	//input
    scalar T,	//input
    scalar Tcsp,//input
    scalar Pcsp,//input
    scalar omegasp,//input

    scalar &alphaz,
    scalar &az,
    scalar &bz,
    scalar &dadtz
) const
{

    scalar Trsp = T/Tcsp; 
    scalar m_omega_z = 0.0; 
       
    if (omegasp > 0.49) 
    {
        m_omega_z = 0.379642 + 1.48503*omegasp - 0.164423*sqr(omegasp) + 0.016666*pow3(omegasp);
    }
    else
    {
        m_omega_z = 0.37464 + 1.54226*omegasp - 0.26992*sqr(omegasp);
    }

    alphaz = sqr(1.0 + omegasp * (1.0 - sqrt(Trsp)));
    az = 0.45724 * alphaz * sqr(RR*1.0e-03*Tcsp)/Pcsp;
    bz = 0.07780 * RR * 1.0e-03* Tcsp/Pcsp;
    dadtz = ( 2.0*0.45725*sqr(RR*1.0e-03*Tcsp)/Pcsp ) * sqrt(alphaz) * (-m_omega_z * 1.0/(2.0*sqrt(T*Tcsp)));  
}




template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::siemix_phase // internal energy with unit of J/kg
(
    int phasei, //input
    scalar p,
    scalar T,
    scalar *comp_phase
) const
{
    scalar amix  = 0.0;
    scalar bmix  = 0.0;
    scalar dadt  = 0.0;
    scalar dadt2 = 0.0;
    scalar mw_phase = this->mwmix(comp_phase);//kg*mol-1
    scalar rho_phase = this->rhomix_phase(phasei, p, T, comp_phase);//kg*m-3

    const scalar cc1 = 1.0 + sqrt(2.0);
    const scalar cc2 = 1.0 - sqrt(2.0);

    scalar VV = mw_phase/rho_phase; //m3*mol-1, phase specific volume

    this->abmix(T, comp_phase, amix, bmix, dadt, dadt2); 
    scalar xxx = ((VV+cc1*bmix)/(VV+cc2*bmix));
    scalar dsie_phase = 0.0;
    if (xxx > 0)
    {
        dsie_phase = ((T*dadt-amix)/(2.0*sqrt(2.0)*bmix)) * log((VV+cc1*bmix)/(VV+cc2*bmix));
    }
    else
    {
        dsie_phase = 0.0;
    }

    dsie_phase = dsie_phase/mw_phase; // J/kg
    // Based on Eq.12, Internal energy is estimated based on ideal(close to absolute 0 T) + departure
    scalar sieReal_phase = 0.0;
    scalar sieig = this->sieig_phase(T, comp_phase);
    sieReal_phase = dsie_phase + sieig;

    return sieReal_phase;
}



template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::cvmix_phase // PY:phase mixture heat capability
(
    int phasei, //input
    scalar p,
    scalar T,
    scalar *comp_phase
) const
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
	cv_phase = delta_sie/deltat;
    }
    return cv_phase;
}


template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::cpmix_phase // PY:phase mixture heat capability
(
    int phasei, //input
    scalar p,
    scalar T,
    scalar *comp_phase
) const
{
    scalar cv_phase = this->cvmix_phase(phasei, p, T, comp_phase);

    scalar amix = 0.0;
    scalar bmix = 0.0;
    scalar dadt = 0.0;
    scalar dadt2 = 0.0;

    this->abmix(T, comp_phase, amix, bmix, dadt, dadt2); 
    scalar mw_phase = this->mwmix(comp_phase);
    scalar rho_phase = this->rhomix_phase(phasei, p, T, comp_phase);
    scalar VV = mw_phase/rho_phase;
    scalar den  = sqr(VV) + 2.0*bmix*VV-sqr(bmix);
    scalar dpdt = RR*1.0e-03/(VV-bmix) - dadt/den;
    scalar dpdv = -RR*1.0e-03*T/sqr(VV-bmix) + 2.0*amix*(VV+bmix)/sqr(den); 
    scalar ct2  = -dpdv*mw_phase/sqr(rho_phase);
    scalar cp_phase = cv_phase + T*sqr(dpdt)/(ct2*sqr(rho_phase));
    //printf("VV,cp_phase= %+5.8f %+5.8f %+5.8f %+5.8f\n",VV,rho_phase,cv_phase,cp_phase);  
    //scalar cp_phase = cv_phase;

    if (cp_phase <= 100)
    {
        cp_phase = 100;  
    }

    return cp_phase;
}



template<class Specie>
inline void Foam::PengRobinsonMix<Specie>::soundspeedmix // m/s
(
    int phasei, //input
    scalar p,
    scalar T,
    scalar *comp_inp,
    scalar& ct_phase, //output
    scalar& cs_phase
) const
{
    scalar cv_phase = this->cvmix_phase(phasei, p, T, comp_inp);
    scalar mw_phase = this->mwmix(comp_inp);
    scalar rho_phase = this->rhomix_phase(phasei, p, T, comp_inp);
    scalar VV = mw_phase/rho_phase;
    scalar amix = 0.0, bmix = 0.0;
    scalar dadt = 0.0, dadt2 = 0.0;
    this->abmix(T, comp_inp, amix, bmix, dadt, dadt2); 
    scalar den  = sqr(VV)+2.0*bmix*VV-sqr(bmix);
    scalar dpdt = RR*1.0e-03/(VV-bmix) - dadt/den;
    scalar dpdv = -RR*1.0e-03*T/sqr(VV-bmix) + 2.0*amix*(VV+bmix)/sqr(den); 
    scalar ct2  = -dpdv*mw_phase/sqr(rho_phase);

    // PY: Deal with the negative sound speed
    if (ct2<0.0)
    {
        ct2 = 1.0e-04;
    }
    scalar cs2  = ct2+T*sqr(dpdt)/(cv_phase*sqr(rho_phase));
    
    cs_phase = sqrt(cs2);
    ct_phase = sqrt(ct2);
}


template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::OAD217_function 
(
    scalar XX,
    scalar YY
) const
{
    int test0 = 1;
    int OAD217 = 0;

    scalar dumdif = (XX-YY);
    scalar dumsum = (XX+YY)/2.0;
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
   
    if (dumdif_abs < (test0*dumsum_abs) || dumdif == 0.0)
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


template<class Specie>
inline void Foam::PengRobinsonMix<Specie>::OAD013_function 
(
    scalar *comp_in,  // input
    scalar *comp_out,  // output
    int& ierr  //outpur
) const
{
    ierr = 0;
    scalar sum = 0.0;
    for (label spid = 0; spid < num_sp; spid++)
    {
        sum += comp_in[spid];
    }
        
    int OAD217_val = this->OAD217_function(sum, 0.0);    
    if (OAD217_val == 0) 
    {
        ierr = 1;
        return;
    }

    for (label spid = 0; spid < num_sp; spid++)
    {
        comp_out[spid] = comp_in[spid]/sum;
    }
}




template<class Specie>
inline void Foam::PengRobinsonMix<Specie>::OAD012_function
(
    scalar p,		//input
    scalar T,		//input
    scalar *comp_inp,	//input
    scalar *equalconstant,//input
    scalar *comp_liq,	//output
    scalar *comp_gas,	//output
    scalar& vaporf,	//output
    int& ierr		//output
) const
{
    int miter = 1;

    ierr = 0;
    scalar crita1 = -1.234567e+10;   
    scalar crita2 = 1.234567e+10;   
    scalar tolerance = 1.0e-14;  
    scalar coeff1 = 0.1; 
    scalar coeff2 = 0.5; //if the output vaporf=0.5, try to enlarge crita2 and give more time for tpn convergence;
    int nitmax = 500; //100
    scalar dumm = 0.0;
    scalar comp_liq0[num_sp] = {0.0};
    scalar comp_gas0[num_sp] = {0.0};

    vaporf = 0.9;
    ierr = 0;

    scalar ad = -crita1, ag = crita1, tetad = -crita1, tetag = crita1;
    //which can be replaced by input comp_gas or comp_liq
    for (label spid = 0; spid < num_sp; spid++)
    {
        comp_liq[spid] = equalconstant[spid] - 1.0;
//printf("OAD012function0******= %+5.8f %+5.8f \n", comp_liq[spid],equalconstant[spid]);  
        if (comp_inp[spid] > 0.0)
        {
            if (comp_liq[spid] > 0.0)
            {
                ag = max(ag, -1.0/comp_liq[spid]);
                dumm = (comp_inp[spid]-1.0)/comp_liq[spid] + comp_inp[spid];
                tetag = max(tetag, dumm);
            } 
            else if (comp_liq[spid] < 0.0)
            {
                ad = min(ad, -1.0/comp_liq[spid]);
                dumm = (comp_inp[spid]-1.0)/comp_liq[spid];
                tetad = min(tetad, dumm);
            }  
            // Remark !!!!! what we will do for comp_liq[spid] = 0.0  
        }
    } 
//printf("OAD012function1******= %+5.8f %+5.8f %+5.8f %+5.8f %+5.8f\n", ag,ad,dumm,tetag,tetad);  
    if (ad > crita2 && ag < -crita2) //Not very sure, no physics, please double check
    {
        vaporf = 0.5;
        for (label spid = 0; spid < num_sp; spid++)
        {
           comp_liq[spid] = comp_inp[spid];
           comp_gas[spid] = comp_inp[spid];
        } 
        ierr = 10; //change to 2
    }
    else if (ad == -crita1)//Not very sure, no physics, please double check
    {
        vaporf = 0.9999999;
        for (label spid = 0; spid < num_sp; spid++)
        {
            comp_gas[spid] = comp_inp[spid];
            comp_liq0[spid] = comp_gas[spid]/equalconstant[spid];
        }
        this->OAD013_function(comp_liq0, comp_liq, ierr);
        ierr = 4; //change to 0
    }
    else if (ag == crita1)//Not very sure, no physics, please double check
    {
        vaporf = 1.0e-07;
        for (label spid = 0; spid < num_sp; spid++)
        {
            comp_liq[spid] = comp_inp[spid];
            comp_gas0[spid] = comp_liq[spid]*equalconstant[spid];
        }        
        this->OAD013_function(comp_gas0, comp_gas, ierr);
        ierr = 3;  //change to 1
    }
    else 
    {
	// the initialization indicates the state is far from two-phase regime;
        if (tetag < 0.0) 
        {
            scalar f0 = -1.0;
            for (label spid = 0; spid < num_sp; spid++)
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
          vaporf = (tetag + tetad)/2.0;
          scalar dx = tetad - tetag;
//printf("OAD012function_00000= %+5.8f %+5.8f\n", vaporf,dx);  
	  //solve Rachford-rice equation 
          int niter = 0;  
	  scalar fonc = 0.0, deri = 0.0, test = 0.0;
          while (niter < nitmax)
          {   
	      fonc = 0.0, deri = 0.0;
              for (label spid = 0; spid < num_sp; spid++)
	      {
                  dumm = comp_liq[spid] / (1.0+comp_liq[spid]*vaporf);
                  fonc = fonc + dumm*comp_inp[spid];
                  deri = deri - sqr(dumm)*comp_inp[spid];
              }
	      dumm =(vaporf - ag) * (ad - vaporf);
              deri = dumm*deri - (2.0*vaporf-ag-ad)*fonc;
              fonc = fonc * dumm;
              test = ((vaporf-tetad)*deri-fonc) * ((vaporf-tetag)*deri-fonc);

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
	          dx = coeff1 * (tetad-tetag); //if TPn flash converged not very well, reduce dx
                  vaporf = tetag + dx;
              }
              else
	      { 
                  dx = coeff2 * fonc/deri;
                  vaporf = vaporf-dx;
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

        for (label spid = 0; spid < num_sp; spid++)
        {   
            comp_liq[spid] = comp_inp[spid]/(1.0 + comp_liq[spid]*vaporf);
            comp_gas[spid] = equalconstant[spid]*comp_liq[spid];
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


//specially for single phase
template<class Specie>
inline void Foam::PengRobinsonMix<Specie>::OAD012_singlephase
(
    scalar p,		//input
    scalar T,		//input
    scalar *comp_inp,	//input
    scalar *equalconstant,//input
    scalar  tpdvalue,
    scalar &vaporfra,
    scalar *comp_liq,	//output
    scalar *comp_gas	//output
) const
{
    if (tpdvalue == 0.0)//Not very sure, no physics, please double check
    {
	vaporfra = 0.99;

	if (system_ == 1 || system_ == 5)
	{
	    comp_liq[1] = 0.1, comp_gas[1] = 0.999999;
	    comp_liq[0] = 0.9, comp_gas[0] = 1.0E-06;
	}
//printf("OAD012_singlephase0 = %+5.8f %+5.8f %+5.8f\n",vaporfra,comp_liq[0], comp_gas[0]);
    }
    else if (tpdvalue == 1.0)//Not very sure, no physics, please double check
    {
	vaporfra = 0.01;
	comp_liq[0] = 0.999999, comp_gas[0] = 0.1;
	comp_liq[1] = 1.0e-06, comp_gas[1] = 0.9;
//printf("OAD012_singlephase1 = %+5.8f %+5.8f %+5.8f\n",vaporfra,comp_liq[0], comp_gas[0]);
    }
}
 

template<class Specie>
inline void Foam::PengRobinsonMix<Specie>::GDEM_SS // PY:accelaration
(
    scalar *vt,
    scalar *vt1,
    scalar *vt2,
    scalar *vout 
) const
{
    // compute error
    scalar err_sum1 = 0.0, err_sum0 = 0.0;
    scalar error_t1[num_sp] = {0.0}, error_t2[num_sp] = {0.0};
    for (label spid = 0; spid < num_sp; spid++)
    {
        error_t1[spid] = vt1[spid] - vt[spid];
        error_t2[spid] = vt2[spid] - vt1[spid];
        err_sum1 += error_t2[spid]*error_t1[spid];
        err_sum0 += error_t1[spid]*error_t1[spid];
    }

    // compute eigenvalue
    scalar lamda = err_sum1/err_sum0;

    // ideal solution
    for (label spid = 0; spid < num_sp; spid++)
    {
        vout[spid] = vt2[spid] + lamda*err_sum1/(1.0-lamda);
    }
}






template<class Specie>
inline void Foam::PengRobinsonMix<Specie>::TPn_flash // PY:with a course initialization
(
    scalar p,	//input
    scalar T,	//input
    scalar *comp,//input
    scalar *comp_liq,		//output
    scalar *comp_gas,		//output
    scalar& vaporfra,		//output
    scalar *equalconstant	//output
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

    equalconstant[num_sp] = {0.0};
    scalar gibbsenergy = 0.0; 
    scalar fugcoef_liq[num_sp] = {0.0}, fugcoef_gas[num_sp] = {0.0};//, fugcoef_mix[num_sp] = {0.0};
    scalar comp_liq0[num_sp] = {0.0}, comp_gas0[num_sp] = {0.0};//, fugcoef_mix[num_sp] = {0.0};

    for (label spid = 0; spid < num_sp; spid++)
    {
        comp_liq[spid] = comp[spid];//use this in PT_x flash, cancle this in hn flash
	comp_gas[spid] = comp[spid];//use this in PT_x flash, cancle this in hn flash
        comp_liq0[spid] = comp_liq[spid];
	comp_gas0[spid] = comp_gas[spid];
    }

    //scalar tpdvalue = 0;
    //this->TPD(p, T, comp, comp_liq, comp_gas, tpdvalue);
    this->fugcoef_phase(iliquid, p, T, comp_liq, fugcoef_liq);
    this->fugcoef_phase(igas,    p, T, comp_gas, fugcoef_gas); 
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

    scalar re_iter1[num_sp] = {1.0};
    scalar gibbs_iter = gibbsenergy;
    scalar vaporfra_iter1 = 0.0;
    scalar vaporfra_iter0 =0.0;
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
	    for (label spid = 0; spid < num_sp; spid++)
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

        if (vaporfra > 0.99999999 && flag_tp > 0) 
	{
	    vaporfra = 0.99999999;
	    for (label spid = 0; spid < num_sp; spid++)
	    {
            	comp_liq[spid] = comp_liq0[spid];
	 	comp_gas[spid] = comp_gas0[spid];
	    }
	    break;
	}
        if (vaporfra < 1.0e-08 && flag_tp > 0) 
	{
	    vaporfra = 1.0e-08;
	    for (label spid = 0; spid < num_sp; spid++)
	    {
            	comp_liq[spid] = comp_liq0[spid];
	 	comp_gas[spid] = comp_gas0[spid];
	    }
	    break;
	}

	vaporfra_iter1 = vaporfra;
        if ((ierr == 10 ) && iter > 2) 
	{
    	    //this->TPD(p, T, comp, comp_liq, comp_gas, tpdvalue);
	    this->fugcoef_phase(iliquid, p, T, comp_liq, fugcoef_liq);
	    this->fugcoef_phase(igas,    p, T, comp_gas, fugcoef_gas);
	    //this->fugcoef_mixture(p, T, comp, fugcoef_mix, gibbsenergy);
	    for (label spid = 0; spid < num_sp; spid++)
	    {
        	equalconstant[spid] = fugcoef_liq[spid]/fugcoef_gas[spid];
	    }
	    this->OAD012_function(p, T, comp, equalconstant, comp_liq, comp_gas, vaporfra, ierr);
	    // Here, for tp diagram, please reconsider it
	    for (label spid = 0; spid < num_sp; spid++)
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
	this->fugcoef_phase(igas,    p, T, comp_gas, fugcoef_gas);
        //printf("TPflash_OAD012_function3 = %+5.8f %+5.8f %+5.8f %+5.8f %+5.8f %+5.8f\n",fugcoef_liq[0], fugcoef_gas[0],equalconstant[0],equalconstant[1],comp_liq[0],comp_gas[0]);
        scalar errr = 0.0;
	scalar re_iter0[num_sp] = {0.0};
	for (label spid = 0; spid < num_sp; spid++)
	{
	    re_iter0[spid] = re_iter1[spid];
            fugcoef_gas[spid] = comp_gas[spid]*fugcoef_gas[spid]*p; 
            fugcoef_liq[spid] = comp_liq[spid]*fugcoef_liq[spid]*p; 
            errr = errr + sqr(fugcoef_liq[spid]/fugcoef_gas[spid]-1.0); 
            equalconstant[spid] = equalconstant[spid]*((fugcoef_liq[spid]/fugcoef_gas[spid])*coeff_tpn);
            re_iter1[spid] = fugcoef_liq[spid]/fugcoef_gas[spid];
	}
        epsilon1  = errr; 
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
	scalar epsilon2  = 1.0e+10;
	if (errr >= 1.0e+05 && iter > 10) 
	{
	    this->fugcoef_phase(iliquid, p, T, comp_liq, fugcoef_liq);
	    this->fugcoef_phase(igas,    p, T, comp_gas, fugcoef_gas);
	    //this->fugcoef_mixture(p, T, comp, fugcoef_mix, gibbsenergy);
	    for (label spid = 0; spid < num_sp; spid++)
	    {
        	equalconstant[spid] = fugcoef_liq[spid]/fugcoef_gas[spid];
	    }

	    this->OAD012_function(p, T, comp, equalconstant, comp_liq, comp_gas, vaporfra, ierr);
	    scalar errr2 = 0.0;
	    for (label spid = 0; spid < num_sp; spid++)
	    {
	    	re_iter0[spid] = re_iter1[spid];
            	fugcoef_gas[spid] = comp_gas[spid]*fugcoef_gas[spid]*p; 
            	fugcoef_liq[spid] = comp_liq[spid]*fugcoef_liq[spid]*p; 
            	errr2 = errr2 + sqr(fugcoef_liq[spid]/fugcoef_gas[spid]-1.0); 
            	equalconstant[spid] = equalconstant[spid]*(fugcoef_liq[spid]/fugcoef_gas[spid]);
            	re_iter1[spid] = fugcoef_liq[spid]/fugcoef_gas[spid];
	    }
	    epsilon2  = errr2; 
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
                scalar lamda_coef[num_sp] = {0.0};
		for (label spid = 0; spid < num_sp; spid++)
		{
		    sum_iter0 += sqr(re_iter0[spid]-1.0); 
	            sum_iter1 += sqr(re_iter1[spid]-1.0); 
		}

		scalar ratioiter = sum_iter1/sum_iter0;
		scalar detiter   = vaporfra_iter1-vaporfra_iter0;
		if (ratioiter>0.8 && detiter<0.1 && sum_iter1>1.0e-5 && sum_iter1<1.0e-3 && vaporfra_iter1>1.0e-07 && vaporfra_iter1<0.9999999) 
		{
		    for (label spid = 0; spid < num_sp; spid++)
		    {
			lamda_coef[spid] = (re_iter0[spid]-1.0)/(re_iter0[spid]-re_iter1[spid]);
			equalconstant[spid] = equalconstant[spid]*pow(re_iter1[spid],lamda_coef[spid]);
		    }
		}
	    }
        }

	scalar gibbsenergyl = 0.0, gibbsenergyg = 0.0;
	this->gibbs_phase(iliquid, p, T, comp_liq, gibbsenergyl);
	this->gibbs_phase(igas,    p, T, comp_gas, gibbsenergyg);

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





template<class Specie>
inline void Foam::PengRobinsonMix<Specie>::TPn_simple // PY:with a course initialization
(
    scalar p,	//input
    scalar T,	//input
    scalar *comp,//input
    scalar *comp_liq,		//output
    scalar *comp_gas,		//output
    scalar& vaporfra,		//output
    scalar *equalconstant	//output
) const
{
    int ierr = 0;
    equalconstant[num_sp] = {0.0};
    this->wilsoneq(p, T, equalconstant);
    scalar vaporf0 = vaporfra;
    this->OAD012_function(p, T, comp, equalconstant, comp_liq, comp_gas, vaporfra, ierr);
//printf("TPn_simple************ = %5.6f %5.6f %5.3i\n",vaporf0,vaporfra,ierr); 
    scalar detvap = vaporfra - vaporf0;
    if (detvap < 0.0)
    {
	detvap = -detvap;
    }
    if (detvap > 0.1)
    {
	//vaporfra = vaporf0 + 0.01*(vaporfra);
    }
}




template<class Specie>
inline void Foam::PengRobinsonMix<Specie>::TPn_output 
(
    scalar p, //input
    scalar T, 
    scalar *comp, 
    scalar *comp_liq, 
    scalar *comp_gas,
    scalar *svphase, //output
    scalar *rhophase,
    scalar *zphase,
    scalar *cpphase,
    scalar *cvphase,
    scalar *ctphase,
    scalar *csphase,
    scalar *siephase,
    scalar *gamphase,
    scalar *dpdtphase,
    scalar *dvdtphase
) const
{
    // Remark
    // Initialization for TPn flash
    // how to get the good vaporfra and comp_liq and comp_gas
    // call TP_flash_ini(press,temp,z_spfre0,z_spfue0,comp_liq,comp_gas,vaporfra,tpdtest0,ierr);
    // or call tpdss(press,temp,comp,xxx,yyy,tpdtest0);  
    const int num_phase = 2;
    const int iliquid = 1, igas = 0;
    int phasei;
    scalar vaporfra_virtual = 1.0e-10;

    scalar vaporfra = 0.5;
    scalar equalconstant[num_sp] = 0.0;
    comp_liq[num_sp] = {comp[num_sp]};
    comp_gas[num_sp] = {comp[num_sp]};
    this->TPn_flash(p, T, comp, comp_liq, comp_gas, vaporfra, equalconstant);

    for (label spid = 0; spid < num_sp; spid++)
    {
        if (vaporfra <= 0.0)
        {
            vaporfra=1.0e-10;
            comp_liq[spid] = comp[spid];
            comp_gas[spid] = comp[spid];
        }
        else if (vaporfra >= 1.0)
        {
            vaporfra = 0.999999999;
            comp_liq[spid] = comp[spid];
            comp_gas[spid] = comp[spid];
        }          
    } 

    scalar mw_gas = this->mwmix(comp_gas);
    scalar mw_liq = this->mwmix(comp_liq);

    scalar sv_gas = this->volmmix_phase(0, p, T, comp_gas);  //gas phase specific volume
    scalar sv_liq = this->volmmix_phase(1, p, T, comp_liq);  //liquid phase specific volume

    scalar rho_gas = mw_gas/sv_gas;
    scalar rho_liq = mw_liq/sv_liq;   

    scalar z_gas = p*sv_gas/(RR*1.0e-03*T); //why donot use the this->Zmix
    scalar z_liq = p*sv_liq/(RR*1.0e-03*T);

    scalar alphagaz_out = vaporfra*sv_gas/(vaporfra*sv_gas+(1.0-vaporfra)* sv_liq); // gas phase volume fraction
    if (alphagaz_out > 1.0  || alphagaz_out == 1.0)
    {
        alphagaz_out = 1.0-vaporfra_virtual;
    }
    else if (alphagaz_out < 0.0 || alphagaz_out == 0.0)
    {
        alphagaz_out = vaporfra_virtual;
    }
    scalar alphaliq_out = 1.0- alphagaz_out;

    //scalar csgas_old=0.0, csliq_old=0.0;
    scalar cv_gas = 0.0, cv_liq = 0.0;   
    scalar cp_gas = 0.0, cp_liq = 0.0;     
    scalar ct_gas = 0.0, ct_liq = 0.0;     
    scalar cs_gas = 0.0, cs_liq = 0.0;     

    scalar amix_gas = 0.0, bmix_gas = 0.0, dadt_gas = 0.0, dadt2 = 0.0;
    scalar amix_liq = 0.0, bmix_liq = 0.0, dadt_liq = 0.0;
    this->abmix(T, comp_gas, amix_gas, bmix_gas, dadt_gas, dadt2); 
    this->abmix(T, comp_liq, amix_liq, bmix_liq, dadt_liq, dadt2); 

    this->cvmix_phase(0, p, T, comp, cv_gas);
    this->cpmix_phase(1, p, T, comp, cp_liq);
    this->soundspeedmix(0, p, T, comp_gas, ct_gas, cs_gas);
    this->soundspeedmix(1, p, T, comp_liq, ct_liq, cs_liq);
    scalar sie_gas = this->siemix_phase(0, p, T, comp_gas);
    scalar sie_liq = this->siemix_phase(1, p, T, comp_liq);
    
    scalar gam_gas = cp_gas/cv_gas;
    scalar gam_liq = cp_liq/cv_liq;

    scalar dpdt_gas = RR*1.0e-03/(sv_gas-bmix_gas)-dadt_gas/(sqr(sv_gas)+2.0*bmix_gas*sv_gas-sqr(bmix_gas));
    scalar dpdt_liq = RR*1.0e-03/(sv_liq-bmix_liq)-dadt_liq/(sqr(sv_liq)+2.0*bmix_liq*sv_liq-sqr(bmix_liq));

    scalar dvdt_p_gas = alphagaz_out*dpdt_gas/(rho_gas*sqr(ct_gas));  
    scalar dvdt_p_liq = alphaliq_out*dpdt_liq/(rho_liq*sqr(ct_liq));

    svphase[0]   = sv_gas,     svphase[1]   = sv_liq;
    rhophase[0]  = rho_gas,    rhophase[1]  = rho_liq;
    cpphase[0]   = cp_gas,     cpphase[1]   = cp_liq;
    cvphase[0]   = cv_gas,     cvphase[1]   = cv_liq;
    ctphase[0]   = ct_gas,     ctphase[1]   = ct_liq;
    csphase[0]   = cs_gas,     csphase[1]   = cs_liq;
    siephase[0]  = sv_gas,     siephase[1]  = sv_liq;
    gamphase[0]  = gam_gas,    gamphase[1]  = gam_liq;
    dpdtphase[0] = dpdt_gas,   dpdtphase[1] = dpdt_liq;
    dvdtphase[0] = dvdt_p_gas, dvdtphase[1] = dvdt_p_liq;
}
// end TPn_output




// TPD
template<class Specie>
inline void Foam::PengRobinsonMix<Specie>::Checkphasesplit 
(
    scalar p,		//input pressure
    scalar T,		//input temperature
    scalar *comp_check,	//input phase compositions
    scalar *comp,	//input overall compositions
    scalar& tpd0,		//output phase state
    scalar& trifile,		//output
    scalar *comp_normal		//output
) const
{
    label iter = 0;
    int itmax = 500;
    scalar tolerance_tpd = 1.0e-05;
    scalar tolerance_new = 1.0e-06;
    scalar eps_ss = 1.0;
    scalar fugcoef_mix0[num_sp] = {0.0};
    scalar gibbsenergy0 = 0.0;
    this->fugcoef_mixture(p, T, comp, fugcoef_mix0, gibbsenergy0); //fugcoef_overall
    //printf("Checkphasesplit_0000 %+5.8f %+5.8f\n",comp[0],comp[1]);  	
    scalar comp_phase0[num_sp] = {0.0}; //YYY
    scalar dz[num_sp] = {0.0};
    comp_normal[num_sp] = {0.0}; //yy
    for (label spid = 0; spid < num_sp; spid++)
    {
	comp_phase0[spid] = comp_check[spid];
        dz[spid] = log(comp[spid])+log(fugcoef_mix0[spid]);
    }

    for (iter = 0; iter < itmax; iter++)
    {
    	scalar ysum = 0.0;
    	for (label spid = 0; spid < num_sp; spid++)
    	{
	    ysum += comp_phase0[spid];
    	}
    	for (label spid = 0; spid < num_sp; spid++)
    	{
	    comp_normal[spid] = comp_phase0[spid]/ysum;
    	}
    	
    	for (label spid = 0; spid < num_sp; spid++)
    	{
	    if (comp_normal[spid] < 0.0)
	    {
                comp_normal[spid] = 1.0e-07;
	    }
	}

	scalar fugcoef_mix1[num_sp] = {0.0};
	scalar detfug_dz[num_sp] = {0.0};
	scalar gibbsenergy1 = 0.0;
        //printf("Checkphasesplit_1000 %+5.8f %+5.8f %+5.8f\n",comp_normal[0],comp_normal[1],ysum);
	this->fugcoef_mixture(p, T, comp_normal, fugcoef_mix1, gibbsenergy1);
	
    	for (label spid = 0; spid < num_sp; spid++)
    	{
	    if (fugcoef_mix1[spid] < 1.0e-15)
	    {
                fugcoef_mix1[spid] = 1.0e-15;
	    }
	}	
        //printf("Checkphasesplit_2222 %+5.2i %+5.8f %+5.8f\n",iter, fugcoef_mix1[0],fugcoef_mix1[1]);

	scalar gstar_y = 0.0;
	scalar f_old[num_sp] = {0.0};
    	for (label spid = 0; spid < num_sp; spid++)
    	{
	    detfug_dz[spid] = fugcoef_mix1[spid]-dz[spid];
	    if (detfug_dz[spid] < 0.0)
	    {
		detfug_dz[spid] = -detfug_dz[spid];
	    }
	    f_old[spid] = sqrt(comp_normal[spid])*(log(comp_normal[spid])+log(detfug_dz[spid]));
	    gstar_y += comp_normal[spid]*(log(comp_normal[spid])+log(fugcoef_mix1[spid])-dz[spid]-1.0);
    	}
        //printf("Checkphasesplit_3333 %+5.2i %+5.8f\n",iter, gstar_y);
	gstar_y = gstar_y +1.0;
	scalar maxcomp = 0.0;
	scalar comp_det[num_sp] = {0.0};
    	for (label spid = 0; spid < num_sp; spid++)
    	{
	    comp_det[spid] = comp_phase0[spid]-comp[spid];
	    if (comp_det[spid] < 0.0)
	    {
		comp_det[spid] = -comp_det[spid];
	    }
	    maxcomp = max(maxcomp, comp_det[spid]);
    	}

	if (maxcomp < tolerance_tpd)
	{
	    trifile = 1.0;
	    tpd0 = 0.0;
	    return;
	}

	scalar ysumnew = 0.0;
	scalar comp_phase[num_sp] = {0.0}; //YYY0
	eps_ss = 0.0;
    	for (label spid = 0; spid < num_sp; spid++)
    	{
	    comp_phase[spid] = comp[spid]*fugcoef_mix0[spid]/fugcoef_mix1[spid]; 
	    ysumnew += comp_phase[spid];
	    eps_ss += sqr(comp_phase[spid] - comp_phase0[spid]);
    	}
        //printf("Checkphasesplit_4444 %+5.2i %+5.8f %+5.8f\n",iter, comp_phase[0],comp_phase[1]);
	if (eps_ss < tolerance_new)
	{
	    trifile = 0.0;
	    if (gstar_y < 0.0)
	    {
     	    	tpd0 = 1.0;
	    }
	    else
	    {
     	    	tpd0 = 0.0;
	    }
	    break;
	}

	scalar fugcoef_mixnew[num_sp] = {0.0};
	scalar comp_normalnew[num_sp] = {0.0}; //yy0
	scalar comp_ref[num_sp] = {0.0};
	scalar g[num_sp] = {0.0};
	scalar s[num_sp] = {0.0};

	scalar sssggg[num_sp] = {0.0};
	scalar gibbsenergynew = 0.0;
	scalar f_new[num_sp] = {0.0};
	if (iter > 20)  //newton method 
	{
	    for (label spid = 0; spid < num_sp; spid++)
	    {
            	comp_normalnew[spid]=comp_phase[spid]/ysumnew; 
	    }
	    this->fugcoef_mixture(p, T, comp_normalnew, fugcoef_mixnew, gibbsenergynew); //fugcoef_yy0 
	    for (label spid = 0; spid < num_sp; spid++)
	    {
		f_new[spid] = sqrt(comp_normalnew[spid])*(log(comp_normalnew[spid])+log(fugcoef_mixnew[spid])-dz[spid]);
	        g[spid] = f_new[spid]-f_old[spid];
	        s[spid] = 2.0*sqrt(comp_normalnew[spid])-2.0*sqrt(comp_normal[spid]);
	    }
	    scalar dotproductsg = 0.0, dotproductgg = 0.0, dotproductsf = 0.0, dotproductsgf = 0.0;
	    for (label spid = 0; spid < num_sp; spid++)
	    {
		dotproductsg += s[spid] * g[spid];
		dotproductgg += g[spid] * g[spid];
		dotproductsf += s[spid] * f_old[spid];
		dotproductsgf += (s[spid]-g[spid]) * f_old[spid];
		sssggg[spid] = dotproductsf*g[spid]/dotproductsg;
	    }
	    scalar sg = 1.0/dotproductsg;
	    scalar ssgg = dotproductsgf + dotproductgg*dotproductsf/dotproductsg;
	    
	    // update comp_phase
	    for (label spid = 0; spid < num_sp; spid++)
	    {
		comp_ref[spid] = 2.0*sqrt(comp_normal[spid]);
		comp_phase[spid] = comp_ref[spid] - (f_old[spid]+sg*ssgg*s[spid]-sssggg[spid]);
		comp_phase[spid] = 0.25*sqr(comp_phase[spid]);
	    }
	}

	for (label spid = 0; spid < num_sp; spid++)
	{
	    comp_phase0[spid] = comp_phase[spid];
	}
	if (iter > 99)
	{
	    scalar sumcompphase0 = 0.0;
	    for (label spid = 0; spid < num_sp; spid++)
	    {
	        sumcompphase0 += comp_phase0[spid];
	    }
	    if (sumcompphase0 > 1.0)
	    {
            	tpd0 = 1.0;
	    }
	    else
	    {
		tpd0 = 0.0;
	    }
	    break;
	}         
    }

//printf("Checkphasesplit_6666 %+5.3i,%+5.8f,%+5.8f,%+5.8f\n",iter,eps_ss,tpd0,comp_normal[0]);  
    trifile = 1.0;
    if (iter >= itmax ) 
    {
	printf("Checkphasesplit() did not converge %5.10f\n",p*1.e-05); 
    }
//FatalErrorInFunction
//<< exit(FatalError); 
}
   




// TPD
template<class Specie>
inline void Foam::PengRobinsonMix<Specie>::TPD 
(
    scalar p,		//input
    scalar T,		//input
    scalar *comp,	//input
    scalar *comp_liq,	//output
    scalar *comp_gas,	//output
    scalar& tpdtest	//output
) const
{
    scalar equalconstant[num_sp] = {0.0};
    scalar tpdliq = 0.0, tpdgas = 0.0, trifileliq = 0.0, trifilegas = 0.0;
    scalar comp_liqout[num_sp] = {0.0};
    scalar comp_gasout[num_sp] = {0.0};

    //this->TPn_flash(p, T, comp, comp_liq, comp_gas, vaporfra, equalconstant);
    this->wilsoneq(p, T, equalconstant);
    for (label spid = 0; spid < num_sp; spid++)
    {
	comp_liq[spid]=comp[spid]*(1.0/equalconstant[spid]);
	comp_gas[spid]=equalconstant[spid]*comp[spid];
    }

    //printf("TPD1111 %5.8f  %5.8f\n",comp[0],comp[1]); 
    this-> Checkphasesplit(p,T,comp_liq,comp,tpdliq,trifileliq,comp_liqout);
    this-> Checkphasesplit(p,T,comp_gas,comp,tpdgas,trifilegas,comp_gasout);

    //check phase state
    scalar sumcompmK = 0.0, sumcompoK = 0.0;
    if (tpdliq == 1.0 || tpdgas == 1.0)
    {
	tpdtest = 2.0; //two phase
    }
    else 
    {
    	for (label spid = 0; spid < num_sp; spid++)
    	{
            sumcompmK += comp[spid]*equalconstant[spid]; //sum(k(i)*z(i))<1, in liquid (Rachford_rice)
            sumcompoK += comp[spid]/equalconstant[spid]; //sum(z(i)/K(i))<1, in gas
    	}
	if (sumcompmK < 1.0 || sumcompoK > 1.0)
	{
	    tpdtest = 1.0; //liquid
	}
    	else if (sumcompmK > 1.0 || sumcompoK <1.0)
    	{
	    tpdtest = 0.0; //gas
    	}
    }
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
template<class Specie>
inline void Foam::PengRobinsonMix<Specie>::kappa_phase
(
    int phasei,
    scalar p, 
    scalar T, 
    scalar *comp_phase,
    scalar& thcond_phase,    
    scalar& viscos
) const
{
    scalar w[num_sp][num_sp] = {{0.0}};
    scalar sigma[num_sp][num_sp] = {{0.0}};
    scalar e_k0[num_sp][num_sp] = {{0.0}};
    scalar k[num_sp] = {0.0}; //need to be sure
    scalar e_k[num_sp] = {0.0};
    scalar kij[num_sp][num_sp] = {{0.0}};
    scalar mwij[num_sp][num_sp] = {{0.0}};
    scalar mw_sp3[num_sp] = {0.0},Vc_sp3[num_sp] = {0.0};

    for (label spid = 0; spid < num_sp; spid++)
    {
	e_k[spid] = Tc_sp[spid]/1.2593;
	mw_sp3[spid] = mw_sp[spid]; //g/mol
	Vc_sp3[spid] = Vc_sp[spid]*1.0e+03; //Vc is m3/kmol, Vc_sp3 is cm3/mol
    }
    for (label spid = 0; spid < num_sp; spid++)
    {
	for (label spjd = 0; spjd < num_sp; spjd++)
	{
	    w[spid][spjd]=0.50*(omega_sp[spid]+omega_sp[spjd]); 
	    sigma[spid][spjd]=sqrt((0.8090*pow(Vc_sp3[spid],1.0/3.0))*(0.8090*pow(Vc_sp3[spjd],1.0/3.0)));
            //sigma(i,j)= bip_kij(i,j)*sqrt((0.809d0*vc(i)**(1.d0/3.d0))*(0.809d0*vc(j)**(1.d0/3.d0)))
	    e_k0[spid][spjd]=sqrt(e_k[spid]*e_k[spjd]);
            //e_k0(i,j)= bip_kij(i,j)*sqrt(e_k(i)*e_k(j))
	    kij[spid][spjd]=sqrt(k[spid]*k[spjd]);
	    mwij[spid][spjd]=2.0*mw_sp3[spid]*mw_sp3[spjd]/(mw_sp3[spid]+mw_sp3[spjd]);
	}
    }

    scalar sigma_mix3=0.0;      
    scalar sigma_mix2=0.0;
    for (label spid = 0; spid < num_sp; spid++)
    {
	for (label spjd = 0; spjd < num_sp; spjd++)
	{
            sigma_mix3 = sigma_mix3 + comp_phase[spid]*comp_phase[spjd]*pow3(sigma[spid][spjd]);
            sigma_mix2 = sigma_mix2 + comp_phase[spid]*comp_phase[spjd]*sqr(sigma[spid][spjd]);
	}
    }
    //scalar sigma_mix_23 = sigma_mix2/sigma_mix3;
    scalar e_k_mix=0.0;
    scalar w_mix=0.0;
    scalar k_mix=0.0;
    scalar mu4=0.0;
    scalar mu[num_sp] = {0.0};
    int num_sp0 = 0, num_spend = num_sp;
		    //co2   ch4	   o2	n2    h2o
    // mu information is in Reid et al., Apendix A
    scalar mu0[5] = {0.0,  0.0,  0.0,  0.0,  1.8};

    if (system_ == 1.0 || system_ == 5.0)
    {
	num_sp0 = 3, num_spend = 4;
    }
    else if (system_ == 6.0)
    {
	num_sp0 = 1, num_spend = 4;
    }
    else if (system_ == 7.0)
    {
	num_sp0 = 0, num_spend = 3;
    }
    else if (system_ == 8.0 || system_ == 10.0)
    {
	num_sp0 = 0, num_spend = 4;
    }
    else if (system_ == 9.0)
    {
	num_sp0 = 0, num_spend = 5;
    }

    for (label spid = num_sp0; spid < num_spend; spid++)
    {
        mu[spid]  = mu0[spid];
    }    
			 //C12H26 C7H16
    scalar mu0fuel[2] = { 0.0,  0.0};

    if (system_ == 1.0 || system_ == 5.0)
    {
	for (label spid = 0; spid < 1; spid++)
	{
	    mu[spid]  = mu0fuel[spid]; 
	}    
    }
    for (label spid = 0; spid < num_sp; spid++)
    {
	for (label spjd = 0; spjd < num_sp; spjd++)
	{
	    e_k_mix= e_k_mix+comp_phase[spid]*comp_phase[spjd]*e_k0[spid][spjd]*pow3(sigma[spid][spjd])/sigma_mix3;
	    w_mix= w_mix+comp_phase[spid]*comp_phase[spjd]*w[spid][spjd]*pow3(sigma[spid][spjd])/sigma_mix3;
	    mu4=mu4+comp_phase[spid]*comp_phase[spjd]*sqr(mu[spid])*sqr(mu[spjd])/pow3(sigma[spid][spjd]);
	    k_mix=k_mix + comp_phase[spid]*comp_phase[spjd]*kij[spid][spjd];
	}
    }

    scalar mw_binary=0.0;
    for (label spid = 0; spid < num_sp; spid++)
    {
	for (label spjd = 0; spjd < num_sp; spjd++)
	{
	    mw_binary=mw_binary+(comp_phase[spid]*comp_phase[spjd]*e_k0[spid][spjd]*sqr(sigma[spid][spjd])*sqrt(mwij[spid][spjd]))/(e_k_mix*sigma_mix2);  
	}
    }
    mu4=mu4*sigma_mix3;
    mw_binary=sqr(mw_binary);
    scalar volu_mix= sigma_mix3/pow3(0.8090);//Eq.(9-5.43);Eq.(9-4.8) bug??
    //scalar volu_mix= pow3(sigma_mix3/0.8090);
    scalar tempc_mix=1.25930*e_k_mix;   
    scalar mur_mix4=pow4(131.30)*mu4/sqr(tempc_mix*volu_mix);//volu_mix in cm3/mol, mu4 in debyes
    scalar t_mix=T/e_k_mix;
    scalar Fcm=1.0-0.2756*w_mix+0.0590350*mur_mix4+k_mix;
    scalar omegav=1.161450*(pow(t_mix,-0.148740)+0.524870*exp(-0.773200*t_mix)+2.161780*exp(-2.437870*t_mix));

    //1. mixture phases viscosity 40.785 ? 26.69 there is error in 9-5.24
    viscos=1e-07*26.69*Fcm*sqrt(mw_binary*T)/(sigma_mix2*omegav);

    scalar mw_phase = this->mwmix(comp_phase);
    scalar Cv_phase = this->cvmix_phase(phasei, p, T, comp_phase); //J/kgK
    scalar alpha=Cv_phase*mw_phase/(RR*1.0e-03)-1.50;
    scalar tempr=T/tempc_mix;

    //compute specific volume, m3/mol
    scalar amix_phase = 0.0, bmix_phase = 0.0, dadt_phase = 0.0, dadt2 = 0.0;
    this->abmix(T, comp_phase, amix_phase, bmix_phase, dadt_phase, dadt2); 
    scalar sv_phase = this->volmmix_phase(phasei, p, T, comp_phase);  //m3/mol

    //2, heat conductivity,W/m.K   
    scalar beta=0.78620-0.71090*w_mix+1.31680*sqr(w_mix);     
    scalar mw_binary0=mw_binary/1.0e+3;
    scalar zeta=2.00+10.50*sqr(tempr);
    scalar psi=1.0+alpha*((0.2150+0.282880*alpha-1.0610*beta+0.266650*zeta)/(0.63660+beta*zeta+1.0610*alpha*beta));

    scalar qq=3.5860*1e-03*sqrt(tempc_mix/mw_binary0)/pow(volu_mix,2.0/3.0);   
    scalar yy=volu_mix/(6.0e+06*sv_phase);      
    scalar G1=(1.0-0.50*yy)/pow3(1.0-yy);

    // Table 10-3      
    scalar a[7] = {2.41660,    -5.0924e-01, 6.61070,    1.4543e+01, 7.9274e-01,-5.86340,    9.1098e+01}; 
    scalar bb[7]= {7.4824e-01, -1.50940,    5.62070,   -8.91390,    8.2019e-01, 1.2801e+01, 1.2811e+02}; 
    scalar c[7] = {-9.1858e-01,-4.9991e+01, 6.4760e+01,-5.63790,   -6.9369e-01, 9.58930,   -5.4217e+01}; 
    scalar d[7] = {1.2712e+02,  6.9983e+01, 2.7039e+01, 7.4344e+01, 6.31730,    6.5529e+01, 5.2381e+02}; 

    scalar B[8] = {0.0};
    for (label ii = 0; ii < 7; ii++)
    {
         B[ii]=a[ii]+bb[ii]*w_mix+c[ii]*mu4+d[ii]*k_mix;
    }
            
    scalar G2=((B[1]/yy)*(1.0-exp(-B[4]*yy))+B[2]*G1*exp(B[5]*yy)+B[3]*G1)/(B[1]*B[4]+B[2]+B[3]);
    //Info<< "G2* * * *"<<G2<<"* * "<< yy <<"* * "<< volu_mix <<"**"<<sv_phase << endl;
    //Info<< "comp_phase* * * *"<<comp_phase[0] <<"* * "<< comp_phase[1] << endl;    
            
    thcond_phase = 31.20*viscos*psi*(1.0/G2+B[6]*yy)/mw_binary0+qq*B[7]*sqr(yy)*sqrt(tempr)*G2;
    //printf("kappa_phase_thcond_phase %5.8f\n\n",thcond_phase);         
}



// ======================================================================//
// binary diffusivity for LOW pressure based on Ping's critical point	 //
// Harstad and Bellan Ref. Ind. Eng. Chem. Res. 2004, 43, 645-654        //
// which is derived from Reid, 2001, Eq. 11-3.1				 //
// Input: feed, p, T            					 //
// Output: Dij under low pressure			 		 //
// Range: 0.1 - 10 bar                                                   //
//									 //
// History :   Creation.	     P. YI		2019.07		 //
// ======================================================================//

template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::Dij_lowP //for gas
(
    scalar p,
    scalar T,
    scalar *comp_inp
) const
{
    // Tcij and Vcmix are interpolated;
    scalar Tcij = this->TcBi(comp_inp);
    scalar Pcij = this->PcBi(comp_inp); 
    scalar Vcij = 1.0e+06*this->volmmix_phase(0, Pcij, Tcij, comp_inp);
    // Calculate s  
    scalar Trij = T/Tcij;
    scalar sumamT = 0.0;
    scalar am[6] = {-0.84211, -0.32643, -0.10053, 0.07747, 0.0127, -0.00995};
    for (label m = 0; m < 6; m++)
    {
         sumamT=am[m]*pow(log(Trij),m);
    }
    scalar sval = 0.0;
    //if (Trij >= 0.2 && Trij <= 10) // default criterion in Bellan
    if (Trij <= 10.0)//original is 10
    {
	sval = exp(sumamT);
    }
    else
    {
    	sval = 0.2304;
	/*scalar sumamT10 = 0.0;
        for (label m = 0; m < 6; m++)
        {
            sumamT10=am[m]*pow(log(10.0),m);
        }
	sval = exp(sumamT10);*/
    }   

    scalar Mij = 1.0/mw_sp[0]+1.0/mw_sp[1];
    //rD is constant parameter of O(1), an empirical adjustment for collistional interaction between i and j
    //for CO2-H2O=0.93, O2-H2O=0.82; 
    scalar rD = 1.5; 
    scalar Vmol = 1.0e+06*this->volmmix_phase(0, p, T, comp_inp);
    scalar Dij_value = 2.81e-05*Vmol*sqrt(Mij*T)*pow(Trij,sval)/(rD*pow(Vcij,2.0/3));
//printf("Vmol sval Tcij Vcij Dij_value, %5.8f, %5.8f, %5.8f\n",Vmol, sval, Dij_value);    
    return Dij_value;
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

template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::Dij_correct //for gas
(
    scalar p,
    scalar T,
    scalar *comp_inp
) const
{
    // Tcij and Vcmix are interpolated;
    scalar Tcij = this->TcBi(comp_inp);
    scalar Pcij = this->PcBi(comp_inp);
    scalar Vcij = 1.0e+06*this->volmmix_phase(0, Pcij, Tcij, comp_inp);
    scalar Vmol = 1.0e+06*this->volmmix_phase(0, p, T, comp_inp);
    // Reduced density
    scalar rho_reduced = Vcij/Vmol;
    // Fit parameters
    scalar c_adjust = 0.58;
    scalar b_adjust = 3.0/2;    
    // Scaling function
    scalar theta_D = c_adjust * pow(rho_reduced,b_adjust);
    return theta_D;
}


template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::Dij_highP //for dense gas
(
    scalar p,
    scalar T,
    scalar *comp_inp
) const
{
    scalar DijL = this->Dij_lowP(p, T, comp_inp);
    scalar WD   = this->Dij_correct(p, T, comp_inp) + 1.0;
    return DijL / WD;
}



// ======================================================================//
// binary diffusivity for LOW pressure based on pseudo's critical point	 //
// Harstad and Bellan Ref. Ind. Eng. Chem. Res. 2004, 43, 645-654        //
// which is derived from Reid, 2001, Eq. 11-3.1				 //
// Input: feed, p, T            					 //
// Output: Dij under low pressure       		 		 //
//									 //
// History :   Creation.	     P. YI		2019.07		 //
// ======================================================================//

template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::Dij_lowP_pesudo //for gas
(
    scalar p,
    scalar T,
    scalar *comp_inp
) const
{
    scalar bip_val = 0.0;
    scalar bipc12n2 = bipc12n2_;   
    scalar bipco2ch4 = bipco2ch4_;   
    scalar bipco2h2o = bipco2h2o_;   
    scalar bipch4o2 = bipch4o2_;  
    scalar bipc7n2 = bipc7n2_;  

	if (system_ == 1.0)
	{
	    bip_val = bipc12n2;
	}
	if (system_ == 2.0)
	{
	    bip_val = bipco2ch4;
	}
	if (system_ == 3.0)
	{
	    bip_val = bipco2h2o;
	}
	if (system_ == 4.0)
	{
	    bip_val = bipch4o2;
	}
	if (system_ == 5.0)
	{
	    bip_val = bipc7n2;
	}
	if (system_ == 12.0)//just used for validation Dij
	{
	    bip_val = 0.0305;
	}
	if (system_ == 13.0)//just used for validation Dij
	{
	    bip_val = 0.103;
	}
	if (system_ == 14.0)//just used for validation Dij
	{
	    bip_val = 0.08;
	}


    // Tcij and Vcmix are interpolated;
    scalar Vc_sp3[num_sp] = {0.0};
    for (label spid = 0; spid < num_sp; spid++)
    {
	Vc_sp3[spid] = Vc_sp[spid]*1.0e+03;//cm3/mol
    }
    scalar sumVij = pow(Vc_sp3[0],1.0/3) + pow(Vc_sp3[1],1.0/3);
    scalar Vcij = pow(sumVij,3.0)/8.0;

    scalar bip_kijB = 1.0-(1.0-bip_val)*pow(Vc_sp3[0]*Vc_sp3[1],1.0/2)/Vcij;
    scalar Tcij = (1.0-bip_kijB)*sqrt(Tc_sp[0]*Tc_sp[1]);

    scalar Trij = T/Tcij;
    scalar sumamT = 0.0;
    scalar am[6] = {-0.84211, -0.32643, -0.10053, 0.07747, 0.0127, -0.00995};
    for (label m = 0; m < 6; m++)
    {
         sumamT=am[m]*pow(log(Trij),m);
    }

    scalar sval = 0.0;
    //if (Trij >= 0.2 && Trij <= 10) // default criterion in Bellan
    if (Trij <= 10)
    {
	sval = exp(sumamT);
    }
    else
    {
	sval = 0.2304;
	/*scalar sumamT10 = 0.0;
        for (label m = 0; m < 6; m++)
        {
            sumamT10=am[m]*pow(log(10.0),m);
        }
	sval = exp(sumamT10);*/
    }

    scalar Mij = 1.0/mw_sp[0]+1.0/mw_sp[1];
    //rD is constant parameter of O(1), an empirical adjustment for collistional interaction between i and j
    //for CO2-H2O=0.93, O2-H2O=0.82; 
    scalar rD = 1.5; 
    scalar Vmol = RR*1.0e+03*T/p;//cm3/mol
    scalar Dij_value = 2.81e-05*Vmol*sqrt(Mij*T)*pow(Trij,sval)/(rD*pow(Vcij,2.0/3));
    printf("Vmol Tcij Vcij pesudo, %5.8f, %5.8f, %5.8f\n",Vmol, Tcij, Vcij);          
    return Dij_value;
}


template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::Dij_correct_pesudo //for gas
(
    scalar p,
    scalar T,
    scalar *comp_inp
) const
{
    // Tcij and Vcmix are interpolated;
    scalar Vc_sp3[num_sp] = {0.0};
    for (label spid = 0; spid < num_sp; spid++)
    {
	Vc_sp3[spid] = Vc_sp[spid]*1.0e+03;//cm3/mol
    }
    scalar sumVij = pow(Vc_sp3[0],1.0/3) + pow(Vc_sp3[1],1.0/3);
    scalar Vcij = pow(sumVij,3.0)/8.0;
    scalar Vmol = RR*1.0e+03*T/p;//cm3/mol
    // Reduced density
    scalar rho_reduced = Vcij/Vmol;
    // Fit parameters
    scalar c_adjust = 0.58;
    scalar b_adjust = 3.0/2;    
    // Scaling function
    scalar theta_D = c_adjust * pow(rho_reduced,b_adjust);
    return theta_D;
}


template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::Dij_highP_pesudo //for dense gas
(
    scalar p,
    scalar T,
    scalar *comp_inp
) const
{
    scalar DijL = this->Dij_lowP_pesudo(p, T, comp_inp);
    scalar WD   = this->Dij_correct_pesudo(p, T, comp_inp) + 1.0;
    
    return DijL / WD;
}




// ======================================================================//
// binary critical Temperature based on PengRobinson			 //
// Using interpolation method proposed by Ping			         //
// Output: critical Temperature				 		 //
//									 //
// History :   Creation.	     P. YI		2019.07		 //
// ======================================================================//

template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::TcBi
(
    scalar *comp_inp
) const
{
    // Tcij and Vcmix are interpolated;
    if (comp_inp[0] < 0.001)
    {
	comp_inp[0] = 0.001;
    }
    else if (comp_inp[0] > 0.999)
    {
	comp_inp[0] = 0.999;
    }

    scalar Tcij_value = 0.0;

    scalar Ac = 0.0, Bc = 0.0, Cc = 0.0, Dc = 0.0, Ec = 0.0;
    if (system_ == 1)
    {
	Ac = -4599.2;
	Bc = 7069.5;
	Cc = -3556.6;
	Dc = 583.23;
	Ec = 639.92;
	Tcij_value = Ac*pow(comp_inp[0],4.0) + Bc*pow(comp_inp[0],3.0) + Cc*pow(comp_inp[0],2.0) + Dc*comp_inp[0] + Ec;
    }
    else if (system_ == 2)
    {
	Ac = 0.0;
	Bc = 2.64300;
	Cc = -42.057;
	Dc = 153.940;
	Ec = 189.330;
	Tcij_value = Ac*pow(comp_inp[0],4.0) + Bc*pow(comp_inp[0],3.0) + Cc*pow(comp_inp[0],2.0) + Dc*comp_inp[0] + Ec;
    }
    else if (system_ == 3)
    {
	Ac = -913.97;
	Bc =  1158.70;
	Cc = -408.49;
	Dc = -181.21;
	Ec =  643.46;
	Tcij_value = Ac*pow(comp_inp[0],4.0) + Bc*pow(comp_inp[0],3.0) + Cc*pow(comp_inp[0],2.0) + Dc*comp_inp[0] + Ec;
    }
    else if (system_ == 12)
    {
	Ac =  0.0;
	Bc = -18.5;
	Cc =  6.536;
	Dc =  75.447;
	Ec =  145.72;
	Tcij_value = 1.0*(Ac*pow(comp_inp[0],4.0) + Bc*pow(comp_inp[0],3.0) + Cc*pow(comp_inp[0],2.0) + Dc*comp_inp[0] + Ec);
    }
    else if (system_ == 13)
    {
	Ac =  0.0;
	Bc = -148.57;
	Cc =  191.49;
	Dc =  -116.12;
	Ec =  128.73;
	Tcij_value = 1.0*(Ac*pow(comp_inp[0],4.0) + Bc*pow(comp_inp[0],3.0) + Cc*pow(comp_inp[0],2.0) + Dc*comp_inp[0] + Ec);
    }
    else if (system_ == 14)
    {
	Ac =  0.0;
	Bc = -199.37;
	Cc =  241.17;
	Dc =  -230.99;
	Ec =  378.66;
	Tcij_value = 1.0*(Ac*pow(comp_inp[0],4.0) + Bc*pow(comp_inp[0],3.0) + Cc*pow(comp_inp[0],2.0) + Dc*comp_inp[0] + Ec);
    }

    return Tcij_value;
}



// ======================================================================//
// binary critical Pressure based on PengRobinson			 //
// Using interpolation method proposed by Ping			         //
// Output: critical Pressure				 		 //
//									 //
// History :   Creation.	     P. YI		2019.07		 //
// ======================================================================//

template<class Specie>
inline Foam::scalar Foam::PengRobinsonMix<Specie>::PcBi
(
    scalar *comp_inp
) const
{
    // Tcij and Vcmix are interpolated;
    if (comp_inp[0] < 0.001)
    {
	comp_inp[0] = 0.001;
    }
    else if (comp_inp[0] > 0.999)
    {
	comp_inp[0] = 0.999;
    }

    scalar Pcij_value = 0.0;
    scalar Ac = 0.0, Bc = 0.0, Cc = 0.0, Dc = 0.0, Ec = 0.0;
    if (system_ == 1)
    {
	if (comp_inp[0] <= 0.8)
	{
	   Ac = 0.0;
	   Bc = 1832.6;
	   Cc = -1398.8;
	   Dc = 366.08;
	   Ec = 14.286;
	}
	else if (comp_inp[0] > 0.8 && comp_inp[0] < 0.97)
	{
	   Ac = 0.0;
	   Bc = 3139.6;
	   Cc = 15764.0;
	   Dc = -28261.0;
	   Ec = 11303.0;
	}
	else if (comp_inp[0] >= 0.97)
	{
	   Ac = 0.0;
	   Bc = 0.0;
	   Cc = 0.0;
	   Dc = -69237.0;
	   Ec = 69588.0;
	}

	Pcij_value = Ac*pow(comp_inp[0],4.0) + Bc*pow(comp_inp[0],3.0) + Cc*pow(comp_inp[0],2.0) + Dc*comp_inp[0] + Ec;
    }
    else if (system_ == 2)
    {
	Ac = 0.0;
	Bc = -40.041;
	Cc = -54.88;
	Dc = 122.82;
	Ec = 42.713;
	Pcij_value = Ac*pow(comp_inp[0],4.0) + Bc*pow(comp_inp[0],3.0) + Cc*pow(comp_inp[0],2.0) + Dc*comp_inp[0] + Ec;
    }
    else if (system_ == 3)
    {
	Ac = 0.0;
	Bc = 102.83;
	Cc = -1778.6;
	Dc = 1519.8;
	Ec = 193.16;
	Pcij_value = Ac*pow(comp_inp[0],4.0) + Bc*pow(comp_inp[0],3.0) + Cc*pow(comp_inp[0],2.0) + Dc*comp_inp[0] + Ec;
    }
    else if (system_ == 12)
    {
	Ac = 0.0;
	Bc = -26.118;
	Cc = -14.933;
	Dc = 50.865;
	Ec = 35.264;
	Pcij_value = 1.0*(Ac*pow(comp_inp[0],4.0) + Bc*pow(comp_inp[0],3.0) + Cc*pow(comp_inp[0],2.0) + Dc*comp_inp[0] + Ec);
    }
    else if (system_ == 13)
    {
	Ac = 0.0;
	Bc = -443.93;
	Cc = 400.32;
	Dc = 54.892;
	Ec = 33.305;
	Pcij_value = 1.0*(Ac*pow(comp_inp[0],4.0) + Bc*pow(comp_inp[0],3.0) + Cc*pow(comp_inp[0],2.0) + Dc*comp_inp[0] + Ec);
    }
    else if (system_ == 14)
    {
	Ac = 0.0;
	Bc = -99.669;
	Cc = -151.52;
	Dc = 203.11;
	Ec = 84.846;
	Pcij_value = 1.0*(Ac*pow(comp_inp[0],4.0) + Bc*pow(comp_inp[0],3.0) + Cc*pow(comp_inp[0],2.0) + Dc*comp_inp[0] + Ec);
    }

    return Pcij_value*1.0e+05;
}








// ======================================================================//
// Mixture Critical Point used in Bellan: 				 //
// Input: feed, species criticals point					 //
// Output: mixture C.P. 				 		 //
//									 //
// Source: L. Sciacovelli and J. Bellan, 2019				 //
//									 //
// History :   Creation.	     P. YI		2019.06		 //
// ======================================================================//
// ** Notice: if use Bellan abmix, the comp_inp is feed rather than each phase compositions
template<class Specie>
inline void Foam::PengRobinsonMix<Specie>::Mixcritic_Bellan 
(
    scalar &Tcmix,
    scalar &Pcmix,
    scalar &Vcmix,
    scalar &Zcmix
) const
{
    //In Bellan's work, single mixture phase is considered
    scalar bipc12n2 = bipc12n2_;   
    scalar bipco2ch4 = bipco2ch4_;   
    scalar bipco2h2o = bipco2h2o_;   
    scalar bipch4o2 = bipch4o2_;  
    scalar bipc7n2 = bipc7n2_;  
    scalar bip_kij[num_sp][num_sp] = { {0.0} };

    if (system_ < 6.0)
    {
	for (label spid = 0; spid < num_sp; spid++)
	{
	    for (label spjd = 0; spjd < num_sp; spjd++)
            {	  
		if (system_ == 1)
		{
		    bip_kij[spid][spjd] = bipc12n2;
		}
		if (system_ == 2)
		{
		    bip_kij[spid][spjd] = bipco2ch4;
		}
		if (system_ == 3)
		{
		    bip_kij[spid][spjd] = bipco2h2o;
		}
		if (system_ == 4)
		{
		    bip_kij[spid][spjd] = bipch4o2;
		}
		if (system_ == 5)
		{
		    bip_kij[spid][spjd] = bipc7n2;
		}

	        if (spid == spjd) 
	        {
		    bip_kij[spid][spjd] = 0.0;
		}
	    }	  	
	}    
    }  
    if (system_ == 6.0)
    {
	string bipkij3_CH4O2N2;
	ifstream input; //infile;
	input.open("bipkij3_CH4O2N2.txt");
	for (label spid = 0; spid < num_sp; spid++)
	{
	    for (label spjd = 0; spjd < num_sp; spjd++)
            {	  
		input >> bip_kij[spid][spjd];
	    }
	}
	input.close();
    }  
    else if (system_ == 7.0)
    {
	//string bipkij3;
	//ifstream input; //infile;
	//input.open("bipkij3.txt");//co2+ch4+o2
	string bipkij3_CO2CH4H2O;
	ifstream input; //infile;	
	input.open("bipkij3_CO2CH4H2O.txt");//co2+ch4+H2O

	for (label spid = 0; spid < num_sp; spid++)
	{
	    for (label spjd = 0; spjd < num_sp; spjd++)
            {	  
		input >> bip_kij[spid][spjd];
	    }
	}
	input.close();
    }  
    else if (system_ == 8.0)
    {
	string bipkij4;
	ifstream input; //infile;
	input.open("bipkij4.txt");
	for (label spid = 0; spid < num_sp; spid++)
	{
	    for (label spjd = 0; spjd < num_sp; spjd++)
            {	  
		input >> bip_kij[spid][spjd];
	    }
	}
	input.close();
    }    
    else if (system_ == 9.0)
    {
	string bipkij5;
	ifstream input; //infile;
	input.open("bipkij5.txt");
	for (label spid = 0; spid < num_sp; spid++)
	{
	    for (label spjd = 0; spjd < num_sp; spjd++)
            {	  
		input >> bip_kij[spid][spjd];
	    }
	}
	input.close();
    }  
    else if (system_ == 10.0)
    {
	string bipkij4_h2o;
	ifstream input; //infile;
	input.open("bipkij4_h2o.txt");
	for (label spid = 0; spid < num_sp; spid++)
	{
	    for (label spjd = 0; spjd < num_sp; spjd++)
            {	  
		input >> bip_kij[spid][spjd];
	    }
	}
	input.close();
    }  
    else if (system_ == 11.0)
    {
	string bipkij_ch;
	ifstream input; //infile;
	input.open("bipkij_ch.txt");
	for (label spid = 0; spid < num_sp; spid++)
	{
	    for (label spjd = 0; spjd < num_sp; spjd++)
            {	  
		input >> bip_kij[spid][spjd];
	    }
	}
	input.close();
    }   
    if (system_ == 12.0)//CH4+n2
    {
	for (label spid = 0; spid < num_sp; spid++)
	{
	    for (label spjd = 0; spjd < num_sp; spjd++)
            {	  
		bip_kij[spid][spjd] = 0.0305;
	        if (spid == spjd) 
	        {
		    bip_kij[spid][spjd] = 0.0;
		}
	    }
	}
    }
    if (system_ == 13.0)//H2+n2
    {
	for (label spid = 0; spid < num_sp; spid++)
	{
	    for (label spjd = 0; spjd < num_sp; spjd++)
            {	  
		bip_kij[spid][spjd] = 0.103;//Bellan2019
	        if (spid == spjd) 
	        {
		    bip_kij[spid][spjd] = 0.0;
		}
	    }
	}
    }
    if (system_ == 14.0)//CH4-H2S
    {
	for (label spid = 0; spid < num_sp; spid++)
	{
	    for (label spjd = 0; spjd < num_sp; spjd++)
            {	  
		bip_kij[spid][spjd] = 0.08;//Bellan2019
	        if (spid == spjd) 
	        {
		    bip_kij[spid][spjd] = 0.0;
		}
	    }
	}
    }
     
    scalar Tc_binary[num_sp][num_sp] = {{0.0}};
    scalar Vc_binary[num_sp][num_sp] = {{0.0}};
    scalar Zc_binary[num_sp][num_sp] = {{0.0}};
    scalar Pc_binary[num_sp][num_sp] = {{0.0}};
    scalar bip_kijB[num_sp][num_sp] = {{0.0}};
    scalar Zc_sp[num_sp] = {0.0};
    for (label spid = 0; spid < num_sp; spid++)
    {
    	for (label spjd = 0; spjd < num_sp; spjd++)
    	{
	    scalar sumvolm3 = pow(Vc_sp[spid],1.0/3)+pow(Vc_sp[spjd],1.0/3);
            Vc_binary[spid][spjd] = (1.0/8)*pow(sumvolm3,3.0);
            bip_kijB[spid][spjd] = 1.0-(1.0-bip_kij[spid][spjd])*pow(Vc_sp[spid]*Vc_sp[spjd],1.0/2)/Vc_binary[spid][spjd];
 	    //printf("Vc_sp %5.3i %5.8f  %5.8f\n", spjd,Vc_binary[spid][spjd], Vc_sp[spjd]);          
    	}
    }    

    for (label spid = 0; spid < num_sp; spid++)
    {
        Zc_sp[spid] = Pc_sp[spid]*Vc_sp[spid]/(RR*1.0e-03*Tc_sp[spid]);
    }

    for (label spid = 0; spid < num_sp; spid++)
    {
    	for (label spjd = 0; spjd < num_sp; spjd++)
    	{
            Tc_binary[spid][spjd] = (1.0-bip_kijB[spid][spjd])*sqrt(Tc_sp[spid]*Tc_sp[spjd]);
            Zc_binary[spid][spjd] = (1.0/2)*(Zc_sp[spid]+Zc_sp[spjd]);
            Pc_binary[spid][spjd] = RR*1.0e-03*Tc_binary[spid][spjd]*Zc_binary[spid][spjd]/Vc_binary[spid][spjd];
 	    //printf("Tc_sp %5.3i %5.8f  %5.8f %5.8f\n", spjd,bip_kijB[spid][spjd], Tc_sp[spjd],Tc_binary[spid][spjd]);          
    	}
    }    

    Tcmix = 0.0, Pcmix = 0.0, Vcmix = 0.0, Zcmix = 0.0;
    if (num_sp == 2)
    {
	Tcmix = Tc_binary[0][1];
	Vcmix = Vc_binary[0][1];
	Pcmix = Pc_binary[0][1];
	Zcmix = Zc_binary[0][1];
 	printf("Tcmix %5.8f %5.8f %5.8f\n", Tcmix,Pcmix, Vcmix);        
    }
FatalErrorInFunction
<< exit(FatalError);  
}







// ======================================================================//
// No phase change: From sie to Temp					 //
// Input: feed, p, T, and sie, ro					 //
// Output: relaxed p and T 				 		 //
//									 //
// History :   Creation.	     P. YI		2019.06		 //
// ======================================================================//

template<class Specie>
inline void Foam::PengRobinsonMix<Specie>::from_sie_ro_to_PT 
(
    scalar sie_inp,
    scalar ro_inp,
    scalar *comp_inp,	//input:feed
    scalar &press,
    scalar &temp
) const
{
    scalar amix = 0.0, bmix = 0.0, dadt = 0.0, dadt2 = 0.0;

    scalar mw_mixture = this->mwmix(comp_inp);
    scalar sv_mixture = mw_mixture/ro_inp;
    int iter = 0;
    int niter_max = 100;
    scalar delta_sie = 1.0e-08;
    while (iter < niter_max)
    {
    	this->abmix(temp, comp_inp, amix, bmix, dadt, dadt2); 
        scalar del = sqr(sv_mixture)+2.0*bmix*sv_mixture-sqr(bmix);
	const scalar cc0 = 2.0*sqrt(2.0);
	const scalar cc1 = 1.0 + sqrt(2.0);
	const scalar cc2 = 1.0 - sqrt(2.0);

	scalar sie_d = ((temp*dadt-amix)/(cc0*bmix))*log((sv_mixture+cc1*bmix)/(sv_mixture+cc2*bmix));

	scalar sieig = this->sieig_phase(temp, comp_inp);
        //call sie_ig(temp0,comp,spdfre(i,1:nspfre),spdfue(i,1:nbffuel),sieig,cv_ig,h_ig,cp_ig)
  	// change the units 
	scalar sie_test=sie_d/mw_mixture + sieig;
	scalar epsilon = sqr((sie_inp-sie_test)/sie_inp);
		     
	if ( epsilon < delta_sie ) 
	{        	  
            press=RR*1.0e-03*temp/(sv_mixture-bmix)- amix/del;
	    break;
	}
	else
	{
	    // update T,P, using NEWton Rapson
            scalar press = RR*1.0e-03*temp/(sv_mixture-bmix) - amix/del;
   	    // compute dudt_v, dudt_v=cv_real
	    // need another subrotine to determin Liq or gas;
            scalar cv_real = this->cvmix_phase(0, press, temp, comp_inp);
	    // Rewton Rapson
	    temp = temp - (sie_inp - sie_test)/(-cv_real); 
	}
	iter ++;
    }
}





    

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::PengRobinsonMix<Specie>::operator+=
(
    const PengRobinsonMix<Specie>& pg
)
{
    scalar Y1 = this->Y();  //YP: Y() is mole number of component in mixture
    Specie::operator+=(pg); //YP: is it mean that sum all Y?  pg means all species Y, so after this, Y() means sumY

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();   // Y1 means mole fraction
        const scalar Y2 = pg.Y()/this->Y(); // pg.Y() means all species corresponding value; Y2 means others mole fraction

        Tc_ = Y1*Tc_ + Y2*pg.Tc_;
        Vc_ = Y1*Vc_ + Y2*pg.Vc_;
        Zc_ = Y1*Zc_ + Y2*pg.Zc_;
        Pc_ = RR*Zc_*Tc_/Vc_;
        omega_ = Y1*omega_ + Y2*pg.omega_;
    }
}


template<class Specie>
inline void Foam::PengRobinsonMix<Specie>::operator*=(const scalar s)
{
     Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


template<class Specie>
Foam::PengRobinsonMix<Specie> Foam::operator+
(
    const PengRobinsonMix<Specie>& pg1,
    const PengRobinsonMix<Specie>& pg2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pg1)
      + static_cast<const Specie&>(pg2)
    );

    if (mag(sp.Y()) < small)
    {
        return PengRobinsonMix<Specie>
        (
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
        const scalar Y1 = pg1.Y()/sp.Y();
        const scalar Y2 = pg2.Y()/sp.Y();

        const scalar Tc = Y1*pg1.Tc_ + Y2*pg2.Tc_;
        const scalar Vc = Y1*pg1.Vc_ + Y2*pg2.Vc_;
        const scalar Zc = Y1*pg1.Zc_ + Y2*pg2.Zc_;

        return PengRobinsonMix<Specie>
        (
            sp,
            Tc,
            Vc,
            Zc,
            RR*Zc*Tc/Vc,
            Y1*pg1.omega_ + Y2*pg2.omega_
        );
    }
}


template<class Specie>
Foam::PengRobinsonMix<Specie> Foam::operator*
(
    const scalar s,
    const PengRobinsonMix<Specie>& pg
)
{
    return PengRobinsonMix<Specie>
    (
        s*static_cast<const Specie&>(pg),
        pg.Tc_,
        pg.Vc_,
        pg.Zc_,
        pg.Pc_,
        pg.omega_
    );
}


template<class Specie>
Foam::PengRobinsonMix<Specie> Foam::operator==
(
    const PengRobinsonMix<Specie>& pg1,
    const PengRobinsonMix<Specie>& pg2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pg1)
     == static_cast<const Specie&>(pg2)
    );

    const scalar Y1 = pg1.Y()/sp.Y();
    const scalar Y2 = pg2.Y()/sp.Y();

    const scalar Tc = Y2*pg2.Tc_ - Y1*pg1.Tc_;
    const scalar Vc = Y2*pg2.Vc_ - Y1*pg1.Vc_;
    const scalar Zc = Y2*pg2.Zc_ - Y1*pg1.Zc_;

    return PengRobinsonMix<Specie>
    (
        sp,
        Tc,
        Vc,
        Zc,
        RR*Zc*Tc/Vc,
        Y2*pg2.omega_ - Y1*pg1.omega_
    );
}


// ************************************************************************* //
