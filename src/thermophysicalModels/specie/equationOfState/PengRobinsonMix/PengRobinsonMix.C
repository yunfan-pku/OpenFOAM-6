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

#include "PengRobinsonMix.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::PengRobinsonMix<Specie>::PengRobinsonMix
(
    const dictionary& dict
)
:
    Specie(dict),
    Tc_(readScalar(dict.subDict("equationOfState").lookup("Tc"))),
    Vc_(readScalar(dict.subDict("equationOfState").lookup("Vc"))),
    Zc_(1.0),
    Pc_(readScalar(dict.subDict("equationOfState").lookup("Pc"))),
    omega_(readScalar(dict.subDict("equationOfState").lookup("omega"))),
    mwsp_(readScalar(dict.subDict("equationOfState").lookup("mwsp"))),

    Tcc12_(readScalar(dict.subDict("equationOfState").lookup("Tcc12"))),
    Vcc12_(readScalar(dict.subDict("equationOfState").lookup("Vcc12"))),
    Pcc12_(readScalar(dict.subDict("equationOfState").lookup("Pcc12"))),
    omegac12_(readScalar(dict.subDict("equationOfState").lookup("omegac12"))),
    mwspc12_(readScalar(dict.subDict("equationOfState").lookup("mwspc12"))),

    Tcc7_(readScalar(dict.subDict("equationOfState").lookup("Tcc7"))),
    Vcc7_(readScalar(dict.subDict("equationOfState").lookup("Vcc7"))),
    Pcc7_(readScalar(dict.subDict("equationOfState").lookup("Pcc7"))),
    omegac7_(readScalar(dict.subDict("equationOfState").lookup("omegac7"))),
    mwspc7_(readScalar(dict.subDict("equationOfState").lookup("mwspc7"))),

    Tcco2_(readScalar(dict.subDict("equationOfState").lookup("Tcco2"))),
    Vcco2_(readScalar(dict.subDict("equationOfState").lookup("Vcco2"))),
    Pcco2_(readScalar(dict.subDict("equationOfState").lookup("Pcco2"))),
    omegaco2_(readScalar(dict.subDict("equationOfState").lookup("omegaco2"))),
    mwspco2_(readScalar(dict.subDict("equationOfState").lookup("mwspco2"))),

    Tcch4_(readScalar(dict.subDict("equationOfState").lookup("Tcch4"))),
    Vcch4_(readScalar(dict.subDict("equationOfState").lookup("Vcch4"))),
    Pcch4_(readScalar(dict.subDict("equationOfState").lookup("Pcch4"))),
    omegach4_(readScalar(dict.subDict("equationOfState").lookup("omegach4"))),
    mwspch4_(readScalar(dict.subDict("equationOfState").lookup("mwspch4"))),

    Tco2_(readScalar(dict.subDict("equationOfState").lookup("Tco2"))),
    Vco2_(readScalar(dict.subDict("equationOfState").lookup("Vco2"))),
    Pco2_(readScalar(dict.subDict("equationOfState").lookup("Pco2"))),
    omegao2_(readScalar(dict.subDict("equationOfState").lookup("omegao2"))),
    mwspo2_(readScalar(dict.subDict("equationOfState").lookup("mwspo2"))),

    Tch2o_(readScalar(dict.subDict("equationOfState").lookup("Tch2o"))),
    Vch2o_(readScalar(dict.subDict("equationOfState").lookup("Vch2o"))),
    Pch2o_(readScalar(dict.subDict("equationOfState").lookup("Pch2o"))),
    omegah2o_(readScalar(dict.subDict("equationOfState").lookup("omegah2o"))),
    mwsph2o_(readScalar(dict.subDict("equationOfState").lookup("mwsph2o"))),

    TcH2_(readScalar(dict.subDict("equationOfState").lookup("TcH2"))),
    VcH2_(readScalar(dict.subDict("equationOfState").lookup("VcH2"))),
    PcH2_(readScalar(dict.subDict("equationOfState").lookup("PcH2"))),
    omegaH2_(readScalar(dict.subDict("equationOfState").lookup("omegaH2"))),
    mwspH2_(readScalar(dict.subDict("equationOfState").lookup("mwspH2"))),

    TcH2S_(readScalar(dict.subDict("equationOfState").lookup("TcH2S"))),
    VcH2S_(readScalar(dict.subDict("equationOfState").lookup("VcH2S"))),
    PcH2S_(readScalar(dict.subDict("equationOfState").lookup("PcH2S"))),
    omegaH2S_(readScalar(dict.subDict("equationOfState").lookup("omegaH2S"))),
    mwspH2S_(readScalar(dict.subDict("equationOfState").lookup("mwspH2S"))),
    
    system_(readScalar(dict.subDict("equationOfState").lookup("system"))),
    num_sp(readScalar(dict.subDict("equationOfState").lookup("num_sp"))),

    bipc12n2_(readScalar(dict.subDict("equationOfState").lookup("bipc12n2"))),
    bipco2ch4_(readScalar(dict.subDict("equationOfState").lookup("bipco2ch4"))),
    bipco2h2o_(readScalar(dict.subDict("equationOfState").lookup("bipco2h2o"))),
    bipch4o2_(readScalar(dict.subDict("equationOfState").lookup("bipch4o2"))),
    bipc7n2_(readScalar(dict.subDict("equationOfState").lookup("bipc7n2"))),

    temp0(readScalar(dict.subDict("equationOfState").lookup("temp0"))),
    press0(readScalar(dict.subDict("equationOfState").lookup("press0"))),

    det_com1(readScalar(dict.subDict("equationOfState").lookup("det_com1"))),
    det_pres(readScalar(dict.subDict("equationOfState").lookup("det_pres"))),
    det_temp(readScalar(dict.subDict("equationOfState").lookup("det_temp"))),
    num_comp(readScalar(dict.subDict("equationOfState").lookup("num_comp"))),
    num_pres(readScalar(dict.subDict("equationOfState").lookup("num_pres"))),
    num_temp(readScalar(dict.subDict("equationOfState").lookup("num_temp"))),

    compinp1(readScalar(dict.subDict("equationOfState").lookup("compinp1"))),
    compinp2(readScalar(dict.subDict("equationOfState").lookup("compinp2"))),
    compinp3(readScalar(dict.subDict("equationOfState").lookup("compinp3"))),
    compinp4(readScalar(dict.subDict("equationOfState").lookup("compinp4"))),
    compinp5(readScalar(dict.subDict("equationOfState").lookup("compinp5"))),
    
    x_fuel0(readScalar(dict.subDict("equationOfState").lookup("x_fuel0"))),
    
    TPnFlag(readScalar(dict.subDict("equationOfState").lookup("TPnFlag"))),    
    TabFlag(readScalar(dict.subDict("equationOfState").lookup("TabFlag"))),    
    HnFlag(readScalar(dict.subDict("equationOfState").lookup("HnFlag"))),    
    TPDFlag(readScalar(dict.subDict("equationOfState").lookup("TPDFlag"))),              
    
    Neg_rootFlag(readScalar(dict.subDict("equationOfState").lookup("Neg_rootFlag"))),    
    PCFlag(readScalar(dict.subDict("equationOfState").lookup("PCFlag"))),
    flag_tp(readScalar(dict.subDict("equationOfState").lookup("flag_tp"))),
    flag_px(readScalar(dict.subDict("equationOfState").lookup("flag_px")))
{
    evap.p_PengRobinsonMix = this;
    //FatalErrorInFunction
    //<< "pppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp000000"
    //<< exit(FatalError);   
    Zc_ = Pc_*Vc_/(RR*Tc_); 
    if (system_ == 1.0) //ch12+n2
    { 
	Pc_sp[0] = Pcc12_, Pc_sp[1] = Pc_; 
	Vc_sp[0] = Vcc12_, Vc_sp[1] = Vc_; 
	Tc_sp[0] = Tcc12_, Tc_sp[1] = Tc_; 
	mw_sp[0] = mwspc12_, mw_sp[1] = mwsp_;  
	omega_sp[0] = omegac12_, omega_sp[1] = omega_;  
    } 
    else if (system_ == 2.0) //co2+ch4
    {
	Pc_sp[0] = Pcco2_, Pc_sp[1] = Pcch4_; 
	Vc_sp[0] = Vcco2_, Vc_sp[1] = Vcch4_; 
	Tc_sp[0] = Tcco2_, Tc_sp[1] = Tcch4_; 
	mw_sp[0] = mwspco2_, mw_sp[1] = mwspch4_; 
	omega_sp[0] = omegaco2_, omega_sp[1] = omegach4_;  
    }
    else if (system_ == 3.0) //co2+h2o
    {
	Pc_sp[0] = Pcco2_, Pc_sp[1] = Pch2o_; 
	Vc_sp[0] = Vcco2_, Vc_sp[1] = Vch2o_; 
	Tc_sp[0] = Tcco2_, Tc_sp[1] = Tch2o_; 
	mw_sp[0] = mwspco2_, mw_sp[1] = mwsph2o_;
	omega_sp[0] = omegaco2_, omega_sp[1] = omegah2o_;
    }
    else if (system_ == 4.0) //ch4+o2
    {
	Pc_sp[0] = Pcch4_, Pc_sp[1] = Pco2_; 
	Vc_sp[0] = Vcch4_, Vc_sp[1] = Vco2_; 
	Tc_sp[0] = Tcch4_, Tc_sp[1] = Tco2_; 
	mw_sp[0] = mwspch4_, mw_sp[1] = mwspo2_; 
	omega_sp[0] = omegach4_, omega_sp[1] = omegao2_; 
    }
    else if (system_ == 5.0) //c7+n2
    {
	Pc_sp[0] = Pcc7_, Pc_sp[1] = Pc_; 
	Vc_sp[0] = Vcc7_, Vc_sp[1] = Vc_; 
	Tc_sp[0] = Tcc7_, Tc_sp[1] = Tc_; 
	mw_sp[0] = mwspc7_, mw_sp[1] = mwsp_;  
	omega_sp[0] = omegac7_, omega_sp[1] = omega_; 
    }
    else if (system_ == 6.0) //ch4+o2+n2
    {
	Pc_sp[0] = Pcch4_, Pc_sp[1] = Pco2_, Pc_sp[2] = Pc_; 
	Vc_sp[0] = Vcch4_, Vc_sp[1] = Vco2_, Vc_sp[2] = Vc_; 
	Tc_sp[0] = Tcch4_, Tc_sp[1] = Tco2_, Tc_sp[2] = Tc_;  
	mw_sp[0] = mwspch4_, mw_sp[1] = mwspo2_, mw_sp[2] = mwsp_; 
	omega_sp[0] = omegach4_, omega_sp[1] = omegao2_, omega_sp[2] = omega_;   
    }

    else if (system_ == 7.0)  //co2+ch4+o2 change to //co2+ch4+H2O: 07/28/2019
    {
	Pc_sp[0] = Pcco2_, Pc_sp[1] = Pcch4_, Pc_sp[2] = Pch2o_; //Pco2_; 
	Vc_sp[0] = Vcco2_, Vc_sp[1] = Vcch4_, Vc_sp[2] = Vch2o_; //Vco2_; 
	Tc_sp[0] = Tcco2_, Tc_sp[1] = Tcch4_, Tc_sp[2] = Tch2o_; //Tco2_;  
	mw_sp[0] = mwspco2_, mw_sp[1] = mwspch4_, mw_sp[2] = mwsph2o_; //mwspo2_; 
	omega_sp[0] = omegaco2_, omega_sp[1] = omegach4_, omega_sp[2] = omegah2o_;//omegao2_;  
    }
    else if (system_ == 8.0)  //co2+ch4+o2+n2
    {
	Pc_sp[0] = Pcco2_, Pc_sp[1] = Pcch4_, Pc_sp[2] = Pco2_, Pc_sp[3] = Pc_; 
	Vc_sp[0] = Vcco2_, Vc_sp[1] = Vcch4_, Vc_sp[2] = Vco2_, Vc_sp[3] = Vc_; 
	Tc_sp[0] = Tcco2_, Tc_sp[1] = Tcch4_, Tc_sp[2] = Tco2_, Tc_sp[3] = Tc_;  
	mw_sp[0] = mwspco2_, mw_sp[1] = mwspch4_, mw_sp[2] = mwspo2_, mw_sp[3] = mwsp_; 
	omega_sp[0] = omegaco2_, omega_sp[1] = omegach4_, omega_sp[2] = omegao2_, omega_sp[3] = omega_;  
    }
    else if (system_ == 9.0) //co2+ch4+o2+n2+ho2
    {
	Pc_sp[0] = Pcco2_, Pc_sp[1] = Pcch4_, Pc_sp[2] = Pco2_, Pc_sp[3] = Pc_, Pc_sp[4] = Pch2o_; 
	Vc_sp[0] = Vcco2_, Vc_sp[1] = Vcch4_, Vc_sp[2] = Vco2_, Vc_sp[3] = Vc_, Vc_sp[4] = Vch2o_; 
	Tc_sp[0] = Tcco2_, Tc_sp[1] = Tcch4_, Tc_sp[2] = Tco2_, Tc_sp[3] = Tc_, Tc_sp[4] = Tch2o_;  
	mw_sp[0] = mwspco2_, mw_sp[1] = mwspch4_, mw_sp[2] = mwspo2_, mw_sp[3] = mwsp_, mw_sp[4] = mwsph2o_; 
	omega_sp[0] = omegaco2_, omega_sp[1] = omegach4_, omega_sp[2] = omegao2_, omega_sp[3] = omega_, omega_sp[4] = omegah2o_;   
    }
    else if (system_ == 10.0) //co2+ch4+o2+ho2
    {
	Pc_sp[0] = Pcco2_, Pc_sp[1] = Pcch4_, Pc_sp[2] = Pco2_, Pc_sp[3] = Pch2o_; 
	Vc_sp[0] = Vcco2_, Vc_sp[1] = Vcch4_, Vc_sp[2] = Vco2_, Vc_sp[3] = Vch2o_; 
	Tc_sp[0] = Tcco2_, Tc_sp[1] = Tcch4_, Tc_sp[2] = Tco2_, Tc_sp[3] = Tch2o_;  
	mw_sp[0] = mwspco2_, mw_sp[1] = mwspch4_, mw_sp[2] = mwspo2_, mw_sp[3] = mwsph2o_; 
	omega_sp[0] = omegaco2_, omega_sp[1] = omegach4_, omega_sp[2] = omegao2_, omega_sp[3] = omegah2o_;   
    }
    else if (system_ == 11.0)   //c12h26+c7h16
    {
	Pc_sp[0] = Pcc12_, Pc_sp[1] = Pcc7_; 
	Vc_sp[0] = Vcc12_, Vc_sp[1] = Vcc7_; 
	Tc_sp[0] = Tcc12_, Tc_sp[1] = Tcc7_;  
	mw_sp[0] = mwspc12_, mw_sp[1] = mwspc7_; 
	omega_sp[0] = omegac12_, omega_sp[1] = omegac7_;   
    }
    else if (system_ == 12.0)   //ch4+n2
    {
	Pc_sp[0] = Pcch4_, Pc_sp[1] = Pc_; 
	Vc_sp[0] = Vcch4_, Vc_sp[1] = Vc_; 
	Tc_sp[0] = Tcch4_, Tc_sp[1] = Tc_; 
	mw_sp[0] = mwspch4_, mw_sp[1] = mwsp_; 
	omega_sp[0] = omegach4_, omega_sp[1] = omega_; 
    }
    else if (system_ == 13.0)   //H2+n2
    {
	Pc_sp[0] = PcH2_, Pc_sp[1] = Pc_; 
	Vc_sp[0] = VcH2_, Vc_sp[1] = Vc_; 
	Tc_sp[0] = TcH2_, Tc_sp[1] = Tc_; 
	mw_sp[0] = mwspH2_, mw_sp[1] = mwsp_; 
	omega_sp[0] = omegaH2_, omega_sp[1] = omega_; 
    }
    else if (system_ == 14.0)   //CH4+H2S
    {
	Pc_sp[0] = Pcch4_, Pc_sp[1] = PcH2S_; 
	Vc_sp[0] = Vcch4_, Vc_sp[1] = VcH2S_; 
	Tc_sp[0] = Tcch4_, Tc_sp[1] = TcH2S_; 
	mw_sp[0] = mwspch4_, mw_sp[1] = mwspH2S_; 
	omega_sp[0] = omegach4_, omega_sp[1] = omegaH2S_; 
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Specie>
void Foam::PengRobinsonMix<Specie>::write(Ostream& os) const
{
    Specie::write(os);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const PengRobinsonMix<Specie>& pg
)
{
    pg.write(os);
    return os;
}


// ************************************************************************* //
