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

Class
    Foam::PengRobinsonMix

Description
    PengRobinsonMix mixture equation of state.

SourceFiles
    PengRobinsonMixI.H
    PengRobinsonMix.C

\*---------------------------------------------------------------------------*/

#ifndef PengRobinsonMix_H
#define PengRobinsonMix_H
#include "evaporation.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class PengRobinsonMix;

template<class Specie>
inline PengRobinsonMix<Specie> operator+
(
    const PengRobinsonMix<Specie>&,
    const PengRobinsonMix<Specie>&
);

template<class Specie>
inline PengRobinsonMix<Specie> operator*
(
    const scalar,
    const PengRobinsonMix<Specie>&
);

template<class Specie>
inline PengRobinsonMix<Specie> operator==
(
    const PengRobinsonMix<Specie>&,
    const PengRobinsonMix<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const PengRobinsonMix<Specie>&
);



/*---------------------------------------------------------------------------*\
                           Class PengRobinsonMix Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class PengRobinsonMix
:
    public Specie
{
    // Private data //PY: must follow the reading order in .C file
    public:
        evaporation<PengRobinsonMix<Specie>> evap;
    private:    
        //- Critical Temperature [K]
        scalar Tc_;
        //- Critical volume [m^3/kmol]
        scalar Vc_;
        //- Critical compression factor [-]
        scalar Zc_;
        //- Critical Pressure [Pa]
        scalar Pc_;
        //- Acentric factor [-]
        scalar omega_;

        //- Molecular weight [-]
        scalar mwsp_;


        //- Critical Temperature [K]
        scalar Tcc12_;
        //- Critical volume [m^3/kmol]
        scalar Vcc12_;
        //- Critical Pressure [Pa]
        scalar Pcc12_;
        //- Acentric factor [-]
        scalar omegac12_;
        //- Molecular weight [-]
        scalar mwspc12_;

        //- Critical Temperature [K]
        scalar Tcc7_;
        //- Critical volume [m^3/kmol]
        scalar Vcc7_;
        //- Critical Pressure [Pa]
        scalar Pcc7_;
        //- Acentric factor [-]
        scalar omegac7_;
        //- Molecular weight [-]
        scalar mwspc7_;

	//- Critical Temperature [K]
        scalar Tcco2_;
	//- Critical volume [m^3/kmol]
        scalar Vcco2_;
        //- Critical Pressure [Pa]
        scalar Pcco2_;
        //- Acentric factor [-]
        scalar omegaco2_;
        //- Molecular weight [-]
        scalar mwspco2_;


	//- Critical Temperature [K]
        scalar Tcch4_;
	//- Critical volume [m^3/kmol]
        scalar Vcch4_;
        //- Critical Pressure [Pa]
        scalar Pcch4_;
        //- Acentric factor [-]
        scalar omegach4_;
        //- Molecular weight [-]
        scalar mwspch4_;

	//- Critical Temperature [K]
        scalar Tco2_;
	//- Critical volume [m^3/kmol]
        scalar Vco2_;
        //- Critical Pressure [Pa]
        scalar Pco2_;
        //- Acentric factor [-]
        scalar omegao2_;
        //- Molecular weight [-]
        scalar mwspo2_;

	//- Critical Temperature [K]
        scalar Tch2o_;
	//- Critical volume [m^3/kmol]
        scalar Vch2o_;
        //- Critical Pressure [Pa]
        scalar Pch2o_;
        //- Acentric factor [-]
        scalar omegah2o_;
        //- Molecular weight [-]
        scalar mwsph2o_;

	//- Critical Temperature [K]
        scalar TcH2_;
	//- Critical volume [m^3/kmol]
        scalar VcH2_;
        //- Critical Pressure [Pa]
        scalar PcH2_;
        //- Acentric factor [-]
        scalar omegaH2_;
        //- Molecular weight [-]
        scalar mwspH2_;

	//- Critical Temperature [K]
        scalar TcH2S_;
	//- Critical volume [m^3/kmol]
        scalar VcH2S_;
        //- Critical Pressure [Pa]
        scalar PcH2S_;
        //- Acentric factor [-]
        scalar omegaH2S_;
        //- Molecular weight [-]
        scalar mwspH2S_;


        //- System index [-]  
        scalar system_;
	//- Number of speices
	int num_sp;

        //- binary coefficient [-]  
        scalar bipc12n2_;   
        scalar bipco2ch4_;   
    	scalar bipco2h2o_;  
    	scalar bipch4o2_;  
    	scalar bipc7n2_;


	//- Critical pressure [Pa]  
        scalar Pc_sp[10];
	//- Critical temperature [K]  
        scalar Tc_sp[10];
	//- Critical volume [m^3/kmol]  
        scalar Vc_sp[10];
	//- Acentric factor [-]  
        scalar omega_sp[10];
	//- Acentric factor [-]  
        scalar mw_sp[10];

        // Initial temperature for thermal table
	scalar temp0;
        // Initial pressure for thermal table
	scalar press0;

        // Variation of compsition for thermal table
	scalar det_com1;
        // Variation of pressure for thermal table
	scalar det_pres;
        // Variation of temperature for thermal table
	scalar det_temp;
        // Number of composition for thermal table
	int num_comp;
        // Number of pressure for thermal table
	int num_pres;
        // Number of temperature for thermal table
	int num_temp;


        // Initialization of composition for thermal table
	scalar compinp1;
	scalar compinp2;
	scalar compinp3;
	scalar compinp4;
	scalar compinp5;
	scalar x_fuel0;

	// Flag of negative root for cubic EOS
	scalar TPnFlag;

	// Flag of negative root for cubic EOS
	scalar TabFlag;

	// Flag of negative root for cubic EOS
	scalar HnFlag;

	// Flag of negative root for cubic EOS
	scalar TPDFlag;
	
	// Flag of negative root for cubic EOS
	scalar Neg_rootFlag;
        // Flag of phase change
	scalar PCFlag;
        // Flag of TP diagram and PX diagram
	scalar flag_tp;
	// flag_px = 1: px diagram//
	// flag_px = 2: Tx diagram
	scalar flag_px;


public:

    // Constructors

        //- Construct from components
        inline PengRobinsonMix
        (
            const Specie& sp,
            const scalar& Tc,
            const scalar& Vc,
            const scalar& Zc,
            const scalar& Pc,
            const scalar& omega,
            const scalar& mwsp,


            const scalar& Tcc12,
            const scalar& Vcc12,
            const scalar& Pcc12,
            const scalar& omegac12,
            const scalar& mwspc12,

            const scalar& Tcc7,
            const scalar& Vcc7,
            const scalar& Pcc7,
            const scalar& omegac7,
            const scalar& mwspc7,

            const scalar& Tcco2,
            const scalar& Vcco2,
            const scalar& Pcco2,
            const scalar& omegaco2,
            const scalar& mwspco2,

            const scalar& Tcch4,
            const scalar& Vcch4,
            const scalar& Pcch4,
            const scalar& omegach4,
            const scalar& mwspch4,

            const scalar& Tco2,
            const scalar& Vco2,
            const scalar& Pco2,
            const scalar& omegao2,
            const scalar& mwspo2,

            const scalar& Tch2o,
            const scalar& Vch2o,
            const scalar& Pch2o,
            const scalar& omegah2o,
            const scalar& mwsph2o,

            const scalar& TcH2,
            const scalar& VcH2,
            const scalar& PcH2,
            const scalar& omegaH2,
            const scalar& mwspH2,

            const scalar& TcH2S,
            const scalar& VcH2S,
            const scalar& PcH2S,
            const scalar& omegaH2S,
            const scalar& mwspH2S,

            const scalar& system,
         

    	    const scalar& bipc12n2,   
	    const scalar& bipco2ch4,   
	    const scalar& bipco2h2o,   
	    const scalar& bipch4o2,  
	    const scalar& bipc7n2 
        );

        //- Construct from dictionary
        PengRobinsonMix(const dictionary& dict);

        //- Construct as named copy
        inline PengRobinsonMix(const word& name, const PengRobinsonMix&);

        //- Construct and return a clone
        inline autoPtr<PengRobinsonMix> clone() const;

        // Selector from dictionary
        inline static autoPtr<PengRobinsonMix> New
        (
            const dictionary& dict
        );


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "PengRobinsonMix<" + word(Specie::typeName_()) + '>';
        }

        // Fundamental properties

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return enthalpy departure [J/kg]
            inline scalar H(const scalar p, const scalar T) const;

            //- Return Cp departure [J/(kg K)]
            inline scalar Cp(scalar p, scalar T) const;

            //- Return entropy [J/(kg K)]
            inline scalar S(const scalar p, const scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;


            //- Return mass fraction of 0 species, temporary
            inline void initialization(scalar p, scalar T) const;

	    //- Return mass fraction of 0 species, temporary
            inline void initPTXdiagram(scalar press, scalar temp, scalar *comp, scalar *comp_liq, scalar *comp_gas) const;

	    //- Generate Thermo table
            inline void genTabulation(scalar press, scalar temp, scalar *comp) const;
            
	    //- Generate Diffusion table
            inline void genDiffusion(scalar press, scalar temp, scalar *comp) const;

            

	    //- Generate Hx flash
            inline void Hn_flash(scalar press, scalar temp, scalar *comp) const;

	    //- Generate Miture Critical Point
            inline void initMix_critical(scalar *comp) const;


	    //- Return mass fraction of 0 species, temporary
            inline Foam::scalar Evaluate_alpha(scalar press, scalar temp, scalar vaporfra, scalar *comp_liq, scalar *comp_gas,scalar *comp) const;

	    //- Return phase compositions and phase state
            inline void initTPD(scalar press, scalar temp, scalar *comp, scalar *comp_liq, scalar *comp_gas, scalar &tpdtest) const;

            //- Return phase mixture ideal internal energy [J/kg]
            inline Foam::scalar sieig_phase(scalar T, scalar *comp_phase) const;

            //- Return phase mixture ideal entropy [J/(kg K)]
            inline Foam::scalar Sig_phase(scalar T, scalar *comp) const;

            //- Return phase mixture ideal enthalpy [J/kg]
            inline Foam::scalar Hig_phase(scalar T, scalar *comp_phase) const;

            //- Return phase mixture ideal enthalpy [J/kg]
            inline Foam::scalar Hig2_phase(scalar T, scalar *comp_phase) const;

            //- Return phase mixture ideal enthalpy [J/kg]
            inline void cvcpig_phase(scalar T, scalar *comp_phase, scalar& cpideal, scalar& cvideal) const;

            //- Return composition []
            inline void calccomp(scalar *comp) const;
            
	    //- Return molecular weight [kg/mol]
            inline scalar mwmix(scalar *comp_inp) const;

	    //- Return phase compostions [kg/mol]
            inline scalar yTox_phase(scalar yGasInp, scalar *comp_gas, scalar *comp_liq) const;

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rhomix_phase(int phasei, scalar p, scalar T, scalar *comp_inp) const;

            //- Return specific volume m^3/mol]
            inline scalar volmmix_phase(int phasei, scalar p, scalar T, scalar *comp_inp) const;

            //- Return speices enthalpy departure [J/kg]
            inline scalar H_phase(int phasei, const scalar p, const scalar T, scalar *comp) const;

            //- Return speices entropy [J/(kg K)]
            inline scalar S_phase(int phasei, const scalar p, const scalar T, scalar *comp) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi_phase(int phasei, scalar p, scalar T, scalar *comp) const;

            //- Return (Cp - Cv) [J/(kg K]
            inline scalar CpMCv(scalar p, scalar T) const;

            //- Return ABZmix: Amixz, Bmixz, Amix2, Bz []
            inline void ABZmix(scalar p, scalar T, scalar *comp_inp, scalar& Amixz, scalar& Bmixz, scalar *Amix2, scalar *Bz) const;

            //- Return abmix: amix, bmix, dadt, dadt2
            inline void abmix(scalar T, scalar *comp_inp, scalar& amix, scalar& bmix, scalar &dadt, scalar &dadt2) const;

            //- Return compression factor []
            inline scalar Zmix(int phasei, scalar Amixz, scalar Bmixz) const;

	    //- Return compression factor root []
            inline void Root_Z(int phasei, scalar Amixz, scalar Bmixz, scalar *rootz) const;
            
            //- Return species compressibility coefficient in each phase []
            inline void fugcoef_phase(int phasei, scalar p, scalar T, scalar *comp_inp, scalar *fugcoef) const;

            //- Return species compressibility coefficient in mixture []
            inline void fugcoef_mixture(scalar p, scalar T, scalar *comp_inp, scalar *fugcoef_mix, scalar& gibbsenergy) const;

            //- Return phase gibbs energy [J/kg]
            inline void gibbs_phase(int phasei, scalar p, scalar T, scalar *comp_inp, scalar &gibbsenergy) const;
            
            //- Return equl_constant [J/kg]
            inline void wilsoneq(scalar p, scalar T, scalar *equl_constant) const;

            //- Return density of each species [kg/m3]
            inline Foam::scalar rho_spd(int phasei, label spid, scalar p, scalar T) const;

            //- Return density of each species [kg/m3]
            inline Foam::scalar vol_spd(int phasei, label spid, scalar p, scalar T) const;

            //- Return specific heat capibility for mixture [J/kgK]
            inline Foam::scalar cvmix_phase(int phasei, scalar p, scalar T, scalar *comp_phase) const;
            inline Foam::scalar cpmix_phase(int phasei, scalar p, scalar T, scalar *comp_phase) const;

            //- Return internal energy for mixture [kg/m3]
            inline Foam::scalar siemix_phase(int phasei, scalar p, scalar T, scalar *comp_phase) const;

            //- Return enthalpy for mixture [J/kg]
            inline Foam::scalar Hmix(int phasei, scalar p, scalar T, scalar *comp_phase) const;

            //- Return enthalpy for mixture [J/kg]
            inline Foam::scalar Hmix2(int phasei, scalar p, scalar T, scalar *comp_phase) const;

	    //- Return pure component EOS parameters
            inline void ab_spd(scalar p,scalar T, scalar Tcsp, scalar Pcsp, scalar omegasp, scalar &alphaz, scalar &az, scalar &bz, scalar &dadtz) const;

	    //- Return enthalpy for spd mixture [J/kg]
            inline Foam::scalar Hspd_mix(scalar p, scalar T, scalar *comp_inp, int spflag) const;

            //- Return soundspeed for mixture [kg/m3]
            inline void soundspeedmix(int phasei, scalar p, scalar T, scalar *comp, scalar& ct_phase, scalar& cs_phase) const;

            //- Return OAD217_function []
            inline Foam::scalar OAD217_function(scalar XX, scalar YY) const;
            
            //- Return OAD013_function []
            inline void OAD013_function(scalar *comp_in, scalar *comp_out, int& ierr) const;

	    //- Return OAD012_function []
            inline void OAD012_function(scalar p, scalar T, scalar *comp_inp, scalar *equalconstant, scalar *comp_liq, scalar *comp_gas, scalar& vaporf, int& ierr) const;

	    //- Return OAD012_function []
            inline void OAD012_singlephase(scalar p, scalar T, scalar *comp_inp, scalar *equalconstant, scalar tpdvalue, scalar &vaporfra, scalar *comp_liq, scalar *comp_gas) const;

            //- acceleration
            inline void GDEM_SS(scalar *vt, scalar *vt1, scalar *v2, scalar *vout) const;

            //- TPn_flash
            inline void TPn_flash(scalar p, scalar T, scalar *comp, scalar *comp_liq, scalar *comp_gas, scalar& vaporfra, scalar *equalconstant) const;
 
            //- TPn_flash
            inline void TPn_simple(scalar p, scalar T, scalar *comp, scalar *comp_liq, scalar *comp_gas, scalar& vaporfra, scalar *equalconstant) const;

           
	    //- TPn_output
            inline void TPn_output(scalar p, scalar T, scalar *comp, scalar *comp_liq, scalar *comp_gas, scalar *svphase, scalar *rhophase, scalar *zphase, scalar *cpphase, scalar *cvphase, scalar *ctphase, scalar *csphase, scalar *siephase, scalar *gamphase, scalar *dpdtphase, scalar *dvdtphase) const;

	    //- Checkphasesplit
	    inline void Checkphasesplit(scalar p, scalar T, scalar *comp_check, scalar *comp, scalar& tpd0, scalar& trifile, scalar *comp_normal) const;

	    //- TPD
	    inline void TPD(scalar p, scalar T, scalar *comp, scalar *comp_liq, scalar *comp_gas, scalar &tpdtest) const;

	    //- Thermal conductivity
            inline void kappa_phase(int phasei, scalar p, scalar T, scalar *comp_phase, scalar& thcond_phase, scalar& viscos) const;

	    //- Mass diffusivity under low pressure, based on Ping's C.P.
            inline scalar Dij_lowP(scalar p, scalar T, scalar *comp_inp) const;

	    //- Mass diffusivity corrections for high pressure, based on Ping's C.P.
            inline scalar Dij_correct(scalar p, scalar T, scalar *comp_inp) const;

	    //- Mass diffusivity under high pressure, based on Ping's C.P.
            inline scalar Dij_highP(scalar p, scalar T, scalar *comp_inp) const;

	    //- Mass diffusivity under low pressure, based on Bellan's binary crtical point
            inline scalar Dij_lowP_pesudo(scalar p, scalar T, scalar *comp_inp) const;

	    //- Mass diffusivity corrections for high pressure, based on Bellan's C.P.
            inline scalar Dij_correct_pesudo(scalar p, scalar T, scalar *comp_inp) const;

	    //- Mass diffusivity under high pressure, based on Bellan's binary crtical point
            inline scalar Dij_highP_pesudo(scalar p, scalar T, scalar *comp_inp) const;

	    //- Binary critical Temperature
            inline scalar TcBi(scalar *comp_inp) const;

	    //- Binary critical Pressure
            inline scalar PcBi(scalar *comp_inp) const;

	    //- Thermal conductivity
            inline void Mixcritic_Bellan(scalar &Tcmix, scalar &Pcmix, scalar &Vcmix, scalar &Zcmix) const;

	    //- Thermal conductivity
            inline void from_sie_ro_to_PT(scalar sie_inp, scalar ro_inp, scalar *comp_inp, scalar &pres, scalar &temp) const;

    // IO

            //- Write to Ostream
            void write(Ostream& os) const;

    // Member operators

        inline void operator+=(const PengRobinsonMix&);
        inline void operator*=(const scalar);

 
    // Friend operators

        friend PengRobinsonMix operator+ <Specie>
        (
            const PengRobinsonMix&,
            const PengRobinsonMix&
        );

        friend PengRobinsonMix operator* <Specie>
        (
            const scalar s,
            const PengRobinsonMix&
        );

        friend PengRobinsonMix operator== <Specie>
        (
            const PengRobinsonMix&,
            const PengRobinsonMix&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const PengRobinsonMix&
        );
};
// * * * * * * * * * * * * * * * * Global data  * * * * * * * * * * * * * * //

// num_spe for every element in the periodic table
//extern PengRobinsonMix num_spe;
//extern PengRobinsonMix pr_spe, Tr_spe;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PengRobinsonMixI.H"

#ifdef NoRepository
    #include "PengRobinsonMix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
