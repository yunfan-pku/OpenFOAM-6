#include "thermotable.H"

template <class ThermoType>
inline Foam::thermotable<ThermoType>::thermotable(
    HYmultiComponentMixture<ThermoType> &multi_mixture_in,
    const word &name,
    const ThermoType &st)
    : 
    ThermoType(name, st),
    table(multi_mixture_in.table),
    multi_mixture_(multi_mixture_in),
    Y_(multi_mixture_in.Y().size(), 0.0)
{
}

template <class ThermoType>
inline void Foam::thermotable<ThermoType>::from_Y_to_z
(
    const scalarList& Yinp,
    scalarList& zi
) const
{
    scalarList Ymw(Yinp.size(),Zero);
    scalar sumYmw = 0.0;
    forAll(Yinp,i)
    {
        Ymw[i] = Yinp[i] / this->multi_mixture_.getLocalThermo(i).W();
        sumYmw += Yinp[i] / this->multi_mixture_.getLocalThermo(i).W();
    }
    forAll(Yinp,i)
    {
        zi[i] = Ymw[i] / sumYmw;
    }
}


// ======================================================================//
// Reverse interpolate tpn_tabulation_reverse.dat               	 //
// Input:                                                                //
//      Internal energy/enthaphy, temperature0, pressure, celli  	 //
// Output: target temperature   			 		 //
//									 //
// History :   Creation.	     P. YI			2019.07  //
// ======================================================================//
template<class ThermoType>
Foam:: scalar Foam::thermotable<ThermoType>::from_tpn_reverse_to_T
(
    const scalar propval,
    const scalar p,
    const scalar T0
)const
{
    scalar eps_tpn = 1.0e-04;
    scalarList zi(Y_.size(),Zero);
    this->from_Y_to_z(Y_, zi);

    scalar temp_guess = T0;//this->clip(temp_guess, ueos_tmin, ueos_tmax);
    scalar guess_value = this->from_tpn_to_property(temp_guess, p, zi[0], 5);
        
    scalar temp_b;
    
    // Table line and column information
    //int rown = table.rownum;
    int rowT = table.dnum[0];//table.indv[1];
    //int rowp = table.dnum[1];//table.indv[2];
    //int rowY = table.dnum[2];//table.indv[3];    
    
    // Find the minimum and maximum temperature in Table
    //point_find posi  = table.find(temp_guess, table.indv[0], 0, table.dnum[0]-1);   
    scalar ueos_tmin = table.indv[0][0];
    scalar ueos_tmax = table.indv[0][table.dnum[0]-1];
        
    scalar ueos_tstep = 0.1*(ueos_tmax-ueos_tmin)/rowT;  
    /*if (celli < 2)
    {
        Info<<"ueos_tmin= "<<ueos_tmin<<", ueos_tmax= "<< ueos_tmax<<", ueos_tstep= "<<ueos_tstep<<", rowT= "<<rowT<<endl;   
        Info<<"guess_value= "<<guess_value<<", temp_guess= "<< temp_guess<<endl;           
    }*/
    
    scalar temp_extra = 0.0; 
    if(guess_value > propval+eps_tpn) //propval is energy input value
    {
        while(guess_value > propval)
        {
            temp_guess -= ueos_tstep;
            if(temp_guess >= ueos_tmax)
            {
                temp_extra = ueos_tmax; //extrapolate_temperature_tpn(temp_guess, pres, &ymas, propval);
                return temp_extra;
            }
            else if(temp_guess <= ueos_tmin)
            {
	        temp_extra = ueos_tmin; //extrapolate_temperature_tpn(temp_guess, pres, &ymas, propval);
                return temp_extra;
            }
            else
            {
                guess_value = this->from_tpn_to_property(temp_guess, p, zi[0], 5);
            }
            /*if (celli < 2)
            {
                Info<<"temp_extra>= "<<temp_extra<<", ueos_tmax="<<ueos_tmax<<", ueos_tmin="<<ueos_tmin<<endl;   
                Info<<"guess_value= "<<guess_value<<", propval= "<<propval<<endl;    
                FatalErrorInFunction
                << exit(FatalError);         
            }*/
        }
    }
    else if(guess_value < propval-eps_tpn)
    {
        while(guess_value < propval)
        {
            temp_guess += ueos_tstep; 

            if(temp_guess >= ueos_tmax)
            {
                temp_extra = ueos_tmax; //extrapolate_temperature_tpn(temp_guess, pres, &ymas, propval);
                return temp_extra;
            }
            else if(temp_guess <= ueos_tmin)
            {
                temp_extra = ueos_tmin; //extrapolate_temperature_tpn(temp_guess, pres, &ymas, propval);
                return temp_extra;
            }
            else
            {
                guess_value = this->from_tpn_to_property(temp_guess, p, zi[0], 5);
            }
            /*if (celli < 2)
            {
                Info<<"temp_extra<= "<<temp_extra<<", ueos_tmax="<<ueos_tmax<<", ueos_tmin="<<ueos_tmin<<endl;   
                Info<<"guess_value= "<<guess_value<<", propval= "<<propval<<endl;  
                FatalErrorInFunction
                << exit(FatalError);                                  
            }*/           
        }
    }
    else
    {
        return temp_guess;
    }

    //posi is a struct; posi is the number index; indv[posi.prev] is the corresponding value;
    point_find posi=table.find(temp_guess, table.indv[0], 0, table.dnum[0]-1); 

    int low = posi.prev; // the lower discreate point correspoing to the input temperature;
    int up  = posi.next; // the uper discreate point correspoing to the input temperature;
    scalar temp_low = table.indv[0][low];
    scalar temp_up = table.indv[0][up];    
    
    scalar low_value  = from_tpn_to_property(temp_low, p, zi[0], 5);//TPn_table_lookup(ueos_temp[low], pres, massfrac);
    scalar high_value = from_tpn_to_property(temp_up, p, zi[0], 5);//TPn_table_lookup(ueos_temp[up], pres, massfrac);
    /*if (celli < 2)
    {
        Info<<"low_value= "<<low_value<<endl;   
        Info<<"high_value= "<<high_value<<endl;           
    }*/
    
    int iT1 = 0;
    if(propval<low_value)
    {
        iT1 = low-1;
        temp_guess = table.indv[0][iT1];
    }
    else if(low_value <= propval && propval <= high_value)
    {
        iT1 = low;
    }
    else
    {
        iT1 = up;
        temp_guess = table.indv[0][iT1+1];
    }
    int iT2 = iT1 + 1;
    
    /*if (celli < 2)
    {
        Info<<"iT1= "<<iT1<<endl;   
    }*/

    if(iT1 <= 0 || iT2 >= rowT + 1)
    {
        //warning_print(916, "3_temperature is out of tabulation bound in reverse lookup. Using the extrapolated value \n");
        temp_extra = this->extrapolate_temperature_tpn(propval,temp_guess, p);
        return temp_extra;
    }
    else
    {
        scalar T1 = table.indv[0][iT1];
        scalar T2 = table.indv[0][iT2];

        low_value =  from_tpn_to_property(T1, p, zi[0], 5);
        high_value = from_tpn_to_property(T2, p, zi[0], 5);

        temp_b = (propval - low_value)*(T2-T1)/(high_value - low_value) + T1;
    }
    /*if (celli < 2)
    {
        Info<<"temp_b= "<<temp_b<<endl;   
    }*/

    return temp_b;
}


// ======================================================================//
// Extrapolate temperature out of range					 //
// Input: guess-temperature, target energy, press			 //
// Output: target temperature   			 		 //
//									 //
// History :   Creation.	     P. YI			2019.07  //
// ======================================================================//

template<class ThermoType>
Foam:: scalar Foam::thermotable<ThermoType>::extrapolate_temperature_tpn
(
    const scalar propval,
    const scalar temp_guess,
    const scalar p
)const
{
    scalarList zi(Y_.size(),Zero);

    this->from_Y_to_z(Y_, zi);

    scalar ueos_tmin = table.indv[0][0];
    scalar ueos_tmax = table.indv[0][table.dnum[0]-1];
    
    int rowT = table.dnum[0];//table.indv[1];
    //int rowp = table.dnum[1];//table.indv[2];
    //int rowY = table.dnum[2];//table.indv[3];  
    
    int it1 = 0, it2 = 0;
    if(temp_guess >= ueos_tmax)
    {
        it1 = rowT-1;
        it2 = rowT;
    }
    else if (temp_guess <= ueos_tmin)
    {
        it1 = 0;//1;
        it2 = 1;//2;
    }
    //point_find posi  = table.find(temp_guess, table.indv[0], 0, table.dnum[0]-1);   
    scalar T1 = table.indv[0][it1];
    scalar T2 = table.indv[0][it2];

    scalar low_value  = from_tpn_to_property(T1, p, zi[0], 5);
    scalar high_value = from_tpn_to_property(T1, p, zi[0], 5);

    scalar temp = (propval - high_value)*(T2 - T1)/(high_value - low_value) + T2;
    return temp;
}



template <class ThermoType>
double Foam::thermotable<ThermoType>::from_tpn_to_property
(
    const double T,
    const double p,
    const double zi,
    int index
) const //index starts from 1
{
    // std::cout<<T<<" "<<p<<" "<<zi<<"\n";
    double input[3] = {T, p, zi};
    if (index == 1)
        return T;
    if (index == 2)
        return p;
    if (index == 3)
        return zi;

    return table.interpolation(input, index - 3);
}

//flag_XYY
template <class ThermoType>
double Foam::thermotable<ThermoType>::get_property_from_tpn_reverse
(
    const double rho,
    const double sie,
    const double zi,
    int index
) const
{
    double t;
    //std::cout<<std::setiosflags(std::ios::fixed)<<std::setprecision(8)<<rho<<" "<<sie<<" "<<zi<<"\n";
    if (index == 1)
    {
        t = table.interpolation_tab(rho, sie, zi, 1);
        // std::cout<<t<<"\n";
        //if(t<100)
        // t=100;
        return t;
    }
    if (index == 2)
        return table.interpolation_tab(rho, sie, zi, 2);
    if (index == 3)
        return zi;
    if (index == 4)
        return rho;
    if (index == 5)
        return sie;
    return from_tpn_to_property(table.interpolation_tab(rho, sie, zi, 1), table.interpolation_tab(rho, sie, zi, 2), zi, index);
}

template<class ThermoType>
Foam:: scalar Foam::thermotable<ThermoType>::mw//g/mol
(
    const scalar p,
    const scalar T
)const
{
    scalarList zi(Y_.size(),Zero);

    this->from_Y_to_z(Y_, zi);
    
    scalar mw_avg = 0.0;
    forAll(Y_, i)
    {
	mw_avg += zi[i] * this->multi_mixture_.getLocalThermo(i).W();
    }

    return mw_avg;
}

template <class ThermoType>
Foam::scalar Foam::thermotable<ThermoType>::mu(
    const scalar p,
    const scalar T) const
{
    scalarList zi(Y_.size(),Zero);
    this->from_Y_to_z(Y_, zi);
    scalar mixture_mu = this->from_tpn_to_property(T, p, zi[0], 10);

    return mixture_mu;
}

template <class ThermoType>
Foam::scalar Foam::thermotable<ThermoType>::rho //kg/m3
    (
        const scalar p,
        const scalar T) const
{
 
    scalarList zi(Y_.size(),Zero);
    this->from_Y_to_z(Y_, zi);
    scalar mixture_rho = this->from_tpn_to_property(T, p, zi[0], 4);
    //  if (celli < 2)
    //  {
    //  Info<< "PY_multiComponentMixture_rho= "<< mixture_rho<<endl;
    //  }

    return mixture_rho;
}

template <class ThermoType>
Foam::scalar Foam::thermotable<ThermoType>::volspecific //m3/mol
    (
        const scalar p,
        const scalar T) const
{
    scalar mw_avg = 1.0e-03 * this->mw(p, T); //kg/mol
    scalar mixture_volspecific = mw_avg / this->rho(p, T);

    return mixture_volspecific;
}
// ======================================================================//
// From energy to T                                                    	 //
// Input:                                                                //
//      Internal energy/enthaphy, temperature0, pressure  	 //
// Output:     			 		                         //
//      target temperature      			 		 //
// Called by : 			 	                        	 //
//      multiComponentMixtureI.H/hePsiThermo.C		        	 //
//									 //
// History :   Creation.	     P. YI			2019.07  //
// ======================================================================//

template <class ThermoType>
Foam::scalar Foam::thermotable<ThermoType>::THE(
    const scalar h,
    const scalar p,
    const scalar T0) const
{
    scalarList zi(Y_.size(),Zero);
    this->from_Y_to_z(Y_, zi);

    
        //double test_rho=this->from_tpn_to_property(T0, p, zi[0], 4);
        //double test_sie=this->from_tpn_to_property(T0, p, zi[0], 5);
        //double test_T=this->get_property_from_tpn_reverse( test_rho, test_sie, zi[0], 1);
        //double test_p=this->get_property_from_tpn_reverse( test_rho, test_sie, zi[0], 2);

    

    //h-hig = (sie-sieig)+pv-RT
    scalar vol = this->volspecific(p, T0);
    scalar rho = this->rho(p, T0);
    scalar sie = h - p * vol; // + 10000.0;
    //scalar cell_T2;
    scalar cell_T1;
    //flag_XYY
    //   Info<<"XYYdebug_T_my "<<rho<<" "<<sie<<" "<<zi[0]<<" "<<cell_T1<<" "<<cell_T2<<" "<<T0<<" "<<cell_T1-cell_T2<<endl;
    cell_T1 = this->get_property_from_tpn_reverse(rho, sie, zi[0], 1);
    //cell_p = this-> get_property_from_tpn_reverse(rho,sie,zi[0],2);
    //if(true||cell_T1==1e308)
    //Info<<"XYYdebug_T "<<sie<<" "<<p<<" "<<T0<<endl;

/*
    if (p < 0)
        cell_T2 = this->from_tpn_reverse_to_T(sie, 1, T0);
    else
        cell_T2 = this->from_tpn_reverse_to_T(sie, p, T0);

    if (cell_T1 - cell_T2 >= 0.1)
        //   Info<<"XYYdebug_T_my "<<rho<<" "<<sie<<" "<<zi[0]<<" "<<cell_T1<<" "<<cell_T2<<" "<<T0<<" "<<cell_T1-cell_T2<<endl;

        //Info<<"cell_T3^^T0= "<<T0<<", T= "<<cell_T<<", Y= "<<Yinp[0]<<", rho= "<<rho<<", sie= "<<sie<<", h= "<<h<<endl;
        //  if (this->cell_ < 2)
        //  {
        //Info<<"cell_T3^^T0= "<<T0<<", T= "<<cell_T<<", Y= "<<Yinp[0]<<", rho= "<<rho<<", sie= "<<sie<<", h= "<<h<<endl;
        // }
        // if(cell_T1==1e308)
        return cell_T2;
      else
      */
       return cell_T1;
}

template <class ThermoType>
Foam::scalar Foam::thermotable<ThermoType>::psi //s2/m2 = rho/p
    (
        const scalar p,
        const scalar T) const
{
    scalar mixture_psi = p * this->volspecific(p, T);
    scalar mw_avg = 1.0e-03 * this->mw(p, T); //kg/mol
    /*
    if (this->cell_ < 2)
    {
        Info << "PY_multiComponentMixture_psi= " << T << ", p=" << p << ", mixture_psi=" << mixture_psi << endl;
    }
    */

    return mw_avg / mixture_psi;
}


template<class ThermoType>
Foam:: scalar Foam::thermotable<ThermoType>::alphah//make sure the unit
(
    const scalar p,
    const scalar T
)const
{
    scalar mixture_alphah = this->kappa(p,T)/this->Cp(p,T); 
    return mixture_alphah;
}

template<class ThermoType>
Foam:: scalar Foam::thermotable<ThermoType>::kappa
(
    const scalar p,
    const scalar T
)const
{
    scalarList zi(Y_.size(),Zero);
    this->from_Y_to_z(Y_, zi);
    scalar mixture_kappa = this->from_tpn_to_property(T,p,zi[0],9);
    return mixture_kappa;
}

template<class ThermoType>
Foam:: scalar Foam::thermotable<ThermoType>::Cp
(
    const scalar p,
    const scalar T
)const
{

    scalarList zi(Y_.size(),Zero);
    this->from_Y_to_z(Y_, zi);
    scalar mixture_cp = this->from_tpn_to_property(T,p,zi[0],6);
    
    return mixture_cp;
}

template<class ThermoType>
Foam:: scalar Foam::thermotable<ThermoType>::Cpv
(
    const scalar p,
    const scalar T
)const
{

    scalarList zi(Y_.size(),Zero);
    this->from_Y_to_z(Y_, zi);
    scalar mixture_cp = this->from_tpn_to_property(T,p,zi[0],6);
    
    return mixture_cp;
}

template<class ThermoType>
Foam:: scalar Foam::thermotable<ThermoType>::HE
(
    const scalar p,
    const scalar T
) const
{
    scalarList zi(Y_.size(),Zero);
    this->from_Y_to_z(Y_, zi);
    
    scalar mixture_vol = this->volspecific(p,T);
    scalar mixture_sie = this->from_tpn_to_property(T,p,zi[0],5);
    scalar mixture_He = mixture_sie + p*mixture_vol;   
/*
    if (celli < 2)
    {
        Info<< "PY_multiComponentMixture___T= " <<T<<", p= "<<p<<", zi="<<zi[0]<<", Yinp="<<Yinp[0]<<endl;      
        Info<< "PY_multiComponentMixture_sie= "<<mixture_sie<<", He="<<mixture_He<< endl;
    }
    */
    return mixture_He;
}

// ======================================================================//
// Why in sensibleEnthalpy or absolute enthalpy, gamma = 1?		 //
// Now I use the normal gamma rather than 1				 //
// History :   Creation.	     P. YI		2019.07		 //
// ======================================================================//
template<class ThermoType>
Foam:: scalar Foam::thermotable<ThermoType>::gamma
(
    const scalar p,
    const scalar T
)const
{
    scalarList zi(Y_.size(),Zero);
    this->from_Y_to_z(Y_, zi);

    scalar mixture_cp = this->from_tpn_to_property(T,p,zi[0],6);
    scalar mixture_cv = this->from_tpn_to_property(T,p,zi[0],7);
    scalar mixture_gamma = mixture_cp/mixture_cv;

    return mixture_gamma;
}
template<class ThermoType>
Foam:: scalar Foam::thermotable<ThermoType>::CpByCpv
(
    const scalar p,
    const scalar T
)const
{
    return this->gamma(p,T);
}