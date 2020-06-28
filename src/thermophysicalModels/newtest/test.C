#include "dictionary.H"
#include "IFstream.H"
#include "specie.H"
#include "multispecie.H"
#include "perfectGas.H"
#include "hConstThermo.H"
#include "speciesTable.H"
#include "sensibleEnthalpy.H"
#include "thermo.H"
#include "constTransport.H"
#include "PengRobinson.H"
#include "PengRobinsonMixture.H"
#include "VLE.H"
#include "janafThermo.H"
#include "chungTransport.H"
#include "chungTransportMixture.H"
#include "HashPtrTable.H"
#include "multithermo.H"
using namespace Foam;
int main()
{
    dictionary dict(IFstream("system/thermotableDict")());
    speciesTable species;
    wordList s(dict.lookup("species"));
    species.transfer(s);
    Info << species << endl;

    typedef chungTransport<species::thermo<janafThermo<PengRobinson<specie>>, sensibleEnthalpy>> Stype;
    dictionary thermoDict(IFstream("system/thermo")());
    HashPtrTable<Stype> speciesThermo(thermoDict);
    typedef species::multithermo<VLE<chungTransportMixture<PengRobinsonMixture<multispecie<Stype>>>>, sensibleEnthalpy> Mtype;

    PtrList<Stype> speciesData(species.size());
    forAll(species, i)
    {
        speciesData.set(
            i,
            new Stype(*speciesThermo[species[i]]));
    }

    dictionary thermoDictM(IFstream("system/thermoMixture")());
    Mtype thermo("test", speciesData, species, thermoDictM);

    scalarList X(2);
    X[0] = 0.7;
    X[1] = 0.3;
    thermo.setX(X);
    //Info << thermo.rho(23000000, 500) << endl;
    // Info<<thermo.Cp(23000000,500)<<endl;
    // Info<<thermo.Hs(23000000,500)<<endl;
    /*
    for(scalar i=400;i<=550;i+=10)
    Info<<thermo.Hs(23000000, i)<<endl;
    */
    scalar dp = 10;
    //Info << "SS=" << ::sqrt(dp / (thermo.rho(5000000 + dp, 500) - thermo.rho(5000000, 500))) << endl;
    //Info << thermo.Hs(23000000, 500) << endl;
    //Info << thermo.Hs(22000000, 500) << endl;
    //Info << thermo.Hs(21000000, 500) << endl;
    //Info << " psi=" << thermo.psi(1000000, 500) << endl;
    //Info << " psi=" << thermo.psi(100000, 520) << endl;
    
    //double hs = thermo.Hs(23000000, 500);
    //Info << hs << endl;
    //double nt = thermo.THE(-32739.9,100000, 576.628);
    //Info << nt << endl;
    //double hn = thermo.Ha(23000000, nt);
    //Info << hn << endl;
    /*
    double te = 505.1362181469;
    Info << te << "," << thermo.Hs(1.50416e+07, te) - 128817 << endl;
    Info << "---------------------------------------" << endl;
    //te=505.1362181469;
    //Info<<te<<","<<thermo.Hs(1.50416e+07, te)-128817<<endl;
    te = 505.1362181470;
    Info << te << "," << thermo.Hs(1.50416e+07, te) - 128817 << endl;
*/
    //Info << thermo.rho(15149200, 567.096) << endl; //15149200	567.096	129.742
/*
    for(scalar i=341.06;i<=341.07;i+=0.0001)
    Info<<i<<","<<thermo.Hs(341.064, i)+32739.9<<endl;
*/
    //Info<<110000<<","<<180<<","<<(thermo.TPn_flash(110000,180))().vaporfra<<endl;
    ofstream fout("de.csv");
    for(double p=10000;p<10000000;p+=100000)
    for(double T=100;T<600;T+=10)
    fout<<p<<","<<T<<","<<(thermo.TPn_flash(p,T))().vaporfra<<std::endl;
    //Info<<"!_!_!_!_!_!_!"<<(thermo.TPn_flash(100000,341.066))().vaporfra<<endl;
    
    //Info<<thermo.Hs(100000,341.064)<<endl;
    //Info<<thermo.Hs(100000,341.066)<<endl;
    //Info<<thermo.Hs(100000,500.066)<<endl;
    //Info<<thermo.Hs(100000,240.066)<<endl;

    //Info<<speciesThermo<<endl;

    //mixture_("mixture", speciesData_,specieNames,thermoDict)
    /*
    specie water_t(1, 18);
    PengRobinson<specie> water(water_t, 647.096, 0.0559, 22.06e+06, 0.344);
    Info << water.rho(100000, 1000) << endl;
    dictionary dict(IFstream("system/thermo")());
    multispecie<PengRobinson<specie>> mst(dict);
    scalarList X(2);
    X[0]=0.5;
    X[1]=0.5;
    mst.setX(X);
    for (int i = 0; i < mst.N_; i++)
        Info << mst.species()[i] << " pc" << mst[i].Pc_ << " Tc" << mst[i].Tc_ << endl;
    Info << mst[0].rho(100000, 1000) << endl;
    Info << mst[1].rho(10, 2000) << endl;
    Info<<mst.X()<<endl;
    Info<<mst.Y()<<endl;
    PengRobinsonMixture<multispecie<janafThermo<PengRobinson<specie>>>> pms(dict);
    pms.setX(X);
    Info<<"Ha= "<<pms.Ha(10000000,1000)<<endl;
    Info<<pms.rho(100000, 1000)<<endl;
    Info<<pms.fugacityCoefficient(100000, 1000)()<<endl;
    typedef VLE<chungTransportMixture<PengRobinsonMixture<multispecie<chungTransport<janafThermo<PengRobinson<specie>>>>>>> vle;
    species::thermo<vle,sensibleEnthalpy> vlet(dict);
    //sensibleEnthalpy<thermo<vle, sensibleEnthalpy>> vlet(dict);
    X[0]=0.5;
    X[1]=0.5;
    vlet.setX(X);
   // autoPtr<vle::solution> sol(vlet.TPn_flash(10000000,450)); 
   // Info<<sol().vaporfra<<endl;
    Info<<vlet.rho(1000000,450)<<endl;
    Info<<vlet.Ha(1000000,450)<<endl;
    Info<<vlet.psi(1000000,450)<<endl;
    Info<<vlet.rho(1000000,450)<<endl;
    Info<<vlet.Ea(100000,450)<<endl;
    
    dictionary dict(IFstream("system/thermo")());
    janafThermo<perfectGas<specie>> h2o(dict.subDict("H2O"));
    Info<<h2o.Ha(10000,2000)<<" "<<h2o.R()<<endl;
    */

    //dictionary dict(IFstream("system/thermo")());
    //chungTransport<species::thermo<janafThermo<PengRobinson<specie>>,sensibleEnthalpy>> water(dict.subDict("H2O"));
    //Info<<water<<endl;
}