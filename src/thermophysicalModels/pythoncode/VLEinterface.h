#include "dictionary.H"
#include "IFstream.H"
#include "specie.H"
#include "perfectGas.H"
#include "hConstThermo.H"
#include "speciesTable.H"
#include "sensibleEnthalpy.H"
#include "thermo.H"
#include "constTransport.H"
#include "PengRobinsonS.H"
#include "PengRobinsonM.H"


#include <vector>
#include <string>
#include <iostream>
#include <cmath>


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
#include "IOmanip.H"

struct solver
{
    std::string path;
    double P;
    double T;
    std::vector<double> X;
    std::vector<std::string> specie;
    Foam::speciesTable* s;
    Foam::dictionary* dict;
    void update();
    Foam::PengRobinsonM<Foam::specie>* PR;
    solver(std::string);
    ~solver();
    void solve(bool flag);
    double vaporfra;
    double density();
    double Z();
    bool twophase(double rho, double lt, double rt);
    std::vector<double> equalconstant;
    std::vector<double> comp_liq;
    std::vector<double> comp_gas;
};


struct solver_new
{
    typedef chungTransport<species::thermo<janafThermo<PengRobinson<specie>>, sensibleEnthalpy>> Stype;
    typedef species::multithermo<VLE<chungTransportMixture<PengRobinsonMixture<multispecie<Stype>>>>, sensibleEnthalpy> Mtype;
    std::string path;

    dictionary dict;
    speciesTable species;
    dictionary thermoDict;
    dictionary thermoDictM;
    Mtype* thermo;
    PtrList<Stype> speciesData;
    double P;
    double T;
    //std::vector<double> X;
    std::vector<std::string> m_specie;
    //Foam::speciesTable* s;
    //Foam::dictionary* dict;
    //void update();
    //Foam::PengRobinsonM<Foam::specie>* PR;
    solver_new(std::string);
    void reset();
    ~solver_new();
    //void solve(bool flag);
    double vaporfra;
    //double density();
    double rho();
    double drhodT();
    double drhodP();
    double drhodXi(int);
    double dZdT();
    double dZdXi(int);
    double drhodP_HP();
    double drhoPdH_HP();
    double drhoPdP_HP();
    double drhoPdXi_HP(int);
    double drhoPdH_HsP();
    double drhoPdP_HsP();
    double drhoPdXi_HsP(int);
    void drhoPdXHP_HsP();
    double dTdXi_HsP(int);
    double dvfdH_HP();
    double dvfdP_HP();
    double dvfdXi_HP(int);

    double dvfdH_HsP();
    double dvfdP_HsP();
    double dvfdXi_HsP(int);
    double drhodXi_HP(int);
    double Z();
    double Gibbs_single();
    double G_departure_Mole();
    double Gideal_Mole();
    double G_Mole();
    bool solveTPD_BFGS();
    double A_single();
    double z_single();
    double muideal_Mole(int);
    double Gideal();
    double Z(int);
    double TPn_flash();
    double TPn_flash_New();
    double TPn_flash_New_TPD();
    double dZdT(int);
    double A();
    double dAdT();
    double B();
    double dBdT();
    double Ha();
    double Hs();
    double Ha_singlePhase(int, std::vector<double>&);
    double dHadT_singlePhase(int, std::vector<double>&);
    double Hideal(std::vector<double>&);
    double dHidealdT(std::vector<double>&);
    double Cp();
    double dHadP();
    double dHadXi(int);
    double dHsdXi(int);
    double W();
    double W(std::vector<double>&);
    double T_HsP(double h, double p, double T0);
    void Ln_fugacityCoefficient();
    void Ln_fugacityCoefficient(int);
    void ddT_Ln_fugacityCoefficient(int);
    void ddxi_Ln_fugacityCoefficient(int, int);
    void dvidT();
    void dvidP();
    void dvidXi(int);
    void  fugacityCoefficient(int, std::vector<double>&);

    double dTdP_HP();
    double dTdH_HP();
    double dTdXi_HP(int);
    //bool twophase(double rho, double lt, double rt);
    std::vector<double> equalconstant;
    std::vector<double> comp_liq;
    std::vector<double> comp_gas;
    std::vector<double> comp;
    std::vector<double> ret;
};





int fun_my(int n);
