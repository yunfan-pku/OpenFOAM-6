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


struct solver
{
    std::string path;
    double P;
    double T;
    std::vector<double> X;
    std::vector<std::string> specie;
    Foam::speciesTable* s;
    Foam::dictionary *dict;
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