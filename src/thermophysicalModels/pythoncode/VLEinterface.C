#include "VLEinterface.h"


void solver::update()
{

    dict = new Foam::dictionary(Foam::IFstream(path.c_str())());
    int n = specie.size();

    char** slist = new char* [n];
    for (int i = 0;i < n;i++)
    {
        slist[i] = new char[20];
        strcpy(slist[i], specie[i].c_str());
        //std::cout << slist[i] << std::endl;
    }
    s = new Foam::speciesTable(n, (const char**)slist);
    PR = new Foam::PengRobinsonM<Foam::specie>(*s, *dict);
    comp_liq.resize(n);
    comp_gas.resize(n);
    equalconstant.resize(n);
    delete slist;

}
solver::solver(std::string path_i) :path(path_i), dict(nullptr), PR(nullptr) {}
solver::~solver() {
    if (dict)
        delete dict;
    if (s)
        delete s;
    if (PR)
        delete PR;
}
void solver::solve(bool flag)
{
    vaporfra = 0.9;
    Foam::scalarList equalconstant_of(X.size(), Foam::Zero);
    Foam::scalarList comp_of(X.size(), Foam::Zero);
    Foam::scalarList comp_liq_of(X.size(), Foam::Zero);
    Foam::scalarList comp_gas_of(X.size(), Foam::Zero);
    for (int i = 0;i < X.size();i++)
    {
        comp_of[i] = X[i];
        equalconstant_of[i] = equalconstant[i];
    }
    PR->TPn_flash(P, T, comp_of, comp_liq_of, comp_gas_of, vaporfra, equalconstant_of, flag);
    for (int i = 0;i < X.size();i++)
    {
        comp_liq[i] = comp_liq_of[i];
        comp_gas[i] = comp_gas_of[i];
        equalconstant[i] = equalconstant_of[i];
    }

}
double solver::Z()
{
    Foam::scalarList comp_liq_of(X.size(), Foam::Zero);
    Foam::scalarList comp_gas_of(X.size(), Foam::Zero);
    Foam::scalarList comp_of(X.size(), Foam::Zero);
    for (int i = 0;i < X.size();i++)
    {
        comp_liq_of[i] = comp_liq[i];
        comp_gas_of[i] = comp_gas[i];
        comp_of[i] = X[i];
    }
    double alphagas = PR->Evaluate_alpha(P, T, vaporfra, comp_liq_of, comp_gas_of, comp_of);
    double VspecificGas = PR->volmmix_phase(0, P, T, comp_gas_of); //m3/mol
    double VspecificLiq = PR->volmmix_phase(1, P, T, comp_liq_of);
    double ZGas = P * VspecificGas / (RR * 1.0e-03 * T);
    double ZLiq = P * VspecificLiq / (RR * 1.0e-03 * T);
    double ZMixture = ZGas * alphagas + ZLiq * (1.0 - alphagas);
    return  ZMixture;

}
double solver::density()
{
    Foam::scalarList comp_liq_of(X.size(), Foam::Zero);
    Foam::scalarList comp_gas_of(X.size(), Foam::Zero);
    Foam::scalarList comp_of(X.size(), Foam::Zero);
    for (int i = 0;i < X.size();i++)
    {
        comp_liq_of[i] = comp_liq[i];
        comp_gas_of[i] = comp_gas[i];
        comp_of[i] = X[i];
    }

    double mw_gas = PR->mwmix(comp_gas_of); //kg/mol
    double mw_liq = PR->mwmix(comp_liq_of); //kg/mol
    double mw_mixture = mw_gas * vaporfra + mw_liq * (1.0 - vaporfra);
    double rhoMixture = P * mw_mixture / (Z() * RR * 1.0e-03 * T);
    return rhoMixture;

}
bool solver::twophase(double rho, double lt, double rt)
{
    double tm = (lt + rt) / 2;

    Foam::scalarList comp_liq_of(X.size(), Foam::Zero);
    Foam::scalarList comp_gas_of(X.size(), Foam::Zero);
    Foam::scalarList comp_of(X.size(), Foam::Zero);
    Foam::scalarList equalconstant_of(X.size(), Foam::Zero);
    for (int i = 0;i < X.size();i++)
    {
        comp_of[i] = X[i];
    }
    while (rt - lt > 1e-4)
    {
        tm = (lt + rt) / 2;
        PR->TPn_flash(P, tm, comp_of, comp_liq_of, comp_gas_of, vaporfra, equalconstant_of, true);

        double rho_gas = PR->rhomix_phase(0, P, tm, comp_gas_of); //kg/m3
        double rho_liq = PR->rhomix_phase(1, P, tm, comp_liq_of);
        if (rho_gas < rho && rho < rho_liq && rho_liq - rho_gas>1e-3)
            return true;
        if (rho_gas < rho && rho > rho_liq)
            rt = tm;
        else
        {
            lt = tm;
        }

    }
    T=tm;
    return false;

}



