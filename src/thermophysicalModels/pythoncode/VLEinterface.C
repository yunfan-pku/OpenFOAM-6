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

