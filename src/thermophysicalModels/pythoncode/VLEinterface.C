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
    T = tm;
    return false;

}



solver_new::solver_new(std::string path_i) :path(path_i), dict(IFstream(path + "system/thermotableDict")()), thermoDict(IFstream(path + "system/thermo")()), thermoDictM(IFstream(path + "system/thermoMixture")()), thermo(nullptr)
{
    wordList s(dict.lookup("species"));
    species.transfer(s);
    Info << species << endl;
    HashPtrTable<Stype> speciesThermo(thermoDict);
    speciesData.resize(species.size());
    forAll(species, i)
    {
        speciesData.set(
            i,
            new Stype(*speciesThermo[species[i]]));
    }
    thermo = new Mtype("test", speciesData, species, thermoDictM);
}

void solver_new::reset()
{
    if (thermo)
        delete thermo;
    wordList s(dict.lookup("species"));
    s.resize(m_specie.size());
    for (int i = 0;i < m_specie.size();i++)
    {
        s[i] = m_specie[i];
    }
    species.transfer(s);
    Info << species << endl;
    HashPtrTable<Stype> speciesThermo(thermoDict);
    speciesData.resize(species.size());
    forAll(species, i)
    {
        speciesData.set(
            i,
            new Stype(*speciesThermo[species[i]]));
    }
    thermo = new Mtype("test", speciesData, species, thermoDictM);
}
/*
void solver_new::update()
{

}
*/


solver_new::~solver_new()
{
    if (thermo)
        delete thermo;
}


double solver_new::Z(int flag)
{
    //Foam::scalarList comp_liq_of(X.size(), Foam::Zero);
    //Foam::scalarList comp_gas_of(X.size(), Foam::Zero);
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        if (flag == 1)
            comp_of[i] = comp_gas[i];
        else if (flag == 0)
            comp_of[i] = comp_liq[i];
    }

    return ((PengRobinsonMixture<multispecie<Stype>>*)thermo)->Z(P, T, flag, &comp_of);
}

double solver_new::dZdT(int flag)
{
    //Foam::scalarList comp_liq_of(X.size(), Foam::Zero);
    //Foam::scalarList comp_gas_of(X.size(), Foam::Zero);
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }

    return ((PengRobinsonMixture<multispecie<Stype>>*)thermo)->dZdT(P, T, flag, &comp_of);
}

double solver_new::TPn_flash()
{
    Foam::scalarList comp_liq_of(comp.size(), Foam::Zero);
    Foam::scalarList comp_gas_of(comp.size(), Foam::Zero);
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution ret(thermo->Mtype::TPn_flash(P, T)());
    comp_liq.resize(comp.size());
    comp_gas.resize(comp.size());
    equalconstant.resize(comp.size());
    for (int i = 0;i < comp.size();i++)
    {
        comp_liq[i] = ret.X_liq()[i];
        comp_gas[i] = ret.X_gas()[i];
        equalconstant[i] = ret.equalconstant()[i];
    }
    vaporfra = ret.vaporfra;
}

double solver_new::TPn_flash_New_TPD()
{
    Foam::scalarList comp_liq_of(comp.size(), Foam::Zero);
    Foam::scalarList comp_gas_of(comp.size(), Foam::Zero);
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution ret(thermo->Mtype::TPn_flash_New_TPD(P, T)());
    comp_liq.resize(comp.size());
    comp_gas.resize(comp.size());
    equalconstant.resize(comp.size());
    for (int i = 0;i < comp.size();i++)
    {
        comp_liq[i] = ret.X_liq()[i];
        comp_gas[i] = ret.X_gas()[i];
        equalconstant[i] = ret.equalconstant()[i];
    }
    vaporfra = ret.vaporfra;
}

double solver_new::TPn_flash_New()
{
    Foam::scalarList comp_liq_of(comp.size(), Foam::Zero);
    Foam::scalarList comp_gas_of(comp.size(), Foam::Zero);
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution ret(thermo->Mtype::TPn_flash_New(P, T)());
    comp_liq.resize(comp.size());
    comp_gas.resize(comp.size());
    equalconstant.resize(comp.size());
    for (int i = 0;i < comp.size();i++)
    {
        comp_liq[i] = ret.X_liq()[i];
        comp_gas[i] = ret.X_gas()[i];
        equalconstant[i] = ret.equalconstant()[i];
    }
    vaporfra = ret.vaporfra;
}

double solver_new::A()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }

    return ((PengRobinsonMixture<multispecie<Stype>>*)thermo)->A(P, T, &comp_of);
}
double solver_new::dAdT()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }

    return ((PengRobinsonMixture<multispecie<Stype>>*)thermo)->dAdT(P, T, &comp_of);
}

double solver_new::B()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }

    return ((PengRobinsonMixture<multispecie<Stype>>*)thermo)->B(P, T, &comp_of);
}
double solver_new::dBdT()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }

    return ((PengRobinsonMixture<multispecie<Stype>>*)thermo)->dBdT(P, T, &comp_of);
}
double solver_new::W()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    return thermo->W(&comp_of);
}
double solver_new::W(std::vector<double>& in)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = in[i];
    }
    return thermo->W(&comp_of);
}

void solver_new::fugacityCoefficient(int flag, std::vector<double>& in)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = in[i];
    }
    Foam::autoPtr<scalarList> pret(thermo->PengRobinsonMixture<multispecie<Stype>>::fugacityCoefficient(P, T, flag, &comp_of));
    ret.resize(comp.size());
    for (int i = 0;i < comp.size();i++)
    {
        ret[i] = pret()[i];
    }
}
double solver_new::Ha()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);

    return thermo->Ha(P, T);
}

double solver_new::Hs()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    return thermo->Hs(P, T);
}
double solver_new::Ha_singlePhase(int flag, std::vector<double>& in)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = in[i];
    }
    return thermo->PengRobinsonMixture<multispecie<Stype>>::Ha(P, T, flag, &(comp_of));
}
double solver_new::dHadT_singlePhase(int flag, std::vector<double>& in)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = in[i];
    }
    return thermo->PengRobinsonMixture<multispecie<Stype>>::dHadT(P, T, flag, &(comp_of));
}

double solver_new::Hideal(std::vector<double>& in)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = in[i];
    }
    return thermo->PengRobinsonMixture<multispecie<Stype>>::Hideal(P, T, &(comp_of));
}
double solver_new::dHidealdT(std::vector<double>& in)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = in[i];
    }
    return thermo->PengRobinsonMixture<multispecie<Stype>>::dHidealdT(P, T, &(comp_of));
}

double solver_new::Cp()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);

    return thermo->Cp(P, T);
}
double solver_new::dHadP()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution ret(thermo->Mtype::TPn_flash(P, T)());
    return thermo->dHadP(P, T, ret);
}

double solver_new::dHadXi(int di)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution ret(thermo->Mtype::TPn_flash(P, T)());
    return thermo->dHadXi(P, T, di, ret);
}

double solver_new::dHsdXi(int di)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution ret(thermo->Mtype::TPn_flash(P, T)());
    return thermo->dHsdXi(P, T, di, ret);
}
void solver_new::Ln_fugacityCoefficient()
{
    autoPtr<scalarList> fugcoef;
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    double z = thermo->Z_gibbs(P, T, &comp_of);
    fugcoef.reset((thermo->ln_fugacityCoefficient(P, T, z, &comp_of)).ptr());
    ret.resize(comp.size());
    for (int i = 0;i < comp.size();i++)
    {
        ret[i] = fugcoef()[i];
    }
}
void solver_new::Ln_fugacityCoefficient(int flag)
{
    autoPtr<scalarList> fugcoef;
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    fugcoef.reset((thermo->fugacityCoefficient(P, T, flag, &comp_of)).ptr());
    ret.resize(comp.size());
    for (int i = 0;i < comp.size();i++)
    {
        ret[i] = ::log(fugcoef()[i]);
    }
}


void solver_new::ddT_Ln_fugacityCoefficient(int flag)
{
    autoPtr<scalarList> ddT_Ln_fugcoef;
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    ddT_Ln_fugcoef.reset((thermo->ddT_Ln_fugacityCoefficient(P, T, flag, &comp_of)).ptr());
    ret.resize(comp.size());
    for (int i = 0;i < comp.size();i++)
    {
        ret[i] = ddT_Ln_fugcoef()[i];
    }
}

void solver_new::ddxi_Ln_fugacityCoefficient(int di, int flag)
{
    autoPtr<scalarList> ddT_Ln_fugcoef;
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    ddT_Ln_fugcoef.reset((thermo->ddxi_Ln_fugacityCoefficient(P, T, di, flag, &comp_of)).ptr());
    ret.resize(comp.size());
    for (int i = 0;i < comp.size();i++)
    {
        ret[i] = ddT_Ln_fugcoef()[i];
    }
}

void solver_new::dvidT()
{
    Foam::scalarList comp_liq_of(comp.size(), Foam::Zero);
    Foam::scalarList comp_gas_of(comp.size(), Foam::Zero);
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }

    //std::cout << "!!comp";
    //for (int i = 0;i < comp.size();i++)
    //{
    //    std::cout << comp[i] << " ";
    //}
    //std::cout << std::endl;
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P, T)());
    //for (int i = 0;i < comp.size();i++)
    //{
    //    std::cout << "liq=" << so.X_liq()[i] << ",gas=" << so.X_gas()[i] << std::endl;
    //}
    solu.reset((thermo->dvidT(P, T, so)).ptr());

    ret.resize(comp.size());
    for (int i = 0;i < comp.size();i++)
    {
        ret[i] = solu()[i];
    }

}

double solver_new::muideal_Mole(int di)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    return  thermo->muideal_Mole(P, T, di, &comp_of);
}

bool solver_new::solveTPD_BFGS()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    bool rb, isVapor;
    Foam::autoPtr<Foam::scalarList> rl;
    std::tie(rb, rl, isVapor) = thermo->solveTPD_BFGS(P, T);
    return  rb;
}

double solver_new::Gideal()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    return  thermo->Gideal(P, T, &comp_of);
}

double solver_new::Gideal_Mole()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    return  thermo->Gideal_Mole(P, T, &comp_of);
}

double solver_new::G_Mole()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    return  thermo->G_Mole(P, T, &comp_of);
}

double solver_new::G_departure_Mole()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    return  thermo->G_departure_Mole(P, T, &comp_of);
}

double solver_new::Gibbs_single()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    return  thermo->G_TPD(P, T, &comp_of);
}


double solver_new::A_single()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    return  thermo->A_TPD(P, T, &comp_of);
}
double solver_new::z_single()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    return  thermo->Z_gibbs(P, T, &comp_of);
}
void solver_new::dvidP()
{
    Foam::scalarList comp_liq_of(comp.size(), Foam::Zero);
    Foam::scalarList comp_gas_of(comp.size(), Foam::Zero);
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }

    //std::cout << "!!comp";
    //for (int i = 0;i < comp.size();i++)
    //{
    //    std::cout << comp[i] << " ";
    //}
    //std::cout << std::endl;
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P, T)());
    //for (int i = 0;i < comp.size();i++)
    //{
    //    std::cout << "liq=" << so.X_liq()[i] << ",gas=" << so.X_gas()[i] << std::endl;
    //}
    solu.reset((thermo->dvidP(P, T, so)).ptr());

    ret.resize(comp.size());
    for (int i = 0;i < comp.size();i++)
    {
        ret[i] = solu()[i];
    }

}

void solver_new::dvidXi(int di)
{
    Foam::scalarList comp_liq_of(comp.size(), Foam::Zero);
    Foam::scalarList comp_gas_of(comp.size(), Foam::Zero);
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P, T)());

    solu.reset((thermo->dvidXi(P, T, di, so)).ptr());

    ret.resize(comp.size());
    for (int i = 0;i < comp.size();i++)
    {
        ret[i] = solu()[i];
    }

}

double  solver_new::dTdP_HP()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P, T)());
    return thermo->dTdP_HP(P, T, so);

}

double  solver_new::dTdH_HP()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P, T)());
    return thermo->dTdH_HP(P, T, so);

}
double  solver_new::dTdXi_HP(int di)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P, T)());
    return thermo->dTdXi_HP(P, T, di, so);

}

double  solver_new::rho()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P, T)());
    return thermo->rho(P, T);

}

double  solver_new::dZdT()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P, T)());
    return thermo->dZdT(P, T, so);

}

double  solver_new::Z()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P, T)());
    return thermo->Z(P, T, so);

}

double  solver_new::drhodT()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P, T)());
    return thermo->drhodT(P, T, so);

}

double  solver_new::drhodP()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P, T)());
    return thermo->drhodP(P, T, so);

}

double  solver_new::drhodXi(int di)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P, T)());
    return thermo->drhodXi(P, T, di, so);
}

double  solver_new::dZdXi(int di)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P, T)());
    return thermo->dZdXi(P, T, di, so);
}

double  solver_new::drhoPdH_HP()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P, T)());
    return thermo->drhoPdH_HP(P, T, so);
}
double  solver_new::drhoPdP_HP()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P, T)());
    return thermo->drhoPdP_HP(P, T, so);
}

double  solver_new::drhoPdH_HsP()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P, T)());
    return thermo->drhoPdH_HsP(P, T, so);
}

double  solver_new::dTdXi_HsP(int di)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P, T)());
    return thermo->dTdXi_HsP(P, T, di, so);
}
double  solver_new::drhoPdP_HsP()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P, T)());
    return thermo->drhoPdP_HsP(P, T, so);
}

double  solver_new::drhodP_HP()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P, T)());
    return thermo->drhodP_HP(P, T, so);
}
double  solver_new::drhoPdXi_HP(int di)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P, T)());
    return thermo->drhoPdXi_HP(P, T, di, so);
}
double  solver_new::drhoPdXi_HsP(int di)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P, T)());
    return thermo->drhoPdXi_HsP(P, T, di, so);
}
double  solver_new::drhodXi_HP(int di)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P, T)());
    return thermo->drhodXi_HP(P, T, di, so);
}

double  solver_new::dvfdP_HP()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P, T)());
    return thermo->dvfdP_HP(P, T, so);
}
double  solver_new::dvfdH_HP()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P, T)());
    return thermo->dvfdH_HP(P, T, so);
}

double solver_new::dvfdXi_HP(int di)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P, T)());
    return thermo->dvfdXi_HP(P, T, di, so);
}


double  solver_new::dvfdP_HsP()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P, T)());
    return thermo->dvfdP_HsP(P, T, so);
}
double  solver_new::dvfdH_HsP()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P, T)());
    return thermo->dvfdH_HsP(P, T, so);
}

double solver_new::dvfdXi_HsP(int di)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P, T)());
    return thermo->dvfdXi_HsP(P, T, di, so);
}

double solver_new::T_HsP(double h, double P, double T0)
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);

    double Ttemp = T0;
    double htemp = thermo->Hs(P, Ttemp);
    while (fabs(htemp - h) > 1e-5)
    {
        //cout << Ttemp << "," << htemp << std::endl;
        Mtype::solution so(thermo->Mtype::TPn_flash(P, Ttemp)());
        cout << thermo->Cp(P, Ttemp) << "," << thermo->Cp_Hs(P, Ttemp, so) << std::endl;
        Ttemp -= (htemp - h) / thermo->Cp_Hs(P, Ttemp, so);
        htemp = thermo->Hs(P, Ttemp);
    }
    return Ttemp;
}

void solver_new::drhoPdXHP_HsP()
{
    Foam::scalarList comp_of(comp.size(), Foam::Zero);
    autoPtr<scalarList> solu;

    for (int i = 0;i < comp.size();i++)
    {
        comp_of[i] = comp[i];
    }
    thermo->setX(comp_of);
    Mtype::solution so(thermo->Mtype::TPn_flash(P, T)());
    autoPtr<scalarList> grad = thermo->drhoPdXHP_HsP(P, T, so);
    ret.resize(grad().size());
    for (int i = 0;i < ret.size();i++)
    {
        ret[i] = grad()[i];
    }

}


int fun_my(int n)
{
    if (n == 0)
        return 1;
    return n * fun_my(n - 1);
}