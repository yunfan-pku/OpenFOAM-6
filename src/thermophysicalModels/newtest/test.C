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

using namespace Foam;
int main()
{
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
    Info << mst[1].rho(100000, 1000) << endl;
    Info<<mst.X()<<endl;
    Info<<mst.Y()<<endl;
    PengRobinsonMixture<multispecie<PengRobinson<specie>>> pms(dict);
    pms.setX(X);
    Info<<pms.rho(100000, 1000)<<endl;
    Info<<pms.fugacityCoefficient(100000, 1000)()<<endl;
}