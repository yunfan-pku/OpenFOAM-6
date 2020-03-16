#ifndef thermotable_H
#define thermotable_H

#include "volFields.H"
#include "PtrList.H"
#include "specie.H"
#include "HYmultiComponentMixture.H"
#include "table.H"
#include "tabulation.H"

namespace Foam
{

template <class ThermoType>
class HYmultiComponentMixture;

template <class ThermoType>
class thermotable
    : public ThermoType
{
public:
    static tabulation table;
    HYmultiComponentMixture<ThermoType> &multi_mixture_;
    //PtrList<volScalarField>& Y_;
    scalarList Y_;
    /*
    bool boundaryflag;
    label patch_;
    label face_;
    label cell_;
*/
    inline thermotable(
        HYmultiComponentMixture<ThermoType> &multi_mixture_in,
        const word &name,
        const ThermoType &st);

    static word typeName()
    {
        //return "thermotable<" + ThermoType::typeName() + '>';
        return ThermoType::typeName();
    };

    static word heName()
    {
        return "e";
    };

    inline void from_Y_to_z(
        const scalarList &Yinp,
        scalarList &zi) const;

    inline scalar THE(
        const scalar H,
        const scalar p,
        const scalar T0) const;

    inline scalar psi(
        const scalar p,
        const scalar T) const;

    inline scalar mu(
        const scalar p,
        const scalar T) const;

    inline scalar volspecific(
        const scalar p,
        const scalar T) const;

    inline scalar rho //kg/m3
        (
            const scalar p,
            const scalar T) const;

    inline double get_property_from_tpn_reverse(
        const double rho,
        const double sie,
        const double zi,
        int index) const;

    double from_tpn_to_property(
        const double T,
        const double p,
        const double zi,
        int index) const;

    scalar from_tpn_reverse_to_T(
        const scalar propval,
        const scalar p,
        const scalar T0) const;

    scalar mw //g/mol
        (
            const scalar p,
            const scalar T) const;
    scalar extrapolate_temperature_tpn(
        const scalar propval,
        const scalar temp_guess,
        const scalar p) const;
    scalar alphah //make sure the unit
        (
            const scalar p,
            const scalar T) const;
    scalar kappa(
        const scalar p,
        const scalar T) const;
    scalar Cp(
        const scalar p,
        const scalar T) const;
    scalar Cpv(
        const scalar p,
        const scalar T) const;
    scalar HE(
        const scalar p,
        const scalar T) const;
    scalar CpByCpv(
        const scalar p,
        const scalar T) const;
    scalar gamma(
        const scalar p,
        const scalar T) const;
};

template <class ThermoType>
tabulation thermotable<ThermoType>::table("ThermoTable.txt");

} // namespace Foam

#include "thermotableI.H"

#ifdef NoRepository
#include "thermotable.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //