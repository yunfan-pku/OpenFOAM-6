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

using namespace Foam;
void regularize(scalarList &a)
{
  scalar sum = 0;
  forAll(a, i)
  {
    sum += a[i];
  }
  forAll(a, i)
  {
    a[i] /= sum;
  }
}

void tablegen(dictionary &dict)
{
  speciesTable s(dict.lookup("species"));
  PengRobinsonM<specie> PR(s, dict);
  string ThermoTable;
  scalar vaporfra = 0.5;
  scalarList equalconstant(s.size(), Zero);
  scalarList comp_liq(s.size(), Zero);
  scalarList comp_gas(s.size(), Zero);
  scalarList comp_liq0(s.size(), Zero);
  scalarList comp_gas0(s.size(), Zero);

  ofstream output;
  word outputfilename(word(dict.lookup("outputfilename")));
  output.open(outputfilename);
  scalarList comp(dict.lookup("component_ini"));
  scalarList comp1(dict.lookup("component_fin"));
  label num_comp(readLabel(dict.lookup("num_comp")));
  regularize(comp);
  regularize(comp1);
  scalarList det_comp(comp.size());
  scalarList t_comp(comp.size());
  if (num_comp > 0)
    forAll(comp, i)
    {
      det_comp[i] = (comp1[i] - comp[i]) / num_comp;
    }

  forAll(comp, isp)
  {
    comp_liq0[isp] = comp[isp];
    comp_gas0[isp] = comp[isp];
  }

  scalar pressure(readScalar(dict.lookup("pressure_ini")));
  scalar pressure1(readScalar(dict.lookup("pressure_fin")));
  label num_pres(readLabel(dict.lookup("num_pre")));
  scalar det_pres;
  if (num_pres > 0)
    det_pres = (pressure1 - pressure) / num_pres;

  scalar t_pres;

  scalar temperature(readScalar(dict.lookup("temperature_ini")));
  scalar temperature1(readScalar(dict.lookup("temperature_fin")));
  label num_temp(readLabel(dict.lookup("num_temp")));
  scalar det_temp;
  if (num_temp > 0)
    det_temp = (temperature1 - temperature) / num_temp;
  scalar t_temp;
  int state;
  label PCFlag(readLabel(dict.lookup("PCFlag")));
  label HEFlag(readLabel(dict.lookup("HEFlag")));
   output << "A"<<std::endl;
  for (label ncomp = 0; ncomp <= num_comp; ncomp++)
  {
    forAll(comp, i)
    {
      t_comp[i] = comp[i] + det_comp[i] * ncomp;
    }
    scalar mw_mixture = 0.0; //kg/mol
    scalar rhoMixture = 0.0; //kg/m3
    scalar sieMixture = 0.0; //J/kg
    scalar CpMixture = 0.0;  //J/kgK
    scalar CvMixture = 0.0;  //J/kgK
    scalar CsMixture = 0.0;  //cm/s
    scalar ZMixture = 0.0;   //[-]
    //scalar ZMixture2  = 0.0;  //[-]
    //scalar GammaMixture = 0.0; //[-]

    scalar kappaMixture = 0.0;
    scalar muMixture = 0.0;
    scalar rho_gas = 0.0, rho_liq = 0.0;
    scalar Cp_gas = 0.0, Cp_liq = 0.0;
    scalar Cv_gas = 0.0, Cv_liq = 0.0;
    scalar alphagas = 0.0, ygas = 0.0;
    scalar kappa_gas = 0.0, kappa_liq = 0.0;
    scalar mu_gas = 0.0, mu_liq = 0.0;
    scalar Ct_gas = 0.0, Ct_liq = 0.0;
    scalar Cs_gas = 0.0, Cs_liq = 0.0;
    scalar VspecificGas = 0.0, VspecificLiq = 0.0;
    scalar ZGas = 0.0, ZLiq = 0.0;
    scalar Dij_binary_high = 0.0;
    scalar hMixture=0;

    for (label npres = 0; npres <= num_pres; npres++)
    {
      t_pres = pressure + npres * det_pres;

      for (label ntemp = 0; ntemp <= num_temp; ntemp++)
      {
        t_temp = temperature + ntemp * det_temp;

        if (PCFlag == 0)
        {
          forAll(comp, isp)
          {
            comp_liq[isp] = comp_liq0[isp];
            comp_gas[isp] = comp_gas0[isp];
          }

          PR.TPn_flash(t_pres, t_temp, t_comp, comp_liq, comp_gas, vaporfra, equalconstant);
          //Info<<"vaporfra= "<<vaporfra<<endl;

          forAll(comp, isp)
          {
            comp_liq0[isp] = comp_liq[isp];
            comp_gas0[isp] = comp_gas[isp];
          }

          if (vaporfra > 1.0)
          {
            vaporfra = 0.9999999;
          }
          else if (vaporfra < 0.0)
          {
            vaporfra = 1.0e-07;
          }
          if (vaporfra >= 0.9999999)
          {
            state = 1;
          }
          else if (vaporfra <= 1.0e-07)
          {
            state = 0;
          }
          else
          {
            state = 2;
          }
          alphagas = PR.Evaluate_alpha(t_pres, t_temp, vaporfra, comp_liq, comp_gas, t_comp);
          scalar mw_gas = PR.mwmix(comp_gas); //kg/mol
          //mw_mixture	        = PR.mwmix(comp);

          scalar mw_liq = PR.mwmix(comp_liq); //kg/mol
          mw_mixture = mw_gas * vaporfra + mw_liq * (1.0 - vaporfra);
          ygas = vaporfra * mw_gas / mw_mixture;

          //Thermal properties
          VspecificGas = PR.volmmix_phase(0, t_pres, t_temp, comp_gas); //m3/mol
          VspecificLiq = PR.volmmix_phase(1, t_pres, t_temp, comp_liq);
          ZGas = t_pres * VspecificGas / (RR * 1.0e-03 * t_temp);
          ZLiq = t_pres * VspecificLiq / (RR * 1.0e-03 * t_temp);
          ZMixture = ZGas * alphagas + ZLiq * (1.0 - alphagas);
          rhoMixture = t_pres * mw_mixture / (ZMixture * RR * 1.0e-03 * t_temp);

          rho_gas = PR.rhomix_phase(0, t_pres, t_temp, comp_gas); //kg/m3
          rho_liq = PR.rhomix_phase(1, t_pres, t_temp, comp_liq);
          //rhoMixture = rho_gas * alphagas + rho_liq * (1.0 - alphagas);

          scalar sieReal_gas = PR.siemix_phase(0, t_pres, t_temp, comp_gas); //J/kg
          scalar sieReal_liq = PR.siemix_phase(1, t_pres, t_temp, comp_liq);
          sieMixture = sieReal_gas * ygas + sieReal_liq * (1.0 - ygas);

          Cp_gas = PR.cpmix_phase(0, t_pres, t_temp, comp_gas); //J/kgK
          Cp_liq = PR.cpmix_phase(1, t_pres, t_temp, comp_liq);
          CpMixture = Cp_gas * ygas + Cp_liq * (1.0 - ygas);

          Cv_gas = PR.cvmix_phase(0, t_pres, t_temp, comp_gas); //J/kgK
          Cv_liq = PR.cvmix_phase(1, t_pres, t_temp, comp_liq);
          CvMixture = Cv_gas * ygas + Cv_liq * (1.0 - ygas);
          //GammaMixture        = CpMixture/CvMixture;

          PR.soundspeedmix(0, t_pres, t_temp, comp_gas, Ct_gas, Cs_gas);
          PR.soundspeedmix(1, t_pres, t_temp, comp_liq, Ct_liq, Cs_liq);

          CsMixture = ::sqrt(1.0 / (rhoMixture * (alphagas / (rho_gas * sqr(Cs_gas)) +
                                                  (1.0 - alphagas) / (rho_liq * sqr(Cs_liq)))));
          //CtMixture           = CsMixture * sqr(GammaMixture);

          //Transport	Properties
          PR.kappa_phase(0, t_pres, t_temp, comp_gas, kappa_gas, mu_gas); //W/mK
          PR.kappa_phase(1, t_pres, t_temp, comp_liq, kappa_liq, mu_liq);
          kappaMixture = kappa_gas * ygas + kappa_liq * (1.0 - ygas); //W/mK
          

          scalar hReal_gas = PR.hmix_phase(0, t_pres, t_temp, comp_gas); //J/kg
          scalar hReal_liq = PR.hmix_phase(1, t_pres, t_temp, comp_liq);
          hMixture = hReal_gas * ygas + hReal_liq * (1.0 - ygas);
          //hMixture=sieMixture + t_pres*mw_mixture / rhoMixture;
          //muMixture = 1.0e+06 * (mu_gas * ygas + mu_liq * (1.0 - ygas)); //up
          muMixture = (mu_gas * ygas + mu_liq * (1.0 - ygas));

          //Dij_binary_high	= PR.Dij_highP(t_pres, t_temp, comp);//cm2/s

          //ZMixture2	        = t_pres*mw_mixture/(rhoMixture*RR*1.0e-03*t_temp);
        }
        else
        {
          rhoMixture = PR.rhomix_phase(0, t_pres, t_temp, t_comp);
          sieMixture = PR.siemix_phase(0, t_pres, t_temp, t_comp);
          CpMixture = PR.cpmix_phase(0, t_pres, t_temp, t_comp);
          CvMixture = PR.cvmix_phase(0, t_pres, t_temp, t_comp);
          PR.soundspeedmix(0, t_pres, t_temp, t_comp, Ct_gas, CsMixture);
          VspecificGas = PR.volmmix_phase(0, t_pres, t_temp, t_comp);
          ZMixture = t_pres * VspecificGas / (RR * 1.0e-03 * t_temp);
          PR.kappa_phase(0, t_pres, t_temp, t_comp, kappaMixture, muMixture);
          alphagas = 1;

          //Dij_binary_high	= PR.Dij_highP_pesudo(t_pres, t_temp, comp);
        }

        //output<<t_temp<<' '<<' '<<t_pres*1.0e-05<<' '<<' '<<//comp[0]<<' '<<' '<<

        //rhoMixture;<<' '<<' '<<CpMixture<<' '<<' '<<

        //kappaMixture;//<<' '<<' '<<muMixture<<' '<<' '<<

        //ZMixture;<<' '<<' '<<Dij_binary_high*1.0e+06;
        //Info << t_temp << "  " << t_pres << "  " << t_comp[0] << "  "<<muMixture<<endl;
        if (HEFlag == 1)
        {
          output << t_temp << "  " << t_pres << "  " << t_comp[0] << "  " <<

              rhoMixture << "  " << sieMixture << "  " <<

              CpMixture << "  " << CvMixture << "  " << CsMixture << "  " <<

              kappaMixture << "  " << muMixture << "  " <<

              ZMixture << "  " << Dij_binary_high << "  " <<

              comp_liq[0] << "  " << comp_gas[0] << "  " << alphagas;
        }
        else if (HEFlag == 0)
        {
          output << t_temp << "  " << t_pres << "  " << t_comp[0] << "  " << hMixture << "  " << t_pres  << "  " << CpMixture << "  " << CvMixture << "  " << CsMixture << "  " <<kappaMixture << "  " << muMixture << "  " <<ZMixture << "  " << Dij_binary_high << "  " <<comp_liq[0] << "  " << comp_gas[0] << "  " << alphagas<<" "<<rhoMixture<<" "<<vaporfra;
          //          1               2                 3                     4                    5                        6                 7                         8                9                     10                  11                   12                          13                      14                  15       16                17
        }

        Info << "t=" << t_temp << "  p=" << t_pres << " vaporfra=" << vaporfra << endl;

        //output << t_temp << ',' << ' ' << t_pres <<" ," << state<<","<<t_pres/(rhoMixture* t_temp);

        output << "\n";
        //Info<<"PY:Screen_Tabulation"<<", kappa_gas="<<kappa_gas<<", ygas="<<ygas<<", kappa_liq="<<kappa_liq<<", kappaMixture="<<kappaMixture<<"\n"<< endl;
      }
    }
  }
  output.close();
  return;
}

void TPdiagram(dictionary &dict)
{
  speciesTable s(dict.lookup("species"));
  PengRobinsonM<specie> PR(s, dict);
  string ThermoTable;
  scalar vaporfra = 0.9;
  scalarList equalconstant(s.size(), Zero);
  scalarList comp_liq(s.size(), Zero);
  scalarList comp_gas(s.size(), Zero);
  scalarList comp_liq0(s.size(), Zero);
  scalarList comp_gas0(s.size(), Zero);

  ofstream output;
  word outputfilename(word(dict.lookup("outputfilename")));
  output.open(outputfilename);
  scalarList comp(dict.lookup("component"));
  regularize(comp);

  forAll(comp, isp)
  {
    comp_liq0[isp] = comp[isp];
    comp_gas0[isp] = comp[isp];
  }

  scalar pressure(readScalar(dict.lookup("pressure_ini")));
  scalar pressure1(readScalar(dict.lookup("pressure_fin")));
  label num_pres(readLabel(dict.lookup("num_pre")));
  scalar det_pres;
  if (num_pres > 0)
    det_pres = (pressure1 - pressure) / num_pres;

  scalar t_pres;

  scalar temperature(readScalar(dict.lookup("temperature_ini")));
  scalar temperature1(readScalar(dict.lookup("temperature_fin")));
  label num_temp(readLabel(dict.lookup("num_temp")));
  scalar det_temp;
  if (num_temp > 0)
    det_temp = (temperature1 - temperature) / num_temp;
  scalar t_temp;

  scalar tpdtest = 0.0;
  scalar tempphase_pt[10] = {0.0};
  scalar presphase = 0.0;
  scalar tpdlast = 1.0;
  scalar dettpd = 0.0;
  int ndettpd = 0;

  for (label npres = 0; npres < num_pres; npres++)
  {
    tpdlast = 1.0;
    t_pres = pressure + npres * det_pres;
    ndettpd = 0;

    for (label ntemp = 0; ntemp < num_temp; ntemp++)
    {
      t_temp = temperature + ntemp * det_temp;

      PR.TPn_flash(t_pres, t_temp, comp, comp_liq, comp_gas, vaporfra, equalconstant);
      //alphagas = PR.Evaluate_alpha(press, temp, vaporfra, comp_liq, comp_gas,comp);

      scalar detCompPhase = 0.0;

      tpdtest = phasestate(vaporfra);
      dettpd = tpdtest - tpdlast; //if tpd change first time, abs(dettpd) > 0;
      scalar maxtpd = max(tpdtest, tpdlast);
      /*
		    if (dettpd != 0.0 && maxtpd == 2.0)
		    {
			ndettpd = ndettpd + 1;
			int ndettpd0 = ndettpd-1;
			for (label nwrite = ndettpd0; nwrite < ndettpd; nwrite++)
			{
			    tempphase_pt[nwrite] = temp;
			}
		    }
            */
      if (dettpd != 0.0)
      {
        double templ = t_temp - det_temp;
        double tempr = t_temp;
        double tempm = (templ + tempr) / 2;
        double tpdtest_t = 0;
        while (tempr - templ > 1e-3)
        {
          tempm = (templ + tempr) / 2;
          PR.TPn_flash(t_pres, tempm, comp, comp_liq, comp_gas, vaporfra, equalconstant);
          tpdtest_t = phasestate(vaporfra);
          if (tpdtest_t == tpdlast)
          {
            templ = tempm;
          }
          else //(tpdtest_t == tpdtest)
          {
            tempr = tempm;
          }
          Info << templ << " " << tempr << endl;
        }
        tempphase_pt[ndettpd++] = tempm;
      }
      presphase = t_pres;

      //if (dettpd != 0.0 && ndettpd == 2)
      if (ntemp == num_temp - 1)
      {
        for (int i = 0; i < ndettpd; i++)
          output << tempphase_pt[i] << ',';
        output << presphase * 1.0e-05;
        output << "\n";
      }
      tpdlast = tpdtest;

      Info << "PY:Screen_TPnFlash"
           << ", x=" << comp_liq[0] << ", y=" << comp_gas[0] << ", vap=" << vaporfra << "\n"
           << endl;
    }
  }
  //}
  output.close();
}

void PXdiagram(dictionary &dict)
{
  speciesTable s(dict.lookup("species"));
  PengRobinsonM<specie> PR(s, dict);
  string ThermoTable;
  scalar vaporfra = 0.9;
  scalarList equalconstant(s.size(), Zero);
  scalarList comp_liq(s.size(), Zero);
  scalarList comp_gas(s.size(), Zero);
  scalarList comp_liq0(s.size(), Zero);
  scalarList comp_gas0(s.size(), Zero);

  ofstream output;
  word outputfilename(word(dict.lookup("outputfilename")));
  output.open(outputfilename);
  scalarList comp(dict.lookup("component_ini"));
  scalarList comp1(dict.lookup("component_fin"));
  label num_comp(readLabel(dict.lookup("num_comp")));
  regularize(comp);
  regularize(comp1);
  scalarList det_comp(comp.size());
  scalarList t_comp(comp.size());
  if (num_comp > 0)
    forAll(comp, i)
    {
      det_comp[i] = (comp1[i] - comp[i]) / num_comp;
    }

  forAll(comp, isp)
  {
    comp_liq0[isp] = comp[isp];
    comp_gas0[isp] = comp[isp];
  }

  scalar pressure(readScalar(dict.lookup("pressure_ini")));
  scalar pressure1(readScalar(dict.lookup("pressure_fin")));
  label num_pres(readLabel(dict.lookup("num_pre")));
  scalar det_pres;
  if (num_pres > 0)
    det_pres = (pressure1 - pressure) / num_pres;

  scalar t_pres;

  scalar temperature(readScalar(dict.lookup("temperature")));

  scalar tpdtest = 0.0;
  scalar compphase_pt[10] = {0.0};
  scalar presphase = 0.0;
  scalar tpdlast = 1.0;
  scalar dettpd = 0.0;
  int ndettpd = 0;

  scalarList comple(s.size(), Zero);
  scalarList compr(s.size(), Zero);
  scalarList compm(s.size(), Zero);

  for (label npres = 0; npres < num_pres; npres++)
  {
    tpdlast = 0;
    t_pres = pressure + npres * det_pres;
    ndettpd = 0;

    for (label ncomp = 0; ncomp <= num_comp; ncomp++)
    {
      forAll(comp, i)
      {
        t_comp[i] = comp[i] + det_comp[i] * ncomp;
      }

      PR.TPn_flash(t_pres, temperature, t_comp, comp_liq, comp_gas, vaporfra, equalconstant);
      //alphagas = PR.Evaluate_alpha(press, temp, vaporfra, comp_liq, comp_gas,comp);

      scalar detCompPhase = 0.0;

      tpdtest = phasestate(vaporfra);
      dettpd = tpdtest - tpdlast; //if tpd change first time, abs(dettpd) > 0;
      scalar maxtpd = max(tpdtest, tpdlast);
      /*
		    if (dettpd != 0.0 && maxtpd == 2.0)
		    {
			ndettpd = ndettpd + 1;
			int ndettpd0 = ndettpd-1;
			for (label nwrite = ndettpd0; nwrite < ndettpd; nwrite++)
			{
			    tempphase_pt[nwrite] = temp;
			}
		    }
            */

      if (dettpd != 0.0)
      {

        forAll(comp, i)
        {
          comple[i] = t_comp[i] - det_comp[i];
          compr[i] = t_comp[i];
          compm[i] = (comple[i] + compr[i]) / 2;
        }
        double tpdtest_t = 0;
        while (compr[0] - comple[0] > 1e-3)
        {
          forAll(comp, i)
          {
            compm[i] = (comple[i] + compr[i]) / 2;
          }
          PR.TPn_flash(t_pres, temperature, compm, comp_liq, comp_gas, vaporfra, equalconstant);
          tpdtest_t = phasestate(vaporfra);
          if (tpdtest_t == tpdlast)
          {
            forAll(comp, i)
            {
              comple[i] = compm[i];
            }
          }
          if (tpdtest_t == tpdtest)
          {
            forAll(comp, i)
            {
              compr[i] = compm[i];
            }
          }
        }
        compphase_pt[ndettpd++] = compm[0];
      }
      presphase = t_pres;

      //if (dettpd != 0.0 && ndettpd == 2)
      if (ncomp == num_comp - 1)
      {
        for (int i = 0; i < ndettpd; i++)
          output << compphase_pt[i] << ' ' << ' ' << ' ';
        output << presphase * 1.0e-05;
        output << "\n";
      }
      tpdlast = tpdtest;

      Info << "PY:Screen_TPnFlash"
           << ", x=" << comp_liq[0] << ", y=" << comp_gas[0] << ", vap=" << vaporfra << "\n"
           << endl;
    }
  }
  //}
  output.close();
}

void TXdiagram(dictionary &dict)
{
  speciesTable s(dict.lookup("species"));
  PengRobinsonM<specie> PR(s, dict);
  string ThermoTable;
  scalar vaporfra = 0.9;
  scalarList equalconstant(s.size(), Zero);
  scalarList comp_liq(s.size(), Zero);
  scalarList comp_gas(s.size(), Zero);
  scalarList comp_liq0(s.size(), Zero);
  scalarList comp_gas0(s.size(), Zero);

  ofstream output;
  word outputfilename(word(dict.lookup("outputfilename")));
  output.open(outputfilename);
  scalarList comp(dict.lookup("component_ini"));
  scalarList comp1(dict.lookup("component_fin"));
  label num_comp(readLabel(dict.lookup("num_comp")));
  regularize(comp);
  regularize(comp1);
  scalarList det_comp(comp.size());
  scalarList t_comp(comp.size());
  if (num_comp > 0)
    forAll(comp, i)
    {
      det_comp[i] = (comp1[i] - comp[i]) / num_comp;
    }

  forAll(comp, isp)
  {
    comp_liq0[isp] = comp[isp];
    comp_gas0[isp] = comp[isp];
  }

  scalar temperature(readScalar(dict.lookup("temperature_ini")));
  scalar temperature1(readScalar(dict.lookup("temperature_fin")));
  label num_temp(readLabel(dict.lookup("num_temp")));
  scalar det_temp;
  if (num_temp > 0)
    det_temp= (temperature1 - temperature) / num_temp;

  scalar t_temp;

  scalar pressure(readScalar(dict.lookup("pressure")));

  scalar tpdtest = 0.0;
  scalar compphase_pt[10] = {0.0};
  scalar tempphase = 0.0;
  scalar tpdlast = 1.0;
  scalar dettpd = 0.0;
  int ndettpd = 0;

  scalarList comple(s.size(), Zero);
  scalarList compr(s.size(), Zero);
  scalarList compm(s.size(), Zero);

  for (label ntemp = 0; ntemp < num_temp; ntemp++)
  {
    tpdlast = 0;
    t_temp = temperature + ntemp * det_temp;
    ndettpd = 0;

    for (label ncomp = 0; ncomp <= num_comp; ncomp++)
    {
      forAll(comp, i)
      {
        t_comp[i] = comp[i] + det_comp[i] * ncomp;
      }

      PR.TPn_flash(pressure, t_temp, t_comp, comp_liq, comp_gas, vaporfra, equalconstant);
      //alphagas = PR.Evaluate_alpha(press, temp, vaporfra, comp_liq, comp_gas,comp);

      scalar detCompPhase = 0.0;

      tpdtest = phasestate(vaporfra);
      dettpd = tpdtest - tpdlast; //if tpd change first time, abs(dettpd) > 0;
      scalar maxtpd = max(tpdtest, tpdlast);
      /*
		    if (dettpd != 0.0 && maxtpd == 2.0)
		    {
			ndettpd = ndettpd + 1;
			int ndettpd0 = ndettpd-1;
			for (label nwrite = ndettpd0; nwrite < ndettpd; nwrite++)
			{
			    tempphase_pt[nwrite] = temp;
			}
		    }
            */

      if (dettpd != 0.0)
      {

        forAll(comp, i)
        {
          comple[i] = t_comp[i] - det_comp[i];
          compr[i] = t_comp[i];
          compm[i] = (comple[i] + compr[i]) / 2;
        }
        double tpdtest_t = 0;
        while (compr[0] - comple[0] > 1e-3)
        {
          forAll(comp, i)
          {
            compm[i] = (comple[i] + compr[i]) / 2;
          }
          PR.TPn_flash(pressure, t_temp, compm, comp_liq, comp_gas, vaporfra, equalconstant);
          tpdtest_t = phasestate(vaporfra);
          if (tpdtest_t == tpdlast)
          {
            forAll(comp, i)
            {
              comple[i] = compm[i];
            }
          }
          if (tpdtest_t == tpdtest)
          {
            forAll(comp, i)
            {
              compr[i] = compm[i];
            }
          }
        }
        compphase_pt[ndettpd++] = compm[0];
      }
      tempphase = t_temp;

      //if (dettpd != 0.0 && ndettpd == 2)
      if (ncomp == num_comp - 1)
      {
        for (int i = 0; i < ndettpd; i++)
          output << compphase_pt[i] << ',';
        output << tempphase;
        output << "\n";
      }
      tpdlast = tpdtest;

      Info << "PY:Screen_TPnFlash"
           << ", x=" << comp_liq[0] << ", y=" << comp_gas[0] << ", vap=" << vaporfra << "\n"
           << endl;
    }
  }
  //}
  output.close();
}

void phase(dictionary &dict)
{
  speciesTable s(dict.lookup("species"));
  PengRobinsonM<specie> PR(s, dict);

  ofstream output;
  word outputfilename(word(dict.lookup("outputfilename")));
  output.open(outputfilename);

  ifstream input;
  word inputfilename(word(dict.lookup("inputfilename")));
  input.open(inputfilename);
  string str;
  input >> str;
  Info << str << endl;
  double c1, c2, p, t, x1, x2, x3;
  char C;
  scalarList t_comp(2);
  scalarList comp_liq(2);
  scalarList comp_gas(2);
  double vaporfra;
  scalarList equalconstant(2);
  double state, st;
  t_comp[0] = c1;
  t_comp[1] = c2;
  while (!input.eof())
  {
    input >> c1 >> C >> p >> C >> t >> C >> c2 >> C >> x1 >> C >> x2 >> C >> x3 >> C >> st;
    Info << c1 << " " << p << " " << t << " " << c2 << " " << x1 << " " << x2 << " " << x3 << " " << endl;
    c1 = c1 / 44 / (c1 / 44 + c2 / 18);
    c2 = 1 - c1;
    p = 10000000;
    t = 350;
    t_comp[0] = c1;
    t_comp[1] = c2;

    PR.TPn_flash(p, t, t_comp, comp_liq, comp_gas, vaporfra, equalconstant);
    if (vaporfra > 1.0)
    {
      vaporfra = 0.9999999;
    }
    else if (vaporfra < 0.0)
    {
      vaporfra = 1.0e-07;
    }
    if (vaporfra >= 0.9999999)
    {
      state = 1;
    }
    else if (vaporfra <= 1.0e-07)
    {
      state = 0;
    }
    else
    {
      state = 2;
    }
    output << x1 << "," << x2 << "," << x3 << "," << state << std::endl;
  }
}

int main()
{

  /*
  dictionary dict(IFstream("system/thermo")());
  PengRobinsonS<specie> A(dict.subDict("N2"));
  Info << A.name() << endl;
  Info << A.W() << endl;
  Info << A.Y() << endl;

  Info << A.Tc() << endl;
  Info << A.Pc() << endl;
  Info << A.Vc() << endl;
  Info << A.omega() << endl;

  // wordList s(dict.lookup("species"));
  speciesTable s(dict.lookup("species"));
  //Info << s << endl;

  PengRobinsonM<specie> B(s, dict);
  scalarList Y(2);
  scalarList liq(2), gas(2),ec(2);
  scalar vf;
  Y[0] = 1;
  Y[1] = 1;
 
  B.setY(Y);
  //Info << B.Y_[0] << B.Y_[1]<< endl;
  B.TPn_flash(1000000, 400, Y, liq, gas, vf, ec);
  */
  dictionary dict(IFstream("system/thermotableDict")());
  if (readLabel(dict.lookup("TableFlag")) == 1)
    tablegen(dict.subDict("table"));
  if (readLabel(dict.lookup("PTFlag")) == 1)
    TPdiagram(dict.subDict("PT"));
  if (readLabel(dict.lookup("PXFlag")) == 1)
    PXdiagram(dict.subDict("PX"));
  if (readLabel(dict.lookup("phaseFlag")) == 1)
    phase(dict.subDict("phase"));
  if (readLabel(dict.lookup("TXFlag")) == 1)
    TXdiagram(dict.subDict("TX"));
}
