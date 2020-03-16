#include<math.h>
namespace XYYthermo{
double a;
double b;
double R;
double dadT;
double d2adT2;
void getThermo(double T)
{
    double MW =28.0134e-3;
    double Tc = 126.1;
    double pc  = 3.39e+06;
    double rhoc = 310.915;
    double omega  = 0.0377;
    R = 8.314/MW;
    double c = 0.37464 + 1.54226*omega - 0.26992*omega*omega;
    a= 0.457236*(R*Tc)*(R*Tc)/pc*(1+c*(1-std::sqrt(T/Tc)))*(1+c*(1-std::sqrt(T/Tc)));
    b= 0.077796*R*Tc/pc;
    double G = c*std::sqrt(T/Tc)/(1+c*(1-std::sqrt(T/Tc)));
    dadT = -1/T*a*G;
    d2adT2 = 0.457236*R*R/T/2*c*(1+c)*Tc/pc*std::sqrt(Tc/T);
}
double getPfromTandRho(double T,double rho)
{
    double v = 1/rho;
    return R*T/(v-b) - a/(v*v+2*v*b-b*b);
}
double getTfromPandRho(double p,double rho )
{
    double CRIT = 1.0e-6;
    double v=1.0/rho;
    double T = 300;
    double p_n = getPfromTandRho(T,rho);
    double diff = p_n - p;

    while (fabs(diff) > CRIT)
    {
      getThermo(T);
      double dpdT = R/(v-b) - 1/(v*v+2*v*b-b*b)*dadT;
      T = T - diff/dpdT;
      p_n = getPfromTandRho(T,rho);
      diff = p_n - p;
    }
    return T;
}

double getsos(double rho,double p)
{
   double an[7]={3.531005280E+00,-1.236609870E-04,-5.029994370E-07,2.435306120E-09,-1.408812350E-12,-1.046976280E+03,2.967474680E+00};
   double v=1/rho;
   double T=getTfromPandRho(p,rho);
   getThermo(T);
   double cp_ideal = R*(an[0] + an[1]*T + an[2]*T*T + an[3]*T*T*T + an[4]*T*T*T*T);
   double cv_ideal = cp_ideal - R;
   double dpdT = R/(v-b) - dadT/(v*v+2*v*b-b*b);
   double dpdv = -R*T/((v-b)*(v-b))*(1-2*a*(1/(R*T*(v+b)*((v*v+2*v*b-b*b)/(v*v-b*b))*((v*v+2*v*b-b*b)/(v*v-b*b)))));
   double K1 = 1/std::sqrt(8)/b*std::log((v+(1-std::sqrt(2))*b)/(v+(1+std::sqrt(2))*b));
   double cp = cv_ideal - T*(dpdT)*(dpdT)/dpdv - K1*T*d2adT2;
   double kT = -1/v/dpdv;
   double av = -dpdT/v/dpdv;
   double ks = kT - v*T*av*av/cp;
   return  std::sqrt(1/rho/ks);
}
}