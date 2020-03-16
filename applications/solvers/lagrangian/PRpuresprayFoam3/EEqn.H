{
  
    volScalarField& he = thermo.he();
    //volScalarField& sos = thermo.sos();

    Info<<"I'm 2nd!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";

    //volScalarField gammastar=rho*sos*sos/p;
    //volScalarField estar=he-p/rho/(gammastar-1);
    

    /*
    std::stringstream sstream;
    sstream<<timei<<"s2"<<".csv";
    std::string s;
    sstream>>s;
    std::ofstream fout(s);

    forAll(he,i)
    {
      fout<<p[i]<<",";
    } 
    fout<<std::endl;
    */


    fvScalarMatrix EEqn
    (
    //    fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
    fvm::ddt(rho, he) +  fvc::div(phi, he)//fvc::div_doubleflux(phi, p,rho,gammastar,estar)//fvc::div(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
      - fvm::laplacian(turbulence->alphaEff(), he)
     ==
        rho*(U&g)
      + parcels.Sh(he)
      + radiation->Sh(thermo, he)
      + Qdot
      + fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();
    radiation->correct();

/*
    forAll(he,i)
    {
      fout<<p[i]<<",";
    } 
    fout<<std::endl;
    fout.close();
    */

    Info<< "T gas min/max   " << min(T).value() << ", "
        << max(T).value() << endl;
    

/*
    p=(gammastar-1)*rho*(he-estar);
    he.write();
    p.write();
    U.write();
    rho.write();
    */

}
