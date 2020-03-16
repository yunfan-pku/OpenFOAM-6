{
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      /*+ (
            he.name() == "he"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )*/
      - fvm::laplacian(turbulence->alphaEff(), he)
     ==
        //parcels.Sh(he)
        rho*(U&g)
      //+ radiation->Sh(thermo, he)
      + Qdot
      + fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();
    
    fvOptions.correct(he);

    //thermo.correct();
    //Info<< "T_gas_min/max = " << min(T).value() << ", "<< max(T).value() << endl;
    //Info<< "he_gas_min/max = " << min(he).value() << ", "<< max(he).value() << endl;
    //radiation->correct();

}
