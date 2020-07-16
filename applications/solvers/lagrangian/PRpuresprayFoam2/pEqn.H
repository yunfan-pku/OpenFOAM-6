/*
Info<< "rho gas min/max!!!   " << min(rho).value() << ", "
        << max(rho).value() << endl;
Info<< "p max/min!!! : " << max(p).value()
    << " " << min(p).value() << endl;
Info<< "psi max/min!!! : " << max(psi).value()
    << " " << min(psi).value() << endl;
    */
rho = thermo.rho();
//Info<< "pEqn_rho#########################################3\n" << endl;
//which is from psiThermo.C
Info<< "rho gas min/max   " << min(rho).value() << ", "
        << max(rho).value() << endl;
Info<< "p max/min01 : " << max(p).value()
    << " " << min(p).value() << endl;
//thermo.correct();
rho = max(rho, rhoMin);
rho = min(rho, rhoMax);
rho.relax();


Info<< "rho gas min/max   " << min(rho).value() << ", "
        << max(rho).value() << endl;


volScalarField rAU(1.0/UEqn.A());
/*
forAll(p,i)
{
fout<<rAU[i]<<",";
} 
fout<<std::endl;
*/
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}
//////////////////////////
if (pimple.transonic())
{
    
    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(psi)
       *(
            fvc::flux(HbyA)
          + MRF.zeroFilter
            (
                rhorAUf*fvc::ddtCorr(rho, U, phi)/fvc::interpolate(rho)
            )
        )
    );

    MRF.makeRelative(fvc::interpolate(psi), phid);
    
    while (pimple.correctNonOrthogonal())
    {
        //csvwrite(filename, p, (runTime.timeName()+s+std::string("p_b_pEqn")).data());
        //csvwrite(filename, psi, (runTime.timeName()+s+std::string("psi_b_pEqn")).data());
        //csvwrite(filename, phid, (runTime.timeName()+s+std::string("phid_b_pEqn")).data());
        //csvwrite(filename, rhorAUf, (runTime.timeName()+s+std::string("rhorAUf_b_pEqn")).data());
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvm::div(phid, p)
          - fvm::laplacian(rhorAUf, p)
         ==
            parcels.Srho()
          + fvOptions(psi, p, rho.name())
        );

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi == pEqn.flux();
        }
        //csvwrite(filename, p, (runTime.timeName()+s+std::string("p_a_pEqn")).data());
    }
    /////////////
}
else
{
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
            fvc::flux(rho*HbyA)
          + rhorAUf*fvc::ddtCorr(rho, U, phi)
        )
    );

    MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, rho, U, phiHbyA, rhorAUf, MRF);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvc::div(phiHbyA)
          - fvm::laplacian(rhorAUf, p)
         ==
            parcels.Srho()
          + fvOptions(psi, p, rho.name())
        );

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
        }
    }
}
/////////////////////////////
#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"
/////////////////////////////////////
//p=(gammastar-1)*rho*(he-estar);//XYY

// Explicitly relax pressure for momentum corrector
p.relax();


/////////////////////////////////////////////////////////////
// Recalculate density from the relaxed pressure
rho = thermo.rho();
//Info<< "Solving_dpdt##########"<<rho<<"##############################3\n" << endl; 
rho = max(rho, rhoMin);
rho = min(rho, rhoMax);
rho.relax();
Info<< "rho max/min : " << max(rho).value()
    << " " << min(rho).value() << endl;

U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);
K = 0.5*magSqr(U);

if (thermo.dpdt())//dpdt is in basicThermo.H
{
    dpdt = fvc::ddt(p);
    //Info<< "Solving_dpdt########################################3\n" << endl; 
}

