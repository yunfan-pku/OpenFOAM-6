rho = thermo.rho();
rho = max(rho, rhoMin);
rho = min(rho, rhoMax);
rho.relax();
// PY 0906
// Thermodynamic density needs to be updated by psi*d(p) after the
// pressure solution
const volScalarField psip0(psi*p);
const volScalarField p0(p);

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

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
    }
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
        /*fvScalarMatrix pEqn
        (
           fvm::ddt(psi, p)
         + fvc::div(phiHbyA)
         - fvm::laplacian(rhorAUf, p)
         //== parcels.Srho() + fvOptions(psi, p, rho.name())
        );*/
        
        fvScalarMatrix pEqn
        (
           fvc::ddt(rho) 
        +  psi*correction(fvm::ddt(p))
        +  fvc::div(phiHbyA)
        -  fvm::laplacian(rhorAUf, p)
        == 
           //parcels.Srho() 
           fvOptions(psi, p, rho.name())
        );
        
        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));
        
        const volScalarField detp = p0 - p;
        scalar maxdep = fabs(max(detp).value()); 
        if (maxdep > 1.0e+04)
        {
            //Info<< "p2= " << maxdep << ", "<< max(p).value() << endl;
            /*forAll(detp, celli)
            {
                if (detp[celli] > maxdep)
                p[celli] = p0[celli] + (p[celli] - p0[celli])*0.1;
            }*/
            //const volScalarField rho2 = 0.5*thermo.rho();
        }

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
        }
    }
}

// PY 0906
// Thermodynamic density update
//thermo.correctRho(psi*p - psip0);

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();
// PY 0906
// Correct pressure and density
pressureControl.limit(p);
p.correctBoundaryConditions();

U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);

K = 0.5*magSqr(U);

//Info<< "p_gas_min/max = " << min(p).value() << ", "<< max(p).value() << endl;

thermo.correct();

if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);
}

