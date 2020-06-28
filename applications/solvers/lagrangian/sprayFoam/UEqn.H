// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);
Info<<"!!!!!!!!!!!!!!!!!!!!!!!!!"<<Y[0][0]<<" "<<p[0]<<" "<<T[0]<<endl;
tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(rho, U) + fvm::div(phi, U)
  + MRF.DDt(rho, U)
  + turbulence->divDevRhoReff(U)
 ==
    rho()*g
  + parcels.SU(U)
  + fvOptions(rho, U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
    K = 0.5*magSqr(U);
}
