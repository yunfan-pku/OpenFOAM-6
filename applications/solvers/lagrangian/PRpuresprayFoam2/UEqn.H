// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

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
//Info<< "UEqn_#########################################4\n" << endl;

/*
foutl1<<"Up1"<<",";
forAll(p,i)
{
  foutl1<<U[i][0]<<",";
} 
foutl1<<std::endl;
foutl1<<"Up2"<<",";
forAll(p,i)
{
  foutl1<<U[i][1]<<",";
} 
foutl1<<std::endl;
foutl1<<"Up3"<<",";
forAll(p,i)
{
  foutl1<<U[i][2]<<",";
} 
foutl1<<std::endl;

foutl1<<"p"<<",";
forAll(p,i)
{
  foutl1<<p[i]<<",";
} 
foutl1<<std::endl;
tmp<volVectorField> gp=-fvc::grad(p);
foutl1<<"gp"<<",";
forAll(p,i)
{
  foutl1<<gp()[i][0]<<",";
} 
foutl1<<std::endl;
*/