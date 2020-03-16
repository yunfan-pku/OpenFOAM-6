// Solve the Momentum equation
//Info<< "UEqn_#########################################1\n" << endl;
MRF.correctBoundaryVelocity(U);
//Info<< "UEqn_#########################################2\n" << endl;
tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(rho, U) + fvm::div(phi, U)
  + MRF.DDt(rho, U)
 // + turbulence->divDevRhoReff(U)
 ==
    rho()*g
  + parcels.SU(U)
  + fvOptions(rho, U)
);
/*
foutl1<<"rho"<<",";
forAll(p,i)
{
  foutl1<<rho[i]<<",";
} 
foutl1<<std::endl;

foutl1<<"U1"<<",";
forAll(p,i)
{
  foutl1<<U[i][0]<<",";
} 
foutl1<<std::endl;
foutl1<<"U2"<<",";
forAll(p,i)
{
  foutl1<<U[i][1]<<",";
} 
foutl1<<std::endl;
foutl1<<"U3"<<",";
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
*/
//Info<< "UEqn_#########################################3\n" << endl;
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