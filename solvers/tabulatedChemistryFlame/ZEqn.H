//-  Solving equatio of mixture fraction
fvScalarMatrix ZEqn
(
	fvm::ddt(rho,Z)
      + fvm::div(phi,Z) 
        == 
	fvm::laplacian(rho*alpha,Z) 
      + fvOptions(rho,Z)
);

// Solve
ZEqn.relax();
fvOptions.constrain(ZEqn);
ZEqn.solve();
fvOptions.correct(Z);

