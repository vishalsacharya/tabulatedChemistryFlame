volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

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

if (pimple.transonic())
{
    surfaceScalarField phid
    (
        "phid",
        (fvc::interpolate(psi)/fvc::interpolate(rho))*phiHbyA
    );
    phiHbyA -= fvc::interpolate(p)*phid;

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvc::div(phiHbyA)
          + fvm::div(phid, p)
          - fvm::laplacian(rhorAUf, p)
         ==
            fvOptions(psi, p, rho.name())
        );

        // Relax the pressure equation to ensure diagonal-dominance
        pEqn.relax();

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
        }
    }
}
else
{
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvc::div(phiHbyA)
          - fvm::laplacian(rhorAUf, p)
         ==
            fvOptions(psi, p, rho.name())
        );

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
        }
    }
}

#include "rhoEqn.H"
{
    dimensionedScalar totalMass = fvc::domainIntegrate(rho);

    scalar sumLocalContErr =
        (fvc::domainIntegrate(mag(rho - Rho(p,psi)))/totalMass).value();

    scalar globalContErr =
        (fvc::domainIntegrate(rho - Rho(p,psi))/totalMass).value();

    cumulativeContErr += globalContErr;

    Info << "time step continuity errors : sum local = " << sumLocalContErr
         << ", global = " << globalContErr
         << ", cumulative = " << cumulativeContErr
         << endl;
}

// Explicitly relax pressure for momentum corrector
p.relax();

U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);
K = 0.5*magSqr(U);

pressureControl.limit(p);
p.correctBoundaryConditions();
rho = Rho(p,psi);

if (!pimple.transonic())
{
    rho.relax();
}
Info<< "rho max/min : " << max(rho).value() << " " << min(rho).value() << endl;

// Be careful
// In case of energy balance equation here you can calculate the DpDt contribution
// The present solver does not need it
