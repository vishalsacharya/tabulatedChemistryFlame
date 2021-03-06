#include "createRDeltaT.H"

Info<< "Reading Z field\n" << endl;
volScalarField Z
(
    IOobject
    (
        "Z",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Creating field rho\n" << endl;
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(1, -3, 0, 0, 0)
);

Info<< "Creating field psi\n" << endl;
volScalarField psi
(
    IOobject
    (
        "psi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionSet(0, -2, 2, 0, 0)
);

#include "createFieldRefs.H"

#include "readingTabulation.H"
#include "updateFields.H"

#include "compressibleCreatePhi.H"

pressureControl pressureControl(p, rho, pimple.dict(), false);

mesh.setFluxRequired(p.name());

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

#include "createMRF.H"

// Check the gravity
if ( (mag(g)).value() != 0.)
{
	Info << endl;
	Info << "------------------------------------------------------------------" << endl;
	Info << " ERROR:   gravity is turned on                                    " << endl;
	Info << "          Please consider that this solver is conceived           " << endl;
	Info << "          for flows iwithout gravity.                             " << endl;
	Info << "          If you think gravity cannot be neglegted, consider the  " << endl;
	Info << "          possibility to use the                                  " << endl;
	Info << "          tabulatedChemistryFlameBuoyantPimple solver             " << endl;
	Info << "------------------------------------------------------------------" << endl;
	Info << endl;
	abort();
}
