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

#include "readhRef.H"
#include "gh.H"
#if STEADYSTATE != 1
	#include "readpRef.H"
#endif

volScalarField p_rgh
(
    IOobject
    (
	"p_rgh",
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

mesh.setFluxRequired(p_rgh.name());

// Force p_rgh to be consistent with p
p_rgh = p - rho*gh;

#if STEADYSTATE == 1

	label pRefCell = 0;
	scalar pRefValue = 0.0;

	setRefCell
	(
		p,
		p_rgh,
		simple.dict(),
		pRefCell,
		pRefValue
	);

#endif

#include "createMRF.H"

// Check the gravity
if ( (mag(g)).value() == 0.)
{
	Info << endl;
	Info << "------------------------------------------------------------------" << endl;
	Info << " ERROR:   gravity is turned off                                   " << endl;
	Info << "          Please consider that this solver is conceived           " << endl;
	Info << "          for flows in which gravity play a major role.           " << endl;
	Info << "          If you think gravity plays a negligible role, consider  " << endl;
	Info << "          the possibility to use the                              " << endl;
	Info << "          tabulatedChemistryFlamePimple solver                    " << endl;
	Info << "------------------------------------------------------------------" << endl;
	Info << endl;
	abort();
}
