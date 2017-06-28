Info<< "Creating field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(0, 0, 0, 0, 1)
);

//- Allocating memory for dynamic viscosity [kg/m/s]
Info<< "Creating field mu\n" << endl;
volScalarField mu
(
    IOobject
    (
        "gas_mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(1, -1, -1, 0, 0)
);

//- Allocating memory for thermal diffusion [m2/s]
Info<< "Creating field alpha\n" << endl;
volScalarField alpha
(
    IOobject
    (
        "gas_alpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(0, 2, -1, 0, 0)
);
