//- Read options
Info<< "Reading Solver options\n" << endl;
IOdictionary solverOptions
(
	IOobject
	(
		"solverOptions",
		U.time().constant(),
		U.db(),
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	)
);

//- Tabulated chemistry
OpenSMOKE::NonAdiabaticFlameletLibraryReader tabulatedChemistry;

//- Importing tabulated chemistry
{
	Foam::string library_path = solverOptions.lookup("libraryPath");
	boost::filesystem::path library_folder = library_path;
	tabulatedChemistry.SetLibraryFolder(library_folder);
}

//- Set minimum value of scalar dissipation rate
dimensionedScalar chiMinimum( "chiMinimum", dimensionSet(0,0,-1,0,0,0,0), readScalar(solverOptions.lookup("chiMinimum")) );

//- List of available species
{
	std::vector<std::string> list_of_species;
	tabulatedChemistry.GetSpeciesAvailable(list_of_species);
	tabulatedChemistry.SetSpeciesToExtract(list_of_species);
}

//- Set adiabatic mode
tabulatedChemistry.SetAdiabaticMode();

//- Read all flamelet
tabulatedChemistry.Read();
tabulatedChemistry.Summary();

//- Names of species in tabulation
Info << "List of available species..." << endl;
for (int j=0;j<tabulatedChemistry.number_of_species();j++)
	Info << j+1 << " " << tabulatedChemistry.species()[j] << endl;


