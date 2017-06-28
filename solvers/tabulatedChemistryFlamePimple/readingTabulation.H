//- Tabulated chemistry
OpenSMOKE::NonAdiabaticFlameletLibraryReader tabulatedChemistry;

boost::filesystem::path library_folder = "/home/chimica2/cuoci/OpenFOAM/cuoci-3.0.x/mySolvers/tabulatedChemistryFlame/run/tutorials/coflowFlame/library/OutputXML";
tabulatedChemistry.SetLibraryFolder(library_folder);
//tabulatedChemistry.SetSpeciesToExtract(list_of_species);

//- Set adiabatic mode
tabulatedChemistry.SetAdiabaticMode();

//- Read all flamelet
tabulatedChemistry.Read();
tabulatedChemistry.Summary();
