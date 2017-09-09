if (runTime.outputTime() == true)
{
	const double small_eps = 1.e-6;
	const double small_chi_st = 1.e-8;
	std::vector<double> extracted(tabulatedChemistry.number_of_species()+1);

	double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

	Info<< "Reconstruct fields... " ;

	// Update scalar dissipation rate
	chi = 2.0 * alpha * magSqr(fvc::grad(Z));
	chi.correctBoundaryConditions();

	//- Mass fraction fields
	PtrList<volScalarField> Y(tabulatedChemistry.number_of_species());
	for(int j=0;j<tabulatedChemistry.number_of_species();j++)
	{
		Y.set
		(
		        j,
		        new volScalarField
		        (
		                IOobject
		                (
				        tabulatedChemistry.species()[j],
				        mesh.time().timeName(),
				        mesh,
				        IOobject::NO_READ,
				        IOobject::AUTO_WRITE
		                ),
		                mesh,
		                dimensionSet(0, 0, 0, 0, 0)
		        )
		);
	}

	// Adiabatic library
	const double defect = 0.;

	// Internal fields
	{
		const scalarField& ZCells   = Z.internalField();
		const scalarField& pCells   = p.internalField(); 
		const scalarField& chiCells = chi.internalField(); 

		forAll(ZCells, celli)
		{
			double max_chi = max(small_chi_st, chiCells[celli]);

			//- Pure oxidizer
			if (ZCells[celli]<=small_eps)
			{
				tabulatedChemistry.ExtractMeanValues(0., 0., max_chi, defect, extracted);
			}

			//- Pure fuel
			else if (ZCells[celli]>=(1.-small_eps))
			{
				tabulatedChemistry.ExtractMeanValues(1., 0., max_chi, defect, extracted);
			}

			//- Mixture
			else
			{
				tabulatedChemistry.ExtractMeanValues(ZCells[celli], 0., max_chi, defect, extracted);
			}
			
			for(int j=0;j<tabulatedChemistry.number_of_species();j++)
				Y[j].ref()[celli] = extracted[j+1];
		}
    	}

	// Boundaries
	forAll(Z.boundaryField(), patchi)
	{
		fvPatchScalarField& pZ     = Z.boundaryFieldRef()[patchi];
		fvPatchScalarField& pp     = p.boundaryFieldRef()[patchi];
		fvPatchScalarField& pchi   = chi.boundaryFieldRef()[patchi];

		fvPatchScalarField& prho   = rho.boundaryFieldRef()[patchi];
		fvPatchScalarField& pT     = T.boundaryFieldRef()[patchi];
		fvPatchScalarField& ppsi   = psi.boundaryFieldRef()[patchi];
		fvPatchScalarField& pmu    = mu.boundaryFieldRef()[patchi];
		fvPatchScalarField& palpha = alpha.boundaryFieldRef()[patchi];
		
		forAll(pZ, facei)
		{
			double max_chi = max(small_chi_st, pchi[facei]);

			//- Pure oxidizer
			if (pZ[facei]<=0.)
			{
				tabulatedChemistry.ExtractMeanValues(0., 0., max_chi, defect, extracted);
			}

			//- Pure fuel
			else if (pZ[facei]>=1.)
			{
				tabulatedChemistry.ExtractMeanValues(1., 0., max_chi, defect, extracted);
			}

			//- Mixture
			else
			{
				tabulatedChemistry.ExtractMeanValues(pZ[facei],0., max_chi,defect, extracted);
			}

			for(int j=0;j<tabulatedChemistry.number_of_species();j++)
				Y[j].boundaryFieldRef()[patchi][facei] = extracted[j+1];
		}
	}

	for(int j=0;j<tabulatedChemistry.number_of_species();j++)
		Y[j].write();
	
	double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

	Info << "done in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per cell)" << endl;
}
