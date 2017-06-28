{
	const double small_eps = 1.e-6;
	const double small_chi_st = 1.e-8;
	std::vector<double> extracted(7);

	double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

	Info<< "Update fields... " ;

	// Adiabatic library
	const double defect = 0.;

	// Internal fields
	{
		const scalarField& ZCells = Z.internalField();
		const scalarField& pCells = p.internalField(); 

		scalarField& TCells      = T.ref();
		scalarField& rhoCells    = rho.ref();
		scalarField& psiCells    = psi.ref();
		scalarField& muCells     = mu.ref();
		scalarField& alphaCells  = alpha.ref();

		forAll(ZCells, celli)
		{
			//double max_chi = max(small_chi_st,chi_st[celli]);
			const double max_chi = 1.;	//TODO

			//- Pure oxidizer
			if (ZCells[celli]<=small_eps)
			{
				tabulatedChemistry.GetMeanValues(0., 0., max_chi, defect, extracted);
			}

			//- Pure fuel
			else if (ZCells[celli]>=(1.-small_eps))
			{
				tabulatedChemistry.GetMeanValues(1., 0., max_chi, defect, extracted);
			}

			//- Mixture
			else
			{
				tabulatedChemistry.GetMeanValues(ZCells[celli], 0., max_chi, defect, extracted);
			}

			TCells[celli] 		= extracted[1];
			rhoCells[celli]		= extracted[2];
			muCells[celli]		= extracted[4];
			alphaCells[celli]	= extracted[5];
			psiCells[celli]  	= rhoCells[celli]/pCells[celli];
		}
    	}

	// Boundaries
	forAll(Z.boundaryField(), patchi)
	{
		fvPatchScalarField& pZ     = Z.boundaryFieldRef()[patchi];
		fvPatchScalarField& pp     = p.boundaryFieldRef()[patchi];

		fvPatchScalarField& prho   = rho.boundaryFieldRef()[patchi];
		fvPatchScalarField& pT     = T.boundaryFieldRef()[patchi];
		fvPatchScalarField& ppsi   = psi.boundaryFieldRef()[patchi];
		fvPatchScalarField& pmu    = mu.boundaryFieldRef()[patchi];
		fvPatchScalarField& palpha = alpha.boundaryFieldRef()[patchi];
		
		forAll(pZ, facei)
		{
			//double max_chi = max(small_chi_st, pchi_st[facei]);
			const double max_chi = 1.;	//TODO

			//- Pure oxidizer
			if (pZ[facei]<=0.)
			{
				tabulatedChemistry.GetMeanValues(0., 0., max_chi, defect, extracted);
			}

			//- Pure fuel
			else if (pZ[facei]>=1.)
			{
				tabulatedChemistry.GetMeanValues(1., 0., max_chi, defect, extracted);
			}

			//- Mixture
			else
			{
				tabulatedChemistry.GetMeanValues(pZ[facei],0., max_chi,defect, extracted);
			}

			pT[facei] 	= extracted[1];
			prho[facei]	= extracted[2];
			pmu[facei]	= extracted[4];
			palpha[facei]	= extracted[5];
			ppsi[facei]  	= prho[facei]/pp[facei];
		}
	}
	
	double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

	Info << "done in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per cell)" << endl;
}
