{
	word alpharScheme("div(phirb,alpha)");
	word YiScheme("div(phi,Yi)");

    // Standard face-flux compression coefficient
    surfaceScalarField phic(speciesMixture.cYi()*mag(phi/mesh.magSf()));

    surfaceScalarField phir(phic*mixture.nHatf());

    const speciesTable& solutionSpecies = speciesMixture.species();

    forAll(solutionSpecies, i)
	{
        volScalarField& Yi = speciesMixture.Y(i);

        scalar maxYi = max(gMax(Yi),gMax(Yi.boundaryField()))+1e-30;

        Yi.oldTime() == Yi.oldTime()/maxYi;
        Yi == Yi/maxYi;

		surfaceScalarField phiComp = fvc::flux
        (
            -fvc::flux(-phir, alpha2, alpharScheme),
            alpha1,
            alpharScheme
        );

        tmp<surfaceScalarField> tYiPhi1Un
        (
            fvc::flux
            (
                phi,
                Yi,
                YiScheme
            )
		+   phiComp*speciesMixture.compressionCoeff(i)
        );

        {
            surfaceScalarField YiPhi10 = tYiPhi1Un;

            MULES::explicitSolve
            (
                geometricOneField(),
                Yi,
                phi,
                YiPhi10,
                zeroField(),
                zeroField(),
                oneField(),
                zeroField()
            );
        }

        Yi.oldTime() == Yi.oldTime()*maxYi;
        Yi == Yi*maxYi;
    }
 }
