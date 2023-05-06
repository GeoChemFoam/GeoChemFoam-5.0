{
    if (correctPhi)
    {
        rAU.ref() = 1.0/URelaxEqn.A();
    }
    else
    {
        rAU = 1.0/URelaxEqn.A();
    }

    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU()));
    volVectorField HbyA(constrainHbyA(rAU()*URelaxEqn.H(), Ucr, pcr));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + MRF.zeroFilter(fvc::interpolate(rho_cr*rAU())*fvc::ddtCorr(Ucr, phicr, Uf))
    );
    MRF.makeRelative(phiHbyA);

    if (p_rgh.needReference())
    {
        fvc::makeRelative(phiHbyA, Ucr);
        adjustPhi(phiHbyA, Ucr, pcr);
        fvc::makeAbsolute(phiHbyA, Ucr);
    }

    surfaceScalarField phig
    (
        rAUf*mixture.phic()
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(pcr, Ucr, phiHbyA, rAUf, MRF);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, pcr) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phicr = phiHbyA - p_rghEqn.flux();

            pcr.relax();

            Ucr = HbyA + rAU()*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
            Ucr.correctBoundaryConditions();
            fvOptions.correct(Ucr);
        }
    }

    #include "continuityErrsRelax.H"

    // Correct Uf if the mesh is moving
    fvc::correctUf(Uf, Ucr, phicr);

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phicr, Ucr);

    if (!correctPhi)
    {
        rAU.clear();
    }
}
