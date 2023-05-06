{
    if (correctPhi)
    {
        rAU.ref() = 1.0/UEqn.A();
    }
    else
    {
        rAU = 1.0/UEqn.A();
    }

    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU()));
    volVectorField HbyA(constrainHbyA(rAU()*UEqn.H(), Uvd, pvd));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + MRF.zeroFilter(fvc::interpolate(rho_vd*rAU())*fvc::ddtCorr(Uvd, phivd, Uf))
    );
    MRF.makeRelative(phiHbyA);

    if (p_rgh.needReference())
    {
        fvc::makeRelative(phiHbyA, Uvd);
        adjustPhi(phiHbyA, Uvd, pvd);
        fvc::makeAbsolute(phiHbyA, Uvd);
    }

    surfaceScalarField phig
    (
        (
          //  mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rho_vd)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(pvd, Uvd, phiHbyA, rAUf, MRF);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, pvd) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phivd = phiHbyA - p_rghEqn.flux();

            pvd.relax();

            Uvd = HbyA + rAU()*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
            Uvd.correctBoundaryConditions();
            fvOptions.correct(Uvd);
        }
    }

    #include "continuityErrs.H"

    // Correct Uf if the mesh is moving
    fvc::correctUf(Uf, Uvd, phivd);

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phivd, Uvd);

    p_rgh = pvd + pcr;
    p == p_rgh + mixture.pc() + rho_vd*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - mixture.pc() - rho_vd*gh;
    }

    if (!correctPhi)
    {
        rAU.clear();
    }
}
