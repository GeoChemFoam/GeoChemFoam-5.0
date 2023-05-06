if (nPhaseFieldSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();
    surfaceScalarField rhoPhiSum
    (
        IOobject
        (
            "rhoPhiSum",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(rhoPhi.dimensions(), Zero)
    );

    tmp<volScalarField> trSubDeltaT;

    if (LTS)
    {
        trSubDeltaT =
            fv::localEulerDdt::localRSubDeltaT(mesh, nPhaseFieldSubCycles);
    }

    for
    (
        subCycle<volScalarField> phaseFieldSubCycle(C, nPhaseFieldSubCycles);
        !(++phaseFieldSubCycle).end();
    )
    {
        #include "phaseFieldEqn.H"
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
    }

    rhoPhi = rhoPhiSum;
}
else
{
    #include "phaseFieldEqn.H"
}

rho == alpha1*rho1 + alpha2*rho2;
