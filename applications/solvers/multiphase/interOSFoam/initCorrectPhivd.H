tmp<volScalarField> rAU;

if (correctPhi)
{
    rAU = new volScalarField
    (
        IOobject
        (
            "rAU",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rAU", dimTime/dimDensity, 1)
    );

    #include "correctPhi.H"
}
else
{
    CorrectPhi
    (
        Uvd,
        phivd,
        pvd,
        dimensionedScalar("rAUf", dimTime/rho_vd.dimensions(), 1),
        geometricZeroField(),
        pimple
    );

    #include "continuityErrs.H"
}
