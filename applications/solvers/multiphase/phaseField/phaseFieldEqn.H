{
    volScalarField Kc = C*scalar(0)/runTime.deltaT();

    Info<< "Solving C  RK4 Equations: ";

    for (int i=0; i<=3; i++)
    {
        Info << " " << scalar(i);
        scalar T_Multiplier = scalar(0.5) + scalar(i/2)*scalar(0.5);
        scalar K_Multiplier = scalar(1)/(scalar(3) + scalar(3)*mag(scalar(1) - scalar((i + 1)/2)));


        surfaceScalarField phiC =
            fvc::flux
            (
                phi,
                C,
                "div(phi,C)"
            );

        surfaceScalarField mob = mixture.mob()*(scalar(1)-sqr(fvc::interpolate(C)));

        surfaceScalarField difFlux = -mob*fvc::snGrad(mixture.eta())*mesh.magSf();

        volScalarField Cflux
        (
             -fvc::div(phi, C)
             -fvc::div(difFlux)
        );

        C = C.oldTime() + runTime.deltaT()*T_Multiplier*Cflux;
        C.correctBoundaryConditions();

        rhoPhi = 0.5*phiC*(rho1 - rho2) + 0.5*phi*(rho1+rho2);

        mixture.correct();

        Kc += K_Multiplier*Cflux;
    }

    C = C.oldTime() + runTime.deltaT()*Kc;
    C.correctBoundaryConditions();

    alpha1 = 0.5*(1+C);
    alpha2 = 0.5*(1-C);

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alpha1.name() << ") = " << gMin(alpha1.internalField())
        << "  Max(" << alpha1.name() << ") = " << gMax(alpha1.internalField())
        << endl;
}
