    // Momentum predictor

    MRF.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> tUEqn
    (
        (1/eps/eps)*fvm::div(phi, U)
      + 1/eps*MRF.DDt(U)
      - fvm::laplacian(nu,U)//+ turbulence->divDevReff(U)
      + fvm::Sp(nu*Kinv,U)
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (steadyState.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));

        fvOptions.correct(U);
    }
