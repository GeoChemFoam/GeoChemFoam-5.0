    // Momentum predictor

    MRF.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> tUEqn
    (
        (1/eps/eps)*fvm::div(phi, U)
      + (1/eps/eps)*MRF.DDt(U)
      + turbulence->divDevReff(U)
      + fvm::Sp(nu*Kinv,U)
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (simple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));

        fvOptions.correct(U);
    }
