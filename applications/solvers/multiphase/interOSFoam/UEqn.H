    MRF.correctBoundaryVelocity(Uvd);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho_vd, Uvd) + fvm::div(rhoPhi, Uvd,"div(rho*phi,U)")
      + MRF.DDt(rho_vd, Uvd)
      + turbulence->divDevRhoReff(rho_vd, Uvd)
     ==
        fvOptions(rho_vd, Uvd)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                  //  mixture.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho_vd)
                  - fvc::snGrad(pvd)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(Uvd);
    }
