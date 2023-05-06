CorrectPhi
(
    U,
    phi,
    p_rgh,
    surfaceScalarField("rAUf", fvc::interpolate(rAU())),
    -speciesMixture.Mflux()*(1/rho1-1/rho2),
    pimple
);

#include "continuityErrs.H"
