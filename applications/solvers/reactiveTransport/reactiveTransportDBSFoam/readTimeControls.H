
adjustTimeStep =
    runTime.controlDict().getOrDefault("adjustTimeStep", false);

maxDeltaEps =
    runTime.controlDict().getOrDefault<scalar>("maxDeltaEps", 0.1);

maxDeltaT =
    runTime.controlDict().getOrDefault<scalar>("maxDeltaT", GREAT);

// ************************************************************************* //
