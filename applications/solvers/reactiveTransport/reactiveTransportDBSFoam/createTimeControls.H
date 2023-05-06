
bool adjustTimeStep =
    runTime.controlDict().getOrDefault("adjustTimeStep", false);

scalar maxDeltaEps =
    runTime.controlDict().getOrDefault<scalar>("maxDeltaEps", 0.1);

scalar maxDeltaT =
    runTime.controlDict().getOrDefault<scalar>("maxDeltaT", GREAT);

// ************************************************************************* //
