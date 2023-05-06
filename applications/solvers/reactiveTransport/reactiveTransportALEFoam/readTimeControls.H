
adjustTimeStep =
    runTime.controlDict().getOrDefault("adjustTimeStep", false);

maxMeshCo =
    runTime.controlDict().getOrDefault<scalar>("maxMeshCo", 1);

maxDeltaT =
    runTime.controlDict().getOrDefault<scalar>("maxDeltaT", GREAT);

label checkFrequency =          runTime.controlDict().lookupOrDefault<label>("checkFrequency", 1.0);

label remeshFrequency =         runTime.controlDict().lookupOrDefault<label>("remeshFrequency", 10000);

// ************************************************************************* //
