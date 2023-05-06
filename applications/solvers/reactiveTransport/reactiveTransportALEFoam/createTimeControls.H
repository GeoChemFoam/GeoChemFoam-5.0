
bool adjustTimeStep =
    runTime.controlDict().getOrDefault("adjustTimeStep", false);

scalar maxMeshCo =
    runTime.controlDict().getOrDefault<scalar>("maxMeshCo", 1);

scalar maxDeltaT =
    runTime.controlDict().getOrDefault<scalar>("maxDeltaT", GREAT);

// ************************************************************************* //
