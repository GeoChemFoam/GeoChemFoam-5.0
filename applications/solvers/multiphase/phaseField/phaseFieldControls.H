const dictionary& phaseFieldControls = mesh.solverDict(C.name());

label nPhaseFieldSubCycles(phaseFieldControls.get<label>("nPhaseFieldSubCycles"));
