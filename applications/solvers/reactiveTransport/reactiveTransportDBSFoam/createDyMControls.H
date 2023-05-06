#include "createControl.H"
#include "createTimeControls.H"

word VoS(mesh.solutionDict().subDict("STEADYSTATE").lookup("VoS"));

scalar psiCoeff(mesh.solutionDict().subDict("STEADYSTATE").lookupOrDefault("psiCoeff",4.0));

bool adaptPsiCoeff(mesh.solutionDict().subDict("STEADYSTATE").lookupOrDefault("adaptPsiCoeff",false));
