/*---------------------------------------------------------------------------*\

License
    This file is part of GeoChemFoam, an Open source software using OpenFOAM
    for multiphase multicomponent reactive transport simulation in pore-scale
    geological domain.

    GeoChemFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version. See <http://www.gnu.org/licenses/>.

    The code was developed by Dr Julien Maes as part of his research work for
    the GeoChemFoam Group at Heriot-Watt University. Please visit our
    website for more information <https://github.com/GeoChemFoam>.

Application
    reactiveTransportDBSFoam

Description
    Solves reactive transport equation with microcontinuum DBS for multi-species flow

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "reactiveMixture.H"
#include "multiComponentTransportMixture.H"
#include "steadyStateControl.H"
#include "dynamicFvMesh.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Steady-state solver for reactive transport with dissolving solid surface."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"

    #include "initContinuityErrs.H"

    #include "createDyMControls.H"
    #include "createFields.H"

    const speciesTable& solutionSpecies = speciesMixture.species();
    const wordList& kineticPhases = speciesMixture.kineticPhases();
    const wordList& kineticPhaseReactions = speciesMixture.kineticPhaseReactions();

    turbulence->validate();
    #include "deltaEpsMax.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "deltaEpsMax.H"
        
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Do any mesh changes
        mesh.controlledUpdate();

        if (mesh.changing())
        {
            MRF.update();

            phi = mesh.Sf() & fvc::interpolate(U);
        }

        #include "epsEqn.H"

        //permeability
        Kinv = kf*pow(1-eps,2)/pow(eps,3);

        volVectorField gradEps = fvc::grad(eps);
        surfaceVectorField gradEpsf = fvc::interpolate(gradEps);
        surfaceVectorField nEpsv = -gradEpsf/(mag(gradEpsf) + deltaN);
        nEpsf = nEpsv & mesh.Sf();

        volScalarField a = mag(fvc::grad(eps));

        scalar lambda = psiCoeff;

        if (VoS=="VoS-psi")
        {
            scalar As = a.weightedAverage(mesh.V()).value();

            a = a*(1-eps)*(1e-3+eps);

            if (adaptPsiCoeff) lambda = As/a.weightedAverage(mesh.V()).value();

            a = lambda*a;


            Info << "psiCoeff=" << lambda << endl;
       }
    	
        steadyStateControl steadyState(mesh);
        while (steadyState.loop())
        {
            // --- Pressure-velocity SIMPLE corrector
            {
                #include "UEqn.H"
                #include "pEqn.H"
            }

            laminarTransport.correct();
            turbulence->correct();

            // Concentration solver
            {
                #include "YiEqn.H"
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

   	runTime.writeNow();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
