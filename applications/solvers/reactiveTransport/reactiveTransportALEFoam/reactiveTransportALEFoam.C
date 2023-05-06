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
    reactiveTransportALEFoam

Description
    Solves reactive transport equation with ALE mesh for multi-species flow

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

    simpleControl simple(mesh);

    #include "createFields.H"
    
    const speciesTable& solutionSpecies = speciesMixture.species();

    #include "createTimeControls.H"

    turbulence->validate();
    #include "meshCourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "meshCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        mesh.update();

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
        if (!steadyState.criteriaSatisfied())
        {
            U == U.oldTime();
            p == p.oldTime();
            forAll(solutionSpecies,i)
            {
                volScalarField Yi = speciesMixture.Y(i);
                Yi == Yi.oldTime(); 
            }
        }

        if (runTime.timeIndex() > remeshFrequency-1 )
        {
            runTime.writeNow();
            return 0;
        }

		if (runTime.timeIndex() % checkFrequency == 0)
		{
			if (mesh.checkMesh(false))
			{
				runTime.writeNow();
				return 0;
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
