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
    interFoam

Group
    grpMultiphaseSolvers

Description
    Solver for two incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//define monitor time indexes
#define ifMonitor  if( runTime.timeIndex()%10== 0 )

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two incompressible, isothermal immiscible fluids"
        " using VOF phase-fraction based interface capturing.\n"
        "With optional mesh motion and mesh topology changes including"
        " adaptive re-meshing."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createAlphaFluxes.H"
    #include "initCorrectPhivd.H"
    #include "initCorrectPhicr.H"
    #include "createUfIfPresent.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "alphaCourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                mesh.update();

                if (mesh.changing())
                {
                    // Do not apply previous time-step mesh compression flux
                    // if the mesh topology changed
                    if (mesh.topoChanging())
                    {
                        talphaPhi1Corr0.clear();
                    }

                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);

                        mixture.correct();
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            U = 0*U;
            phi = 0*phi;

            // Initialize alpha for capillary relaxation step
            alpha1.oldTime() = alpha1;
            int nRelax = ::floor(runTime.deltaTValue()/deltaTCR-1e-9)+1;
            runTime.setDeltaT(runTime.deltaTValue()/nRelax,false);

            nSteps++;
            nRelaxTot += nRelax;


            for (int i=0;i<nRelax;i++)
            { 
                nRelaxSteps++;

                Info << "Time=" << runTime.timeName() << ", relaxation iteration:" << i+1 << " out of " << nRelax << endl;

                #include "alphaRelaxEqnSubCycle.H"

                mixture.correct();
                
                if (pimple.frozenFlow())
                {
                    continue;
                }

                rhoPhi = rhoPhivd + rhoPhicr;

                #include "URelaxEqn.H"

                // --- PISO loop
                while (pimple.correct())
                {
                    #include "pRelaxEqn.H"
                }

                #include "continuityErrsRelax.H"

                Info << "\n  Ucrmax = " << max(mag(Ucr)).value() << " m/s  " ;


                alpha1.oldTime() = alpha1;
                Ucr.oldTime() = Ucr;
                rho_cr.oldTime() = rho_cr;
                phicr.oldTime() = phicr;

                U += Ucr/nRelax;
                phi += phicr/nRelax;


                if (i>=relaxMin-1 && i> 0.9*(nIterLast-1) && mixture.eqnResidual() < residual) 
                {
                    Info << "\n         Time=" << runTime.timeName() << ", capillary relaxation finished in " << 1+i << " steps out of " << nRelax << endl;
                    nIterLast = i;
                    break;
                }
                if (i==nRelax-1)
                {
                    Info << "\n         Time=" << runTime.timeName() << ", capillary relaxation finished in " << 1+i << " steps out of " << nRelax << endl;
                    nIterLast = i;
                }
            }

            Info << "\n         Time=" << runTime.timeName() << ", end of capillary relaxation" << endl; 
            Info << "\n         Average number of relaxation steps="<< nRelaxSteps/nSteps << endl;
            Info << "\n         Percent of relaxation steps performed="<< nRelaxSteps/nRelaxTot*100 << "%" << endl;
            Info << "\n" << endl;

            runTime.setDeltaT(nRelax*runTime.deltaTValue(),false);

            rho_vd = rho_cr;


            if (pimple.frozenFlow())
            {
                continue;
            }

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            U += Uvd;
            phi += phivd;

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();

        //monitor average and max velocity
        ifMonitor
        {
                        Info << "\n         Umax = " << max(mag(U)).value() << " m/s  "
                        << "Uavg = " << mag(average(U)).value() << " m/s";
        }


        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
