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

\*---------------------------------------------------------------------------*/

#include "steadyStateControl.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(steadyStateControl, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::steadyStateControl::read()
{
    solutionControl::read(true);
    const dictionary& steadyStateDict =dict();
    nCorr_ = steadyStateDict.lookupOrDefault<label>("nCorrectors",1000);
    return true;
}


bool Foam::steadyStateControl::criteriaSatisfied()
{
    if (residualControl_.empty())
    {
        return false;
    }

    bool achieved = true;
    bool checked = false;    // safety that some checks were indeed performed

    const dictionary& solverDict = mesh_.solverPerformanceDict();
    forAllConstIters(solverDict, iter)
    {
        const entry& solverPerfDictEntry = *iter;

        const word& fieldName = solverPerfDictEntry.keyword();
        label fieldi = applyToField(fieldName);
       
        if (fieldi == -1)
        {
            forAll(residualControl_, i)
            {
                if (residualControl_[i].name == fieldName)
                {
                    fieldi = i;
                }
            }
        }
        if (fieldi != -1)
        {
            Pair<scalar> residuals = maxResidual(solverPerfDictEntry);
            checked = true;

            const bool absCheck =
                (residuals.last() < residualControl_[fieldi].absTol);

            achieved = achieved && absCheck;

            if (debug)
            {
                Info<< algorithmName_ << " solution statistics:" << endl;

                Info<< "    " << fieldName << ": tolerance = "
                    << residuals.last()
                    << " (" << residualControl_[fieldi].absTol << ")"
                    << endl;
            }
        }
    }

    return checked && achieved;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::steadyStateControl::steadyStateControl
(
    fvMesh& mesh,
    const word& dictName,
    const bool verbose
)
:
    solutionControl(mesh, dictName),
    initialised_(false),
    iter_counter(0)
{
    read();

    if (verbose)
    {
        Info<< nl << algorithmName_;

        if (residualControl_.empty())
        {
            const scalar duration =
                mesh_.time().endTime().value()
              - mesh_.time().startTime().value();

            Info<< ": no convergence criteria found. "
                << "Calculations will run for " << duration << " steps."
                << nl;
        }
        else
        {
            Info<< ": convergence criteria" << nl;
            for (const fieldData& ctrl : residualControl_)
            {
                Info<< "    field " << ctrl.name << token::TAB
                    << " tolerance " << ctrl.absTol
                    << nl;
            }
        }
        Info<< endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::steadyStateControl::loop()
{
    solutionControl::setFirstIterFlag(true, true);

    read();

    if (initialised_)
    {
        iter_counter++;
        if (iter_counter==nCorr_ + 1)
        {
            Info<< nl
            << algorithmName_
            << " solution not converged within "
            << nCorr_ << " iterations" << nl << endl;

            // Set to finalise calculation
            return false;
        }
        if (criteriaSatisfied())
        {
            Info<< nl
            << algorithmName_
            << " solution converged in "
            << iter_counter << " iterations" << nl << endl;
            // Set to finalise calculation
            return false;
        }

    }
    else
    {
        initialised_ = true;
        storePrevIterFields();
    }

    return true;
}


// ************************************************************************* //
