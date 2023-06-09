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

#include "solutionSurfaceMultiComponentMixture.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solutionSurfaceMultiComponentMixture::solutionSurfaceMultiComponentMixture
(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const fvMesh& mesh
)
:
    basicMultiComponentMixture(thermoDict, specieNames, mesh),
    surfaceMixture_(thermoDict,thermoDict.subDict("surfaceSpecies").toc(),mesh),
    surfaceSpecies_(surfaceMixture_.species()),
	surfaceMasters_(thermoDict.subDict("surfaceMasters").toc())
{
}

Foam::solutionSurfaceMultiComponentMixture::solutionSurfaceMultiComponentMixture
(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const fvMesh& mesh,
    const word& phaseName
)
:
    basicMultiComponentMixture(thermoDict, specieNames, mesh, phaseName),
    surfaceMixture_(thermoDict,thermoDict.subDict("surfaceSpecies").toc(),mesh),
    surfaceSpecies_(surfaceMixture_.species()),
	surfaceMasters_(thermoDict.subDict("surfaceMasters").toc())
{
}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ************************************************************************* //
