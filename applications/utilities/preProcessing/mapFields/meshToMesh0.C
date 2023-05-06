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

#include "meshToMesh0.H"
#include "processorFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(meshToMesh0, 0);
}

const Foam::scalar Foam::meshToMesh0::directHitTol = 1e-10;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshToMesh0::meshToMesh0
(
    const fvMesh& meshFrom,
    const fvMesh& meshTo,
    const HashTable<word>& patchMap,
    const wordList& cuttingPatchNames
)
:
    fromMesh_(meshFrom),
    toMesh_(meshTo),
    patchMap_(patchMap),
    cellAddressing_(toMesh_.nCells()),
    boundaryAddressing_(toMesh_.boundaryMesh().size()),
    inverseDistanceWeightsPtr_(nullptr),
    inverseVolumeWeightsPtr_(nullptr),
    cellToCellAddressingPtr_(nullptr),
    V_(0)
{
    forAll(fromMesh_.boundaryMesh(), patchi)
    {
        fromMeshPatches_.insert
        (
            fromMesh_.boundaryMesh()[patchi].name(),
            patchi
        );
    }

    forAll(toMesh_.boundaryMesh(), patchi)
    {
        toMeshPatches_.insert
        (
            toMesh_.boundaryMesh()[patchi].name(),
            patchi
        );
    }

    forAll(cuttingPatchNames, i)
    {
        if (toMeshPatches_.found(cuttingPatchNames[i]))
        {
            cuttingPatches_.insert
            (
                cuttingPatchNames[i],
                toMeshPatches_.find(cuttingPatchNames[i])()
            );
        }
        else
        {
            WarningInFunction
                << "Cannot find cutting-patch " << cuttingPatchNames[i]
                << " in destination mesh" << endl;
        }
    }

    forAll(toMesh_.boundaryMesh(), patchi)
    {
        // Add the processor patches in the toMesh to the cuttingPatches list
        if (isA<processorPolyPatch>(toMesh_.boundaryMesh()[patchi]))
        {
            cuttingPatches_.insert
            (
                toMesh_.boundaryMesh()[patchi].name(),
                patchi
            );
        }
    }

    calcAddressing();
}


Foam::meshToMesh0::meshToMesh0
(
    const fvMesh& meshFrom,
    const fvMesh& meshTo
)
:
    fromMesh_(meshFrom),
    toMesh_(meshTo),
    cellAddressing_(toMesh_.nCells()),
    boundaryAddressing_(toMesh_.boundaryMesh().size()),
    inverseDistanceWeightsPtr_(nullptr),
    inverseVolumeWeightsPtr_(nullptr),
    cellToCellAddressingPtr_(nullptr),
    V_(0)
{
    // check whether both meshes have got the same number
    // of boundary patches
    if (fromMesh_.boundary().size() != toMesh_.boundary().size())
    {
        FatalErrorInFunction
            << "Incompatible meshes: different number of patches, "
            << "fromMesh = " << fromMesh_.boundary().size()
            << ", toMesh = " << toMesh_.boundary().size()
            << exit(FatalError);
    }

    forAll(fromMesh_.boundaryMesh(), patchi)
    {
        if
        (
            fromMesh_.boundaryMesh()[patchi].name()
         != toMesh_.boundaryMesh()[patchi].name()
        )
        {
            FatalErrorInFunction
                << "Incompatible meshes: different patch names for patch "
                << patchi
                << ", fromMesh = " << fromMesh_.boundary()[patchi].name()
                << ", toMesh = " << toMesh_.boundary()[patchi].name()
                << exit(FatalError);
        }

        if
        (
            fromMesh_.boundaryMesh()[patchi].type()
         != toMesh_.boundaryMesh()[patchi].type()
        )
        {
            FatalErrorInFunction
                << "Incompatible meshes: different patch types for patch "
                << patchi
                << ", fromMesh = " << fromMesh_.boundary()[patchi].type()
                << ", toMesh = " << toMesh_.boundary()[patchi].type()
                << exit(FatalError);
        }

        fromMeshPatches_.insert
        (
            fromMesh_.boundaryMesh()[patchi].name(),
            patchi
        );

        toMeshPatches_.insert
        (
            toMesh_.boundaryMesh()[patchi].name(),
            patchi
        );

        patchMap_.insert
        (
            toMesh_.boundaryMesh()[patchi].name(),
            fromMesh_.boundaryMesh()[patchi].name()
        );
    }

    calcAddressing();
}


// ************************************************************************* //
