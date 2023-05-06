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

#ifndef MapConsistentVolFields_H
#define MapConsistentVolFields_H

#include "GeometricField.H"
#include "meshToMesh0.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type, class CombineOp>
void MapConsistentVolFields
(
    const IOobjectList& objects,
    const meshToMesh0& meshToMesh0Interp,
    const meshToMesh0::order& mapOrder,
    const CombineOp& cop
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    const fvMesh& meshSource = meshToMesh0Interp.fromMesh();
    const fvMesh& meshTarget = meshToMesh0Interp.toMesh();

    IOobjectList fields = objects.lookupClass(fieldType::typeName);

    forAllConstIters(fields, fieldIter)
    {
        Info<< "    interpolating " << (*fieldIter)->name()
            << endl;

        // Read field. Do not auto-load old-time field
        fieldType fieldSource(*fieldIter(), meshSource, false);

        IOobject fieldTargetIOobject
        (
            (*fieldIter)->name(),
            meshTarget.time().timeName(),
            meshTarget,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        if (fieldTargetIOobject.typeHeaderOk<fieldType>(true))
        {
            // Read fieldTarget
            fieldType fieldTarget(fieldTargetIOobject, meshTarget, false);

            // Interpolate field
            meshToMesh0Interp.interpolate
            (
                fieldTarget,
                fieldSource,
                mapOrder,
                cop
            );

            // Write field
            fieldTarget.write();
        }
        else
        {
            fieldTargetIOobject.readOpt(IOobject::NO_READ);

            // Interpolate field
            fieldType fieldTarget
            (
                fieldTargetIOobject,
                meshToMesh0Interp.interpolate
                (
                    fieldSource,
                    mapOrder,
                    cop
                )
            );

            // Write field
            fieldTarget.write();
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
