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

InNamespace
    Foam

Description
    Gets the indices of (source)particles that have been appended to the
    target cloud and maps the lagrangian fields accordingly.

\*---------------------------------------------------------------------------*/

#ifndef MapLagrangianFields_H
#define MapLagrangianFields_H

#include "cloud.H"
#include "GeometricField.H"
#include "meshToMesh0.H"
#include "IOobjectList.H"
#include "CompactIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Gets the indices of (source)particles that have been appended to the
//  target cloud and maps the lagrangian fields accordingly.
template<class Type>
void MapLagrangianFields
(
    const string& cloudName,
    const IOobjectList& objects,
    const meshToMesh0& meshToMesh0Interp,
    const labelList& addParticles
)
{
    const fvMesh& meshTarget = meshToMesh0Interp.toMesh();

    {
        IOobjectList fields = objects.lookupClass(IOField<Type>::typeName);

        forAllConstIters(fields, fieldIter)
        {
            Info<< "    mapping lagrangian field "
                << (*fieldIter)->name() << endl;

            // Read field (does not need mesh)
            IOField<Type> fieldSource(*fieldIter());

            // Map
            IOField<Type> fieldTarget
            (
                IOobject
                (
                    (*fieldIter)->name(),
                    meshTarget.time().timeName(),
                    cloud::prefix/cloudName,
                    meshTarget,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                addParticles.size()
            );

            forAll(addParticles, i)
            {
                fieldTarget[i] = fieldSource[addParticles[i]];
            }

            // Write field
            fieldTarget.write();
        }
    }

    {
        IOobjectList fieldFields =
            objects.lookupClass(IOField<Field<Type>>::typeName);

        forAllConstIters(fieldFields, fieldIter)
        {
            Info<< "    mapping lagrangian fieldField "
                << (*fieldIter)->name() << endl;

            // Read field (does not need mesh)
            IOField<Field<Type>> fieldSource(*fieldIter());

            // Map - use CompactIOField to automatically write in
            // compact form for binary format.
            CompactIOField<Field<Type>, Type> fieldTarget
            (
                IOobject
                (
                    (*fieldIter)->name(),
                    meshTarget.time().timeName(),
                    cloud::prefix/cloudName,
                    meshTarget,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                addParticles.size()
            );

            forAll(addParticles, i)
            {
                fieldTarget[i] = fieldSource[addParticles[i]];
            }

            // Write field
            fieldTarget.write();
        }
    }

    {
        IOobjectList fieldFields =
            objects.lookupClass(CompactIOField<Field<Type>, Type>::typeName);

        forAllConstIters(fieldFields, fieldIter)
        {
            Info<< "    mapping lagrangian fieldField "
                << (*fieldIter)->name() << endl;

            // Read field (does not need mesh)
            CompactIOField<Field<Type>, Type> fieldSource(*fieldIter());

            // Map
            CompactIOField<Field<Type>, Type> fieldTarget
            (
                IOobject
                (
                    (*fieldIter)->name(),
                    meshTarget.time().timeName(),
                    cloud::prefix/cloudName,
                    meshTarget,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                addParticles.size()
            );

            forAll(addParticles, i)
            {
                fieldTarget[i] = fieldSource[addParticles[i]];
            }

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
