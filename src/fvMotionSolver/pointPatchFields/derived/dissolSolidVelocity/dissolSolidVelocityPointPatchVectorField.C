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

#include "dissolSolidVelocityPointPatchVectorField.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "polyMesh.H"
#include "primitivePatchInterpolation.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dissolSolidVelocityPointPatchVectorField::
dissolSolidVelocityPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    CName_("C"),
    k_(0.0),
    rhos_(0.0),
    Mw_(0.0)
{}


dissolSolidVelocityPointPatchVectorField::
dissolSolidVelocityPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    CName_(dict.lookup("CName")),
    k_(dict.get<scalar>("k")),
    rhos_(dict.get<scalar>("rhos")),
    Mw_(dict.get<scalar>("Mw"))
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }
}


dissolSolidVelocityPointPatchVectorField::
dissolSolidVelocityPointPatchVectorField
(
    const dissolSolidVelocityPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    CName_(ptf.CName_),
    k_(ptf.k_),
    rhos_(ptf.rhos_),
    Mw_(ptf.Mw_)
{}


dissolSolidVelocityPointPatchVectorField::
dissolSolidVelocityPointPatchVectorField
(
    const dissolSolidVelocityPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    CName_(ptf.CName_),
    k_(ptf.k_),
    rhos_(ptf.rhos_),
    Mw_(ptf.Mw_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dissolSolidVelocityPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const polyMesh& mesh = this->internalField().mesh()();

    vectorField nd = patch().pointNormals();

    //- set-up interpolator
    primitivePatchInterpolation patchInterpolator
    (
        mesh.boundaryMesh()[patch().index()]
    );


    const volScalarField& C = this->db().objectRegistry::lookupObject<volScalarField>(CName_);
    const fvPatchField<scalar>& Cp = C.boundaryField()[patch().index()];
    scalarField pointValues =
    patchInterpolator.faceToPointInterpolate<scalar>
    (
        Cp
    );


    Field<vector>::operator=
    (
        nd*k_*max(pointValues,0.0)*Mw_/rhos_
    );

    fixedValuePointPatchField<vector>::updateCoeffs();
}


void dissolSolidVelocityPointPatchVectorField::write(Ostream& os) const
{
    pointPatchField<vector>::write(os);
    os.writeEntry("CName", CName_);
    os.writeEntry("k",k_);
    os.writeEntry("rhos",rhos_);
    os.writeEntry("Mw",Mw_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    dissolSolidVelocityPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
