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

#include "reactiveSurfaceConcentrationMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reactiveSurfaceConcentrationMixedFvPatchScalarField::
reactiveSurfaceConcentrationMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    k_(0.0),
    scoeff_(0.0)
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


Foam::reactiveSurfaceConcentrationMixedFvPatchScalarField::
reactiveSurfaceConcentrationMixedFvPatchScalarField
(
    const reactiveSurfaceConcentrationMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    k_(ptf.k_),
    scoeff_(ptf.scoeff_)
{}


Foam::reactiveSurfaceConcentrationMixedFvPatchScalarField::
reactiveSurfaceConcentrationMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    k_(readScalar(dict.lookup("k"))),
    scoeff_(readScalar(dict.lookup("scoeff")))
{
    this->refValue() = Zero;

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(refValue());
    }

    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


Foam::reactiveSurfaceConcentrationMixedFvPatchScalarField::
reactiveSurfaceConcentrationMixedFvPatchScalarField
(
    const reactiveSurfaceConcentrationMixedFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf),
    k_(tppsf.k_),
    scoeff_(tppsf.scoeff_)
{}


Foam::reactiveSurfaceConcentrationMixedFvPatchScalarField::
reactiveSurfaceConcentrationMixedFvPatchScalarField
(
    const reactiveSurfaceConcentrationMixedFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF),
    k_(tppsf.k_),
    scoeff_(tppsf.scoeff_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::reactiveSurfaceConcentrationMixedFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    const dictionary& speciesTransportProperties = db().lookupObject<IOdictionary> ("thermoPhysicalProperties");

    const dictionary& solutionSpeciesDict = speciesTransportProperties.subDict("solutionSpecies");

    const dictionary& subdict = solutionSpeciesDict.subDict(this->internalField().name());

    dimensionedScalar D("D",subdict);


    scalarField lambda = scoeff_*k_/D.value()/(this->patch().deltaCoeffs());

    valueFraction() = lambda/(lambda + 1);

    mixedFvPatchScalarField::updateCoeffs();
}

void Foam::reactiveSurfaceConcentrationMixedFvPatchScalarField::write(Ostream& os) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("k") << k_ << token::END_STATEMENT << nl;
    os.writeKeyword("scoeff") << scoeff_ << token::END_STATEMENT << nl;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        reactiveSurfaceConcentrationMixedFvPatchScalarField
    );
}

// ************************************************************************* //
