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

#include "globalConcentrationMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::globalConcentrationMixedFvPatchScalarField::
globalConcentrationMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    alphaName_("none"),
	H_(0.0)
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


Foam::globalConcentrationMixedFvPatchScalarField::
globalConcentrationMixedFvPatchScalarField
(
    const globalConcentrationMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    alphaName_(ptf.alphaName_),
	H_(ptf.H_)
{}


Foam::globalConcentrationMixedFvPatchScalarField::
globalConcentrationMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    alphaName_("none"),
	H_(0.0)
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

	IOdictionary transportProperties 
	(
		IOobject
		(
			"transportProperties",
			this->db().time().constant(),
			this->db(),
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	);

    alphaName_ = IOobject::groupName("alpha", transportProperties.get<wordList>("phases")[0]);
	IOdictionary thermoPhysicalProperties 
	(
		IOobject
		(
			"thermoPhysicalProperties",
			this->db().time().constant(),
			this->db(),
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	);

	const dictionary& solutionSpeciesDict = thermoPhysicalProperties.subDict("solutionSpecies");

    word speciesName(iF.name());

	const dictionary& subdict = solutionSpeciesDict.subDict(speciesName);

   	H_ = dimensionedScalar("H",subdict).value();
}


Foam::globalConcentrationMixedFvPatchScalarField::
globalConcentrationMixedFvPatchScalarField
(
    const globalConcentrationMixedFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf),
    alphaName_(tppsf.alphaName_),
    H_(tppsf.H_)
{}


Foam::globalConcentrationMixedFvPatchScalarField::
globalConcentrationMixedFvPatchScalarField
(
    const globalConcentrationMixedFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF),
    alphaName_(tppsf.alphaName_),
    H_(tppsf.H_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::globalConcentrationMixedFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
	const volScalarField& alpha1 = this->db().objectRegistry::lookupObject<volScalarField>(alphaName_);

	const fvPatchField<scalar>& ap = patch().patchField<volScalarField, scalar>(alpha1);
	scalarField api=ap.patchInternalField(); 

    scalarField lambda = -(1-H_)/(ap+H_*(1-ap))*(ap-api);

    valueFraction() = lambda/(lambda + 1);

    mixedFvPatchScalarField::updateCoeffs();
}

void Foam::globalConcentrationMixedFvPatchScalarField::write(Ostream& os) const
{
    mixedFvPatchScalarField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        globalConcentrationMixedFvPatchScalarField
    );
}

// ************************************************************************* //
