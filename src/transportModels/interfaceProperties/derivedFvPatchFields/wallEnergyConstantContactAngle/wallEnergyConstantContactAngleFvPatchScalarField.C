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
    the Carbonate Reservoir Group at Heriot-Watt University. Please visit our
    website for more information <https://carbonates.hw.ac.uk>.

\*---------------------------------------------------------------------------*/

#include "wallEnergyConstantContactAngleFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallEnergyConstantContactAngleFvPatchScalarField::
wallEnergyConstantContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    theta0_(0.0),
    eps_(0.0)
{}


Foam::wallEnergyConstantContactAngleFvPatchScalarField::
wallEnergyConstantContactAngleFvPatchScalarField
(
    const wallEnergyConstantContactAngleFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    theta0_(ptf.theta0_),
    eps_(ptf.eps_)
{}


Foam::wallEnergyConstantContactAngleFvPatchScalarField::
wallEnergyConstantContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    theta0_(readScalar(dict.lookup("theta0"))),
    eps_(0.0)
{
    fvPatchField<scalar>::operator=(patchInternalField());

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

    eps_=dimensionedScalar("eps",transportProperties).value();

    updateCoeffs();
}



Foam::wallEnergyConstantContactAngleFvPatchScalarField::
wallEnergyConstantContactAngleFvPatchScalarField
(
    const wallEnergyConstantContactAngleFvPatchScalarField& tppsf
)
:
    fixedGradientFvPatchScalarField(tppsf),
    theta0_(tppsf.theta0_),
    eps_(tppsf.eps_)
{}


Foam::wallEnergyConstantContactAngleFvPatchScalarField::
wallEnergyConstantContactAngleFvPatchScalarField
(
    const wallEnergyConstantContactAngleFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(tppsf, iF),
    theta0_(tppsf.theta0_),
    eps_(tppsf.eps_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wallEnergyConstantContactAngleFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalar theta = degToRad()*theta0_;

    scalar b1 = cos(theta);

    scalarField Cp = *this;

    Info << eps_ << endl;

    gradient() = 0.5*sqrt(2.0)*b1*(1-Cp*Cp)/eps_;

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::wallEnergyConstantContactAngleFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("theta0",theta0_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        wallEnergyConstantContactAngleFvPatchScalarField
    );
}

// ************************************************************************* //
