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

Class
    Foam::wallEnergyConstantContactAngleFvPatchScalarField

Description
    Impose contact angle using
    grad(C).n = sqrt(2)/2 *(1-C^2)/eps 

SourceFiles
    wallEnergyConstantContactAngleFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef wallEnergyConstantContactAngleFvPatchScalarField_H
#define wallEnergyConstantContactAngleFvPatchScalarField_H

#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class wallEnergyConstantContactAngleFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class wallEnergyConstantContactAngleFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{

    // Private data
    scalar theta0_;

    scalar eps_;

public:

    //- Runtime type information
    TypeName("wallEnergyConstantContactAngle");


    // Constructors

        //- Construct from patch and internal field
        wallEnergyConstantContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        wallEnergyConstantContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given wallEnergyConstantContactAngleFvPatchScalarField
        // onto a new patch
        wallEnergyConstantContactAngleFvPatchScalarField
        (
            const wallEnergyConstantContactAngleFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        wallEnergyConstantContactAngleFvPatchScalarField
        (
            const wallEnergyConstantContactAngleFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new wallEnergyConstantContactAngleFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        wallEnergyConstantContactAngleFvPatchScalarField
        (
            const wallEnergyConstantContactAngleFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new wallEnergyConstantContactAngleFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
