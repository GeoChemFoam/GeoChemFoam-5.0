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
    Foam::diffInterfaceProperties

Description
    Contains the interface properties.

    Properties to aid interFoam:
    -# Correct the alpha boundary condition for dynamic contact angle.
    -# Calculate interface curvature.

SourceFiles
    diffInterfaceProperties.C

\*---------------------------------------------------------------------------*/

#ifndef diffInterfaceProperties_H
#define diffInterfaceProperties_H

#include "IOdictionary.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class diffInterfaceProperties Declaration
\*---------------------------------------------------------------------------*/

class diffInterfaceProperties
{
    // Private data

        //- Keep a reference to the transportProperties dictionary
        const dictionary& transportPropertiesDict_;

        /*----------------------------------------------*/
        //- Surface tension
        dimensionedScalar sigma_;

        //- interface mean-field thickness
        dimensionedScalar eps_;

        //- interface mixing force magnitude
        dimensionedScalar lambda_;

        //- interface free mixture energy density
        dimensionedScalar beta_;

        //-interface mobility
        dimensionedScalar mob_;


        const volScalarField& C_;
        const volVectorField& U_;

        //chemical potential
        volScalarField eta_;

    // Private Member Functions

        //- No copy construct
        diffInterfaceProperties(const diffInterfaceProperties&) = delete;

        //- No copy assignment
        void operator=(const diffInterfaceProperties&) = delete;

        //- Re-calculate the interface curvature
        void calculateEta();

public:

    // Constructors

        //- Construct from volume fraction field gamma and IOdictionary
        diffInterfaceProperties
        (
            const volScalarField& C,
            const volVectorField& U,
            const IOdictionary&
        );

        const volScalarField& eta() const
        {
            return eta_;
        }

        const dimensionedScalar& sigma() const
        {
            return sigma_;
        }

        const dimensionedScalar& eps() const
        {
            return eps_;
        }

        const dimensionedScalar& lambda() const
        {
            return lambda_;
        }

        const dimensionedScalar& beta() const
        {
            return beta_;
        }

        const dimensionedScalar& mob() const
        {
            return mob_;
        }

        tmp<volScalarField> betaEta() const
        {
            return beta_*eta_;
        }


    // Member Functions

        void correct();
        
        //- Read transportProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
