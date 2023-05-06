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

#include "diffInterfaceProperties.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcAverage.H"
#include "fvCFD.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::diffInterfaceProperties::calculateEta()
{
    //free mixture energy demsity
    volScalarField psi =C_*C_*C_-C_;

    //chemical potential
    eta_ = psi-sqr(eps_)*fvc::laplacian(C_);
    eta_.correctBoundaryConditions();
}

Foam::diffInterfaceProperties::diffInterfaceProperties
(
    const volScalarField& C,
    const volVectorField& U,
    const IOdictionary& dict
)
:
    transportPropertiesDict_(dict),
    sigma_(dimensionedScalar("sigma",dict)),
    eps_(dimensionedScalar("eps",dict)),
    lambda_(3*sigma_*eps_/std::sqrt(8.0)),
    beta_(3*sigma_/eps_/std::sqrt(8.0)),
    mob_(dimensionedScalar("mob",dict)),
    C_(C),
    U_(U),
    eta_
    (
        IOobject
        (
            "eta",
            C_.time().timeName(),
            C_.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        C_.mesh()
    )
{
    calculateEta();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::diffInterfaceProperties::correct()
{
    calculateEta();
}


bool Foam::diffInterfaceProperties::read()
{
    return true;
}


// ************************************************************************* //
