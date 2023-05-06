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
    Foam::twoPhasePhaseFieldMixture

Description
    A two-phase iphase-field mixture model

SourceFiles
    twoPhasePhaseFieldMixture.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhasePhaseFieldMixture_H
#define twoPhasePhaseFieldMixture_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class twoPhasePhaseFieldMixture Declaration
\*---------------------------------------------------------------------------*/

class twoPhasePhaseFieldMixture
{
protected:

    // Protected data

        word phase1Name_;
        word phase2Name_;

        volScalarField C_;
        volScalarField alpha1_;
	volScalarField alpha2_;


public:

    // Constructors

        //- Construct from components
        twoPhasePhaseFieldMixture
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    ~twoPhasePhaseFieldMixture()
    {}


    // Member Functions

        const word& phase1Name() const
        {
            return phase1Name_;
        }

        const word& phase2Name() const
        {
            return phase2Name_;
        }

        //- Return the phase-field fraction
        const volScalarField& C() const
        {
            return C_;
        }

        //- Return the phase-field fraction 
        volScalarField& C()
        {
            return C_;
        }



        //- Return the phase-fraction of phase 1
        const volScalarField& alpha1() const
        {
            return alpha1_;
        }

        //- Return the phase-fraction of phase 1
        volScalarField& alpha1()
        {
            return alpha1_;
        }

        //- Return the phase-fraction of phase 2
        const volScalarField& alpha2() const
        {
            return alpha2_;
        }

        //- Return the phase-fraction of phase 2
        volScalarField& alpha2()
        {
            return alpha2_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
