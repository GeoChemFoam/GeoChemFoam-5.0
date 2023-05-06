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
    Foam::steadyStateControl

Description
    Steady state control class to supply convergence information/checks for
    the steady state loop.

\*---------------------------------------------------------------------------*/

#ifndef steadyStateControl_H
#define steadyStateControl_H

#include "solutionControl.H"

//- Declare that steady state Control will be used
#define STEADYSTATE_CONTROL

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class steadyStateControl Declaration
\*---------------------------------------------------------------------------*/

class steadyStateControl
:
    public solutionControl
{

protected:

    // Protected Data

        //- Initialised flag
        bool initialised_;

		label iter_counter;

       label nCorr_;


    // Protected Member Functions

        //- Read controls from fvSolution dictionary
        bool read();

private:

    // Private member functions

        //- No copy construct
        steadyStateControl(const steadyStateControl&) = delete;

        //- No copy assignment
        void operator=(const steadyStateControl&) = delete;


public:

    // Static Data Members

        //- Run-time type information
        TypeName("steadyStateControl");


    // Constructors

        //- Construct from mesh and the name of control sub-dictionary
        steadyStateControl
        (
            fvMesh& mesh,
            const word& dictName = "STEADYSTATE",
            const bool verbose = true
        );


    //- Destructor
    virtual ~steadyStateControl() = default;

    //- Return true if all convergence checks are satisfied
    bool criteriaSatisfied();


    // Member Functions

        // Solution control

            //- SIMPLE loop
            virtual bool loop();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
