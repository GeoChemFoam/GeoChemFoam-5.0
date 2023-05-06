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
    Foam::collapsingFace

Description
    A \c topoSetFaceSource to select faces that collapse into internal faces.

    Operands:
    \table
      Operand   | Type    | Location
      output    | faceSet | $FOAM_CASE/constant/polyMesh/sets/\<set\>
    \endtable

Usage
    Minimal example by using \c system/topoSetDict.actions:
    \verbatim
    {
        // Mandatory (inherited) entries
        name        <name>;
        type        faceSet;
        action      <action>;

        // Mandatory entries
        source      collapsingFace;
        tol      1e-5;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                         | Type | Req'd | Dflt
      name       | Name of faceSet                     | word |  yes  | -
      type       | Type name: faceSet                  | word |  yes  | -
      action     | Action applied on faces - see below | word |  yes  | -
      source     | Source name: collapsingFace           | word |  yes  | -
      tol        | acdeptd distance between two internal faces | scalar | yes | -


    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new faceSet from selected faces
      add      | Add selected faces into this faceSet
      subtract | Remove selected faces from this faceSet
    \endverbatim

See also
    - Foam::topoSetSource
    - Foam::topoSetFaceSource

SourceFiles
    collapsingFace.C

\*---------------------------------------------------------------------------*/

#ifndef collapsingFace_H
#define collapsingFace_H

#include "topoSetFaceSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class collapsingFace Declaration
\*---------------------------------------------------------------------------*/

class collapsingFace
:
    public topoSetFaceSource
{
    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- Tolerance (i.e. cos of angle between normal_ and faceNormal)
        const scalar tol_;


    // Private Member Functions

		bool checkHit(label faceI, label faceJ) const;

		bool checkIntersect(label faceI,label faceJ) const;


public:

    //- Runtime type information
    TypeName("collapsingFace");


    // Constructors

        //- Construct from components
        collapsingFace
        (
            const polyMesh& mesh,
            const scalar tol
        );

        //- Construct from dictionary
        collapsingFace(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        collapsingFace(const polyMesh& mesh, Istream&);


    //- Destructor
    virtual ~collapsingFace() = default;


    // Member Functions

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet& set
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
