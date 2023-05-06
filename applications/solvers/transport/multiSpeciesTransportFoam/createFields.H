Info<< "Reading field U\n" << endl;

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

#include "createFvOptions.H"

Info << "Create species mixture\n" << endl;
multiComponentTransportMixture<basicMultiComponentMixture> speciesMixture(mesh);

