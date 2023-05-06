Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

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

volScalarField eps
(
    IOobject
    (
        "eps",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

surfaceScalarField nEpsf
(
   IOobject
   (
      "nEpsf",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
   ),
   mesh,
   dimensionedScalar("nEpsf",dimArea, 0.0)
);

Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

#include "createPhi.H"

Info << "Create species mixture\n" << endl;
multiComponentTransportMixture<reactiveMixture> speciesMixture(mesh);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

singlePhaseTransportModel laminarTransport(U, phi);


dimensionedScalar nu("nu",dimViscosity,transportProperties);


autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

dimensionedScalar kf("kf",transportProperties);

dimensionedScalar Mw("Mw",transportProperties);

dimensionedScalar rhos("rhos", transportProperties);

//permeability
volScalarField Kinv = kf*pow(1-eps,2)/pow(eps,3);

Info << "Reading reaction rate if present\n" << endl;
volScalarField R
(
    IOobject
    (
        "R",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("R",dimless/dimTime, 0.0)
);

dimensionedScalar deltaN = 1e-13/pow(average(mesh.V()), 1.0/3.0);

#include "createMRF.H"
#include "createFvOptions.H"
