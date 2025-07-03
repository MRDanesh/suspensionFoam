#include "createRDeltaT.H"

Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",						// Reading properties from transportProperties file
        runTime.constant(),							// It means the transportProperties file is in constant Foulder
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

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

// Temperature Field
volScalarField T
(
   IOobject
  (
       "T",
       runTime.timeName(),
       mesh,
       IOobject::MUST_READ,
       IOobject::AUTO_WRITE
  ),
   mesh
);
//

// Thermal Diffusivity
dimensionedScalar kappa
(
    transportProperties.lookup("kappa")
);
//

// Thermal Expansion 
dimensionedScalar beta
(
    transportProperties.lookup("beta")
);
//

// Reference Temperature 
dimensionedScalar Tref
(
    transportProperties.lookup("Tref")
);
//

// Fluid density 
dimensionedScalar rho_f
(
    transportProperties.lookup("rho_f")
);

// Particle density
dimensionedScalar rho_p
(
    transportProperties.lookup("rho_p")
);

// Particle radius
dimensionedScalar a
(
    transportProperties.lookup("a")
);

// fluid's dynamic viscosity
dimensionedScalar eta_f
(
    transportProperties.lookup("eta_f")
);

volScalarField eta_s
(
  IOobject
  (
  "eta_s",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar("eta_s", dimensionSet(1, -1, -1, 0, 0, 0, 0), (eta_f).value())  // replace with actual value
);

// initial suspension volume fraction
dimensionedScalar gl0
(
    transportProperties.lookup("gl0")
);

#include "readGravitationalAcceleration.H"
#include "createPhi.H"

volScalarField gl
(
  IOobject
  (
  "gl",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar("gl", dimless, (gl0).value())  // replace with actual value
);

// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject("rho", runTime.timeName(), mesh),
    mesh,
    dimensionedScalar("rho", dimDensity, (rho_p).value())  // replace with actual value
);

rho.oldTime();


// Mass flux
// Initialisation does not matter because rhoPhi is reset after the
// gl solution, before it is used in the U equation.
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho_f*phi
);




label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvOptions.H"
