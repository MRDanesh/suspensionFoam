// Conservation equation of non-colloidal particles

// Parameters
scalar Kc = 0.41;
scalar Keta = 0.62;

// Compute strain rate magnitude gammaDot = sqrt(2*S:S)
volSymmTensorField S(symm(fvc::grad(U)));
volScalarField gammaDot("gammaDot", sqrt(2.0)*mag(S));

// Compute log(eta_s), where eta_s is suspension viscosity
volScalarField logEta("logEta", log(eta_s)); // eta must be declared somewhere


fvScalarMatrix glEqn
(
    fvm::ddt(gl)+fvm::div(phi,gl)+fvm::div(kappa,gl)
);


rho

rhoPhi