
rho = gl * rho_p + (1.0-gl) * rho_f;
surfaceScalarField glf = fvc::interpolate(gl);

rhoPhi = glf * phi * rho_p + (scalar(1.0) - glf) * phi * rho_f;
eta_s = eta_f;

// Solve the Momentum equation

// Ensure boundary values match internal logic
//eta_s.correctBoundaryConditions();

// Optional: Smooth the field numerically for stability
//eta_s.relax();

MRF.correctBoundaryVelocity(U);

tmp<volTensorField> tgradU = fvc::grad(U);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
  - fvc::div(eta_s * dev2(2*symm(fvc::grad(U))-fvc::grad(U))) - fvm::laplacian(eta_s, U)
  + rho*g*beta*(T-Tref)
);
// Note: if there is confusion, check https://www.cfd-online.com/Forums/openfoam-solving/58214-calculating-divdevreff.html
//+ turbulence->divDevRhoReff(rho, U)
fvVectorMatrix& UEqn = tUEqn.ref();
UEqn.relax();
fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == fvc::grad(p));

    fvOptions.correct(U);
}

