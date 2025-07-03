// Adding Temperature Equation

volScalarField kappaRho = kappa * rho;

fvScalarMatrix TEqn
(
    fvm::ddt(rho, T)+fvm::div(rhoPhi,T)-fvm::laplacian(kappaRho,T) 
);