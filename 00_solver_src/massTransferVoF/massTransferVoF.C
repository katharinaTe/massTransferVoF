/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "massTransferVoF.H"
#include "localEulerDdtScheme.H"
#include "fvCorrectPhi.H"
#include "geometricZeroField.H"
#include "addToRunTimeSelectionTable.H"

// new code
#include "fvmLaplacian.H"
#include "fvcSnGrad.H"
#include "fvmDiv.H"


// end new code

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(massTransferVoF, 0);
    addToRunTimeSelectionTable(solver, massTransferVoF, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::massTransferVoF::massTransferVoF(fvMesh& mesh)
:
    twoPhaseVoFSolver
    (
        mesh,
        autoPtr<twoPhaseVoFMixture>(new incompressibleTwoPhaseVoFMixture(mesh))
    ),

    mixture
    (
        refCast<incompressibleTwoPhaseVoFMixture>(twoPhaseVoFSolver::mixture)
    ),

    p
    (
        IOobject
        (
            "p",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rho*buoyancy.gh
    ),
    
    // new code
        C
    (
        IOobject
        (
            "C",
            runTime.name(),
            mesh,
            //IOobject::MUST_READ_IF_MODIFIED,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ), 
    
    
    transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    Dwater
    (
        transportProperties.lookup("Dwater")
    ),

    Dair
    (
        transportProperties.lookup("Dair")
    ),

    Schmidtnumber
    (
        transportProperties.lookup("Schmidtnumber")
    ),
    
    Temp
    (
        transportProperties.lookup("Temp")
    ),
   
    DT
    (
        IOobject
        (
            "DT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Dwater * Dair / (alpha1 * Dair + (1 - alpha1) * Dwater)
    ),
    

    He("He",dimless, 1),
    

    // end new code
    pressureReference_
    (
        p,
        p_rgh,
        pimple.dict()
    ),

    momentumTransport
    (
        U,
        phi,
        alphaPhi1,
        mixture
    )
{
    // Read the controls
    readControls();
    
    // start new code
    
    if (transportProperties.found("H2S"))
    {
    Info<< "Using H2S Henry-law correlation" << nl;
    scalar Hcp = transportProperties.subDict("H2S").lookupOrDefault<scalar>("Hcp", 1e-3);
    scalar Ctemp = transportProperties.subDict("H2S").lookupOrDefault<scalar>("Ctemp", 2200);

    scalar Hcp_Temp = 1.0/(Hcp*Foam::exp(Ctemp*(1.0/Temp.value()-1.0/298.15))*8.314* Temp.value());
    He = dimensionedScalar("He", dimless, Hcp_Temp);
    }
    else if (transportProperties.found("O2"))
    {
    Info<< "Using O2 Henry-law correlation" << nl;
    scalar Hcp = transportProperties.subDict("O2").lookupOrDefault<scalar>("Hcp", 1.3e-5);
    scalar Ctemp = transportProperties.subDict("O2").lookupOrDefault<scalar>("Ctemp", 1700);

    scalar Hcp_Temp = 1.0/(Hcp*Foam::exp(Ctemp*(1.0/Temp.value()-1.0/298.15))*8.314* Temp.value());
    He = dimensionedScalar("He", dimless, Hcp_Temp);
    }
    else        // user-specified constant Hcpgen - general value for Hcp (not temperature-dependent)
    {
    Info<< "Using custom constant Henry coefficient" << nl;
    scalar Hcpgen = transportProperties.lookup<scalar>("HenryConstant");
    He = dimensionedScalar("He", dimless, Hcpgen);
    }
    
    // end new code

    if (correctPhi || mesh.topoChanging())
    {
        rAU = new volScalarField
        (
            IOobject
            (
                "rAU",
                runTime.name(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimTime/dimDensity, 1)
        );
    }

    if (!runTime.restart() || !divergent())
    {
        correctUphiBCs(U_, phi_, true);

        fv::correctPhi
        (
            phi_,
            U,
            p_rgh,
            rAU,
            autoPtr<volScalarField>(),
            pressureReference(),
            pimple
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::massTransferVoF::~massTransferVoF()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::massTransferVoF::prePredictor()
{
    twoPhaseVoFSolver::prePredictor();

    const dimensionedScalar& rho1 = mixture.rho1();
    const dimensionedScalar& rho2 = mixture.rho2();

    // Calculate the mass-flux from the accumulated alphaPhi1
    rhoPhi = (alphaPhi1*(rho1 - rho2) + phi*rho2);

    if (pimple.predictTransport())
    {
        momentumTransport.predict();
    }
}


void Foam::solvers::massTransferVoF::pressureCorrector()
{
    incompressiblePressureCorrector(p);
}


void Foam::solvers::massTransferVoF::thermophysicalPredictor()
{}


void Foam::solvers::massTransferVoF::postCorrector()
{
    if (pimple.correctTransport())
    {
        momentumTransport.correct();
    }
}

void Foam::solvers::massTransferVoF::postSolve()
{    
    // new code from here
        
           
 
fvScalarMatrix CEqn =
        fvm::ddt(C)
      + fvm::div(phi, C, "div(phi,C)")
      - fvm::laplacian(fvc::interpolate(DT), C, "laplacian(C)");    

if (mesh_.foundObject<volScalarField>("nut"))
{
        
    const volScalarField& nutField =
        mesh_.lookupObject<volScalarField>("nut");   // <- the key line
        
    // effective diffusivity for the species
    surfaceScalarField gammaEff =
        fvc::interpolate(nutField)/Schmidtnumber;

    surfaceScalarField phiCiTurb =
        (
            (
                (fvc::interpolate(DT)+gammaEff) * (1-He) / (fvc::interpolate(alpha1)+(1-fvc::interpolate(alpha1))*He)
            )
            * fvc::snGrad(alpha1)
        ) 
        * mesh.magSf();
        
    CEqn -= fvm::laplacian(gammaEff, C);   
    CEqn += fvm::div(phiCiTurb, C, "div(phi,C)");     

}
else
{
    surfaceScalarField phiCiLam =
        (
            (
                fvc::interpolate(DT) * (1-He) / (fvc::interpolate(alpha1)+(1-fvc::interpolate(alpha1))*He)
            )
            * fvc::snGrad(alpha1)
        ) 
        * mesh.magSf();
    CEqn += fvm::div(phiCiLam, C, "div(phi,C)");
}


    


    CEqn.solve();
        
    // end new code*/
}

// ************************************************************************* //
