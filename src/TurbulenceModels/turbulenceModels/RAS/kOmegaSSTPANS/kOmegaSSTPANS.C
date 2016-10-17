/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

#include "kOmegaSSTPANS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTPANS<BasicTurbulenceModel>::kOmegaSSTPANS::F1
(
    const volScalarField& CDkOmega
) const
{
    tmp<volScalarField> CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar("1.0e-10", dimless/sqr(dimTime), 1.0e-10)
    );

    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/this->betaStar_)*sqrt(kU_)/(omegaU_*this->y_),
                scalar(500)*(this->mu()/this->rho_)/(sqr(this->y_)*omegaU_)
            ),
            (4*this->alphaOmega2_*(fK_/fOmega_))*kU_
            /(CDkOmegaPlus*sqr(this->y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
kOmegaSSTPANS<BasicTurbulenceModel>::kOmegaSSTPANS::F2() const
{
    tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/this->betaStar_)*sqrt(kU_)/(omegaU_*this->y_),
            scalar(500)*(this->mu()/this->rho_)/(sqr(this->y_)*omegaU_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
kOmegaSSTPANS<BasicTurbulenceModel>::kOmegaSSTPANS::F3() const
{
    tmp<volScalarField> arg3 = min
    (
        150*(this->mu()/this->rho_)/(omegaU_*sqr(this->y_)),
        scalar(10)
    );

    return 1 - tanh(pow4(arg3));
}


template<class BasicTurbulenceModel>
void kOmegaSSTPANS<BasicTurbulenceModel>::correctNut
(
    const volScalarField& S2,
    const volScalarField& F2
)
{
    this->nut_ = this->a1_*kU_/max(this->a1_*omegaU_, this->b1_*F2*sqrt(S2));
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kOmegaSSTPANS<BasicTurbulenceModel>::correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))), this->F23());
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> kOmegaSSTPANS<BasicTurbulenceModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            kU_,
            dimVolume*this->rho_.dimensions()*kU_.dimensions()/dimTime
        )
    );
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> kOmegaSSTPANS<BasicTurbulenceModel>::omegaSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omegaU_,
            dimVolume*this->rho_.dimensions()*omegaU_.dimensions()/dimTime
        )
    );
}

template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> kOmegaSSTPANS<BasicTurbulenceModel>::Qsas
(
    const volScalarField& S2,
    const volScalarField& gamma,
    const volScalarField& beta
) const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omegaU_,
            dimVolume*this->rho_.dimensions()*omegaU_.dimensions()/dimTime
        )
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegaSSTPANS<BasicTurbulenceModel>::kOmegaSSTPANS
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    kOmegaSST<BasicTurbulenceModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    fEpsilon_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "fEpsilon",
            this->coeffDict_,
            1.0
        )
    ),
    uLim_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "fKupperLimit",
            this->coeffDict_,
            1.0
        )
    ),
    loLim_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "fKlowerLimit",
            this->coeffDict_,
            0.1
        )
    ),

    fK_
    (
        IOobject
        (
            IOobject::groupName("fK", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zero", loLim_)
    ),

    fOmega_
    (
        IOobject
        (
            "fOmega",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fEpsilon_/fK_
    ),

    delta_
    (
        LESdelta::New
        (
            IOobject::groupName("delta", U.group()),
            *this,
            this->coeffDict_
        )
    ),

    kU_
    (
        IOobject
        (
            IOobject::groupName("kU", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->k_*fK_,
        this->k_.boundaryField().types()
    ),

    omegaU_
    (
        IOobject
        (
            IOobject::groupName("omegaU", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->omega_*fOmega_,
        this->omega_.boundaryField().types()
    )
{
    bound(kU_, min(fK_)*this->kMin_);
    bound(omegaU_, min(fOmega_)*this->omegaMin_);

    if (type == typeName)
    {
        correctNut();
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kOmegaSSTPANS<BasicTurbulenceModel>::read()
{
    if (kOmegaSST<BasicTurbulenceModel>::read())
    {
        fEpsilon_.readIfPresent(this->coeffDict());
        uLim_.readIfPresent(this->coeffDict());
        loLim_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}

template<class BasicTurbulenceModel>
void kOmegaSSTPANS<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    //eddyViscosity<RASModel<BasicTurbulenceModel> >::correct();
    BasicTurbulenceModel::correct();
    
    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField S2(2*magSqr(symm(tgradU())));
    volScalarField GbyNu((tgradU() && dev(twoSymm(tgradU()))));
    volScalarField G(this->GName(), nut*GbyNu);
    tgradU.clear();

    // Update omegaU and G at the wall
    omegaU_.boundaryFieldRef().updateCoeffs();

    volScalarField CDkOmega
    (
        (2*this->alphaOmega2_*(fK_/fOmega_))*
        (fvc::grad(kU_) & fvc::grad(omegaU_))/omegaU_
    );

    volScalarField F1(this->F1(CDkOmega));
    volScalarField F23(this->F23());

    {
        volScalarField gamma(this->gamma(F1));
        volScalarField beta(this->beta(F1));
        volScalarField betaL
        (
            gamma*this->betaStar_ - (gamma *this->betaStar_/fOmega_)
            + (beta/fOmega_)
        );


        // Unresolved Turbulent frequency equation
        tmp<fvScalarMatrix> omegaUEqn
        (
            fvm::ddt(alpha, rho, omegaU_)
          + fvm::div(alphaRhoPhi, omegaU_)
          - fvm::laplacian(alpha*rho*DomegaUEff(F1), omegaU_)
         ==
            alpha*rho*gamma
           *min
            (
                GbyNu,
                (this->c1_/this->a1_)*this->betaStar_*omegaU_
                *max(this->a1_*omegaU_, this->b1_*F23*sqrt(S2))
            )
          - fvm::SuSp((2.0/3.0)*alpha*rho*gamma*divU, omegaU_)
          - fvm::Sp(alpha*rho*betaL*omegaU_, omegaU_)
          - fvm::SuSp
            (
                alpha*rho*(F1 - scalar(1))*CDkOmega/omegaU_,
                omegaU_
            )
          + Qsas(S2, gamma, beta)
          + omegaSource()
          + fvOptions(alpha, rho, omegaU_)
        );

        omegaUEqn.ref().relax();
        fvOptions.constrain(omegaUEqn.ref());
        omegaUEqn.ref().boundaryManipulate(omegaU_.boundaryFieldRef());
        solve(omegaUEqn);
        fvOptions.correct(omegaU_);
        bound(omegaU_, min(fOmega_)*this->omegaMin_);
    }

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kUEqn
    (
        fvm::ddt(alpha, rho, kU_)
      + fvm::div(alphaRhoPhi, kU_)
      - fvm::laplacian(alpha*rho*DkUEff(F1), kU_)
     ==
        min(alpha*rho*G, (this->c1_*this->betaStar_)*alpha*rho*kU_*omegaU_)
      - fvm::SuSp((2.0/3.0)*alpha*rho*divU, kU_)
      - fvm::Sp(alpha*rho*this->betaStar_*omegaU_, kU_)
      + kSource()
      + fvOptions(alpha, rho, kU_)
    );

    kUEqn.ref().relax();
    fvOptions.constrain(kUEqn.ref());
    solve(kUEqn);
    fvOptions.correct(kU_);
    bound(kU_, min(fK_)*this->kMin_);


    // Calculation of Turbulent kinetic energy and Frequency
    this->k_ = kU_/fK_;
    this->k_.correctBoundaryConditions();

    this->omega_ = omegaU_/fOmega_;
    this->omega_.correctBoundaryConditions();

    bound(this->k_, this->kMin_);
    bound(this->omega_, this->omegaMin_);

    correctNut(S2, F23);

    // Recalculate fK with new kU and epsilonU
    
    // Calculate the unresolved turbulence length scale
    volScalarField lu_
    (
        sqrt(kU_)/(this->betaStar_*omegaU_)
    );

    // update fK
    fK_.primitiveFieldRef() = min
    (
        max(sqrt(1.0/this->betaStar_)*pow(delta()/lu_,2.0/3.0),loLim_),
        uLim_
    );

    // update fOmega
    fOmega_ = fEpsilon_/fK_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
