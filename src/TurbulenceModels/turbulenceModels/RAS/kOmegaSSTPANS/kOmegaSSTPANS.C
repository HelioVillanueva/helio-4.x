/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kOmegaSSTPANS<BasicTurbulenceModel>::correctPANSCoeffs()
{
    fK_ = min
    (
        max
        (
            sqrt(betaStar_.value())*(pow(pow(cellVolume,1.0/3.0)*
                (betaStar_.value()*omegaU_)/sqrt(kU_),2.0/3.0)),
            loLimVec
        ),
        uLimVec
    );

    fOmega_ = fEpsilon_/fK_;

}

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
                (scalar(1)/betaStar_)*sqrt(kU_)/(omegaU_*y_),
                scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omegaU_)
            ),
            (4*alphaOmega2_*(fK_/fOmega_))*kU_/(CDkOmegaPlus*sqr(y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTPANS<BasicTurbulenceModel>::kOmegaSSTPANS::F2() const
{
    tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(kU_)/(omegaU_*y_),
            scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omegaU_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTPANS<BasicTurbulenceModel>::kOmegaSSTPANS::F3() const
{
    tmp<volScalarField> arg3 = min
    (
        150*(this->mu()/this->rho_)/(omegaU_*sqr(y_)),
        scalar(10)
    );

    return 1 - tanh(pow4(arg3));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTPANS<BasicTurbulenceModel>::kOmegaSSTPANS::F23() const
{
    tmp<volScalarField> f23(F2());

    if (F3_)
    {
        f23() *= F3();
    }

    return f23;
}


template<class BasicTurbulenceModel>
void kOmegaSSTPANS<BasicTurbulenceModel>::correctNut(const volScalarField& S2)
{
    this->nut_ = a1_*kU_/max(a1_*omegaU_, b1_*F23()*sqrt(S2));
    this->nut_.correctBoundaryConditions();

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kOmegaSSTPANS<BasicTurbulenceModel>::correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))));
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
    eddyViscosity<RASModel<BasicTurbulenceModel> >
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    alphaK1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK1",
            this->coeffDict_,
            0.85
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK2",
            this->coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega1",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega2",
            this->coeffDict_,
            0.856
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma1",
            this->coeffDict_,
            5.0/9.0
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma2",
            this->coeffDict_,
            0.44
        )
    ),
    beta1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta1",
            this->coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta2",
            this->coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a1",
            this->coeffDict_,
            0.31
        )
    ),
    b1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "b1",
            this->coeffDict_,
            1.0
        )
    ),
    c1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "c1",
            this->coeffDict_,
            10.0
        )
    ),
    F3_
    (
        Switch::lookupOrAddToDict
        (
            "F3",
            this->coeffDict_,
            false
        )
    ),

    y_(wallDist::New(this->mesh_).y()),

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
    uLimVec
    (
        dimensionedScalar("uLimVec", uLim_)
    ),
    loLimVec
    (
        dimensionedScalar("loLimVec", loLim_)
    ),
    cellVolume
    (
        IOobject
        (
            "cellVolume",
            this->runTime_.timeName(),
            this->mesh_
        ),
        this->mesh_,
        dimensionedScalar("zero", dimVolume, 0.0)
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

    k_
    (
        IOobject
        (
            IOobject::groupName("k", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    kU_
    (
        IOobject
        (
            IOobject::groupName("kU", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omegaU_
    (
        IOobject
        (
            IOobject::groupName("omegaU", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
    //Initialize variable cellVolume
    cellVolume.internalField() = this->mesh_.V();

    kU_ = k_*fK_;
    omegaU_ = omega_*fOmega_;

    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);
    bound(kU_, this->kMin_);
    bound(omegaU_, this->omegaMin_);

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
    if (eddyViscosity<RASModel<BasicTurbulenceModel> >::read())
    {
        alphaK1_.readIfPresent(this->coeffDict());
        alphaK2_.readIfPresent(this->coeffDict());
        alphaOmega1_.readIfPresent(this->coeffDict());
        alphaOmega2_.readIfPresent(this->coeffDict());
        gamma1_.readIfPresent(this->coeffDict());
        gamma2_.readIfPresent(this->coeffDict());
        beta1_.readIfPresent(this->coeffDict());
        beta2_.readIfPresent(this->coeffDict());
        betaStar_.readIfPresent(this->coeffDict());
        a1_.readIfPresent(this->coeffDict());
        b1_.readIfPresent(this->coeffDict());
        c1_.readIfPresent(this->coeffDict());
        F3_.readIfPresent("F3", this->coeffDict());
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

    eddyViscosity<RASModel<BasicTurbulenceModel> >::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField S2(2*magSqr(symm(tgradU())));
    volScalarField GbyNu((tgradU() && dev(twoSymm(tgradU()))));
    volScalarField G(this->GName(), nut*GbyNu);
    tgradU.clear();

    // Update omega and G at the wall
    omegaU_.boundaryField().updateCoeffs();

    volScalarField CDkOmega
    (
        (2*alphaOmega2_*(fK_/fOmega_))*(fvc::grad(kU_) & fvc::grad(omegaU_))
        /omegaU_
    );

    volScalarField F1(this->F1(CDkOmega));

    {
        volScalarField gamma(this->gamma(F1));
        volScalarField beta(this->beta(F1));
        volScalarField betaL
        (
            gamma*betaStar_ - (gamma *betaStar_/fOmega_) +(beta/fOmega_)
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
                (c1_/a1_)*betaStar_*omegaU_*max(a1_*omegaU_, b1_*F23()*sqrt(S2))
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
        );

        omegaUEqn().relax();

        omegaUEqn().boundaryManipulate(omegaU_.boundaryField());

        solve(omegaUEqn);
        bound(omegaU_, min(fOmega_)*this->omegaMin_);
    }

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kUEqn
    (
        fvm::ddt(alpha, rho, kU_)
      + fvm::div(alphaRhoPhi, kU_)
      - fvm::laplacian(alpha*rho*DkUEff(F1), kU_)
     ==
        min(alpha*rho*G, (c1_*betaStar_)*alpha*rho*kU_*omegaU_)
      - fvm::SuSp((2.0/3.0)*alpha*rho*divU, kU_)
      - fvm::Sp(alpha*rho*betaStar_*omegaU_, kU_)
      + kSource()
    );

    kUEqn().relax();
    solve(kUEqn);
    bound(kU_, min(fK_)*this->kMin_);

    // Calculation of Turbulent kinetic energy and Frequency
    k_ = kU_/fK_;
    omega_ = omegaU_/fOmega_;

    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    correctNut(S2);
    
    //Info << "Recalculating fK with new kU and omegaU" << endl;

    // Recalculate fK with new kU and epsilonU
    correctPANSCoeffs();

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
