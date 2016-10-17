/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenFOAM Foundation
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

Class
    Foam::RASModels::kOmegaSSTPANS

Group
    grpLESTurbulence

Description
    Partially-Average Navier-Stokes model based on the k-omega-SST RAS model.

    References:
    \verbatin
        Luo, D.; Yan, C.; Liu, H. & Zhao, R. (2014).
        Comparative assessment of PANS and DES for simulation of flow past a
        circular cylinder.
        Journal of wind engineering and industrial aerodynamics, Vol 134, 65-77
    \endverbatim

    The model coefficients are
    \verbatim
        kOmegaSSTPANSCoeffs
        {
            // Default SST coefficients
            alphaK1     0.85;
            alphaK2     1.0;
            alphaOmega1 0.5;
            alphaOmega2 0.856;
            beta1       0.075;
            beta2       0.0828;
            betaStar    0.09;
            gamma1      5/9;
            gamma2      0.44;
            a1          0.31;
            b1          1.0;
            c1          10.0;
            F3          no;

            // Default PANS coefficients
            Cs          0.11;
            kappa       0.41;
            zeta2       3.51;
            sigmaPhi    2.0/3.0;
            C           2;

            fEpsilon    1.0;
            fKupperLimit  1.0;
            fKlowerLimit  0.1;

            // Delta must be specified for SAS e.g.
            delta cubeRootVol;

            cubeRootVolCoeffs
            {}
        }
    \endverbatim

SourceFiles
    kOmegaSSTPANS.C

\*---------------------------------------------------------------------------*/

#ifndef kOmegaSSTPANS_H
#define kOmegaSSTPANS_H

#include "kOmegaSST.H"
#include "LESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                     Class kOmegaSSTPANS Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kOmegaSSTPANS
:
    public kOmegaSST<BasicTurbulenceModel>
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        kOmegaSSTPANS(const kOmegaSSTPANS&);
        void operator=(const kOmegaSSTPANS&);


protected:

    // Protected data

        // Model constants

            dimensionedScalar fEpsilon_;
            dimensionedScalar uLim_;
            dimensionedScalar loLim_;


        // Fields

            volScalarField fK_;
            volScalarField fOmega_;
            //- Run-time selectable delta model
            autoPtr<Foam::LESdelta> delta_;
            volScalarField kU_;
            volScalarField omegaU_;

    // Private Member Functions   
    
        tmp<volScalarField> F1(const volScalarField& CDkOmega) const;
        tmp<volScalarField> F2() const;
        tmp<volScalarField> F3() const;

        void correctNut(const volScalarField& S2, const volScalarField& F2);
    
    // Protected Member Functions

        virtual void correctNut();

        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> omegaSource() const;
        virtual tmp<fvScalarMatrix> Qsas
        (
            const volScalarField& S2,
            const volScalarField& gamma,
            const volScalarField& beta
        ) const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kOmegaSSTPANS");


    // Constructors

        //- Construct from components
        kOmegaSSTPANS
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~kOmegaSSTPANS()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Access function to filter width
        inline const volScalarField& delta() const
        {
            return delta_();
        }

        //- Return the effective diffusivity for unresolved k
        tmp<volScalarField> DkUEff(const volScalarField& F1) const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DkUEff",
                    (fK_/fOmega_)*this->alphaK(F1)*this->nut_ + this->nu()
                )
            );
        }

        //- Return the effective diffusivity for unresolved omega
        tmp<volScalarField> DomegaUEff(const volScalarField& F1) const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DomegaUEff",
                    (fK_/fOmega_)*this->alphaOmega(F1)*this->nut_ + this->nu()
                )
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return this->k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
              return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "epsilon",
                        this->mesh_.time().timeName(),
                        this->mesh_
                    ),
                    this->betaStar_*this->k_*this->omega_,
                    this->omega_.boundaryField().types()
                )
            );
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return this->omega_;
        }


        //- Return the unresolved turbulence kinetic energy
        virtual tmp<volScalarField> kU() const
        {
            return kU_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> omegaU() const
        {
            return omegaU_;
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "kOmegaSSTPANS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
