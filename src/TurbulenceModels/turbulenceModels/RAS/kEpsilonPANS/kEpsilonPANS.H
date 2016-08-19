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

Class
    Foam::RASModels::kEpsilonPANS

Group
    grpRASTurbulence

Description
    PANS implementation based on the Standard k-epsilon turbulence model for 
    incompressible and compressible flows including rapid distortion theory 
    (RDT) based compression term. The fK parameter, responsible for switching 
    from RAS to DNS is calculated based on [Girimaji,2005]:

    fK_ = min( max( sqrt(Cmu_.value())*(pow(pow(cellVolume,1.0/3.0)/
           (pow(kU_,1.5)/epsilonU_),2.0/3.0)) , loLimVec ), uLimVec );

    bounded between loLim (DNS -> 0.0) and uLim (RAS -> 1.0). The lower the 
    uLim, defined in turbulenceProperties dict, more scales are resolved.


    Reference:
    \verbatim
        Standard model:
            Launder, B. E., & Spalding, D. B. (1972).
            Lectures in mathematical models of turbulence.

            Launder, B. E., & Spalding, D. B. (1974).
            The numerical computation of turbulent flows.
            Computer methods in applied mechanics and engineering,
            3(2), 269-289.

        For the RDT-based compression term:
            El Tahry, S. H. (1983).
            k-epsilon equation for compressible reciprocating engine flows.
            Journal of Energy, 7(4), 345-353.

        PANS:
            Girimaji, S. S. & Abdol-Hamid K. S. (2005).
            Partially-averaged Navier-Stokes model for turbulence: 
            Implementation and Validation.
            43rd AIAA Aerospace Science Meeting and Exhibit.

            Girimaji, S. S. (2006).
            Partially-averaged Navier-Stokes method for turbulence: A Reynolds
            averaged Navier-Stokes to direct numerical simulation bridging
            method.
            Journal of Applied Mechanics, Vol 73, 413-421.

            Girimaji, S. & Jeong, E. & Srinivasan, R. (2006).
            Partially-averaged Navier-Stokes method for turbulence: Fixed point
            analysis and comparison with unsteady partially averaged 
            Navier-Stokes.
            Journal of Applied Mechanics, Vol 73, 422-429.

    \endverbatim

    The default model coefficients are
    \verbatim
        kEpsilonPANSCoeffs
        {
            Cmu         0.09;
            C1          1.44;
            C2          1.92;
            C3          -0.33;
            sigmak      1.0;
            sigmaEps    1.3;
            fEpsilon    1.0;
            fKupperLimit  1.0;
            fKlowerLimit  0.1;
        }
    \endverbatim

SourceFiles
    kEpsilonPANS.C

\*---------------------------------------------------------------------------*/

#ifndef kEpsilonPANS_H
#define kEpsilonPANS_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class kEpsilonPANS Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kEpsilonPANS
:
    public eddyViscosity<RASModel<BasicTurbulenceModel> >
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        kEpsilonPANS(const kEpsilonPANS&);
        kEpsilonPANS& operator=(const kEpsilonPANS&);


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar C3_;
            dimensionedScalar sigmak_;
            dimensionedScalar sigmaEps_;
            dimensionedScalar fEpsilon_;
            dimensionedScalar uLim_;
            dimensionedScalar loLim_;
            dimensionedScalar uLimVec;
            dimensionedScalar loLimVec;

        // Fields

            volScalarField cellVolume;
            volScalarField fK_;
            volScalarField C2U;
            volScalarField k_;
            volScalarField kU_;
            volScalarField epsilon_;
            volScalarField epsilonU_;


    // Protected Member Functions

        virtual void correctNut();
        virtual void correctPANSCoeffs();
        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> epsilonSource() const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kEpsilonPANS");


    // Constructors

        //- Construct from components
        kEpsilonPANS
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
    virtual ~kEpsilonPANS()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for unresolved k
        tmp<volScalarField> DkUEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DkUEff",
                    (this->nut_/(fK_*fK_*sigmak_/fEpsilon_)
                     + this->nu())
                )
            );
        }

        //- Return the effective diffusivity for unresolved epsilon
        tmp<volScalarField> DepsilonUEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DepsilonUEff",
                    (this->nut_/(fK_*fK_*sigmaEps_/fEpsilon_)
                     + this->nu())
                )
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "kEpsilonPANS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
