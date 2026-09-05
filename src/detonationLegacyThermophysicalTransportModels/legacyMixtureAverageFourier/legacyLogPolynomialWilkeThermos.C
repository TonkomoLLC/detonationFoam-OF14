/*---------------------------------------------------------------------------*\
  Register the one OpenFOAM-14 thermo package combination needed by the
  detonationFoam legacy NS_mixtureAverage compatibility layer:

    coefficientWilkeMulticomponentMixture
    + logPolynomialTransport<...,8>
    + janafThermo
    + sensibleInternalEnergy
    + perfectGas

  All component models are native OpenFOAM-14 classes.  This file only adds
  the missing run-time registration combination; it does not duplicate them.
\*---------------------------------------------------------------------------*/

#include "psiMulticomponentThermo.H"
#include "coefficientWilkeMulticomponentMixture.H"
#include "logPolynomialTransport.H"
#include "janafThermo.H"
#include "sensibleInternalEnergy.H"
#include "perfectGas.H"
#include "specie.H"
#include "thermo.H"
#include "forThermo.H"
#include "makeFluidMulticomponentThermo.H"

namespace Foam
{

forThermo
(
    logPolynomialTransport,
    sensibleInternalEnergy,
    janafThermo,
    perfectGas,
    specie,
    makeFluidMulticomponentThermos,
    psiThermo,
    psiMulticomponentThermo,
    coefficientWilkeMulticomponentMixture
);

} // End namespace Foam

// ************************************************************************* //
