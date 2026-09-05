/*---------------------------------------------------------------------------*\
  legacyBinaryDiffusionCoefficient
\*---------------------------------------------------------------------------*/
#include "legacyBinaryDiffusionCoefficient.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace Function2s
{
    addScalarFunction2(legacyBinaryDiffusionCoefficient);
}
}

Foam::Function2s::legacyBinaryDiffusionCoefficient::
legacyBinaryDiffusionCoefficient
(
    const word& name,
    const unitSets&,
    const dictionary& dict
)
:
    FieldFunction2<scalar, legacyBinaryDiffusionCoefficient>(name),
    coeffs_(dict.lookup("coeffs"))
{}

Foam::scalar
Foam::Function2s::legacyBinaryDiffusionCoefficient::value
(
    const scalar p,
    const scalar T
) const
{
    const scalar lnT = log(T);

    return
        1e-4
       *exp
        (
            coeffs_[0]
          + lnT
           *(
                coeffs_[1]
              + lnT*(coeffs_[2] + lnT*coeffs_[3])
            )
        )
       /(p/101325.0);
}

void Foam::Function2s::legacyBinaryDiffusionCoefficient::write
(
    Ostream& os,
    const unitSets&
) const
{
    writeEntry(os, "coeffs", coeffs_);
}

// ************************************************************************* //
