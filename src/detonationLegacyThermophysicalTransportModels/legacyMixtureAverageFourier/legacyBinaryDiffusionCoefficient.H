/*---------------------------------------------------------------------------*\
  legacyBinaryDiffusionCoefficient
  OpenFOAM Foundation 14 Function2 reproducing the detonationFoam OF8
  binary-diffusion fit used by solverTypeNS_mixtureAverage.
\*---------------------------------------------------------------------------*/
#ifndef legacyBinaryDiffusionCoefficient_H
#define legacyBinaryDiffusionCoefficient_H

#include "Function2.H"
#include "FixedList.H"

namespace Foam
{
namespace Function2s
{

class legacyBinaryDiffusionCoefficient
:
    public FieldFunction2<scalar, legacyBinaryDiffusionCoefficient>
{
    // Diff1..Diff4 from the legacy constant/binaryDiff dictionary.
    FixedList<scalar, 4> coeffs_;

public:

    TypeName("legacyBinaryDiffusionCoefficient");

    legacyBinaryDiffusionCoefficient
    (
        const word& name,
        const unitSets& units,
        const dictionary& dict
    );

    virtual tmp<Function2<scalar>> clone() const
    {
        return tmp<Function2<scalar>>
        (
            new legacyBinaryDiffusionCoefficient(*this)
        );
    }

    // Exact published OF8 fit, with p in Pa and T in K:
    // D = 1e-4*exp(c0 + ln(T)*(c1 + ln(T)*(c2 + ln(T)*c3)))
    //     /(p/101325)
    virtual scalar value(const scalar p, const scalar T) const;

    virtual void write(Ostream& os, const unitSets& units) const;
};

} // End namespace Function2s
} // End namespace Foam

#endif

// ************************************************************************* //
