/*---------------------------------------------------------------------------*\
  Density-based face reconstruction/flux coefficients for detonationFluid.

  Kurganov and Tadmor follow the OpenFOAM-14 shockFluid formulation.
  Gate B1 restored HLL using the original detonationFoam wave-speed estimate.
  Gate B2 restored HLLC using the same signal speeds plus the contact wave.
  Gate B3 restored HLLCP, including its pressure-ratio sensor, low-Mach pressure
  blending, and pressure-difference mass/momentum/energy correction.
  Gate B4 restores AUSM+ with the original fourth-order Mach split and
  fifth-order pressure split (alpha=3/16, beta=1/8).
  Gate B5 restores AUSM+up, including its all-speed pressure and velocity
  corrections (Kp=0.25, Ku=0.75, sigma=1, fa=1).

  HLLC/HLLCP/AUSM+/AUSM+up are represented by conservative split advection plus explicit
  face pressure/work and (for HLLCP) correction fluxes.  This keeps the OF14
  foamRun module, chemistry, transport and mesh-update path unified.
\*---------------------------------------------------------------------------*/

#include "detonationFluid.H"

void Foam::solvers::detonationFluid::fluxPredictor()
{
    if (!pos_.valid())
    {
        pos_ = surfaceScalarField::New
        (
            "pos",
            mesh,
            dimensionedScalar(dimless, 1)
        );
        neg_ = surfaceScalarField::New
        (
            "neg",
            mesh,
            dimensionedScalar(dimless, -1)
        );
    }

    rhoPos_ = interpolate(rho_, pos_());
    rhoNeg_ = interpolate(rho_, neg_());

    const volVectorField rhoU(rho_*U_);
    rhoUPos_ = interpolate(rhoU, pos_(), U_.name());
    rhoUNeg_ = interpolate(rhoU, neg_(), U_.name());

    UPos_ = surfaceVectorField::New("U_pos", rhoUPos_()/rhoPos_());
    UNeg_ = surfaceVectorField::New("U_neg", rhoUNeg_()/rhoNeg_());

    const volScalarField& T = thermo_.T();
    const volScalarField rPsi("rPsi", 1.0/thermo_.psi());
    const surfaceScalarField rPsiPos(interpolate(rPsi, pos_(), T.name()));
    const surfaceScalarField rPsiNeg(interpolate(rPsi, neg_(), T.name()));

    pPos_ = surfaceScalarField::New("p_pos", rhoPos_()*rPsiPos);
    pNeg_ = surfaceScalarField::New("p_neg", rhoNeg_()*rPsiNeg);

    surfaceScalarField phivPos("phiv_pos", UPos_() & mesh.Sf());
    surfaceScalarField phivNeg("phiv_neg", UNeg_() & mesh.Sf());

    if (mesh.moving())
    {
        phivPos -= mesh.phi();
        phivNeg -= mesh.phi();
    }

    const volScalarField c("c", sqrt(thermo_.Cp()/thermo_.Cv()*rPsi));
    const surfaceScalarField cPos(interpolate(c, pos_(), T.name()));
    const surfaceScalarField cNeg(interpolate(c, neg_(), T.name()));
    const surfaceScalarField cSfPos("cSf_pos", cPos*mesh.magSf());
    const surfaceScalarField cSfNeg("cSf_neg", cNeg*mesh.magSf());

    const dimensionedScalar vZero("vZero", dimVolume/dimTime, 0);

    // Native OF14/Kurganov one-sided signal bounds are the default.
    surfaceScalarField ap
    (
        "ap",
        max(max(phivPos + cSfPos, phivNeg + cSfNeg), vZero)
    );
    surfaceScalarField am
    (
        "am",
        min(min(phivPos - cSfPos, phivNeg - cSfNeg), vZero)
    );

    // HLL and HLLC use the original detonationFoam density-square-root
    // weighted signal-speed estimate.
    tmp<surfaceScalarField> tUvPos;
    tmp<surfaceScalarField> tUvNeg;
    tmp<surfaceScalarField> tSLeft;
    tmp<surfaceScalarField> tSRight;

    if (fluxScheme_ == "HLL" || fluxScheme_ == "HLLC")
    {
        const surfaceScalarField sqrtRhoPos
        (
            "HLL_sqrtRho_pos",
            sqrt(rhoPos_())
        );
        const surfaceScalarField sqrtRhoNeg
        (
            "HLL_sqrtRho_neg",
            sqrt(rhoNeg_())
        );
        const surfaceScalarField wPos
        (
            "HLL_w_pos",
            sqrtRhoPos/(sqrtRhoPos + sqrtRhoNeg)
        );

        tUvPos = surfaceScalarField::New
        (
            "HLL_Uv_pos",
            phivPos/mesh.magSf()
        );
        tUvNeg = surfaceScalarField::New
        (
            "HLL_Uv_neg",
            phivNeg/mesh.magSf()
        );

        const surfaceScalarField cTilde
        (
            "HLL_cTilde",
            wPos*cPos + (1.0 - wPos)*cNeg
        );
        const surfaceScalarField UvTilde
        (
            "HLL_UvTilde",
            wPos*tUvPos() + (1.0 - wPos)*tUvNeg()
        );

        tSLeft = surfaceScalarField::New
        (
            "HLL_sLeft",
            min(tUvPos() - cPos, UvTilde - cTilde)
        );
        tSRight = surfaceScalarField::New
        (
            "HLL_sRight",
            max(tUvNeg() + cNeg, UvTilde + cTilde)
        );

        if (fluxScheme_ == "HLL")
        {
            ap = max(tSRight()*mesh.magSf(), vZero);
            am = min(tSLeft()*mesh.magSf(), vZero);
        }
    }

    if (fluxScheme_ == "HLLC")
    {
        // HLLC cannot in general be represented by the single HLL
        // dissipation coefficient.  Build its four standard branches directly
        // into the same cached split-advection/pressure-work fields consumed
        // by the OF14 momentum, species and energy equations.
        aPos_ = surfaceScalarField::New
        (
            "a_pos",
            mesh,
            dimensionedScalar(dimless, 0)
        );
        aNeg_ = surfaceScalarField::New
        (
            "a_neg",
            mesh,
            dimensionedScalar(dimless, 0)
        );
        // Preserve the orientation metadata carried by reconstructed face
        // volumetric fluxes.  This is important on modern OpenFOAM releases:
        // constructing these fields from dimensioned zero would make them
        // unoriented even though they are later used as oriented fluxes.
        aSf_ = surfaceScalarField::New("aSf", 0.0*phivPos);
        aphivPos_ = surfaceScalarField::New("aphiv_pos", 0.0*phivPos);
        aphivNeg_ = surfaceScalarField::New("aphiv_neg", 0.0*phivNeg);
        pFlux_ = surfaceScalarField::New
        (
            "pFlux",
            mesh,
            dimensionedScalar(dimPressure, 0)
        );
        pEnergyFlux_ = surfaceScalarField::New
        (
            "pEnergyFlux",
            0.0*pPos_()*phivPos
        );
        pWork_ = surfaceScalarField::New
        (
            "pWork",
            mesh,
            dimensionedScalar(dimPressure, 0)
        );
        Uf_ = surfaceVectorField::New
        (
            "Uf",
            mesh,
            dimensionedVector(dimVelocity, Zero)
        );
        // CFL signal magnitude is intentionally unoriented.
        amaxSf_ = surfaceScalarField::New("amaxSf", 0.0*mag(phivPos));

        auto setHLLCFace = []
        (
            const scalar rhoL,
            const scalar rhoR,
            const vector& UL,
            const vector& UR,
            const scalar pL,
            const scalar pR,
            const scalar cL,
            const scalar cR,
            const vector& Sf,
            const scalar phivL,
            const scalar phivR,
            scalar& aP,
            scalar& aN,
            scalar& aS,
            scalar& qP,
            scalar& qN,
            scalar& pF,
            scalar& pEF,
            scalar& pW,
            vector& Uf,
            scalar& amax
        )
        {
            const scalar area = mag(Sf);
            const vector n = Sf/area;
            const scalar uL = phivL/area;
            const scalar uR = phivR/area;

            const scalar sqrtRhoL = sqrt(rhoL);
            const scalar sqrtRhoR = sqrt(rhoR);
            const scalar wL = sqrtRhoL/(sqrtRhoL + sqrtRhoR);
            const scalar uTilde = wL*uL + (1.0 - wL)*uR;
            const scalar cTilde = wL*cL + (1.0 - wL)*cR;

            const scalar SL = min(uL - cL, uTilde - cTilde);
            const scalar SR = max(uR + cR, uTilde + cTilde);

            amax = max(mag(SL), mag(SR))*area;

            auto setHLLFallback = [&]()
            {
                const scalar apf = max(SR*area, scalar(0));
                const scalar amf = min(SL*area, scalar(0));
                const scalar den = apf - amf;

                aP = apf/den;
                aN = 1.0 - aP;
                aS = amf*aP;
                qP = aP*phivL - aS;
                qN = aN*phivR + aS;
                pF = aP*pL + aN*pR;
                pEF = qP*pL + qN*pR + aS*(pL - pR);
                pW = pF;
                Uf = aP*UL + aN*UR;
            };

            const scalar starDen =
                rhoL*(SL - uL) - rhoR*(SR - uR);

            if (mag(starDen) < VSMALL)
            {
                setHLLFallback();
                return;
            }

            const scalar SStar =
            (
                pR - pL
              + rhoL*uL*(SL - uL)
              - rhoR*uR*(SR - uR)
            )/starDen;

            const scalar pStarL = pL + rhoL*(SL - uL)*(SStar - uL);
            const scalar pStarR = pR + rhoR*(SR - uR)*(SStar - uR);

            if
            (
                (SStar > 0 && mag(SL - SStar) < VSMALL)
             || (SStar <= 0 && SR > 0 && mag(SR - SStar) < VSMALL)
            )
            {
                setHLLFallback();
                return;
            }

            aS = 0.0;

            if (SL > 0)
            {
                aP = 1.0;
                aN = 0.0;
                qP = phivL;
                qN = 0.0;
                pF = pL;
                pEF = pL*phivL;
                pW = pL;
                Uf = UL;
            }
            else if (SStar > 0)
            {
                const scalar dS = SL - SStar;
                aP = 1.0;
                aN = 0.0;
                qP = SStar*(SL - uL)/dS*area;
                qN = 0.0;
                pF = (SL*pStarL - SStar*pL)/dS;
                pEF = SStar*(SL*pStarL - uL*pL)/dS*area;
                pW = 0.5*(pStarL + pStarR);
                Uf = UL + (pStarL - pL)/(rhoL*(SL - uL))*n;
            }
            else if (SR > 0)
            {
                const scalar dS = SR - SStar;
                aP = 0.0;
                aN = 1.0;
                qP = 0.0;
                qN = SStar*(SR - uR)/dS*area;
                pF = (SR*pStarR - SStar*pR)/dS;
                pEF = SStar*(SR*pStarR - uR*pR)/dS*area;
                pW = 0.5*(pStarL + pStarR);
                Uf = UR + (pStarR - pR)/(rhoR*(SR - uR))*n;
            }
            else
            {
                aP = 0.0;
                aN = 1.0;
                qP = 0.0;
                qN = phivR;
                pF = pR;
                pEF = pR*phivR;
                pW = pR;
                Uf = UR;
            }
        };

        scalarField& aP = aPos_.ref().primitiveFieldRef();
        scalarField& aN = aNeg_.ref().primitiveFieldRef();
        scalarField& aS = aSf_.ref().primitiveFieldRef();
        scalarField& qP = aphivPos_.ref().primitiveFieldRef();
        scalarField& qN = aphivNeg_.ref().primitiveFieldRef();
        scalarField& pF = pFlux_.ref().primitiveFieldRef();
        scalarField& pEF = pEnergyFlux_.ref().primitiveFieldRef();
        scalarField& pW = pWork_.ref().primitiveFieldRef();
        vectorField& Uf = Uf_.ref().primitiveFieldRef();
        scalarField& amax = amaxSf_.ref().primitiveFieldRef();

        const scalarField& rhoP = rhoPos_().primitiveField();
        const scalarField& rhoN = rhoNeg_().primitiveField();
        const vectorField& UP = UPos_().primitiveField();
        const vectorField& UN = UNeg_().primitiveField();
        const scalarField& pP = pPos_().primitiveField();
        const scalarField& pN = pNeg_().primitiveField();
        const scalarField& cP = cPos.primitiveField();
        const scalarField& cN = cNeg.primitiveField();
        const vectorField& Sf = mesh.Sf().primitiveField();
        const scalarField& phP = phivPos.primitiveField();
        const scalarField& phN = phivNeg.primitiveField();

        forAll(aP, facei)
        {
            setHLLCFace
            (
                rhoP[facei], rhoN[facei],
                UP[facei], UN[facei],
                pP[facei], pN[facei],
                cP[facei], cN[facei],
                Sf[facei], phP[facei], phN[facei],
                aP[facei], aN[facei], aS[facei],
                qP[facei], qN[facei], pF[facei], pEF[facei], pW[facei],
                Uf[facei], amax[facei]
            );
        }

        forAll(mesh.boundary(), patchi)
        {
            fvsPatchScalarField& aPp = aPos_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& aNp = aNeg_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& aSp = aSf_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& qPp = aphivPos_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& qNp = aphivNeg_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& pFp = pFlux_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& pEFp = pEnergyFlux_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& pWp = pWork_.ref().boundaryFieldRef()[patchi];
            fvsPatchVectorField& Ufp = Uf_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& amaxp = amaxSf_.ref().boundaryFieldRef()[patchi];

            const fvsPatchScalarField& rhoPp = rhoPos_().boundaryField()[patchi];
            const fvsPatchScalarField& rhoNp = rhoNeg_().boundaryField()[patchi];
            const fvsPatchVectorField& UPp = UPos_().boundaryField()[patchi];
            const fvsPatchVectorField& UNp = UNeg_().boundaryField()[patchi];
            const fvsPatchScalarField& pPp = pPos_().boundaryField()[patchi];
            const fvsPatchScalarField& pNp = pNeg_().boundaryField()[patchi];
            const fvsPatchScalarField& cPp = cPos.boundaryField()[patchi];
            const fvsPatchScalarField& cNp = cNeg.boundaryField()[patchi];
            const fvsPatchVectorField& Sfp = mesh.Sf().boundaryField()[patchi];
            const fvsPatchScalarField& phPp = phivPos.boundaryField()[patchi];
            const fvsPatchScalarField& phNp = phivNeg.boundaryField()[patchi];

            forAll(aPp, facei)
            {
                setHLLCFace
                (
                    rhoPp[facei], rhoNp[facei],
                    UPp[facei], UNp[facei],
                    pPp[facei], pNp[facei],
                    cPp[facei], cNp[facei],
                    Sfp[facei], phPp[facei], phNp[facei],
                    aPp[facei], aNp[facei], aSp[facei],
                    qPp[facei], qNp[facei], pFp[facei], pEFp[facei],
                    pWp[facei], Ufp[facei], amaxp[facei]
                );
            }
        }

        phi_ = aphivPos_()*rhoPos_() + aphivNeg_()*rhoNeg_();
        return;
    }

    if (fluxScheme_ == "HLLCP")
    {
        // Original HLLCP pressure sensor.  Each cell receives the minimum
        // neighbour pressure ratio, the cell sensor is cubed, and the result
        // is interpolated to faces.  f=1 recovers the HLLC-like star pressure;
        // smaller f activates the pressure correction around strong jumps.
        volScalarField fCells
        (
            IOobject
            (
                "HLLCP_fCells",
                runTime.name(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar(dimless, 1.0)
        );

        scalarField& fCell = fCells.primitiveFieldRef();
        const scalarField& pCell = p_.primitiveField();
        const labelListList& cellCells = mesh.cellCells();

        forAll(fCell, celli)
        {
            forAll(cellCells[celli], cellj)
            {
                const label nbr = cellCells[celli][cellj];
                fCell[celli] = min
                (
                    fCell[celli],
                    min(pCell[celli]/pCell[nbr], pCell[nbr]/pCell[celli])
                );
            }
        }

        forAll(fCell, celli)
        {
            const scalar f = fCell[celli];
            fCell[celli] = f*f*f;
        }
        fCells.correctBoundaryConditions();

        const tmp<surfaceScalarField> tfFace = fvc::interpolate(fCells);
        const surfaceScalarField& fFace = tfFace();

        aPos_ = surfaceScalarField::New
        (
            "a_pos",
            mesh,
            dimensionedScalar(dimless, 0)
        );
        aNeg_ = surfaceScalarField::New
        (
            "a_neg",
            mesh,
            dimensionedScalar(dimless, 0)
        );
        aSf_ = surfaceScalarField::New("aSf", 0.0*phivPos);
        aphivPos_ = surfaceScalarField::New("aphiv_pos", 0.0*phivPos);
        aphivNeg_ = surfaceScalarField::New("aphiv_neg", 0.0*phivNeg);
        pFlux_ = surfaceScalarField::New
        (
            "pFlux",
            mesh,
            dimensionedScalar(dimPressure, 0)
        );
        pEnergyFlux_ = surfaceScalarField::New
        (
            "pEnergyFlux",
            0.0*pPos_()*phivPos
        );
        pWork_ = surfaceScalarField::New
        (
            "pWork",
            mesh,
            dimensionedScalar(dimPressure, 0)
        );
        Uf_ = surfaceVectorField::New
        (
            "Uf",
            mesh,
            dimensionedVector(dimVelocity, Zero)
        );
        amaxSf_ = surfaceScalarField::New("amaxSf", 0.0*mag(phivPos));

        // HLLCP adds a pressure-difference mass flux that is not expressible
        // as qL*rhoL + qR*rhoR.  Preserve it, and its corresponding tangential
        // momentum and kinetic-energy transport, as explicit conservative
        // correction fluxes.
        rhoPhiCorrection_ = surfaceScalarField::New
        (
            "rhoPhiCorrection",
            0.0*(rhoPos_()*phivPos)
        );
        rhoUPhiCorrection_ = surfaceVectorField::New
        (
            "rhoUPhiCorrection",
            0.0*(phivPos*rhoUPos_())
        );
        rhoEPhiCorrection_ = surfaceScalarField::New
        (
            "rhoEPhiCorrection",
            0.0*pPos_()*phivPos
        );

        auto setHLLCPFace = []
        (
            const scalar rhoL,
            const scalar rhoR,
            const vector& UL,
            const vector& UR,
            const scalar pL,
            const scalar pR,
            const scalar cL,
            const scalar cR,
            const scalar fSensor,
            const vector& Sf,
            const scalar phivL,
            const scalar phivR,
            scalar& aP,
            scalar& aN,
            scalar& aS,
            scalar& qP,
            scalar& qN,
            scalar& pF,
            scalar& pEF,
            scalar& pW,
            vector& Uf,
            scalar& amax,
            scalar& mCorr,
            vector& momCorr,
            scalar& eCorr
        )
        {
            const scalar area = mag(Sf);
            const vector n = Sf/area;
            const scalar uL = phivL/area;
            const scalar uR = phivR/area;

            const scalar sqrtRhoL = sqrt(rhoL);
            const scalar sqrtRhoR = sqrt(rhoR);
            const scalar wL = sqrtRhoL/(sqrtRhoL + sqrtRhoR);
            const scalar wR = 1.0 - wL;
            const vector UTilde = wL*UL + wR*UR;
            const scalar uTilde = wL*uL + wR*uR;
            const scalar cTilde = wL*cL + wR*cR;
            const scalar MaL = uL/cL;
            const scalar MaR = uR/cR;
            const scalar MaTilde = wL*MaL + wR*MaR;

            const scalar SL = min(uL - cL, uTilde - cTilde);
            const scalar SR = max(uR + cR, uTilde + cTilde);
            amax = max(mag(SL), mag(SR))*area;

            mCorr = 0.0;
            momCorr = Zero;
            eCorr = 0.0;

            auto setHLLFallback = [&]()
            {
                const scalar apf = max(SR*area, scalar(0));
                const scalar amf = min(SL*area, scalar(0));
                const scalar den = apf - amf;

                aP = apf/den;
                aN = 1.0 - aP;
                aS = amf*aP;
                qP = aP*phivL - aS;
                qN = aN*phivR + aS;
                pF = aP*pL + aN*pR;
                pEF = qP*pL + qN*pR + aS*(pL - pR);
                pW = pF;
                Uf = aP*UL + aN*UR;
                mCorr = 0.0;
                momCorr = Zero;
                eCorr = 0.0;
            };

            const scalar aL = rhoL*(SL - uL);
            const scalar aR = rhoR*(SR - uR);
            const scalar starDen = aR - aL;

            if
            (
                mag(starDen) < VSMALL
             || mag(SR - SL) < VSMALL
             || mag(cTilde) < VSMALL
            )
            {
                setHLLFallback();
                return;
            }

            // Original HLLCP contact speed and common star pressure.
            const scalar SStar =
                (aR*uR - aL*uL + pL - pR)/starDen;
            const scalar pStar =
                (aR*pL - aL*pR - aL*aR*(uL - uR))/starDen;

            const scalar theta = min(max(mag(MaL), mag(MaR)), scalar(1));
            const scalar pAvg = 0.5*(pL + pR);
            const scalar pStarStar = theta*pStar + (1.0 - theta)*pAvg;
            const scalar pStarStarStar =
                fSensor*pStarStar + (1.0 - fSensor)*pStar;

            const scalar phip =
                (fSensor - 1.0)
               *SL*SR/(SR - SL)
               /(1.0 + mag(MaTilde))
               *(pR - pL)/(cTilde*cTilde);

            if
            (
                (SStar > 0 && mag(SL - SStar) < VSMALL)
             || (SStar <= 0 && SR > 0 && mag(SR - SStar) < VSMALL)
            )
            {
                setHLLFallback();
                return;
            }

            aS = 0.0;

            if (SL > 0)
            {
                aP = 1.0;
                aN = 0.0;
                qP = phivL;
                qN = 0.0;
                pF = pL;
                pEF = pL*phivL;
                pW = pL;
                Uf = UL;
            }
            else if (SStar > 0)
            {
                const scalar dS = SL - SStar;
                aP = 1.0;
                aN = 0.0;
                qP = SStar*(SL - uL)/dS*area;
                qN = 0.0;
                pF = (SL*pStarStarStar - SStar*pL)/dS;
                pEF = SStar*(SL*pStar - uL*pL)/dS*area;
                pW = pStar;
                Uf = UL
                   + (pStarStarStar - pL)/(rhoL*(SL - uL))*n;
                mCorr = phip*area;
                momCorr = mCorr*UTilde;
                eCorr = 0.5*mCorr*magSqr(UTilde);
            }
            else if (SR > 0)
            {
                const scalar dS = SR - SStar;
                aP = 0.0;
                aN = 1.0;
                qP = 0.0;
                qN = SStar*(SR - uR)/dS*area;
                pF = (SR*pStarStarStar - SStar*pR)/dS;
                pEF = SStar*(SR*pStar - uR*pR)/dS*area;
                pW = pStar;
                Uf = UR
                   + (pStarStarStar - pR)/(rhoR*(SR - uR))*n;
                mCorr = phip*area;
                momCorr = mCorr*UTilde;
                eCorr = 0.5*mCorr*magSqr(UTilde);
            }
            else
            {
                aP = 0.0;
                aN = 1.0;
                qP = 0.0;
                qN = phivR;
                pF = pR;
                pEF = pR*phivR;
                pW = pR;
                Uf = UR;
            }
        };

        scalarField& aP = aPos_.ref().primitiveFieldRef();
        scalarField& aN = aNeg_.ref().primitiveFieldRef();
        scalarField& aS = aSf_.ref().primitiveFieldRef();
        scalarField& qP = aphivPos_.ref().primitiveFieldRef();
        scalarField& qN = aphivNeg_.ref().primitiveFieldRef();
        scalarField& pF = pFlux_.ref().primitiveFieldRef();
        scalarField& pEF = pEnergyFlux_.ref().primitiveFieldRef();
        scalarField& pW = pWork_.ref().primitiveFieldRef();
        vectorField& Uf = Uf_.ref().primitiveFieldRef();
        scalarField& amax = amaxSf_.ref().primitiveFieldRef();
        scalarField& mCorr = rhoPhiCorrection_.ref().primitiveFieldRef();
        vectorField& momCorr = rhoUPhiCorrection_.ref().primitiveFieldRef();
        scalarField& eCorr = rhoEPhiCorrection_.ref().primitiveFieldRef();

        const scalarField& rhoP = rhoPos_().primitiveField();
        const scalarField& rhoN = rhoNeg_().primitiveField();
        const vectorField& UP = UPos_().primitiveField();
        const vectorField& UN = UNeg_().primitiveField();
        const scalarField& pP = pPos_().primitiveField();
        const scalarField& pN = pNeg_().primitiveField();
        const scalarField& cP = cPos.primitiveField();
        const scalarField& cN = cNeg.primitiveField();
        const scalarField& fP = fFace.primitiveField();
        const vectorField& Sf = mesh.Sf().primitiveField();
        const scalarField& phP = phivPos.primitiveField();
        const scalarField& phN = phivNeg.primitiveField();

        forAll(aP, facei)
        {
            setHLLCPFace
            (
                rhoP[facei], rhoN[facei],
                UP[facei], UN[facei],
                pP[facei], pN[facei],
                cP[facei], cN[facei], fP[facei],
                Sf[facei], phP[facei], phN[facei],
                aP[facei], aN[facei], aS[facei],
                qP[facei], qN[facei], pF[facei], pEF[facei], pW[facei],
                Uf[facei], amax[facei],
                mCorr[facei], momCorr[facei], eCorr[facei]
            );
        }

        forAll(mesh.boundary(), patchi)
        {
            fvsPatchScalarField& aPp = aPos_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& aNp = aNeg_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& aSp = aSf_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& qPp = aphivPos_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& qNp = aphivNeg_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& pFp = pFlux_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& pEFp = pEnergyFlux_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& pWp = pWork_.ref().boundaryFieldRef()[patchi];
            fvsPatchVectorField& Ufp = Uf_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& amaxp = amaxSf_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& mCorrp =
                rhoPhiCorrection_.ref().boundaryFieldRef()[patchi];
            fvsPatchVectorField& momCorrp =
                rhoUPhiCorrection_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& eCorrp =
                rhoEPhiCorrection_.ref().boundaryFieldRef()[patchi];

            const fvsPatchScalarField& rhoPp = rhoPos_().boundaryField()[patchi];
            const fvsPatchScalarField& rhoNp = rhoNeg_().boundaryField()[patchi];
            const fvsPatchVectorField& UPp = UPos_().boundaryField()[patchi];
            const fvsPatchVectorField& UNp = UNeg_().boundaryField()[patchi];
            const fvsPatchScalarField& pPp = pPos_().boundaryField()[patchi];
            const fvsPatchScalarField& pNp = pNeg_().boundaryField()[patchi];
            const fvsPatchScalarField& cPp = cPos.boundaryField()[patchi];
            const fvsPatchScalarField& cNp = cNeg.boundaryField()[patchi];
            const fvsPatchScalarField& fPp = fFace.boundaryField()[patchi];
            const fvsPatchVectorField& Sfp = mesh.Sf().boundaryField()[patchi];
            const fvsPatchScalarField& phPp = phivPos.boundaryField()[patchi];
            const fvsPatchScalarField& phNp = phivNeg.boundaryField()[patchi];

            forAll(aPp, facei)
            {
                setHLLCPFace
                (
                    rhoPp[facei], rhoNp[facei],
                    UPp[facei], UNp[facei],
                    pPp[facei], pNp[facei],
                    cPp[facei], cNp[facei], fPp[facei],
                    Sfp[facei], phPp[facei], phNp[facei],
                    aPp[facei], aNp[facei], aSp[facei],
                    qPp[facei], qNp[facei], pFp[facei], pEFp[facei],
                    pWp[facei], Ufp[facei], amaxp[facei],
                    mCorrp[facei], momCorrp[facei], eCorrp[facei]
                );
            }
        }

        phi_ =
            aphivPos_()*rhoPos_()
          + aphivNeg_()*rhoNeg_()
          + rhoPhiCorrection_();
        return;
    }

    if (fluxScheme_ == "AUSM+")
    {
        // Original detonationFoam AUSM+ (Liou 1996).  AUSM+ separates the
        // advective Mach-number flux from the pressure flux.  Represent the
        // selected upwind advective state with qP/qN and store the pressure
        // split explicitly, preserving the same mass flux for continuity and
        // species transport.
        aPos_ = surfaceScalarField::New
        (
            "a_pos",
            mesh,
            dimensionedScalar(dimless, 0)
        );
        aNeg_ = surfaceScalarField::New
        (
            "a_neg",
            mesh,
            dimensionedScalar(dimless, 0)
        );
        aSf_ = surfaceScalarField::New("aSf", 0.0*phivPos);
        aphivPos_ = surfaceScalarField::New("aphiv_pos", 0.0*phivPos);
        aphivNeg_ = surfaceScalarField::New("aphiv_neg", 0.0*phivNeg);
        pFlux_ = surfaceScalarField::New
        (
            "pFlux",
            mesh,
            dimensionedScalar(dimPressure, 0)
        );
        pEnergyFlux_ = surfaceScalarField::New
        (
            "pEnergyFlux",
            0.0*pPos_()*phivPos
        );
        pWork_ = surfaceScalarField::New
        (
            "pWork",
            mesh,
            dimensionedScalar(dimPressure, 0)
        );
        Uf_ = surfaceVectorField::New
        (
            "Uf",
            mesh,
            dimensionedVector(dimVelocity, Zero)
        );
        // Use the physical acoustic characteristic bound for time stepping.
        // The original AUSM+ flux itself does not require HLL-style signal
        // speeds, but the OF14 module still needs a conservative CFL estimate.
        amaxSf_ = surfaceScalarField::New("amaxSf", 0.0*mag(phivPos));

        auto setAUSMPlusFace = []
        (
            const scalar rhoL,
            const scalar rhoR,
            const vector& UL,
            const vector& UR,
            const scalar pL,
            const scalar pR,
            const scalar cL,
            const scalar cR,
            const vector& Sf,
            const scalar phivL,
            const scalar phivR,
            scalar& aP,
            scalar& aN,
            scalar& aS,
            scalar& qP,
            scalar& qN,
            scalar& pF,
            scalar& pEF,
            scalar& pW,
            vector& Uf,
            scalar& amax
        )
        {
            const scalar alpha = 3.0/16.0;
            const scalar beta = 0.125;
            const scalar area = mag(Sf);
            const scalar uL = phivL/area;
            const scalar uR = phivR/area;
            const scalar c12 = 0.5*(cL + cR);

            // Guard only a degenerate acoustic state.  Positive-temperature
            // reacting-gas states are well away from this path.
            if (mag(c12) < VSMALL)
            {
                aP = 0.5;
                aN = 0.5;
                aS = 0.0;
                qP = 0.5*phivL;
                qN = 0.5*phivR;
                pF = 0.5*(pL + pR);
                pEF = qP*pL + qN*pR;
                pW = pF;
                Uf = 0.5*(UL + UR);
                amax = max(mag(uL) + cL, mag(uR) + cR)*area;
                return;
            }

            const scalar MaL = uL/c12;
            const scalar MaR = uR/c12;
            const scalar absMaL = mag(MaL);
            const scalar absMaR = mag(MaR);

            scalar Ma4L = max(MaL, scalar(0));
            scalar P5L = MaL >= 0 ? 1.0 : 0.0;

            if (absMaL < 1.0)
            {
                const scalar mm1 = sqr(MaL) - 1.0;
                Ma4L = 0.25*sqr(MaL + 1.0) + beta*sqr(mm1);
                P5L =
                    0.25*sqr(MaL + 1.0)*(2.0 - MaL)
                  + alpha*MaL*sqr(mm1);
            }

            scalar Ma4R = min(MaR, scalar(0));
            scalar P5R = MaR < 0 ? 1.0 : 0.0;

            if (absMaR < 1.0)
            {
                const scalar mm1 = sqr(MaR) - 1.0;
                Ma4R = -0.25*sqr(MaR - 1.0) - beta*sqr(mm1);
                P5R =
                    0.25*sqr(MaR - 1.0)*(2.0 + MaR)
                  - alpha*MaR*sqr(mm1);
            }

            const scalar Ma12 = Ma4L + Ma4R;
            const scalar P12 = P5L*pL + P5R*pR;
            const scalar phiVol = area*c12*Ma12;

            aS = 0.0;
            pF = P12;
            amax = max(mag(uL) + cL, mag(uR) + cR)*area;

            if (Ma12 >= 0.0)
            {
                aP = 1.0;
                aN = 0.0;
                qP = phiVol;
                qN = 0.0;
                pEF = phiVol*pL;
                pW = pL;
                Uf = UL;
            }
            else
            {
                aP = 0.0;
                aN = 1.0;
                qP = 0.0;
                qN = phiVol;
                pEF = phiVol*pR;
                pW = pR;
                Uf = UR;
            }
        };

        scalarField& aP = aPos_.ref().primitiveFieldRef();
        scalarField& aN = aNeg_.ref().primitiveFieldRef();
        scalarField& aS = aSf_.ref().primitiveFieldRef();
        scalarField& qP = aphivPos_.ref().primitiveFieldRef();
        scalarField& qN = aphivNeg_.ref().primitiveFieldRef();
        scalarField& pF = pFlux_.ref().primitiveFieldRef();
        scalarField& pEF = pEnergyFlux_.ref().primitiveFieldRef();
        scalarField& pW = pWork_.ref().primitiveFieldRef();
        vectorField& Uf = Uf_.ref().primitiveFieldRef();
        scalarField& amax = amaxSf_.ref().primitiveFieldRef();

        const scalarField& rhoP = rhoPos_().primitiveField();
        const scalarField& rhoN = rhoNeg_().primitiveField();
        const vectorField& UP = UPos_().primitiveField();
        const vectorField& UN = UNeg_().primitiveField();
        const scalarField& pP = pPos_().primitiveField();
        const scalarField& pN = pNeg_().primitiveField();
        const scalarField& cP = cPos.primitiveField();
        const scalarField& cN = cNeg.primitiveField();
        const vectorField& Sf = mesh.Sf().primitiveField();
        const scalarField& phP = phivPos.primitiveField();
        const scalarField& phN = phivNeg.primitiveField();

        forAll(aP, facei)
        {
            setAUSMPlusFace
            (
                rhoP[facei], rhoN[facei],
                UP[facei], UN[facei],
                pP[facei], pN[facei],
                cP[facei], cN[facei],
                Sf[facei], phP[facei], phN[facei],
                aP[facei], aN[facei], aS[facei],
                qP[facei], qN[facei], pF[facei], pEF[facei], pW[facei],
                Uf[facei], amax[facei]
            );
        }

        forAll(mesh.boundary(), patchi)
        {
            fvsPatchScalarField& aPp = aPos_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& aNp = aNeg_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& aSp = aSf_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& qPp = aphivPos_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& qNp = aphivNeg_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& pFp = pFlux_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& pEFp = pEnergyFlux_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& pWp = pWork_.ref().boundaryFieldRef()[patchi];
            fvsPatchVectorField& Ufp = Uf_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& amaxp = amaxSf_.ref().boundaryFieldRef()[patchi];

            const fvsPatchScalarField& rhoPp = rhoPos_().boundaryField()[patchi];
            const fvsPatchScalarField& rhoNp = rhoNeg_().boundaryField()[patchi];
            const fvsPatchVectorField& UPp = UPos_().boundaryField()[patchi];
            const fvsPatchVectorField& UNp = UNeg_().boundaryField()[patchi];
            const fvsPatchScalarField& pPp = pPos_().boundaryField()[patchi];
            const fvsPatchScalarField& pNp = pNeg_().boundaryField()[patchi];
            const fvsPatchScalarField& cPp = cPos.boundaryField()[patchi];
            const fvsPatchScalarField& cNp = cNeg.boundaryField()[patchi];
            const fvsPatchVectorField& Sfp = mesh.Sf().boundaryField()[patchi];
            const fvsPatchScalarField& phPp = phivPos.boundaryField()[patchi];
            const fvsPatchScalarField& phNp = phivNeg.boundaryField()[patchi];

            forAll(aPp, facei)
            {
                setAUSMPlusFace
                (
                    rhoPp[facei], rhoNp[facei],
                    UPp[facei], UNp[facei],
                    pPp[facei], pNp[facei],
                    cPp[facei], cNp[facei],
                    Sfp[facei], phPp[facei], phNp[facei],
                    aPp[facei], aNp[facei], aSp[facei],
                    qPp[facei], qNp[facei], pFp[facei], pEFp[facei],
                    pWp[facei], Ufp[facei], amaxp[facei]
                );
            }
        }

        phi_ = aphivPos_()*rhoPos_() + aphivNeg_()*rhoNeg_();
        return;
    }

    if (fluxScheme_ == "AUSM+up")
    {
        // Original detonationFoam AUSM+up (Liou 2006).  In addition to the
        // AUSM+ polynomial splitting, AUSM+up adds a pressure-difference
        // correction to the interface Mach number and a velocity-difference
        // correction to the split pressure.  Keep the corrected volumetric
        // mass flux as the single continuity/species convection flux.
        aPos_ = surfaceScalarField::New
        (
            "a_pos",
            mesh,
            dimensionedScalar(dimless, 0)
        );
        aNeg_ = surfaceScalarField::New
        (
            "a_neg",
            mesh,
            dimensionedScalar(dimless, 0)
        );
        aSf_ = surfaceScalarField::New("aSf", 0.0*phivPos);
        aphivPos_ = surfaceScalarField::New("aphiv_pos", 0.0*phivPos);
        aphivNeg_ = surfaceScalarField::New("aphiv_neg", 0.0*phivNeg);
        pFlux_ = surfaceScalarField::New
        (
            "pFlux",
            mesh,
            dimensionedScalar(dimPressure, 0)
        );
        pEnergyFlux_ = surfaceScalarField::New
        (
            "pEnergyFlux",
            0.0*pPos_()*phivPos
        );
        pWork_ = surfaceScalarField::New
        (
            "pWork",
            mesh,
            dimensionedScalar(dimPressure, 0)
        );
        Uf_ = surfaceVectorField::New
        (
            "Uf",
            mesh,
            dimensionedVector(dimVelocity, Zero)
        );
        amaxSf_ = surfaceScalarField::New("amaxSf", 0.0*mag(phivPos));

        auto setAUSMPlusUpFace = []
        (
            const scalar rhoL,
            const scalar rhoR,
            const vector& UL,
            const vector& UR,
            const scalar pL,
            const scalar pR,
            const scalar cL,
            const scalar cR,
            const vector& Sf,
            const scalar phivL,
            const scalar phivR,
            scalar& aP,
            scalar& aN,
            scalar& aS,
            scalar& qP,
            scalar& qN,
            scalar& pF,
            scalar& pEF,
            scalar& pW,
            vector& Uf,
            scalar& amax
        )
        {
            const scalar beta = 0.125;
            const scalar Kp = 0.25;
            const scalar Ku = 0.75;
            const scalar sigma = 1.0;
            const scalar fa = 1.0;
            const scalar alpha = 3.0/16.0*(5.0*sqr(fa) - 4.0);

            const scalar area = mag(Sf);
            const scalar uL = phivL/area;
            const scalar uR = phivR/area;
            const scalar c12 = sqrt(0.5*(sqr(cL) + sqr(cR)));

            if (mag(c12) < VSMALL || rhoL + rhoR < VSMALL)
            {
                aP = 0.5;
                aN = 0.5;
                aS = 0.0;
                qP = 0.5*phivL;
                qN = 0.5*phivR;
                pF = 0.5*(pL + pR);
                pEF = qP*pL + qN*pR;
                pW = pF;
                Uf = 0.5*(UL + UR);
                amax = max(mag(uL) + cL, mag(uR) + cR)*area;
                return;
            }

            const scalar MaL = uL/c12;
            const scalar MaR = uR/c12;
            const scalar absMaL = mag(MaL);
            const scalar absMaR = mag(MaR);

            scalar Ma4L = max(MaL, scalar(0));
            scalar P5L = MaL >= 0 ? 1.0 : 0.0;
            if (absMaL < 1.0)
            {
                const scalar mm1 = sqr(MaL) - 1.0;
                Ma4L = 0.25*sqr(MaL + 1.0) + beta*sqr(mm1);
                P5L =
                    0.25*sqr(MaL + 1.0)*(2.0 - MaL)
                  + alpha*MaL*sqr(mm1);
            }

            scalar Ma4R = min(MaR, scalar(0));
            scalar P5R = MaR < 0 ? 1.0 : 0.0;
            if (absMaR < 1.0)
            {
                const scalar mm1 = sqr(MaR) - 1.0;
                Ma4R = -0.25*sqr(MaR - 1.0) - beta*sqr(mm1);
                P5R =
                    0.25*sqr(MaR - 1.0)*(2.0 + MaR)
                  - alpha*MaR*sqr(mm1);
            }

            const scalar MaBarSqr =
                (sqr(uL) + sqr(uR))/(2.0*sqr(c12));
            const scalar Ma12 =
                Ma4L + Ma4R
              - 2.0*Kp/fa*max(1.0 - sigma*MaBarSqr, scalar(0))
               *(pR - pL)/((rhoL + rhoR)*sqr(c12));

            const scalar P12 =
                P5L*pL + P5R*pR
              - Ku*fa*c12*P5L*P5R*(rhoL + rhoR)*(uR - uL);

            const scalar phiVol = area*c12*Ma12;

            aS = 0.0;
            pF = P12;
            amax = max(mag(uL) + cL, mag(uR) + cR)*area;

            if (Ma12 >= 0.0)
            {
                aP = 1.0;
                aN = 0.0;
                qP = phiVol;
                qN = 0.0;
                pEF = phiVol*pL;
                pW = pL;
                Uf = UL;
            }
            else
            {
                aP = 0.0;
                aN = 1.0;
                qP = 0.0;
                qN = phiVol;
                pEF = phiVol*pR;
                pW = pR;
                Uf = UR;
            }
        };

        scalarField& aP = aPos_.ref().primitiveFieldRef();
        scalarField& aN = aNeg_.ref().primitiveFieldRef();
        scalarField& aS = aSf_.ref().primitiveFieldRef();
        scalarField& qP = aphivPos_.ref().primitiveFieldRef();
        scalarField& qN = aphivNeg_.ref().primitiveFieldRef();
        scalarField& pF = pFlux_.ref().primitiveFieldRef();
        scalarField& pEF = pEnergyFlux_.ref().primitiveFieldRef();
        scalarField& pW = pWork_.ref().primitiveFieldRef();
        vectorField& Uf = Uf_.ref().primitiveFieldRef();
        scalarField& amax = amaxSf_.ref().primitiveFieldRef();

        const scalarField& rhoP = rhoPos_().primitiveField();
        const scalarField& rhoN = rhoNeg_().primitiveField();
        const vectorField& UP = UPos_().primitiveField();
        const vectorField& UN = UNeg_().primitiveField();
        const scalarField& pP = pPos_().primitiveField();
        const scalarField& pN = pNeg_().primitiveField();
        const scalarField& cP = cPos.primitiveField();
        const scalarField& cN = cNeg.primitiveField();
        const vectorField& Sf = mesh.Sf().primitiveField();
        const scalarField& phP = phivPos.primitiveField();
        const scalarField& phN = phivNeg.primitiveField();

        forAll(aP, facei)
        {
            setAUSMPlusUpFace
            (
                rhoP[facei], rhoN[facei],
                UP[facei], UN[facei],
                pP[facei], pN[facei],
                cP[facei], cN[facei],
                Sf[facei], phP[facei], phN[facei],
                aP[facei], aN[facei], aS[facei],
                qP[facei], qN[facei], pF[facei], pEF[facei], pW[facei],
                Uf[facei], amax[facei]
            );
        }

        forAll(mesh.boundary(), patchi)
        {
            fvsPatchScalarField& aPp = aPos_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& aNp = aNeg_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& aSp = aSf_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& qPp = aphivPos_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& qNp = aphivNeg_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& pFp = pFlux_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& pEFp = pEnergyFlux_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& pWp = pWork_.ref().boundaryFieldRef()[patchi];
            fvsPatchVectorField& Ufp = Uf_.ref().boundaryFieldRef()[patchi];
            fvsPatchScalarField& amaxp = amaxSf_.ref().boundaryFieldRef()[patchi];

            const fvsPatchScalarField& rhoPp = rhoPos_().boundaryField()[patchi];
            const fvsPatchScalarField& rhoNp = rhoNeg_().boundaryField()[patchi];
            const fvsPatchVectorField& UPp = UPos_().boundaryField()[patchi];
            const fvsPatchVectorField& UNp = UNeg_().boundaryField()[patchi];
            const fvsPatchScalarField& pPp = pPos_().boundaryField()[patchi];
            const fvsPatchScalarField& pNp = pNeg_().boundaryField()[patchi];
            const fvsPatchScalarField& cPp = cPos.boundaryField()[patchi];
            const fvsPatchScalarField& cNp = cNeg.boundaryField()[patchi];
            const fvsPatchVectorField& Sfp = mesh.Sf().boundaryField()[patchi];
            const fvsPatchScalarField& phPp = phivPos.boundaryField()[patchi];
            const fvsPatchScalarField& phNp = phivNeg.boundaryField()[patchi];

            forAll(aPp, facei)
            {
                setAUSMPlusUpFace
                (
                    rhoPp[facei], rhoNp[facei],
                    UPp[facei], UNp[facei],
                    pPp[facei], pNp[facei],
                    cPp[facei], cNp[facei],
                    Sfp[facei], phPp[facei], phNp[facei],
                    aPp[facei], aNp[facei], aSp[facei],
                    qPp[facei], qNp[facei], pFp[facei], pEFp[facei],
                    pWp[facei], Ufp[facei], amaxp[facei]
                );
            }
        }

        phi_ = aphivPos_()*rhoPos_() + aphivNeg_()*rhoNeg_();
        return;
    }

    aPos_ = surfaceScalarField::New
    (
        "a_pos",
        fluxScheme_ == "Tadmor"
      ? surfaceScalarField::New("a_pos", mesh, 0.5)
      : ap/(ap - am)
    );
    aNeg_ = surfaceScalarField::New("a_neg", 1.0 - aPos_());

    phivPos *= aPos_();
    phivNeg *= aNeg_();

    aSf_ = surfaceScalarField::New
    (
        "aSf",
        fluxScheme_ == "Tadmor"
      ? -0.5*max(mag(am), mag(ap))
      : am*aPos_()
    );

    aphivPos_ = surfaceScalarField::New("aphiv_pos", phivPos - aSf_());
    aphivNeg_ = surfaceScalarField::New("aphiv_neg", phivNeg + aSf_());

    pFlux_ = surfaceScalarField::New
    (
        "pFlux",
        aPos_()*pPos_() + aNeg_()*pNeg_()
    );
    pEnergyFlux_ = surfaceScalarField::New
    (
        "pEnergyFlux",
        aphivPos_()*pPos_()
      + aphivNeg_()*pNeg_()
      + aSf_()*(pPos_() - pNeg_())
    );
    pWork_ = surfaceScalarField::New
    (
        "pWork",
        aPos_()*pPos_() + aNeg_()*pNeg_()
    );
    Uf_ = surfaceVectorField::New
    (
        "Uf",
        aPos_()*UPos_() + aNeg_()*UNeg_()
    );
    amaxSf_ = surfaceScalarField::New
    (
        "amaxSf",
        max(mag(aphivPos_()), mag(aphivNeg_()))
    );

    phi_ = aphivPos_()*rhoPos_() + aphivNeg_()*rhoNeg_();
}

// ************************************************************************* //
