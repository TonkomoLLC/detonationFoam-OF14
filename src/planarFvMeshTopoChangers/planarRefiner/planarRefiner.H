#ifndef planarRefiner_fvMeshTopoChanger_H
#define planarRefiner_fvMeshTopoChanger_H

#include "fvMeshTopoChanger.H"
#include "IOdictionary.H"
#include "word.H"
#include "volFields.H"
#include "autoPtr.H"

namespace Foam
{

class undoableMeshCutter;

namespace fvMeshTopoChangers
{

class planarRefiner : public fvMeshTopoChanger
{
    dictionary dict_;
    label refineInterval_;
    word fieldName_;
    scalar lowerRefineLevel_;
    scalar upperRefineLevel_;
    label maxCells_;

    // Legacy refinement-only global event cap. Retained unchanged when
    // automaticUnrefinement=false.
    label maxRefinementIterations_;
    label nRefinementIterations_;

    // Optional reversible AMR controls. One planar refinement comprises two
    // binary directional cuts; maxRefinementLevel is therefore represented
    // internally by a binary-cut depth of 2*maxRefinementLevel.
    Switch automaticUnrefinement_;
    label unrefineInterval_;
    scalar lowerUnrefineLevel_;
    scalar upperUnrefineLevel_;
    label maxUnrefinementPassesPerUpdate_;
    label maxRefinementLevel_;
    label nUnrefinementPasses_;

    word geometry_;
    label timeIndex_;
    mutable IOdictionary state_;

    // The OpenFOAM undo tree is intentionally runtime-only. Version 1.1.0 does
    // not attempt to serialize/reconstruct undoableMeshCutter history.
    autoPtr<undoableMeshCutter> reversibleCutter_;

    // Runtime binary-cut depth. This is a registered field so fvMesh maps it
    // through the two directional topology changes performed by one planar
    // refinement. It is deliberately not written because restart of active
    // reversible history is prohibited in this candidate.
    volScalarField refinementLevel_;

    word detectedGeometry() const;
    void validateGeometry();
    void validateAutomaticUnrefinement();
    dictionary directionalRefineDict() const;
    labelList selectCells() const;
    labelList selectSplitFacesForUnrefinement() const;
    bool unrefineOnePass();

public:
    TypeName("planarRefiner");

    planarRefiner(fvMesh& mesh, const dictionary& dict);
    planarRefiner(const planarRefiner&) = delete;
    virtual ~planarRefiner();

    virtual bool update();
    virtual void topoChange(const polyTopoChangeMap&);
    virtual void mapMesh(const polyMeshMap&);
    virtual void distribute(const polyDistributionMap&);
    virtual bool write(const bool write=true) const;

    void operator=(const planarRefiner&) = delete;
};

}
}

#endif
