#ifndef planarRefiner_fvMeshTopoChanger_H
#define planarRefiner_fvMeshTopoChanger_H

#include "fvMeshTopoChanger.H"
#include "IOdictionary.H"
#include "word.H"

namespace Foam
{
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
    label maxRefinementIterations_;
    label nRefinementIterations_;
    word geometry_;
    label timeIndex_;
    mutable IOdictionary state_;

    word detectedGeometry() const;
    void validateGeometry();
    dictionary directionalRefineDict() const;
    labelList selectCells() const;

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
