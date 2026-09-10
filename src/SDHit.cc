#include "SDHit.hh"

G4ThreadLocal G4Allocator<SDHit>* SDHitAllocator = nullptr;

SDHit::SDHit(G4int volID)
    : volumeID(volID), edep(0.0), tmin(DBL_MAX) {}
    
void SDHit::AddEdep(G4double val) {
    edep += val;
}

G4double SDHit::GetEdep() const {
    return edep;
}

void SDHit::UpdateTmin(G4double t) {
    if (t < tmin) {
        tmin = t;
    }
}

