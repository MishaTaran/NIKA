#ifndef SIZES_HH
#define SIZES_HH

#include <G4SystemOfUnits.hh>
#include <algorithm>

namespace Sizes
{
    const G4bool viewMode = false;

    namespace Instrument {
        const G4double sizeX = 81.0 * mm;
        const G4double sizeY = 83.8 * mm;
        const G4double sizeZ = 83.0 * mm;

        const G4double halfX = sizeX / 2.0;
        const G4double halfY = sizeY / 2.0;
        const G4double halfZ = sizeZ / 2.0;

        const G4double bottomZ = 0.0 * mm;

        inline G4double centerZ() { return sizeZ / 2.0; }
    }


    namespace Trigger {

        const G4int nLayers = 3;
        const G4double panelSizeX = 45.0 * mm;
        const G4double panelSizeY = 45.0 * mm;
        const G4double halfX = panelSizeX / 2.0;
        const G4double halfY = panelSizeY / 2.0;
        const G4double polkaSmallXY = 9.4 * mm;
        const G4double polkaSmallZ = 1.5 * mm;
        const G4double standXY = 3.0 * mm;
        const G4double standrubXY = 2.3 * mm;
        const G4double standZ = 13.0 * mm;
        const G4double spmXY = 5.0 * mm;
        const G4double spmZ = 0.6 * mm;
        const G4double SiPMFrame = 0.5 * mm;
        const G4double SiPMWindowThick = 0.21 * mm;
        const G4double plataTOFXY = 1.5 * mm;
        const G4double standAlXY = 1.0 * mm;
        const G4double standAlZ = 8.0 * mm;
        const G4double thickness = 5.0 * mm;
        const G4double tyvek = 0.5 * mm;
        const G4double thicknessTrig = 12.0 * mm;
        const G4double alSizeX = 66.0 * mm;
        const G4double alSizeY = 66.0 * mm;
    }


    namespace VetoAC {
        //Боковая часть
        const G4double VetoHeight = 73.0 * mm;
        const G4double thickness = 5.0 * mm;
        const G4double alTSide = 1.0 * mm;
        const G4double alTSideDown = 3.75 * mm;
        const G4double PlataSide = 2.9 * mm;
        const G4double Plata = 3.9 * mm;
        const G4double alTSideUp = 1.5 * mm;
        const G4double sideVetoXY = 33.4 * mm;
        const G4double polkaBigY = 2.0 * mm;
        const G4double polkaBigX = 66.0 * mm;
        const G4double polkaBigZ = 75.75 * mm;

        //Нижняя часть
        const G4double downVetoXY = 30.9 * mm;
        const G4double VetoDownXY = 30.9 * mm;
        const G4double alTDownXY = 15.0 * mm;
        const G4double alTDownZ = 3.8 * mm;
        const G4double alTDown = 1.0 * mm;
        const G4double polkaXY = 1.1 * mm;
        const G4double polkaZ = 13.0 * mm;
        const G4double DownVeloAll = alTDownZ + alTDown * 2.0 + Plata + thickness + 2.0 * Trigger::tyvek - 0.05 * mm;
    }

    namespace Envelope {
        const G4double halfX = std::max(Instrument::halfX, Instrument::sizeX / 2.0);
        const G4double halfY = std::max(Instrument::halfY, Instrument::sizeY / 2.0);
    }

    namespace CubeInner {
        const G4double halfX = Envelope::halfX;
        const G4double halfY = Envelope::halfY;
    }
}
#endif