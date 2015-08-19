//
// Created by lichon on 18.08.15.
//

#ifndef SRC_AAA_H
#define SRC_AAA_H

#include "G4LogicalVolumeStore.hh"

class LogicalVolumeStoreFix {
public:
    static void copyToSingleton(G4LogicalVolumeStore* theStore);
};


#endif //SRC_AAA_H
