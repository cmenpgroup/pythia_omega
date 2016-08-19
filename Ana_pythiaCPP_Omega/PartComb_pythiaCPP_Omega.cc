#include "PartComb_pythiaCPP_Omega.h"

PartComb_pythiaCPP_Omega::PartComb_pythiaCPP_Omega()
{
    nCtr = 6; // combination counter
    // 1 : Number of pi+
    // 2 : Number of pi-
    // 3 : Number of pi0
    // 4 : Number of 2 photons
    // 5 : Number of omega mesons
    // 6 : Number of particle combinations listed below
    
    nCombPhoton = 8; // number of particle combination with 2 photons
    
    LabelPartComb.push_back("FS 1"); // Final State 1 = nPip>=1 && nPim>=1 && nPi0>=1
    LabelPartComb.push_back("FS 2"); // Final State 2 = nPip==1 && nPim==1 && nPi0==1
    LabelPartComb.push_back("FS 3"); // Final State 3 = nPip>=1 && nPim>=1 && nPhoton>=2
    LabelPartComb.push_back("FS 4"); // Final State 4 = nPip==1 && nPim==1 && nPhoton==2
    LabelPartComb.push_back("FS 1 && nOmega>=1");
    LabelPartComb.push_back("FS 2 && nOmega==1");
    LabelPartComb.push_back("FS 3 && nOmega>=1");
    LabelPartComb.push_back("FS 4 && nOmega==1");
    LabelPartComb.push_back("FS 1 && nOmega==0");
    LabelPartComb.push_back("FS 2 && nOmega==0");
    LabelPartComb.push_back("FS 3 && nOmega==0");
    LabelPartComb.push_back("FS 4 && nOmega==0");
    LabelPartComb.push_back("FS 1 && nOmega==0 && nEta==0");
    LabelPartComb.push_back("FS 2 && nOmega==0 && nEta==0");
    LabelPartComb.push_back("FS 3 && nOmega==0 && nEta==0");
    LabelPartComb.push_back("FS 4 && nOmega==0 && nEta==0");
    LabelPartComb.push_back("FS 1 && no other");
    LabelPartComb.push_back("FS 2 && no other");
    LabelPartComb.push_back("FS 3 && no other");
    LabelPartComb.push_back("FS 4 && no other");
}


