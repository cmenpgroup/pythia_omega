#include <iostream>
#include <fstream>
#include <string>
#include "HistManager_pythiaCPP_Omega.h"
#include "PartComb_pythiaCPP_Omega.h"
#include "PDG_pythiaCPP_Omega.h"
#include "Cuts_pythiaCPP_Omega.h"
#include "Decay1_pythiaCPP_Omega.h"
#include "Decay2_pythiaCPP_Omega.h"
#include "Decay3_pythiaCPP_Omega.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TVector3.h"
#include "TFile.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TH1D.h"
#include "TH2D.h"

using namespace std;

#define MAX_ANA 2
#define MAX_TRACKS 50

HistManager_pythiaCPP_Omega myHistManager; // declare the histogram manager

int process (string inFile, int iAna, int MaxEvents, int dEvents);
bool Cut_CLAS6_Theta_EC(double theta);
bool Cut_CLAS6_Theta_TOF(double theta);
bool Cut_CLAS6_Theta_Ana(TLorentzVector V4, int iAna, string detName);
void FillHists_omega(TLorentzVector pip, TLorentzVector pim, TLorentzVector pi0, double nu, int iPC);