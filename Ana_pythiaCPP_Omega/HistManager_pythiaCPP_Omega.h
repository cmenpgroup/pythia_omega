#ifndef HISTMANAGER_H
#define HISTMANAGER_H
#include <vector>
#include <string>
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "PartComb_pythiaCPP_Omega.h"
#include "Decay1_pythiaCPP_Omega.h"
#include "Decay2_pythiaCPP_Omega.h"
#include "Decay3_pythiaCPP_Omega.h"

class HistManager_pythiaCPP_Omega
{
private:
    int MAX_SECTORS; // max. number of CLAS sectors
    double LIGHTSPEED; // speed of light in cm/ns
    double BEAM_ENERGY; // electron beam energy in GeV
    
    TH1D *hIntermedPart;
    TH2D *hPartPerEvt;
    TH2D *hIMomega;
    TH2D *hIMomega_CutDalitz;
    TH2D *hMMsq_X;
    TH2D *hMMsq_Xrecoil;
    TH2D *hMMsq_Xpip;
    TH2D *hMMsq_Xpim;
    TH2D *hMMsq_Xpi0;
    TH2D *hOpAng_PipPim;
    TH2D *hOpAng_PipPimPi0;
    TH1D *hQsq_NoCuts;
    TH1D *hNu_NoCuts;
    TH1D *hW_NoCuts;
    TH1D *hQsq;
    TH1D *hNu;
    TH1D *hW;
    TH2D *hz_fracEnergy;
    TH2D *hOpAng_TwoPhoton;
    TH2D *hIMomega_VS_IMPipPim[20];  // one histogram for each particle combination
    TH2D *hDalitz_pip[20];  // one Dalitz histogram for each particle combination
    TH2D *hDalitz_pip_CutDalitz[20];  // one Dalitz histogram for each particle combination with Dalitz cut
    TH2D *hIM_NotOmega_Decay1[4]; // pi+ pi- pi0 inv. mass for events that are not omega, 1 part. decay
    TH2D *hIM_NotOmega_Decay2[4]; // pi+ pi- pi0 inv. mass for events that are not omega, 2 part. decay
    TH2D *hIM_NotOmega_Decay3[4]; // pi+ pi- pi0 inv. mass for events that are not omega, 3 part. decay
    TH2D *hIM_Unknown_Decay1[4]; // pi+ pi- pi0 inv. mass for events with unknown parent, 1 part. decay
    TH2D *hIM_Unknown_Decay2[4]; // pi+ pi- pi0 inv. mass for events with unknown parent, 2 part. decay
    TH2D *hIM_Unknown_Decay3[4]; // pi+ pi- pi0 inv. mass for events with unknown parent, 3 part. decay
public:
    HistManager_pythiaCPP_Omega();
    void BookHist();
    void WriteHist(std::string outFile);
    TH1D* Get_hIntermedPart() { return hIntermedPart; };
    TH1D* Get_hQsq_NoCuts() { return hQsq_NoCuts; };
    TH1D* Get_hNu_NoCuts() { return hNu_NoCuts; };
    TH1D* Get_hW_NoCuts() { return hW_NoCuts; };
    TH1D* Get_hQsq() { return hQsq; };
    TH1D* Get_hNu() { return hNu; };
    TH1D* Get_hW() { return hW; };
    
    TH2D* Get_hPartPerEvt() { return hPartPerEvt; };
    TH2D* Get_hIMomega() { return hIMomega; };
    TH2D* Get_hIMomega_CutDalitz() { return hIMomega_CutDalitz; };
    TH2D* Get_hMMsq_X() { return hMMsq_X; };
    TH2D* Get_hMMsq_Xrecoil() { return hMMsq_Xrecoil; };
    TH2D* Get_hMMsq_Xpip() { return hMMsq_Xpip; };
    TH2D* Get_hMMsq_Xpim() { return hMMsq_Xpim; };
    TH2D* Get_hMMsq_Xpi0() { return hMMsq_Xpi0; };
    TH2D* Get_hOpAng_PipPim() { return hOpAng_PipPim; };
    TH2D* Get_hOpAng_PipPimPi0() { return hOpAng_PipPimPi0; };
    TH2D* Get_hz_fracEnergy() { return hz_fracEnergy; };
    TH2D* Get_hOpAng_TwoPhoton() { return hOpAng_TwoPhoton; };
    
    TH2D* Get_hIMomega_VS_IMPipPim(int index) { return hIMomega_VS_IMPipPim[index]; };
    TH2D* Get_hDalitz_pip(int index) { return hDalitz_pip[index]; };
    TH2D* Get_hDalitz_pip_CutDalitz(int index) { return hDalitz_pip_CutDalitz[index]; };

    TH2D* Get_hIM_NotOmega_Decay1(int index) { return hIM_NotOmega_Decay1[index]; };
    TH2D* Get_hIM_NotOmega_Decay2(int index) { return hIM_NotOmega_Decay2[index]; };
    TH2D* Get_hIM_NotOmega_Decay3(int index) { return hIM_NotOmega_Decay3[index]; };

    TH2D* Get_hIM_Unknown_Decay1(int index) { return hIM_Unknown_Decay1[index]; };
    TH2D* Get_hIM_Unknown_Decay2(int index) { return hIM_Unknown_Decay2[index]; };
    TH2D* Get_hIM_Unknown_Decay3(int index) { return hIM_Unknown_Decay3[index]; };
    
    int GetMAX_SECTORS() { return MAX_SECTORS; };
    double GetLIGHTSPEED() { return LIGHTSPEED; };
    double GetBEAM_ENERGY() { return BEAM_ENERGY; };
    
};
#endif