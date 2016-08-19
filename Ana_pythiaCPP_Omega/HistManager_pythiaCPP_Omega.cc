#include <vector>
#include <string>
#include "HistManager_pythiaCPP_Omega.h"
#include <iostream>

HistManager_pythiaCPP_Omega::HistManager_pythiaCPP_Omega()
{
    MAX_SECTORS = 6;
    LIGHTSPEED = 30.0;
    BEAM_ENERGY = 5.01; // 4.5
}

//
// BookHist - routine to set up histograms
//
void HistManager_pythiaCPP_Omega::BookHist(){

    int j, k;

    char hname[100];
    char htitle[100];
    
    int nTheta = 180;
    double ThetaLo = 0.0;
    double ThetaHi = 180.0;

    int nIM2pi = 200;
    double IM2piLo = 0.0;
    double IM2piHi = 1.0;

    int nIMomega = 125;
    double IMomegaLo = 0.0;
    double IMomegaHi = 2.5;

    int nMMsq = 200;
    double MMsqLo = 0.0;
    double MMsqHi = 5.0;
    
    int nDalitzX = 100;
    double DalitzX_Lo = 0.0;
    double DalitzX_Hi = 1.0;
    
    int nDalitzY = 100;
    double DalitzY_Lo = 0.0;
    double DalitzY_Hi = 1.0;

    PartComb_pythiaCPP_Omega myComb;
    int nComb = myComb.Get_nPartComb();
    int maxCtr = myComb.Get_nCtr();
    int maxCombPhoton = myComb.Get_nCombPhoton();
    
    Decay1_pythiaCPP_Omega myDecay1; // declare the object to handle the one body decay
    int nDecay1 = myDecay1.Get_nDecays();

    Decay2_pythiaCPP_Omega myDecay2; // declare the object to handle the one body decay
    int nDecay2 = myDecay2.Get_nDecays();

    Decay3_pythiaCPP_Omega myDecay3; // declare the object to handle the one body decay
    int nDecay3 = myDecay3.Get_nDecays();
    
    hIntermedPart = new TH1D("hIntermedPart","Intermediate Particle ID", 3500, 0.5, 3500.5);
    hPartPerEvt = new TH2D("hPartPerEvt", "Part. Comb. per Event",nComb+1,-0.5,nComb+0.5,maxCtr,0.5,maxCtr+0.5);

    hIMomega = new TH2D("hIMomega", "Mass(#omega)",nIMomega,IMomegaLo,IMomegaHi,nComb,0.5,nComb+0.5);
    hIMomega_CutDalitz = new TH2D("hIMomega_CutDalitz", "Mass(#omega) w/ Dalitz CUt",nIMomega,IMomegaLo,IMomegaHi,nComb,0.5,nComb+0.5);
    
    hMMsq_X = new TH2D("hMMsq_X", "MM^{2}, total",nMMsq,MMsqLo,MMsqHi+5.0,nComb,0.5,nComb+0.5);
    hMMsq_Xrecoil = new TH2D("hMMsq_Xrecoil", "MM^{2} off proton",nMMsq,MMsqLo,MMsqHi,nComb,0.5,nComb+0.5);
    hMMsq_Xpip = new TH2D("hMMsq_Xpip", "MM^{2} off pi^{+}",nMMsq,MMsqLo,MMsqHi+2.0,nComb,0.5,nComb+0.5);
    hMMsq_Xpim = new TH2D("hMMsq_Xpim", "MM^{2} off pi^{-}",nMMsq,MMsqLo,MMsqHi+2.0,nComb,0.5,nComb+0.5);
    hMMsq_Xpi0 = new TH2D("hMMsq_Xpi0", "MM^{2} off pi^{0}",nMMsq,MMsqLo,MMsqHi+2.0,nComb,0.5,nComb+0.5);
    
    hOpAng_PipPim = new TH2D("hOpAng_PipPim","Opening Angle(#pi^{+},#pi^{-})",nTheta,ThetaLo,ThetaHi,nComb,0.5,nComb+0.5);

    hOpAng_PipPimPi0 = new TH2D("hOpAng_PipPimPi0","Opening Angle(#pi^{+}#pi^{-},#pi^{0})",nTheta,ThetaLo,ThetaHi,nComb,0.5,nComb+0.5);

    hQsq_NoCuts = new TH1D("hQsq_NoCuts","Q^{2}", 100, 0., 4.);
    hNu_NoCuts = new TH1D("hNu_NoCuts","#nu", 100, 0., 5.);
    hW_NoCuts = new TH1D("hW_NoCuts", "W", 250, 0., 5.);

    hQsq = new TH1D("hQsq","Q^{2}", 100, 0., 4.);
    hNu = new TH1D("hNu","#nu", 100, 0., 5.);
    hW = new TH1D("hW", "W", 250, 0., 5.);

    hz_fracEnergy = new TH2D("hz_fracEnergy", "Fractional Energy", 150, 0, 1.5,nComb,0.5,nComb+0.5);
    hOpAng_TwoPhoton = new TH2D("hOpAng_TwoPhoton","Opening Angle(#gamma,#gamma)",nTheta, ThetaLo, ThetaHi,maxCombPhoton,0.5,maxCombPhoton+0.5);
    
    for(j=0; j<nComb; j++){
        (nComb<10) ? sprintf(hname,"hIMomega_VS_IMPipPim_0%i",j) : sprintf(hname,"hIMomega_VS_IMPipPim_%i",j);
        sprintf(htitle,"Mass(#omega) vs. Mass(#pi^{+} #pi^{-}), %s",myComb.Get_PartComb(j).c_str());
        hIMomega_VS_IMPipPim[j] = new TH2D(hname,htitle,nIM2pi,IM2piLo,IM2piHi,nIMomega,IMomegaLo,IMomegaHi);

        (nComb<10) ? sprintf(hname,"hDalitz_pip_0%i",j) : sprintf(hname,"hDalitz_pip_%i",j);
        sprintf(htitle,"Dalitz for #omega study, case %i",j);
        hDalitz_pip[j] = new TH2D(hname, htitle, nDalitzX, DalitzX_Lo, DalitzX_Hi,nDalitzY, DalitzY_Lo, DalitzY_Hi);

        (nComb<10) ? sprintf(hname,"hDalitz_pip_CutDalitz_0%i",j) : sprintf(hname,"hDalitz_pip_CutDalitz_%i",j);
        sprintf(htitle,"Dalitz for #omega study, case %i",j);
        hDalitz_pip_CutDalitz[j] = new TH2D(hname, htitle, nDalitzX, DalitzX_Lo, DalitzX_Hi,nDalitzY, DalitzY_Lo, DalitzY_Hi);
    }
    
    for(k=0; k<4; k++){
        sprintf(hname,"hIM_NotOmega_Decay1_%i",k+1);
        sprintf(htitle,"One particle decay, Case %i",k+1);
        hIM_NotOmega_Decay1[k] = new TH2D(hname,htitle,nIMomega,IMomegaLo,IMomegaHi,nDecay1,0.5,nDecay1+0.5);

        sprintf(hname,"hIM_NotOmega_Decay2_%i",k+1);
        sprintf(htitle,"Two particle decay, Case %i",k+1);
        hIM_NotOmega_Decay2[k] = new TH2D(hname,htitle,nIMomega,IMomegaLo,IMomegaHi,nDecay2,0.5,nDecay2+0.5);
        
        sprintf(hname,"hIM_NotOmega_Decay3_%i",k+1);
        sprintf(htitle,"Three particle decay, Case %i",k+1);
        hIM_NotOmega_Decay3[k] = new TH2D(hname,htitle,nIMomega,IMomegaLo,IMomegaHi,nDecay3,0.5,nDecay3+0.5);
    }
}

//
// WriteHist - routine to write histograms to the output file
//
void HistManager_pythiaCPP_Omega::WriteHist(string outFile){

    PartComb_pythiaCPP_Omega myComb;
    int nComb = myComb.Get_nPartComb();
    
    int j, k;
    
    TFile *out = new TFile(outFile.c_str(), "recreate");
    out->cd();

    hIntermedPart->GetXaxis()->SetTitle("Intermediate Particle ID");
    hIntermedPart->GetYaxis()->SetTitle("Counts");
    hIntermedPart->Write();

    hPartPerEvt->GetXaxis()->SetTitle("Number of Particles per Event");
    hPartPerEvt->GetYaxis()->SetTitle("Particle Combination");
    hPartPerEvt->Write();

    hIMomega->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
    hIMomega->GetYaxis()->SetTitle("Particle Combination");
    hIMomega->Write();

    
    hIMomega_CutDalitz->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
    hIMomega_CutDalitz->GetYaxis()->SetTitle("Particle Combination");
    hIMomega_CutDalitz->Write();
    
    hMMsq_X->GetXaxis()->SetTitle("MM^{2} (GeV/c^{2})^{2}");
    hMMsq_X->GetYaxis()->SetTitle("Particle Combination");
    hMMsq_X->Write();

    hMMsq_Xrecoil->GetXaxis()->SetTitle("MM^{2} (GeV/c^{2})^{2}");
    hMMsq_Xrecoil->GetYaxis()->SetTitle("Particle Combination");
    hMMsq_Xrecoil->Write();
    
    hMMsq_Xpip->GetXaxis()->SetTitle("MM^{2} (GeV/c^{2})^{2}");
    hMMsq_Xpip->GetYaxis()->SetTitle("Particle Combination");
    hMMsq_Xpip->Write();
    
    hMMsq_Xpim->GetXaxis()->SetTitle("MM^{2} (GeV/c^{2})^{2}");
    hMMsq_Xpim->GetYaxis()->SetTitle("Particle Combination");
    hMMsq_Xpim->Write();
    
    hMMsq_Xpi0->GetXaxis()->SetTitle("MM^{2} (GeV/c^{2})^{2}");
    hMMsq_Xpi0->GetYaxis()->SetTitle("Particle Combination");
    hMMsq_Xpi0->Write();
    
    hQsq_NoCuts->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
    hQsq_NoCuts->GetYaxis()->SetTitle("Counts");
    hQsq_NoCuts->Write();

    hNu_NoCuts->GetXaxis()->SetTitle("\nu (GeV)");
    hNu_NoCuts->GetYaxis()->SetTitle("Counts");
    hNu_NoCuts->Write();

    hW_NoCuts->GetXaxis()->SetTitle("W (GeV)");
    hW_NoCuts->GetYaxis()->SetTitle("Counts");
    hW_NoCuts->Write();

    hQsq->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
    hQsq->GetYaxis()->SetTitle("Counts");
    hQsq->Write();
    
    hNu->GetXaxis()->SetTitle("\nu (GeV)");
    hNu->GetYaxis()->SetTitle("Counts");
    hNu->Write();
    
    hW->GetXaxis()->SetTitle("W (GeV)");
    hW->GetYaxis()->SetTitle("Counts");
    hW->Write();
    
    hz_fracEnergy->GetXaxis()->SetTitle("z");
    hz_fracEnergy->GetYaxis()->SetTitle("Counts");
    hz_fracEnergy->Write();

    hOpAng_PipPim->GetXaxis()->SetTitle("#theta_{#pi^{+},#pi^{-}}(deg.)");
    hOpAng_PipPim->GetYaxis()->SetTitle("Particle Combination");
    hOpAng_PipPim->Write();

    hOpAng_PipPimPi0->GetXaxis()->SetTitle("#theta_{#pi^{+}#pi^{-},#pi^{0}}(deg.)");
    hOpAng_PipPimPi0->GetYaxis()->SetTitle("Particle Combination");
    hOpAng_PipPimPi0->Write();

    hOpAng_TwoPhoton->GetXaxis()->SetTitle("#theta_{#gamma,#gamma}(deg.)");
    hOpAng_TwoPhoton->GetYaxis()->SetTitle("Particle Combination with 2 photons");
    hOpAng_TwoPhoton->Write();
    
    for(j=0; j<nComb; j++){
        hIMomega_VS_IMPipPim[j]->GetXaxis()->SetTitle("#pi^{+} #pi^{-} Mass (GeV/c^{2})");
        hIMomega_VS_IMPipPim[j]->GetYaxis()->SetTitle("#pi^{+} #pi^{-} #pi^{0} Mass (GeV/c^{2})");
        hIMomega_VS_IMPipPim[j]->Write();
    }
    
    for(j=0; j<nComb; j++){
        hDalitz_pip[j]->GetXaxis()->SetTitle("M^{2} (#pi^{0} #pi^{+}) (GeV/c^{2})^{2}");
        hDalitz_pip[j]->GetYaxis()->SetTitle("M^{2} (#pi^{-} #pi^{+}) (GeV/c^{2})^{2}");
        hDalitz_pip[j]->Write();
    }
    
    for(j=0; j<nComb; j++){
        hDalitz_pip_CutDalitz[j]->GetXaxis()->SetTitle("M^{2} (#pi^{0} #pi^{+}) (GeV/c^{2})^{2}");
        hDalitz_pip_CutDalitz[j]->GetYaxis()->SetTitle("M^{2} (#pi^{-} #pi^{+}) (GeV/c^{2})^{2}");
        hDalitz_pip_CutDalitz[j]->Write();
    }
    
    for(k=0; k<4; k++){
        hIM_NotOmega_Decay1[k]->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
        hIM_NotOmega_Decay1[k]->GetYaxis()->SetTitle("Particle Combination");
        hIM_NotOmega_Decay1[k]->Write();
    
        hIM_NotOmega_Decay2[k]->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
        hIM_NotOmega_Decay2[k]->GetYaxis()->SetTitle("Particle Combination");
        hIM_NotOmega_Decay2[k]->Write();
        
        hIM_NotOmega_Decay3[k]->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
        hIM_NotOmega_Decay3[k]->GetYaxis()->SetTitle("Particle Combination");
        hIM_NotOmega_Decay3[k]->Write();
    }
    
    out->Close();
}
