#include "Decay3_pythiaCPP_Omega.h"

Decay3_pythiaCPP_Omega::Decay3_pythiaCPP_Omega(){
    this->Init();

    // fill the map with the histogram numbers
    int nhist = 1; // start the histogram list at 1
    for (std::map<string,bool>::iterator it=evtDecay.begin(); it!=evtDecay.end(); ++it){
        evtDecayNum[it->first] = nhist;
        nhist++;
    }
}

void Decay3_pythiaCPP_Omega::Init(){
    evtDecay["RhoMPipProton"]=false;
    evtDecay["RhoPPimProton"]=false;
    evtDecay["Rho0Pi0Proton"]=false;
    evtDecay["RhoPRho0Neutron"]=false;
    evtDecay["RhoMPipDeltaP"]=false;
    evtDecay["RhoPEtaNeutron"]=false;
    evtDecay["RhoPRhoMProton"]=false;
    evtDecay["Rho0PipDelta0"]=false;
    evtDecay["RhoMPi0DeltaPP"]=false;
    evtDecay["OmegaPi0Proton"]=false;
    evtDecay["OmegaPimDeltaPP"]=false;
    evtDecay["OmegaRhoMDeltaPP"]=false;
    evtDecay["OmegaPipNeutron"]=false;
    evtDecay["OmegaOmegaProton"]=false;
    evtDecay["OmegaKStar892PLambda"]=false;
    evtDecay["OmegaK0SigmaP"]=false;
    evtDecay["OmegaK0Sigma1385"]=false;
    evtDecay["OmegaPi0Eta"]=false;
    evtDecay["OmegaPi0RhoP"]=false;
    evtDecay["Rho0EtaProton"]=false;
    evtDecay["RhoMRhoPDeltaP"]=false;
    evtDecay["OmegaRhoPNeutron"]=false;
    evtDecay["OmegaPipDelta0"]=false;
    evtDecay["PimPi0DeltaPP"]=false;
    evtDecay["Pi0EtaProton"]=false;
    evtDecay["Pi0EtaDelta0"]=false;
    evtDecay["PipEtaNeutron"]=false;
    evtDecay["PimEtaDeltaPP"]=false;
    evtDecay["RhoKPSigma0"]=false;
    evtDecay["Pi0Rho0DeltaP"]=false;
    evtDecay["PimRhoPDeltaP"]=false;
    evtDecay["Pi0RhoPDelta0"]=false;
    evtDecay["PipPimDeltaP"]=false;
    evtDecay["KStar892PKStar892MProton"]=false;
    evtDecay["OmegaOmegaDeltaP"]=false;
    evtDecay["OmegaPi0DeltaP"]=false;
    evtDecay["OmegaPi0Pip"]=false;
    evtDecay["OmegaPimPip"]=false;
    evtDecay["OmegaEtaProton"]=false;
    evtDecay["OmegaEtaDeltaP"]=false;
    evtDecay["OmegaEtaPrimeProton"]=false;
    evtDecay["OmegaRho0Proton"]=false;
    evtDecay["OmegaRho0DeltaP"]=false;
    evtDecay["OmegaRhoPDelta0"]=false;
    evtDecay["OmegaKPLambda"]=false;
    evtDecay["OmegaKPSigma0"]=false;
    evtDecay["OmegaKPSigmaStar0"]=false;
    evtDecay["Pi0PipDelta0"]=false;
    evtDecay["Pi0PimRhoP"]=false;
    evtDecay["Pi0EtaDeltaP"]=false;
    evtDecay["Pi0EtaPrimeProton"]=false;
    evtDecay["Pi0EtaPrimeDeltaP"]=false;
    evtDecay["Pi0K0Sigma1385"]=false;
    evtDecay["Pi0KStar892PLambda"]=false;
    evtDecay["Pi0KStar892SigmaP"]=false;
    evtDecay["Pi0KStar892PSigma0"]=false;
    evtDecay["Pi0KStar892PSigmaStar0"]=false;
    evtDecay["Pi0PipPim"]=false;
    evtDecay["Pi0PipEta"]=false;
    evtDecay["Pi0PipRhoM"]=false;
    evtDecay["Pi0K0SigmaP"]=false;
    evtDecay["Pi0KPSigmaStar0"]=false;
    evtDecay["PimPipRhoP"]=false;
    evtDecay["PimKStar892PSigmaP"]=false;
    evtDecay["PimEtaPrimeDeltaPP"]=false;
    evtDecay["PimKPSigma1385"]=false;
    evtDecay["PimKStar892PSigma1385"]=false;
    evtDecay["PipEtaPrimeNeutron"]=false;
    evtDecay["PipEtaPrimeDelta0"]=false;
    evtDecay["PipK0Lambda"]=false;
    evtDecay["PipK0Sigma0"]=false;
    evtDecay["PipK0SigmaStar0"]=false;
    evtDecay["PipKStar892Lambda"]=false;
    evtDecay["PipKStar892Sigma0"]=false;
    evtDecay["PipKStar892SigmaStar0"]=false;
    evtDecay["EtaRho0DeltaP"]=false;
    evtDecay["EtaPrimeRho0Proton"]=false;
    evtDecay["EtaPrimeRho0DeltaP"]=false;
    evtDecay["EtaRhoPDelta0"]=false;
    evtDecay["EtaPrimeRhoPNeutron"]=false;
    evtDecay["EtaPrimeRhoPDelta0"]=false;
    evtDecay["EtaRhoMDeltaPP"]=false;
    evtDecay["EtaKPLambda"]=false;
    evtDecay["EtaK0SigmaP"]=false;
    evtDecay["EtaK0Sigma1385"]=false;
    evtDecay["EtaKStar892PSigma0"]=false;
    evtDecay["EtaEtaPrimeDeltaP"]=false;
    evtDecay["EtaPrimeKPLambda"]=false;
    evtDecay["EtaKStar892PLambda"]=false;
    evtDecay["EtaKPSigma0"]=false;
    evtDecay["EtaEtaDeltaP"]=false;
    evtDecay["EtaEtaProton"]=false;
    evtDecay["EtaEtaPrimeProton"]=false;
    evtDecay["Rho0RhoPDelta0"]=false;
    evtDecay["Rho0Rho0DeltaP"]=false;
    evtDecay["Rho0RhoMDeltaPP"]=false;
    evtDecay["Rho0KPLambda"]=false;
    evtDecay["Rho0KPSigmaStar0"]=false;
    evtDecay["Rho0KStar892PLambda"]=false;
    evtDecay["RhoPK0Lambda"]=false;
    evtDecay["RhoPKStar892Lambda"]=false;
    evtDecay["RhoPK0Sigma0"]=false;
    evtDecay["RhoMKPSigmaP"]=false;
    evtDecay["RhoMKPSigma1385"]=false;
    evtDecay["KPAntiKStar892Neutron"]=false;
    evtDecay["KPAntiKStar892Delta0"]=false;
    evtDecay["KPAntiK0Delta0"]=false;
    evtDecay["KStar892PAntiK0Neutron"]=false;
    evtDecay["KStar892AntiK0DeltaP"]=false;
    evtDecay["K0AntiK0DeltaP"]=false;
    evtDecay["KStar892PAntiK0Delta0"]=false;
    evtDecay["KStar892PAntiKStar892Neutron"]=false;
    evtDecay["RhoPK0SigmaStar0"]=false;
    evtDecay["EtaKPSigmaStar0"]=false;
    evtDecay["KMKStar892DeltaPP"]=false;
    evtDecay["KPKStar892MDeltaP"]=false;
    evtDecay["KPPhi1020Sigma0"]=false;
    evtDecay["K0AntiKStar892Proton"]=false;
    evtDecay["K0AntiKStar892DeltaP"]=false;
    evtDecay["K0KStar892MDeltaPP"]=false;
    evtDecay["AntiKStar892K0Proton"]=false;
    evtDecay["RhoMKStar892PSigmaP"]=false;
}

// return the number of decay channels
int Decay3_pythiaCPP_Omega::Get_nDecays(){
    int ret = 0 ;
    if(evtDecay.size()==evtDecayNum.size()){
        ret = evtDecay.size();
    }else{
        cout<<"Decay3_pythiaCPP_Omega::Get_nDecays() - mismatch in maps evtDecay and evtDecayNum"<<endl;
    }
    return ret;
}

bool Decay3_pythiaCPP_Omega::Check_DecayPID(vector<int> DecayPID, string partName1, string partName2, string partName3){
    int k;
    bool ret = false;
    bool Check = false;
    PDG_pythiaCPP_Omega tempPDG;
    
//    cout<<partName1<<"\t"<<partName2<<"\t"<<partName3<<"\t";
    if(DecayPID.size()==3){
        for(k=0; k<6; k++){
            Check = (tempPDG.Get_id2name(DecayPID[TriCycle[k][0]]).compare(partName1.c_str())==0 && tempPDG.Get_id2name(DecayPID[TriCycle[k][1]]).compare(partName2.c_str())==0 && tempPDG.Get_id2name(DecayPID[TriCycle[k][2]]).compare(partName3.c_str())==0);
            ret = (ret || Check);
        }
    }
//    cout<<ret<<endl;
    return ret;
}

void Decay3_pythiaCPP_Omega::Set_All(vector<int> DecayPID){
    bool foundEvt = false;
    
    if(DecayPID.size()==3){
        foundEvt = evtDecay["RhoMPipProton"] = this->Check_DecayPID(DecayPID,"rho-","pi+","p");
        if(!foundEvt) foundEvt = evtDecay["RhoPPimProton"] = this->Check_DecayPID(DecayPID,"rho+","pi-","p");
        if(!foundEvt) foundEvt = evtDecay["Rho0Pi0Proton"] = this->Check_DecayPID(DecayPID,"rho0","pi0","p");
        if(!foundEvt) foundEvt = evtDecay["RhoPRho0Neutron"] = this->Check_DecayPID(DecayPID,"rho+","rho0","n");
        if(!foundEvt) foundEvt = evtDecay["RhoMPipDeltaP"] = this->Check_DecayPID(DecayPID,"rho-","pi+","Delta+");
        if(!foundEvt) foundEvt = evtDecay["RhoPEtaNeutron"] = this->Check_DecayPID(DecayPID,"rho+","eta","n");
        if(!foundEvt) foundEvt = evtDecay["RhoPRhoMProton"] = this->Check_DecayPID(DecayPID,"rho-","rho+","p");
        if(!foundEvt) foundEvt = evtDecay["Rho0PipDelta0"] = this->Check_DecayPID(DecayPID,"rho0","pi+","Delta0");
        if(!foundEvt) foundEvt = evtDecay["RhoMPi0DeltaPP"] = this->Check_DecayPID(DecayPID,"rho-","pi0","Delta++");
        if(!foundEvt) foundEvt = evtDecay["OmegaPi0Proton"] = this->Check_DecayPID(DecayPID,"omega","pi0","p");
        if(!foundEvt) foundEvt = evtDecay["OmegaPimDeltaPP"] = this->Check_DecayPID(DecayPID,"omega","pi-","Delta++");
        if(!foundEvt) foundEvt = evtDecay["OmegaRhoMDeltaPP"] = this->Check_DecayPID(DecayPID,"omega","rho-","Delta++");
        if(!foundEvt) foundEvt = evtDecay["OmegaPipNeutron"] = this->Check_DecayPID(DecayPID,"omega","pi+","n");
        if(!foundEvt) foundEvt = evtDecay["OmegaOmegaProton"] = this->Check_DecayPID(DecayPID,"omega","omega","p");
        if(!foundEvt) foundEvt = evtDecay["Rho0EtaProton"] = this->Check_DecayPID(DecayPID,"p","rho0","eta");
        if(!foundEvt) foundEvt = evtDecay["RhoMRhoPDeltaP"] = this->Check_DecayPID(DecayPID,"rho+","rho-","Delta+");
        if(!foundEvt) foundEvt = evtDecay["OmegaRhoPNeutron"] = this->Check_DecayPID(DecayPID,"n","rho+","omega");
        if(!foundEvt) foundEvt = evtDecay["OmegaPipDelta0"] = this->Check_DecayPID(DecayPID,"pi+","omega","Delta0");
        if(!foundEvt) foundEvt = evtDecay["PimPi0DeltaPP"] = this->Check_DecayPID(DecayPID,"pi0","pi-","Delta++");
        if(!foundEvt) foundEvt = evtDecay["Pi0EtaProton"] = this->Check_DecayPID(DecayPID,"eta","pi0","p");
        if(!foundEvt) foundEvt = evtDecay["Pi0EtaDelta0"] = this->Check_DecayPID(DecayPID,"Delta0","eta","pi+");
        if(!foundEvt) foundEvt = evtDecay["PipEtaNeutron"] = this->Check_DecayPID(DecayPID,"eta","pi+","n");
        if(!foundEvt) foundEvt = evtDecay["PimEtaDeltaPP"] = this->Check_DecayPID(DecayPID,"eta","pi-","Delta++");
        if(!foundEvt) foundEvt = evtDecay["RhoKPSigma0"] = this->Check_DecayPID(DecayPID,"Sigma0","K+","rho0");
        if(!foundEvt) foundEvt = evtDecay["Pi0Rho0DeltaP"] = this->Check_DecayPID(DecayPID,"pi0","rho0","Delta+");
        if(!foundEvt) foundEvt = evtDecay["PimRhoPDeltaP"] = this->Check_DecayPID(DecayPID,"pi-","rho+","Delta+");
        if(!foundEvt) foundEvt = evtDecay["Pi0RhoPDelta0"] = this->Check_DecayPID(DecayPID,"pi0","rho+","Delta0");
        if(!foundEvt) foundEvt = evtDecay["PipPimDeltaP"] = this->Check_DecayPID(DecayPID,"pi+","pi-","Delta+");
        if(!foundEvt) foundEvt = evtDecay["KStar892PKStar892MProton"] = this->Check_DecayPID(DecayPID,"p","K*(892)-","K*(892)+");
        
        if(!foundEvt) foundEvt = evtDecay["OmegaOmegaDeltaP"] = this->Check_DecayPID(DecayPID,"omega","omega","Delta+");
        if(!foundEvt) foundEvt = evtDecay["OmegaPi0DeltaP"] = this->Check_DecayPID(DecayPID,"omega","pi0","Delta+");
        if(!foundEvt) foundEvt = evtDecay["OmegaPi0Pip"] = this->Check_DecayPID(DecayPID,"pi+","omega","pi0");
        if(!foundEvt) foundEvt = evtDecay["OmegaPimPip"] = this->Check_DecayPID(DecayPID,"pi-","omega","pi+");
        if(!foundEvt) foundEvt = evtDecay["OmegaEtaProton"] = this->Check_DecayPID(DecayPID,"omega","eta","p");
        if(!foundEvt) foundEvt = evtDecay["OmegaEtaDeltaP"] = this->Check_DecayPID(DecayPID,"omega","eta","Delta+");
        if(!foundEvt) foundEvt = evtDecay["OmegaEtaPrimeProton"] = this->Check_DecayPID(DecayPID,"p","eta'(958)","omega");
        if(!foundEvt) foundEvt = evtDecay["OmegaRho0Proton"] = this->Check_DecayPID(DecayPID,"omega","rho0","p");
        if(!foundEvt) foundEvt = evtDecay["OmegaRho0DeltaP"] = this->Check_DecayPID(DecayPID,"omega","rho0","Delta+");
        if(!foundEvt) foundEvt = evtDecay["OmegaRhoPDelta0"] = this->Check_DecayPID(DecayPID,"omega","rho+","Delta0");
        if(!foundEvt) foundEvt = evtDecay["OmegaKPLambda"] = this->Check_DecayPID(DecayPID,"omega","K+","Lambda");
        if(!foundEvt) foundEvt = evtDecay["OmegaKPSigma0"] = this->Check_DecayPID(DecayPID,"omega","K+","Sigma0");
        if(!foundEvt) foundEvt = evtDecay["OmegaKPSigmaStar0"] = this->Check_DecayPID(DecayPID,"omega","K+","Sigma*0");
        if(!foundEvt) foundEvt = evtDecay["OmegaKStar892PLambda"] = this->Check_DecayPID(DecayPID,"omega","K*(892)+","Lambda");
        if(!foundEvt) foundEvt = evtDecay["OmegaK0SigmaP"] = this->Check_DecayPID(DecayPID,"omega","K0","Sigma+");
        if(!foundEvt) foundEvt = evtDecay["OmegaK0Sigma1385"] = this->Check_DecayPID(DecayPID,"omega","K0","Sigma(1385)");
        if(!foundEvt) foundEvt = evtDecay["OmegaPi0Eta"] = this->Check_DecayPID(DecayPID,"omega","pi0","eta");
        if(!foundEvt) foundEvt = evtDecay["OmegaPi0RhoP"] = this->Check_DecayPID(DecayPID,"pi0","omega","rho+");
        if(!foundEvt) foundEvt = evtDecay["Pi0PipDelta0"] = this->Check_DecayPID(DecayPID,"pi+","pi0","Delta0");
        if(!foundEvt) foundEvt = evtDecay["Pi0PimRhoP"] = this->Check_DecayPID(DecayPID,"pi-","pi0","rho+");
        if(!foundEvt) foundEvt = evtDecay["Pi0EtaDeltaP"] = this->Check_DecayPID(DecayPID,"pi0","eta","Delta+");
        if(!foundEvt) foundEvt = evtDecay["Pi0EtaPrimeProton"] = this->Check_DecayPID(DecayPID,"pi0","eta'(958)","p");
        if(!foundEvt) foundEvt = evtDecay["Pi0EtaPrimeDeltaP"] = this->Check_DecayPID(DecayPID,"eta'(958)","pi0","Delta+");
        if(!foundEvt) foundEvt = evtDecay["Pi0K0Sigma1385"] = this->Check_DecayPID(DecayPID,"pi0","K0","Sigma(1385)");
        if(!foundEvt) foundEvt = evtDecay["Pi0KStar892PLambda"] = this->Check_DecayPID(DecayPID,"pi0","K*(892)+","Lambda");
        if(!foundEvt) foundEvt = evtDecay["Pi0KStar892SigmaP"] = this->Check_DecayPID(DecayPID,"pi0","K*(892)0","Sigma+");
        if(!foundEvt) foundEvt = evtDecay["Pi0KStar892PSigma0"] = this->Check_DecayPID(DecayPID,"pi0","K*(892)+","Sigma0");
        if(!foundEvt) foundEvt = evtDecay["Pi0KStar892PSigmaStar0"] = this->Check_DecayPID(DecayPID,"pi0","K*(892)+","Sigma*0");
        if(!foundEvt) foundEvt = evtDecay["Pi0PipPim"] = this->Check_DecayPID(DecayPID,"pi0","pi+","pi-");
        if(!foundEvt) foundEvt = evtDecay["Pi0PipEta"] = this->Check_DecayPID(DecayPID,"pi0","pi+","eta");
        if(!foundEvt) foundEvt = evtDecay["Pi0PipRhoM"] = this->Check_DecayPID(DecayPID,"pi0","pi+","rho-");
        if(!foundEvt) foundEvt = evtDecay["Pi0K0SigmaP"] = this->Check_DecayPID(DecayPID,"pi0","K0","Sigma+");
        if(!foundEvt) foundEvt = evtDecay["Pi0KPSigmaStar0"] = this->Check_DecayPID(DecayPID,"pi0","K+","Sigma*0");
        if(!foundEvt) foundEvt = evtDecay["PimPipRhoP"] = this->Check_DecayPID(DecayPID,"pi+","pi-","rho+");
        if(!foundEvt) foundEvt = evtDecay["PimKStar892PSigmaP"] = this->Check_DecayPID(DecayPID,"pi-","K*(892)+","Sigma+");
        if(!foundEvt) foundEvt = evtDecay["PimEtaPrimeDeltaPP"] = this->Check_DecayPID(DecayPID,"pi-","eta'(958)","Delta++");
        if(!foundEvt) foundEvt = evtDecay["PimKPSigma1385"] = this->Check_DecayPID(DecayPID,"pi-","K+","Sigma(1385)");
        if(!foundEvt) foundEvt = evtDecay["PimKStar892PSigma1385"] = this->Check_DecayPID(DecayPID,"pi-","K*(892)+","Sigma(1385)");
        if(!foundEvt) foundEvt = evtDecay["PipEtaPrimeNeutron"] = this->Check_DecayPID(DecayPID,"pi+","eta'(958)","n");
        if(!foundEvt) foundEvt = evtDecay["PipEtaPrimeDelta0"] = this->Check_DecayPID(DecayPID,"pi+","eta'(958)","Delta0");
        if(!foundEvt) foundEvt = evtDecay["PipK0Lambda"] = this->Check_DecayPID(DecayPID,"pi+","K0","Lambda");
        if(!foundEvt) foundEvt = evtDecay["PipK0Sigma0"] = this->Check_DecayPID(DecayPID,"pi+","K0","Sigma0");
        if(!foundEvt) foundEvt = evtDecay["PipK0SigmaStar0"] = this->Check_DecayPID(DecayPID,"pi+","K0","Sigma*0");
        if(!foundEvt) foundEvt = evtDecay["PipKStar892Lambda"] = this->Check_DecayPID(DecayPID,"pi+","K*(892)0","Lambda");
        if(!foundEvt) foundEvt = evtDecay["PipKStar892Sigma0"] = this->Check_DecayPID(DecayPID,"pi+","K*(892)0","Sigma0");
        if(!foundEvt) foundEvt = evtDecay["PipKStar892SigmaStar0"] = this->Check_DecayPID(DecayPID,"pi+","K*(892)0","Sigma*0");
        if(!foundEvt) foundEvt = evtDecay["EtaRho0DeltaP"] = this->Check_DecayPID(DecayPID,"eta","rho0","Delta+");
        if(!foundEvt) foundEvt = evtDecay["EtaPrimeRho0Proton"] = this->Check_DecayPID(DecayPID,"eta'(958)","rho0","p");
        if(!foundEvt) foundEvt = evtDecay["EtaPrimeRho0DeltaP"] = this->Check_DecayPID(DecayPID,"Delta+","rho0","eta'(958)");
        if(!foundEvt) foundEvt = evtDecay["EtaK0SigmaP"] = this->Check_DecayPID(DecayPID,"eta","K0","Sigma+");
        if(!foundEvt) foundEvt = evtDecay["EtaK0Sigma1385"] = this->Check_DecayPID(DecayPID,"eta","K0","Sigma(1385)");
        if(!foundEvt) foundEvt = evtDecay["EtaKPSigmaStar0"] = this->Check_DecayPID(DecayPID,"eta","K+","Sigma*0");
        if(!foundEvt) foundEvt = evtDecay["EtaKStar892PSigma0"] = this->Check_DecayPID(DecayPID,"eta","K*(892)+","Sigma0");
        if(!foundEvt) foundEvt = evtDecay["EtaEtaPrimeDeltaP"] = this->Check_DecayPID(DecayPID,"eta'(958)","eta","Delta+");
        if(!foundEvt) foundEvt = evtDecay["EtaRhoPDelta0"] = this->Check_DecayPID(DecayPID,"eta","rho+","Delta0");
        if(!foundEvt) foundEvt = evtDecay["EtaPrimeRhoPNeutron"] = this->Check_DecayPID(DecayPID,"eta'(958)","rho+","n");
        if(!foundEvt) foundEvt = evtDecay["EtaPrimeRhoPDelta0"] = this->Check_DecayPID(DecayPID,"Delta0","rho+","eta'(958)");
        if(!foundEvt) foundEvt = evtDecay["EtaRhoMDeltaPP"] = this->Check_DecayPID(DecayPID,"eta","rho-","Delta++");
        if(!foundEvt) foundEvt = evtDecay["EtaKPLambda"] = this->Check_DecayPID(DecayPID,"eta","K+","Lambda");
        if(!foundEvt) foundEvt = evtDecay["EtaPrimeKPLambda"] = this->Check_DecayPID(DecayPID,"K+","eta'(958)","Lambda");
        if(!foundEvt) foundEvt = evtDecay["EtaKStar892PLambda"] = this->Check_DecayPID(DecayPID,"Lambda","eta","K*(892)+");
        if(!foundEvt) foundEvt = evtDecay["EtaKPSigma0"] = this->Check_DecayPID(DecayPID,"eta","K+","Sigma0");
        if(!foundEvt) foundEvt = evtDecay["EtaEtaDeltaP"] = this->Check_DecayPID(DecayPID,"eta","eta","Delta+");
        if(!foundEvt) foundEvt = evtDecay["EtaEtaProton"] = this->Check_DecayPID(DecayPID,"eta","eta","p");
        if(!foundEvt) foundEvt = evtDecay["EtaEtaPrimeProton"] = this->Check_DecayPID(DecayPID,"eta","eta'(958)","p");
        if(!foundEvt) foundEvt = evtDecay["Rho0RhoPDelta0"] = this->Check_DecayPID(DecayPID,"rho+","rho0","Delta0");
        if(!foundEvt) foundEvt = evtDecay["Rho0Rho0DeltaP"] = this->Check_DecayPID(DecayPID,"rho0","rho0","Delta+");
        if(!foundEvt) foundEvt = evtDecay["Rho0RhoMDeltaPP"] = this->Check_DecayPID(DecayPID,"rho0","rho-","Delta++");
        if(!foundEvt) foundEvt = evtDecay["Rho0KPLambda"] = this->Check_DecayPID(DecayPID,"Lambda","K+","rho0");
        if(!foundEvt) foundEvt = evtDecay["Rho0KPSigmaStar0"] = this->Check_DecayPID(DecayPID,"Sigma*0","K+","rho0");
        if(!foundEvt) foundEvt = evtDecay["Rho0KStar892PLambda"] = this->Check_DecayPID(DecayPID,"Lambda","K*(892)+","rho0");
        if(!foundEvt) foundEvt = evtDecay["RhoPK0Lambda"] = this->Check_DecayPID(DecayPID,"rho+","K0","Lambda");
        if(!foundEvt) foundEvt = evtDecay["RhoPKStar892Lambda"] = this->Check_DecayPID(DecayPID,"rho+","K*(892)0","Lambda");
        if(!foundEvt) foundEvt = evtDecay["RhoPK0Sigma0"] = this->Check_DecayPID(DecayPID,"rho+","K0","Sigma0");
        if(!foundEvt) foundEvt = evtDecay["RhoPK0SigmaStar0"] = this->Check_DecayPID(DecayPID,"rho+","K0","Sigma*0");
        if(!foundEvt) foundEvt = evtDecay["RhoMKPSigmaP"] = this->Check_DecayPID(DecayPID,"rho-","K+","Sigma+");
        if(!foundEvt) foundEvt = evtDecay["RhoMKPSigma1385"] = this->Check_DecayPID(DecayPID,"rho-","K+","Sigma(1385)");
        if(!foundEvt) foundEvt = evtDecay["KPAntiKStar892Neutron"] = this->Check_DecayPID(DecayPID,"K+","anti-K*(892)0","n");
        if(!foundEvt) foundEvt = evtDecay["KPAntiKStar892Delta0"] = this->Check_DecayPID(DecayPID,"K+","anti-K*(892)0","Delta0");
        if(!foundEvt) foundEvt = evtDecay["KPAntiK0Delta0"] = this->Check_DecayPID(DecayPID,"K+","anti-K0","Delta0");
        if(!foundEvt) foundEvt = evtDecay["KStar892PAntiK0Neutron"] = this->Check_DecayPID(DecayPID,"K*(892)+","anti-K0","n");
        if(!foundEvt) foundEvt = evtDecay["KStar892AntiK0DeltaP"] = this->Check_DecayPID(DecayPID,"K*(892)0","anti-K0","Delta+");
        if(!foundEvt) foundEvt = evtDecay["K0AntiK0DeltaP"] = this->Check_DecayPID(DecayPID,"K0","anti-K0","Delta+");
        if(!foundEvt) foundEvt = evtDecay["KStar892PAntiK0Delta0"] = this->Check_DecayPID(DecayPID,"Delta0","anti-K0","K*(892)+");
        if(!foundEvt) foundEvt = evtDecay["KStar892PAntiKStar892Neutron"] = this->Check_DecayPID(DecayPID,"K*(892)+","anti-K*(892)0","n");
        if(!foundEvt) foundEvt = evtDecay["KMKStar892DeltaPP"] = this->Check_DecayPID(DecayPID,"K*(892)0","K-","Delta++");
        if(!foundEvt) foundEvt = evtDecay["KPKStar892MDeltaP"] = this->Check_DecayPID(DecayPID,"K+","K*(892)-","Delta+");
        if(!foundEvt) foundEvt = evtDecay["KPPhi1020Sigma0"] = this->Check_DecayPID(DecayPID,"K+","phi(1020)","Sigma0");
        if(!foundEvt) foundEvt = evtDecay["K0AntiKStar892Proton"] = this->Check_DecayPID(DecayPID,"K0","anti-K*(892)0","p");
        if(!foundEvt) foundEvt = evtDecay["K0AntiKStar892DeltaP"] = this->Check_DecayPID(DecayPID,"K0","anti-K*(892)0","Delta+");
        if(!foundEvt) foundEvt = evtDecay["K0KStar892MDeltaPP"] = this->Check_DecayPID(DecayPID,"K0","K*(892)-","Delta++");
        if(!foundEvt) foundEvt = evtDecay["AntiKStar892K0Proton"] = this->Check_DecayPID(DecayPID,"K*(892)0","anti-K0","p");
        if(!foundEvt) foundEvt = evtDecay["RhoMKStar892PSigmaP"] = this->Check_DecayPID(DecayPID,"Sigma+","K*(892)+","rho-");
    }else{
        cout<<"Decay3_pythiaCPP_Omega::Set_All - Wrong number of particles.  Must be 3."<<endl;
    }
}

bool Decay3_pythiaCPP_Omega::Get_All(){
    bool ret = false;
    
    for (std::map<string,bool>::iterator it=evtDecay.begin(); it!=evtDecay.end(); ++it){
        if(it->second){
            ret = true;
            break;
        }
    }
    return ret;
}

bool Decay3_pythiaCPP_Omega::Get_evtDecay(string evtName){
    bool ret = false;
    
    if(evtDecay.count(evtName)>0)
        ret = evtDecay[evtName];
    else
        cout<<"Decay3_pythiaCPP_Omega::Get_evtDecay() - no match for "<<evtName<<endl;
    
    return ret;
}

// returns the name of the decay topology that triggered the event
string Decay3_pythiaCPP_Omega::Get_trueDecay(){
    string ret;
    
    for (std::map<string,bool>::iterator it=evtDecay.begin(); it!=evtDecay.end(); ++it){
        if(it->second){
            ret = it->first;
            break;
        }
    }
    return ret;
}

// returns the histogram number for the decay topology that triggered the event
int Decay3_pythiaCPP_Omega::Get_trueDecayNum(){
    int ret = 0;
    
    for (std::map<string,bool>::iterator it=evtDecay.begin(); it!=evtDecay.end(); ++it){
        if(it->second){
            ret = this->Get_evtDecayNum(it->first);
            break;
        }
    }
    return ret;
}

// returns the histogram number given the decay name
int Decay3_pythiaCPP_Omega::Get_evtDecayNum(string evtName){
    int ret = 0;
    
    if(evtDecayNum.count(evtName)>0)
        ret = evtDecayNum[evtName];
    else
        cout<<"Decay3_pythiaCPP_Omega::Get_evtDecayNum() - no match for "<<evtName<<endl;
    
    return ret;
}

void Decay3_pythiaCPP_Omega::Print_HistNumbers(){
    
    cout<<"Decay3_pythiaCPP_Omega::Print_HistNumbers()"<<endl;
    for (std::map<string,int>::iterator it=evtDecayNum.begin(); it!=evtDecayNum.end(); ++it){
        cout<<it->first<<"\t"<<it->second<<endl;
    }
}



