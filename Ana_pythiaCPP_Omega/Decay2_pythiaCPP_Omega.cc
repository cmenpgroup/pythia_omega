#include "Decay2_pythiaCPP_Omega.h"

Decay2_pythiaCPP_Omega::Decay2_pythiaCPP_Omega(){
    this->Init();

    // fill the map with the histogram numbers
    int nhist = 1; // start the histogram list at 1
    for (std::map<string,bool>::iterator it=evtDecay.begin(); it!=evtDecay.end(); ++it){
        evtDecayNum[it->first] = nhist;
        nhist++;
    }
}

void Decay2_pythiaCPP_Omega::Init(){
    evtDecay["RhoPDelta0"]=false;
    evtDecay["Rho0DeltaP"]=false;
    evtDecay["RhoMDeltaPP"]=false;
    evtDecay["Rho0Pi0"]=false;
    evtDecay["RhoPPim"]=false;
    evtDecay["Pi0Eta"]=false;
    evtDecay["KStar892Sigma1385"]=false;
    evtDecay["KStar892PLambda"]=false;
    evtDecay["EtaProton"]=false;
    evtDecay["EtaPrimeProton"]=false;
    evtDecay["OmegaProton"]=false;
    evtDecay["OmegaDeltaP"]=false;
    evtDecay["RhoMPip"]=false;
    evtDecay["EtaPrimeDeltaP"]=false;
    evtDecay["EtaDeltaP"]=false;
    evtDecay["EtaEta"]=false;
    evtDecay["KStar892PSigma0"]=false;
    evtDecay["KStar892PSigmaStar0"]=false;
    evtDecay["OmegaPi0"]=false;
    evtDecay["OmegaEta"]=false;
    evtDecay["OmegaPip"]=false;
    evtDecay["Rho0Eta"]=false;
    evtDecay["RhoPK0"]=false;
    evtDecay["KStar892SigmaP"]=false;
    evtDecay["K0Sigma1385"]=false;
    evtDecay["OmegaRho0"]=false;
    evtDecay["K0SigmaP"]=false;
    evtDecay["AntiKStar892Neutron"]=false;
    evtDecay["AntiKStar892K0"]=false;
    evtDecay["RhoPRhoM"]=false;
    evtDecay["EtaEtaPrime"]=false;
    evtDecay["Pi0EtaPrime"]=false;
    evtDecay["PipEtaPrime"]=false;
    evtDecay["RhoPEta"]=false;
    evtDecay["RhoPRho0"]=false;
    evtDecay["OmegaKP"]=false;
    evtDecay["OmegaRhoP"]=false;
    evtDecay["AntiK0KStar892P"]=false;
    evtDecay["PipKStar892"]=false;
    evtDecay["Rho0EtaPrime"]=false;
    evtDecay["OmegaOmega"]=false;
    evtDecay["PipEta"]=false;
    evtDecay["EtaSigma1385"]=false;
    evtDecay["RhoMPi0"]=false;
    evtDecay["OmegaEtaPrime"]=false;
    evtDecay["AntiK0Delta0"]=false;
    evtDecay["AntiK0DeltaP"]=false;
    evtDecay["Rho0Rho0"]=false;
    evtDecay["AntiK0Proton"]=false;
    evtDecay["Pi0DeltaP"]=false;
    evtDecay["KStar892MProton"]=false;
    evtDecay["KStar892MDeltaPP"]=false;
    evtDecay["KStar892MKStar892P"]=false;
    evtDecay["AntiKStar892KP"]=false;
    evtDecay["AntiK0K0"]=false;
    evtDecay["PimEta"]=false;
    evtDecay["KStar892PPi0"]=false;
    evtDecay["RhoPNeutron"]=false;
    evtDecay["EtaKP"]=false;
    evtDecay["AntiKStar892DeltaP"]=false;
    evtDecay["EtaLambda"]=false;
    evtDecay["PipDelta0"]=false;
    evtDecay["OmegaPim"]=false;
    evtDecay["Pi0Proton"]=false;
    evtDecay["PipK0"]=false;
    evtDecay["PipNeutron"]=false;
    evtDecay["PimDeltaPP"]=false;
    evtDecay["AntiK0KStar892"]=false;
    evtDecay["Rho0Proton"]=false;
    evtDecay["KPLambda"]=false;
}

// return the number of decay channels
int Decay2_pythiaCPP_Omega::Get_nDecays(){
    int ret = 0 ;
    if(evtDecay.size()==evtDecayNum.size()){
        ret = evtDecay.size();
    }else{
        cout<<"Decay2_pythiaCPP_Omega::Get_nDecays() - mismatch in maps evtDecay and evtDecayNum"<<endl;
    }
    return ret;
}

bool Decay2_pythiaCPP_Omega::Check_DecayPID(vector<int> DecayPID, string partName1, string partName2){
    int k;
    bool ret = false;
    bool Check = false;
    PDG_pythiaCPP_Omega tempPDG;
    
    if(DecayPID.size()==2){
        for(k=0; k<2; k++){
            Check = (tempPDG.Get_id2name(DecayPID[k]).compare(partName1.c_str())==0 && tempPDG.Get_id2name(DecayPID[1-k]).compare(partName2.c_str())==0);
            ret = (ret || Check);
        }
    }
    return ret;
}

void Decay2_pythiaCPP_Omega::Set_All(vector<int> DecayPID){
    bool foundEvt = false;
    
    if(DecayPID.size()==2){
        foundEvt = evtDecay["AntiKStar892K0"] = this->Check_DecayPID(DecayPID,"K0","anti-K*(892)0");
        if(!foundEvt) foundEvt = evtDecay["AntiKStar892Neutron"] = this->Check_DecayPID(DecayPID,"n","anti-K*(892)0");
        if(!foundEvt) foundEvt = evtDecay["K0Sigma1385"] = this->Check_DecayPID(DecayPID,"K0","Sigma(1385)");
        if(!foundEvt) foundEvt = evtDecay["K0SigmaP"] = this->Check_DecayPID(DecayPID,"K0","Sigma+");
        if(!foundEvt) foundEvt = evtDecay["KStar892Sigma1385"] = this->Check_DecayPID(DecayPID,"K*(892)0","Sigma(1385)");
        if(!foundEvt) foundEvt = evtDecay["KStar892SigmaP"] = this->Check_DecayPID(DecayPID,"K*(892)0","Sigma+");
        if(!foundEvt) foundEvt = evtDecay["KStar892PSigma0"] = this->Check_DecayPID(DecayPID,"K*(892)+","Sigma0");
        if(!foundEvt) foundEvt = evtDecay["KStar892PSigmaStar0"] = this->Check_DecayPID(DecayPID,"K*(892)+","Sigma*0");
        if(!foundEvt) foundEvt = evtDecay["Rho0Pi0"] = this->Check_DecayPID(DecayPID,"rho0","pi0");
        if(!foundEvt) foundEvt = evtDecay["RhoPK0"] = this->Check_DecayPID(DecayPID,"rho+","K0");
        if(!foundEvt) foundEvt = evtDecay["RhoPPim"] = this->Check_DecayPID(DecayPID,"rho+","pi-");
        if(!foundEvt) foundEvt = evtDecay["RhoMPip"] = this->Check_DecayPID(DecayPID,"rho-","pi+");
        if(!foundEvt) foundEvt = evtDecay["RhoMDeltaPP"] = this->Check_DecayPID(DecayPID,"rho-","Delta++");
        if(!foundEvt) foundEvt = evtDecay["Rho0DeltaP"] = this->Check_DecayPID(DecayPID,"rho0","Delta+");
        if(!foundEvt) foundEvt = evtDecay["RhoPDelta0"] = this->Check_DecayPID(DecayPID,"rho+","Delta0");
        if(!foundEvt) foundEvt = evtDecay["Pi0Eta"] = this->Check_DecayPID(DecayPID,"pi0","eta");
        if(!foundEvt) foundEvt = evtDecay["KStar892PLambda"] = this->Check_DecayPID(DecayPID,"K*(892)+","Lambda");
        if(!foundEvt) foundEvt = evtDecay["EtaEta"] = this->Check_DecayPID(DecayPID,"eta","eta");
        if(!foundEvt) foundEvt = evtDecay["EtaProton"] = this->Check_DecayPID(DecayPID,"p","eta");
        if(!foundEvt) foundEvt = evtDecay["EtaDeltaP"] = this->Check_DecayPID(DecayPID,"eta","Delta+");
        if(!foundEvt) foundEvt = evtDecay["EtaPrimeProton"] = this->Check_DecayPID(DecayPID,"p","eta'(958)");
        if(!foundEvt) foundEvt = evtDecay["EtaPrimeDeltaP"] = this->Check_DecayPID(DecayPID,"Delta+","eta'(958)");
        if(!foundEvt) foundEvt = evtDecay["OmegaDeltaP"] = this->Check_DecayPID(DecayPID,"omega","Delta+");
        if(!foundEvt) foundEvt = evtDecay["OmegaProton"] = this->Check_DecayPID(DecayPID,"omega","p");
        if(!foundEvt) foundEvt = evtDecay["OmegaPi0"] = this->Check_DecayPID(DecayPID,"omega","pi0");
        if(!foundEvt) foundEvt = evtDecay["OmegaPip"] = this->Check_DecayPID(DecayPID,"omega","pi+");
        if(!foundEvt) foundEvt = evtDecay["OmegaPim"] = this->Check_DecayPID(DecayPID,"omega","pi-");
        if(!foundEvt) foundEvt = evtDecay["OmegaEta"] = this->Check_DecayPID(DecayPID,"omega","eta");
        if(!foundEvt) foundEvt = evtDecay["OmegaEtaPrime"] = this->Check_DecayPID(DecayPID,"omega","eta'(958)");
        if(!foundEvt) foundEvt = evtDecay["OmegaRho0"] = this->Check_DecayPID(DecayPID,"omega","rho0");
        if(!foundEvt) foundEvt = evtDecay["OmegaKP"] = this->Check_DecayPID(DecayPID,"omega","K+");
        if(!foundEvt) foundEvt = evtDecay["OmegaRhoP"] = this->Check_DecayPID(DecayPID,"omega","rho+");
        if(!foundEvt) foundEvt = evtDecay["OmegaOmega"] = this->Check_DecayPID(DecayPID,"omega","omega");
        if(!foundEvt) foundEvt = evtDecay["Rho0Eta"] = this->Check_DecayPID(DecayPID,"rho0","eta");
        if(!foundEvt) foundEvt = evtDecay["RhoPRhoM"] = this->Check_DecayPID(DecayPID,"rho+","rho-");
        if(!foundEvt) foundEvt = evtDecay["RhoPRho0"] = this->Check_DecayPID(DecayPID,"rho+","rho0");
        if(!foundEvt) foundEvt = evtDecay["RhoPEta"] = this->Check_DecayPID(DecayPID,"rho+","eta");
        if(!foundEvt) foundEvt = evtDecay["RhoPNeutron"] = this->Check_DecayPID(DecayPID,"rho+","n");
        if(!foundEvt) foundEvt = evtDecay["RhoMPi0"] = this->Check_DecayPID(DecayPID,"rho-","pi0");
        if(!foundEvt) foundEvt = evtDecay["Rho0EtaPrime"] = this->Check_DecayPID(DecayPID,"rho0","eta'(958)");
        if(!foundEvt) foundEvt = evtDecay["Rho0Rho0"] = this->Check_DecayPID(DecayPID,"rho0","rho0");
        if(!foundEvt) foundEvt = evtDecay["EtaEtaPrime"] = this->Check_DecayPID(DecayPID,"eta","eta'(958)");
        if(!foundEvt) foundEvt = evtDecay["Pi0EtaPrime"] = this->Check_DecayPID(DecayPID,"pi0","eta'(958)");
        if(!foundEvt) foundEvt = evtDecay["PipEtaPrime"] = this->Check_DecayPID(DecayPID,"pi+","eta'(958)");
        if(!foundEvt) foundEvt = evtDecay["PipEta"] = this->Check_DecayPID(DecayPID,"pi+","eta");
        if(!foundEvt) foundEvt = evtDecay["PimEta"] = this->Check_DecayPID(DecayPID,"pi-","eta");
        if(!foundEvt) foundEvt = evtDecay["AntiK0Proton"] = this->Check_DecayPID(DecayPID,"anti-K0","p");
        if(!foundEvt) foundEvt = evtDecay["AntiK0Delta0"] = this->Check_DecayPID(DecayPID,"anti-K0","Delta0");
        if(!foundEvt) foundEvt = evtDecay["AntiK0DeltaP"] = this->Check_DecayPID(DecayPID,"anti-K0","Delta+");
        if(!foundEvt) foundEvt = evtDecay["AntiK0KStar892P"] = this->Check_DecayPID(DecayPID,"anti-K0","K*(892)+");
        if(!foundEvt) foundEvt = evtDecay["AntiK0K0"] = this->Check_DecayPID(DecayPID,"anti-K0","K0");
        if(!foundEvt) foundEvt = evtDecay["PipKStar892"] = this->Check_DecayPID(DecayPID,"pi+","K*(892)0");
        if(!foundEvt) foundEvt = evtDecay["EtaSigma1385"] = this->Check_DecayPID(DecayPID,"eta","Sigma(1385)");
        if(!foundEvt) foundEvt = evtDecay["KStar892MProton"] = this->Check_DecayPID(DecayPID,"K*(892)-","p");
        if(!foundEvt) foundEvt = evtDecay["KStar892MDeltaPP"] = this->Check_DecayPID(DecayPID,"K*(892)-","Delta++");
        if(!foundEvt) foundEvt = evtDecay["KStar892MKStar892P"] = this->Check_DecayPID(DecayPID,"K*(892)-","K*(892)+");
        if(!foundEvt) foundEvt = evtDecay["KStar892PPi0"] = this->Check_DecayPID(DecayPID,"K*(892)+","pi0");
        if(!foundEvt) foundEvt = evtDecay["AntiKStar892KP"] = this->Check_DecayPID(DecayPID,"anti-K*(892)0","K+");
        if(!foundEvt) foundEvt = evtDecay["AntiKStar892DeltaP"] = this->Check_DecayPID(DecayPID,"anti-K*(892)0","Delta+");
        if(!foundEvt) foundEvt = evtDecay["Pi0DeltaP"] = this->Check_DecayPID(DecayPID,"pi0","Delta+");
        if(!foundEvt) foundEvt = evtDecay["PipDelta0"] = this->Check_DecayPID(DecayPID,"pi+","Delta0");
        if(!foundEvt) foundEvt = evtDecay["PipNeutron"] = this->Check_DecayPID(DecayPID,"pi+","n");
        if(!foundEvt) foundEvt = evtDecay["PimDeltaPP"] = this->Check_DecayPID(DecayPID,"pi-","Delta++");
        if(!foundEvt) foundEvt = evtDecay["PipK0"] = this->Check_DecayPID(DecayPID,"pi+","K0");
        if(!foundEvt) foundEvt = evtDecay["Pi0Proton"] = this->Check_DecayPID(DecayPID,"pi0","p");
        if(!foundEvt) foundEvt = evtDecay["EtaKP"] = this->Check_DecayPID(DecayPID,"eta","K+");
        if(!foundEvt) foundEvt = evtDecay["EtaLambda"] = this->Check_DecayPID(DecayPID,"eta","Lambda");
        if(!foundEvt) foundEvt = evtDecay["AntiK0KStar892"] = this->Check_DecayPID(DecayPID,"anti-K0","K*(892)0");
        if(!foundEvt) foundEvt = evtDecay["Rho0Proton"] = this->Check_DecayPID(DecayPID,"rho0","p");
        if(!foundEvt) foundEvt = evtDecay["KPLambda"] = this->Check_DecayPID(DecayPID,"K+","Lambda");
    }else{
        cout<<"Decay2_pythiaCPP_Omega::Set_All - Wrong number of particles.  Must be 2."<<endl;
    }
}

bool Decay2_pythiaCPP_Omega::Get_All(){
    bool ret = false;
    
    for (std::map<string,bool>::iterator it=evtDecay.begin(); it!=evtDecay.end(); ++it){
        if(it->second){
            ret = true;
            break;
        }
    }
    return ret;
}

bool Decay2_pythiaCPP_Omega::Get_evtDecay(string evtName){
    bool ret = false;
    
    if(evtDecay.count(evtName)>0)
        ret = evtDecay[evtName];
    else
        cout<<"Decay2_pythiaCPP_Omega::Get_evtDecay() - no match for "<<evtName<<endl;
                                          
    return ret;
}

// returns the name of the decay topology that triggered the event
string Decay2_pythiaCPP_Omega::Get_trueDecay(){
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
int Decay2_pythiaCPP_Omega::Get_trueDecayNum(){
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
int Decay2_pythiaCPP_Omega::Get_evtDecayNum(string evtName){
    int ret = 0;
    
    if(evtDecayNum.count(evtName)>0)
        ret = evtDecayNum[evtName];
    else
        cout<<"Decay2_pythiaCPP_Omega::Get_evtDecayNum() - no match for "<<evtName<<endl;
    
    return ret;
}

void Decay2_pythiaCPP_Omega::Print_HistNumbers(){
    
    cout<<"Decay2_pythiaCPP_Omega::Print_HistNumbers()"<<endl;
    for (std::map<string,int>::iterator it=evtDecayNum.begin(); it!=evtDecayNum.end(); ++it){
        cout<<it->first<<"\t"<<it->second<<endl;
    }
}
