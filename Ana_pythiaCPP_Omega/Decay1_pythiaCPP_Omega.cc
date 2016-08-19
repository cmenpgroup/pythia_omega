#include "Decay1_pythiaCPP_Omega.h"

Decay1_pythiaCPP_Omega::Decay1_pythiaCPP_Omega(){
    this->Init();
    
    // fill the map with the histogram numbers
    int nhist = 1; // start the histogram list at 1
    for (std::map<string,bool>::iterator it=evtDecay.begin(); it!=evtDecay.end(); ++it){
        evtDecayNum[it->first] = nhist;
        nhist++;
    }
}

// initialize all decay topology flags to false
void Decay1_pythiaCPP_Omega::Init(){
    evtDecay["RhoM"]=false;
    evtDecay["Eta"]=false;
    evtDecay["EtaPrime"]=false;
    evtDecay["Omega"]=false;
    evtDecay["KStar892P"]=false;
    evtDecay["Lambda"]=false;
    evtDecay["SigmaP"]=false;
    evtDecay["Sigma1385"]=false;
}

// return the number of decay channels
int Decay1_pythiaCPP_Omega::Get_nDecays(){
    int ret = 0 ;
    if(evtDecay.size()==evtDecayNum.size()){
        ret = evtDecay.size();
    }else{
        cout<<"Decay1_pythiaCPP_Omega::Get_nDecays() - mismatch in maps evtDecay and evtDecayNum"<<endl;
    }
    return ret;
}

// check whether decay particle id. for the event
bool Decay1_pythiaCPP_Omega::Check_DecayPID(vector<int> DecayPID, string partName1){
    bool ret = false;
    PDG_pythiaCPP_Omega tempPDG;
    
    if(DecayPID.size()==1){
        ret = (tempPDG.Get_id2name(DecayPID[0]).compare(partName1.c_str())==0);
    }

    return ret;
}

// set the flags for each single-particle decay topology
void Decay1_pythiaCPP_Omega::Set_All(vector<int> DecayPID){
    bool foundEvt = false;
    
    if(DecayPID.size()==1){
        foundEvt = evtDecay["RhoM"] = this->Check_DecayPID(DecayPID,"rho-");
        if(!foundEvt) foundEvt = evtDecay["Eta"] = this->Check_DecayPID(DecayPID,"eta");
        if(!foundEvt) foundEvt = evtDecay["EtaPrime"] = this->Check_DecayPID(DecayPID,"eta'(958)");
        if(!foundEvt) foundEvt = evtDecay["Omega"] = this->Check_DecayPID(DecayPID,"omega");
        if(!foundEvt) foundEvt = evtDecay["KStar892P"] = this->Check_DecayPID(DecayPID,"K*(892)+");
        if(!foundEvt) foundEvt = evtDecay["Lambda"] = this->Check_DecayPID(DecayPID,"Lambda");
        if(!foundEvt) foundEvt = evtDecay["SigmaP"] = this->Check_DecayPID(DecayPID,"Sigma+");
        if(!foundEvt) foundEvt = evtDecay["Sigma1385"] = this->Check_DecayPID(DecayPID,"Sigma(1385)");
    }else{
        cout<<"Decay1_pythiaCPP_Omega::Set_All - Wrong number of particles.  Must be 1."<<endl;
    }
}

// returns true or false whether the event contained a one particle decay
bool Decay1_pythiaCPP_Omega::Get_All(){
    bool ret = false;
    
    for (std::map<string,bool>::iterator it=evtDecay.begin(); it!=evtDecay.end(); ++it){
        if(it->second){
            ret = true;
            break;
        }
    }
    return ret;
}

// returns true or false whether given the decay was the event
bool Decay1_pythiaCPP_Omega::Get_evtDecay(string evtName){
    bool ret = false;
    
    if(evtDecay.count(evtName)>0)
        ret = evtDecay[evtName];
    else
        cout<<"Decay1_pythiaCPP_Omega::Get_evtDecay() - no match for "<<evtName<<endl;
    
    return ret;
}

// returns the name of the decay topology that triggered the event
string Decay1_pythiaCPP_Omega::Get_trueDecay(){
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
int Decay1_pythiaCPP_Omega::Get_trueDecayNum(){
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
int Decay1_pythiaCPP_Omega::Get_evtDecayNum(string evtName){
    int ret = 0;
    
    if(evtDecayNum.count(evtName)>0)
        ret = evtDecayNum[evtName];
    else
        cout<<"Decay1_pythiaCPP_Omega::Get_evtDecayNum() - no match for "<<evtName<<endl;
    
    return ret;
}

void Decay1_pythiaCPP_Omega::Print_HistNumbers(){
    
    cout<<"Decay1_pythiaCPP_Omega::Print_HistNumbers()"<<endl;
    for (std::map<string,int>::iterator it=evtDecayNum.begin(); it!=evtDecayNum.end(); ++it){
        cout<<it->first<<"\t"<<it->second<<endl;
    }
}

