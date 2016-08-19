#ifndef DECAY2_PYTHIACPP_OMEGA_H
#define DECAY2_PYTHIACPP_OMEGA_H
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include "PDG_pythiaCPP_Omega.h"

using namespace std;

class Decay2_pythiaCPP_Omega
{
    map<string,bool> evtDecay;
    map<string,int> evtDecayNum;
    
public:
    Decay2_pythiaCPP_Omega();
    int Get_nDecays();
    bool Check_DecayPID(vector<int> DecayPID, string partName1, string partName2);
    void Init();
    void Set_All(vector<int> DecayPID);

    bool Get_All();
    bool Get_evtDecay(string evtName);
    int Get_evtDecayNum(string evtName);
    void Print_HistNumbers();
    
    string Get_trueDecay();
    int Get_trueDecayNum();
};

#endif

