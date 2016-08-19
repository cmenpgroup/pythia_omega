#ifndef DECAY3_PYTHIACPP_OMEGA_H
#define DECAY3_PYTHIACPP_OMEGA_H
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include "PDG_pythiaCPP_Omega.h"

using namespace std;

class Decay3_pythiaCPP_Omega
{
    int TriCycle[6][3] = {{0,1,2},{1,2,0},{2,0,1},{1,0,2},{0,2,1},{2,1,0}};
    
    map<string,bool> evtDecay;
    map<string,int> evtDecayNum;

public:
    Decay3_pythiaCPP_Omega();
    int Get_nDecays();
    bool Check_DecayPID(vector<int> DecayPID, string partName1, string partName2, string partName3);
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

