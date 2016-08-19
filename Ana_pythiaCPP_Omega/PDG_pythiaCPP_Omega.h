#ifndef PDG_PYTHIACPP_OMEGA_H
#define PDG_PYTHIACPP_OMEGA_H
#include <vector>
#include <string>
#include <iostream>

using namespace std;

class PDG_pythiaCPP_Omega
{
    vector<int> PDGid;
    vector<string> PDGname;
    vector<double> PDGmass;
    
public:
    PDG_pythiaCPP_Omega();
    bool Check_nPDG();
    int Get_nPDGid() {return PDGid.size(); };
    int Get_nPDGname() {return PDGname.size(); };
    int Get_nPDGmass() {return PDGmass.size(); };
    
    string Get_PDGname(int num) {return PDGname[num]; };
    int Get_PDGid(int num) {return PDGid[num]; };
    double Get_PDGmass(int num) {return PDGmass[num]; };
    
    string Get_id2name(int PDGcode);
    int Get_name2id(string PDGstring);
    double Get_name2mass(string PDGstring);
    
};

#endif

