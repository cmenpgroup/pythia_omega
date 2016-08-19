#ifndef CUTS_PYTHIACPP_OMEGA_H
#define CUTS_PYTHIACPP_OMEGA_H
#include <vector>
#include <string>

using namespace std;

class Cuts_pythiaCPP
{
    vector<string> CutsLabel;
    vector<double> RangeQSquared;
    vector<double> RangeWcut;
    
    int topo_nelec;
    int topo_npim;
    int topo_npip;
    int topo_ngam;
    
    double dalitzParentMass;
    double dalitzDaughter1Mass;
    double dalitzDaughter2Mass;
    double dalitzDaughter3Mass;
    
    bool cuts_omega_Q2;
    bool cuts_omega_W;
    bool cuts_omega_NumDetPart;
    bool cuts_omega_dalitz;
    bool cuts_omega_All;
    
public:
    Cuts_pythiaCPP();
    int Get_nCuts() {return CutsLabel.size();};
	string Get_CutsLabel(int num) {return CutsLabel[num];};
    double Get_QSquared_lo() {return RangeQSquared[0];};
    double Get_QSquared_hi() {return RangeQSquared[1];};
    double Get_Wcut_lo() {return RangeWcut[0];};
    double Get_Wcut_hi() {return RangeWcut[1];};
    int Get_Topo_nElec() {return topo_nelec;};
    int Get_Topo_nPim() {return topo_npim;};
    int Get_Topo_nPip() {return topo_npip;};
    int Get_Topo_nGam() {return topo_ngam;};
    void Set_Topo_nElec(int nElec) {topo_nelec = nElec;};
    void Set_Topo_nPim(int nPim) {topo_npim = nPim;};
    void Set_Topo_nPip(int nPip) {topo_npip = nPip;};
    void Set_Topo_nGam(int nGam) {topo_ngam = nGam;};
    bool Check_QSquared(double Qsq);
    bool Check_Wcut(double W);
    bool Check_NumDetPart(int nElec, int nPim, int nPip, int nGam);
    void Print_Cuts();
    
    void InitCuts();
    void SetCut_QSquared(double Qsq);
    bool GetCut_QSquared() {return cuts_omega_Q2;};
    void SetCut_Wcut(double W);
    bool GetCut_Wcut() {return cuts_omega_W;};
    void SetCut_NumDetPart(int nElec, int nPim, int nPip, int nGam);
    bool GetCut_NumDetPart() {return cuts_omega_NumDetPart;};
    void SetCut_Dalitz(double Msq12, double Msq23);
    bool GetCut_Dalitz() {return cuts_omega_dalitz;};
    
    void SetDalitz_Daughter(double mass, int index);
    double GetDalitz_Daughter(int index);
    void SetDalitz_Parent(double mass) {dalitzParentMass = mass;};
    double GetDalitz_Parent() {return dalitzParentMass;};
    bool Check_Dalitz(double Msq12, double Msq23);
};
#endif
