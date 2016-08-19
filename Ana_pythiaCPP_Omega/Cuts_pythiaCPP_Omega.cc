#include <vector>
#include <string>
#include "Cuts_pythiaCPP_Omega.h"
#include <iostream>
#include "math.h"

Cuts_pythiaCPP::Cuts_pythiaCPP()
{
    CutsLabel.push_back("NoCuts");
    CutsLabel.push_back("QSquared");
    CutsLabel.push_back("Wcut");
    CutsLabel.push_back("NumDetPart");
    CutsLabel.push_back("Dalitz");
    
    topo_nelec = 1;
    topo_npim = 1;
    topo_npip = 1;
    topo_ngam = 2;
    
    RangeQSquared.push_back(1.0); // Lower limit on Q^2 (in Gev^2)
    RangeQSquared.push_back(100000.0); // Upper limit on Q^2 (in Gev^2)

    RangeWcut.push_back(2.0); // Lower limit on W (in Gev)
    RangeWcut.push_back(100000.0); // Upper limit on W (in Gev)
    
    this->InitCuts();
}

// initialize the cuts
void Cuts_pythiaCPP::InitCuts()
{
    cuts_omega_Q2 = false;
    cuts_omega_W = false;
    cuts_omega_NumDetPart = false;
    cuts_omega_All = false;
    cuts_omega_dalitz = false;
}

// check the cut on Q^2
bool Cuts_pythiaCPP::Check_QSquared(double Qsq)
{
	bool ret = (Qsq >= this->Get_QSquared_lo() && Qsq < this->Get_QSquared_hi()) ? true : false;
	
	return ret;
}

// set the value of the Q2 cut
void Cuts_pythiaCPP::SetCut_QSquared(double Qsq)
{
    cuts_omega_Q2 = this->Check_QSquared(Qsq);
}

// check the cut on W
bool Cuts_pythiaCPP::Check_Wcut(double W)
{
    bool ret = (W >= this->Get_Wcut_lo() && W < this->Get_Wcut_hi()) ? true : false;
    
    return ret;
}

// set the value of the pi0 mass cut
void Cuts_pythiaCPP::SetCut_Wcut(double W)
{
    cuts_omega_W = this->Check_Wcut(W);
}

// check the cut on particle topology
bool Cuts_pythiaCPP::Check_NumDetPart(int nElec, int nPim, int nPip, int nGam)
{
    bool check_elec = (nElec == this->Get_Topo_nElec());
    bool check_pim = (nPim == this->Get_Topo_nPim());
    bool check_pip = (nPip == this->Get_Topo_nPip());
    bool check_gam = (nGam == this->Get_Topo_nGam());
    
    bool ret = (check_elec && check_pim && check_pip && check_gam) ? true : false;
    
    return ret;
}

// set the value of the particle topology cut
void Cuts_pythiaCPP::SetCut_NumDetPart(int nElec, int nPim, int nPip, int nGam)
{
    cuts_omega_NumDetPart = this->Check_NumDetPart(nElec,nPim,nPip,nGam);
}

void Cuts_pythiaCPP::SetDalitz_Daughter(double mass, int index)
{
    switch(index){
        case 1: this->dalitzDaughter1Mass = mass; break;
        case 2: this->dalitzDaughter2Mass = mass; break;
        case 3: this->dalitzDaughter3Mass = mass; break;
        default:
            cout<<"Cuts_pythiaCPP::SetDalitz_Daughter - Incorrect index "<<index<<endl;
            exit(0);
            break;
    }
}

double Cuts_pythiaCPP::GetDalitz_Daughter(int index)
{
    double ret;
    switch(index){
        case 1: ret = this->dalitzDaughter1Mass; break;
        case 2: ret = this->dalitzDaughter2Mass; break;
        case 3: ret = this->dalitzDaughter3Mass; break;
        default:
            cout<<"Cuts_pythiaCPP::GetDalitz_Daughter - Incorrect index "<<index<<endl;
            exit(0);
            break;
    }
    return ret;
}

// check the cut on the Dalitz plot
bool Cuts_pythiaCPP::Check_Dalitz(double Msq12, double Msq23)
{
    bool ret;
    
    double Mparent = this->GetDalitz_Parent();
    double Mdaughter1 = this->GetDalitz_Daughter(1);
    double Mdaughter2 = this->GetDalitz_Daughter(2);
    double Mdaughter3 = this->GetDalitz_Daughter(3);
 
    double Msq12_Lo = (Mdaughter1 + Mdaughter2)*(Mdaughter1 + Mdaughter2);
    double Msq12_Hi = (Mparent - Mdaughter3)*(Mparent - Mdaughter3);
    double Msq23_Lo = (Mdaughter2 + Mdaughter3)*(Mdaughter2 + Mdaughter3);
    double Msq23_Hi = (Mparent - Mdaughter1)*(Mparent - Mdaughter1);
    
    bool Msq12_Limits = ((Msq12>= Msq12_Lo) && (Msq12< Msq12_Hi));
    bool Msq23_Limits = ((Msq23>= Msq23_Lo) && (Msq23< Msq23_Hi));
    
    if (Msq12_Limits && Msq23_Limits) {
        double Msq = Mparent*Mparent;
        double Msq1 = Mdaughter1*Mdaughter1;
        double Msq2 = Mdaughter2*Mdaughter2;
        double Msq3 = Mdaughter3*Mdaughter3;
    
        double E2 = (Msq12 - Msq1 + Msq2)/(2.0*sqrt(Msq12));
        double E3 = (Msq - Msq12 - Msq3)/(2.0*sqrt(Msq12));
    
        double Esum = E2 + E3;
    
        double EdiffMin = sqrt(E2*E2 - Msq2) + sqrt(E3*E3 - Msq3);
        double EdiffMax = sqrt(E2*E2 - Msq2) - sqrt(E3*E3 - Msq3);
    
        bool min = (Msq23 >= (Esum*Esum -  EdiffMin*EdiffMin)) ? true : false;
        bool max = (Msq23 < (Esum*Esum -  EdiffMax*EdiffMax)) ? true : false;
        ret = (min && max);
    }else{
        ret = false;
    }
    
    return ret;
}

// set the value of the Dalitz cut
void Cuts_pythiaCPP::SetCut_Dalitz(double Msq12, double Msq23)
{
    cuts_omega_dalitz = this->Check_Dalitz(Msq12,Msq23);
}

// print the cut information
void Cuts_pythiaCPP::Print_Cuts()
{
	int ii;
    cout<<"Cut Info (pythiaCPP Omega)"<<endl;
    cout<<"========================="<<endl;
    
    for(ii=0;ii<this->Get_nCuts();ii++){
        cout << this->Get_CutsLabel(ii) << "\t";
        if (this->Get_CutsLabel(ii).compare("QSquared")==0) {
            cout << "[" << this->Get_QSquared_lo() << "," << this->Get_QSquared_hi() << "] (GeV^2)" << endl;
        }else if (this->Get_CutsLabel(ii).compare("Wcut")==0) {
            cout << "[" << this->Get_Wcut_lo() << "," << this->Get_Wcut_hi() << "] (GeV)" << endl;
        }else{
            cout << endl;
        }
    }
    cout << endl;
}
