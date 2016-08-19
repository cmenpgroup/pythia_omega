#include "PDG_pythiaCPP_Omega.h"

PDG_pythiaCPP_Omega::PDG_pythiaCPP_Omega(){
    PDGid.push_back(1);     // down
    PDGid.push_back(-1);    // anti-down
    PDGid.push_back(2);     // up
    PDGid.push_back(-2);    // anti-up
    PDGid.push_back(3);     // strange
    PDGid.push_back(-3);    // anti-strange
    PDGid.push_back(11);    // electron
    PDGid.push_back(-11);   // anti-electron
    PDGid.push_back(21);    // gluon
    PDGid.push_back(22);    // gamma
    PDGid.push_back(111);   // pi0
    PDGid.push_back(211);   // pi+
    PDGid.push_back(-211);  // pi-
    PDGid.push_back(221);   // eta
    PDGid.push_back(113);   // rho0
    PDGid.push_back(213);   // rho+
    PDGid.push_back(-213);  // rho-
    PDGid.push_back(223);   // omega
    PDGid.push_back(130);   // K0-long
    PDGid.push_back(310);   // K0-short
    PDGid.push_back(311);   // K0
    PDGid.push_back(-311);  // anti-K0
    PDGid.push_back(313);   // K*(892)0
    PDGid.push_back(-313);  // anti-K*(892)0
    PDGid.push_back(323);   // K*(892)+
    PDGid.push_back(-323);  // K*(892)-
    PDGid.push_back(321);   // K+
    PDGid.push_back(-321);  // K-
    PDGid.push_back(331);   // eta'(958)
    PDGid.push_back(333);   // phi(1020)
    PDGid.push_back(2101);  // (ud)0
    PDGid.push_back(2103);  // (ud)1
    PDGid.push_back(2203);  // (uu)1
    PDGid.push_back(2112);  // neutron
    PDGid.push_back(2212);  // proton
    PDGid.push_back(1114);  // Delta-
    PDGid.push_back(2114);  // Delta0
    PDGid.push_back(2214);  // Delta+
    PDGid.push_back(2224);  // Delta++
    PDGid.push_back(3114);  // Sigma*-
    PDGid.push_back(3122);  // Lambda
    PDGid.push_back(3112);  // Sigma-
    PDGid.push_back(3212);  // Sigma0
    PDGid.push_back(3214);  // Sigma*0
    PDGid.push_back(3222);  // Sigma+
    PDGid.push_back(3224);  // Sigma(1385)

    PDGname.push_back("down");     // down
    PDGname.push_back("anti-down");    // anti-down
    PDGname.push_back("up");     // up
    PDGname.push_back("anti-up");    // anti-up
    PDGname.push_back("strange");     // strange
    PDGname.push_back("anti-strange");    // anti-strange
    PDGname.push_back("e-");    // electron
    PDGname.push_back("e+");   // anti-electron
    PDGname.push_back("gluon");    // gluon
    PDGname.push_back("gamma");    // gamma
    PDGname.push_back("pi0");   // pi0
    PDGname.push_back("pi+");   // pi+
    PDGname.push_back("pi-");  // pi-
    PDGname.push_back("eta");   // eta
    PDGname.push_back("rho0");   // rho0
    PDGname.push_back("rho+");   // rho+
    PDGname.push_back("rho-");  // rho-
    PDGname.push_back("omega");  // omega
    PDGname.push_back("KL0");   // K0-long
    PDGname.push_back("KS0");   // K0-short
    PDGname.push_back("K0");   // K0
    PDGname.push_back("anti-K0");  // anti-K0
    PDGname.push_back("K*(892)0");   // K*(892)0
    PDGname.push_back("anti-K*(892)0");  // anti-K*(892)0
    PDGname.push_back("K*(892)+");   // K*(892)+
    PDGname.push_back("K*(892)-");  // K*(892)-
    PDGname.push_back("K+");   // K+
    PDGname.push_back("K-");  // K-
    PDGname.push_back("eta'(958)");   // eta'(958)
    PDGname.push_back("phi(1020)");   // phi(1020)
    PDGname.push_back("(ud)0");  // (ud)0
    PDGname.push_back("(ud)1");  // (ud)1
    PDGname.push_back("(uu)1");  // (uu)1
    PDGname.push_back("n");  // neutron
    PDGname.push_back("p");  // proton
    PDGname.push_back("Delta-");  // Delta-
    PDGname.push_back("Delta0");  // Delta0
    PDGname.push_back("Delta+");  // Delta+
    PDGname.push_back("Delta++");  // Delta++
    PDGname.push_back("Sigma*-");  // Sigma*-
    PDGname.push_back("Lambda");  // Lambda
    PDGname.push_back("Sigma-");  // Sigma-
    PDGname.push_back("Sigma0");  // Sigma0
    PDGname.push_back("Sigma*0");  // Sigma*0
    PDGname.push_back("Sigma+");  // Sigma+
    PDGname.push_back("Sigma(1385)");  // Sigma(1385)
    
    PDGmass.push_back(0.0048);     // down
    PDGmass.push_back(0.0048);    // anti-down
    PDGmass.push_back(0.0023);     // up
    PDGmass.push_back(0.0023);    // anti-up
    PDGmass.push_back(0.095);     // strange
    PDGmass.push_back(0.095);    // anti-strange
    PDGmass.push_back(0.000511);    // electron
    PDGmass.push_back(0.000511);   // anti-electron
    PDGmass.push_back(0.0);    // gluon
    PDGmass.push_back(0.0);    // gamma
    PDGmass.push_back(0.135);   // pi0
    PDGmass.push_back(0.138);   // pi+
    PDGmass.push_back(0.138);  // pi-
    PDGmass.push_back(0.547862);   // eta
    PDGmass.push_back(0.77562);   // rho0
    PDGmass.push_back(0.77562);   // rho+
    PDGmass.push_back(0.77562);  // rho-
    PDGmass.push_back(0.78265);  // omega
    PDGmass.push_back(0.497611);   // K0-long
    PDGmass.push_back(0.497611);   // K0-short
    PDGmass.push_back(0.497611);   // K0
    PDGmass.push_back(0.497611);  // anti-K0
    PDGmass.push_back(0.892);   // K*(892)0
    PDGmass.push_back(0.892);  // anti-K*(892)0
    PDGmass.push_back(0.892);   // K*(892)+
    PDGmass.push_back(0.892);  // K*(892)-
    PDGmass.push_back(0.493677);   // K+
    PDGmass.push_back(0.493677);  // K-
    PDGmass.push_back(0.958);   // eta'(958)
    PDGmass.push_back(1.020);   // phi(1020)
    PDGmass.push_back(0.01);  // (ud)0
    PDGmass.push_back(0.01);  // (ud)1
    PDGmass.push_back(0.01);  // (uu)1
    PDGmass.push_back(0.935);  // neutron
    PDGmass.push_back(0.938);  // proton
    PDGmass.push_back(1.232);  // Delta-
    PDGmass.push_back(1.232);  // Delta0
    PDGmass.push_back(1.232);  // Delta+
    PDGmass.push_back(1.232);  // Delta++
    PDGmass.push_back(-1.0);  // Sigma*-
    PDGmass.push_back(1.115683);  // Lambda
    PDGmass.push_back(1.197449);  // Sigma-
    PDGmass.push_back(1.192642);  // Sigma0
    PDGmass.push_back(-1.0);  // Sigma*0
    PDGmass.push_back(1.18937);  // Sigma+
    PDGmass.push_back(1.385);  // Sigma(1385)
    
    if(this->Check_nPDG()==false){
        cout<<"PDG_pythiaCPP_Omega::Init - Mismatch in the PDG vectors."<<endl;
        exit(0);
    }
}

bool PDG_pythiaCPP_Omega::Check_nPDG(){
    bool CheckID = (this->Get_nPDGname()==this->Get_nPDGid()) ? true : false;
    bool CheckMass = (this->Get_nPDGname()==this->Get_nPDGmass()) ? true : false;
    return (CheckID && CheckMass);
}

int PDG_pythiaCPP_Omega::Get_name2id(string PDGstring){
    int ii;
    int ret = 0;
    
    for(ii=0;ii<this->Get_nPDGname();ii++){
        if (this->Get_PDGname(ii).compare(PDGstring)==0){
            ret = this->Get_PDGid(ii);
            break;
        }
    }
    return ret;
}

double PDG_pythiaCPP_Omega::Get_name2mass(string PDGstring){
    int ii;
    double ret = 0.0;
    
    for(ii=0;ii<this->Get_nPDGname();ii++){
        if (this->Get_PDGname(ii).compare(PDGstring)==0){
            ret = this->Get_PDGmass(ii);
            break;
        }
    }
    return ret;
}

string PDG_pythiaCPP_Omega::Get_id2name(int PDGcode){
    string ret;
    
    switch (PDGcode) {
        case 1: ret = "down"; break;
        case -1: ret = "anti-down"; break;
        case 2: ret = "up"; break;
        case -2: ret = "anti-up"; break;
        case 3: ret = "strange"; break;
        case -3: ret = "anti-strange"; break;
        case 11: ret = "e-"; break;
        case -11: ret = "e+"; break;
        case 21: ret = "gluon"; break;
        case 22: ret = "gamma"; break;
        case 111: ret = "pi0"; break;
        case 211: ret = "pi+"; break;
        case -211: ret = "pi-"; break;
        case 221: ret = "eta"; break;
        case 113: ret = "rho0"; break;
        case 213: ret = "rho+"; break;
        case -213: ret = "rho-"; break;
        case 223: ret = "omega"; break;
        case 130: ret = "KL0"; break;
        case 310: ret = "KS0"; break;
        case 311: ret = "K0"; break;
        case -311: ret = "anti-K0"; break;
        case -313: ret = "anti-K*(892)0"; break;
        case 313: ret = "K*(892)0"; break;
        case 323: ret = "K*(892)+"; break;
        case -323: ret = "K*(892)-"; break;
        case 321: ret = "K+"; break;
        case -321: ret = "K-"; break;
        case 331: ret = "eta'(958)"; break;
        case 333: ret = "phi(1020)"; break;
        case 2101: ret = "(ud)0"; break;
        case 2103: ret = "(ud)1"; break;
        case 2203: ret = "(uu)1"; break;
        case 2112: ret = "n"; break;
        case 2212: ret = "p"; break;
        case 1114: ret = "Delta-"; break;
        case 2114: ret = "Delta0"; break;
        case 2214: ret = "Delta+"; break;
        case 2224: ret = "Delta++"; break;
        case 3114: ret = "Sigma*-"; break;
        case 3122: ret = "Lambda"; break;
        case 3112: ret = "Sigma-"; break;
        case 3212: ret = "Sigma0"; break;
        case 3214: ret = "Sigma*0"; break;
        case 3222: ret = "Sigma+"; break;
        case 3224: ret = "Sigma(1385)"; break;
        default: ret = "unknown"; break;
    }
    return ret;
}

