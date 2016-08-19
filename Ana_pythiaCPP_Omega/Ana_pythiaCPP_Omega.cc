#include "Ana_pythiaCPP_Omega.h"

void PrintUsage(char *processName)
{
    cerr << processName << " <options> <filename>\n";
    cerr << "\toptions are:\n";
    cerr << "\t-o<filename>\tROOT output file (def. = ctProcess_omega.root).\n";
    cerr << "\t-M#\t\tprocess maximum # of events.\n";
    cerr << "\t-D#\t\tinform user when # of events have been processed (def. = 1000).\n";
    cerr << "\t-A#\t\tAnalysis type (def. = 0)\n";
    cerr << "\t-h\t\tprint the above" << endl;
}

void PrintAnalysisTime(float tStart, float tStop){
    //time to complete function
    float minutes = 0;
    float seconds = 0;
    minutes = (tStop - tStart)/1000000;
    minutes = (minutes)/60;
    seconds = fmod(minutes,1);
    minutes = minutes-seconds;
    seconds = seconds*60;
    
    if (minutes==0){
        cout<<endl<<"Completed in "<<seconds<<" seconds."<<endl<<endl;
    }
    else{
        cout<<endl<<"Completed in "<<minutes<<" minutes and "<<seconds<<" seconds."<<endl<<endl;
    }
}

int main (int argc, char **argv) {
    
    extern char *optarg;
    int c;
    extern int optind;
    
    int i;
    
    int dEvents = 1000; // increment of events for processing print statement
    int MaxEvents = 0; // max. number of events to process
    int iAna = 0; // analysis type
    int totOmegas = 0;
    
    string inFile;
    string outFile = "Ana_pythiaCPP_Omega.root";
    
    float timeStart = clock(); // start time
    
    for (i = 0; i < argc; ++i) cerr << argv[i] << " "; cerr << endl;
    while ((c = getopt(argc,argv, "o:M:D:T:ih")) != -1 ) {
        switch (c) {
            case 'o': outFile = optarg; break;
            case 'M': MaxEvents = atoi(optarg); break;
            case 'D': dEvents = atoi(optarg); break;
            case 'A': iAna = atoi(optarg); break;
            case 'h':
                PrintUsage(argv[0]);
                exit(0);
                break;
                
            default:
                cerr << "Unrecognized argument: " << optarg << endl;
                PrintUsage(argv[0]);
                exit(0);
                break;
        }
    }
    
    myHistManager.BookHist();
    
    for (i = optind; i < argc; ++i) {
        inFile = argv[i]; // process all arguments on command line.
        if (inFile != '-') { // we have a file to process
            cout << "Analyzing file " << inFile << endl; // let user know which file is being processed
            // process the root file and return number of processed events
            totOmegas = process(inFile,iAna,MaxEvents,dEvents);
            cout<<"=========================="<<endl;
            cout<<totOmegas<<" omegas found"<<endl; // print out stats
            cout<<"=========================="<<endl;
        }
    }
    
    myHistManager.WriteHist(outFile); // write the histograms to file
    
    float timeStop = clock();
    PrintAnalysisTime(timeStart,timeStop);
}

int process (string inFile, int iAna, int MaxEvents, int dEvents) {
    
    int i, j, k;
    int ii, jj;
    int iPC; // particle combination index
    int iPhoton;
    
    int totOmegas = 0;
    int nPip, nPim, nPi0, nPhoton, nOmega, nEta, nEtaPrime, nPhi; // particle counters per event
    int nUnknown = 0;
    
    int unknownParent; // index for the unknown particle
    
    TLorentzVector beam;
    TLorentzVector target;
    TLorentzVector electron;
    TLorentzVector recoil;
    TLorentzVector BeamMinusElectron;
    TLorentzVector missing;
    TLorentzVector W_TLV;
    TLorentzVector pip;
    TLorentzVector pim;
    TLorentzVector pi0;
    TLorentzVector photon[2];
    TLorentzVector twoPhotons;
    TLorentzVector omega;
    TLorentzVector PipPim;
    TLorentzVector PipPimPi0;
    TLorentzVector X;
    TLorentzVector Xrecoil;
    TLorentzVector Xpip;
    TLorentzVector Xpim;
    TLorentzVector Xpi0;
    
    Float_t Qsq, W;
    
    int ntNpart, ntTargetA, ntTargetZ, ntProcess;
    int ntks[MAX_TRACKS], ntPID[MAX_TRACKS], ntParent[MAX_TRACKS];
    double ntEbeam, ntNu;
    double ntPx[MAX_TRACKS], ntPy[MAX_TRACKS], ntPz[MAX_TRACKS], ntP[MAX_TRACKS], ntE[MAX_TRACKS];
    
    double ChargedPionAngle;
    double ThreePionAngle;
    double TwoPhotonAngle;
    
    vector<int> tempPid;
    vector<int> tempParent;
    vector<int> tempKs;
    vector<int> unknownDecayList;
    
    Cuts_pythiaCPP myCuts;
    PDG_pythiaCPP_Omega myPDG; // declare the PDG object
    Decay1_pythiaCPP_Omega myDecay1; // declare the object to handle the one body decay
    Decay2_pythiaCPP_Omega myDecay2; // declare the object to handle the two body decay
    Decay3_pythiaCPP_Omega myDecay3; // declare the object to handle the three body decay
    
    bool cuts_Qsq;
    bool cuts_W;
    bool cuts_beam;
    bool cuts_target;
    bool cuts_electron;
    bool cuts_recoil;
    bool evtAll;
    bool evtDecay1;
    bool evtDecay2;
    bool evtDecay3;
    bool evtDecay4;
    bool evtOneBodyDecay;
    bool evtTwoBodyDecay;
    bool evtThreeBodyDecay;
    bool evtFourBodyDecay;
    bool no_other_mesons;
    bool unknownCheck;
    
    // open text files for writing events with certain particle combinations
    char OutFile[100];
    sprintf(OutFile,"PartList09.dat");
    std::ofstream fout09(OutFile);
    
    sprintf(OutFile,"PartList10.dat");
    std::ofstream fout10(OutFile);
    
    sprintf(OutFile,"PartList11.dat");
    std::ofstream fout11(OutFile);
    
    sprintf(OutFile,"PartList12.dat");
    std::ofstream fout12(OutFile);
    
    TFile *fm = new TFile(inFile.c_str(),"READ"); // data file containing the trees
    
    TTree *myTree = (TTree*)fm->Get("MC");
    
    int nEntries = (int)myTree->GetEntries();
    
    cout<<"Number of entries = "<<nEntries<<endl;
    
    if(MaxEvents==0) MaxEvents = nEntries;
    
    myTree->SetBranchAddress("Ntracks",&ntNpart);
    myTree->SetBranchAddress("Eb",&ntEbeam);
    myTree->SetBranchAddress("tarA",&ntTargetA);
    myTree->SetBranchAddress("tarZ",&ntTargetZ);
    myTree->SetBranchAddress("process",&ntProcess);
    myTree->SetBranchAddress("nu",&ntNu);
    
    myTree->SetBranchAddress("ks",&ntks);
    myTree->SetBranchAddress("type",&ntPID);
    myTree->SetBranchAddress("parent",&ntParent);
    myTree->SetBranchAddress("px",&ntPx);
    myTree->SetBranchAddress("py",&ntPy);
    myTree->SetBranchAddress("pz",&ntPz);
    myTree->SetBranchAddress("p",&ntP);
    myTree->SetBranchAddress("E",&ntE);
    
    for(ii=0; ii<MaxEvents; ii++){

        myTree->GetEntry(ii); // retrieve the event from the ntuple

        if(!(ii % dEvents)) cout<<ii<<endl;
        
        cuts_beam = false;
        cuts_target = false;
        cuts_electron = false;
        cuts_recoil = false;
        cuts_Qsq = false;
        cuts_W =false;
        
        evtAll=false;
        evtDecay1=false;
        evtDecay2=false;
        evtDecay3=false;
        evtDecay4=false;
        evtOneBodyDecay=false;
        evtTwoBodyDecay=false;
        evtThreeBodyDecay=false;
        evtFourBodyDecay=false;
        no_other_mesons=false;
        
        unknownCheck=false;
        
        iPhoton = 0;
        nPip = 0;
        nPim = 0;
        nPi0 = 0;
        nPhoton = 0;
        nOmega = 0;
        nEta = 0;
        nEtaPrime = 0;
        nPhi = 0;
        
        unknownParent = -1; // init. unknown particle index
        unknownDecayList.clear();
        
        tempKs.clear();
        tempPid.clear();
        tempParent.clear();
        
        for (i=0; i<ntNpart; i++) {
            tempKs.push_back(ntks[i]);
            tempPid.push_back(ntPID[i]);
            tempParent.push_back(ntParent[i]);
            
            if(myPDG.Get_id2name(ntPID[i]).compare("e-")==0 && ntParent[i]==0){

                beam.SetPxPyPzE(ntPx[i],ntPy[i],ntPz[i],ntE[i]);
                cuts_beam = true;
            }
            
            if(myPDG.Get_id2name(ntPID[i]).compare("p")==0 && ntParent[i]==0){
                target.SetPxPyPzE(ntPx[i],ntPy[i],ntPz[i],ntE[i]);
                cuts_target =true;
            }
            
            if(myPDG.Get_id2name(ntPID[i]).compare("p")==0 && ntks[i]==1){
                recoil.SetPxPyPzE(ntPx[i],ntPy[i],ntPz[i],ntE[i]);
                cuts_recoil = true;
            }
            
            if(myPDG.Get_id2name(ntPID[i]).compare("e-")==0 && ntks[i]==1){
                electron.SetPxPyPzE(ntPx[i],ntPy[i],ntPz[i],ntE[i]);
                if(Cut_CLAS6_Theta_Ana(electron, iAna, "EC")) cuts_electron = true;
            }
            
            if(myPDG.Get_id2name(ntPID[i]).compare("pi+")==0 && ntks[i]==1){
                pip.SetPxPyPzE(ntPx[i],ntPy[i],ntPz[i],ntE[i]);
                if(Cut_CLAS6_Theta_Ana(pip, iAna, "TOF")) nPip++;
            }
            if(myPDG.Get_id2name(ntPID[i]).compare("pi-")==0 && ntks[i]==1){
                pim.SetPxPyPzE(ntPx[i],ntPy[i],ntPz[i],ntE[i]);
                if(Cut_CLAS6_Theta_Ana(pim, iAna, "TOF")) nPim++;
            }
            if(myPDG.Get_id2name(ntPID[i]).compare("pi0")==0 && ntks[i]==11){
                pi0.SetPxPyPzE(ntPx[i],ntPy[i],ntPz[i],ntE[i]);
                if(Cut_CLAS6_Theta_Ana(pi0, iAna, "TOF")) nPi0++;
            }
            if(myPDG.Get_id2name(ntPID[i]).compare("gamma")==0 && ntks[i]==1){
                if(iPhoton<2){
                    photon[iPhoton].SetPxPyPzE(ntPx[i],ntPy[i],ntPz[i],ntE[i]);
                    if(Cut_CLAS6_Theta_Ana(photon[iPhoton], iAna, "EC")) nPhoton++;
                }
                iPhoton++;
            }
            if(myPDG.Get_id2name(ntPID[i]).compare("omega")==0 && ntks[i]==11){
                totOmegas++;
                nOmega++;
                omega.SetPxPyPzE(ntPx[i],ntPy[i],ntPz[i],ntE[i]);
            }

            if(myPDG.Get_id2name(ntPID[i]).compare("eta")==0 && ntks[i]==11){
                nEta++;
            }

            if(myPDG.Get_id2name(ntPID[i]).compare("phi(1020)")==0 && ntks[i]==11){
                nPhi++;
            }

            if(myPDG.Get_id2name(ntPID[i]).compare("eta'(958)")==0 && ntks[i]==11){
                nEtaPrime++;
            }

            if(myPDG.Get_id2name(ntPID[i]).compare("unknown")==0 && (ntks[i]==11 || ntks[i]==21)){
                nUnknown++;
                unknownParent = i+1; // parent id of the unknown particle
//                unknownDecayList.clear();
            }
            
            if(ntParent[i]==unknownParent){
                unknownDecayList.push_back(ntPID[i]);
            }

            if(ntks[i]==11) myHistManager.Get_hIntermedPart()->Fill(ntPID[i]);
        }

/*        if(unknownDecayList.size()==0){
            cout<<endl<<"Event: "<<ii<<"\tProcess: "<<ntProcess<<endl;
            for(jj=0; jj<ntNpart; jj++){
                cout<<jj+1<<"\t"<<tempKs[jj]<<"\t"<<tempPid[jj]<<"\t"<<tempParent[jj]<<"\t"<<myPDG.Get_id2name(tempPid[jj])<<endl;
            }
        }
*/
        if(unknownDecayList.size()==1){
            myDecay1.Init();
            myDecay1.Set_All(unknownDecayList);
            evtOneBodyDecay = myDecay1.Get_All();
            if(!evtOneBodyDecay){
                cout<<myPDG.Get_id2name(unknownDecayList[0])<<endl;
            }
        }
        
        if(unknownDecayList.size()==2){
            myDecay2.Init();
            myDecay2.Set_All(unknownDecayList);
            evtTwoBodyDecay = myDecay2.Get_All();
            if(!evtTwoBodyDecay){
                cout<<myPDG.Get_id2name(unknownDecayList[0])<<"\t"<<myPDG.Get_id2name(unknownDecayList[1])<<endl;
            }
        }
        
        if(unknownDecayList.size()==3){
            myDecay3.Init();
            myDecay3.Set_All(unknownDecayList);
            evtThreeBodyDecay = myDecay3.Get_All();
            if(!evtThreeBodyDecay){
                cout<<myPDG.Get_id2name(unknownDecayList[0])<<"\t"<<myPDG.Get_id2name(unknownDecayList[1])<<"\t"<<myPDG.Get_id2name(unknownDecayList[2])<<endl;
            }
        }

        if(unknownDecayList.size()==4){
            evtFourBodyDecay = true;
        }
        
        if(nPip>=1 && nPim>=1 && nPi0>=1) evtDecay1 = true;
        if(nPip==1 && nPim==1 && nPi0==1) evtDecay2 = true;
        if(nPip>=1 && nPim>=1 && nPhoton>=2) evtDecay3 = true;
        if(nPip==1 && nPim==1 && nPhoton==2) evtDecay4 = true;
        
        if(nEta==0 && nEtaPrime==0 && nPhi==0) no_other_mesons = true; // set if these mesons are not in the event
        
        if(evtOneBodyDecay || evtTwoBodyDecay || evtThreeBodyDecay || evtFourBodyDecay) evtAll = true;
        
        if(cuts_beam && cuts_target && cuts_electron){

            BeamMinusElectron = beam - electron;
            
            Qsq = -1.0*BeamMinusElectron.Mag2();
            myHistManager.Get_hQsq_NoCuts()->Fill(Qsq);
        
            myHistManager.Get_hNu_NoCuts()->Fill(ntNu);
        
            W_TLV = BeamMinusElectron + target;
            W = W_TLV.M();
            myHistManager.Get_hW_NoCuts()->Fill(W);
            
            myCuts.SetCut_QSquared(Qsq);
            cuts_Qsq = myCuts.GetCut_QSquared();
            
            myCuts.SetCut_Wcut(W);
            cuts_W = myCuts.GetCut_Wcut();
        
            if(cuts_Qsq && cuts_W){

                myHistManager.Get_hQsq()->Fill(Qsq);
                myHistManager.Get_hNu()->Fill(ntNu);
                myHistManager.Get_hW()->Fill(W);
                
                myHistManager.Get_hPartPerEvt()->Fill(nPip,1);
                myHistManager.Get_hPartPerEvt()->Fill(nPim,2);
                myHistManager.Get_hPartPerEvt()->Fill(nPi0,3);
                myHistManager.Get_hPartPerEvt()->Fill(nPhoton,4);
                myHistManager.Get_hPartPerEvt()->Fill(nOmega,5);
        
                if(evtDecay1) myHistManager.Get_hPartPerEvt()->Fill(1,6);
                if(evtDecay2) myHistManager.Get_hPartPerEvt()->Fill(2,6);
                if(evtDecay3) myHistManager.Get_hPartPerEvt()->Fill(3,6);
                if(evtDecay4) myHistManager.Get_hPartPerEvt()->Fill(4,6);
                if(evtDecay1 && nOmega>=1) myHistManager.Get_hPartPerEvt()->Fill(5,6);
                if(evtDecay2 && nOmega==1) myHistManager.Get_hPartPerEvt()->Fill(6,6);
                if(evtDecay3 && nOmega>=1) myHistManager.Get_hPartPerEvt()->Fill(7,6);
                if(evtDecay4 && nOmega==1) myHistManager.Get_hPartPerEvt()->Fill(8,6);
                if(evtDecay1 && nOmega==0) myHistManager.Get_hPartPerEvt()->Fill(9,6);
                if(evtDecay2 && nOmega==0) myHistManager.Get_hPartPerEvt()->Fill(10,6);
                if(evtDecay3 && nOmega==0) myHistManager.Get_hPartPerEvt()->Fill(11,6);
                if(evtDecay4 && nOmega==0) myHistManager.Get_hPartPerEvt()->Fill(12,6);
                if(evtDecay1 && nOmega==0 && no_other_mesons) myHistManager.Get_hPartPerEvt()->Fill(13,6);
                if(evtDecay2 && nOmega==0 && no_other_mesons) myHistManager.Get_hPartPerEvt()->Fill(14,6);
                if(evtDecay3 && nOmega==0 && no_other_mesons) myHistManager.Get_hPartPerEvt()->Fill(15,6);
                if(evtDecay4 && nOmega==0 && no_other_mesons) myHistManager.Get_hPartPerEvt()->Fill(16,6);
            
                if(cuts_recoil){
                    X = BeamMinusElectron + target;
                    Xrecoil = X - recoil;
                    Xpip = X - pip;
                    Xpim = X - pim;
                    Xpi0 = X - pi0;
                }
                
                PipPim = pip + pim; // TLVector of the pi+ pi- pair
                if(evtDecay1 || evtDecay2) PipPimPi0 = PipPim + pi0; // TLVector of pi+ pi- pi0
                if(evtDecay3 || evtDecay4){
                    twoPhotons = photon[0] + photon[1]; // 2 photon recon. TLVector
                    PipPimPi0 = PipPim + twoPhotons; // TLVector of pi+ pi- pi0, with pi0 recon from 2 photons
                    TwoPhotonAngle = TMath::RadToDeg()*photon[0].Angle(photon[1].Vect());
                }
                
                
                if(evtDecay1){
                    FillHists_omega(pip,pim,pi0,ntNu,1);
                    if(cuts_recoil){
                        myHistManager.Get_hMMsq_X()->Fill(X.Mag2(),1);
                        myHistManager.Get_hMMsq_Xrecoil()->Fill(Xrecoil.Mag2(),1);
                        myHistManager.Get_hMMsq_Xpip()->Fill(Xpip.Mag2(),1);
                        myHistManager.Get_hMMsq_Xpim()->Fill(Xpim.Mag2(),1);
                        myHistManager.Get_hMMsq_Xpi0()->Fill(Xpi0.Mag2(),1);
                    }
                    if(nOmega){
                        FillHists_omega(pip,pim,pi0,ntNu,5);
                        if(cuts_recoil){
                            myHistManager.Get_hMMsq_X()->Fill(X.Mag2(),5);
                            myHistManager.Get_hMMsq_Xrecoil()->Fill(Xrecoil.Mag2(),5);
                            myHistManager.Get_hMMsq_Xpip()->Fill(Xpip.Mag2(),5);
                            myHistManager.Get_hMMsq_Xpim()->Fill(Xpim.Mag2(),5);
                            myHistManager.Get_hMMsq_Xpi0()->Fill(Xpi0.Mag2(),5);
                        }
                    }else{
                        FillHists_omega(pip,pim,pi0,ntNu,9);
                        fout09<<endl<<"Event: "<<ii<<"\tProcess: "<<ntProcess<<endl;
                        for(jj=0; jj<ntNpart; jj++){
                            fout09<<jj+1<<"\t"<<tempKs[jj]<<"\t"<<tempPid[jj]<<"\t"<<tempParent[jj]<<"\t"<<myPDG.Get_id2name(tempPid[jj])<<endl;
                        }
                        if(cuts_recoil){
                            myHistManager.Get_hMMsq_X()->Fill(X.Mag2(),9);
                            myHistManager.Get_hMMsq_Xrecoil()->Fill(Xrecoil.Mag2(),9);
                            myHistManager.Get_hMMsq_Xpip()->Fill(Xpip.Mag2(),9);
                            myHistManager.Get_hMMsq_Xpim()->Fill(Xpim.Mag2(),9);
                            myHistManager.Get_hMMsq_Xpi0()->Fill(Xpi0.Mag2(),9);
                        }
                        
                        if(no_other_mesons) FillHists_omega(pip,pim,pi0,ntNu,13);
                        if(!evtAll) FillHists_omega(pip,pim,pi0,ntNu,17);
                    }
                    if(evtOneBodyDecay) myHistManager.Get_hIM_NotOmega_Decay1(0)->Fill(PipPimPi0.M(),myDecay1.Get_trueDecayNum());
                    if(evtTwoBodyDecay) myHistManager.Get_hIM_NotOmega_Decay2(0)->Fill(PipPimPi0.M(),myDecay2.Get_trueDecayNum());
                    if(evtThreeBodyDecay) myHistManager.Get_hIM_NotOmega_Decay3(0)->Fill(PipPimPi0.M(),myDecay3.Get_trueDecayNum());
                }

                if(evtDecay2){
                    FillHists_omega(pip,pim,pi0,ntNu,2);
                    if(cuts_recoil){
                        myHistManager.Get_hMMsq_X()->Fill(X.Mag2(),2);
                        myHistManager.Get_hMMsq_Xrecoil()->Fill(Xrecoil.Mag2(),2);
                        myHistManager.Get_hMMsq_Xpip()->Fill(Xpip.Mag2(),2);
                        myHistManager.Get_hMMsq_Xpim()->Fill(Xpim.Mag2(),2);
                        myHistManager.Get_hMMsq_Xpi0()->Fill(Xpi0.Mag2(),2);
                    }
                
                    if(nOmega){
                        FillHists_omega(pip,pim,pi0,ntNu,6);
                        if(cuts_recoil){
                            myHistManager.Get_hMMsq_X()->Fill(X.Mag2(),6);
                            myHistManager.Get_hMMsq_Xrecoil()->Fill(Xrecoil.Mag2(),6);
                            myHistManager.Get_hMMsq_Xpip()->Fill(Xpip.Mag2(),6);
                            myHistManager.Get_hMMsq_Xpim()->Fill(Xpim.Mag2(),6);
                            myHistManager.Get_hMMsq_Xpi0()->Fill(Xpi0.Mag2(),6);
                        }
                    }else{
                        FillHists_omega(pip,pim,pi0,ntNu,10);
                        fout10<<endl<<"Event: "<<ii<<"\tProcess: "<<ntProcess<<endl;
                        for(jj=0; jj<ntNpart; jj++){
                            fout10<<jj+1<<"\t"<<tempKs[jj]<<"\t"<<tempPid[jj]<<"\t"<<tempParent[jj]<<"\t"<<myPDG.Get_id2name(tempPid[jj])<<endl;
                        }
                        if(cuts_recoil){
                            myHistManager.Get_hMMsq_X()->Fill(X.Mag2(),10);
                            myHistManager.Get_hMMsq_Xrecoil()->Fill(Xrecoil.Mag2(),10);
                            myHistManager.Get_hMMsq_Xpip()->Fill(Xpip.Mag2(),10);
                            myHistManager.Get_hMMsq_Xpim()->Fill(Xpim.Mag2(),10);
                            myHistManager.Get_hMMsq_Xpi0()->Fill(Xpi0.Mag2(),10);
                        }
                        
                        if(no_other_mesons) FillHists_omega(pip,pim,pi0,ntNu,14);
                        if(!evtAll) FillHists_omega(pip,pim,pi0,ntNu,18);
                    }
                    if(evtOneBodyDecay) myHistManager.Get_hIM_NotOmega_Decay1(1)->Fill(PipPimPi0.M(),myDecay1.Get_trueDecayNum());
                    if(evtTwoBodyDecay) myHistManager.Get_hIM_NotOmega_Decay2(1)->Fill(PipPimPi0.M(),myDecay2.Get_trueDecayNum());
                    if(evtThreeBodyDecay) myHistManager.Get_hIM_NotOmega_Decay3(1)->Fill(PipPimPi0.M(),myDecay3.Get_trueDecayNum());
                }
                
                if(evtDecay3){
                    FillHists_omega(pip,pim,twoPhotons,ntNu,3);
                    myHistManager.Get_hOpAng_TwoPhoton()->Fill(TwoPhotonAngle,0);
                    if(cuts_recoil){
                        myHistManager.Get_hMMsq_X()->Fill(X.Mag2(),3);
                        myHistManager.Get_hMMsq_Xrecoil()->Fill(Xrecoil.Mag2(),3);
                        myHistManager.Get_hMMsq_Xpip()->Fill(Xpip.Mag2(),3);
                        myHistManager.Get_hMMsq_Xpim()->Fill(Xpim.Mag2(),3);
                        myHistManager.Get_hMMsq_Xpi0()->Fill(Xpi0.Mag2(),3);
                    }
                    
                    if(nOmega){
                        FillHists_omega(pip,pim,twoPhotons,ntNu,7);
                        myHistManager.Get_hOpAng_TwoPhoton()->Fill(TwoPhotonAngle,2);
                        if(cuts_recoil){
                            myHistManager.Get_hMMsq_X()->Fill(X.Mag2(),7);
                            myHistManager.Get_hMMsq_Xrecoil()->Fill(Xrecoil.Mag2(),7);
                            myHistManager.Get_hMMsq_Xpip()->Fill(Xpip.Mag2(),7);
                            myHistManager.Get_hMMsq_Xpim()->Fill(Xpim.Mag2(),7);
                            myHistManager.Get_hMMsq_Xpi0()->Fill(Xpi0.Mag2(),7);
                        }
                    }else{
                        FillHists_omega(pip,pim,twoPhotons,ntNu,11);
                        myHistManager.Get_hOpAng_TwoPhoton()->Fill(TwoPhotonAngle,4);
                        fout11<<endl<<"Event: "<<ii<<"\tProcess: "<<ntProcess<<endl;
                        for(jj=0; jj<ntNpart; jj++){
                            fout11<<jj+1<<"\t"<<tempKs[jj]<<"\t"<<tempPid[jj]<<"\t"<<tempParent[jj]<<"\t"<<myPDG.Get_id2name(tempPid[jj])<<endl;
                        }
                        if(cuts_recoil){
                            myHistManager.Get_hMMsq_X()->Fill(X.Mag2(),11);
                            myHistManager.Get_hMMsq_Xrecoil()->Fill(Xrecoil.Mag2(),11);
                            myHistManager.Get_hMMsq_Xpip()->Fill(Xpip.Mag2(),11);
                            myHistManager.Get_hMMsq_Xpim()->Fill(Xpim.Mag2(),11);
                            myHistManager.Get_hMMsq_Xpi0()->Fill(Xpi0.Mag2(),11);
                        }
                        
                        if(no_other_mesons) FillHists_omega(pip,pim,twoPhotons,ntNu,15);
                        if(!evtAll) FillHists_omega(pip,pim,pi0,ntNu,19);
                    }
                    if(evtOneBodyDecay) myHistManager.Get_hIM_NotOmega_Decay1(2)->Fill(PipPimPi0.M(),myDecay1.Get_trueDecayNum());
                    if(evtTwoBodyDecay) myHistManager.Get_hIM_NotOmega_Decay2(2)->Fill(PipPimPi0.M(),myDecay2.Get_trueDecayNum());
                    if(evtThreeBodyDecay) myHistManager.Get_hIM_NotOmega_Decay3(2)->Fill(PipPimPi0.M(),myDecay3.Get_trueDecayNum());
                }
                
                if(evtDecay4){
                    FillHists_omega(pip,pim,twoPhotons,ntNu,4);
                    myHistManager.Get_hOpAng_TwoPhoton()->Fill(TwoPhotonAngle,1);
                    if(cuts_recoil){
                        myHistManager.Get_hMMsq_X()->Fill(X.Mag2(),4);
                        myHistManager.Get_hMMsq_Xrecoil()->Fill(Xrecoil.Mag2(),4);
                        myHistManager.Get_hMMsq_Xpip()->Fill(Xpip.Mag2(),4);
                        myHistManager.Get_hMMsq_Xpim()->Fill(Xpim.Mag2(),4);
                        myHistManager.Get_hMMsq_Xpi0()->Fill(Xpi0.Mag2(),4);
                    }
                    
                    if(nOmega){
                        FillHists_omega(pip,pim,twoPhotons,ntNu,8);
                        myHistManager.Get_hOpAng_TwoPhoton()->Fill(TwoPhotonAngle,3);
                        if(cuts_recoil){
                            myHistManager.Get_hMMsq_X()->Fill(X.Mag2(),8);
                            myHistManager.Get_hMMsq_Xrecoil()->Fill(Xrecoil.Mag2(),8);
                            myHistManager.Get_hMMsq_Xpip()->Fill(Xpip.Mag2(),8);
                            myHistManager.Get_hMMsq_Xpim()->Fill(Xpim.Mag2(),8);
                            myHistManager.Get_hMMsq_Xpi0()->Fill(Xpi0.Mag2(),8);
                        }
                    }else{
                        FillHists_omega(pip,pim,twoPhotons,ntNu,12);
                        myHistManager.Get_hOpAng_TwoPhoton()->Fill(TwoPhotonAngle,5);
                        fout12<<endl<<"Event: "<<ii<<"\tProcess: "<<ntProcess<<endl;
                        for(jj=0; jj<ntNpart; jj++){
                            fout12<<jj+1<<"\t"<<tempKs[jj]<<"\t"<<tempPid[jj]<<"\t"<<tempParent[jj]<<"\t"<<myPDG.Get_id2name(tempPid[jj])<<endl;
                        }
                        if(cuts_recoil){
                            myHistManager.Get_hMMsq_X()->Fill(X.Mag2(),12);
                            myHistManager.Get_hMMsq_Xrecoil()->Fill(Xrecoil.Mag2(),12);
                            myHistManager.Get_hMMsq_Xpip()->Fill(Xpip.Mag2(),12);
                            myHistManager.Get_hMMsq_Xpim()->Fill(Xpim.Mag2(),12);
                            myHistManager.Get_hMMsq_Xpi0()->Fill(Xpi0.Mag2(),12);
                        }
                        
                        if(no_other_mesons)FillHists_omega(pip,pim,twoPhotons,ntNu,16);
                        if(!evtAll) FillHists_omega(pip,pim,pi0,ntNu,20);
                        if(!evtAll){
                            cout<<endl<<"Event: "<<ii<<"\tProcess: "<<ntProcess<<endl;
                            for(jj=0; jj<ntNpart; jj++){
                                cout<<jj+1<<"\t"<<tempKs[jj]<<"\t"<<tempPid[jj]<<"\t"<<tempParent[jj]<<"\t"<<myPDG.Get_id2name(tempPid[jj])<<endl;
                            }
                        }
                    }
                    if(evtOneBodyDecay) myHistManager.Get_hIM_NotOmega_Decay1(3)->Fill(PipPimPi0.M(),myDecay1.Get_trueDecayNum());
                    if(evtTwoBodyDecay) myHistManager.Get_hIM_NotOmega_Decay2(3)->Fill(PipPimPi0.M(),myDecay2.Get_trueDecayNum());
                    if(evtThreeBodyDecay) myHistManager.Get_hIM_NotOmega_Decay3(3)->Fill(PipPimPi0.M(),myDecay3.Get_trueDecayNum());
                }
            }
        }
    }
    
    fout09.close();
    fout10.close();
    fout11.close();
    fout12.close();
    
    cout<<"Unknowns "<<nUnknown<<endl;
    return totOmegas;
}

void FillHists_omega(TLorentzVector pip, TLorentzVector pim, TLorentzVector pi0, double nu, int iPC){
    TLorentzVector PipPim = pip + pim;
    TLorentzVector PipPi0 = pip + pi0;    
    TLorentzVector PipPimPi0 = pip + pim + pi0;
    
    double ChargedPionAngle = TMath::RadToDeg()*pip.Angle(pim.Vect());
    double ThreePionAngle = TMath::RadToDeg()*PipPim.Angle(pi0.Vect());
    
    myHistManager.Get_hIMomega()->Fill(PipPimPi0.M(),iPC);
    myHistManager.Get_hIMomega_VS_IMPipPim(iPC-1)->Fill(PipPim.M(),PipPimPi0.M());
    myHistManager.Get_hOpAng_PipPim()->Fill(ChargedPionAngle,iPC);
    myHistManager.Get_hOpAng_PipPimPi0()->Fill(ThreePionAngle,iPC);
    myHistManager.Get_hz_fracEnergy()->Fill(PipPimPi0.E()/nu,iPC);
    myHistManager.Get_hDalitz_pip(iPC-1)->Fill(PipPi0.M2(),PipPim.M2());
    
    PDG_pythiaCPP_Omega myPDG; // create a PDG object to extract the particle masses
    Cuts_pythiaCPP myCuts; // creat the Cuts object to test the Dalitz cut
    myCuts.SetDalitz_Parent(myPDG.Get_name2mass("omega"));
    myCuts.SetDalitz_Daughter(myPDG.Get_name2mass("pi0"),1);
    myCuts.SetDalitz_Daughter(myPDG.Get_name2mass("pi+"),2);
    myCuts.SetDalitz_Daughter(myPDG.Get_name2mass("pi-"),3);
    myCuts.SetCut_Dalitz(PipPi0.M2(),PipPim.M2());
    
    if(myCuts.GetCut_Dalitz()){
        myHistManager.Get_hIMomega_CutDalitz()->Fill(PipPimPi0.M(),iPC);
        myHistManager.Get_hDalitz_pip_CutDalitz(iPC-1)->Fill(PipPi0.M2(),PipPim.M2());
    }
}

// Check that polar angle is within the geometrical acceptance of the EC.  Input angle in degrees.
bool Cut_CLAS6_Theta_EC(double theta){
    double theta_min = 10.0; // lower limit on polar angle, in degrees
    double theta_max = 60.0; // lower limit on polar angle, in degrees

    return (theta >= theta_min && theta < theta_max) ? true : false;
}

// Check that polar angle is within the geometrical acceptance of the TOF.  Input angle in degrees.
bool Cut_CLAS6_Theta_TOF(double theta){
    double theta_min = 10.0; // lower limit on polar angle, in degrees
    double theta_max = 160.0; // lower limit on polar angle, in degrees
    
    return (theta >= theta_min && theta < theta_max) ? true : false;
}

bool Cut_CLAS6_Theta_Ana(TLorentzVector V4, int iAna, string detName){
    bool ret = false;
    
    if(iAna==1){
        if(detName.compare("EC")==0){
            if(Cut_CLAS6_Theta_EC(V4.Theta()* TMath::RadToDeg())) ret = true;
        }
        if(detName.compare("TOF")==0){
            if(Cut_CLAS6_Theta_TOF(V4.Theta()* TMath::RadToDeg())) ret = true;
        }
    }else{
        ret = true;
    }
    return ret;
}

void PrintTLorentzVector(TLorentzVector TLV){
    cout <<"Px "<<TLV.Px()<<"\t";
    cout <<"Py "<<TLV.Py()<<"\t";
    cout <<"Pz "<<TLV.Pz()<<"\t";
    cout <<"E " <<TLV.E() <<"\t";
    cout <<"M " <<TLV.M() <<"\t";
    cout <<"M^2 " <<TLV.M2() <<endl;
}
