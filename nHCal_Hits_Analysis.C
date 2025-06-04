#include "MyConstants.h"

void nHCal_Hits_Analysis(){
//const char strang[]="podio_output"){

  gSystem->Exec("date");
  TString flavor = "nHCal_Hits"; 

  // >>>>> If streaming a runlist from SDCC JLab:
  
  ///////////
  // >>>>> If using local runlist (reading locally stored files):
  TString strang = "nhcal_only_tile5cm_absorber3cm_scintillator0.8cm_11layers_neutron_p1gev_phi45_theta170_10events"; // local runlist
  TString runlist = TString("local_runlists/") + strang + TString("_runlist.txt");
  ///////////

  ///////////
  // streaming runlist (default): 
  //TString runlist = TString("runlists/") + strang + TString("_runlist.txt");  
  ///////////
  
  TString outfile = TString("out.") + strang + TString("-") + flavor + TString(".root");
  TFile *ofile = TFile::Open(outfile,"RECREATE"); // RECREATE overwrites an existing file of the same name

  cout << "Analyzed data is of the type: \n " << strang << " \n";
  cout << "+ Runlist: " << runlist << " \n";
  cout << "+ Output file: " << outfile << " \n";

  TChain *mychain = new TChain("events");

  // if reading a single local file: (locally or streaming)
  //mychain->Add(infile);

  // if reading from a run list (locally or streaming):
  std::ifstream in(runlist);
  std::string file("");
  while (in >> file) mychain->Add(file.data());

  //////////////////////////////////////  
  // Initialize reader
  TTreeReader tree_reader(mychain);
  
  // Get generated particle information (after GEANT; before GEANT is in "GeneratedParticles"):
  TTreeReaderArray<int> partGenStat(tree_reader, "MCParticles.generatorStatus");
  TTreeReaderArray<int> partSimStat(tree_reader, "MCParticles.simulatorStatus");
  TTreeReaderArray<double> partMomX(tree_reader, "MCParticles.momentum.x");
  TTreeReaderArray<double> partMomY(tree_reader, "MCParticles.momentum.y");
  TTreeReaderArray<double> partMomZ(tree_reader, "MCParticles.momentum.z");
  //
  TTreeReaderArray<int> partPdg(tree_reader, "MCParticles.PDG");
  TTreeReaderArray<double> partMass(tree_reader, "MCParticles.mass");
  TTreeReaderArray<double> partEndpointX(tree_reader, "MCParticles.endpoint.x");
  TTreeReaderArray<double> partEndpointY(tree_reader, "MCParticles.endpoint.y");
  TTreeReaderArray<double> partEndpointZ(tree_reader, "MCParticles.endpoint.z");
  TTreeReaderArray<double> partVertexX(tree_reader, "MCParticles.vertex.x");
  TTreeReaderArray<double> partVertexY(tree_reader, "MCParticles.vertex.y");
  TTreeReaderArray<double> partVertexZ(tree_reader, "MCParticles.vertex.z");

  // Get parent and daugther information:
  TTreeReaderArray<int> parents_index(tree_reader, "_MCParticles_parents.index");
  TTreeReaderArray<unsigned int> parents_begin(tree_reader, "MCParticles.parents_begin");
  TTreeReaderArray<unsigned int> parents_end(tree_reader, "MCParticles.parents_end");
  
  TTreeReaderArray<int> daughters_index(tree_reader, "_MCParticles_daughters.index");
  TTreeReaderArray<unsigned int> daughters_begin(tree_reader, "MCParticles.daughters_begin");
  TTreeReaderArray<unsigned int> daughters_end(tree_reader, "MCParticles.daughters_end");  
  
  //// Define Histograms

  //generatorStatus, simulatorStatus, number of daughters
  TH1D *generatorStatus = new TH1D("generatorStatus","Status of generated particles, all; generatorStatus",101,0,100);
  TH1D *simulatorStatus = new TH1D("simulatorStatus","Status of simulated particles, all; simulatorStatus",1000,500000,200000000);
  TH1D *kpmfromphi_simulatorStatus = new TH1D("kpmfromphi_simulatorStatus","Status of simulated K^{#pm} from #phi(1020) decay; simulatorStatus",1000,500000,200000000);
  TH1D *h_idaughters = new TH1D("h_idaughters","Number of daughters, all; Number of daughters",11,0,10);

  // eta (pseudorapidity)
  TH1D *partEta = new TH1D("partEta","Eta of thrown particles; #eta",120,-6.,6.);
  TH1D *recEta = new TH1D("recEta","Eta of reconstructed tracks that have matching thrown particle; #eta",120,-6.,6.);

  TH1D *electronEta = new TH1D("electronEta", "Eta of thrown e; #eta",120,-6.,6.);
  TH1D *protonEta = new TH1D("protonEta","Eta of thrown p;#eta",120,-6.,6.);
  TH1D *muonEta = new TH1D("muonEta","Eta of thrown #mu;#eta",120,-6.,6.);  
  TH1D *pionEta = new TH1D("pionEta","Eta of thrown #pi;#eta",120,-6.,6.);
  TH1D *pi0Eta = new TH1D("pi0Eta","Eta of thrown #pi;#eta",120,-6.,6.);
  TH1D *kaonEta = new TH1D("kaonEta","Eta of thrown K;#eta",120,-6.,6.);
  TH1D *h_idaughters_kaon = new TH1D("h_idaughters_kaon","Number of daughters, K^{#pm}; Number of daughters for K^{#pm}",11,0,10);
  TH1D *EndpointZ_kpm_mupm = new TH1D("EndpointZ_kpm_mupm","generated endpoint.z of K^{#pm} decaying to #mu^{#pm}; endpoint.z [cm]",150,-4000.,4000.);
  TH1D *EndpointZ_kpm_pipm = new TH1D("EndpointZ_kpm_pipm","generated endpoint.z of K^{#pm} decaying to #pi^{#pm}; endpoint.z [cm]",150,-4000.,4000.);
  
  TH1D *rho0Eta = new TH1D("rho0Eta","Eta of thrown #rho^{0};#eta",120,-6.,6.);
  TH1D *pipmfromrho0Eta = new TH1D("pipmfromrho0Eta","generated #eta of pi^{#pm} from rho^{0} decay;#eta",120,-6.,6.);
  
  // eta phi(1020)
  TH1D *phiEta = new TH1D("phiEta","Eta of thrown #phi(1020);#eta",120,-6.,6.);
  TH1D *kpmfromphiEta = new TH1D("kpmfromphiEta","generated #eta of K^{#pm} from #phi(1020) decay;#eta",120,-6.,6.);
  
  TH1D *jpsiEta = new TH1D("jpsiEta","Eta of thrown J/#psi;#eta",120,-6.,6.);
  TH1D *epmfromjpsiEta = new TH1D("epmfromjpsiEta","generated #eta of e^{#pm} from J/#psi decay;#eta",120,-6.,6.);
  
  // momentum
  TH1D *partMom = new TH1D("partMom","Mom of thrown particles; P [GeV]",150,0.,40.);
  TH1D *kpmfromphiMom = new TH1D("kpmfromphiMom","Momentum of thrown K^{#pm} from #phi(1020) decay; p [GeV]",150,0.,40.);

  // endpoint z, decay length, and related studies - for decay kaons from phi(1020)
  TH1D *kpmfromphiEndpointZ = new TH1D("kpmfromphiEndpointZ","generated endpoint.z of K^{#pm} from #phi(1020) decay; endpoint.z [cm]",150,-4000.,4000.);
  TH1D *kpmfromphiEndpointZ_nHCal = new TH1D("kpmfromphiEndpointZ_nHCal","generated endpoint.z of K^{#pm} from #phi(1020) decay in the nHCal acceptance; endpoint.z [cm]",150,-4000.,4000.);
  TH2D *kpmfromphiSimstatus_EndpointZ = new TH2D("kpmfromphiSimstatus_EndpointZ","generated endpoint.z of K^{#pm} from #phi(1020) decay vs. its simulator status; endpoint.z [cm]; simulator status", 150,-4000.,4000., 1000,50000000,200000000);

  // theta (polar angle)
  TH1D *partTheta = new TH1D("partTheta","Theta of thrown charged particles; #theta [rad]",150,0.,3.2);
  // phi (azimuthal angle)
  TH1D *partPhi = new TH1D("partPhi","Phi of thrown charged particles; #phi [rad]",150,-3.2,3.2);

  //// Reset global counters : 
  // count events:
  int ievgen = 0;
  // count total generated particles:
  int ngen_electrons = 0; //+-
  int ngen_protons = 0; //+-
  int ngen_muons = 0; //+-
  int ngen_pions = 0; //+-
  int ngen_pi0=0;
  int ngen_kaons = 0; //+-
  int ngen_rho0 = 0;
  int ngen_rhop = 0;
  int ngen_phi = 0;
  int ngen_omega = 0;
  int ngen_jpsi = 0;
  int ngen_upsilon = 0;
  int ngen_d0 = 0;
  int ngen_b0 = 0;
  // count total number of decays (generated level):
  int ndecay_pi0_gg = 0;
  int ndecay_kpm_mupm = 0;
  int ndecay_kpm_pipm = 0;
  int ndecay_kpm_baryon = 0;
  int ndecay_kpm_em = 0;
  int ndecay_rho0_pp = 0;
  int ndecay_rho0_mumu = 0;
  int ndecay_rho0_ee = 0;
  int ndecay_phi_kk = 0;
  int ndecay_jpsi_mumu = 0;
  int ndecay_jpsi_ee = 0;
  // couunt total number of decay daughters in calo (eta) acceptance:
  int decay_phi_kaonpm_0_nHCal = 0;
  int decay_phi_kaonpm_1_nHCal = 0;
  int decay_phi_kaonpm_2_nHCal = 0;
  // count number of particles with non-zero daughters (gen level):
  int nonzero_daughters_kpm = 0;
  // tag decays on generated particle level ( 0 or 1 for a given generated particle ):
  int is_kpmdecay_mupm = 0;
  int is_kpmdecay_pipm = 0;
  int is_kpmdecay_em = 0;
  int is_kpmdecay_baryon = 0;
  int is_rho0decay_pp = 0; 
  int is_rho0decay_mumu = 0;
  int is_rho0decay_ee = 0;
  int is_phidecay_kk = 0;
  int is_phidecay_kk_reco = 0; // tag that checks if both daughters are reco
  int is_jpsidecay_mumu = 0;
  int is_jpsidecay_ee = 0;
  // count number of each (!) decay daughters in HCalo (eta) acceptances on generated particle (!) level:
  int n_this_decay_phi_kaonpm_k1_rec_nHCal = 0;
  int n_this_decay_phi_kaonpm_k1_rec_bHCal = 0;
  int n_this_decay_phi_kaonpm_k1_rec_lfHCal = 0;
  int n_this_decay_phi_kaonpm_k2_rec_nHCal = 0;
  int n_this_decay_phi_kaonpm_k2_rec_bHCal = 0;
  int n_this_decay_phi_kaonpm_k2_rec_lfHCal = 0;
  for (int k1 = 0; k1 <=3; k1++)
    {
      for (int k2 = 0; k2 <=3; k2++)
	{
	  HCalMatrixphi_kaonpm_rec[k1][k2]=0.;  
	}
    }
  // count total number of decay particles (reco level):
  int ndecay_kpm_mupm_rec = 0; // not yet used
  int ndecay_rho0_pionpm_rec = 0; // not yet used
  int ndecay_phi_kaonpm_rec = 0; 
  int ndecay_jpsi_epm_rec = 0; // not yet used
  int ndecay_jpsi_muonpm_rec = 0; // not yet used
  // count total number of decay particles (reco level) in nHCal eta acceptance:
  int ndecay_kpm_mupm_nHCal = 0;
  int ndecay_rho0_pionpm_nHCal = 0;
  int ndecay_phi_kaonpm_nHCal = 0;
  int ndecay_jpsi_epm_nHCal = 0;
  int ndecay_jpsi_muonpm_nHCal = 0;
  
  cout << "+ Ready to run over events... \n"; 
  
  while(tree_reader.Next()) { // Loop over events

    ievgen++;
    
    //cout << "+ Entering event #: " << ievgen << " \n";
    if(ievgen % 10000 == 0 && ievgen != 0 ){
      cout << "+ processed " << ievgen << " events \n";
      }

    //cout << "Event #: " << ievgen << ", " << partGenStat.GetSize() << " gen particles, " << parents_index.GetSize() << " parent particles, " << daughters_index.GetSize() << " daughter particles \n";   // parent_index and daughter_index must be of the same length since they are in the same tree (is that what pushback does?)

    // start a gigantic loop over the generated particles:
    for(unsigned int i=0; i<partGenStat.GetSize(); i++) // Loop over generated particles
      {
	//cout << "++ Entering generated particle #: " << i << " \n";
	
	int pdg = TMath::Abs(partPdg[i]);
	TVector3 trueMom(partMomX[i],partMomY[i],partMomZ[i]);
	float trueEta = trueMom.PseudoRapidity();
	float truePhi = trueMom.Phi();
	float trueTheta = trueMom.Theta();

	// according to Wouter, "end" is exclusive for daughters - so no "+1" and stop 1 before "end" (can't do it for parents, otherwise we get negative indices)
	int i_parents_begin = parents_begin[i];
        int i_parents_end = parents_end[i];
	int i_parents = parents_end[i] - parents_begin[i];

	int i_daughters_begin = daughters_begin[i]; 
        int i_daughters_end = daughters_end[i] - 1;	
	int i_daughters = daughters_end[i] - daughters_begin[i];

	generatorStatus->Fill(partGenStat[i]);
	simulatorStatus->Fill(partSimStat[i]);
	h_idaughters->Fill(i_daughters);

	// reset particle decays on generated level for each new generated particle:
	is_kpmdecay_mupm = 0;
	is_kpmdecay_pipm = 0;
	is_kpmdecay_em = 0;
	is_kpmdecay_baryon = 0;
	is_rho0decay_pp = 0;
	is_rho0decay_mumu = 0;
	is_rho0decay_ee = 0;
	is_phidecay_kk = 0;
	is_phidecay_kk_reco = 0; // were daughters reco?
	is_jpsidecay_mumu = 0;
	is_jpsidecay_ee = 0;
	n_this_decay_phi_kaonpm_k1_rec_nHCal = 0;
	n_this_decay_phi_kaonpm_k1_rec_bHCal = 0;
	n_this_decay_phi_kaonpm_k1_rec_lfHCal = 0;
	n_this_decay_phi_kaonpm_k2_rec_nHCal = 0;
	n_this_decay_phi_kaonpm_k2_rec_bHCal = 0;
	n_this_decay_phi_kaonpm_k2_rec_lfHCal = 0;
	// tag if for a given decay, all daughters were reconstructed:
	

	// partGenStat[i]== 1: stable; 2: decay; 4: beam particle;  21, 23, 61, 62, 63, 71, ... are [di]quark-related)
    //generatorStatus is set in DD4Hep -> Geant4InputAction::setGeneratorStatus. It is possible that it stays in its initial value, 0, "empty". 
	//Consider only selected generated particles:
	//if( (partGenStat[i] == 1) || (partGenStat[i] == 2) )  // Select only stable or decay particles (I trust the beam particles now... (partGenStat[i] == 4) count them separately?)
	//{
	//if(ievgen==11741)
	//{
	//  cout << "Ev#: " << ievgen << ", P-index: " << i <<", PDG: " << partPdg[i] << ", GenStatus:" << partGenStat[i] << ", z-momentum: " << partMomZ[i] << ", y-momentum: " << partMomY[i] << ", x-momentum: " << partMomX[i] << ", i_parents: "<< i_parents<<", i_daughters: " << i_daughters << ", pb: " << parents_index[i_parents_begin] << ", pe: " << parents_index[i_parents_end] <<  ", db: " << daughters_index[i_daughters_begin] << ", de: " << daughters_index[i_daughters_end] << " \n";
	//}	  

	  // Bookkeeping of decay particles:
	  
	  // *** charged-kaon decays (do NOT request partGenStat[i] == 2 - they are tagged as "stable", yet may have daughters...):
	if( partPdg[i] == 321 || partPdg[i] == -321 ) 
	    {
	      //if( partGenStat[i] != 1 )
	      //{
	      //  cout << "Event " << ievgen << ", generated kaon with GenStatus not 1: " << partGenStat[i] << "  \n";
	      //}
	      if( i_daughters > 0 )
		{
		  nonzero_daughters_kpm++;
		  //cout << "Event " << ievgen << " with gen decaying kpm #: " << partPdg[i] << ", daughter 1:" << partPdg[daughters_index[i_daughters_begin]] << ", i_daughters:" << i_daughters << "  \n";
		  
		  // count the kpm to mupm and pipm decays - for now require exactly 1 daughter (which may mean some decays are escaping):
		  if( i_daughters == 1 && ( partPdg[daughters_index[i_daughters_begin]] == 13 || partPdg[daughters_index[i_daughters_begin]] == -13 ))
		    {
		      ndecay_kpm_mupm++;
		      is_kpmdecay_mupm = 1;
		      EndpointZ_kpm_mupm->Fill(partEndpointZ[i]);
		    } // end of K+-->mu+- decay
		  else if( i_daughters == 1 && ( partPdg[daughters_index[i_daughters_begin]] == 211 || partPdg[daughters_index[i_daughters_begin]] == -211 ))
		    {
		      ndecay_kpm_pipm++;
		      is_kpmdecay_pipm = 1;
		      EndpointZ_kpm_pipm->Fill(partEndpointZ[i]);
		    }// end of K+-->pi+- decay
		  // count the appreances that the first daughter is a baryon (whatever that means...), an electron or a photon:
		  if( partPdg[daughters_index[i_daughters_begin]] == 2112 || partPdg[daughters_index[i_daughters_begin]] == 2212 )
		    {
		      ndecay_kpm_baryon++;
		      is_kpmdecay_baryon = 1;
		    } // end "first daughter is a baryon"
		  else if( partPdg[daughters_index[i_daughters_begin]] == 11 || partPdg[daughters_index[i_daughters_begin]] == -11 || partPdg[daughters_index[i_daughters_begin]] == 22 )
		    {
		      ndecay_kpm_em++;
		      is_kpmdecay_em = 1;
		    }// end "first daughter is electron or photon"
		}// end of non-zero number of daughters

	      //cout << "Event " << ievgen << ", generated kaon with SimStat: " <<  partSimStat[i] << ", ismupm: " << is_kpmdecay_mupm << ", ispipm: " << is_kpmdecay_pipm << ", isem: " << is_kpmdecay_em << ", isbaryon: " << is_kpmdecay_baryon << "  \n";
	      
	    } // end of charged-kaon decays
	// *** rho decays: 
	else if( partPdg[i] == 113 )
	  {
	    //cout << "Event " << ievgen << " with gen rho0 #: " << partPdg[i] << ", daughter 1:" << partPdg[daughters_index[i_daughters_begin]] << ", daughter 2:" << partPdg[daughters_index[i_daughters_begin]+1] << ", i_daughters:" << i_daughters << "  \n";
	    
	    // count the 2-body rho0 decays:
	    if( i_daughters == 2 )
	      {
		if( (partPdg[daughters_index[i_daughters_begin]] == 211 ) && ( partPdg[daughters_index[i_daughters_begin]+1] == -211 ) || (partPdg[daughters_index[i_daughters_begin]] == -211 ) && ( partPdg[daughters_index[i_daughters_begin]+1] == 211) )
		  {
		    //cout << "-> Event " << ievgen << " found rho0 decayed into pi+ pi-: "  << partPdg[daughters_index[i_daughters_begin]] << " and " << partPdg[daughters_index[i_daughters_begin]+1] << "  \n";
		    
		    // count it:
		    ndecay_rho0_pp++;
		    // tag it:
		    is_rho0decay_pp = 1;
		    
		    TVector3 trueMom_rho0_pi1(partMomX[daughters_index[i_daughters_begin]],partMomY[daughters_index[i_daughters_begin]],partMomZ[daughters_index[i_daughters_begin]]);
		    TVector3 trueMom_rho0_pi2(partMomX[daughters_index[i_daughters_begin]+1],partMomY[daughters_index[i_daughters_begin]+1],partMomZ[daughters_index[i_daughters_begin]+1]);
		    TVector3 trueMom_rho0_pi12=trueMom_rho0_pi1 + trueMom_rho0_pi2;
		    
		    float trueEta_rho0_pi1 = trueMom_rho0_pi1.PseudoRapidity();
		    float trueEta_rho0_pi2 = trueMom_rho0_pi2.PseudoRapidity();
		    pipmfromrho0Eta->Fill(trueEta_rho0_pi1);
		    pipmfromrho0Eta->Fill(trueEta_rho0_pi2);
		    
		    //	      cout << "--> Event " << ievgen << " rho0 decay to 2pi: generated rho0 eta: " << trueEta << ", pi1: " << trueEta_rho0_pi1 << ", pi2: " << trueEta_rho0_pi2 << "  \n";
		    //  cout << "            trueMomrho0 X: " << trueMom.X() << ", trueMomrho0 Y: " << trueMom.Y() <<", trueMomrho0 Z: " << trueMom.Z() << "  \n";
		    // cout << "            trueMom_rho0_pi12 X: " << trueMom_rho0_pi12.X() << ", trueMom_rho0_pi12 Y: " << trueMompi12.Y() <<", trueMom_rho0_pi12 Z: " << trueMom_rho0_pi12.Z() << "  \n";
		    
		  } // end of rho0 to pi+pi- decays
		else if( (partPdg[daughters_index[i_daughters_begin]] == 13 ) && ( partPdg[daughters_index[i_daughters_begin]+1] == -13 ) || (partPdg[daughters_index[i_daughters_begin]] == -13 ) && ( partPdg[daughters_index[i_daughters_begin]+1] == 13) )
		  {
		    //count it:
		    ndecay_rho0_mumu++;
		    // tag it:
		    is_rho0decay_mumu = 1;
		  } // end of rho0 to mumu decays
		else if( (partPdg[daughters_index[i_daughters_begin]] == 11 ) && ( partPdg[daughters_index[i_daughters_begin]+1] == -11 ) || (partPdg[daughters_index[i_daughters_begin]] == -11 ) && ( partPdg[daughters_index[i_daughters_begin]] == 11) )
		  {
		    // count it:
		    ndecay_rho0_ee++;
		    // tag it:
		    is_rho0decay_ee = 1;
		  } // end of ee decays
	      } // end of 2-body decays
	  } // end of rho0 decays
	// *** pi0 decays
	else if( partPdg[i] == 111 )
	  {
	    //cout << "Event with gen pi0 #: " << ievgen << ", daughter 1:" << partPdg[daughters_index[i_daughters_begin]] << ", daughter 2:" << partPdg[daughters_index[i_daughters_begin]+1] << ", i_daughters:" << i_daughters << "  \n";
	    if( (partPdg[daughters_index[i_daughters_begin]] == 22) && (partPdg[daughters_index[i_daughters_begin]+1] == 22) )
	      {
		ndecay_pi0_gg++;
	      } // end of gg decays
	  }// end of pi0 decay
	// *** phi(1020) decays:
	else if( partPdg[i] == 333 )
	  {
	    //cout << "Event " << ievgen << " with gen phi(1020): " << partPdg[i] << ", daughter 1:" << partPdg[daughters_index[i_daughters_begin]] << ", daughter 2:" << partPdg[daughters_index[i_daughters_begin]+1] << ", i_daughters:" << i_daughters << "  \n";
	    if( i_daughters != 2 )
	      {
		cout << "Event " << ievgen << " with gen phi(1020), i_daughters= " << i_daughters << " \n";
	      }
	    if( ( partPdg[daughters_index[i_daughters_begin]] != 321 && partPdg[daughters_index[i_daughters_begin]] != -321 ) || ( partPdg[daughters_index[i_daughters_begin]+1] != 321 && partPdg[daughters_index[i_daughters_begin]+1] != -321 ) )
	      {
		cout << "Event " << ievgen << " with gen phi(1020), pdg of daughter1= " << partPdg[daughters_index[i_daughters_begin]] << ", pdg daughter2: " << partPdg[daughters_index[i_daughters_begin]+1] << " \n";
	      }
	    
	    // count the 2-body phi(1020) decays (checked Sartre - they all have exactly 2 daughters, and they ALL decay to kaons, even when all genstatus types are allowed):
	    if( i_daughters == 2 )
	      {
		if( (partPdg[daughters_index[i_daughters_begin]] == 321 ) && ( partPdg[daughters_index[i_daughters_begin]+1] == -321 ) || (partPdg[daughters_index[i_daughters_begin]] == -321 ) && ( partPdg[daughters_index[i_daughters_begin]+1] == 321) )
		  {
		    //cout << "-> Event " << ievgen << " found phi(1020) decayed into K+ K-: "  << partPdg[daughters_index[i_daughters_begin]] << " and " << partPdg[daughters_index[i_daughters_begin]+1] << "  \n";
		    
		    // count it:
		    ndecay_phi_kk++;
		    // tag it:
		    is_phidecay_kk = 1;
		    
		    TVector3 trueMom_phi_k1(partMomX[daughters_index[i_daughters_begin]],partMomY[daughters_index[i_daughters_begin]],partMomZ[daughters_index[i_daughters_begin]]);
		    TVector3 trueMom_phi_k2(partMomX[daughters_index[i_daughters_begin]+1],partMomY[daughters_index[i_daughters_begin]+1],partMomZ[daughters_index[i_daughters_begin]+1]);
		    
		    // phi meson level:
		    TVector3 trueMom_phi_k12=trueMom_phi_k1 + trueMom_phi_k2;
		    
		    // kaon level:
		    float trueEta_phi_k1 = trueMom_phi_k1.PseudoRapidity();
		    float trueEta_phi_k2 = trueMom_phi_k2.PseudoRapidity();
		    kpmfromphiEta->Fill(trueEta_phi_k1);
		    kpmfromphiEta->Fill(trueEta_phi_k2);
		    kpmfromphiMom->Fill(trueMom_phi_k1.Mag());
		    kpmfromphiMom->Fill(trueMom_phi_k2.Mag());
		    kpmfromphi_simulatorStatus->Fill(partSimStat[daughters_index[i_daughters_begin]]);
		    kpmfromphi_simulatorStatus->Fill(partSimStat[daughters_index[i_daughters_begin]+1]);
		    kpmfromphiEndpointZ->Fill(partEndpointZ[daughters_index[i_daughters_begin]]);
		    kpmfromphiEndpointZ->Fill(partEndpointZ[daughters_index[i_daughters_begin]+1]);
		    kpmfromphiSimstatus_EndpointZ->Fill(partEndpointZ[daughters_index[i_daughters_begin]], partSimStat[daughters_index[i_daughters_begin]]);
		    kpmfromphiSimstatus_EndpointZ->Fill(partEndpointZ[daughters_index[i_daughters_begin]+1], partSimStat[daughters_index[i_daughters_begin]+1]);
		    
		    // kaon 1 in nHCal eta acceptance:
		    if(calo_eta_acceptance("nhcal", trueEta_phi_k1))
		      {
			kpmfromphiEndpointZ_nHCal->Fill(partEndpointZ[daughters_index[i_daughters_begin]]);
		      }
		    // kaon 2 in nHCal eta acceptance:
		    if(calo_eta_acceptance("nhcal", trueEta_phi_k2))
		      {
			kpmfromphiEndpointZ_nHCal->Fill(partEndpointZ[daughters_index[i_daughters_begin]+1]);
		      }
		    
		    //cout << "--> Event " << ievgen << " phi(1020) decay to 2K: generated phi eta: " << trueEta << ", K1: " << trueEta_phi_k1 << ", K2: " << trueEta_phi_k2 << "  \n";
		    //cout << "            trueMomphi X: " << trueMom.X() << ", trueMomphi Y: " << trueMom.Y() <<", trueMomphi Z: " << trueMom.Z() << "  \n";
		    //cout << "            trueMom_phi_k12 X: " << trueMom_phi_k12.X() << ", trueMom_phi_k12 Y: " << trueMom_phi_k12.Y() <<", trueMom_phi_k12 Z: " << trueMom_phi_k12.Z() << "  \n";
		    //cout << "           endpoint Z K1: " << partEndpointZ[daughters_index[i_daughters_begin]] << ", endpoint Z K2: " << partEndpointZ[daughters_index[i_daughters_begin]+1] << "  \n";
		    //cout << "           generator status K1: " <<  partGenStat[daughters_index[i_daughters_begin]] << ", simulator status K1: " << partSimStat[daughters_index[i_daughters_begin]]  << ", generator status K2: " << partGenStat[daughters_index[i_daughters_begin]+1] << ", simulator status K2: " << partSimStat[daughters_index[i_daughters_begin]+1] << "  \n";
		    
		    
		  } // end of phi to k+k- decays
	      } // end of 2-body decays
	  } // end of phi meson decays
	// *** jpsi decays:
	else if( partPdg[i] == 443 )
	  {
	    cout << "Event " << ievgen << " with gen J/Psi #: " << partPdg[i] << ", daughter 1:" << partPdg[daughters_index[i_daughters_begin]] << ", daughter 2:" << partPdg[daughters_index[i_daughters_begin]+1] << ", i_daughters:" << i_daughters << "  \n";
	    // count the 2-body jpsi decays:
	    if( i_daughters == 2 )
	      {
		if( (partPdg[daughters_index[i_daughters_begin]] == 11 ) && ( partPdg[daughters_index[i_daughters_begin]+1] == -11 ) || (partPdg[daughters_index[i_daughters_begin]] == -11 ) && ( partPdg[daughters_index[i_daughters_begin]+1] == 11) )
		  {
		    //cout << "-> Event " << ievgen << " found J/Psi decayed into e+ e-: "  << partPdg[daughters_index[i_daughters_begin]] << " and " << partPdg[daughters_index[i_daughters_begin]+1] << "  \n";
		    ndecay_jpsi_ee++;
		    is_jpsidecay_ee = 1;
		    
		    TVector3 trueMom_jpsi_e1(partMomX[daughters_index[i_daughters_begin]],partMomY[daughters_index[i_daughters_begin]],partMomZ[daughters_index[i_daughters_begin]]);
		    TVector3 trueMom_jpsi_e2(partMomX[daughters_index[i_daughters_begin]+1],partMomY[daughters_index[i_daughters_begin]+1],partMomZ[daughters_index[i_daughters_begin]+1]);
		    TVector3 trueMom_jpsi_e12=trueMom_jpsi_e1 + trueMom_jpsi_e2;
		    
		    float trueEta_jpsi_e1 = trueMom_jpsi_e1.PseudoRapidity();
		    float trueEta_jpsi_e2 = trueMom_jpsi_e2.PseudoRapidity();
		    epmfromjpsiEta->Fill(trueEta_jpsi_e1);
		    epmfromjpsiEta->Fill(trueEta_jpsi_e2);
		    
		    //cout << "--> Event " << ievgen << " J/Psi decay to 2e: generated J/Psi eta: " << trueEta << ", e1: " << trueEta_jpsi_e1 << ", e2: " << trueEta_jpsi_e2 << "  \n";
		    //cout << "            trueMomjpsi X: " << trueMom.X() << ", trueMomjpsi Y: " << trueMom.Y() <<", trueMomjpsi Z: " << trueMom.Z() << "  \n";
		    //cout << "            trueMom_jpsi_e12 X: " << trueMom_jpsi_e12.X() << ", trueMom_jpsi_e12 Y: " << trueMom_jpsi_e12.Y() <<", trueMom_jpsi_e12 Z: " << trueMom_jpsi_e12.Z() << "  \n";
		    
		  } // end of jpsi to e+e- decays
		else if( (partPdg[daughters_index[i_daughters_begin]] == 13 ) && ( partPdg[daughters_index[i_daughters_begin]+1] == -13 ) || (partPdg[daughters_index[i_daughters_begin]] == -13 ) && ( partPdg[daughters_index[i_daughters_begin]+1] == 13) )
		  {
		    // count it:
		    ndecay_jpsi_mumu++;
		    // tag it:
		    is_jpsidecay_mumu = 1;
		  } // end of jpsi to mumu decays
	      } // end of 2-body decays
	  } // end of jpsi decays
	
	// end of decay bookkeeping on generated particle level 
	
	
	// count any generated particles (within the type selection we made above):
	if( pdg == 11){
	  ngen_electrons++;
	  electronEta->Fill(trueEta);
	}// electrons                                                                                                     
	else if( pdg == 13){
	  ngen_muons++;
	  muonEta->Fill(trueEta);
	  //cout << "*** Event " << ievgen << ", generated muon \n";    
	}// muons                                                                                                         
	else if( pdg == 211){
	  ngen_pions++;
	  pionEta->Fill(trueEta);
	  //cout << "*** Event " << ievgen << ", generated pi+- \n";    
	}//pions_pm
	else if( pdg == 111){
	  ngen_pi0++;
	  pi0Eta->Fill(trueEta);
	}//pi0  
	else if( pdg == 321 ){
	  ngen_kaons++;
	  kaonEta->Fill(trueEta);
	  h_idaughters_kaon->Fill(i_daughters);
	} // kaons_pm                                                                                                     
	else if( pdg == 113){
	  ngen_rho0++;
	  rho0Eta->Fill(trueEta);
	} // rho(770)    
	else if( pdg == 443){
	  ngen_jpsi++;
	  jpsiEta->Fill(trueEta);
	} // J/Psi(1S)                                                                                                    
	else if( pdg == 2212){
	  ngen_protons++;
	  protonEta->Fill(trueEta);
	}// protons
	else if( pdg == 213){
	  ngen_rhop++;
	}// rhop
	else if( pdg == 333){
	  ngen_phi++;
	  phiEta->Fill(trueEta);
	}// phi(1020)
	else if( pdg == 223){
	  ngen_omega++;
	}// omega(982)
	else if( pdg == 553){
	  ngen_upsilon++;
	}// Upsilon(1S)
	else if( pdg == 421){
	  ngen_d0++;
	}// D0
	else if( pdg == 511){
	  ngen_b0++;
	}// B0
	
	
	//Fill all true eta: 
	partEta->Fill(trueEta);
	
	// Fill all true momentum:
	partMom->Fill(trueMom.Mag());
	
	// Fill all true phi: 
	partPhi->Fill(truePhi);
	
	// Fill all true theta: 
	partTheta->Fill(trueTheta);
	
	//} // End stable or decay particles condition
      } // End loop over thrown particles, within that event

    // XXX Do some things on event level after looping over all generated particles and getting all the info:
    
    
    //for(unsigned int k=0; k<trackMomX.GetSize(); k++){ // Loop over all reconstructed tracks, thrown or not                                 
    
    //  TVector3 recMom(trackMomX[k], trackMomY[k], trackMomZ[k]);
    
    //} // End loop over all reconstructed tracks, within that event 
    
    // now go to next event
    
  } // End loop over events
  
  // Calculate fractions of decay particles:
  double fraction_rho0_pionpm_nHCal = 0.;
  fraction_rho0_pionpm_nHCal = ndecay_rho0_pp?(double(ndecay_rho0_pionpm_nHCal)/(2*double(ndecay_rho0_pp))):0;
  double fraction_phi_kaonpm_rec = 0.;
  fraction_phi_kaonpm_rec = ndecay_phi_kk?(double(ndecay_phi_kaonpm_rec)/(2*double(ndecay_phi_kk))):0;
  double fraction_phi_kaonpm_nHCal = 0.;
  fraction_phi_kaonpm_nHCal = ndecay_phi_kk?(double(ndecay_phi_kaonpm_nHCal)/(2*double(ndecay_phi_kk))):0;
  double fraction_jpsi_epm_nHCal = 0.;
  fraction_jpsi_epm_nHCal = ndecay_jpsi_ee?(double(ndecay_jpsi_epm_nHCal)/(2*double(ndecay_jpsi_ee))):0;
  //
  float HCalMatrixphi_kaonpm_rec_fraction[4][4];
  for (int k1 = 0; k1 <=3; k1++)
    {
      for (int k2 = 0; k2 <=3; k2++)
	{
	  HCalMatrixphi_kaonpm_rec_fraction[k1][k2]=HCalMatrixphi_kaonpm_rec[3][3]?(100*double(HCalMatrixphi_kaonpm_rec[k1][k2])/double(HCalMatrixphi_kaonpm_rec[3][3])):0;  
	}
    }

  
  
  cout << "Number of generated events: " << ievgen << " \n\n";
  cout << "Number of generated electrons +-: " << ngen_electrons << " \n";
  cout << "Number of generated protons +-: " << ngen_protons << " \n";
  cout << "Number of generated muons +-: " << ngen_muons << " \n";
  cout << "Number of generated pions +-: " << ngen_pions << " \n";
  cout << "Number of generated pi0: " << ngen_pi0 << ", of which decay into 2 gamma: " << ndecay_pi0_gg <<  " \n";
  cout << "Number of generated kaons +-: " << ngen_kaons << ", of which have non-zero daughters: " << nonzero_daughters_kpm << ", of which decay into mu+- (1 daughter): " << ndecay_kpm_mupm << ", into pi+- (1 daughter): " << ndecay_kpm_pipm << ", involve an electron or photon: " << ndecay_kpm_em << ", end up as baryon: " << ndecay_kpm_baryon <<" \n";
  cout << "Number of generated rho0: " << ngen_rho0 << ", of which decay into pi+ pi-: " << ndecay_rho0_pp << ", into mu+ mu-: " << ndecay_rho0_mumu << ", into e+ e-: " << ndecay_rho0_ee << " \n";
  cout << "        " << ndecay_rho0_pionpm_nHCal << " reconstructed pi+ pi- make it into the nHCal acceptance, with corresponds to a fraction " << fraction_rho0_pionpm_nHCal << " \n";
  cout << "Number of generated rho+: " << ngen_rhop << " \n";
  cout << "Number of generated phi: " << ngen_phi <<", of which decay into K+ K-: " << ndecay_phi_kk << " \n";
  cout << "         Of these " << 2*ndecay_phi_kk << " decay K+ K-, " << ndecay_phi_kaonpm_rec << " are reconstructed by ePIC (fraction " << fraction_phi_kaonpm_rec << "), and \n ";
  cout << "        " << ndecay_phi_kaonpm_nHCal << " reconstructed K+ K- make it into the nHCal acceptance, with corresponds to a fraction (of generated decay K+ K-) " << fraction_phi_kaonpm_nHCal << " \n";
  cout << "  --> 0 kaons: " << decay_phi_kaonpm_0_nHCal << "\n";
  cout << "  --> 1 kaon:  " << decay_phi_kaonpm_1_nHCal << "\n";
  cout << "  --> 2 kaons: " << decay_phi_kaonpm_2_nHCal << "\n";
  cout << "Distributions in nHCal (0), bHCal (1), lfHCal (2), and all HCals (3), for (K1, K2): \n";
  for (int k1 = 0; k1 <=3; k1++)
    {
      for (int k2 = 0; k2 <=3; k2++)
	{
	  cout << "( K1, K2 ) = ( " << k1 << ", " << k2 << " ) = " << HCalMatrixphi_kaonpm_rec[k1][k2] << " , fraction [%]: " << HCalMatrixphi_kaonpm_rec_fraction[k1][k2] <<"\n";  
	}
    }
  
  cout << "Number of generated omega: " << ngen_omega << " \n";
  cout << "Number of generated J/Psi: " << ngen_jpsi << " , of which decay into e+ e-: " << ndecay_jpsi_ee << ", into mu+ mu-: " << ndecay_jpsi_mumu << " \n";
  cout << "        " << ndecay_jpsi_epm_nHCal << " reconstructed e+ e- make it into the nHCal acceptance, with corresponds to a fraction " << fraction_jpsi_epm_nHCal << " \n";
  cout << "Number of generated Upsilon: " << ngen_upsilon << " \n";
  cout << "Number of generated D0: " << ngen_d0 << " \n";
  cout << "Number of generated B0: " << ngen_b0 << " \n\n";
    
  ofile->Write(); // Write histograms to file
  ofile->Close(); // Close output file

  cout << "Output histograms written in: " << outfile << " \n";
  
  cout << "Thank you for running Caro's macro.\n";
  gSystem->Exec("date");
}
