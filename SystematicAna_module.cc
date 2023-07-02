//  Code developed for estimating charge from SpacePoint object with "reco3d" module but this object is not availble in MUSUN sample so it is not is use for the time being // So now going to developo code for calculating the tracklength using Space Point object : 
//  Date - 24/April/2020/ 
// LArSoft includes
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/Simulation/LArG4Parameters.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/Exception.h"
#include "dune/AnaUtils/DUNEAnaSpacePointUtils.h"
// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib/pow.h" // cet::sum_of_squares()

// ROOT includes. Note: To look up the properties of the ROOT classes,
// use the ROOT web site; e.g.,
// <https://root.cern.ch/doc/master/annotated.html>
#include "TH1.h"
#include "TF2.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TVector3.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TGraph2D.h"
#include "TPolyLine3D.h"
#include"TFile.h"
#include"TCanvas.h"
#include"TLegend.h"
#include "TStyle.h"
#include <cassert>
// C++ inclu
//ides
#include <map>
#include <cmath>

#include <Eigen/Dense>

using namespace Eigen;
#include "Utilsfunction.cc"
using namespace std;

namespace lar {
namespace example {


  
  class SystematicAna : public art::EDAnalyzer
  {
  public:

  
    //
    struct Config {

      // Save some typing:
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      // One Atom for each parameter
      fhicl::Atom<art::InputTag> SimulationLabel {
        Name("SimulationLabel"),
        Comment("tag of the input data product with the detector simulation information")
        };

      fhicl::Atom<art::InputTag> HitLabel {
        Name("HitLabel"),
        Comment("tag of the input data product with reconstructed hits")
        };

      fhicl::Atom<art::InputTag> ClusterLabel {
        Name("ClusterLabel"),
        Comment("tag of the input data product with reconstructed clusters")
        };

      fhicl::Atom<int> PDGcode1 {
        Name("PDGcode1"),
        Comment("particle type (PDG ID) of the primary particle to be selected")
        };

      fhicl::Atom<int> PDGcode2 {
        Name("PDGcode2"),
        Comment("particle type (PDG ID) of the primary particle to be selected")
        };


    }; // Config

  
    //
    using Parameters = art::EDAnalyzer::Table<Config>;

    explicit SystematicAna(Parameters const& config);
    virtual void beginJob() override;
   virtual void beginRun(const art::Run& run) override;
    virtual void analyze (const art::Event& event) override;

  private:

    // The parameters we'll read from the .fcl file.
    art::InputTag fSimulationProducerLabel; ///< The name of the producer that tracked simulated particles through the detector
    art::InputTag fHitProducerLabel;        ///< The name of the producer that created hits
    art::InputTag fClusterProducerLabel;    ///< The name of the producer that created clusters
     art::InputTag fTrackProducerLabel;
     int fSelectedPDG1;          
     int fSelectedPDG2;             ///< PDG code of particle we'll focus on
     double fBinSize;                        ///< For dE/dx work: the value of dx.

    // Pointers to the histograms we'll create.
    TH1D* fPDGCodeHist;     ///< PDG code of all particles
    TH1D* fMomentumHist;    ///< momentum [GeV] of all selected particles
    TH1D* fTrackLengthHist; ///< true length [cm] of all selected particles
    TH1D* fHitIntegralHist; ///< Hit ADC Integral

    // The n-tuples we'll create.
    TTree* fSimulationNtuple;     ///< tuple with simulated data
    TTree* fSimulationTrackTree;
    TTree* fReconstructionNtuple; ///< tuple with reconstructed data


       geo::GeometryCore const *fGeometry;  
   
    // The comment lines with the @ symbols define groups in doxygen.
    /// @name The variables that will go into both n-tuples.
    /// @{
    int fEvent;        ///< number of the event being processed
    int fRun;          ///< number of the run being processed
    int fSubRun;       ///< number of the sub-run being processed
    /// @}

    /// @name The variables that will go into the simulation n-tuple.
    /// @{
    int fSimPDG;       ///< PDG ID of the particle being processed
    int fSimTrackID;   ///< GEANT ID of the particle being processed
     int wireNumber;

    // Arrays for 4-vectors: (x,y,z,t) and (Px,Py,Pz,E).
    // Note: old-style C++ arrays are considered obsolete. However,
    // to create simple n-tuples, we still need to use them.
    double fStartXYZT[4]; ///< (x,y,z,t) of the true start of the particle
    double fEndXYZT[4];   ///< (x,y,z,t) of the true end of the particle
    double fStartPE[4];   ///< (Px,Py,Pz,E) at the true start of the particle
    double fEndPE[4];     ///< (Px,Py,Pz,E) at the true end of the particle


           const double lineXM = -720;   const double lineYM = -600;    const double lineZM = 0;
           const double lineXP = +720;   const double lineYP = +600;    const double lineZP = 6500;    /// Detector planes  position on XYZ axis 

        int tmx, tmy, tmz, tpx, tpy, tpz;
        bool cmx, dmx, cmy, dmy, cmz, dmz;
        bool cpx, dpx, cpy, dpy, cpz, dpz;

           /// Number of dE/dx bins in a given track.
    int fSimNdEdxBins;
    int fRecoPDG;       ///< PDG ID of the particle being processed
    int fRecoTrackID;   

     double RecTrackExitXposmx=0.0, RecTrackExitYposmy=0.0, RecTrackExitZposmz=0.0;
     double RecTrackExitXpospx=0.0, RecTrackExitYpospy=0.0, RecTrackExitZpospz=0.0;
     
     double  Extpl_track_length1=0.0, Stp_muon_track_len=0.0, Extpl_track_length=0.0, trackLength=0.0, SPHitCharge =0.0; 
      std::vector<double> fSimWireEnergy; std::vector<int> fSimWireNumber; std::vector<double> fSimWireCharge; std::vector<double>  fTrackLength;std::vector<double> fStp_muon_Trj_TrackLength; 
      std::vector<double> fStp_muon_TrackLen; std::vector<double> fNon_Stp_Trjtl; std::vector<double> fTrackStartXpos; std::vector<double> fTrackStartYpos; std::vector<double> fTrackStartZpos; 
      std::vector<double> fTrackEndXpos; std::vector<double> fTrackEndYpos; std::vector<double> fTrackEndZpos; std::vector<double >fTrackExitXpos;std::vector<double >fTrackExitYpos;
      std::vector<double >fTrackExitZpos;std::vector<double>fExtpl_track_length;std::vector<double>fNon_Stp_Extpltl;std::vector<double >fRTrackExitXpos; std::vector<double >fRTrackExitYpos; 
      std::vector<double >fRTrackExitZpos;
     std::vector<double> fRecodEdxBins; std::vector<double> fRecoChargeInt; std::vector<double> fRecoChargeADC; std::vector<double> fRecoWireNumber;std::vector<double> fRecoHitCharge;
    std::vector<double> fRecoPeakTime;  std::vector<double> fRec_SpacePoint_X;std::vector<double> fRec_SpacePoint_Y;  std::vector<double> fRec_SpacePoint_Z; std::vector<double> SPLPd_X;
    std::vector<double> SPLPd_Y;  std::vector<double> SPLPd_Z; std::vector<double> fHitCharge_SP;std::vector<double> fSimEnergyDepo; std::vector<double> fHitdQ_track_segment;
     std::vector<double>muonDir; std::vector<double> fMuon_mom;
  
        geo::GeometryCore const* fGeometryService;   ///< pointer to Geometry provider
      //  detinfo::DetectorClocks const* fTimeService; ///< pointer to detector clock time service provider
        double                   fElectronsToGeV;    ///< conversion factor
        int                      fTriggerOffset;     ///< (units of ticks) time of expected neutrino event

  }; // class EnergyAna

 
  SystematicAna::SystematicAna(Parameters const& config)
    : EDAnalyzer(config)
    , fSimulationProducerLabel(config().SimulationLabel())
    , fHitProducerLabel       (config().HitLabel())
    , fClusterProducerLabel   (config().ClusterLabel())
    , fSelectedPDG1            (config().PDGcode1())
    , fSelectedPDG2            (config().PDGcode2())
  {
    fGeometryService = lar::providerFrom<geo::Geometry>();

    consumes<std::vector<simb::MCParticle>>(fSimulationProducerLabel);
    consumes<std::vector<sim::SimChannel>>(fSimulationProducerLabel);
    consumes<art::Assns<simb::MCTruth, simb::MCParticle>>(fSimulationProducerLabel);
    consumes<std::vector<recob::Hit>>(fHitProducerLabel);
    consumes<std::vector<recob::Cluster>>(fClusterProducerLabel);
    consumes<art::Assns<recob::Cluster, recob::Hit>>(fHitProducerLabel);
  }


  //-----------------------------------------------------------------------
  void SystematicAna::beginJob()
  {

    art::ServiceHandle<art::TFileService const> tfs;
    fPDGCodeHist     = tfs->make<TH1D>("pdgcodes",";PDG Code;",                  5000, -2500, 2500);
    fMomentumHist    = tfs->make<TH1D>("mom",     ";particle Momentum (GeV);",    100, 0.,    10.);
    fTrackLengthHist = tfs->make<TH1D>("length",  ";particle track length (cm);", 200, 0, 5000);
    fHitIntegralHist = tfs->make<TH1D>("hitintegral",  ";Hit Integral (sumadc);", 200, 0, 1000.);
   // hdl_TrueVsRec   = tfs->make<TH2F>("hdl_TrueVsRec", "TrueVsRec", 100, 0, 100, 100, 0, 100);
  
    fSimulationNtuple     = tfs->make<TTree>("EnergyAnaSimulation",    "EnergyAna Simulation");
    fSimulationTrackTree     = tfs->make<TTree>("TrackExtrapolation",    "EnergyAna Extrapolation");
    fReconstructionNtuple = tfs->make<TTree>("EnergyAnaReconstruction","EnergyAna Reconstruction");

    
            fSimulationNtuple->Branch("Event",      			 &fEvent,          "Event/I");
            fSimulationNtuple->Branch("SubRun",    		         &fSubRun,         "SubRun/I");
            fSimulationNtuple->Branch("Run",        			 &fRun,            "Run/I");
            fSimulationNtuple->Branch("TrackID",   		         &fSimTrackID,     "TrackID/I");
            fSimulationNtuple->Branch("PDG",        			 &fSimPDG,         "PDG/I");
            fSimulationNtuple->Branch("StartXYZT",  			 fStartXYZT,       "StartXYZT[4]/D");
            fSimulationNtuple->Branch("EndXYZT",    			 fEndXYZT,         "EndXYZT[4]/D");
            fSimulationNtuple->Branch("StartPE",    			 fStartPE,         "StartPE[4]/D");
            fSimulationNtuple->Branch("EndPE",     		         fEndPE,           "EndPE[4]/D");
            fSimulationNtuple->Branch("WireEnergy",   			 &fSimWireEnergy);
            fSimulationNtuple->Branch("WireNumber", 			 &fSimWireNumber); 
            fSimulationNtuple->Branch("WireCharge",                      &fSimWireCharge);

            fSimulationTrackTree->Branch("TrackLength",                  &fTrackLength);
            fSimulationTrackTree->Branch("Stp_muon_TrackLength",         &fStp_muon_Trj_TrackLength);  
            fSimulationTrackTree->Branch("Non_Stp_Trjtl",                &fNon_Stp_Trjtl); 
            fSimulationTrackTree->Branch("TrackExitXpos",                &fTrackExitXpos);
            fSimulationTrackTree->Branch("TrackExitYpos",   		 &fTrackExitYpos);
            fSimulationTrackTree->Branch("TrackExitZpos",  	         &fTrackExitZpos);
            fSimulationTrackTree->Branch("SpacePoint_X",		 &fRec_SpacePoint_X);
            fSimulationTrackTree->Branch("SpacePoint_Y",		 &fRec_SpacePoint_Y);
            fSimulationTrackTree->Branch("SpacePoint_Z",  		 &fRec_SpacePoint_Z);
            fSimulationTrackTree->Branch("SPExtpl_track_length",         &fExtpl_track_length);
            fSimulationTrackTree->Branch("Stp_muon_tracklen",	        &fStp_muon_TrackLen);
           fSimulationTrackTree->Branch("Non_Stp_Extpltl",              &fNon_Stp_Extpltl);
            fSimulationTrackTree->Branch("RTrackExitXpos",  		 &fRTrackExitXpos);
            fSimulationTrackTree->Branch("RTrackExitYpos",  	         &fRTrackExitYpos);
    	    fSimulationTrackTree->Branch("RTrackExitZpos",   		 &fRTrackExitZpos);
            fSimulationTrackTree->Branch("RecoChargeInt",   		 &fRecoChargeInt);
            fSimulationTrackTree->Branch("RecoChargeADC",   		 &fRecoChargeADC);
            fSimulationTrackTree->Branch("RecoWireNumber",   		 &fRecoWireNumber);
            fSimulationTrackTree->Branch("RecoHitCharge",		 &fRecoHitCharge);
            fSimulationTrackTree->Branch("RecoPeakTime",   		 &fRecoPeakTime);
            fSimulationTrackTree->Branch("WireEnergy",   		 &fSimWireEnergy);
            fSimulationTrackTree->Branch("WireNumber", 			 &fSimWireNumber);
            fSimulationTrackTree->Branch("WireCharge",                   &fSimWireCharge);
            fSimulationTrackTree->Branch("SimEnergyDepo",	         &fSimEnergyDepo);
            fSimulationTrackTree->Branch("hitdQ_track_segment",		 &fHitdQ_track_segment);
            fSimulationTrackTree->Branch("muon_directon",		 &muonDir);
            fSimulationTrackTree->Branch("Muon_momentum",                &fMuon_mom);


            fReconstructionNtuple->Branch("Event", 			  &fEvent,          "Event/I");
            fReconstructionNtuple->Branch("SubRun", 			  &fSubRun,         "SubRun/I");
            fReconstructionNtuple->Branch("Run",    			  &fRun,            "Run/I");
            fReconstructionNtuple->Branch("TrackID",			  &fRecoTrackID,    "TrackID/I");
            fReconstructionNtuple->Branch("PDG",    			  &fRecoPDG,        "PDG/I");
        
    }



  void SystematicAna::beginRun(const art::Run& /*run*/)
  {
    art::ServiceHandle<sim::LArG4Parameters const> larParameters;
    fElectronsToGeV = 1./larParameters->GeVToElectrons();
  }

  //-----------------------------------------------------------------------
  void SystematicAna::analyze(const art::Event& event)
  {
    // Start by fetching some basic event information for our n-tuple.
    fEvent  = event.id().event();
    fRun    = event.run();
    fSubRun = event.subRun();

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(event);
     auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(event, clockData);
     art::Handle< std::vector<simb::MCParticle> > particleHandle;
    if (!event.getByLabel(fSimulationProducerLabel, particleHandle))
      {
        throw cet::exception("SystematicAna")
          << " No simb::MCParticle objects in this event - "
          << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
      }

       art::Handle< std::vector<recob::Hit> > hitHandle;
      if (!event.getByLabel(fHitProducerLabel, hitHandle)) return;

         
      fTrackLength.clear();fStp_muon_Trj_TrackLength.clear();fTrackExitXpos.clear();fTrackExitYpos.clear();fTrackExitZpos.clear();fRec_SpacePoint_X.clear();fRec_SpacePoint_Y.clear(); 
      fRec_SpacePoint_Z.clear();SPLPd_X.clear();SPLPd_Y.clear();SPLPd_Z.clear();fExtpl_track_length.clear();fStp_muon_TrackLen.clear(); fNon_Stp_Trjtl.clear(); fNon_Stp_Extpltl.clear();
      fRTrackExitXpos.clear(); fRTrackExitYpos.clear(); fRTrackExitZpos.clear();fSimWireEnergy.clear();fSimWireNumber.clear();fSimWireCharge.clear();fRecoChargeInt.clear();fRecoChargeADC.clear();
     fTrackExitXpos.clear(); fTrackExitYpos.clear(); fTrackExitZpos.clear();
      fRecoPeakTime.clear(); fRecoHitCharge.clear(); fSimEnergyDepo.clear();fHitdQ_track_segment.clear(); fHitCharge_SP.clear();  // clearing hits charge associated with Spacepoint 
      muonDir.clear();fMuon_mom.clear();
//  Clearing Vectors End here :
   
    auto simChannelHandle = event.getValidHandle<std::vector<sim::SimChannel>>(fSimulationProducerLabel);
    std::map< int, const simb::MCParticle* > particleMap;

    for ( auto const& particle : (*particleHandle) )
      {
	fSimTrackID = particle.TrackId();
	particleMap[fSimTrackID] = &particle;
	fSimPDG = particle.PdgCode();
	fPDGCodeHist->Fill( fSimPDG );
       trackLength=0.0;

	// only with information from the primary particles in the
	// event, whose PDG codes match a value supplied in the .fcl file.
	if( particle.Process() != "primary"  &&  ( fSimPDG != fSelectedPDG1 || fSimPDG != fSelectedPDG2 ))
	  continue;
	// A particle has a trajectory, consisting of a set of
	// 4-positions and 4-mommenta.
	const size_t numberTrajectoryPoints = particle.NumberTrajectoryPoints();

//	const int last = numberTrajectoryPoints - 1;
//	const TLorentzVector& positionStart = particle.Position(0);
//	const TLorentzVector& positionEnd   = particle.Position(last);
	const TLorentzVector& momentumStart = particle.Momentum(0);
//	const TLorentzVector& momentumEnd   = particle.Momentum(last);
//	fMomentumHist->Fill( momentumStart.P() );
//	positionStart.GetXYZT( fStartXYZT );
//	positionEnd.GetXYZT( fEndXYZT );
	momentumStart.GetXYZT( fStartPE );
//	momentumEnd.GetXYZT( fEndPE );
    //   fMuon_mom.push_back(momentumStart.Mag());
      TVector3 dirVMC = momentumStart.Vect();
      double dmv = dirVMC.Mag();
     if (dmv==0) continue;
     TVector3 udir = dirVMC*(1.0/dmv);
        TVector3 lastpos(0,0,0);
         bool first = true;
        int tpn = 0;
        for(size_t tp=1; tp<numberTrajectoryPoints; ++tp)
             {
                TVector3 pos_trj = particle.Position(tp).Vect();
               if (pos_trj.X()<lineXP && pos_trj.X()>lineXM && pos_trj.Y()<lineYP && pos_trj.Y()>lineYM &&  pos_trj.Z()<lineZP && pos_trj.Z()>lineZM)
                  { if (!first)
                     {
                      trackLength +=(pos_trj-lastpos).Mag();
                      //Traj_gr->SetPoint(tpn,pos_trj.X(),pos_trj.Y(),pos_trj.Z());
                      tpn++;
                      }
                    else
		      {
			first = false;
		      }
		       lastpos = pos_trj; 
                  }
                }
//	const double trackLength = ( positionEnd - positionStart ).Rho();
          fTrackLength.push_back(trackLength);
         fTrackExitXpos.push_back(lastpos.X()); fTrackExitYpos.push_back(lastpos.Y()); fTrackExitZpos.push_back(lastpos.Z());
	//cout<< "Trajectory track length of Particle is: "<< trackLength << " cm\t"<<"End point => \t"<<lastpos.X()<<",\t"<<lastpos.Y()<<",\t"<<lastpos.Z()<<endl;
         UtilsPrint(trackLength);
	fTrackLengthHist->Fill(trackLength);
 }  // MC particle for loop new      
 // i ##############################################    END of TRUE  TRACK EXTRAPOLATION ############################################            
            TVector3  P0,P1,P2,PEnd, p1, p2, p;
           double dqdl=0.0,dqdlADC=0.0;  
           art::ServiceHandle<geo::Geometry> geom;
            for( auto const& hit : (*hitHandle) )
                    {
                     if ( fGeometryService->SignalType( hit.Channel() ) != geo::kCollection) continue;
                     p.SetXYZ(detProp.ConvertTicksToX(hit.PeakTime(), hit.WireID().Plane, hit.WireID().TPC, hit.WireID().Cryostat ), geom->Wire(hit.WireID()).GetCenter().Z(), 0 );
                //   if(((p-p1).Dot(p2-p1) > 0.0) && ((p-p2).Dot(p1-p2) > 0.0)) { // condition to get hits along a track segments // 
                      dqdl=dqdl+ hit.Integral(); dqdlADC=dqdlADC + hit.SummedADC(); // here * is Dot product 
                 // if(p1.Dot(p2)){cout<<"Dot product syntex "<<"\n";}      
                        }
                     cout<< "Totale Charge from Stopping muon Track Segment : "<<dqdl<<"\t ADC charge =>\t "<<dqdlADC<<"\n";
                       if(dqdl!=0) { fHitdQ_track_segment.push_back(dqdl);}
                         //  P0 = P0 + tl*V;
                           dqdl=0.0;
                           dqdlADC=0.0;
//   *********************************8    Charge and Energy Deposite Analysis  **********************************//  
      double ChargeDepo = 0.0;
      double energyDepositB = 0.0;
      for( auto const& channel : (*simChannelHandle) )
          {
            double energyDepositA = 0.0; //jd
            wireNumber = (int)channel.Channel();//jd ///< Set wire number

            if ( fGeometryService->SignalType( channel.Channel() ) != geo::kCollection )
              continue;
            auto const& timeSlices = channel.TDCIDEMap();
           for ( auto const& timeSlice : timeSlices )
              {
                   auto const& chargeDeposits = timeSlice.second;   //jd  
                 for (auto const& chargeDeposit : chargeDeposits)  //jd
                   {  
                     energyDepositA += chargeDeposit.numElectrons * fElectronsToGeV; // jd
                     energyDepositB += chargeDeposit.numElectrons * fElectronsToGeV; // jd
                     ChargeDepo += chargeDeposit.numElectrons;
                   } // End chargeDeposit loop
              } // For each time slice
                  fSimWireEnergy.push_back(energyDepositA); //jd  // old filling place
                  fSimWireNumber.push_back(wireNumber); //jd     // old filling place 
               //   fSimWireCharge.push_back(ChargeDepo);
                      //    fSimulationNtuple->Fill();
          } // For each SimChannel

        fSimWireCharge.push_back(ChargeDepo);
        fSimEnergyDepo.push_back(energyDepositB);
        double hitCharge = 0; //jd
    // double LifeTimeCurrection = 0;
    
     for ( auto const& hit : (*hitHandle) )
      {
       //  fHitIntegralHist->Fill(hit.Integral());
  
         if ( fGeometryService->SignalType( hit.Channel() ) != geo::kCollection) continue; 
          MF_LOG_DEBUG("SystematicAna") << "Hit in collection plane"<<"JD \t "<<hitCharge<< std::endl;
           hitCharge += (hit.Integral()); 

     MF_LOG_DEBUG("SystematicAna")
              << std::endl;
                 fRecoChargeInt.push_back(hit.Integral());  // We should work with this quantity , 
                 fRecoChargeADC.push_back(hit.SummedADC());
             //     fRecoChargeTrue.push_back(hit.PeakAmplitude());  // this have noise .....
                 fRecoWireNumber.push_back(hit.Channel());
                 fRecoPeakTime.push_back(hit.PeakTime());
        
     MF_LOG_DEBUG("SystematicAna")
                  << "Hit index = " << hit.LocalIndex()
                  << " channel number = " << hit.Channel()
                  << " start TDC tick = " << hit.StartTick()
                  << " end TDC tick = " << hit.EndTick()
                  << " peak TDC tick = " << hit.PeakTime()
                  << " sigma peak time = " << hit.SigmaPeakTime()
              //    << " adjusted start TDC tick = " << fTimeService->TPCTick2TDC(hit.StartTick())
               //   << " adjusted end TDC tick = " << fTimeService->TPCTick2TDC(hit.EndTick())
              //    << " adjusted peak TDC tick = " << fTimeService->TPCTick2TDC(hit.PeakTime())
                  << " time = " << time
                  << std::endl;

      } //// for each Hit
 //  cout<< "Total hit Charge from Raw Hit = "<< TMath::Log(hitCharge) <<"\n";
   fRecoHitCharge.push_back(hitCharge); // old 
  // fRecoChargeTrue.push_back();
// ***************************************     Charge and Energy Deposite Analysis **************************// 
 // }  // if stopping final...............

    fSimulationTrackTree->Fill();
 
//   } // if(SpacePoint)
} // EnergyAna::analyze()

DEFINE_ART_MODULE(SystematicAna)

} // namespace example
} // namespace lar



