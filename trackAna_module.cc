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

#include "Utilsfunction.cc"
using namespace std;

namespace lar {
namespace example {


  
  class trackAna : public art::EDAnalyzer
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

           fhicl::Atom<art::InputTag> SpacePointLabel {
        Name("SpacePointLabel"),
        Comment("tag of the input data product with reconstructed Space Point ")
        };

      fhicl::Atom<art::InputTag> TrackLabel {
        Name("TrackLabel"),
        Comment("tag of the input data product with reconstructed hits")
        };   
 
      fhicl::Atom<int> PDGcode1 {
        Name("PDGcode1"),
        Comment("particle type (PDG ID) of the primary particle to be selected")
        };

      fhicl::Atom<int> PDGcode2 {
        Name("PDGcode2"),
        Comment("particle type (PDG ID) of the primary particle to be selected")
        };


   //   fhicl::Atom<double> BinSize {
   //     Name("BinSize"),
     //   Comment("dx [cm] used for the dE/dx calculation")
      //  };

    }; // Config

  
    //
    using Parameters = art::EDAnalyzer::Table<Config>;

    explicit trackAna(Parameters const& config);
    virtual void beginJob() override;
   virtual void beginRun(const art::Run& run) override;
    virtual void analyze (const art::Event& event) override;

  private:

    // The parameters we'll read from the .fcl file.
    art::InputTag fSimulationProducerLabel; ///< The name of the producer that tracked simulated particles through the detector
    art::InputTag fHitProducerLabel;        ///< The name of the producer that created hits
    art::InputTag fClusterProducerLabel;    ///< The name of the producer that created clusters
     art::InputTag fSpacePointProducerLabel;   
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
     std::vector<double> fRecodEdxBins; std::vector<double> fRecoChargeInt; std::vector<double> fRecoChargeTrue; std::vector<double> fRecoWireNumber;std::vector<double> fRecoHitCharge;
    std::vector<double> fRecoPeakTime;  std::vector<double> fRec_SpacePoint_X;std::vector<double> fRec_SpacePoint_Y;  std::vector<double> fRec_SpacePoint_Z; std::vector<double> SPLPd_X;
    std::vector<double> SPLPd_Y;  std::vector<double> SPLPd_Z; std::vector<double> fHitCharge_SP;std::vector<double> fSimEnergyDepo; std::vector<double> fHitdQ_track_segment;
     std::vector<double>muonDir; std::vector<double> fMuon_mom;
  
        geo::GeometryCore const* fGeometryService;   ///< pointer to Geometry provider
      //  detinfo::DetectorClocks const* fTimeService; ///< pointer to detector clock time service provider
        double                   fElectronsToGeV;    ///< conversion factor
        int                      fTriggerOffset;     ///< (units of ticks) time of expected neutrino event

  }; // class EnergyAna

 
  trackAna::trackAna(Parameters const& config)
    : EDAnalyzer(config)
    , fSimulationProducerLabel(config().SimulationLabel())
    , fHitProducerLabel       (config().HitLabel())
    , fClusterProducerLabel   (config().ClusterLabel())
    , fSpacePointProducerLabel (config().SpacePointLabel())
    , fTrackProducerLabel       (config().TrackLabel())
    , fSelectedPDG1            (config().PDGcode1())
    , fSelectedPDG2            (config().PDGcode2())
   // , fBinSize                (config().BinSize())
  {
    // Get a pointer to the geometry service provider.
    fGeometryService = lar::providerFrom<geo::Geometry>();
    // The same for detector TDC clock services.
   // fTimeService = lar::providerFrom<detinfo::DetectorClocksService>();
    // Access to detector properties.
  //  const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  //  fTriggerOffset = detprop->TimeOffsetU();

    consumes<std::vector<simb::MCParticle>>(fSimulationProducerLabel);
    consumes<std::vector<sim::SimChannel>>(fSimulationProducerLabel);
    consumes<art::Assns<simb::MCTruth, simb::MCParticle>>(fSimulationProducerLabel);
    consumes<std::vector<recob::Hit>>(fHitProducerLabel);
    consumes<std::vector<recob::Cluster>>(fClusterProducerLabel);
    consumes<art::Assns<recob::Cluster, recob::Hit>>(fHitProducerLabel);
    consumes<std::vector<recob::SpacePoint>>(fSpacePointProducerLabel);
    consumes<std::vector<recob::Track>>(fTrackProducerLabel);
  }


  //-----------------------------------------------------------------------
  void trackAna::beginJob()
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
            fSimulationTrackTree->Branch("RecoChargeTrue",   		 &fRecoChargeTrue);
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



  void trackAna::beginRun(const art::Run& /*run*/)
  {
    // How to convert from number of electrons to GeV. The ultimate
    // source of this conversion factor is
    // ${LARCOREOBJ_INC}/larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h.
    // But sim::LArG4Parameters might in principle ask a database for it.
    art::ServiceHandle<sim::LArG4Parameters const> larParameters;
    fElectronsToGeV = 1./larParameters->GeVToElectrons();
  }

  //-----------------------------------------------------------------------
  void trackAna::analyze(const art::Event& event)
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
        throw cet::exception("trackAna")
          << " No simb::MCParticle objects in this event - "
          << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
      }

       art::Handle< std::vector<recob::Hit> > hitHandle;
      if (!event.getByLabel(fHitProducerLabel, hitHandle)) return;

       TGraph2D * SP_gr = new TGraph2D();
       TGraph2D * Traj_gr = new TGraph2D();
   // TGraph2D * SP_AvePos_gr = new TGraph2D();
    // TGraph2D *SP_TracePoint_gr = new TGraph2D();
   //  TGraph2D *SPLPd_gr = new TGraph2D();
   //   Clearing All the vectors here 
         
      fTrackLength.clear();fStp_muon_Trj_TrackLength.clear();fTrackExitXpos.clear();fTrackExitYpos.clear();fTrackExitZpos.clear();fRec_SpacePoint_X.clear();fRec_SpacePoint_Y.clear(); 
      fRec_SpacePoint_Z.clear();SPLPd_X.clear();SPLPd_Y.clear();SPLPd_Z.clear();fExtpl_track_length.clear();fStp_muon_TrackLen.clear(); fNon_Stp_Trjtl.clear(); fNon_Stp_Extpltl.clear();
      fRTrackExitXpos.clear(); fRTrackExitYpos.clear(); fRTrackExitZpos.clear();fSimWireEnergy.clear();fSimWireNumber.clear();fSimWireCharge.clear();fRecoChargeInt.clear();fRecoChargeTrue.clear();
     fTrackExitXpos.clear(); fTrackExitYpos.clear(); fTrackExitZpos.clear();
      fRecoPeakTime.clear(); fRecoHitCharge.clear(); fSimEnergyDepo.clear();fHitdQ_track_segment.clear(); fHitCharge_SP.clear();  // clearing hits charge associated with Spacepoint 
      muonDir.clear();fMuon_mom.clear();
//  Clearing Vectors End here :
   
      art::Handle< std::vector<recob::SpacePoint> > recobspacepoints;
    event.getByLabel(fSpacePointProducerLabel, recobspacepoints);
    if (!recobspacepoints->empty())
   {
    
    auto simChannelHandle = event.getValidHandle<std::vector<sim::SimChannel>>(fSimulationProducerLabel);
    std::map< int, const simb::MCParticle* > particleMap;

    for ( auto const& particle : (*particleHandle) )
      {
	fSimTrackID = particle.TrackId();
	particleMap[fSimTrackID] = &particle;
	fSimPDG = particle.PdgCode();
	fPDGCodeHist->Fill( fSimPDG );
       trackLength=0.0;

	// For this example, we want to fill the n-tuples and histograms
	// only with information from the primary particles in the
	// event, whose PDG codes match a value supplied in the .fcl file.
	if( particle.Process() != "primary"  &&  ( fSimPDG != fSelectedPDG1 || fSimPDG != fSelectedPDG2 ))
	  continue;
	// A particle has a trajectory, consisting of a set of
	// 4-positions and 4-mommenta.
	const size_t numberTrajectoryPoints = particle.NumberTrajectoryPoints();

	// For trajectories, as for vectors and arrays, the first
	// point is #0, not #1.
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
  //    cout<< "Start Momentum of the muons \t"<< momentumEnd.Mag()<<"\n";   
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
                      Traj_gr->SetPoint(tpn,pos_trj.X(),pos_trj.Y(),pos_trj.Z());
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
 // ##########   Reconstructed Space Point and Track Length Estimation ###############
            art::ServiceHandle<geo::Geometry> geom;
            dune_ana::DUNEAnaSpacePointUtils SPcharge;
            art::Handle< std::vector<recob::SpacePoint> > recobspacepoints;
            event.getByLabel(fSpacePointProducerLabel, recobspacepoints);
              
         for(size_t isp=0;isp<recobspacepoints->size(); ++isp)
	      {  
                 art::Ptr<recob::SpacePoint> recSPoint(recobspacepoints, isp);
                 const recob::SpacePoint& RecSP = *recSPoint; 
             
    // std::vector<art::Ptr<recob::Hit>> spacePointHits = ::dune_ana::DUNEAnaSpacePointUtils::GetHits(recSPoint, event, "pandora");
        std::vector<art::Ptr<recob::Hit>> spacePointHits  =  SPcharge.GetHits(recSPoint, event, "pandora");       
         
       for ( auto const& SPhits : (spacePointHits) )
            {
                  //  auto hitChannelNumber = hit.Channel();
              auto HitCharge = SPhits->Integral();  
              SPHitCharge = SPHitCharge + HitCharge; 
             if ( fGeometryService->SignalType( SPhits->Channel() ) != geo::kCollection)          
                  continue;	 
                fHitCharge_SP.push_back(SPhits->Integral());
     //  double SP_hit_x = detProp.ConvertTicksToX(SPhits->PeakTime(), SPhits->WireID().Plane, SPhits->WireID().TPC, SPhits->WireID().Cryostat );
     //  double SP_hit_z = geom->Wire(SPhits->WireID()).GetCenter().Z();
    //   cout<<"Hits and Spacepoint Z and X coordinate is : "<< SP_hit_z<<"\t"<<RecSP.XYZ()[2]<<"\t X => "<<SP_hit_x<<"\t"<<RecSP.XYZ()[0]<<"\n";
              }
               // fHitCharge_SP.push_back(spacePointHits.Integral());
               SP_gr->SetPoint(isp,RecSP.XYZ()[0],RecSP.XYZ()[1],RecSP.XYZ()[2]);		
               
                       fRec_SpacePoint_X.push_back(RecSP.XYZ()[0]);
                       fRec_SpacePoint_Y.push_back(RecSP.XYZ()[1]);
                       fRec_SpacePoint_Z.push_back(RecSP.XYZ()[2]);
                } // for () 
     //  fHitCharge_SP.push_back(SPHitCharge);
      SPHitCharge = 0.0;

       //graphplot( SP_gr);  // check Utilsfunction.cc  
 // $$i$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  PCA   ##########################################

      double* xyz[3];
     TVector3 pos(0,0,0);
    TVector3 dir;
      xyz[0] = fRec_SpacePoint_X.data();
      xyz[1] = fRec_SpacePoint_Y.data();
      xyz[2] = fRec_SpacePoint_Z.data();
     
    size_t npts = fRec_SpacePoint_X.size();
    
         if (npts < 2)
            {
               throw cet::exception("tpcvechitfinder2_module.cc: too few TPCClusters to fit a line in linefit");
            }
      
      TMatrixDSym covmat(3);  // covariance matrix (use symmetric version)
     // position is just the average of the coordinates
      double psum[3] = {0,0,0};
        for (size_t ipoint=0; ipoint<npts; ++ipoint)
              { 
               for(size_t j=0; j<3; j++) 
                  {
                    psum[j] += xyz[j][ipoint];
                   }
                 }
       for (size_t j=0; j<3; ++j)
        {
          psum[j] /= npts;
        }
       pos.SetXYZ(psum[0],psum[1],psum[2]);

     for(size_t i=0; i<3; ++i)
        {
          for (size_t j=0; j<= i; ++j)
            {
              double csum=0;
              for (size_t ipoint=0; ipoint<npts; ++ipoint)
                {
                  csum += (xyz[i][ipoint] - psum[i]) * (xyz[j][ipoint] - psum[j]);
                }
              csum /= (npts-1);
              covmat[i][j] = csum;
              covmat[j][i] = csum;
            }
        }
      TVectorD eigenvalues(3);
      TMatrixD eigenvectors = covmat.EigenVectors(eigenvalues);
      
     double dirv[3] = {0,0,0};
      for (size_t i=0; i<3; ++i)
        {
          dirv[i]=eigenvectors[i][0]; 
         }
      dir.SetXYZ(dirv[0],dirv[1],dirv[2]);
  
  //   cout<<"Ave Space Point Position P  =(\t "<<pos[0]<<",\t"<<pos[1]<<",\t"<<pos[2]<<"\t)"<<endl;
 //  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$   PCA  END ################################################               
 //  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$         Track Length Estimation from Space Point  $$$$$  
  
     double dm = dir.Mag(); 
   //  if (dm == 0) continue;
    dir *= (1.0/dm);

 tmx=(lineXM - pos[0])/(dir[0]);  tmy=(lineYM- pos[1])/(dir[1]); tmz=(lineZM- pos[2])/(dir[2]);
 tpx=(lineXP - pos[0])/(dir[0]);  tpy=(lineYP- pos[1])/(dir[1]); tpz=(lineZP- pos[2])/(dir[2]);

 RecTrackExitXposmx  = lineXM;   RecTrackExitYposmy  = lineYM;   RecTrackExitZposmz  = lineZM;
 RecTrackExitXpospx  = lineXP;   RecTrackExitYpospy  = lineYP;   RecTrackExitZpospz  = lineZP;


double RecTrackExitYposmx=pos[1]+dir[1]*tmx; double RecTrackExitXposmy=pos[0]+dir[0]*tmy; double RecTrackExitXposmz=pos[0]+dir[0]*tmz;
double RecTrackExitZposmx=pos[2]+dir[2]*tmx; double RecTrackExitZposmy=pos[2]+dir[2]*tmy; double RecTrackExitYposmz=pos[1]+dir[1]*tmz;
double RecTrackExitYpospx=pos[1]+dir[1]*tpx; double RecTrackExitXpospy=pos[0]+dir[0]*tpy; double RecTrackExitXpospz=pos[0]+dir[0]*tpz;
double RecTrackExitZpospx=pos[2]+dir[2]*tpx; double RecTrackExitZpospy=pos[2]+dir[2]*tpy; double RecTrackExitYpospz=pos[1]+dir[1]*tpz;

 if((RecTrackExitYposmx >= -600. && RecTrackExitYposmx <= 600.) && (RecTrackExitZposmx >= 0. && RecTrackExitZposmx <= 6500.) )
    { fRTrackExitXpos.push_back( RecTrackExitXposmx ); fRTrackExitYpos.push_back( RecTrackExitYposmx );fRTrackExitZpos.push_back(RecTrackExitZposmx);
     }
if((RecTrackExitXposmy >= -720. && RecTrackExitXposmy <= 720.) && (RecTrackExitZposmy >= 0. && RecTrackExitZposmy <= 6500.))
    {fRTrackExitXpos.push_back( RecTrackExitXposmy);fRTrackExitYpos.push_back(RecTrackExitYposmy);fRTrackExitZpos.push_back(RecTrackExitZposmy);
    }
if((RecTrackExitXposmz >= -720. && RecTrackExitXposmz <= 720.) && (RecTrackExitYposmz >= -600. && RecTrackExitYposmz <= 600.))
 {fRTrackExitXpos.push_back(RecTrackExitXposmz);fRTrackExitYpos.push_back(RecTrackExitYposmz);fRTrackExitZpos.push_back(RecTrackExitZposmz);
 }
 if((RecTrackExitYpospx >= -600. && RecTrackExitYpospx <= 600.) && (RecTrackExitZpospx >= 0. && RecTrackExitZpospx <= 6500.))
 {fRTrackExitXpos.push_back( RecTrackExitXpospx );fRTrackExitYpos.push_back( RecTrackExitYpospx );fRTrackExitZpos.push_back( RecTrackExitZpospx );
 }
 if((RecTrackExitXpospy >= -720. && RecTrackExitXpospy <= 720.) && (RecTrackExitZpospy >= 0. && RecTrackExitZpospy <= 6500.))
                                                                                                                                                    
 {fRTrackExitXpos.push_back( RecTrackExitXpospy );fRTrackExitYpos.push_back(RecTrackExitYpospy);fRTrackExitZpos.push_back(RecTrackExitZpospy); 
   }
 if( (RecTrackExitXpospz >= -720. && RecTrackExitXpospz <= 720.) && (RecTrackExitYpospz >= -600. && RecTrackExitYpospz <= 600.) )
  { fRTrackExitXpos.push_back(RecTrackExitXpospz);fRTrackExitYpos.push_back( RecTrackExitYpospz );fRTrackExitZpos.push_back(RecTrackExitZpospz);
   }
          /*
  for(size_t i=0; i<2; ++i)
  { std::cout<<"End point of the PCA_Tracks = (\t "<<fRTrackExitXpos[i]<< ",\t"<<fRTrackExitYpos[i]<<",\t"<<fRTrackExitZpos[i]<<")\t"<<std::endl;
        }*/ 

 Extpl_track_length1 =sqrt(pow(fRTrackExitXpos[0]-fRTrackExitXpos[1], 2)+pow(fRTrackExitYpos[0]-fRTrackExitYpos[1], 2)+pow(fRTrackExitZpos[0]-fRTrackExitZpos[1], 2) * 1.0);

    fExtpl_track_length.push_back(Extpl_track_length1);
    std::cout<<"Extrapolated Track Length from Space-Point is = \t "<< Extpl_track_length1<< "\t cm "<<std::endl;
   
  //  SP_AvePos_gr->SetPoint(0,pos[0],pos[1],pos[2]);
  //  SP_AvePos_gr->SetPoint(1,fRTrackExitXpos[0],fRTrackExitYpos[0],fRTrackExitZpos[0]);
  //  SP_AvePos_gr->SetPoint(2,fRTrackExitXpos[1],fRTrackExitYpos[1],fRTrackExitZpos[1]);


 
    //    int const tl(100.);   // Track segment size set here ex 100 cm, events have track length then 100 cm will be removed  
        TVector3 V = dir.Unit();
        TVector3  P0,P1,P2,PEnd, p1, p2, p,Ptrj,U,SP1,SP2;
        SP1.SetXYZ(fRec_SpacePoint_X.front(), fRec_SpacePoint_Y.front(), fRec_SpacePoint_Z.front());
        SP2.SetXYZ(fRec_SpacePoint_X.back(), fRec_SpacePoint_Y.back(), fRec_SpacePoint_Z.back());
        P0.SetXYZ(fRTrackExitXpos[0],fRTrackExitYpos[0],fRTrackExitZpos[0]);
        PEnd.SetXYZ(fRTrackExitXpos[1],fRTrackExitYpos[1],fRTrackExitZpos[1]);  // any end point can be used either 
        Ptrj.SetXYZ(fTrackExitXpos.back(),fTrackExitXpos.back(),fTrackExitXpos.back());
       bool stopping_muon_E1 = false;  // CASE 1 stopping muons events 
       bool stopping_muon_E2 = false;
// Block for start MC momentum 
        for ( auto const& particle : (*particleHandle) )
        {   
            if( particle.Process() != "primary"  &&  ( fSimPDG != fSelectedPDG1 || fSimPDG != fSelectedPDG2 ))
          continue;
         const size_t numberTrajectoryPoints = particle.NumberTrajectoryPoints();
         const int last = numberTrajectoryPoints - 1;     //  code bloc for stopping muons based on truth info 
          const TLorentzVector& positionStart = particle.Position(0);
          const TLorentzVector& positionEnd   = particle.Position(last);
     if (positionStart.X()<lineXP && positionStart.X()>lineXM && positionStart.Y()<lineYP && positionStart.Y()>lineYM &&  positionStart.Z()<lineZP && positionStart.Z()>lineZM) {stopping_muon_E1 = true;}
     if (positionEnd.X()<lineXP && positionEnd.X()>lineXM && positionEnd.Y()<lineYP && positionEnd.Y()>lineYM &&  positionEnd.Z()<lineZP && positionEnd.Z()>lineZM) {stopping_muon_E2 = true;}

              for(size_t tp=1; tp<numberTrajectoryPoints; ++tp)
                {
                   TVector3 pos_tpc = particle.Position(tp).Vect();
                       if (pos_tpc.X()<lineXP && pos_tpc.X()>lineXM && pos_tpc.Y()<lineYP && pos_tpc.Y()>lineYM &&  pos_tpc.Z()<lineZP && pos_tpc.Z()>lineZM)
                      { 
                        fMuon_mom.push_back(particle.Momentum(tp).Vect().Mag());                 
                 //     cout<<"Mpmentur of MC particle  =>\t "<< particle.Momentum(tp).Vect().Mag()<<"\n";
                     break;
                      }
                 } // for loop
         }  // particle 

     StoppingMuon_track( stopping_muon_E1, stopping_muon_E2);
  // cout<< "Trajectory track end point of Particle is: cm\t"<<"End point => \t"<<fTrackExitXpos.back()<<",\t"<<fTrackExitYpos.back()<<",\t"<<fTrackExitZpos.back()<<endl;
// Block for Stopping Muons Analysis $$$

     stopping_muon_E1 = false;  // CASE 2 Stopping Muons events 
     stopping_muon_E2 = false;
    int N_sp1(0), N_sp2(0);   
      U = PEnd - P0;    //  3D projection algo source : https://stackoverflow.com/questions/9368436/3d-perpendicular-point-on-line-from-3d-point
     
       for(size_t isp=0;isp<recobspacepoints->size(); ++isp)
            { art::Ptr<recob::SpacePoint> recSPoint(recobspacepoints, isp);
             const recob::SpacePoint& RecSP = *recSPoint;
            TVector3 Q1 = RecSP.XYZ();
            TVector3 P0Q = Q1 - P0;
            TVector3 W2_1 = P0Q - U*(P0Q.Dot(U)/(U.Mag()*U.Mag()));
            TVector3 PL1 = Q1 - W2_1;
         //   SP_gr->SetPoint(isp,PL1.X(),PL1.Y(), PL1.Z());
          //  SP_AvePos_gr->SetPoint(isp, RecSP.XYZ()[0], RecSP.XYZ()[1], RecSP.XYZ()[2]);  
             double SP_dis1 = (PL1 - P0).Mag();
             if(SP_dis1 <= 5) { N_sp1++; if(N_sp1 == 4 ) {stopping_muon_E1 = true; break;}} // look for 5 space points from end points 
            } // for(); 

        for(size_t isp= 0; isp<recobspacepoints->size();++isp)
            { art::Ptr<recob::SpacePoint> recSPoint(recobspacepoints, isp);
             const recob::SpacePoint& RecSP = *recSPoint;
            TVector3 Q2 = RecSP.XYZ();
            TVector3 PEndQ = Q2 - PEnd;
            TVector3 W2_2 = PEndQ - U*(PEndQ.Dot(U)/(U.Mag()*U.Mag()));
            TVector3 PL2 = Q2 - W2_2;  
           //  SP_gr->SetPoint(isp,PL2.X(),PL2.Y(), PL2.Z());
            double SP_dis2 = (PL2 - PEnd).Mag();
            if( SP_dis2 <= 5 ) { N_sp2++; if(N_sp2 == 4) {stopping_muon_E2 = true;break;}} // look for 5 space points from second end points      
         //  SP_AvePos_gr->SetPoint(isp, RecSP.XYZ()[0], RecSP.XYZ()[1], RecSP.XYZ()[2]); 
         } // for();                    
  
        //StoppingMuon_track( stopping_muon_E1, stopping_muon_E2);
        stopping_muon_E1 = false;   // CASE 3 Stopping muons events
          stopping_muon_E2 = false;
          if((Ptrj - P0).Mag()>=10) {stopping_muon_E1 = true;}   // Stopping point 10cm  away from Endpoint
             if((Ptrj - PEnd).Mag()>=10) {stopping_muon_E2 = true;}
            StoppingMuon_track( stopping_muon_E1, stopping_muon_E2);
          
           stopping_muon_E1 = false;  // CASE 4 Stopping muons events
           stopping_muon_E2 = false;
         if (SP1.X()<lineXP && SP1.X()>lineXM && SP1.Y()<lineYP && SP1.Y()>lineYM &&  SP1.Z()<lineZP && SP1.Z()>lineZM) {stopping_muon_E1 = true;}  
         if (SP2.X()<lineXP && SP2.X()>lineXM && SP2.Y()<lineYP && SP2.Y()>lineYM &&  SP2.Z()<lineZP && SP2.Z()>lineZM) {stopping_muon_E2 = true;}  
           StoppingMuon_track( stopping_muon_E1, stopping_muon_E2);

    /*
            P1 = P0 + tl*V;
            P2 = P0 - tl*V;
          double d01 = (PEnd - P1).Mag();
          double d02 = (PEnd - P2).Mag();
          double tl_trace = (P1 - P0).Mag();double dqdl=0.0,dqdlADC=0.0;
         // double d0 = tl_trace;
          int n= 0;




   if(StoppingMuon_track( stopping_muon_E1, stopping_muon_E2))  // XOR logic to pass the first level cuts   
     { // cout<<"This is a stopping events :"<<"Trajectory Track length = \t"<<fTrackLength.back()<<"\t Reconstructed track length =\t"<<Extpl_track_length1<<" \t mometum \t"<<fMuon_mom.back()<<"\n";
                 
                      fStp_muon_TrackLen.push_back(Extpl_track_length1);
                      fStp_muon_Trj_TrackLength.push_back(fTrackLength.back()); 

        if(d01 < d02 )
           {
            while( tl_trace < Extpl_track_length1)
                  { //  SP_TracePoint_gr->SetPoint(n,P0.X(),P0.Y(),P0.Z());
                    //  cout<< "Point at the line (Forward Dir ) is : (\t"<<P0(0)<<",\t"<<P0(1)<<",\t"<<P0(2)<<"\t"<<"distance = \t"<<tl_trace<<"\n";
                    p1.SetXYZ(P0.X(),P0.Z(),0);
                    p2.SetXYZ((P0+tl*V).X(),(P0+tl*V).Z(),0);

                   for( auto const& hit : (*hitHandle) )
                       {
                     if ( fGeometryService->SignalType( hit.Channel() ) != geo::kCollection) continue;
                     p.SetXYZ(detProp.ConvertTicksToX(hit.PeakTime(), hit.WireID().Plane, hit.WireID().TPC, hit.WireID().Cryostat ), geom->Wire(hit.WireID()).GetCenter().Z(), 0 );
                     
                    if(((p-p1).Dot(p2-p1) > 0.0) && ((p-p2).Dot(p1-p2) > 0.0)) { dqdl=dqdl+ hit.Integral(); dqdlADC=dqdlADC + hit.SummedADC();} // here * is Dot product 
                 // if(p1.Dot(p2)){cout<<"Dot product syntex "<<"\n";}      
               }
             cout<< "Totale Charge from Stopping muon Track Segment : "<<dqdl/tl<<"\t ADC charge =>\t "<<dqdlADC/tl<<"\n";
           if(dqdl!=0) { fHitdQ_track_segment.push_back(dqdl/tl);}
            P0 = P0 + tl*V;
            dqdl=0.0;
            dqdlADC=0.0;
            n++;
            tl_trace = tl_trace + 100;
           }    
        }  // if (d01 or d02)
   else // (d01 not less than d02)
      {
     while( tl_trace < Extpl_track_length1)
      {  // SP_TracePoint_gr->SetPoint(n,P0.X(),P0.Y(),P0.Z());
         //  cout<< "Point at the line (Reverse dir ) is : (\t"<<P0(0)<<",\t"<<P0(1)<<",\t"<<P0(2)<<")"<<"distance = \t"<<tl_trace<<"\n";

          p1.SetXYZ(P0.X(),P0.Z(),0);
          p2.SetXYZ((P0-tl*V).X(),(P0-tl*V).Z(),0);
              for( auto const& hit : (*hitHandle) )
               {
                   if (fGeometryService->SignalType( hit.Channel() ) != geo::kCollection)  continue;
              p.SetXYZ(detProp.ConvertTicksToX(hit.PeakTime(), hit.WireID().Plane, hit.WireID().TPC, hit.WireID().Cryostat ), geom->Wire(hit.WireID()).GetCenter().Z(),0 );
                 if(((p-p1).Dot(p2-p1) > 0.0) && ((p-p2).Dot(p1-p2) > 0.0)) { dqdl=dqdl+ hit.Integral(); dqdlADC=dqdlADC + hit.SummedADC();}
                 }
               cout<< "Totale Charge from Stopping muon Track Segment : "<<dqdl/tl<<" \t ADC charge =>\t"<<dqdlADC/tl<<"\n";
             //  dqdl=dqdl+HitCharge;                                                                    
             //               //  break;           
             if(dqdl!=0)  { fHitdQ_track_segment.push_back(dqdl/tl);}
                P0 = P0 - tl*V;
                dqdl=0.0;
                dqdlADC=0.0;
                n++;
                tl_trace = tl_trace + 100;
                }  //  while loop                                                                                                                                                       
            }  // else d01 < d02      

    }  // if Stopping
    else  // block for through-going muons 
        {
         // cout<<"This is a Through-going muons events :"<<"Trajectory Track length = \t"<<fTrackLength.back()<<"\t Reconstructed track length =\t"<<Extpl_track_length1<<"=tl_trace >\t"<<tl_trace<<"\n"; 
    
            fNon_Stp_Extpltl.push_back(Extpl_track_length1);
            fNon_Stp_Trjtl.push_back(fTrackLength.back());
  }  this braket not needed for next steps       
    if(d01 < d02 )
      {   
          while( tl_trace < Extpl_track_length1)  //  have issue with track smaller then 500 cm
            { //  SP_TracePoint_gr->SetPoint(n,P0.X(),P0.Y(),P0.Z());
              //  cout<< "Point at the line (Forward Dir ) is : (\t"<<P0(0)<<",\t"<<P0(1)<<",\t"<<P0(2)<<"\t"<<"distance = \t"<<tl_trace<<"\n";
               p1.SetXYZ(P0.X(),P0.Z(),0);
               p2.SetXYZ((P0+tl*V).X(),(P0+tl*V).Z(),0);

           for( auto const& hit : (*hitHandle) )
             {
               if ( fGeometryService->SignalType( hit.Channel() ) != geo::kCollection) continue;
                     p.SetXYZ(detProp.ConvertTicksToX(hit.PeakTime(), hit.WireID().Plane, hit.WireID().TPC, hit.WireID().Cryostat ), geom->Wire(hit.WireID()).GetCenter().Z(), 0 );

             if(((p-p1).Dot(p2-p1) > 0.0) && ((p-p2).Dot(p1-p2) > 0.0)) { dqdl=dqdl+ hit.Integral(); dqdlADC=dqdlADC + hit.SummedADC();} // here * is Dot product 
            // if(p1.Dot(p2)){cout<<"Dot product syntex "<<"\n";}      
             }
           cout<< "Totale Charge from Track Segment : "<<dqdl/tl<<"\t ADC charge =>\t"<<dqdlADC/tl<<"\n";
       //  break;
       if(dqdl!=0) { fHitdQ_track_segment.push_back(dqdl/tl);}
            P0 = P0 + tl*V; 
            dqdl=0.0;
            dqdlADC=0.0;
            n++;
            tl_trace = tl_trace + 100;
     }  //  while loop       
   }  // if (d01 or d02)
   else // (d01 not less than d02)
      {
     while( tl_trace < Extpl_track_length1)
      {  // SP_TracePoint_gr->SetPoint(n,P0.X(),P0.Y(),P0.Z());
         //  cout<< "Point at the line (Reverse dir ) is : (\t"<<P0(0)<<",\t"<<P0(1)<<",\t"<<P0(2)<<")"<<"distance = \t"<<tl_trace<<"\n";
       
            p1.SetXYZ(P0.X(),P0.Z(),0);
            p2.SetXYZ((P0-tl*V).X(),(P0-tl*V).Z(),0);     
              for( auto const& hit : (*hitHandle) )
               {  
                   if (fGeometryService->SignalType( hit.Channel() ) != geo::kCollection)  continue;
              p.SetXYZ(detProp.ConvertTicksToX(hit.PeakTime(), hit.WireID().Plane, hit.WireID().TPC, hit.WireID().Cryostat ), geom->Wire(hit.WireID()).GetCenter().Z(),0 );
                 if(((p-p1).Dot(p2-p1) > 0.0) && ((p-p2).Dot(p1-p2) > 0.0)) { dqdl=dqdl+ hit.Integral(); dqdlADC=dqdlADC + hit.SummedADC();}
                 }
               cout<< "Totale Charge from Track Segment : "<<dqdl/tl<< "\tADC Charge => \t"<<dqdlADC/tl<<"\n";
             //  dqdl=dqdl+HitCharge;                                                                    
             //  break;           
             if(dqdl!=0)  { fHitdQ_track_segment.push_back(dqdl/tl);}                                                                                              
                 P0 = P0 - tl*V;                                                                                     
                 dqdl=0.0;
                 dqdlADC=0.0;                                          
                 n++;    
                 tl_trace = tl_trace + 100;
              }  //  while loop                                                                                                                                                       
             } // else //  //   End estimation dQ/dl here ...
    } /// else  through going muons events 
// END Stopping Muons Analysis $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 
*/  
    //    if(Extpl_track_length1 >= 1200 && trackLength <= 500 )
     //   {
  /*         TCanvas *tc1 = new TCanvas("tc1","Rec Space-Point",0,20,800,800);
            TH2D *h = new TH2D("h","Space",100,-720,720,100,-600,600);
           // TH2D *h1 = new TH2D("h1","Space",100,-720,720,100,-600,600);
             gStyle->SetOptStat(0);
             h->GetZaxis()->SetLimits(0.0,6500);
             SP_gr->SetHistogram(h); 
             SP_gr->Draw("p0");
             SP_gr->SetTitle("DUNE-FD ; X (cm); Y (cm); Z (cm)");
             SP_gr->GetXaxis()->CenterTitle();SP_gr->GetYaxis()->CenterTitle();SP_gr->GetZaxis()->CenterTitle();
             SP_gr->GetXaxis()->SetTitleOffset(1.5);SP_gr->GetYaxis()->SetTitleOffset(1.5);SP_gr->GetZaxis()->SetTitleOffset(1.3);
          //   SP_gr->SetMinimum(0.0);SP_gr->SetMaximum(6500);
             SP_gr->GetXaxis()->SetTitleSize(0.04);SP_gr->GetYaxis()->SetTitleSize(0.04);SP_gr->GetZaxis()->SetTitleSize(0.04);         
             SP_gr->GetXaxis()->SetLabelSize(0.025); SP_gr->GetYaxis()->SetLabelSize(0.025); SP_gr->GetZaxis()->SetLabelSize(0.025);            
              
     
         
          // Traj_gr->SetMarkerColor(kGreen);
          //  Traj_gr->SetMarkerStyle(10);
          //  Traj_gr->SetMarkerSize(0.5);
          //  Traj_gr->Draw("p0 same");
          

          
        //    SP_AvePos_gr->SetPoint(0,pos[0],pos[1],pos[2]);
        //   SP_AvePos_gr->SetPoint(1,fRTrackExitXpos[0],fRTrackExitYpos[0],fRTrackExitZpos[0]);
        //    SP_AvePos_gr->SetPoint(2,fRTrackExitXpos[1],fRTrackExitYpos[1],fRTrackExitZpos[1]);
            SP_AvePos_gr->SetMarkerColor(kRed);
            SP_AvePos_gr->SetMarkerStyle(20);
            SP_AvePos_gr->SetMarkerSize(1.5);
            SP_AvePos_gr->Draw("pLINE same");
            SP_AvePos_gr->SetLineColor(kRed);
       

            // gStyle->SetMinimumZero(kFALSE); 
       
          SP_TracePoint_gr->SetMarkerColor(kBlue);
             SP_TracePoint_gr->SetMarkerStyle(20);
             SP_TracePoint_gr->SetMarkerSize(1.5);
             SP_TracePoint_gr->Draw("p same");
          //   h1->GetZaxis()->SetRangeUser(0.0,6500);
             SP_TracePoint_gr->SetMinimum(0.0);
             SP_TracePoint_gr->SetMaximum(6500);
      
          SPLPd_gr->SetMarkerColor(kYellow); 
          SPLPd_gr->SetMarkerStyle(20);
          SPLPd_gr->SetMarkerSize(1.5);
          SPLPd_gr->Draw("p same");
          SPLPd_gr->SetMinimum(0.0);
          SPLPd_gr->SetMaximum(6500);
         
           TLegend *leg = new TLegend(0.15,0.9,0.35,0.7);
            // leg->SetHeader("Track Length");
            leg->AddEntry(SP_AvePos_gr,"Reco(Line)","l");
         //  leg->AddEntry(SP_TracePoint_gr,"Sphere Centers","p");
          //  leg->AddEntry(Traj_gr,"Traj","p");
             leg->AddEntry(SP_gr,"Space Points","p0");
           //  leg->AddEntry(SPLPd_gr,"Space Points at d = 5 cm from line","p0");
             leg->Draw();
        
           //  tc1->Print("SpacePoint_location.root");
            tc1->SaveAs("/dune/app/users/jdsingh/DUNE_SP2020Work/DUNEWork/MUSUNWork/myplot/SpacePointMC.pdf");
     //    }// if(TGraph2D) 
  //  else continue;
   //   }//stopping  else 
  
*/
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
          MF_LOG_DEBUG("EnergyAna") << "Hit in collection plane"<<"JD \t "<<hitCharge<< std::endl;
           hitCharge += (hit.Integral()); 

     MF_LOG_DEBUG("EnergyAna")
              << std::endl;
                 fRecoChargeInt.push_back(hit.Integral());  // We should work with this quantity , 
              //   fRecoChargeTrue.push_back(hit.PeakAmplitude());  // this have noise .....
                 fRecoWireNumber.push_back(hit.Channel());
                 fRecoPeakTime.push_back(hit.PeakTime());
        
     MF_LOG_DEBUG("EnergyAna")
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
 //  fRecoChargeTrue.push_back(hitCharge*fElectronsToGeV);
// ***************************************     Charge and Energy Deposite Analysis **************************// 
 // }  // if stopping final...............

    fSimulationTrackTree->Fill();
 
   } // if(SpacePoint)
} // EnergyAna::analyze()

DEFINE_ART_MODULE(trackAna)

} // namespace example
} // namespace lar



