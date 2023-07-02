/**
 * @file   hitExample_module.cc
 * @brief  A basic "skeleton" to read in art::Event records from a file,
*/

// LArSoft includes
#include "lardata/RecoBaseProxy/ChargedSpacePoints.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PointCharge.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/Exception.h"

// Utility libraries
#include "cetlib/pow.h" // cet::sum_of_squares()
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT includes. Note: To look up the properties of the ROOT classes,
// use the ROOT web site; e.g.,
// <https://root.cern.ch/doc/master/annotated.html>
#include "TH1.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TVector3.h"

// C++ includes
#include <cmath>
#include <map>
using namespace std;

namespace {


  

} // local namespace



namespace lar {
namespace example {


  class JDSpacePointAna : public art::EDAnalyzer
  {
  public:

 
    struct Config {

      // Save some typing:
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      // One Atom for each parameter
      fhicl::Atom<art::InputTag> SimulationLabel {
        Name("SimulationLabel"),
        Comment("tag of the input data product with the detector simulation information")
        };

    /*  fhicl::Atom<art::InputTag> TrackLabel {
        Name("TrackLabel"),
        Comment("tag of the input data product with reconstructed hits")
        };*/

         fhicl::Atom<art::InputTag> JDSpacePointLabel {
        Name("JDSpacePointLabel"),
        Comment("tag of the input data product with reconstructed Space Point ")
        };

     

      fhicl::Atom<int> PDGcode {
        Name("PDGcode"),
        Comment("particle type (PDG ID) of the primary particle to be selected")
        };

      fhicl::Atom<double> BinSize {
        Name("BinSize"),
        Comment("dx [cm] used for the dE/dx calculation")
        };

    }; // Config

    
    using Parameters = art::EDAnalyzer::Table<Config>;

   

    /// Constructor: configures the module (see the Config structure above)
    explicit JDSpacePointAna(Parameters const& config);

    virtual void beginJob() override;

    // This method is called once, at the start of each run. It's a
    // good place to read databases or files that may have
    // run-dependent information.
    virtual void beginRun(const art::Run& run) override;

    // The analysis routine, called once per event.
    virtual void analyze (const art::Event& event) override;

  private:

    // The stuff below is the part you'll most likely have to change to
    // go from this custom example to your own task.

    // The parameters we'll read from the .fcl file.
    art::InputTag fSimulationProducerLabel; ///< The name of the producer that tracked simulated particles through the detector
    art::InputTag fTrackProducerLabel;        ///< The name of the producer that created hits
    art::InputTag fSpacePointProducerLabel;        ///< The name of the producer that created hits
    art::InputTag fTrkSpptAssocModuleLabel;  
 art::InputTag fClusterProducerLabel;    ///< The name of the producer that created clusters
    int fSelectedPDG;                       ///< PDG code of particle we'll focus on
    double fBinSize;                        ///< For dE/dx work: the value of dx.

    // Pointers to the histograms we'll create.
    TH1D* fPDGCodeHist;     ///< PDG code of all particles
    TH1D* fMomentumHist;    ///< momentum [GeV] of all selected particles
     TTree* fSimulationNtuple;     ///< tuple with simulated data
    TTree* fReconstructionNtuple; ///< tuple with reconstructed data
      

     int fEvent;        ///< number of the event being processed
    int fRun;          ///< number of the run being processed
    int fSubRun;       ///< number of the sub-run being processed
    /// @}

       // Geometry Info
      geo::GeometryCore const *fGeometry;    ///< pointer to geometery provider
                     ///< Conversation factor


    /// @name The variables that will go into the simulation n-tuple.
    /// @{
    int fSimPDG;       ///< PDG ID of the particle being processed
    int fSimTrackID;   ///< GEANT ID of the particle being processed
    int tmx, tmy, tmz, tpx, tpy, tpz;
    bool cmx, dmx, cmy, dmy, cmz, dmz; 
    bool cpx, dpx, cpy, dpy, cpz, dpz; 
  //  double fRecSapcePointX;    //  Reconstructed Space point X position 

    
             double lineXM = -720;    double lineYM = -600;     double lineZM = 0;
             double lineXP = +720;    double lineYP = +600;     double lineZP = 6500;    /// Detector planes  position on XYZ axis 
     double RTrackExitXposmx; 
    double RTrackExitYposmy;
     double RTrackExitZposmz;
     double RTrackExitXpospx;
     double RTrackExitYpospy;
     double RTrackExitZpospz;
     double MCTrackExitXposmx;
     double MCTrackExitYposmy;
     double MCTrackExitZposmz;
     double MCTrackExitXpospx;
     double MCTrackExitYpospy;
     double MCTrackExitZpospz;
     double   IPl_track_length; 
     double   IPl_MCtrack_length;
    // Arrays for 4-vectors: (x,y,z,t) and (Px,Py,Pz,E).
    // Note: old-style C++ arrays are considered obsolete. However,
    // to create simple n-tuples, we still need to use them.
    double fStartXYZT[4]; ///< (x,y,z,t) of the true start of the particle
    double fEndXYZT[4];   ///< (x,y,z,t) of the true end of the particle
    double fStartPE[4];   ///< (Px,Py,Pz,E) at the true start of the particle
    double fEndPE[4];     ///< (Px,Py,Pz,E) at the true end of the particle
    

    std::vector<double >fTrackExitXpos;
    std::vector<double >fTrackExitYpos;
    std::vector<double >fTrackExitZpos;


        int fSimNdEdxBins;                      ///< dE/dx bins in a track   
      std::vector<double> fSimWireEnergy;  ///< vector for dE/dx values 
      std::vector<double>  fTrackLength;
      std::vector<double>  fExtplTrackLength;
      std::vector<int>    fSimWireNumber;
      std::vector<double> fRec_SpacePoint_X;  
      std::vector<double> fRec_SpacePoint_Y;
      std::vector<double> fRec_SpacePoint_Z;
       std::vector<double> fPointCharge;
      std::vector<double> fPointCharge5m;
      std::vector<double> fPointCharge10m;
      std::vector<double> fPointCharge20m;
   
    // Other variables that will be shared between different methods.
    geo::GeometryCore const* fGeometryService;   ///< pointer to Geometry provider
    detinfo::DetectorClocks const* fTimeService; ///< pointer to detector clock time service provider
    double                   fElectronsToGeV;    ///< conversion factor
    int                      fTriggerOffset;     ///< (units of ticks) time of expected neutrino event

  }; // class hitExample

  /// @}
  // END hitExample group -------------------------------------------------


  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class implementation

  //-----------------------------------------------------------------------
  // Constructor
  //
  // Note that config is a Table<Config>, and to access the Config
  // value we need to use an operator: "config()". In the same way,
  // each element in Config is an Atom<Type>, so to access the type we
  // again use the call operator, e.g. "SimulationLabel()".
  //
 JDSpacePointAna::JDSpacePointAna(Parameters const& config)
    : EDAnalyzer(config)
    , fSimulationProducerLabel(config().SimulationLabel())
   // , fTrackProducerLabel       (config().TrackLabel())
    , fSpacePointProducerLabel      (config().JDSpacePointLabel())
   // , fTrkSpptAssocModuleLabel      (config().TrkSpptLabel()) //(pset.get<std::string>("TrkSpptAssocModuleLabel"))
   // , fClusterProducerLabel   (config().ClusterLabel())
    , fSelectedPDG            (config().PDGcode())
    , fBinSize                (config().BinSize())
  {
    // Get a pointer to the geometry service provider.
    fGeometryService = lar::providerFrom<geo::Geometry>();
    // The same for detector TDC clock services.
    fTimeService = lar::providerFrom<detinfo::DetectorClocksService>();
    // Access to detector properties.
    const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    fTriggerOffset = detprop->TriggerOffset();

    // Since art 2.8, you can and should tell beforehand, here in the constructor,
    // all the data the module is going to read ("consumes") or might read
    // ("may_consume"). Diligence here will in the future help the framework
    // execute modules in parallel, making sure the order is correct.
    consumes<std::vector<simb::MCParticle>>(fSimulationProducerLabel);
    consumes<std::vector<sim::SimChannel>>(fSimulationProducerLabel);
    consumes<art::Assns<simb::MCTruth, simb::MCParticle>>(fSimulationProducerLabel);
    //consumes<std::vector<recob::Track>>(fTrackProducerLabel);
     consumes<std::vector<recob::SpacePoint>>(fSpacePointProducerLabel);   //  Space point producer class 
    //consumes<std::vector<recob::Cluster>>(fClusterProducerLabel);
   // consumes<art::Assns<recob::Cluster, recob::Track>>(fTrackProducerLabel);

  }


  //-----------------------------------------------------------------------
  void JDSpacePointAna::beginJob()
  {
    // Get the detector length, to determine the maximum bin edge of one
    // of the histograms.
  //  const double detectorLength = DetectorDiagonal(*fGeometryService);

    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us.
    art::ServiceHandle<art::TFileService const> tfs;

  
    fPDGCodeHist     = tfs->make<TH1D>("pdgcodes",";PDG Code;",                  5000, -2500, 2500);
    fMomentumHist    = tfs->make<TH1D>("MC_momentum",     ";particle Momentum (GeV);",    100, 0.,    10.);
   

    // Define our n-tuples, which are limited forms of ROOT
    // TTrees. Start with the TTree itself.
     fSimulationNtuple     = tfs->make<TTree>("MCAnalysis",    "EnergyOnWireTree");
     fReconstructionNtuple     = tfs->make<TTree>("RECAnalysis",    "Rec_SapcepointSimulation");
    

    // Define the branches (columns) of our simulation n-tuple. To
    // write a variable, we give the address of the variable to
    // TTree::Branch.
    fSimulationNtuple->Branch("Event",       &fEvent,          "Event/I");
    fSimulationNtuple->Branch("SubRun",      &fSubRun,         "SubRun/I");
    fSimulationNtuple->Branch("Run",         &fRun,            "Run/I");
    fSimulationNtuple->Branch("TrackID",     &fSimTrackID,     "TrackID/I");
    fSimulationNtuple->Branch("PDG",         &fSimPDG,         "PDG/I");
    fSimulationNtuple->Branch("StartXYZT",   fStartXYZT,       "StartXYZT[4]/D");
    fSimulationNtuple->Branch("EndXYZT",     fEndXYZT,         "EndXYZT[4]/D");
    fSimulationNtuple->Branch("StartPE",     fStartPE,         "StartPE[4]/D");
    fSimulationNtuple->Branch("EndPE",       fEndPE,           "EndPE[4]/D");
    fSimulationNtuple->Branch("NdEdx",       &fSimNdEdxBins,   "NdEdx/I");
    fSimulationNtuple->Branch("TrackLength",    &fTrackLength);
    fSimulationNtuple->Branch("ExtplTrackLength",    &fExtplTrackLength);
     fSimulationNtuple->Branch("TrackExitXpos",    &fTrackExitXpos);
    fSimulationNtuple->Branch("TrackExitYpos",    &fTrackExitYpos);
    fSimulationNtuple->Branch("TrackExitZpos",    &fTrackExitZpos);
    
  
    // fReconstructionNtuple->Branch("Event",   	     &fEvent,             "Event/I");
   
 
       fReconstructionNtuple->Branch("SpacePoint_X", &fRec_SpacePoint_X);
       fReconstructionNtuple->Branch("SpacePoint_Y", &fRec_SpacePoint_Y);
       fReconstructionNtuple->Branch("SpacePoint_Z", &fRec_SpacePoint_Z);
       fReconstructionNtuple->Branch("pointq", &fPointCharge); 
       fReconstructionNtuple->Branch("pointq5m", &fPointCharge5m); 
       fReconstructionNtuple->Branch("pointq10m", &fPointCharge10m);
       fReconstructionNtuple->Branch("pointq20m", &fPointCharge20m);      
                                                         
  }

  //-----------------------------------------------------------------------
  // art expects this function to have a art::Run argument; C++
  // expects us to use all the arguments we are given, or it will
  // generate an "unused variable" warning. But we don't actually need
  // nor use the art::Run object in this example. The trick to prevent
  // that warning is to omit (or comment out) the name of the
  // parameter.

  void JDSpacePointAna::beginRun(const art::Run& /*run*/)
  {
    
    art::ServiceHandle<sim::LArG4Parameters const> larParameters;
    fElectronsToGeV = 1./larParameters->GeVToElectrons();
  }
 

 
 void JDSpacePointAna::analyze(const art::Event& event)
  {

      fEvent  = event.id().event();
      fRun    = event.run();
     fSubRun = event.subRun(); 

    art::ServiceHandle<geo::Geometry> geom;
   //auto const *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();;

    art::Handle< std::vector<simb::MCParticle> > particleHandle;
    // auto simChannelHandle = event.getValidHandle<std::vector<sim::SimChannel>>(fSimulationProducerLabel);

   if (!event.getByLabel(fSimulationProducerLabel, particleHandle))
      {
         throw cet::exception("hitExample")
	  << " No simb::MCParticle objects in this event - "
	  << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
      }

     fTrackLength.clear(); 
     fExtplTrackLength.clear();
     fTrackExitXpos.clear();
     fTrackExitYpos.clear();
     fTrackExitZpos.clear();

    std::map< int, const simb::MCParticle* > particleMap;

     for ( auto const& particle : (*particleHandle) )
      {
          
          fSimTrackID = particle.TrackId();
          
          particleMap[fSimTrackID] = &particle;
       
        fSimPDG = particle.PdgCode();
	fPDGCodeHist->Fill( fSimPDG );

        const size_t numberTrajectoryPoints = particle.NumberTrajectoryPoints();

          const int last = numberTrajectoryPoints - 1;
	const TLorentzVector& positionStart = particle.Position(0);
	const TLorentzVector& positionEnd   = particle.Position(last);
	const TLorentzVector& momentumStart = particle.Momentum(0);
	const TLorentzVector& momentumEnd   = particle.Momentum(last);

         fMomentumHist->Fill( momentumStart.P() );
     
         positionStart.GetXYZT( fStartXYZT );
	positionEnd.GetXYZT( fEndXYZT );
	momentumStart.GetXYZT( fStartPE );
	momentumEnd.GetXYZT( fEndPE );

       // TVector3 endMC = positionEnd.End<TVector3>();
        TVector3 endMC = positionEnd.Vect();
        TVector3 unitMC =  endMC.Unit();

     // TVector3 pos = particle.Position(numberTrajectoryPoints).Vect();  // Let the particle go inside the detector 
     // double const tmpArray[]={pos.X(),pos.Y(),pos.Z()};
      
     // geo::TPCID tpcid=geom->FindTPCAtPosition(tmpArray);
    // if (!tpcid.isValid) continue;
          const double trackLength = ( positionEnd - positionStart ).Rho();
          fTrackLength.push_back(trackLength);
         


 


 tmx = (lineXM - positionEnd.X())/(unitMC.X());  tmy  = (lineYM- positionEnd.Y())/(unitMC.Y());   tmz  = (lineZM- positionEnd.Z())/(unitMC.Z());
 tpx = (lineXP - positionEnd.X())/(unitMC.X());  tpy = (lineYP- positionEnd.Y())/(unitMC.Y());   tpz  = (lineZP- positionEnd.Z())/(unitMC.Z());
             
 double MCTrackExitYposmx=positionEnd.Y() + tmx;double  MCTrackExitXposmy=positionEnd.X()+ tmy; double  MCTrackExitXposmz=positionEnd.X()+ tmz; 
double  MCTrackExitZposmx=positionEnd.Z() +tmx;double MCTrackExitZposmy=positionEnd.Z() + tmy; double   MCTrackExitYposmz=positionEnd.Y() + tmz;


double  MCTrackExitYpospx=positionEnd.Y() + tpx; double  MCTrackExitXpospy=positionEnd.X() + tpy;  double MCTrackExitXpospz=positionEnd.X()+tpz;
double  MCTrackExitZpospx=positionEnd.Z() + tpx; double  MCTrackExitZpospy=positionEnd.Z() + tpy; double MCTrackExitYpospz=positionEnd.Y()+tpz;

    if((MCTrackExitYposmx >= -600. && MCTrackExitYposmx < 600.) && (MCTrackExitZposmx >= 0. && MCTrackExitZposmx < 6500.) )
          {  
                  MCTrackExitXposmx  = lineXM;
         
         fTrackExitXpos.push_back( MCTrackExitXposmx );
         fTrackExitYpos.push_back( MCTrackExitYposmx );
         fTrackExitZpos.push_back( MCTrackExitZposmx ); 

     IPl_MCtrack_length =  sqrt(pow(MCTrackExitXposmx -  positionEnd.X(), 2) + pow(MCTrackExitYposmx -  positionEnd.Y(), 2) 
                                 + pow(MCTrackExitZposmx -  positionEnd.Z(), 2) * 1.0);  

    
    fExtplTrackLength.push_back(IPl_MCtrack_length);
        }
            else 
            if((MCTrackExitXposmy > -720. && MCTrackExitXposmy < 720.) && (MCTrackExitZposmy > 0. && MCTrackExitZposmy < 6500.))
               {   
                      MCTrackExitYposmy  = lineYM;
                      
                        fTrackExitXpos.push_back( MCTrackExitXposmx );
                        fTrackExitYpos.push_back( MCTrackExitYposmx );
                        fTrackExitZpos.push_back( MCTrackExitZposmx ); 

              IPl_MCtrack_length =  sqrt(pow(MCTrackExitXposmy -  positionEnd.X(), 2) + pow(MCTrackExitYposmy -  positionEnd.Y(), 2) 
                                 + pow(MCTrackExitZposmy -  positionEnd.Z(), 2) * 1.0);  

   
    fExtplTrackLength.push_back(IPl_MCtrack_length);           
          } 
          else
           if((MCTrackExitXposmz >= -720. && MCTrackExitXposmz < 720.) && (MCTrackExitYposmz >= -600. && MCTrackExitYposmz < 600.))
               {     MCTrackExitZposmz  = lineZM;
                      
                        fTrackExitXpos.push_back( MCTrackExitXposmx );
        		 fTrackExitYpos.push_back( MCTrackExitYposmx );
         		fTrackExitZpos.push_back( MCTrackExitZposmx );

              IPl_MCtrack_length =  sqrt(pow(MCTrackExitXposmz -  positionEnd.X(), 2) + pow(MCTrackExitYposmz -  positionEnd.Y(), 2) 
                                 + pow(MCTrackExitZposmz -  positionEnd.Z(), 2) * 1.0);  

    
    fExtplTrackLength.push_back(IPl_MCtrack_length);
      }  
            else 
              if((MCTrackExitYpospx >= -600. && MCTrackExitYpospx < 600.) && (MCTrackExitZpospx >= 0. && MCTrackExitZpospx < 6500.))
                   {   
                      MCTrackExitXpospx  = lineXP;
                     
                        fTrackExitXpos.push_back( MCTrackExitXposmx );
       		       fTrackExitYpos.push_back( MCTrackExitYposmx );
                      fTrackExitZpos.push_back( MCTrackExitZposmx ); 

                  IPl_MCtrack_length =  sqrt(pow(MCTrackExitXpospx -  positionEnd.X(), 2) + pow(MCTrackExitYpospx -  positionEnd.Y(), 2) 
                                 + pow(MCTrackExitZpospx -  positionEnd.Z(), 2) * 1.0);  

            
             fExtplTrackLength.push_back(IPl_MCtrack_length);
                }  
           else  
           if((MCTrackExitXpospy >= -720. && MCTrackExitXpospy < 720.) && (MCTrackExitZpospy >= 0. && MCTrackExitZpospy < 6500.))
                {    MCTrackExitYpospy  = lineYP;
                      
                        fTrackExitXpos.push_back( MCTrackExitXposmx );
		         fTrackExitYpos.push_back( MCTrackExitYposmx );
		         fTrackExitZpos.push_back( MCTrackExitZposmx );

                  IPl_MCtrack_length =  sqrt(pow(MCTrackExitXpospy -  positionEnd.X(), 2) + pow(MCTrackExitYpospy -  positionEnd.Y(), 2) 
                                 + pow(MCTrackExitZpospy -  positionEnd.Z(), 2) * 1.0);  

             
             fExtplTrackLength.push_back(IPl_MCtrack_length);                                                        
                } 
           else
           if((MCTrackExitXpospz >= -720. && MCTrackExitXpospz < 720.) && (MCTrackExitYpospz >= -600. && MCTrackExitYpospz < 600.) )
                {  MCTrackExitZpospz  = lineZP;
                      
  		     fTrackExitXpos.push_back( MCTrackExitXposmx );
         	     fTrackExitYpos.push_back( MCTrackExitYposmx );
                     fTrackExitZpos.push_back( MCTrackExitZposmx );
                  IPl_MCtrack_length =  sqrt(pow(MCTrackExitXpospz -  positionEnd.X(), 2) + pow(MCTrackExitYpospz -  positionEnd.Y(), 2) 
                                 + pow(MCTrackExitZpospz -  positionEnd.Z(), 2) * 1.0);  

            
             fExtplTrackLength.push_back(IPl_MCtrack_length);
          } 
            else { continue;
                   //std::cout << "No, Point is out of the XY Plane   =====" << std::endl;
                   } 
          

      MF_LOG_DEBUG("hitExample")
	  << "Track length: " << trackLength << " cm";

	// Fill a histogram of the track length.
	
	MF_LOG_DEBUG("hitExample")
	  << "track ID=" << fSimTrackID
	  << " (PDG ID: " << fSimPDG << ") "
	  << trackLength << " cm long, momentum "
	  << momentumStart.P() << " GeV/c, has "
	  << numberTrajectoryPoints << " trajectory points";

        

  } // loop over all particles in the event.

     fSimulationNtuple->Fill(); 



     // ##########   Reconstructed Space Point Estimation ###############

            art::Handle< std::vector<recob::SpacePoint> > recobspacepoints;
            event.getByLabel(fSpacePointProducerLabel, recobspacepoints);

             art::Handle< std::vector<recob::PointCharge> > recobCharge;
            event.getByLabel(fSpacePointProducerLabel, recobCharge);

           // auto qHandle = event.getValidHandle< std::vector<recob::PointCharge> >(fSpacePointProducerLabel);
          //  auto const& recobspacepoints = *event.getValidHandle<vector<recob::SpacePoint>>();
       //    auto points = proxy::getChargedSpacePoints(event, recobspacepoints);
            fRec_SpacePoint_X.clear();
            fRec_SpacePoint_Y.clear();
            fRec_SpacePoint_Z.clear();
            fPointCharge.clear();
            fPointCharge5m.clear();
            fPointCharge10m.clear();
            fPointCharge20m.clear();
              
            
             if (!recobspacepoints->empty())
            {
	        
	         double xpt[3] = {-15,400,0};
	         double dcos[3] = {-0.2,-0.15,0.9899};
            
	   // TGraph2D *gri = new TGraph2D();
	  //  TGraph2D *gro = new TGraph2D();
          //  TCanvas *c = new TCanvas("c","TGraph2D Event Display",0,0,800,800);
	    for (size_t isp=0;isp<recobspacepoints->size(); ++isp)
	      {

                  art::Ptr<recob::SpacePoint> recSPoint(recobspacepoints, isp);
                 const recob::SpacePoint& RecSP = *recSPoint;

                 art::Ptr<recob::PointCharge> recCharge(recobCharge, isp);
                 const recob::PointCharge& RecQ = *recCharge;
                
                
		double x = RecSP.XYZ()[0];
		double y = RecSP.XYZ()[1];
		double z = RecSP.XYZ()[2];
                double Q = RecQ.charge();
		//cout << "Space point: " << x << " " << y << " " << z << endl;

		double xd = (x-xpt[0]);
		double yd = (y-xpt[1]);
		double zd = (z-xpt[2]);

		double xc = (yd * dcos[2] - dcos[1]*zd); // signs might be wrong here but we'll square it anyway
		double yc = (xd * dcos[2] - dcos[0]*zd);
		double zc = (xd * dcos[1] - dcos[0]*yd);
		double sdiff = TMath::Sqrt(xc*xc+yc*yc+zc*zc);

		if (sdiff < 0.0) // 0 cm cut from the road
		  {
		    //gri->SetPoint(inroad,z,x,y);
                     cout << "Got here" << endl;
                   // fRecSpacePointXHist->Fill( x );
		   // ++inroad;
		  }
		else
		 {
               
                       fRec_SpacePoint_X.push_back(x);
                       fRec_SpacePoint_Y.push_back(y);
                       fRec_SpacePoint_Z.push_back(z);
                       fPointCharge.push_back(Q);
                       if(z <= 500.00){fPointCharge5m.push_back(Q);}
                       if(z <= 1000.00){fPointCharge10m.push_back(Q);}
                       if(z <= 2000.00){fPointCharge20m.push_back(Q);}
                       
                  
                   
                   
		  }
                   
	      }
              
                
          }
           
  

    fReconstructionNtuple->Fill();  
   
   //  reconstrution analysis  @JAYDIP SINGH  

} // hitExample::analyze()


  // This macro has to be defined for this module to be invoked from a
  // .fcl file; see hitExample.fcl for more information.
  DEFINE_ART_MODULE(JDSpacePointAna)

} // namespace example
} // namespace lar
