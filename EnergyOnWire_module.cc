/* This is a module intended to imitate the AnalysisPracticeAM_moudle.cc
 * There was a problem calculating RMS for energy per wire collected in 
 * that edition, and this is an attempt to verify the results from that 
 * study
 *
 * Everything is copied over, except the names of the module and analysis 
 * are change, also instead of calculating the path length of the module 
 * and dividing by it to the get the energy loss quantity dE/dx
 * --------------------------------
 * Kevin Ingles
 * 20 December 2017
 * EnergyOnWire_module.cc
 * --------------------------------
 */
 
// LArSoft
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// ART Framework
#include "canvas/Utilities/Exception.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"

// Utilities
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib/pow.h"

// ROOT 
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"

// STL 
#include <map>
#include <vector>
#include <string>
#include <cmath>

namespace
{
  double DetectorDiagonal(geo::GeometryCore const& geom);
 
 // bool TDCIDETimeCompare( const sim::TDCIDE&, const sim::TDCIDE& );
}

namespace lar
{
  namespace wireEnergy
  {
    class EnergyOnWire : public art::EDAnalyzer
    {
    public:
      struct Config ///<FHICL Parameter definitions
      {
	using Name = fhicl::Name;
	using Comment = fhicl::Comment;
	using Atom = fhicl::Atom<art::InputTag>;

	Atom SimulationLabel
	{
	  Name("SimulationLabel"),
	    Comment("Tag of the input data product with the detector simulation")
	    };

        fhicl::Atom<art::InputTag> HitLabel {
        Name("HitLabel"),
        Comment("tag of the input data product with reconstructed hits")
        };

      fhicl::Atom<art::InputTag> ClusterLabel {
        Name("ClusterLabel"),
        Comment("tag of the input data product with reconstructed clusters")
        };

	fhicl::Atom<int> PDGCode
	{
	  Name("PDGCode"),
	    Comment("Particle type(ID) of the primary particle to be selected")
	    };
	fhicl::Atom<double> BinSize
	{
	  Name("BinSize"),
	    Comment("dx [cm] used for the dE/dc calculation")
	    };
      };
      // End Config
      // Public Member Functions
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit EnergyOnWire(Parameters const& config);    ///< Constructor

      virtual void beginJob() override;
      virtual void beginRun(const art::Run& run) override;
      virtual void analyze (const art::Event& event) override;

    private:
      // Parameters from .fcl file
      art::InputTag fSimulationProducerLabel; ///< Name of producer that track
       art::InputTag fHitProducerLabel;        ///< The name of the producer that created hits
    art::InputTag fClusterProducerLabel;
      int fSelectedPDG;                       ///< Particle I
      double fBinSize;                        ///< dx value for dE/dx

      // Pointers to histograms created
      TH1D *fPDGCodeHist;
      TH1D *fMomentumHist;     ///< in GeV
      TH1D *fTrackLengthHist;  ///< in cm, true length

      // N-tuples created
      TTree *fSimulationNtuple;
      TTree* fReconstructionNtuple;
           
      // Variables going into n-tuples
      int fEvent;   ///< number of events being processed
      int fRun;     ///< ''     '' runs   ''    ''
      int fSubRun;  ///< ''     '' subruns''    ''                            
      // -----------------------------------------
      //Simulation Tree Info
      int fSimPDG;      ///< particle ID
      int fSimTrackID;  ///< GEANT track ID

      // 4-vector arrays for data
      double fStartXYZT[4];  ///< (x,y,z,t) position start
      double fEndXYZT[4];    ///< position end
      double fStartPE[4];    ///< (Px,Py,Pz,E) at start
      double fEndPE[4];      ///< at end

      int fSimNdEdxBins;                      ///< dE/dx bins in a track   
      std::vector<double> fSimWireEnergy;  ///< vector for dE/dx values 
      std::vector<int>    fSimWireNumber;


        int fRecoPDG;       ///< PDG ID of the particle being processed
    int fRecoTrackID;   ///< GEANT ID of the particle being processed

    /// Number of dE/dx bins in a given track.
    int fRecoNdEdxBins;

    /// The vector that will be used to accumulate dE/dx values as a function of range.
    std::vector<double> fRecodEdxBins;
	  
      // Geometry Info
      geo::GeometryCore const *fGeometry;
      geo::GeometryCore const* fGeometryService;   ///< pointer to Geometry provider
    detinfo::DetectorClocks const* fTimeService;    ///< pointer to geometery provider
      double fElectronsToGeV;                ///< Conversation factor
       int                      fTriggerOffset; 

      // Optional, output .txt file giving particle number, initial and final 
      // position and momentum and range
      // ------------------------------------------
    }; // End EnergyOnWire
    // Class Implementation
    EnergyOnWire::EnergyOnWire(Parameters const& config)
      : EDAnalyzer(config)
      , fSimulationProducerLabel (config().SimulationLabel())
      , fHitProducerLabel       (config().HitLabel())
      , fClusterProducerLabel   (config().ClusterLabel())
      , fSelectedPDG             (config().PDGCode())
      , fBinSize                 (config().BinSize())
    {
      fGeometry = lar::providerFrom<geo::Geometry>();

        fTimeService = lar::providerFrom<detinfo::DetectorClocksService>();
    // Access to detector properties.
    const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    fTriggerOffset = detprop->TriggerOffset();    


      consumes<std::vector<simb::MCParticle>>(fSimulationProducerLabel);
    consumes<std::vector<sim::SimChannel>>(fSimulationProducerLabel);
    consumes<art::Assns<simb::MCTruth, simb::MCParticle>>(fSimulationProducerLabel);
    consumes<std::vector<recob::Hit>>(fHitProducerLabel);
    consumes<std::vector<recob::Cluster>>(fClusterProducerLabel);
    consumes<art::Assns<recob::Cluster, recob::Hit>>(fHitProducerLabel);
    }
    // ------------------------------------------
    void EnergyOnWire::beginJob()
    {
      const double detectorLength = DetectorDiagonal(*fGeometry);
      art::ServiceHandle<art::TFileService> tfs;

      // Make Histograms
      fPDGCodeHist     = tfs->make<TH1D>("pdgcodes",";PDG Code;",                  5000, -2500, 2500);
      fMomentumHist    = tfs->make<TH1D>("mom",     ";particle Momentum (GeV);",    100, 0.,    10.);
      fTrackLengthHist = tfs->make<TH1D>("length",  ";particle track length (cm);", 200, 0, detectorLength);

      // Make Trees   
      fSimulationNtuple     = tfs->make<TTree>("EnergyOnWireTree",    "EnergyOnWireTree");
      fReconstructionNtuple = tfs->make<TTree>("EnergyOnWireReconstruction","EnergyOnWireReconstruction");

	  
      // Make Branches to complete n-tuples
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
      fSimulationNtuple->Branch("WireEnergy",    &fSimWireEnergy);
      fSimulationNtuple->Branch("WireNumber",  &fSimWireNumber);


    fReconstructionNtuple->Branch("Event",   &fEvent,          "Event/I");
    fReconstructionNtuple->Branch("SubRun",  &fSubRun,         "SubRun/I");
    fReconstructionNtuple->Branch("Run",     &fRun,            "Run/I");
    fReconstructionNtuple->Branch("TrackID", &fRecoTrackID,    "TrackID/I");
    fReconstructionNtuple->Branch("PDG",     &fRecoPDG,        "PDG/I");
    fReconstructionNtuple->Branch("NdEdx",   &fRecoNdEdxBins,  "NdEdx/I");
    fReconstructionNtuple->Branch("dEdx",    &fRecodEdxBins);
    }
    // ------------------------------------------
    void EnergyOnWire::beginRun(const art::Run& run)
    {
      art::ServiceHandle<sim::LArG4Parameters> larParameters;
      fElectronsToGeV = 1./larParameters->GeVToElectrons();
    }
    void EnergyOnWire::analyze(const art::Event& event)
    {
      fEvent = event.id().event();
      fRun = event.run();
      fSubRun = event.subRun();

      auto particleHandle = event.getValidHandle<std::vector<simb::MCParticle>>(fSimulationProducerLabel);
      auto simChannelHandle = event.getValidHandle<std::vector<sim::SimChannel>>(fSimulationProducerLabel);
      std::map<int,const simb::MCParticle*> particleMap;
      // Fill other information into output from the simulation
      for (auto const& particle : (*particleHandle))
	{
	  // Fill particleMap
	  fSimTrackID = particle.TrackId();
	  particleMap[fSimTrackID] = &particle;
	  // FIll PDG Histogram
	  fSimPDG = particle.PdgCode();
	  fPDGCodeHist->Fill(fSimPDG);
	  // Get kinetmatics of particle
	  const size_t numberTrajectoryPoints = particle.NumberTrajectoryPoints();
	  const int last = numberTrajectoryPoints - 1;
	  const TLorentzVector& positionStart = particle.Position(0);
	  const TLorentzVector& positionEnd   = particle.Position(last);
	  const TLorentzVector& momentumStart = particle.Momentum(0);
	  const TLorentzVector& momentumEnd   = particle.Momentum(last);
	  positionStart.GetXYZT(fStartXYZT);
	  positionEnd.GetXYZT  (fEndXYZT);
	  momentumStart.GetXYZT(fStartPE);
	  momentumEnd.GetXYZT  (fEndPE);
	  const double trackLength= (positionEnd-positionStart).Rho();
	  // Fill momentum and trackLength Histogram
	  fMomentumHist->Fill(momentumStart.P());
	  fTrackLengthHist->Fill(trackLength);

	} // End particle loop

      for (auto const& channel : (*simChannelHandle))
	{
	  float energyDeposit = 0;                 ///< Intialize dE value 
	  int wireNumber = (int)channel.Channel(); ///< Set wire number
	  fSimWireEnergy.clear();                    ///< Clear vectors at beginning of loop
	  fSimWireNumber.clear();

	  // Check if collection plane
	  auto const channelNumber = channel.Channel();
	  if (fGeometry->SignalType(channelNumber) != geo::kCollection)
	    continue;

	  // Begin looping over timeSlices, steps along track
	  auto const& timeSlices = channel.TDCIDEMap(); 
	  for (auto const& timeSlice : timeSlices)
	    {
	      auto const& chargeDeposits = timeSlice.second;     
	      for (auto const& chargeDeposit : chargeDeposits)
		{
		  energyDeposit += chargeDeposit.numElectrons * fElectronsToGeV; 
		  
		} // End chargeDeposit loop
	    } // End timeSlice loop
	  fSimWireEnergy.push_back(energyDeposit);
	  fSimWireNumber.push_back(wireNumber);
	  fSimulationNtuple->Fill();
	} // End channel loop

  // ###############   Reconstruction    Analysis Start from here #######################  
  /*  art::Handle< std::vector<recob::Hit> > hitHandle;
    if (!event.getByLabel(fHitProducerLabel, hitHandle)) return;

    std::map< int, std::vector<double> > dEdxMap;

    // For every Hit:
    for ( auto const& hit : (*hitHandle) )
      {
	// The channel associated with this hit.
	auto hitChannelNumber = hit.Channel();

	// We have a hit. For this example let's just focus on the
	// hits in the collection plane.
	if ( fGeometryService->SignalType( hitChannelNumber ) != geo::kCollection )
	  continue;

	MF_LOG_DEBUG("EnergyOnWire")
	  << "Hit in collection plane"
	  << std::endl;

	

	typedef sim::SimChannel::StoredTDC_t TDC_t;
	TDC_t start_tdc    = fTimeService->TPCTick2TDC( hit.StartTick() );
	TDC_t end_tdc      = fTimeService->TPCTick2TDC( hit.EndTick()   );
	TDC_t hitStart_tdc = fTimeService->TPCTick2TDC( hit.PeakTime() - 3.*hit.SigmaPeakTime() );
	TDC_t hitEnd_tdc   = fTimeService->TPCTick2TDC( hit.PeakTime() + 3.*hit.SigmaPeakTime() );

	start_tdc = std::max(start_tdc, hitStart_tdc);
	end_tdc   = std::min(end_tdc,   hitEnd_tdc  );

	// In the simulation section, we started with particles to find
	// channels with a matching track ID. Now we search in reverse:
	// search the SimChannels for matching channel number, then look
	// at the tracks inside the channel.

	for ( auto const& channel : (*simChannelHandle) )
	  {
	    auto simChannelNumber = channel.Channel();
	    if ( simChannelNumber != hitChannelNumber ) continue;

	    MF_LOG_DEBUG("EnergyOnWire")
	      << "SimChannel number = " << simChannelNumber
	      << std::endl;

	    // The time slices in this channel.
	    auto const& timeSlices = channel.TDCIDEMap();


	    // We have to create "dummy" time slices for the search.
	    sim::TDCIDE startTime;
	    sim::TDCIDE endTime;
	    startTime.first = start_tdc;
	    endTime.first   = end_tdc;

	    // Here are the fast searches:

	    // Find a pointer to the first channel with time >= start_tdc.
	    auto const startPointer
	      = std::lower_bound( timeSlices.begin(), timeSlices.end(), startTime, TDCIDETimeCompare);

	    // From that time slice, find the last channel with time < end_tdc.
	    auto const endPointer
	      = std::upper_bound( startPointer,       timeSlices.end(), endTime,   TDCIDETimeCompare);

	    // Did we find anything? If not, skip.
	    if ( startPointer == timeSlices.end() || startPointer == endPointer ) continue;
	    MF_LOG_DEBUG("EnergyOnWire")
	      << "Time slice start = " << (*startPointer).first
	      << std::endl;

	    // Loop over the channel times we found that match the hit
	    // times.
	    for ( auto slicePointer = startPointer; slicePointer != endPointer; slicePointer++)
	      {
		auto const timeSlice = *slicePointer;
		auto time = timeSlice.first;


		MF_LOG_DEBUG("EnergyOnWire")
		  << "Hit index = " << hit.LocalIndex()
		  << " channel number = " << hitChannelNumber
		  << " start TDC tick = " << hit.StartTick()
		  << " end TDC tick = " << hit.EndTick()
		  << " peak TDC tick = " << hit.PeakTime()
		  << " sigma peak time = " << hit.SigmaPeakTime()
		  << " adjusted start TDC tick = " << fTimeService->TPCTick2TDC(hit.StartTick())
		  << " adjusted end TDC tick = " << fTimeService->TPCTick2TDC(hit.EndTick())
		  << " adjusted peak TDC tick = " << fTimeService->TPCTick2TDC(hit.PeakTime())
		  << " adjusted start_tdc = " << start_tdc
		  << " adjusted end_tdc = " << end_tdc
		  << " time = " << time
		  << std::endl;

		// Loop over the energy deposits.
		auto const& energyDeposits = timeSlice.second;
		for ( auto const& energyDeposit : energyDeposits )
		  {
		   
		    auto search = particleMap.find( energyDeposit.trackID );

		    if ( search == particleMap.end() ) continue;

		    // "search" points to a pair in the map: <track ID, MCParticle*>
		    int trackID = (*search).first;
		    const simb::MCParticle& particle = *((*search).second);

		    // Is this a primary particle, with a PDG code that
		    // matches the user input?
		    if ( particle.Process() != "primary"
			 || particle.PdgCode() != fSelectedPDG )
		      continue;

		    // Determine the dE/dx of this particle.
		    const TLorentzVector& positionStart = particle.Position(0);
		    TVector3 location( energyDeposit.x,
				       energyDeposit.y,
				       energyDeposit.z );
		    double distance = ( location - positionStart.Vect() ).Mag();
		    unsigned int bin = int( distance / fBinSize );
		    double energy = energyDeposit.numElectrons * fElectronsToGeV;

		    auto& track_dEdX = dEdxMap[trackID];
		    if ( track_dEdX.size() < bin+1 )
		      {
			// Increase the vector size, padding it with
			// zeroes.
			track_dEdX.resize( bin+1, 0 );
		      }

		    // Add the energy to the dE/dx bins for this track.
		    track_dEdX[bin] += energy;

		  } // loop over energy deposits
	      } // loop over time slices
	  } // for each SimChannel
      } // for each Hit

    // We have a map of dE/dx vectors. Write each one of them to the
    // reconstruction n-tuple.
    for ( const auto& dEdxEntry : dEdxMap )
      {
	
	fRecoTrackID = dEdxEntry.first;

	
	fRecoPDG = particleMap[fRecoTrackID]->PdgCode();

	// Get the number of bins for this track.
	const std::vector<double>& dEdx = dEdxEntry.second;
	fRecoNdEdxBins = dEdx.size();

	fRecodEdxBins = dEdx;

	fReconstructionNtuple->Fill();
      }
*/



    } // End void analyze
    DEFINE_ART_MODULE(EnergyOnWire)
  } // End wireEnergy
} // End lar
namespace
{
  double DetectorDiagonal(geo::GeometryCore const& geom)
  {
    const double length = geom.DetLength();
    const double width = 2. * geom.DetHalfWidth();
    const double height = 2. * geom.DetHalfHeight();

    return std::sqrt(cet::sum_of_squares(length,width,height));
  }

 // bool TDCIDETimeCompare( const sim::TDCIDE& lhs, const sim::TDCIDE& rhs )
 // {
 //   return lhs.first < rhs.first;
  //}
}
