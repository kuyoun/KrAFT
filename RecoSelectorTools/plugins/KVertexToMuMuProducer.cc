#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
//#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
//#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/VolumeBasedEngine/interface/VolumeBasedMagneticField.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/GluedGeomDet.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

#include <string>
#include <fstream>

//#define DEBUGPLOT

#ifdef DEBUGPLOT
#include "TH1F.h"
#include "TH2F.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#endif

class KVertexToMuMuProducer : public edm::EDFilter
{
public:
  KVertexToMuMuProducer(const edm::ParameterSet& pset);
  ~KVertexToMuMuProducer() {};

  bool filter(edm::Event& event, const edm::EventSetup& eventSetup);

private:
  bool isGoodTrack(const reco::TrackRef& track, const reco::BeamSpot* beamSpot) const;
  const pat::Muon* matchMuon(const reco::TrackRef trackRef,
                             pat::MuonCollection::const_iterator muonsBegin,
                             pat::MuonCollection::const_iterator muonsEnd);

private:
  constexpr static double muonMass_ = 0.1056583715;
  edm::InputTag muonLabel_;

  unsigned int pdgId_;
  double rawMassMin_, rawMassMax_, massMin_, massMax_;

  double cut_minPt_, cut_maxEta_;
  double cut_trackChi2_, cut_trackSignif_, cut_DCA_;
  int cut_trackNHit_;
  double cut_vertexChi2_, cut_minLxy_, cut_maxLxy_, cut_vtxSignif_;

  double muonDPt_, muonDR_;

  unsigned int minNumber_, maxNumber_;

  //const TrackerGeometry* trackerGeom_;
  const MagneticField* bField_;
#ifdef DEBUGPLOT
  TH1F* hRawMass_, * hFitMass_;
  TH2F* hRawMassVsFitMass_;
#endif
};

KVertexToMuMuProducer::KVertexToMuMuProducer(const edm::ParameterSet& pset)
{
  muonLabel_ = pset.getParameter<edm::InputTag>("src");

  edm::ParameterSet trackPSet = pset.getParameter<edm::ParameterSet>("track");
  cut_minPt_ = trackPSet.getParameter<double>("minPt");
  cut_maxEta_ = trackPSet.getParameter<double>("maxEta");
  cut_trackChi2_ = trackPSet.getParameter<double>("chi2");
  cut_trackNHit_  = trackPSet.getParameter<int>("nHit");
  cut_trackSignif_ = trackPSet.getParameter<double>("signif");
  cut_DCA_ = trackPSet.getParameter<double>("DCA");

  edm::ParameterSet vertexPSet = pset.getParameter<edm::ParameterSet>("vertex");
  cut_vertexChi2_ = vertexPSet.getParameter<double>("chi2");
  cut_minLxy_ = vertexPSet.getParameter<double>("minLxy");
  cut_maxLxy_ = vertexPSet.getParameter<double>("maxLxy");
  cut_vtxSignif_ = vertexPSet.getParameter<double>("signif");

  pdgId_ = pset.getParameter<unsigned int>("pdgId");
  rawMassMin_ = pset.getParameter<double>("rawMassMin");
  rawMassMax_ = pset.getParameter<double>("rawMassMax");
  massMin_ = pset.getParameter<double>("massMin");
  massMax_ = pset.getParameter<double>("massMax");

  minNumber_ = pset.getParameter<unsigned int>("minNumber");
  maxNumber_ = pset.getParameter<unsigned int>("maxNumber");

  produces<reco::VertexCompositeCandidateCollection>();
  produces<std::vector<double> >("lxy");
  produces<std::vector<double> >("l3D");

#ifdef DEBUGPLOT
  edm::Service<TFileService> fs;
  hRawMass_ = fs->make<TH1F>("hRawMass", "raw mass;Raw mass (GeV/c^{2});Entries", 100, 0, 5);
  hFitMass_ = fs->make<TH1F>("hFitMass", "fit mass;Fit mass (GeV/c^{2});Entries", 100, 0, 5);
  hRawMassVsFitMass_ = fs->make<TH2F>("hRawMassVsFitMass", "raw vs fit;Raw mass (GeV/c^{2};Fit mass (GeV/c^{2}", 100, 0, 5, 100, 0, 5);
#endif
}

bool KVertexToMuMuProducer::filter(edm::Event& event, const edm::EventSetup& eventSetup)
{
  using namespace reco;
  using namespace edm;
  using namespace std;

  typedef reco::VertexCompositeCandidateCollection VCCandColl;

  std::auto_ptr<VCCandColl> decayCands(new VCCandColl);
  std::auto_ptr<std::vector<double> > decayLengths(new std::vector<double>);
  std::auto_ptr<std::vector<double> > decayLengths3D(new std::vector<double>);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  event.getByLabel("offlineBeamSpot", beamSpotHandle);
  const double pvx = beamSpotHandle->position().x();
  const double pvy = beamSpotHandle->position().y();
  const double pvz = beamSpotHandle->position().z();

  edm::ESHandle<MagneticField> bFieldHandle;
  eventSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
  bField_ = bFieldHandle.product();

  //edm::ESHandle<TrackerGeometry> trackerGeomHandle;
  //eventSetup.get<TrackerDigiGeometryRecord>().get(trackerGeomHandle);
  //trackerGeom_ = trackerGeomHandle.product();

  edm::ESHandle<GlobalTrackingGeometry> glbTkGeomHandle;
  eventSetup.get<GlobalTrackingGeometryRecord>().get(glbTkGeomHandle);
  //glbTkGeom_ = glbTkGeomHandle.product();

  edm::Handle<std::vector<pat::Muon> > muonHandle;
  event.getByLabel(muonLabel_, muonHandle);
	


  for ( int i=0, n= muonHandle->size() ; i<n; ++i )
  {
    const pat::Muon& muon1 = muonHandle->at(i);
    if ( muon1.charge() >= 0 ) continue;
    TrackRef trackRef1;
    if ( muon1.isGlobalMuon() ) trackRef1 = muon1.globalTrack();
    else if ( muon1.isTrackerMuon() ) trackRef1 = muon1.innerTrack();
    else continue;
    if ( !isGoodTrack(trackRef1, beamSpotHandle.product()) ) continue;

    TransientTrack transTrack1(*trackRef1, bField_, glbTkGeomHandle);
    if ( !transTrack1.impactPointTSCP().isValid() ) continue;
    FreeTrajectoryState ipState1 = transTrack1.impactPointTSCP().theState();

    for ( int j=0; j<n; ++j )
    {
      const pat::Muon& muon2 = muonHandle->at(j);
      if ( muon2.charge() <= 0 ) continue;
      TrackRef trackRef2;
      if ( muon2.isGlobalMuon() ) trackRef2 = muon2.globalTrack();
      else if ( muon2.isTrackerMuon() ) trackRef2 = muon2.innerTrack();
      else continue;
      if ( !isGoodTrack(trackRef2, beamSpotHandle.product()) ) continue;

      TransientTrack transTrack2(*trackRef2, bField_, glbTkGeomHandle);
      if ( !transTrack2.impactPointTSCP().isValid() ) continue;
      FreeTrajectoryState ipState2 = transTrack2.impactPointTSCP().theState();

      // Measure distance between tracks at their closest approach
      ClosestApproachInRPhi cApp;
      cApp.calculate(ipState1, ipState2);
      if ( !cApp.status() ) continue;
      const float dca = fabs(cApp.distance());
      if ( dca < 0. || dca > cut_DCA_ ) continue;
      GlobalPoint cxPt = cApp.crossingPoint();
      if (std::hypot(cxPt.x(), cxPt.y()) > 120. || std::abs(cxPt.z()) > 300.) continue;

      TrajectoryStateClosestToPoint caState1 = transTrack1.trajectoryStateClosestToPoint(cxPt);
      TrajectoryStateClosestToPoint caState2 = transTrack2.trajectoryStateClosestToPoint(cxPt);
      if ( !caState1.isValid() or !caState2.isValid() ) continue;

      const double rawEnergy = std::hypot(caState1.momentum().mag(), muonMass_)
                             + std::hypot(caState2.momentum().mag(), muonMass_);
      const double rawMass = sqrt(rawEnergy*rawEnergy - (caState1.momentum()+caState2.momentum()).mag2());

#ifdef DEBUGPLOT
      hRawMass_->Fill(rawMass);
#endif
      if ( rawMassMin_ > rawMass or rawMassMax_ < rawMass ) continue;

      // Build Vertex
      std::vector<TransientTrack> transTracks;
      transTracks.push_back(transTrack1);
      transTracks.push_back(transTrack2);
      KalmanVertexFitter fitter(true);
      TransientVertex transVertex = fitter.vertex(transTracks);

      if ( !transVertex.isValid() or transVertex.totalChiSquared() < 0. ) continue;

      const reco::Vertex vertex = transVertex;
      if ( vertex.normalizedChi2() > cut_vertexChi2_ ) continue;

      std::vector<TransientTrack> refittedTracks;
      if ( transVertex.hasRefittedTracks() ) refittedTracks = transVertex.refittedTracks();

      typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
      typedef ROOT::Math::SVector<double, 3> SVector3;

      GlobalPoint vtxPos(vertex.x(), vertex.y(), vertex.z());

      SMatrixSym3D totalCov = beamSpotHandle->rotatedCovariance3D() + vertex.covariance();
      SVector3 distanceVectorXY(vertex.x() - pvx, vertex.y() - pvy, 0.);

      double rVtxMag = ROOT::Math::Mag(distanceVectorXY);
      double sigmaRvtxMag = sqrt(ROOT::Math::Similarity(totalCov, distanceVectorXY)) / rVtxMag;
      if( rVtxMag < cut_minLxy_ or rVtxMag > cut_maxLxy_ or rVtxMag / sigmaRvtxMag < cut_vtxSignif_ ) continue;

      SVector3 distanceVector3D(vertex.x() - pvx, vertex.y() - pvy, vertex.z() - pvz);
      const double rVtxMag3D = ROOT::Math::Mag(distanceVector3D);

      // Cuts finished, now we create the candidates and push them back into the collections.

      std::auto_ptr<TrajectoryStateClosestToPoint> traj1;
      std::auto_ptr<TrajectoryStateClosestToPoint> traj2;

      if ( refittedTracks.empty() )
      {
        traj1.reset(new TrajectoryStateClosestToPoint(transTrack1.trajectoryStateClosestToPoint(vtxPos)));
        traj2.reset(new TrajectoryStateClosestToPoint(transTrack2.trajectoryStateClosestToPoint(vtxPos)));
      }
      else
      {
        TransientTrack* refTrack1 = 0, * refTrack2 = 0;
        for ( std::vector<TransientTrack>::iterator refTrack = refittedTracks.begin();
              refTrack != refittedTracks.end(); ++refTrack )
        {
          if ( refTrack->track().charge() < 0 ) refTrack1 = &*refTrack;
          else refTrack2 = &*refTrack;
        }
        if ( refTrack1 == 0 or refTrack2 == 0 ) continue;
        traj1.reset(new TrajectoryStateClosestToPoint(refTrack1->trajectoryStateClosestToPoint(vtxPos)));
        traj2.reset(new TrajectoryStateClosestToPoint(refTrack2->trajectoryStateClosestToPoint(vtxPos)));
      }
      if( !traj1->isValid() or !traj2->isValid() ) continue;

      GlobalVector mom1(traj1->momentum());
      GlobalVector mom2(traj2->momentum());
      GlobalVector mom(mom1+mom2);

      //cleanup stuff we don't need anymore
      traj1.reset();
      traj2.reset();

      Particle::Point vtx(vertex.x(), vertex.y(), vertex.z());
      const Vertex::CovarianceMatrix vtxCov(vertex.covariance());
      double vtxChi2(vertex.chi2());
      double vtxNdof(vertex.ndof());

      const double candE1 = hypot(mom1.mag(), muonMass_);
      const double candE2 = hypot(mom2.mag(), muonMass_);

      Particle::LorentzVector candLVec(mom.x(), mom.y(), mom.z(), candE1+candE2);
#ifdef DEBUGPLOT
      hFitMass_->Fill(candLVec.mass());
      hRawMassVsFitMass_->Fill(rawMass, candLVec.mass());
#endif
      if ( massMin_ > candLVec.mass() or massMax_ < candLVec.mass() ) continue;
//	std::cout<<"all step passed!"<<std::endl;
      // Match to muons
      pat::Muon newMuon1(muon1);
      pat::Muon newMuon2(muon2);
      newMuon1.setP4(reco::Candidate::LorentzVector(mom1.x(), mom1.y(), mom1.z(), candE1));
      newMuon2.setP4(reco::Candidate::LorentzVector(mom2.x(), mom2.y(), mom2.z(), candE2));
      VertexCompositeCandidate* cand = new VertexCompositeCandidate(0, candLVec, vtx, vtxCov, vtxChi2, vtxNdof);
      cand->addDaughter(newMuon1);
      cand->addDaughter(newMuon2);

      cand->setPdgId(pdgId_);
      AddFourMomenta addP4;
      addP4.set(*cand);

      decayCands->push_back(*cand);
      decayLengths->push_back(rVtxMag);
      decayLengths3D->push_back(rVtxMag3D);

    }
  }

  const unsigned int nCands = decayCands->size();
	  //std::cout<<"Event : "<<event.id()<<"   total passed :  "<<nCands<<std::endl;
  event.put(decayCands);
  event.put(decayLengths, "lxy");
  event.put(decayLengths3D, "l3D");

  return (nCands >= minNumber_ and nCands <= maxNumber_);
}



bool KVertexToMuMuProducer::isGoodTrack(const reco::TrackRef& track, const reco::BeamSpot* beamSpot) const
{
  //Turn off quality cuts - we discovered muon tracks does not pass this selection
  //const static reco::TrackBase::TrackQuality trackQual = reco::TrackBase::qualityByName("loose");
  //if ( !track->quality(trackQual) ) return false;
  //if ( track->normalizedChi2() >= cut_trackChi2_ ) return false;
  //if ( track->numberOfValidHits() < cut_trackNHit_ ) return false;
  if ( track->pt() < cut_minPt_ or abs(track->eta()) > cut_maxEta_ ) return false;

  FreeTrajectoryState initialFTS = trajectoryStateTransform::initialFreeState(*track, bField_);
  TSCBLBuilderNoMaterial blsBuilder;
  TrajectoryStateClosestToBeamLine tscb( blsBuilder(initialFTS, *beamSpot) );
  if ( !tscb.isValid() ) return false;
  //if ( tscb.transverseImpactParameter().significance() <= cut_trackSignif_ ) return false;

  return true;
}

const pat::Muon* KVertexToMuMuProducer::matchMuon(const reco::TrackRef trackRef,
                                                  pat::MuonCollection::const_iterator muonsBegin,
                                                  pat::MuonCollection::const_iterator muonsEnd)
{
  const double trackPt = trackRef->pt();
  const double trackEta = trackRef->eta();
  const double trackPhi = trackRef->phi();
  for ( auto mu = muonsBegin; mu != muonsEnd; ++mu )
  {
    if ( !mu->isTrackerMuon() and !mu->isPFMuon() ) continue;
    const reco::TrackRef muonTrackRef = mu->innerTrack();
    if ( trackRef == muonTrackRef ) return &(*mu);

    const double muonPt = mu->pt();
    const double muonEta = mu->eta();
    const double muonPhi = mu->phi();
    if ( std::abs(muonPt-trackPt)/trackPt < muonDPt_ and
         reco::deltaR(trackEta, muonEta, trackPhi, muonPhi) > muonDR_ ) return &(*mu);
  }
  return 0;
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(KVertexToMuMuProducer);
