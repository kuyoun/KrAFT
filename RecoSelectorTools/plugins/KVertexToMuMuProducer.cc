#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
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

class KVertexToMuMuProducer : public edm::EDFilter
{
public:
  KVertexToMuMuProducer(const edm::ParameterSet& pset);
  ~KVertexToMuMuProducer() {};

  bool filter(edm::Event& event, const edm::EventSetup& eventSetup);

private:
  bool isGoodTrack(const reco::TrackRef& track, const reco::BeamSpot* beamSpot) const;
  const reco::Muon* matchMuon(const reco::TrackRef trackRef,
                              reco::MuonCollection::const_iterator muonsBegin,
                              reco::MuonCollection::const_iterator muonsEnd);

private:
  constexpr static double muonMass_ = 0.1056583715; 
  edm::InputTag muonLabel_;
  edm::InputTag trackLabel_;

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
};

KVertexToMuMuProducer::KVertexToMuMuProducer(const edm::ParameterSet& pset)
{
  edm::ParameterSet trackPSet = pset.getParameter<edm::ParameterSet>("track");
  trackLabel_ = trackPSet.getParameter<edm::InputTag>("src");
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

  edm::ParameterSet muonPSet = pset.getParameter<edm::ParameterSet>("muon");
  muonLabel_ = muonPSet.getParameter<edm::InputTag>("src");
  muonDPt_ = muonPSet.getParameter<double>("dPtRel");
  muonDR_  = muonPSet.getParameter<double>("dR");

  pdgId_ = pset.getParameter<unsigned int>("pdgId");
  rawMassMin_ = pset.getParameter<double>("rawMassMin");
  rawMassMax_ = pset.getParameter<double>("rawMassMax");
  massMin_ = pset.getParameter<double>("massMin");
  massMax_ = pset.getParameter<double>("massMax");

  minNumber_ = pset.getParameter<unsigned int>("minNumber");
  maxNumber_ = pset.getParameter<unsigned int>("maxNumber");

  produces<reco::VertexCompositeCandidateCollection>();
  produces<std::vector<double> >("lxy");
}

bool KVertexToMuMuProducer::filter(edm::Event& event, const edm::EventSetup& eventSetup)
{
  using namespace reco;
  using namespace edm;
  using namespace std;

  typedef reco::VertexCompositeCandidateCollection VCCandColl;

  std::auto_ptr<VCCandColl> decayCands(new VCCandColl);
  std::auto_ptr<std::vector<double> > decayLengths(new std::vector<double>);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  event.getByLabel("offlineBeamSpot", beamSpotHandle);

  edm::Handle<reco::TrackCollection> trackHandle;
  event.getByLabel(trackLabel_, trackHandle);

  edm::ESHandle<MagneticField> bFieldHandle;
  eventSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
  bField_ = bFieldHandle.product();

  //edm::ESHandle<TrackerGeometry> trackerGeomHandle;
  //eventSetup.get<TrackerDigiGeometryRecord>().get(trackerGeomHandle);
  //trackerGeom_ = trackerGeomHandle.product();

  edm::ESHandle<GlobalTrackingGeometry> glbTkGeomHandle;
  eventSetup.get<GlobalTrackingGeometryRecord>().get(glbTkGeomHandle);
  //glbTkGeom_ = glbTkGeomHandle.product();

  edm::Handle<reco::MuonCollection> muonHandle;
  event.getByLabel(muonLabel_, muonHandle);

  for ( int i=0, n=trackHandle->size(); i<n; ++i )
  {
    TrackRef trackRef1(trackHandle, i);
    // Positive particle in the 1st index (pi+, proton, K+...)
    if ( trackRef1->charge() < 0 ) continue;
    if ( !isGoodTrack(trackRef1, beamSpotHandle.product()) ) continue;

    TransientTrack transTrack1(*trackRef1, bField_, glbTkGeomHandle);
    if ( !transTrack1.impactPointTSCP().isValid() ) continue;
    FreeTrajectoryState ipState1 = transTrack1.impactPointTSCP().theState();

    for ( int j=0; j<n; ++j )
    {
      TrackRef trackRef2(trackHandle, j);
      // Negative particle in the 2nd index (pi-, anti-proton, K-...)
      if ( trackRef2->charge() > 0 ) continue;
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
      GlobalPoint beamSpotPos(beamSpotHandle->position().x(),
                              beamSpotHandle->position().y(),
                              beamSpotHandle->position().z());

      SMatrixSym3D totalCov = beamSpotHandle->rotatedCovariance3D() + vertex.covariance();
      SVector3 distanceVector(vertex.x() - beamSpotPos.x(), vertex.y() - beamSpotPos.y(), 0.);

      double rVtxMag = ROOT::Math::Mag(distanceVector);
      double sigmaRvtxMag = sqrt(ROOT::Math::Similarity(totalCov, distanceVector)) / rVtxMag;
      if( rVtxMag < cut_minLxy_ or rVtxMag > cut_maxLxy_ or rVtxMag / sigmaRvtxMag < cut_vtxSignif_ ) continue;

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
          if ( refTrack->track().charge() > 0 ) refTrack1 = &*refTrack;
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
      if ( massMin_ > candLVec.mass() or massMax_ < candLVec.mass() ) continue;

      // Match to muons
      const reco::Muon* muon1 = matchMuon(trackRef1, muonHandle->begin(), muonHandle->end());
      const reco::Muon* muon2 = matchMuon(trackRef2, muonHandle->begin(), muonHandle->end());
      if ( !muon1 and !muon2 ) continue;

      RecoChargedCandidate cand1(trackRef1->charge(), Particle::LorentzVector(mom1.x(), mom1.y(), mom1.z(), candE1), vtx);
      RecoChargedCandidate cand2(trackRef2->charge(), Particle::LorentzVector(mom2.x(), mom2.y(), mom2.z(), candE2), vtx);
      cand1.setTrack(trackRef1);
      cand2.setTrack(trackRef2);
      const int pdgId1 = -13*trackRef1->charge();
      const int pdgId2 = -13*trackRef2->charge();
      cand1.setPdgId(pdgId1);
      cand2.setPdgId(pdgId2);
      VertexCompositeCandidate* cand = new VertexCompositeCandidate(0, candLVec, vtx, vtxCov, vtxChi2, vtxNdof);
      cand->addDaughter(cand1);
      cand->addDaughter(cand2);

      cand->setPdgId(pdgId_);
      AddFourMomenta addP4;
      addP4.set(*cand);

      decayCands->push_back(*cand);
      decayLengths->push_back(rVtxMag);
      
    }
  }

  const unsigned int nCands = decayCands->size();
  event.put(decayCands);
  event.put(decayLengths, "lxy");

  return (nCands >= minNumber_ and nCands <= maxNumber_);
}

bool KVertexToMuMuProducer::isGoodTrack(const reco::TrackRef& track, const reco::BeamSpot* beamSpot) const
{
  const static reco::TrackBase::TrackQuality trackQual = reco::TrackBase::qualityByName("loose");
  if ( !track->quality(trackQual) ) return false;
  if ( track->normalizedChi2() >= cut_trackChi2_ ) return false;
  if ( track->numberOfValidHits() < cut_trackNHit_ ) return false;
  if ( track->pt() < cut_minPt_ or abs(track->eta()) > cut_maxEta_ ) return false;

  FreeTrajectoryState initialFTS = trajectoryStateTransform::initialFreeState(*track, bField_);
  TSCBLBuilderNoMaterial blsBuilder;
  TrajectoryStateClosestToBeamLine tscb( blsBuilder(initialFTS, *beamSpot) );
  if ( !tscb.isValid() ) return false;
  if ( tscb.transverseImpactParameter().significance() <= cut_trackSignif_ ) return false;

  return true;
}

const reco::Muon* KVertexToMuMuProducer::matchMuon(const reco::TrackRef trackRef, 
                                                   reco::MuonCollection::const_iterator muonsBegin,
                                                   reco::MuonCollection::const_iterator muonsEnd)
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
