#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include <TVector2.h>

typedef reco::VertexCompositeCandidateCollection VCCandColl;

template<typename T>
class KJpsiProducer : public edm::EDFilter
{
public:
  KJpsiProducer(const edm::ParameterSet& pset);
  ~KJpsiProducer() {};

  bool filter(edm::Event& event, const edm::EventSetup& eventSetup);

private:
  bool isGoodTrack(const reco::TrackRef& track, const GlobalPoint& pvPoint) const;
  bool buildTransientTrack(edm::ESHandle<TransientTrackBuilder>& trackBuilder, const pat::Muon& muon, reco::TransientTrack& transTrack) const;
  bool buildTransientTrack(edm::ESHandle<TransientTrackBuilder>& trackBuilder, const pat::Electron& electron, reco::TransientTrack& transTrack) const;
  std::auto_ptr<edm::ValueMap<double> > getPtrValueMap( edm::OrphanHandle< VCCandColl > outHandle, const std::vector<double>& vmValue) const;

private:
  constexpr static double muonMass_ = 0.1056583715;
  constexpr static double electronMass_ = 0.0005;
  edm::InputTag leptonLabel_;
  edm::InputTag jetLabel_;
  edm::InputTag goodPrimaryVertexLabel_;

  unsigned int pdgId_, leptonId_;
  double leptonMass_;
  double rawMassMin_, rawMassMax_, massMin_, massMax_;

  double cut_minPt_, cut_maxEta_;
  double cut_trackChi2_, cut_trackSignif_, cut_DCA_;
  int cut_trackNHit_;
  double cut_vertexChi2_, cut_minLxy_, cut_maxLxy_, cut_vtxSignif_;

  unsigned int minNumber_, maxNumber_;

};

template<typename T>
KJpsiProducer<T>::KJpsiProducer(const edm::ParameterSet& pset)
{
  pdgId_ = 443;
  if ( typeid(T) == typeid(pat::Muon) )
  {
    leptonId_ = 13;
    leptonMass_ = muonMass_;
  }
  else if ( typeid(T) == typeid(pat::Electron) )
  {
    leptonId_ = 11;
    leptonMass_ = electronMass_;
  }

  leptonLabel_ = pset.getParameter<edm::InputTag>("src");
  jetLabel_ = pset.getParameter<edm::InputTag>("jetSrc");
  goodPrimaryVertexLabel_ = pset.getParameter<edm::InputTag>("primaryVertex");

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

  rawMassMin_ = pset.getParameter<double>("rawMassMin");
  rawMassMax_ = pset.getParameter<double>("rawMassMax");
  massMin_ = pset.getParameter<double>("massMin");
  massMax_ = pset.getParameter<double>("massMax");

  minNumber_ = pset.getParameter<unsigned int>("minNumber");
  maxNumber_ = pset.getParameter<unsigned int>("maxNumber");

  produces<reco::VertexCompositeCandidateCollection>();
  produces<edm::ValueMap<double> >("lxy");
  produces<edm::ValueMap<double> >("l3D");
  produces<edm::ValueMap<double> >("jetDR");
  produces<edm::ValueMap<double> >("vProb");

}

template<typename T>
bool KJpsiProducer<T>::filter(edm::Event& event, const edm::EventSetup& eventSetup)
{
  using namespace reco;
  using namespace edm;
  using namespace std;

  std::auto_ptr<VCCandColl> decayCands(new VCCandColl);
  std::vector<double> decayLengths;
  std::vector<double> decayLengths3D;
  std::vector<double> minJetDR;
  std::vector<double> vProb;

  edm::Handle< reco::VertexCollection >  goodPVHandle;
  event.getByLabel(goodPrimaryVertexLabel_ , goodPVHandle);
  if ( goodPVHandle->empty() ) {
    edm::OrphanHandle< VCCandColl > outHandle = event.put(decayCands); 
    event.put( getPtrValueMap( outHandle, decayLengths), "lxy");
    event.put( getPtrValueMap( outHandle, decayLengths3D), "l3D");
    event.put( getPtrValueMap( outHandle, vProb), "vProb");
    event.put( getPtrValueMap( outHandle, minJetDR), "jetDR");
    return false;
  }
  const reco::Vertex goodPV = goodPVHandle->at(0);

  const double pvx = goodPV.position().x();
  const double pvy = goodPV.position().y();
  const double pvz = goodPV.position().z();

  edm::ESHandle<TransientTrackBuilder> trackBuilder;
  eventSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",trackBuilder);

  edm::Handle<std::vector<T> > leptonHandle;
  event.getByLabel(leptonLabel_, leptonHandle);

  // Build transient tracks and separate by its charge
  std::vector<TransientTrack> transTracks1, transTracks2;
  for ( auto& lepton : *leptonHandle )
  {
    reco::TransientTrack transTrack;
    if ( !buildTransientTrack(trackBuilder, lepton, transTrack) ) continue;
    if ( !transTrack.impactPointTSCP().isValid() ) continue;

    if ( lepton.charge() < 0 ) transTracks1.push_back(transTrack);
    else if ( lepton.charge() > 0 ) transTracks2.push_back(transTrack);
  }

  // Make pairings
  for ( auto& transTrack1 : transTracks1 )
  {
    FreeTrajectoryState ipState1 = transTrack1.impactPointTSCP().theState();
    for ( auto& transTrack2 : transTracks2 )
    {
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
      SMatrixSym3D totalCov = goodPV.covariance() + vertex.covariance();
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
        for ( auto& refTrack : refittedTracks )
        {
          if ( refTrack.track().charge() < 0 ) refTrack1 = &refTrack;
          else refTrack2 = &refTrack;
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

      const double candE1 = hypot(mom1.mag(), leptonMass_);
      const double candE2 = hypot(mom2.mag(), leptonMass_);

      const math::XYZTLorentzVector candLVec(mom.x(), mom.y(), mom.z(), candE1+candE2);
      if ( massMin_ > candLVec.mass() or massMax_ < candLVec.mass() ) continue;

      // Match to muons
      VertexCompositeCandidate* cand = new VertexCompositeCandidate(0, candLVec, vtx, vtxCov, vtxChi2, vtxNdof);
      reco::LeafCandidate newLep1(+leptonId_, math::XYZTLorentzVector(mom1.x(), mom1.y(), mom1.z(), candE1));
      reco::LeafCandidate newLep2(-leptonId_, math::XYZTLorentzVector(mom2.x(), mom2.y(), mom2.z(), candE2));
      cand->addDaughter(newLep1);
      cand->addDaughter(newLep2);

      cand->setPdgId(pdgId_);
      AddFourMomenta addP4;
      addP4.set(*cand);

      edm::Handle<std::vector<pat::Jet> > jetHandle;
      event.getByLabel(jetLabel_, jetHandle);

      double minDR = 9999.;
      for ( auto& jet : *jetHandle )
      {
        const double dR = deltaR(*cand, jet);
        minDR = std::min(minDR, dR);
      }

      decayCands->push_back(*cand);
      LogDebug("KJpsiProducer")<<"minJetDR : "<<minDR<<"at Jet size : "<<jetHandle->size()<<std::endl;
      minJetDR.push_back(minDR);
      vProb.push_back( TMath::Prob( vtxChi2, (int) vtxNdof));
      decayLengths.push_back(rVtxMag);
      decayLengths3D.push_back(rVtxMag3D);
    }
  }

  const unsigned int nCands = decayCands->size();
  edm::OrphanHandle< VCCandColl > outHandle = event.put(decayCands);

  event.put( getPtrValueMap( outHandle, decayLengths), "lxy");
  event.put( getPtrValueMap( outHandle, decayLengths3D), "l3D");
  event.put( getPtrValueMap( outHandle, vProb), "vProb");
  event.put( getPtrValueMap( outHandle, minJetDR), "jetDR");

  return (nCands >= minNumber_ and nCands <= maxNumber_);
}

template<typename T>
std::auto_ptr<edm::ValueMap<double> > KJpsiProducer<T>::getPtrValueMap( edm::OrphanHandle< VCCandColl > outHandle, const std::vector<double>& vmValue) const
{
  std::auto_ptr<edm::ValueMap<double> > temp_valuemap( new edm::ValueMap<double> );
  edm::ValueMap<double>::Filler filler(*temp_valuemap);
  filler.insert(outHandle, vmValue.begin(), vmValue.end());
  filler.fill();
  return temp_valuemap;
}

template<typename T>
bool KJpsiProducer<T>::buildTransientTrack(edm::ESHandle<TransientTrackBuilder>& trackBuilder,
                                        const pat::Muon& muon, reco::TransientTrack& transTrack) const
{
  reco::TrackRef trackRef;
  if ( muon.isGlobalMuon() ) trackRef = muon.globalTrack();
  else if ( muon.isTrackerMuon() ) trackRef = muon.innerTrack();
  else return false;

  transTrack = trackBuilder->build(trackRef);
  return true;
}

template<typename T>
bool KJpsiProducer<T>::buildTransientTrack(edm::ESHandle<TransientTrackBuilder>& trackBuilder,
                                        const pat::Electron& electron, reco::TransientTrack& transTrack) const
{
  if ( electron.gsfTrack().isNull() ) return false;

  transTrack = trackBuilder->build(electron.gsfTrack());
  return true;
}

#include "FWCore/Framework/interface/MakerMacros.h"
typedef KJpsiProducer<pat::Muon> KJpsiMuMuProducer;
typedef KJpsiProducer<pat::Electron> KJpsiElElProducer;
DEFINE_FWK_MODULE(KJpsiMuMuProducer);
DEFINE_FWK_MODULE(KJpsiElElProducer);

