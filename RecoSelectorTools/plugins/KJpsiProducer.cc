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
//#define DEBUGPLOT

#ifdef DEBUGPLOT
#include "TH1F.h"
#include "TH2F.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#endif


typedef reco::VertexCompositeCandidateCollection VCCandColl;

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
  std::auto_ptr<edm::ValueMap<double> > getPtrValueMap( edm::OrphanHandle< VCCandColl > outHandle, std::vector<double> vmValue); 

private:
  constexpr static double muonMass_ = 0.1056583715;
  constexpr static double electronMass_ = 0.0005;
  edm::EDGetTokenT<std::vector<pat::Muon> > muonToken_;
  edm::EDGetTokenT<std::vector<pat::Electron> > electronToken_;
  edm::EDGetTokenT<std::vector<pat::Jet> > jetToken_;
  edm::EDGetTokenT<reco::VertexCollection> goodPrimaryVertexToken_;

  unsigned int pdgId_;
  double rawMassMin_, rawMassMax_, massMin_, massMax_;

  double cut_minPt_, cut_maxEta_;
  double cut_trackChi2_, cut_trackSignif_, cut_DCA_;
  int cut_trackNHit_;
  double cut_vertexChi2_, cut_minLxy_, cut_maxLxy_, cut_vtxSignif_;

  unsigned int minNumber_, maxNumber_;

};
KJpsiProducer::KJpsiProducer(const edm::ParameterSet& pset)
{
  muonToken_ = consumes<std::vector<pat::Muon> >(pset.getParameter<edm::InputTag>("muonSrc"));
  electronToken_ = consumes<std::vector<pat::Electron> >(pset.getParameter<edm::InputTag>("electronSrc"));
  jetToken_ = consumes<std::vector<pat::Jet> >(pset.getParameter<edm::InputTag>("jetSrc"));
  goodPrimaryVertexToken_ = consumes<reco::VertexCollection>(pset.getParameter<edm::InputTag>("primaryVertex"));

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
  produces<edm::ValueMap<double> >("lxy");
  produces<edm::ValueMap<double> >("l3D");
  produces<edm::ValueMap<double> >("jetDR");
  produces<edm::ValueMap<double> >("vProb");

}

bool KJpsiProducer::filter(edm::Event& event, const edm::EventSetup& eventSetup)
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
  event.getByToken(goodPrimaryVertexToken_ , goodPVHandle);
  if ( goodPVHandle->size() ==0 ) {
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

  edm::Handle<std::vector<pat::Muon> > muonHandle;
  event.getByToken(muonToken_, muonHandle);

  edm::Handle<std::vector<pat::Electron> > electronHandle;
  event.getByToken(electronToken_, electronHandle);

  vector<TransientTrack> transTrack1;
  vector<TransientTrack> transTrack2;
  std::vector<reco::LeafCandidate> lep1;
  std::vector<reco::LeafCandidate> lep2;
  vector<int> lep_id;
  for ( int i=0, n=muonHandle->size(); i<n; ++i )
  {
    reco::TransientTrack t1, t2;
    const pat::Muon& muon1 = muonHandle->at(i);
    if ( muon1.charge() >= 0 ) continue;
    if ( !buildTransientTrack(trackBuilder, muon1, t1) ) continue;
    if ( !t1.impactPointTSCP().isValid() ) continue;

    for ( int j=0; j<n; ++j )
    {
      const pat::Muon& muon2 = muonHandle->at(j);
      if ( muon2.charge() <= 0 ) continue;
      if ( !buildTransientTrack(trackBuilder, muon2, t2) ) continue;
      if ( !t2.impactPointTSCP().isValid() ) continue;

      transTrack1.push_back(t1);
      transTrack2.push_back(t2);
      lep1.push_back( reco::LeafCandidate(muon1) );
      lep2.push_back( reco::LeafCandidate(muon2) );
      lep_id.push_back(13);
    }
  }
  for( unsigned int i=0 ; i< transTrack1.size();++i) {
    FreeTrajectoryState ipState1 = transTrack1[i].impactPointTSCP().theState();
    FreeTrajectoryState ipState2 = transTrack2[i].impactPointTSCP().theState();

    // Measure distance between tracks at their closest approach
    ClosestApproachInRPhi cApp;
    cApp.calculate(ipState1, ipState2);
    if ( !cApp.status() ) continue;
    const float dca = fabs(cApp.distance());
    if ( dca < 0. || dca > cut_DCA_ ) continue;
    GlobalPoint cxPt = cApp.crossingPoint();
    if (std::hypot(cxPt.x(), cxPt.y()) > 120. || std::abs(cxPt.z()) > 300.) continue;
    TrajectoryStateClosestToPoint caState1 = transTrack1[i].trajectoryStateClosestToPoint(cxPt);
    TrajectoryStateClosestToPoint caState2 = transTrack2[i].trajectoryStateClosestToPoint(cxPt);
    if ( !caState1.isValid() or !caState2.isValid() ) continue;

    // Build Vertex
    std::vector<TransientTrack> transTracks;
    transTracks.push_back(transTrack1[i]);
    transTracks.push_back(transTrack2[i]);
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

    //SMatrixSym3D totalCov = beamSpotHandle->rotatedCovariance3D() + vertex.covariance();
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
      traj1.reset(new TrajectoryStateClosestToPoint(transTrack1[i].trajectoryStateClosestToPoint(vtxPos)));
      traj2.reset(new TrajectoryStateClosestToPoint(transTrack2[i].trajectoryStateClosestToPoint(vtxPos)));
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

    double lepMass = 0.0;
    if ( lep_id[i] == 13) lepMass = muonMass_; 
    else if ( lep_id[i] ==1 ) lepMass = electronMass_; 
    const double candE1 = hypot(mom1.mag(), lepMass);
    const double candE2 = hypot(mom2.mag(), lepMass);

    Particle::LorentzVector candLVec(mom.x(), mom.y(), mom.z(), candE1+candE2);
    if ( massMin_ > candLVec.mass() or massMax_ < candLVec.mass() ) continue;

    // Match to muons
    reco::LeafCandidate newLep1( lep1[i] );
    reco::LeafCandidate newLep2( lep2[i] );
    newLep1.setP4(reco::Candidate::LorentzVector(mom1.x(), mom1.y(), mom1.z(), candE1));
    newLep2.setP4(reco::Candidate::LorentzVector(mom2.x(), mom2.y(), mom2.z(), candE2));
    VertexCompositeCandidate* cand = new VertexCompositeCandidate(0, candLVec, vtx, vtxCov, vtxChi2, vtxNdof);
    cand->addDaughter(newLep1);
    cand->addDaughter(newLep2);

    cand->setPdgId(pdgId_);
    AddFourMomenta addP4;
    addP4.set(*cand);

    edm::Handle<std::vector<pat::Jet> > jetHandle;
    event.getByToken(jetToken_, jetHandle);
     
    double minDR = 9999.;
    for ( int i=0, n=jetHandle->size(); i<n; ++i ) {
      Double_t deta = cand->eta() - jetHandle->at(i).eta();
      Double_t dphi = TVector2::Phi_mpi_pi( cand->phi()- jetHandle->at(i).phi());
      Double_t dr = TMath::Sqrt( deta*deta+dphi*dphi );
      if ( minDR > dr ) minDR = dr ;
    }
    decayCands->push_back(*cand);
    LogDebug("KJpsiProducer")<<"minJetDR : "<<minDR<<"at Jet size : "<<jetHandle->size()<<std::endl;
    minJetDR.push_back(minDR);
    vProb.push_back( TMath::Prob( vtxChi2, (int) vtxNdof)); 
    decayLengths.push_back(rVtxMag);
    decayLengths3D.push_back(rVtxMag3D);

  }

  const unsigned int nCands = decayCands->size();
  edm::OrphanHandle< VCCandColl > outHandle = event.put(decayCands); 

  event.put( getPtrValueMap( outHandle, decayLengths), "lxy");
  event.put( getPtrValueMap( outHandle, decayLengths3D), "l3D");
  event.put( getPtrValueMap( outHandle, vProb), "vProb");
  event.put( getPtrValueMap( outHandle, minJetDR), "jetDR");

  return (nCands >= minNumber_ and nCands <= maxNumber_);
}


std::auto_ptr<edm::ValueMap<double> > KJpsiProducer::getPtrValueMap( edm::OrphanHandle< VCCandColl > outHandle, std::vector<double> vmValue) 
{
  std::auto_ptr<edm::ValueMap<double> > temp_valuemap( new edm::ValueMap<double> );
  edm::ValueMap<double>::Filler filler(*temp_valuemap);
  filler.insert(outHandle, vmValue.begin(), vmValue.end() );
  filler.fill();
  return temp_valuemap;
}

bool KJpsiProducer::buildTransientTrack(edm::ESHandle<TransientTrackBuilder>& trackBuilder,
                                              const pat::Muon& muon, reco::TransientTrack& transTrack) const
{
  reco::TrackRef trackRef;
  if ( muon.isGlobalMuon() ) trackRef = muon.globalTrack();
  else if ( muon.isTrackerMuon() ) trackRef = muon.innerTrack();
  else return false;

  transTrack = trackBuilder->build(trackRef);
  return true;
}
bool KJpsiProducer::buildTransientTrack(edm::ESHandle<TransientTrackBuilder>& trackBuilder,
                                              const pat::Electron& electron, reco::TransientTrack& transTrack) const
{
  if ( electron.gsfTrack().isNull() ) return false;

  transTrack = trackBuilder->build(electron.gsfTrack());
  return true;
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(KJpsiProducer);

