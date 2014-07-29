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
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include <TVector2.h>
#include <TMath.h>
//#define DEBUGPLOT

#ifdef DEBUGPLOT
#include "TH1F.h"
#include "TH2F.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#endif


class KJpsiToElElProducer : public edm::EDFilter
{
public:
  KJpsiToElElProducer(const edm::ParameterSet& pset);
  ~KJpsiToElElProducer() {};

  bool filter(edm::Event& event, const edm::EventSetup& eventSetup);

private:
  bool isGoodTrack(const reco::TrackRef& track, const GlobalPoint& pvPoint) const;
  bool buildTransientTrack(edm::ESHandle<TransientTrackBuilder>& trackBuilder, const pat::Electron& electron, reco::TransientTrack& transTrack);

private:
  constexpr static double electronMass_ = 0.0005;
  edm::InputTag electronTag_, jetTag_;
  edm::InputTag goodPrimaryVertexTag_;

  unsigned int pdgId_;
  double rawMassMin_, rawMassMax_, massMin_, massMax_;

  double cut_minPt_, cut_maxEta_;
  double cut_trackChi2_, cut_trackSignif_, cut_DCA_;
  int cut_trackNHit_;
  double cut_vertexChi2_, cut_minLxy_, cut_maxLxy_, cut_vtxSignif_;

  unsigned int minNumber_, maxNumber_;

#ifdef DEBUGPLOT
  TH1F* hRawMass_, * hFitMass_;
  TH2F* hRawMassVsFitMass_;
#endif
};

KJpsiToElElProducer::KJpsiToElElProducer(const edm::ParameterSet& pset)
{
  electronTag_ = pset.getParameter<edm::InputTag>("src");
  jetTag_ = pset.getParameter<edm::InputTag>("jet");
  goodPrimaryVertexTag_ = pset.getParameter<edm::InputTag>("primaryVertex");

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
  produces<std::vector<double> >("jetDR");
  produces<std::vector<double> >("vProb");

#ifdef DEBUGPLOT
  edm::Service<TFileService> fs;
  hRawMass_ = fs->make<TH1F>("hRawMass", "raw mass;Raw mass (GeV/c^{2});Entries", 100, 0, 5);
  hFitMass_ = fs->make<TH1F>("hFitMass", "fit mass;Fit mass (GeV/c^{2});Entries", 100, 0, 5);
  hRawMassVsFitMass_ = fs->make<TH2F>("hRawMassVsFitMass", "raw vs fit;Raw mass (GeV/c^{2};Fit mass (GeV/c^{2}", 100, 0, 5, 100, 0, 5);
#endif
}

bool KJpsiToElElProducer::filter(edm::Event& event, const edm::EventSetup& eventSetup)
{
  using namespace reco;
  using namespace edm;
  using namespace std;

  typedef reco::VertexCompositeCandidateCollection VCCandColl;

  std::auto_ptr<VCCandColl> decayCands(new VCCandColl);
  std::auto_ptr<std::vector<double> > decayLengths(new std::vector<double>);
  std::auto_ptr<std::vector<double> > decayLengths3D(new std::vector<double>);
  std::auto_ptr<std::vector<double> > minJetDR(new std::vector<double>);
  std::auto_ptr<std::vector<double> > vProb(new std::vector<double>);

  edm::Handle< reco::VertexCollection >  goodPVHandle;
  event.getByLabel(goodPrimaryVertexTag_ , goodPVHandle);
  const reco::Vertex goodPV = goodPVHandle->at(0);

  const double pvx = goodPV.position().x();
  const double pvy = goodPV.position().y();
  const double pvz = goodPV.position().z();

  edm::ESHandle<TransientTrackBuilder> trackBuilder;
  eventSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilder);

  edm::Handle<std::vector<pat::Electron> > electronHandle;
  event.getByLabel(electronTag_, electronHandle);

  reco::TransientTrack transTrack1, transTrack2;
  for ( int i=0, n=electronHandle->size(); i<n; ++i )
  {
    const pat::Electron& electron1 = electronHandle->at(i);
    if ( electron1.charge() >= 0 ) continue;
    if ( !buildTransientTrack(trackBuilder, electron1, transTrack1) ) continue;
    if ( !transTrack1.impactPointTSCP().isValid() ) continue;
    FreeTrajectoryState ipState1 = transTrack1.impactPointTSCP().theState();

    for ( int j=0; j<n; ++j )
    {
      const pat::Electron& electron2 = electronHandle->at(j);
      if ( electron2.charge() <= 0 ) continue;
      if ( !buildTransientTrack(trackBuilder, electron2, transTrack2) ) continue;
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

      using namespace ROOT::Math;
      ROOT::Math::LorentzVector<PxPyPzE4D<double> > sum_vec = electron1.p4() + electron2.p4();
      const double rawMass = sum_vec.M();
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

      const double candE1 = hypot(mom1.mag(), electronMass_);
      const double candE2 = hypot(mom2.mag(), electronMass_);

      Particle::LorentzVector candLVec(mom.x(), mom.y(), mom.z(), candE1+candE2);
#ifdef DEBUGPLOT
      hFitMass_->Fill(candLVec.mass());
      hRawMassVsFitMass_->Fill(rawMass, candLVec.mass());
#endif
      if ( massMin_ > candLVec.mass() or massMax_ < candLVec.mass() ) continue;

      // Match to electrons
      pat::Electron newLep1(electron1);
      pat::Electron newLep2(electron2);
      newLep1.setP4(reco::Candidate::LorentzVector(mom1.x(), mom1.y(), mom1.z(), candE1));
      newLep2.setP4(reco::Candidate::LorentzVector(mom2.x(), mom2.y(), mom2.z(), candE2));
      VertexCompositeCandidate* cand = new VertexCompositeCandidate(0, candLVec, vtx, vtxCov, vtxChi2, vtxNdof);
      cand->addDaughter(newLep1);
      cand->addDaughter(newLep2);

      cand->setPdgId(pdgId_);
      AddFourMomenta addP4;
      addP4.set(*cand);

      edm::Handle<std::vector<pat::Jet> > jetHandle;
      event.getByLabel(jetTag_, jetHandle);
     
      double minDR = 9999.; 
      for ( int i=0, n=jetHandle->size(); i<n; ++i ) {
        Double_t deta = cand->eta() - jetHandle->at(i).eta();
        Double_t dphi = TVector2::Phi_mpi_pi( cand->phi()-jetHandle->at(i).phi());
        Double_t dr = TMath::Sqrt( deta*deta+dphi*dphi );
        if ( minDR > dr ) minDR = dr ;
      }

      minJetDR->push_back(minDR);
      vProb->push_back( TMath::Prob( vtxChi2, (int) vtxNdof));
      decayCands->push_back(*cand);
      decayLengths->push_back(rVtxMag);
      decayLengths3D->push_back(rVtxMag3D);

    }
  }

  const unsigned int nCands = decayCands->size();
  event.put(decayCands);
  event.put(vProb,"vProb");
  event.put(minJetDR,"jetDR");
  event.put(decayLengths, "lxy");
  event.put(decayLengths3D, "l3D");

  return (nCands >= minNumber_ and nCands <= maxNumber_);
}

bool KJpsiToElElProducer::buildTransientTrack(edm::ESHandle<TransientTrackBuilder>& trackBuilder,
                                              const pat::Electron& electron, reco::TransientTrack& transTrack)
{
  if ( electron.gsfTrack().isNull() ) return false;

  transTrack = trackBuilder->build(electron.gsfTrack());
  return true;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(KJpsiToElElProducer);
