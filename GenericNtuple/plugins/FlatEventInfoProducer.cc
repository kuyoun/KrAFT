#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

//#include "FWCore/Common/interface/TriggerNames.h"
//#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include <memory>
#include <vector>
#include <string>

class FlatEventInfoProducer : public edm::EDProducer
{
public:
  FlatEventInfoProducer(const edm::ParameterSet& pset);
  void produce(edm::Event& event, const edm::EventSetup& eventSetup);
  //void beginRun(edm::Event& event, const edm::EventSetup& eventSetup);

private:
  typedef std::vector<double> doubles;
  typedef std::vector<std::string> strings;

  edm::InputTag vertexLabel_;
  edm::InputTag genInfoLabel_;
  //strings hltPaths_ // TODO : implement HLT info
  //HLTConfigProducer hltConfig_;

};

FlatEventInfoProducer::FlatEventInfoProducer(const edm::ParameterSet& pset)
{
  vertexLabel_ = pset.getParameter<edm::InputTag>("vertex");
  genInfoLabel_ = pset.getParameter<edm::InputTag>("genInfo");
  //hltPaths_ = pset.getParameter<strings>("hltPaths");

  produces<int>("pvN");
  produces<double>("pvX");
  produces<double>("pvY");
  produces<double>("pvZ");

  produces<int>("pdfId1");
  produces<int>("pdfId2");
  produces<double>("pdfX1");
  produces<double>("pdfX2");
  produces<double>("pdfQ");

/*
  for ( auto& hltPath : hltPaths_ )
  {
    produces<bool>(hltPath);
  }
*/
}

void FlatEventInfoProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<reco::VertexCollection> vertexHandle;
  event.getByLabel(vertexLabel_, vertexHandle);

  const int nPV = vertexHandle->size();
  double pvX = 0, pvY = 0, pvZ = 0;
  if ( nPV > 0 )
  {
    pvX = vertexHandle->at(0).x();
    pvY = vertexHandle->at(0).y();
    pvZ = vertexHandle->at(0).z();
  }
  event.put(std::auto_ptr<int>(new int(nPV)), "pvN");
  event.put(std::auto_ptr<double>(new double(pvX)), "pvX");
  event.put(std::auto_ptr<double>(new double(pvY)), "pvY");
  event.put(std::auto_ptr<double>(new double(pvZ)), "pvZ");

  if ( !event.isRealData() )
  {
    edm::Handle<GenEventInfoProduct> genInfoHandle;
    event.getByLabel(genInfoLabel_, genInfoHandle);

    const float q = genInfoHandle->pdf()->scalePDF;
    const int id1 = genInfoHandle->pdf()->id.first;
    const int id2 = genInfoHandle->pdf()->id.second;
    const double x1 = genInfoHandle->pdf()->x.first;
    const double x2 = genInfoHandle->pdf()->x.second;

    event.put(std::auto_ptr<int>(new int(id1)), "pdfId1");
    event.put(std::auto_ptr<int>(new int(id2)), "pdfId2");
    event.put(std::auto_ptr<double>(new double(x1)), "pdfX1");
    event.put(std::auto_ptr<double>(new double(x2)), "pdfX2");
    event.put(std::auto_ptr<double>(new double(q)), "pdfQ");
  }

/*
  edm::Handle<edm::TriggerResults> hltHandle;
  event.getByLabel(edm::InputTag("TriggerResults"), hltHandle);
  for ( auto& hltPath : hltPaths_ )
  {
    
  }
*/

}

DEFINE_FWK_MODULE(FlatEventInfoProducer);
