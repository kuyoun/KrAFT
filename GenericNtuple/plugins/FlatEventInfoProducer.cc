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

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include <memory>
#include <vector>
#include <string>

class FlatEventInfoProducer : public edm::EDProducer
{
public:
  FlatEventInfoProducer(const edm::ParameterSet& pset);
  void produce(edm::Event& event, const edm::EventSetup& eventSetup);
  void beginRun(edm::Run& run, const edm::EventSetup& eventSetup);

private:
  typedef std::vector<double> doubles;
  typedef std::vector<std::string> strings;

  edm::InputTag vertexLabel_;
  edm::InputTag genInfoLabel_;
  edm::InputTag hltLabel_;

  std::string processName_;
  std::map<std::string, strings> hltGroup_;
  HLTConfigProvider hltConfig_;

};

FlatEventInfoProducer::FlatEventInfoProducer(const edm::ParameterSet& pset)
{
  vertexLabel_ = pset.getParameter<edm::InputTag>("vertex");
  genInfoLabel_ = pset.getParameter<edm::InputTag>("genInfo");
  processName_ = pset.getParameter<std::string>("hltProcessName");

  produces<int>("pvN");
  produces<double>("pvX");
  produces<double>("pvY");
  produces<double>("pvZ");

  produces<int>("pdfId1");
  produces<int>("pdfId2");
  produces<double>("pdfX1");
  produces<double>("pdfX2");
  produces<double>("pdfQ");

  edm::ParameterSet hltSet = pset.getParameter<edm::ParameterSet>("HLT");
  for ( auto& hltSetName : hltSet.getParameterNamesForType<strings>() )
  {
    const std::string hltGroupName = hltSetName;
    strings& hltPaths = hltGroup_[hltGroupName];
    hltPaths = hltSet.getParameter<strings>(hltGroupName);
    for ( auto& hltPath : hltPaths )
    {
      hltPath = HLTConfigProvider::removeVersion(hltPath);
    }

    produces<int>("HLT"+hltSetName);
  }
}

void FlatEventInfoProducer::beginRun(edm::Run& run, const edm::EventSetup& eventSetup)
{
  bool changed = true;
  if ( !hltConfig_.init(run, eventSetup, processName_, changed) ) 
  {
    edm::LogError("FlatEventInfoProducer") << "HLT config extraction failure with process name " << processName_;
  }

  //if ( changed ) 
  //{
  //  // The HLT config has actually changed wrt the previous Run, hence rebook your
  //  // histograms or do anything else dependent on the revised HLT config
  //}
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

  edm::Handle<edm::TriggerResults> hltHandle;
  event.getByLabel(edm::InputTag("TriggerResults", "", processName_), hltHandle);
  for ( auto key = hltGroup_.begin(); key != hltGroup_.end(); ++key )
  {
    const std::string& hltGroupName = key->first;
    const strings& hltPaths = key->second;

    bool isPassed = false;
    for ( auto& hltPath : hltPaths )
    {
      const strings hltPathsWithV = HLTConfigProvider::restoreVersion(hltConfig_.triggerNames(), hltPath);
      if ( hltPathsWithV.empty() ) continue;
      const std::string& trigName = hltPathsWithV[0];
      const unsigned int trigIndex = hltConfig_.triggerIndex(trigName);
      if ( trigIndex < hltHandle->size() )
      {
        if ( hltHandle->accept(trigIndex) ) { isPassed = true; break; }
      }
      //const int psValue = hltConfig_.prescaleValue(event, eventSetup, trigName);
    }
    event.put(std::auto_ptr<int>(new int(isPassed)), "HLT"+hltGroupName);
  }

}

DEFINE_FWK_MODULE(FlatEventInfoProducer);
