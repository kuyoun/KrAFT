import FWCore.ParameterSet.Config as cms

pdfWeight = cms.EDProducer("PDFWeightProducer",
    pdfName = cms.string("cteq66.LHgrid"),
    altPdfNames = cms.vstring([]),
)

