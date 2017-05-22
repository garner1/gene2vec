library(RTCGAToolbox)

# Valid aliases
getFirehoseDatasets()

# Valid stddata runs
stddata = getFirehoseRunningDates()
stddata

# Valid analysis running dates (will return 3 recent date)
gisticDate = getFirehoseAnalyzeDates(last=3)
gisticDate

# READ mutation data and clinical data 
brcaData = getFirehoseData (dataset="READ", runDate="20150402",forceDownload = TRUE,
                            Clinic=TRUE, Mutation=TRUE)

data(RTCGASample)
RTCGASample

# Differential gene expression analysis for gene level RNA data.
diffGeneExprs = getDiffExpressedGenes(dataObject=RTCGASample,DrawPlots=TRUE,
                                      adj.method="BH",adj.pval=0.05,raw.pval=0.05,
                                      logFC=2,hmTopUpN=10,hmTopDownN=10)

RTCGASampleCN = getData(RTCGASample,"GISTIC")
RTCGASampleCN[1:3,1:5]
