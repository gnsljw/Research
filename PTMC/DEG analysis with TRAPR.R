source('TRAPR_Code.R')
ls()

PTMC <- TRAPR.Data.ReadExpressionTable('(1) PTMC(36)PTC(450) All gene symbol.txt', Exp1 = c(1:36), Exp2 = c(37:486), Tag = c('PTMC', 'PTC'))
str(PTMC)

TRAPR.DataVisualization(PTMC, 'box', logged = F)
TRAPR.DataVisualization(PTMC, 'DS', logged = F)
TRAPR.DataVisualization(PTMC, 'MA', logged = F)

PTMC <- TRAPR.Filter.ZeroValue(PTMC)# Filtering for zero values
PTMC <- TRAPR.Filter.LowVariance(PTMC) # Filtering for genes with low variance
PTMC <-  TRAPR.Normalize(PTMC, Method = 'UpperQuartile')

TRAPR.DataVisualization(PTMC, 'box', logged = F)
TRAPR.DataVisualization(PTMC, 'DS', logged = F)
TRAPR.DataVisualization(PTMC, 'MA', logged = F)

PTMC <- TRAPR.StatisticalTest(PTMC, Method = 'ttest', FDRControl = 'BH', PvalueThre = 0.05, FCThre = 0.5)
TestMatrix <- PTMC$CurrentMatrix[PTMC$DEGIndex,]
zTestMatrix <- (TestMatrix - rowMeans(TestMatrix)) / apply(TestMatrix, 1, var)
rownames(zTestMatrix) <- PTMC$DEGName
colnames(zTestMatrix) <- PTMC$SampleTag
heatmap(zTestMatrix)
TRAPR.ResultVisualization(PTMC, 'VO')
TRAPR.ResultVisualization(PTMC, 'HM')

TRAPR.Data.DEGResulttoFile(PTMC, FileName = '(1) TRAPR Result.txt')
TRAPR.Data.DEGNameListtoFile(PTMC)
