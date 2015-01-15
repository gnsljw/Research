####################################################################################################
#######################################      TRAPR      ############################################
####################################################################################################

# By Lim Jae Hyun, Soo Youn Lee
# Last modified : 29-August-2012

######################### Main Function ############################################################
TRAPR <- function(File, sep = '\t', Exp1, Exp2, Method = 'ttest', FDRControl = 'none', Pvalue = 0.01)
{
	Data <- TRAPR.Data.ReadExpressionTable(File, sep='\t', Exp1, Exp2)
	Data <- TRAPR.Filter.ZeroValues(Data)
	Data <- TRAPR.Normalize.UpperQuartile(Data)
	Data <- TRAPR.StatisticalTest(Data, Method, FDRControl, Pvalue)
	TRAPR.ResultVisualization(Data, 'ALL')
}

######################### Data Manipulation ########################################################

TRAPR.Data.ReadExpressionTable <- function(File, sep = '\t', Exp1, Exp2, Tag = c('Exp1', 'Exp2'))
{
	cat('opening ', File, '\n')
	cat('It will take a while when your file size is large\n')
	RawData <- read.table(File, header=F, sep=sep)
	RawData.ColNum <- dim(RawData)[2]
	RawData.RowNum <- dim(RawData)[1]
	Data <- list()
	class(Data) <- 'TRAPR'

	if(is.vector(Exp1, mode='numeric')) Data$Exp1 <- Exp1 else stop('Exp1 must be numeric vector')
	if(is.vector(Exp2, mode='numeric')) Data$Exp2 <- Exp2 else stop('Exp2 must be numeric vector')
	
	Data$Tag <- Tag
	Data$SampleTag <- as.character(as.matrix(RawData[1, 2:RawData.ColNum]))
	Data$GeneTag <- as.character(as.matrix(RawData[2:RawData.RowNum, 1]))
	Data$RawMatrix <- matrix(as.numeric(as.matrix(RawData[2:RawData.RowNum, 2:RawData.ColNum])), ncol=RawData.ColNum - 1)
	
	Data$CurrentMatrix <- Data$RawMatrix
	Data$CurrentSample <- Data$SampleTag
	Data$CurrentGene <- Data$GeneTag
	
	Data$Exp1OnlyGene <- 'NA'
	Data$Exp2OnlyGene <- 'NA'
	Data$NonExpressedGene <- 'NA'
	
	Data$pvalues <- 'NA'
	Data$qvalues <- 'NA'
	Data$DEGName <- 'NA'
	Data$DEGIndex <- 'NA'
	Data$FC <- 'NA'
	
	cat(File, ' Loaded\n')
	return(Data)
}

TRAPR.Data.ReadGeneList <- function(File)
{
	cat('Reading ', File, ' as GeneList')
	temp <- read.table(File, header=FALSE)
	if(dim(temp)[2] > 1) stop('Invalid GeneList File')
	temp <- as.character(temp$V1)
	return(temp)
}

TRAPR.Data.ExpressionMatrixtoFile <- function(Data, FileName = 'output.txt')
{
	cat('Save Data$CurrentMatrix as ', FileName)
	TestMatrix <- Data$CurrentMatrix
	Header <- c('GeneName', Data$CurrentSample)
	TestMatrix <- cbind(Data$CurrentGene, TestMatrix)
	TestMatrix <- rbind(Header, TestMatrix)
	write.table(TestMatrix, file = FileName, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
}

TRAPR.Data.DEGNameListtoFile <- function(Data, FileNamePrefix = 'DEGList')
{	
	if(length(Data$DEGName) < 2) stop('DEG Not available')
	cat('Save DEGNames as', FileNamePrefix, '_up.txt which is upregulated in Exp1 and', FileNamePrefix, '_down.txt which is downregulated in Exp1')
	Up.Index <- which(Data$FC[Data$DEGIndex] > 0)
	Down.Index <- which(Data$FC[Data$DEGIndex] < 0)
	write(Data$DEGName[Up.Index], file = paste(FileNamePrefix, '_Up.txt', sep=""), sep = '\n')
	write(Data$DEGName[Down.Index], file = paste(FileNamePrefix, '_Down.txt', sep=""), sep = '\n')
}

TRAPR.Data.DEGResulttoFile <- function(Data, FileName = 'Result.txt')
{
	if(length(Data$DEGName) < 2) stop('DEG Not available')
	DEG.mean <- apply(Data$CurrentMatrix[Data$DEGIndex,], 1, mean)
	Result.Name <- c('GeneName', Data$DEGName)
	Result.FC <- c('FC', Data$FC[Data$DEGIndex])
	Result.pvalue <- c('pvalue', Data$pvalues[Data$DEGIndex])
	Result.qvalue <- c('qvalue', Data$qvalues[Data$DEGIndex])
	Result.RPKM <- c('Relative Expression Levels', DEG.mean)
	Result <- matrix(c(Result.Name, Result.FC, Result.pvalue, Result.qvalue, Result.RPKM), ncol = 5)
	write.table(Result, file = FileName, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
}

TRAPR.Data.ChangeExp <- function(Data, Exp1, Exp2)
{
	if(is.vector(Exp1, mode='numeric')) Data$Exp1 <- Exp1 else stop('Exp1 must be numeric vector')
	if(is.vector(Exp2, mode='numeric')) Data$Exp2 <- Exp2 else stop('Exp2 must be numeric vector')
	return(Data)
}

######################### Filters ########################################################################

TRAPR.Filter.ZeroValue <- function(Data)
{
	cat('Filter for Zero values\n')
	cat('')
	Exp1.Sum <- apply(Data$CurrentMatrix[,Data$Exp1], 1, sum)
	Exp2.Sum <- apply(Data$CurrentMatrix[,Data$Exp2], 1, sum)
	All.Sum <- apply(Data$CurrentMatrix, 1, sum)
	
	Exp1Zero.Index <- which(Exp1.Sum == 0)
	Exp2Zero.Index <- which(Exp2.Sum == 0)
	AllZero.Index <- which(All.Sum == 0)
	
	Exp1Only.Index <- setdiff(Exp2Zero.Index, AllZero.Index)
	Exp2Only.Index <- setdiff(Exp1Zero.Index, AllZero.Index)
	Whole.Index <- union(Exp1Zero.Index, Exp2Zero.Index)
	cat(as.character(length(AllZero.Index)), 'genes have zero values for all samples\n')
	cat(as.character(length(Exp2Only.Index)), 'genes have zero values for Exp1 samples\n')
	cat(as.character(length(Exp1Only.Index)), 'genes have zero values for Exp2 samples\n')
	
	if(length(Whole.Index) < 1) stop('Nothing to filter')
	Data$Exp1OnlyGene <- Data$CurrentGene[Exp1Only.Index]
	Data$Exp2OnlyGene <- Data$CurrentGene[Exp2Only.Index]
	Data$CurrentMatrix <- Data$CurrentMatrix[-Whole.Index,]
	Data$CurrentGene <- Data$CurrentGene[-Whole.Index]
	return(Data)
}

TRAPR.Filter.GeneList <- function(Data, GeneList)
{
	cat('Filter for Gene List\n')
	cat(as.character(length(GeneList)), 'genes are in your GeneList and', as.character(length(Data$CurrentGene)), 'genes are in your Data$CurrentGene\n')
	cat(as.character(length(Data$CurrentGene) - length(intersect(GeneList, Data$CurrentGene))), 'genes are filtered')
	
	temp <- c(1:length(Data$CurrentGene))
	names(temp) <- Data$CurrentGene
	as <- temp[intersect(Data$CurrentGene, CDGList)]
	
	Data$CurrentMatrix <- Data$CurrentMatrix[as,]
	Data$CurrentGene <- Data$CurrentGene[as]
	return(Data)
}

TRAPR.Filter.LowVariance <- function(Data, Thre = 0.1)
{
	cat('Filter for Low Variance Genes, it will stabilize your DEG.\n')
	cat(as.character(Thre * 100), 'Percent of your genes will be filtered, which is', as.character(length(Data$CurrentGene) * Thre))
	TestMatrix <- Data$CurrentMatrix
	TestMatrix.Var <- apply(TestMatrix, 1, var)
	VarThre = quantile(TestMatrix.Var, Thre)
	TestMatrix.In <- which(TestMatrix.Var > VarThre)
	
	Data$CurrentMatrix <- TestMatrix[TestMatrix.In,]
	Data$CurrentGene <- Data$CurrentGene[TestMatrix.In]
	return(Data)
}

TRAPR.Filter.LowExpression <- function(Data, Method = 'mean', Thre = 0.01)
{
	TestMatrix <- Data$CurrentMatrix

	cat('Filter for Low Expressed Genes, which might be read mismatch.\n')
	METHODS <- c('mean', 'min', 'max', 'median')
	Method <- pmatch(Method, METHODS)
	if (is.na(Method))
		stop('invalid expression estimation method')
	if (Method == -1)
		stop('ambiguous expression estimation method')
	if (Method == 1)
		TestMatrix.ex <- apply(TestMatrix, 1, mean)
	if (Method == 2)
		TestMatrix.ex <- apply(TestMatrix, 1, min)
	if (Method == 3)
		TestMatrix.ex <- apply(TestMatrix, 1, max)
	if (Method == 4)
		TestMatrix.ex <- apply(TestMatrix, 1, median)
	TestMatrix.In <- which(TestMatrix.ex > Thre)
	cat(as.character(length(TestMatrix.In)), 'genes are filtered')
	
	Data$CurrentMatrix <- TestMatrix[TestMatrix.In,]
	Data$CurrentGene <- Data$CurrentGene[TestMatrix.In]
	return(Data)
}

TRAPR.Filter.SampleDeletion <- function(Data, Outlier)
{
	Data$CurrentMatrix <- Data$CurrentMatrix[,-Outlier]
	Data$CurrentSample <- Data$CurrentSample[-Outlier]
	Data$Exp1 <- Data$Exp1[-intersect(Data$Exp1, Outlier)]
	Data$Exp2 <- Data$Exp2[-intersect(Data$Exp2, Outlier)]
	return(Data)
}

TRAPR.Filter.GeneDeletion <- function(Data, Outlier)
{
	Data$CurrentMatrix <- Data$CurrentMatrix[-Outlier]
	Data$CurrentGene <- Data$CurrentGene[-Outlier]
	return(Data)
}

######################### Transformation #################################################
TRAPR.Transformation.VSN <- function(Data)
{
	cat('Variance Stabilization Normalization\n')
	library(vsn)
	Data$CurrentMatrix <- vsn2(Data$CurrentMatrix)@hx
	return(Data)
}

TRAPR.Transformation.log2 <- function(Data)
{
	cat('log2 transformation')
	Data$CurrentMatrix <- log2(Data$CurrentMatrix + unique(sort(Data$CurrentMatrix))[2])
	Data$CurrentMatrix <- Data$CurrentMatrix
	return(Data)
}

######################### Normalization ##################################################
TRAPR.Normalize.UpperQuartile <- function(Data)
{
	TestMatrix <- Data$CurrentMatrix
	TestMatrix.Sums <- apply(TestMatrix, 1, sum)
	TestMatrix.NotZero <- which(TestMatrix.Sums != 0)
	TestMatrix.UpperQuartile <- apply(TestMatrix[TestMatrix.NotZero,], 2, function(x) quantile(x, 0.75))
	TestMatrix.UpperQuartile <- TestMatrix.UpperQuartile / mean(TestMatrix.UpperQuartile)
	nMatrix <- TestMatrix[,1] / TestMatrix.UpperQuartile[1]
	for(i in 2:dim(TestMatrix)[2])
	{
		nMatrix <- c(nMatrix, TestMatrix[,i]/TestMatrix.UpperQuartile[i])
	}
	nMatrix <- matrix(nMatrix, ncol=dim(TestMatrix)[2])
	Data$CurrentMatrix <- nMatrix
	return(Data)
}

TRAPR.Normalize.Quantile <- function(Data)
{
	library(preprocessCore)
	Data$CurrentMatrix <- normalize.quantiles(Data$CurrentMatrix)
	return(Data)
}

TRAPR.Normalize.Median <- function(Data)
{
	TestMatrix <- Data$CurrentMatrix
	TestMatrix.Sums <- apply(TestMatrix, 1, sum)
	TestMatrix.NotZero <- which(TestMatrix.Sums != 0)
	TestMatrix.Median <- apply(TestMatrix[TestMatrix.NotZero,], 2, median)
	TestMatrix.Median <- TestMatrix.Median / mean(TestMatrix.Median)
	nMatrix <- TestMatrix[,1] / TestMatrix.Median[1]
	for(i in 2:dim(TestMatrix)[2])
	{
		nMatrix <- c(nMatrix, TestMatrix[,i]/TestMatrix.Median[i])
	}
	nMatrix <- matrix(nMatrix, ncol=dim(TestMatrix)[2])
	Data$CurrentMatrix <- nMatrix
	return(Data)	
}

TRAPR.Normalize.Mean <- function(Data)
{
	TestMatrix <- Data$CurrentMatrix
	TestMatrix.Sums <- apply(TestMatrix, 1, sum)
	TestMatrix.NotZero <- which(TestMatrix.Sums != 0)
	TestMatrix.Mean <- apply(TestMatrix[TestMatrix.NotZero,], 2, mean)
	TestMatrix.Mean <- TestMatrix.Mean / mean(TestMatrix.Mean)
	nMatrix <- TestMatrix[,1] / TestMatrix.Mean[1]
	for(i in 2:dim(TestMatrix)[2])
	{
		nMatrix <- c(nMatrix, TestMatrix[,i]/TestMatrix.Mean[i])
	}
	nMatrix <- matrix(nMatrix, ncol=dim(TestMatrix)[2])
	Data$CurrentMatrix <- nMatrix
	return(Data)
}

TRAPR.Normalize <- function(Data, Method = 'UpperQuartile')
{
	cat(Method, 'Normalization \n')
	METHODS <- c('UpperQuartile', 'Quantile', 'Median', 'Mean')
	Method <- pmatch(Method, METHODS)
	if (is.na(Method))
		stop('invalid statistical test method')
	if (Method == -1)
		stop('ambiguous statistical test method')
	if (Method == 1)
		Data <- TRAPR.Normalize.UpperQuartile(Data)
	if (Method == 2)
		Data <- TRAPR.Normalize.Quantile(Data)
	if (Method == 3)
		Data <- TRAPR.Normalize.Median(Data)
	if (Method == 4)
		Data <- TRAPR.Normalize.Mean(Data)
	return(Data)
}

######################### Statistical Test ###############################################
TRAPR.StatisticalTest.ttest <- function(Data)
{
	Data$pvalues <- apply(Data$CurrentMatrix, 1, function(x) as.numeric(t.test(x[Data$Exp1], x[Data$Exp2], var.equal=FALSE)$p.value))
	return(Data)
}

TRAPR.StatisticalTest.Wilcoxon <- function(Data)
{
	Data$pvalues <- apply(Data$CurrentMatrix, 1, function(x) as.numeric(wilcox.test(x[Data$Exp1], x[Data$Exp2], var.equal=FALSE)$p.value))
	return(Data)
}

TRAPR.StatisticalTest.FC <- function(Data)
{
	Data$pvalues <- c(rep('NA', dim(Data$CurrentMatrix)[1]))
	return(Data)
}

TRAPR.StatisticalTest.EdgeR <- function(Data)
{
	library('edgeR')
	d <- Data$CurrentMatrix
	if(min(d) < 0) stop('EdgeR needs non-negative matrix')
	rownames(d) <- Data$CurrentGene
	group <- c(rep('Exp1', length(Data$Exp1)), rep('Exp2', length(Data$Exp2)))
	d <- DGEList(counts = d, group = group)
	d <- calcNormFactors(d)
	d <- estimateCommonDisp(d)
	de.com <- exactTest(d)
	Data$FC <- de.com$table$logFC
	Data$pvalues <- de.com$table$PValue
	return(Data)
}

TRAPR.StatisticalTest <- function(Data, Method = 'ttest', FDRControl = 'none', PvalueThre = 0.01, FCThre = 0) 
{
	cat('Statistical Analysis by', Method, 'FDR Threshould for selecting DEG is \'FDR <', as.character(PvalueThre), '\', \' Fold Change(log2(Exp1 / Exp2)) >', as.character(FCThre), '\'\n')
	temp <- c(1:length(Data$GeneTag))
	names(temp) <- Data$GeneTag
	FCIndex <- temp[Data$CurrentGene]

	Data$FC <- log2(apply(Data$RawMatrix[FCIndex,Data$Exp1], 1, mean) / apply(Data$RawMatrix[FCIndex,Data$Exp2], 1, mean))

	METHODS <- c('ttest', 'wilcoxon', 'edgeR', 'FC')
	Method <- pmatch(Method, METHODS)
	if (is.na(Method))
		stop('invalid statistical test method')
	if (Method == -1)
		stop('ambiguous statistical test method')
	if (Method == 1)
		Data <- TRAPR.StatisticalTest.ttest(Data)
	if (Method == 2)
		Data <- TRAPR.StatisticalTest.Wilcoxon(Data)
	if (Method == 3)
		Data <- TRAPR.StatisticalTest.EdgeR(Data)
	if (Method == 4)
		Data <- TRAPR.StatisticalTest.FC(Data)
	
	FDRMETHODS <- c('holm', 'hotchberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none')
	Test <- pmatch(FDRControl, FDRMETHODS)
	if (is.na(Test)) stop('invalid FDR Control Method')
	if (Method < 0) stop('ambiguous FDR Control Method')
	if (Method > 0) Data$qvalues <- as.numeric(p.adjust(Data$pvalues, method = FDRControl))
	PvalIndex <- which(Data$qvalues < PvalueThre)
	FCIndex <- which(abs(Data$FC) > FCThre)
	Data$DEGIndex <- intersect(PvalIndex, FCIndex)
	cat(as.character(length(Data$DEGIndex)), ' genes selected \n')
	Data$DEGName <- Data$CurrentGene[Data$DEGIndex]
	return(Data)
}

TRAPR.StatisticalTest.ReThreshould <- function(Data, PvalueThre = 0.01, FCThre = 0)
{
	if(Data$qvalues == 'NA') stop('You need to estimate FC and qvalues by StatisticalTest Function, first')
	PvalIndex <- which(Data$qvalues < PvalueThre)
	FCIndex <- which(Data$FC < FCThre)
	Data$DEGIndex <- intersect(PvalIndex, FCIndex)
	cat(as.character(length(Data$DEGIndex)), 'genes selected')
	Data$DEGName <- Data$CurrentGene[Data$DEGIndex]
	return(Data)
}

######################### Visualization ###########################################

TRAPR.DataVisualization <- function(Data, type_of_plot,logged = TRUE)
{
	cat('Visualization for Data$CurrentMatrix\n')
	library(gridExtra)
	cat('Loading gridExtra package \n')
	library(ggplot2)
	cat('Loading ggplot2 package \n')
	library(reshape2)
	cat('Loading reshape2 package \n')
	
	if (type_of_plot != 'box' && type_of_plot != 'MA'&& type_of_plot != 'MV' && type_of_plot != 'SE'&& type_of_plot != 'DS'&& type_of_plot != 'ALL')
	stop('err please check plot names !!!! \n')
	
	
	#Data <- COPD
	TestMatrix <- Data$CurrentMatrix
	rownames(TestMatrix) <- Data$CurrentGene
	colnames(TestMatrix) <- Data$CurrentSample
	sample_level = Data$Tag
	
	#'Mean' and 'Variance'
	TestMatrix.mean <- (apply(TestMatrix, 1, mean))
	TestMatrix.var <- (apply(TestMatrix, 1, var))
	Mean = TestMatrix.mean
	Variance = TestMatrix.var
		
	#'MA'
	Exp1.Mean <- apply(TestMatrix[ ,Data$Exp1], 1, mean)
	Exp2.Mean <- apply(TestMatrix[ ,Data$Exp2], 1, mean)
	
	
	#Scatter plot for Exp1 vs Exp2
	Exp1.Mean <- apply(TestMatrix[ ,Data$Exp1], 1, mean)
	Exp2.Mean <- apply(TestMatrix[ ,Data$Exp2], 1, mean)
	
	##################################### plot functions########################################

	Boxplot.data <- function()
	{
		tt = melt(TestMatrix)
		group = c( rep(sample_level[1], length(Data$Exp1)*dim(TestMatrix )[1]), rep(sample_level[2], length(Data$Exp2) * dim(TestMatrix)[1]) )
		Data_value = data.frame(tt,group)
		colnames(Data_value) =c("gene","sample","value","category")
		return(Data_value)
	}
	Boxplot <- function(input_data)
	{	
		p = ggplot(input_data, aes(sample, value))
		p = p + geom_boxplot(aes(fill=category ,colour=category ))
		p = p + stat_summary(fun.y=median, geom="crossbar", ymin=-1, ymax=1, color="grey90")
		p = p + labs(x = "sample", y = "RPKM")
		#p = p + theme(plot.margin = unit(c(0,0,1,0), "cm"))
		#p = p + theme(title="boxplot")
		return(p)
	}
	
	MA <- function(M,A)
	{
		MA_TestDataframe <- data.frame(M,A)
		p = ggplot(MA_TestDataframe, aes(A, M)) 
		p = p + geom_point(shape=8, colour = "#FF3399")
		#p = p + geom_hline(aes(yintercept=0),colour = "#6633CC",size = 1)
		#p = p + opts(title="MA plot for Exp1 vs Exp2")
		p = p + geom_smooth(method="loess")
		return(p)
	}
	
	Mean_Variance.data <- function()
	{
		MV_TestDataframe <- data.frame(Mean,Variance )
		return(MV_TestDataframe)
	}
	
	Mean_Variance <- function(input_data)
	{
		p = ggplot(input_data, aes(Mean, Variance)) 
		p = p + geom_point(shape=8, colour = "#FFFF00")
		p = p +  geom_smooth(method=lm,colour="#66FF00") + scale_colour_hue(l=50)
		#p = p + opts(title="Mean vs Variance")
		return(p)

	}
	
	Scatter_plot_Exp1_Exp2.data <- function()
	{
		SE_TestDataframe <- data.frame(Exp1.Mean,Exp2.Mean )
		return(SE_TestDataframe)
	}
	
	Scatter_plot_Exp1_Exp2 <- function(input_data)
	{
		p = ggplot(input_data, aes(Exp1.Mean, Exp2.Mean)) 
		p = p + geom_point(shape=8, colour = "#0000FF")
		#p = p + geom_smooth(method = lm, colour = "#FF33CC",fullrange=T,size = 1)
		p = p + geom_smooth(method="loess", colour = "#FF33CC",fullrange=T,size = 1)
		#plot(p)
		#p = p + opts(title="Scatter plot for Exp1 vs Exp2")
		return(p)
	}
	
	Density_plot.data <- function()
	{
		DS_TestDataframe <- data.frame(TestMatrix)
		m.d = melt(DS_TestDataframe)
		return(m.d)
		#### before Normalization ####
		#a = data.frame(sample$CurrentMatrix)
		#b = melt(a)
		#ggplot(b, aes(value, colour = variable)) +  geom_density()
	}
	
	Density_plot <- function(input_data,les)
	{
		p = ggplot(input_data, aes(value, colour = variable)) 
		p = p + geom_density(size=0.5)
		#p = p + opts(legend.key.size = unit(les, "cm"))
		#p = p + opts(title="Density Plot")
		return(p)

	}
	#####################################################################################
	##################################### boxplot ########################################
	if(type_of_plot == 'box') 
		{
		if (logged) 
			{
				in_data = Boxplot.data()
				plot(Boxplot(in_data))				
			}
		else
			{
				in_data = Boxplot.data()
				in_data$value = log2(in_data$value)
				plot(Boxplot(in_data))
			}
		}
	#################################################################################
	##################################### MA ########################################
	
	if(type_of_plot == 'MA') 
		{
		if (logged)
			{
			#'MA'
				M <- Exp1.Mean - Exp2.Mean
				A <- (Exp1.Mean + Exp2.Mean)/2
				plot(MA(M,A))
			}
		else
			{
				#'MA'
				ExpIndex <- intersect(which(Exp1.Mean > 0), which(Exp2.Mean > 0))
				M <- log2(Exp1.Mean[ExpIndex]) - log2(Exp2.Mean[ExpIndex])
				A <- log2(Exp1.Mean[ExpIndex]) + log2(Exp2.Mean[ExpIndex])
				plot(MA(M,A))
			}
			
		}
	##################################################################################################
	##################################### Mean' and 'Variance ########################################
if(type_of_plot == 'MV') 
		{
		if (logged)
			{							
				in_data = Mean_Variance.data()
				#print(head(in_data))
				plot(Mean_Variance(in_data)) 
			}
		
		else
			{
				in_data = log2(Mean_Variance.data())
				plot(Mean_Variance(in_data)) 
			}
		}

	##################################################################################################
	##################################### Scatter plot for Exp1 vs Exp2 ##############################
if(type_of_plot == 'SE') 
		{
		if (logged)
			{							
				in_data = Scatter_plot_Exp1_Exp2.data()
				plot(Scatter_plot_Exp1_Exp2(in_data)) 
			}
		
		else
			{
				in_data = log2(Scatter_plot_Exp1_Exp2.data())
				plot(Scatter_plot_Exp1_Exp2(in_data))
			}
		}	
	##################################################################################################
	##################################### density plot ##############################
if(type_of_plot == 'DS') 
		{
		if (logged)
			{							
				in_data = Density_plot.data()
				plot(Density_plot(in_data,0.6))	
			}		
		else
			{
				in_data = Density_plot.data()
				in_data$value = log2(in_data$value)
				plot(Density_plot(in_data,0.6))
			}
		}		
	##################################################################################################
	##################################### ALL ########################################################
	if(type_of_plot == 'ALL') 
		{
			if (logged)
			{
				box_in_data = Boxplot.data()
				
				M <- Exp1.Mean - Exp2.Mean
				A <- (Exp1.Mean + Exp2.Mean)/2
				
				MV_in_data = Mean_Variance.data()
				
				SE_in_data = Scatter_plot_Exp1_Exp2.data()
	
				DS_in_data = Density_plot.data()
				
				grid.arrange(MA(M,A),Mean_Variance(MV_in_data),Scatter_plot_Exp1_Exp2(SE_in_data),Density_plot(DS_in_data,0.3),Boxplot(box_in_data),ncol=3)
				#ggsave(p, file="Genre_04_04.PNG", width=10, height=10)
		}
	else 
		{
				box_in_data = Boxplot.data()
				box_in_data$value = log2(box_in_data$value)
				
				ExpIndex <- intersect(which(Exp1.Mean > 0), which(Exp2.Mean > 0))
				M <- log2(Exp1.Mean[ExpIndex]) - log2(Exp2.Mean[ExpIndex])
				A <- log2(Exp1.Mean[ExpIndex]) + log2(Exp2.Mean[ExpIndex])
				
				MV_in_data = log2(Mean_Variance.data())

				SE_in_data = log2(Scatter_plot_Exp1_Exp2.data())
				
				DS_in_data = Density_plot.data()
				DS_in_data$value = log2(DS_in_data$value)
				
				grid.arrange( MA(M,A),Mean_Variance(MV_in_data),Scatter_plot_Exp1_Exp2(SE_in_data),Density_plot(DS_in_data,0.3),Boxplot(box_in_data),ncol=3)
				#ggsave(p, file="Genre_04_04.PNG", width=10, height=10)
		}
	}
	
}
#TRRAP.ResultVisualization(COPD,1.2,0.01,"VO")
#source("TRAPR_Code-2.r")	
	
TRAPR.ResultVisualization <- function(Data,type_of_plot)
{
	if (type_of_plot != 'VO' && type_of_plot != 'HM'& type_of_plot != 'ALL')
	stop('err please check plot names !!!! \n')
	
	library(gridExtra)
	cat('Loading gridExtra package \n')
	library(gplots)
	cat('Loading gplots package \n')
	library(ggplot2)
	cat('Loading ggplot2 package \n')
	
	if(length(Data$DEGName) < 2) stop('DEG not selected')

	TestMatrix <- Data$CurrentMatrix
	colnames(TestMatrix) <- Data$CurrentSample
	rownames(TestMatrix) <- Data$CurrentGene

	x.min <- quantile(Data$FC, probs = 0.01, na.rm=T)
	x.max <- quantile(Data$FC, probs = 0.99, na.rm=T)
	x.lim <- max(abs(x.min), abs(x.max))
	y.max <- quantile(-log10(Data$pvalue), probs=0.9999, na.rm=T)
	y.lim <- abs(y.max)

	if(type_of_plot == 'VO') 
	{	
		deg_index = rep(0,dim(TestMatrix)[1])
		deg_index[Data$DEGIndex]=1
		threshold = as.factor(deg_index)
		#threshold = as.factor(abs(Data$FC) > FC1 & Data$qvalue < FDR)
		Plot_matrix = cbind(TestMatrix,threshold,"FC" = Data$FC,"pvalues" = -log10(Data$pvalues))
		VO_TestDataframe <- data.frame(Plot_matrix)
		
		p = ggplot(VO_TestDataframe, aes(FC,pvalues, colour=threshold))
		p  = p + geom_point(size=1.75) #alpha=0.4,
		p  = p +  xlim(c(-x.lim, x.lim)) + ylim(c(0, y.lim)) 
		#p  = p + opts(legend.position = "none")
		p  = p + xlab("log2 fold change") + ylab("-log10 p-value")
		plot(p)
	}
	
	if(type_of_plot == 'HM') 
	{
		HM_TestMatrix <- TestMatrix[Data$DEGIndex,]
		rc = rainbow(nrow(HM_TestMatrix),start=0,end=.3)
		cc <- rainbow(ncol(HM_TestMatrix), start=0, end=.3)
		heatmap.2(HM_TestMatrix,trace="none",col=topo.colors(16),hclustfun = hclust,distfun = dist, scale ='row', RowSideColors=rc, ColSideColors=cc,main="HeatMap")
		#topo.colors(16),cm.colors(255)
	}
	
	if(type_of_plot == 'ALL') 
	{
		deg_index = rep(0,dim(TestMatrix)[1])
		deg_index[Data$DEGIndex]=1
		threshold = as.factor(deg_index)
		#threshold = as.factor(abs(Data$FC) > FC1 & Data$qvalue < FDR)
		Plot_matrix = cbind(TestMatrix,threshold,"FC" = Data$FC,"pvalues" = -log10(Data$pvalues))
		VO_TestDataframe <- data.frame(Plot_matrix)
		
		p = ggplot(VO_TestDataframe, aes(FC,pvalues, colour=threshold))
		p  = p + geom_point(size=1.75) #alpha=0.4,
		p  = p +  xlim(c(-x.lim, x.lim)) + ylim(c(0, y.lim)) 
		#p  = p + opts(legend.position = "none")
		p  = p + xlab("log2 fold change") + ylab("-log10 p-value")
		
		
		HM_TestMatrix <- TestMatrix[Data$DEGIndex,]
		rc = rainbow(nrow(HM_TestMatrix),start=0,end=.3)
		cc <- rainbow(ncol(HM_TestMatrix), start=0, end=.3)
		dev.new()
		heatmap.2(HM_TestMatrix,trace="none",col=topo.colors(16),hclustfun = hclust,distfun = dist, scale ='row', RowSideColors=rc, ColSideColors=cc,main="HeatMap")
		dev.new()
		plot(p)


	}
	
}