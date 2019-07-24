################################
### BASICS OF R
###

### The "#" indicates that something is a comment, nothing to the right of the # will be executed ###

print("hello") #THIS WILL NOT BE EXECUTED (though the print command will)

### Comments should be used liberally to document your code! ###

################################
## First things first - getting help
# R has a built-in help function that can bring up documentation on any function
?help  
?'+' 

# You can also search for any word (rather than just functions)
??help

################################
## Other places to get help
# http://cran.r-project.org/doc/manuals/r-release/R-intro.html
# http://www.cyclismo.org/tutorial/R/

################################
## R as a calculator

3^3
log10(100)
exp(4)

################################
# Creating objects and assigning values
# Objects are assigned using <- or =
# I'll be using <- 
# Note, R is case sensitive! a and A are different

x <- 5
y <- 2
x * y
test = "hello world"
print(test)


# There's a ton of great stuff out there, from interactive tools like CodeSchool.com 
# to text books and YouTube videos. If you want to explore more, I suggest trying a 
# few different methods to see which one fits your learning style best. The fastest
# way to learn R is to try analyzing your own data with it, though!

# In this brief lecture, we'll be mostly running line-by-line together, so don't
# worry if you don't immediately understand everything. Feel free to take notes
# as comments throughout this script, and check out the pdf documentation if you
# feel lost or want to revist anything.




################################
### MSSTATS
###

## The following analyses are also done in an R-markdown file with the same filename. ##


##
## 0. Prerequisites
##

# install dependency packages for MSstats from CRAN
install.packages(c("BiocManager", "reshape2","ggplot2"))

# install depency packages from BioConductor
BiocManager::install('MSstats')

# load MSstats
library(MSstats)

# use the '?' symbol to get documentation on MSstats
?MSstats

##
## 1. Preparing the data for MSstats input
##

# check the working directory
getwd()

# read in the Skyline export data
raw <- read.csv('MSstats_Input.csv')

# display the first few rows of the dataframe
head(raw)

# display the names of the columns
colnames(raw)

# display details about the dataframe
str(raw)

# format the dataframe so that it plays nicely with MSstats
raw_msstats <- SkylinetoMSstatsFormat(raw, filter_with_Qvalue = FALSE)

# check the details of the new object
str(raw_msstats)

# count how many proteins are in the data
unique(raw_msstats$ProteinName)

# count how many peptides are in the data
unique(raw_msstats$PeptideSequence)

# What else can you summarize about the dataframe?


##
## 2. Data processing with the `dataProcess` function
##


# see help page for more details
?dataProcess

# process the data 
quant_tmp <- dataProcess(raw = raw_msstats, 
                         normalization="globalStandards", 
                         nameStandards="VVLSGSDATLAYSAFK",
                         censoredInt = '0')

# display a summary of the dataProcess() output
str(quant_tmp)

# display just the column "ProcessedData"
head(quant_tmp$ProcessedData)

# display just the column "RunlevelData"
head(quant_tmp$RunlevelData)

# Show which summarization method was used (what do we expect this to say? hint: we mentioned it above!)
quant_tmp$SummaryMethod


### Using different summarization options 

# use the linear model for summarization
quant_nonorm <- dataProcess(raw = raw_msstats, 
                            normalization=FALSE,
                            censoredInt = '0')

# What's the difference between the linear and TMP summary methods? 
head(quant_tmp$RunlevelData)
head(quant_nonorm$RunlevelData)


## 2.2 Visualization of data processing

# display documentation on the dataProcessPlots() function
?dataProcessPlots

# generate QC plots
dataProcessPlots(data = quant_tmp, type = "QCplot", address = 'MSstats_')

# generate QC plots for all proteins only.
dataProcessPlots(data = quant_tmp, type = "QCplot", 
                 width = 7, height = 7,
                 which.Protein = 'allonly',
                 address = 'MSstats_')

# generate profile plots using the data summarized by TMP
dataProcessPlots(data = quant_tmp, type="Profileplot", 
                 width = 7, height = 7, address = "MSstats_")

# generate profile plots using the data with no normalization
dataProcessPlots(data = quant_nonorm, type="Profileplot", 
                 width = 7, height = 7, address = "MSstats_nonorm_")


# generate the profile plot for just one protein, NP_036620
dataProcessPlots(data = quant_tmp, type = "Profileplot", 
                 which.Protein = "NP_036620", 
                 width = 7, height = 7,
                 address = "NP_036620_")

# generate the condition plots
dataProcessPlots(data = quant_tmp, type = "conditionplot", 
                 width = 7, height - 7,
                 address = "MSstats_")

# display documentation for the groupComparison() function
?groupComparison

# check unique conditions and check order of condition information
levels(quant_tmp$ProcessedData$GROUP_ORIGINAL)

# create a contrast matrix for Diseased vs Healthy
comparison <- matrix(c(1, -1), nrow=1)
row.names(comparison) <- c("Diseased-Healthy")
comparison

# perform the group comparison
gpcomp_tmp <- groupComparison(contrast.matrix = comparison, data = quant_tmp)
names(gpcomp_tmp$ComparisonResult)
head(gpcomp_tmp$ComparisonResult)
head(gpcomp_tmp$ModelQC)
head(gpcomp_tmp$fittedmodel)

# pull just the results out of the whole group comparison output
gpcomp_res <- gpcomp_tmp$ComparisonResult

# subset only proteins with adjusted p-value < 0.05 and a FC > 2^2
list_sig <- gpcomp_res[gpcomp_res$adj.pvalue < 0.05 & abs(gpcomp_res$log2FC) > 2 , ]
head(list_sig)
nrow(list_sig)

##
## 3.2 Visualization of differentially abundant proteins
##

# display documnetation on groupComparisonPlots() function
?groupComparisonPlots

# generate a volcano plot of the analyzed data
groupComparisonPlots(data = gpcomp_tmp$ComparisonResult, 
                     type = 'VolcanoPlot',
                     sig = 0.05, FCcutoff = 2^2, 
                     address = 'MSstats_')

# generate a heatmap of the analyzed data
# ! Note that heatmaps can only be generated when multiple comparisons were made.
# In our experiment here, we only did one comparison, but I'm leaving this code here
# so that you can use it on your own data, if you need to!

#groupComparisonPlots(data = gpcomp_tmp$ComparisonResult, 
#                     type = 'Heatmap', 
#                     address = 'TMP_')


# generate a comparison plot
groupComparisonPlots(data = gpcomp_tmp$ComparisonResult,
                     type = 'ComparisonPlot', 
                     address = 'MSstats_')

##
## 4. Planning future experimental designs
##

# display documentation for the designSampleSize() function
?designSampleSize

# calculate the number of samples to achieve a range of fold changes from 1.1-1.5, at a fixed 90% power
design_size <- designSampleSize(data = gpcomp_tmp$fittedmodel,
                                 desiredFC = c(1.1, 1.5), FDR = 0.05,
                                 power = 0.9,
                                 numSample = TRUE)


# generate a plot from the sample size analysis results
designSampleSizePlots(data = design_size)


## 4.3. Calculating statistical power 

# calculate statistical power at various fold-change values with 7 replicates
design_power <- designSampleSize(data = gpcomp_tmp$fittedmodel,
                                  desiredFC = c(1.1, 1.5), FDR = 0.05,
                                  power = TRUE,
                                  numSample = 7)

# generate a plot from the statistical power analysis results
designSampleSizePlots(data = design_power)

##
## 5. Protein subject quantification 
##

# display documentation for the quantification() function
?quantification

# generate quantification results
sampleQuant <- quantification(quant_tmp)
head(sampleQuant)
