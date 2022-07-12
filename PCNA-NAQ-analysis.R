
# # Install the packages ----------------------------------------------------

#Requirements
#For the mapping you will need to install R, bwa and samtools
#For the demultiplexing you will need R and the R ShortRead package
#For the analysis you will need R and the GenomicAlignments package



# Step 1: Get the functions -----------------------------------------------
#Select the following line and click Run:
source("/home/p.krijger_cbs-niob.local/projects/FM/FM_run1/submission/R/functions_PCNA_NAQ_220712.R")



# Step 2: Demultiplex the fastq files -------------------------------------
#This only needs to be done 1x for each experiment

fastqF<-"/home/p.krijger_cbs-niob.local/projects/FM/FM_run1/submission/R/fastq/"
info.file<-"/home/p.krijger_cbs-niob.local/projects/FM/FM_run1/submission/R/info.tsv"

#Run the following line to demultiplex and trim the fastq
demux(info.file, fastqF)



# Step 3: Map the data ----------------------------------------------------

plasmidFile<-'/home/p.krijger_cbs-niob.local/projects/FM/FM_run1/submission/R/bwa/FM_plasmids.fa'
bamFolder<- "/home/p.krijger_cbs-niob.local/projects/FM/FM_run1/submission/R/bam/"

#Make a bwa index based on the plasmid sequences (only need to be done 1x)
makeIndex(plasmidFile)

#Map the reads (only needs to be done 1x per experiment)
mapreads (ncores = 8
          ,index = plasmidFile
          ,fastqF = fastqF
          ,bamF=bamFolder
          ,info.file=info.file
          ,MQ=60) 




# Step 4: Define the bam files you want to analyze. -----------------------

#Define the results folder:
outF<-"/home/p.krijger_cbs-niob.local/projects/FM/FM_run1/submission/R/results/"

#automatically load all bam files stored in 1 folder:
bamFiles<-list.files(bamFolder, pattern = "\\.bam$", full.names = TRUE)


# Step 5: Extract the fragments from the bam files. -----------------------

fragList<-list()  
for (i in seq_along(bamFiles)){
  message("Extracting fragments from:",bamFiles[i])
  fragList[[i]] <- getFragments(bamFiles[i])
}

#You can filter the reads for mapping quality in case you did not already do this during the mapping
# fragListMQ<-list()  
# for (i in seq_along(bamFiles)){
#   message("Extracting fragments from:",bamFiles[i])
#   fragListMQ[[i]] <- getFragmentsV2(bamFiles[i], MQ=60)
#   
# }



# Step 6: Count the number of fragments identified for each plasmid -------
# ou can filter on fragment size

counts<-CountFragments(fragList,minLen=125,maxLen=160)
write.table(x = counts,file = paste0(outF,"fragCounts.txt",quote = FALSE, row.names = FALSE))


# Step 5: Plot the size distribution of the fragments. --------------------
# usage plotSize(fragments, plasmid=NULL/pRS415/pLoxHisCac3/207LoadingControl, plot="density/hist")

#plot size distribution
par(mfrow=c(2,4))
plotSize(fragments=fragList[[1]])
plotSize(fragList[[1]], plasmid="pRS415")
plotSize(fragList[[1]], plasmid="pLoxHisCac3")
plotSize(fragList[[1]], plasmid="207LoadingControl")

plotSize(fragments=fragList[[1]], plot = "hist")
plotSize(fragList[[1]], plasmid="pRS415", plot = "hist")
plotSize(fragList[[1]], plasmid="pLoxHisCac3", plot = "hist")
plotSize(fragList[[1]], plasmid="207LoadingControl", plot = "hist")



#make 1 big pdf for all experiments.  
filetype=pdf(paste0(outF,"fragsize.pdf"),
      width = 30/2.54,
      height = 20/2.54)

par(mfrow=c(4,4)) #4 columns and 4 rows
for(i in 1:length(fragList)){
  plotSize(fragments=fragList[[i]])
  plotSize(fragList[[i]], plasmid="pRS415")
  plotSize(fragList[[i]], plasmid="pLoxHisCac3")
  plotSize(fragList[[i]], plasmid="207LoadingControl")
 
}
dev.off()


# Step 6: Plot the coverage -----------------------------------------------

#Plot the normalized coverage for each plasmid after size selection
#usage: plotCov(fragments, minLen, maxLen,pRS415/pLoxHisCac3/207LoadingControl, col="black", yMax=NULL)


par(mfrow=c(2,3))
minLen = 125
maxLen = 160
plotCov(fragList[[1]], plasmid="pRS415", minLen=minLen, maxLen=maxLen)
plotCov(fragList[[1]], plasmid="pLoxHisCac3", minLen=minLen, maxLen=maxLen)
plotCov(fragList[[1]], plasmid="207LoadingControl", minLen=minLen, maxLen=maxLen)

minLen = 125
maxLen = 160
plotCov(fragList[[2]], plasmid="pRS415", minLen=minLen, maxLen=maxLen)
plotCov(fragList[[2]], plasmid="pLoxHisCac3", minLen=minLen, maxLen=maxLen)
plotCov(fragList[[2]], plasmid="207LoadingControl", minLen=minLen, maxLen=maxLen)


# make an overlay of 2 coverage plots ------------
par(mfrow=c(1,2))
plot2Cov(fragList[[1]],fragList[[2]], minLen=125, maxLen=160,plasmid="pRS415")
plot2Cov(fragList[[1]],fragList[[2]], minLen=125, maxLen=160,plasmid="pLoxHisCac3")

par(mfrow=c(2,2))
plot2Cov(fragList[[1]],fragList[[2]], minLen=125, maxLen=160,plasmid="pRS415")
plot2Cov(fragList[[1]],fragList[[2]], minLen=125, maxLen=160,plasmid="pLoxHisCac3")
plot2Cov(fragListMQ[[1]],fragListMQ[[2]], minLen=125, maxLen=160,plasmid="pRS415")
plot2Cov(fragListMQ[[1]],fragListMQ[[2]], minLen=125, maxLen=160,plasmid="pLoxHisCac3")




#Now make 1 big pdf with 4 rows, 3 columns

minLen = 125
maxLen = 160
filetype=pdf(paste0(outF,"coverage.pdf"),
             width = 20/2.54,
             height = 30/2.54)
par(mfrow=c(4,2))
for(i in 1:length(fragList)){
  message(i)
  plotCov(fragList[[i]], plasmid="pRS415", minLen=minLen, maxLen=maxLen)
  plotCov(fragList[[i]], plasmid="pLoxHisCac3", minLen=minLen, maxLen=maxLen)
}
dev.off()


#Or probally more usefull defining the yMax and only plotting the pRS415 plasmid
minLen = 125
maxLen = 160
filetype=pdf(paste0(outF,"coverage_pRS415.pdf"),
             width = 30/2.54,
             height = 20/2.54)
par(mfrow=c(4,4))
for(i in 1:length(fragList)){
    message(i)
    plotCov(fragList[[i]], plasmid="pRS415", minLen=minLen, maxLen=maxLen, yMax=300)
}
dev.off()



# Step 7: Export coverage to a csv file -----------------------------------
#Save coverage data as csv file to look at the data using Excel
minLen = 125
maxLen = 160
makeTable(fragList, outF, minLen, maxLen)  




 
