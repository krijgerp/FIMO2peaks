getOrientation<-function(peaks,fimo, colname="FIMO", strongestMotif=TRUE, overlappingMotif=TRUE){
  
  
  if(!is.data.frame(fimo)){
    message('fimo file is not a data frame, trying to read it as a file')
    fimo<-read.table(fimo,header=TRUE,sep="\t", stringsAsFactors = F)
  }
  
  if (!class(peaks) == "GRanges"){
    message('peaks file is not a GRanges object, trying to read it as a narrowPeak file')
    peaks<- import(peaks, format = "narrowPeak")
  }
  
  #FIMO sequence name should match peak name
  if(any(!fimo$sequence_name %in% peaks$name)){
    message('ERROR: fimo file contains sequence names that are not in peaks file')
    return(NULL)
  }
  
  
  #check for duplicates in peaks
  if(any(duplicated(peaks$name))){
    message('peaks file contains duplicate sequence names, using only the first of these')
    peaks<-peaks[!duplicated(peaks$name),]
  }
  
  if(strongestMotif){
    #strongest motif
    fimo_strong<-fimo[order(fimo$score,decreasing=TRUE),]
    fimo_strong<-fimo_strong[!duplicated(fimo_strong$sequence_name),]
    rownames(fimo_strong)<-as.character(fimo_strong$sequence_name)
    mcols(peaks)[paste0(colname,'_strongest')]<-fimo_strong[as.character(peaks$name),'strand']
  }
  
  if(overlappingMotif){
    #count overlapping motifs and keep unique motif
    freq.table<-as.data.frame.matrix(table(fimo$sequence_name, fimo$strand))
    colnames(freq.table)<-c("reverse","forward")
    freq.table$strand<-"*"
    freq.table$strand[ which( freq.table$forward > 0 & freq.table$reverse == 0)] <- "+"
    freq.table$strand[ which( freq.table$forward == 0 & freq.table$reverse > 0)] <- "-"
    freq.table$strand[ which( freq.table$forward > 0 & freq.table$reverse > 0)]  <- "both"
    message('Add strand information back to peak file')
    mcols(peaks)[colname]<-freq.table[as.character(peaks$name),'strand'] 
  }
  
  return(peaks)
}
