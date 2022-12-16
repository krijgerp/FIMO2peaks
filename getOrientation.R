getOrientation<-function(peaks,fimo, colname="FIMO_default", strongestMotif=TRUE, overlappingMotif=TRUE){

  #check for duplicates in peaks
  if(any(duplicated(peaks$name))){
    message('peaks file contains duplicate peak names, using only the first of these')
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
