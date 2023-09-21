plotCTCFmotif <- function(peaks, zoom, y1 = -5, y2 = -10, strandCol = "FIMO", Mb=TRUE) {
  
  if (!class(peaks) == "GRanges"){
    message('ERROR: peaks file is not a GRanges object')
    return(NULL)
  }
  
  if (!strandCol %in% names(mcols(peaks))){
    stop(paste('ERROR: peaks file does not contain column',strandCol))
    
  }
  
  
  x.wid <- 2.5e-08 * width(zoom)
  
  if(Mb){
    Xdiv<-1e6
  }
else{
  Xdiv<-1
}
  
  MOI <- subsetByOverlaps(peaks, zoom)
  MOI_withMotif<-MOI[as.vector(!is.na(mcols(MOI)[paste0(strandCol)]))]
  
  
  # plot positive strands
  positive_strand_peaks <- MOI_withMotif[mcols(MOI_withMotif)[[strandCol]] == "+"]
  
 if(length(positive_strand_peaks)>0) {
    triangle <- lapply(start(positive_strand_peaks) / Xdiv, function(x) list(xx = c(x, x, x + x.wid),
                                                                            yy = c(y1, y2, (y1 + y2) / 2)))
    lapply(triangle, function(x) polygon(x$xx, x$yy, col = "red", border = "red"))
  }
  
  # plot negative strands
  negative_strand_peaks <- MOI_withMotif[mcols(MOI_withMotif)[[strandCol]] == "-"]
  if (length(negative_strand_peaks)>0) {
    triangle <- lapply(start(negative_strand_peaks) / Xdiv, function(x) list(xx = c(x, x, x - x.wid),
                                                                            yy = c(y1, y2, (y1 + y2) / 2)))
    lapply(triangle, function(x) polygon(x$xx, x$yy, col = "blue", border = "blue"))
  }
  
  # plot unorientated motifs as grey rectangles
  both_strand_peaks <- MOI_withMotif[mcols(MOI_withMotif)[[strandCol]] == "both"]
  if (length(both_strand_peaks)>0) {
    rect(start(both_strand_peaks) / Xdiv, y1, end(both_strand_peaks) / Xdiv, y2, col = "lightgrey", border = "lightgrey")
  }
  
  # plot strands without a motif as black rectangles
  noMotif_peaks <- MOI[is.na(mcols(MOI)[[strandCol]])]
  if(length(noMotif_peaks)>0) {
    rect(start(noMotif_peaks) / Xdiv, y1, end(noMotif_peaks) / Xdiv, y2, col = "black", border = "black")
  }
}
