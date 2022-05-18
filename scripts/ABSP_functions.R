
### FUNCTIONS



# FIND CG COORDINATES FUNCTION --------------------------------------------

findCGcoordinates <- function(sequence, sequence_reverse, chr, coord_start, coord_end, seq_strand) {
  # Find CG positions on unconverted sequences (both strands)
  CG_pos <- matchPattern("CG", sequence, max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE, algorithm="auto")
  CG_pos_RC <- matchPattern("CG", sequence_reverse, max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE, algorithm="auto")
  CG_pos <- as.data.frame(ranges(CG_pos))[,c(1,2)] ; colnames(CG_pos) <- c("C_pos","G_pos")
  CG_pos_RC <- as.data.frame(ranges(CG_pos_RC))[,c(1,2)] ; colnames(CG_pos_RC) <- c("C_pos_RC","G_pos_RC")
  # Find CG coordinates depending if base sequence in from plus or minus strand (primer design)
  if(seq_strand=="plus") {
    CG_pos <- CG_pos %>% 
      dplyr::mutate(CG_nb=dplyr::row_number(), chr=chr, start= coord_start + C_pos -1, end= coord_start + G_pos -1) %>% 
      .[,c(3,4,5,6,1,2)]
    CG_pos_RC <- CG_pos_RC %>% 
      dplyr::mutate(CG_nb=rev(dplyr::row_number()), chr=chr, start= coord_end - G_pos_RC +1, end= coord_end - C_pos_RC +1) %>% 
      .[,c(3,4,5,6,1,2)]    
  }
  if(seq_strand=="minus") {
    CG_pos <- CG_pos %>% 
      dplyr::mutate(CG_nb=rev(dplyr::row_number()), chr=chr, start= coord_end - G_pos +1, end= coord_end - C_pos +1) %>% 
      .[,c(3,4,5,6,1,2)]
    CG_pos_RC <- CG_pos_RC %>% 
      dplyr::mutate(CG_nb=dplyr::row_number(), chr=chr, start= coord_start + C_pos_RC -1, end= coord_start + G_pos_RC -1) %>% 
      .[,c(3,4,5,6,1,2)]
  }
  CG_coord <- list("CG_pos"=CG_pos,"CG_pos_RC"=CG_pos_RC)
  return(CG_coord)
}



# BISULFITE CONVERSION FUNCTION -------------------------------------------

bisulfiteConversion <- function(sequence) {
  sequence_BS <- str_replace_all(as.character(sequence),"CG","YG") %>% 
    str_replace_all(.,"C","T") %>% 
    str_replace_all(.,"Y","C")
  sequence_BS <- DNAString(sequence_BS)
  sequence_BS_RC <- reverseComplement(sequence_BS)
  seq_BS <- list()
  seq_BS[["seq_BS"]] <- sequence_BS
  seq_BS[["seq_BS_RC"]] <- sequence_BS_RC
  return(seq_BS)
}





# PEAK TABLE OF SEQUENCING ------------------------------------------------

seqtable <- function(seq) {
  # Retrieve Peak Amplitude Matrix A C G T
  signals <- as.data.frame(peakAmpMatrix(seq)) 
  # Rename columns
  colnames(signals) <- c("peak_A","peak_C","peak_G","peak_T")
  # Add position numbers (row numbers)
  signals <- signals %>%
    dplyr::mutate(position = dplyr::row_number())
  # Rearrange columns
  signals <- signals[,c(5,1,4,3,2)]
  # Compute new columns
  signals <- signals %>%
    dplyr::mutate(total=peak_A+peak_C+peak_G+peak_T) %>% # Total of peak values for all bases at each position
    dplyr::mutate(ratio_A=round(peak_A/total,2)) %>% # ratio of A peaks over total peaks at each position
    dplyr::mutate(ratio_T=round(peak_T/total,2)) %>% # ratio of T peaks over total peaks at each position
    dplyr::mutate(ratio_G=round(peak_G/total,2)) %>% # ratio of G peaks over total peaks at each position
    dplyr::mutate(ratio_C=round(peak_C/total,2)) %>% # ratio of C peaks over total peaks at each position
    dplyr::mutate(primary_base=strsplit(primarySeq(seq,string = TRUE),"")[[1]]) %>%
    dplyr::mutate(quality_score=seq@QualityReport@qualityPhredScores) %>%
    rowwise() %>%
    dplyr::mutate(primary_ratio=max(ratio_A,ratio_T,ratio_G,ratio_C)) # ratio of the main peak
  # Output table
  return(signals)
}




# MIXED BASE PEAK REPORT FOR SEQUENCING -----------------------------------

MixedReport <- function(table,consecutive_threshold, th_mixed_position) {
  
  # determine the mixed base peaks status using the defined threshold
  table <- dplyr::mutate(table, non_mixed=ifelse(primary_ratio>=th_mixed_position,T,F))
  
  serial_min <- c() 
  serial_max <- c()
  
  if (!is.na(consecutive_threshold)) {
    mixed <- table$non_mixed
    positions <- c(0:(consecutive_threshold-1)) # depending on the threshold, increment for positions to compare mixed status

    for (n in 1:(length(mixed)-consecutive_threshold)) {
      if (!FALSE %in% mixed[n+positions]) {
        serial_min <- append(serial_min, n) # first positions of 5 consecutives positions series of non-mixed peaks positions
      }
    }
    
    for (n in consecutive_threshold:length(mixed)) {
      if (!FALSE %in% mixed[n-positions]) {
        serial_max <- append(serial_max, n) # last positions of 5 consecutives positions series of non-mixed peaks positions
      }
    }
  }
  
  MixedReport <- list()
  
  MixedReport[["MixedPosCutOff"]] <- th_mixed_position
  MixedReport[["rawSeqLength"]] <- max(table$position)
  MixedReport[["rawPrimaryRatios"]] <- table$primary_ratio
  MixedReport[["rawNonMixedPos"]] <- table$non_mixed
  MixedReport[["rawMeanPrimaryRatio"]] <- round(mean(table$primary_ratio),2)
  MixedReport[["rawPercNonMixedPos"]] <- round(length(which(MixedReport$rawNonMixedPos==T))/length(MixedReport$rawNonMixedPos)*100,2)
  
  
  MixedReport[["trimmedStartPos"]] <- NA
  MixedReport[["trimmedFinishPos"]] <- NA
  MixedReport[["trimmedSeqLength"]] <- NA
  MixedReport[["trimmedPrimaryRatios"]] <- NA
  MixedReport[["trimmedMeanPrimaryRatio"]] <- NA
  MixedReport[["trimmedNonMixedPos"]] <- NA
  MixedReport[["trimmedPercNonMixedPos"]] <- NA
  
  if (is.null(serial_min)==F & is.null(serial_max)==F & !is.na(consecutive_threshold)) {
    if (length(serial_min)==1 & length(serial_max)==1) {
      MixedReport$trimmedStartPos <- serial_min
      MixedReport$trimmedFinishPos <- serial_max
    } else {
      MixedReport$trimmedStartPos <- min(serial_min)
      MixedReport$trimmedFinishPos <- max(serial_max)
    }
    MixedReport$trimmedSeqLength <- MixedReport$trimmedFinishPos - MixedReport$trimmedStartPos +1
    MixedReport$trimmedPrimaryRatios <- MixedReport$rawPrimaryRatios[MixedReport$trimmedStartPos:MixedReport$trimmedFinishPos]
    MixedReport$trimmedMeanPrimaryRatio <- round(mean(MixedReport$trimmedPrimaryRatios),2)
    MixedReport$trimmedNonMixedPos <- MixedReport$rawNonMixedPos[MixedReport$trimmedStartPos:MixedReport$trimmedFinishPos]
    MixedReport$trimmedPercNonMixedPos <- round(length(which(MixedReport$trimmedNonMixedPos==T))/length(MixedReport$trimmedNonMixedPos)*100,2)
  }
  
  return(MixedReport)
}


# TRIMMED SEQUENCE OF SEQUENCING  -----------------------------------------

trimSeq <- function(sangerRead,MixedReport) {
  Trimmed_seq <- list()
  Trimmed_seq[["StartPos"]] <- max(sangerRead@QualityReport@trimmedStartPos, MixedReport$trimmedStartPos)
  Trimmed_seq[["FinishPos"]] <- min(sangerRead@QualityReport@trimmedFinishPos, MixedReport$trimmedFinishPos)
  
  if(is.na(Trimmed_seq$StartPos) & is.na(Trimmed_seq$FinishPos)) {
    
    Trimmed_seq[["SeqLength"]] <- NA
    
    Trimmed_seq[["qualityPhredScores"]] <- NA
    Trimmed_seq[["MeanQualityScore"]] <- NA
    Trimmed_seq[["PrimaryRatios"]] <- NA
    Trimmed_seq[["MeanPrimaryRatios"]] <- NA
    Trimmed_seq[["NonMixedPos"]] <- NA
    Trimmed_seq[["PercNonMixedPos"]] <- NA
    Trimmed_seq[["PrimarySeq"]] <- NA
    
  } else {
    
    Trimmed_seq[["SeqLength"]] <- Trimmed_seq$FinishPos - Trimmed_seq$StartPos +1
    
    if(Trimmed_seq$SeqLength < 0) {
      Trimmed_seq$StartPos <- NA
      Trimmed_seq$FinishPos <- NA
      Trimmed_seq$SeqLength <- NA  
      
      Trimmed_seq[["qualityPhredScores"]] <- NA
      Trimmed_seq[["MeanQualityScore"]] <- NA
      Trimmed_seq[["PrimaryRatios"]] <- NA
      Trimmed_seq[["MeanPrimaryRatios"]] <- NA
      Trimmed_seq[["NonMixedPos"]] <- NA
      Trimmed_seq[["PercNonMixedPos"]] <- NA
      Trimmed_seq[["PrimarySeq"]] <- NA
      
    } else {
      
      Trimmed_seq[["qualityPhredScores"]] <- sangerRead@QualityReport@qualityPhredScores[Trimmed_seq$StartPos:Trimmed_seq$FinishPos]
      Trimmed_seq[["MeanQualityScore"]] <- round(mean(Trimmed_seq$qualityPhredScores),2)
      Trimmed_seq[["PrimaryRatios"]] <- MixedReport$rawPrimaryRatios[Trimmed_seq$StartPos:Trimmed_seq$FinishPos]
      Trimmed_seq[["MeanPrimaryRatios"]] <- round(mean(Trimmed_seq$PrimaryRatios),2)
      Trimmed_seq[["NonMixedPos"]] <- MixedReport$rawNonMixedPos[Trimmed_seq$StartPos:Trimmed_seq$FinishPos]
      Trimmed_seq[["PercNonMixedPos"]] <- round(length(which(Trimmed_seq$NonMixedPos==T))/length(Trimmed_seq$NonMixedPos)*100,2)
      Trimmed_seq[["PrimarySeq"]] <- primarySeq(sangerRead)[Trimmed_seq$StartPos:Trimmed_seq$FinishPos]
    }
  }
  
  
  
  return(Trimmed_seq)
}




# TRIMING PLOT OF SEQUENCING ----------------------------------------------

trimming_plot <- function(sangerRead, seqtable, MixedReport, Trimmed_seq, 
                          th_quality_phred, th_mixed_position, downloadname, colors_tables) {
  
  # Add colors depending on threshold
  seqtable <- dplyr::mutate(seqtable, 
                             quality_color = ifelse(quality_score>=th_quality_phred,colors_tables["correct"],colors_tables["incorrect"]),
                             ratio_color=  ifelse(primary_ratio>=th_mixed_position,colors_tables["correct"],colors_tables["incorrect"]))
  
  # Quality score plot
  quality_plot <- plot_ly(seqtable, x = ~position, y = ~quality_score, type = "scatter", mode="markers", marker = list(color= ~quality_color),
                          hoverinfo = 'text',
                          height=450,
                          text = ~paste("</br> Position:", position,
                                        "</br> Primary base:", primary_base,
                                        "</br> Quality score:", quality_score)) %>% 
    layout(xaxis = list(range = c(0,(MixedReport$rawSeqLength+2))),
           yaxis = list(title = "Quality phred score",range = c(0,65), ticksuffix="  "))
  
  if(!is.na(sangerRead@QualityReport@trimmedStartPos) & !is.na(sangerRead@QualityReport@trimmedFinishPos)) {
    quality_plot <- quality_plot %>% 
      add_segments(x = sangerRead@QualityReport@trimmedStartPos, xend = sangerRead@QualityReport@trimmedStartPos, y = 0, yend = 65, marker= NULL,
                   line=list(width=2,color= "coral"), text = ~paste("</br> Trimmed start position :", sangerRead@QualityReport@trimmedStartPos)) %>% 
      add_segments(x = sangerRead@QualityReport@trimmedFinishPos, xend = sangerRead@QualityReport@trimmedFinishPos, y = 0, yend = 65, marker= NULL,
                   line=list(width=2,color= "coral"),text = ~paste("</br> Trimmed end position :", sangerRead@QualityReport@trimmedFinishPos))            
  }
  
  # Mixed base peak plot
  mixed_plot <- plot_ly(seqtable, x = ~position, y = ~primary_ratio, type = "scatter", mode="markers", marker = list(color= ~ratio_color),
                        hoverinfo = 'text',
                        height=450,
                        text = ~paste("</br> Position:", position,
                                      "</br> Primary base:", primary_base,
                                      "</br> Primary ratio:", primary_ratio)) %>% 
    layout(xaxis = list(range = c(0,(MixedReport$rawSeqLength+2))),
           yaxis = list(title = "Peak ratio of primary base",range = c(0,1.1), ticksuffix="  "))
  
  if(!is.na(MixedReport$trimmedStartPos) & !is.na(MixedReport$trimmedFinishPos)) {
    mixed_plot <- mixed_plot %>% 
      add_segments(x = MixedReport$trimmedStartPos, xend = MixedReport$trimmedStartPos, y = 0, yend = 1.1, marker= NULL,
                   line=list(width=2,color= "turquoise"), text = ~paste("</br> Trimmed start position :", MixedReport$trimmedStartPos)) %>% 
      add_segments(x = MixedReport$trimmedFinishPos, xend = MixedReport$trimmedFinishPos, y = 0, yend = 1.1, marker= NULL,
                   line=list(width=2,color= "turquoise"),text = ~paste("</br> Trimmed end position :", MixedReport$trimmedFinishPos)) 
  }
  
  
  # Add for each position the sequence info
  seqtable <- dplyr::mutate(seqtable,
                             raw = "Raw sequence",
                            quality_trim =
                            ifelse(!is.na(sangerRead@QualityReport@trimmedStartPos) & !is.na(sangerRead@QualityReport@trimmedFinishPos), 
                                   ifelse(position %in% sangerRead@QualityReport@trimmedStartPos:sangerRead@QualityReport@trimmedFinishPos, "Quality score trimming", NA),
                                   NA),
                            mixed_trim =
                              ifelse(!is.na(MixedReport$trimmedStartPos) & !is.na(MixedReport$trimmedFinishPos), 
                                     ifelse(position %in% MixedReport$trimmedStartPos:MixedReport$trimmedFinishPos , "Mixed peak trimming", NA),
                                     NA),
                            trimmed = 
                              ifelse(!is.na(Trimmed_seq$StartPos) & !is.na(Trimmed_seq$FinishPos), 
                              ifelse(position %in% Trimmed_seq$StartPos:Trimmed_seq$FinishPos , "Trimmed sequence", NA), NA))
  
  # Plot sequences on position axis
  seq_plot <- plot_ly(seqtable, x = ~position, y = ~raw, type = "scatter", mode="lines", 
                      line=list(width=10, color = as.vector(colors_tables["incorrect"])),
                      height = 200,
                      hoverinfo = 'text',
                      text = ~paste("</br> <b> Raw sequence </b>",
                                    "</br> Length :",MixedReport$rawSeqLength,"bp",
                                    "</br> Mean of phred quality score :", round(sangerRead@QualityReport@rawMeanQualityScore,2),
                                    "</br> Mean of primary peak ratio :", MixedReport$rawMeanPrimaryRatio,
                                    "</br> Percentage of non-mixed positions :", MixedReport$rawPercNonMixedPos*100,"%")) %>% 
    add_trace(seqtable, x = ~position, y = ~mixed_trim, type = "scatter", mode="lines", 
              line=list(width=10, color = "turquoise"),
              hoverinfo = 'text',
              text = ~paste("</br> <b> Trimmed sequence based on mixed base peak ratios </b>",
                            "</br> Length :", MixedReport$trimmedSeqLength,
                            "</br> Mean of primary peak ratio :", MixedReport$trimmedMeanPrimaryRatio,
                            "</br> Percentage of non-mixed positions :", MixedReport$trimmedPercNonMixedPos,"%")) %>% 
    add_trace(seqtable,  x = ~position, y = ~quality_trim, type = "scatter", mode="lines",
              line=list(width=10, color = "coral"),
              hoverinfo = 'text',
              text = ~paste("</br> <b> Trimmed sequence based on phred quality scores </b>",
                            "</br> Length :", sangerRead@QualityReport@trimmedSeqLength,
                            "</br> Mean of phred quality score :", round(sangerRead@QualityReport@trimmedMeanQualityScore,2))) %>%
    add_trace(seqtable,  x = ~position, y = ~trimmed, type = "scatter", mode="lines",
              line=list(width=10, color = as.vector(colors_tables["correct"])),
              hoverinfo = 'text',
              text = ~paste("</br> <b> Trimmed sequence </b>",
                            "</br> Length :",Trimmed_seq$SeqLength,"bp",
                            "</br> Mean of phred quality score :", Trimmed_seq$MeanQualityScore,
                            "</br> Mean of primary peak ratio :", Trimmed_seq$MeanPrimaryRatio,
                            "</br> Percentage of non-mixed positions :", Trimmed_seq$PercNonMixedPos*100,"%")) %>%
    layout(xaxis = list(title = "Base position",range = c(0,(MixedReport$rawSeqLength+2))), 
           yaxis = list(title = "", ticksuffix="  "))
  
  # Combine plots in one plot
  plot <- subplot(quality_plot,mixed_plot,seq_plot, nrows = 3,  titleY = TRUE, titleX = TRUE, shareX = TRUE,
                  heights = c(0.40,0.40,0.2))%>%
    layout(legend = list(x = 1.05, y = 1), 
           title = "Quality score and primary base peak ratio plots",
           xaxis = list(title = "Base position"), margin=0.15,
           showlegend = FALSE, height=1000, width=800) %>% 
    config(displaylogo = FALSE, toImageButtonOptions = list(format = "png", # one of png, svg, jpeg, webp
                                                             filename = downloadname,
                                                             scale = 3))
  
  return(plot)
}


# BISULFITE CONVERSION RATES FUNCTION -------------------------------------

BisulfiteConvRate <- function (sequences, alg_pos, seqtable_trim_F, seqtable_trim_R, indel, useSeqF, useSeqR) {

  
  # Get sequence without C from CG dinucleotides
  
  ## chosen strand for amplification after bisulfite conversion corresponds to seq$seq
  sequence <- str_replace_all(as.character(sequences$seq),"CG","NG")
  sequence <- DNAString(sequence)
  
  # Find all remaining C positions in sequence
  C_data <- matchPattern("C", sequence, max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE, algorithm="auto")
  C_data <- data.frame(C_nb = seq(1:length(C_data)), position_dna = ranges(C_data)@start)
  
  # Get positions of these C in reverse sequence : G that should be transformed into A
  C_data_RC <- data.frame(C_nb = seq(1:nrow(C_data)), position_dna = length(sequence) - C_data$position_dna +1)
  C_data_RC <- C_data_RC[order(C_data_RC$position_dna),]
  
  # Find base peaks and conversion rate at each C position (so G in reverse sequence)
  
  BisConRate <- list()
  
  ## Forward
  if(useSeqF==T){
    seqtable_trim_F <- as.data.frame(seqtable_trim_F)
    C_data <- C_data %>% 
      # Filter positions comprised in alignment ranges
      filter(position_dna>=alg_pos$subject_F_start & position_dna <= (alg_pos$subject_F_end+length(indel$dna_ref_F))) %>% 
      # Ajust position in case of gaps on DNA reference sequence aligned
      dplyr::mutate(position_dna = position_dna + rowSums(outer(position_dna,indel$dna_ref_F,`>`)))
    
    C_data <- C_data %>% 
      # Find C equivalent position from DNA reference sequence on sequencing result
      dplyr::mutate(position_seq_trim = position_dna - alg_pos$subject_F_start + alg_pos$pattern_F_start) %>%
      # Ajust position in case of gaps on sequencing result aligned
      dplyr::mutate(position_seq_trim = position_seq_trim - rowSums(outer(position_seq_trim,indel$sequencing_F,`>`))) %>% 
      # Get the position on raw sequencing file
      dplyr::mutate(position_seq_raw = seqtable_trim_F[position_seq_trim,"position"]) %>%
      # Retrieve peak values from matrix and compute conversion rate
      dplyr::mutate(peak_A=seqtable_trim_F[position_seq_trim,"peak_A"],
                    peak_T=seqtable_trim_F[position_seq_trim,"peak_T"],
                    peak_G=seqtable_trim_F[position_seq_trim,"peak_G"],
                    peak_C=seqtable_trim_F[position_seq_trim,"peak_C"]) %>% 
      dplyr::mutate(conversion_rate=round((peak_T/(peak_C+peak_T)),2))
    BisConRate[["F"]] <- C_data
  }
  
  ## Reverse
  if(useSeqR==T){
    seqtable_trim_R <- as.data.frame(seqtable_trim_R)
    C_data_RC <- C_data_RC %>% 
      # Filter positions comprised in alignment ranges
      filter(position_dna>=alg_pos$subject_R_start & position_dna <= (alg_pos$subject_R_end+length(indel$dna_ref_R))) %>% 
      # Ajust position in case of gaps on DNA reference sequence aligned
      dplyr::mutate(position_dna = position_dna + rowSums(outer(position_dna,indel$dna_ref_R,`>`)))
    
    C_data_RC <- C_data_RC %>% 
      # Find C equivalent position from DNA reference sequence on sequencing result
      dplyr::mutate(position_seq_trim = position_dna - alg_pos$subject_R_start + alg_pos$pattern_R_start) %>%
      # Ajust position in case of gaps on sequencing result aligned
      dplyr::mutate(position_seq_trim = position_seq_trim - rowSums(outer(position_seq_trim,indel$sequencing_R,`>`))) %>%
      # Get the position on raw sequencing file
      dplyr::mutate(position_seq_raw = seqtable_trim_R[position_seq_trim,"position"]) %>%
      # Retrieve peak values from matrix and compute conversion rate
      dplyr::mutate(peak_A=seqtable_trim_R[position_seq_trim,"peak_A"],
                    peak_T=seqtable_trim_R[position_seq_trim,"peak_T"],
                    peak_G=seqtable_trim_R[position_seq_trim,"peak_G"],
                    peak_C=seqtable_trim_R[position_seq_trim,"peak_C"]) %>% 
      dplyr::mutate(conversion_rate=round((peak_A/(peak_G+peak_A)),2))
    BisConRate[["R"]] <- C_data_RC
  }
  
  return(BisConRate)
}



# GET MAXIMUM ALIGNED SEQUENCE FUNCTION -----------------------------------


getMaxAligned <- function (sequenceBS, seq_strand, alg_coord_start, alg_coord_end, seq_start, seq_end) {
  if (seq_strand=="plus") {
    pos_start = alg_coord_start-seq_start+1
    pos_end = alg_coord_end-seq_start+1
    seq_max_aligned = sequenceBS[pos_start:pos_end]
    return(seq_max_aligned)
  }
  if (seq_strand=="minus") {
    pos_start = seq_end-alg_coord_end+1
    pos_end = seq_end-alg_coord_start+1
    seq_max_aligned = sequenceBS[pos_start:pos_end]
    return(seq_max_aligned)
  }
}



# CHECK CG COVER FUNCTION -------------------------------------------------

checkCGcover <- function(data_pos,type="S",subject_start,subject_end) {
  if (type == "S") {
    data_pos <- data_pos[which(data_pos[,c(5)] >= subject_start & data_pos[,c(5)] <= subject_end),]
  }
  if (type == "RC") {
    data_pos <- data_pos[which(data_pos[,c(5)] >= subject_start & data_pos[,c(5)] <= subject_end),]
  }
  if(nrow(data_pos)==0) {
    NoCGcovered <- TRUE
  } else {
    NoCGcovered <- FALSE
  }
  return(NoCGcovered)
}




# GET BASE PEAKS FUNCTION -------------------------------------------------

getPeaks <- function (x, alg_pos, seqtable_trim, indel, sequencing) {
  # Get peaks of bases function
  seqtable_trim <- as.data.frame(seqtable_trim)
  
  ## FORWARD
  if (sequencing=="forward") {
    CG_data <- x %>% 
      # Filter CG positions comprised in alignment ranges
      filter (C_pos >=alg_pos$subject_F_start & G_pos <= (alg_pos$subject_F_end+length(indel$dna_ref_F))) %>% 
      # Ajust position in case of gaps on DNA reference sequence aligned
      dplyr::mutate(C_pos = C_pos + rowSums(outer(C_pos,indel$dna_ref_F,`>`))) %>%
      dplyr::mutate(G_pos = G_pos + rowSums(outer(G_pos,indel$dna_ref_F,`>`)))
    
    # Retrieve peak values from matrix and compute methylation rate
    CG_data <- CG_data %>% 
      # Get position of C from CG on sequencing result
      dplyr::mutate(position_seq_trim = C_pos - alg_pos$subject_F_start + alg_pos$pattern_F_start) %>% 
      # Ajust position in case of gaps on sequencing result aligned
      dplyr::mutate(position_seq_trim = position_seq_trim - rowSums(outer(position_seq_trim,indel$sequencing_F,`>`))) %>% 
      # Get the position on raw sequencing file
      dplyr::mutate(position_seq_raw = seqtable_trim[position_seq_trim,"position"]) %>%
      # Get base peaks values
      dplyr::mutate(peak_A=seqtable_trim[position_seq_trim,"peak_A"],
                    peak_T=seqtable_trim[position_seq_trim,"peak_T"],
                    peak_G=seqtable_trim[position_seq_trim,"peak_G"],
                    peak_C=seqtable_trim[position_seq_trim,"peak_C"]) %>% 
      # Methylation percentage computation
      dplyr::mutate(meth_F=round((peak_C/(peak_C+peak_T))*100,2))
  }
  
  ## REVERSE
  if (sequencing=="reverse") {
    CG_data <- x %>%
      # Filter CG positions comprised in alignment ranges
      filter (C_pos_RC >=alg_pos$subject_R_start & G_pos_RC <= (alg_pos$subject_R_end+length(indel$dna_ref_R))) %>%
      # Ajust position in case of gaps on DNA reference sequence aligned
      dplyr::mutate(C_pos_RC = C_pos_RC + rowSums(outer(C_pos_RC,indel$dna_ref_R,`>`))) %>%
      dplyr::mutate(G_pos_RC = G_pos_RC + rowSums(outer(G_pos_RC,indel$dna_ref_R,`>`)))
    
    # Retrieve peak values from matrix and compute methylation rate
    CG_data <- CG_data %>% 
      # Get position of C from CG (so G for reverse) on sequencing result
      dplyr::mutate(position_seq_trim = G_pos_RC - alg_pos$subject_R_start + alg_pos$pattern_R_start) %>%
      # Ajust position in case of gaps on sequencing result aligned
      dplyr::mutate(position_seq_trim = position_seq_trim - rowSums(outer(position_seq_trim,indel$sequencing_R,`>`))) %>% 
      # Get the position on raw sequencing file
      dplyr::mutate(position_seq_raw = seqtable_trim[position_seq_trim,"position"]) %>%
      # Get base peaks values
      dplyr::mutate(peak_A=seqtable_trim[position_seq_trim,"peak_A"],
                    peak_T=seqtable_trim[position_seq_trim,"peak_T"],
                    peak_G=seqtable_trim[position_seq_trim,"peak_G"],
                    peak_C=seqtable_trim[position_seq_trim,"peak_C"]) %>% 
      # Methylation percentage computation
      dplyr::mutate(meth_R=round((peak_G/(peak_G+peak_A))*100,2))
  }
  return(CG_data)
}




# INDIVIDUAL METHYLATION PLOT FUNCTION ------------------------------------

individual_meth_plot <- function(meth, start, end, genome, plot_title, bases_colors) {
  # Data track of methylation
  genomeID <- unique(genome@seqinfo@genome)
  track_meth <- DataTrack(range=meth, genome=genomeID, lwd.title = 2, showColorBar = FALSE, cex.axis = 1.5, lineheight = 1,
                          type="heatmap", background.title = "transparent", col = "black", col.border.title="transparent",
                          col.axis="transparent", col.frame="black", ylim = c(0,100), 
                          lwd.border=10, lwd=10, min.height=10, min.width=10, lwd.grid=10, 
                          fontcolor.title="black", name="Methylation",
                          ncolor=5, gradient=c("grey90","black"))
  
  # CG track
  track_CG <- AnnotationTrack(meth, name ="CpG", fontcolor.title="black", col=NULL, fill = "navy", 
                              background.title = "transparent", col.border.title = "transparent", 
                              rotate.title = F)
  # Sequence track
  sTrack <- SequenceTrack(genome, add53=F , noLetters=T, fontcolor=bases_colors)
  # Genomic axis track
  gtrack <- GenomeAxisTrack(cex=2,col="black",fontcolor="black", distFromAxis = 5, lwd=2)
  plot_ticks <- seq(from=round(start,-2), to=round(end+50,-2), by=20)
  ############################################################################################
  plotTracks(list(track_meth,track_CG,sTrack,gtrack),
             sizes = c(6,1,1,2.5), main=plot_title, cex.main = 2.2,
             cex.title = 2, title.width = 1.2,from=start, to=end, add53=TRUE,labelPos="below",ticksAt=plot_ticks)
}



# COLOR BAR OR TILE WITH CONDITION FORMATTABLE FUNCTION -------------------

color_condi <- function(operator, value, colorT, colorF, textcolorT = NULL, textcolorF = NULL, block = "inline-block", fun = NULL) { 
  if(!is.null(fun)) {fun <- match.fun(fun)}
  
  if(!is.null(textcolorT)&!is.null(textcolorF)) {
    formatter("span", 
              style = x ~ formattable::style(display = block, 
                                             padding = "0 4px", 
                                             "border-radius" = "4px", 
                                             "padding-right" = "2px",
                                             direction = "rtl", 
                                             "unicode-bidi" = "plaintext",
                                             "color" = csscolor(
                                               ifelse(is.na(x) | x=="  NA","white",
                                                      if (operator==">=") {ifelse(x >= value, textcolorT, textcolorF)}
                                                      else if (operator=="<=") {ifelse(x <= value, textcolorT, textcolorF)}
                                                      else if (operator=="<") {ifelse(x < value, textcolorT, textcolorF)}
                                                      else if (operator==">") {ifelse(x > value, textcolorT, textcolorF)}
                                                      else if (operator=="==") {ifelse(x == value, textcolorT, textcolorF)}
                                                      else if (operator=="!=") {ifelse(x != value, textcolorT, textcolorF)})),
                                             "background-color" = csscolor(
                                               ifelse(is.na(x) | x=="  NA","white",
                                                      if (operator==">=") {ifelse(x >= value, colorT, colorF)}
                                                      else if (operator=="<=") {ifelse(x <= value, colorT, colorF)}
                                                      else if (operator=="<") {ifelse(x < value, colorT, colorF)}
                                                      else if (operator==">") {ifelse(x > value, colorT, colorF)}
                                                      else if (operator=="==") {ifelse(x == value, colorT, colorF)}
                                                      else if (operator=="!=") {ifelse(x != value, colorT, colorF)})),
                                             width = if(!is.null(fun)) { percent(fun(as.numeric(x)))}))
  }
  
  else {
    formatter("span", 
              style = x ~ formattable::style(display = block, 
                                             padding = "0 4px", 
                                             "border-radius" = "4px", 
                                             "padding-right" = "2px",
                                             direction = "rtl", 
                                             "unicode-bidi" = "plaintext",
                                             "background-color" = csscolor(
                                               ifelse(is.na(x) | x=="  NA","white",
                                                      if (operator==">=") {ifelse(x >= value, colorT, colorF)}
                                                      else if (operator=="<=") {ifelse(x <= value, colorT, colorF)}
                                                      else if (operator=="<") {ifelse(x < value, colorT, colorF)}
                                                      else if (operator==">") {ifelse(x > value, colorT, colorF)}
                                                      else if (operator=="==") {ifelse(x == value, colorT, colorF)}
                                                      else if (operator=="!=") {ifelse(x != value, colorT, colorF)})),
                                             width = if(!is.null(fun)) { percent(fun(as.numeric(x)))}))
  }
}




# COLOR TEXT AND BACKGROUND FORMATTABLE FUNCTION --------------------------

color_format <- function(text_color, background_color) { 
  formatter("span", style = x ~ formattable::style(display = "block",
                                                   padding = "0 4px", 
                                                   "border-radius" = "4px",
                                                   "color" = text_color,
                                                   "background-color" = background_color)) }






# COLOR OF BASE BACKGROUND FORMATTABLE FUNCTION ---------------------------

color_base <- function(text_color, bases_colors) { 
  formatter("span", style = x ~ formattable::style(display = "block",
                                                   padding = "0 4px", 
                                                   "border-radius" = "4px",
                                                   "color" = text_color,
                                                   "background-color" = ifelse(x=="A", bases_colors["A"],
                                                                               ifelse(x=="T", bases_colors["T"],
                                                                                      ifelse(x=="G", bases_colors["G"],
                                                                                             ifelse(x=="C", bases_colors["C"],
                                                                                                    ifelse(x=="N", bases_colors["N"],""))))))) }



# COLORS FOR METHYLATION INTERVALS FORMATTABLE FUNCTION -------------------

color_meth <- function(b_colors=c(grey(0.80),grey(0.5),grey(0.3),"black"),t_colors=c("black","white")) { 
  formatter("span", 
            style = x ~ formattable::style(display = "block", 
                                           padding = "0 4px", 
                                           "border-radius" = "4px", 
                                           "padding-right" = "2px",
                                           direction = "rtl", 
                                           "unicode-bidi" = "plaintext",
                                           "background-color" = csscolor(ifelse(is.na(x),"#F2F2F2",
                                                                                ifelse(x<=25,b_colors[1],
                                                                                       ifelse(x<=50,b_colors[2],
                                                                                              ifelse(x<=75,b_colors[3],
                                                                                                     ifelse(x<=100,b_colors[4],"white")))))),
                                           "color" = csscolor(ifelse(is.na(x),"black",
                                                                     ifelse(x<=50,t_colors[1],t_colors[2])))))
}



# EXPORT FORMATTABLE FUNCTION ---------------------------------------------

export_formattable <- function(f, file, width = "100%", height = NULL, background = "white", delay = 0.2){
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
}





# ADD ALPHA TO COLORS FUNCTION --------------------------------------------

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}






# DISPLAY COLORS FROM A COLOR VECTOR --------------------------------------


displayColors <- function(a) {
  image(1:length(a), 1, as.matrix(1:length(a)), col=a, axes=FALSE , xlab="", ylab="")
  axis(1, at=1:length(a), labels = names(a), srt=45, las=2, lwd.ticks=0, lwd = 0)}





# DISPLAY COLORS AND SHAPES FOR PLOTS -------------------------------------


displayColorsShapes <- function(a,b) {
  
  points <- data.frame(
    group=rep(as.factor(1:min(length(a),length(b))),3), 
    value=c(rep(1,min(length(a),length(b))),rep(2,min(length(a),length(b))),rep(3,min(length(a),length(b))))
  )
  a <- as.vector(a[1:min(length(a),length(b))])
  b <- as.vector(b[1:min(length(a),length(b))])
  
  ggplot(points,aes(x=value, y=group)) +
    geom_line(aes(col=group)) +
    geom_point(aes(col=group, shape=group), fill="white", size=2) +
    scale_y_discrete(limits=rev) +
    scale_fill_manual(values=a) +
    scale_color_manual(values=a) +
    scale_shape_manual(values=b) +
    theme_bw() + 
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.text.x.bottom = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.text = element_text(size=11),
          plot.margin = unit(c(1,1,1,1),"line"))
}









#************************************************************************************









# PREPROCESS DATA FUNCTION ------------------------------------------------

preprocess_data <- function(data_files, GroupOrder, cloning_exp, clone_thresholds) {
  
  
  # create list of dataframes from all file paths :
  data_all <- lapply(data_files, read.table, header=TRUE, check.names = FALSE)
  
  
  # For clones : Convert methylation levels to methylated or unmethylated 
  if (cloning_exp==T) {
    data_all <- lapply(data_all, function(x) {
      
      x <- mutate(x,meth=meth_avg)
      # if more than 20% of CpG have a methylation level between 20% and 80%, it can be a mix of clones and results are not accurates
      # so all methylation are considered as Not Available
      if( length( which(x$meth>clone_thresholds[1] & x$meth<clone_thresholds[2]) ) > (length(x$meth)*clone_thresholds[3]) ) { # if more than 20% of partial positions
        x$meth <- NA
      } 
      x[which(x$meth>clone_thresholds[1] & x$meth<clone_thresholds[2]),"meth"] <- NA
      x[which(x$meth<=clone_thresholds[1]),"meth"] <- 0
      x[which(x$meth>=clone_thresholds[2]),"meth"] <- 100
      return(x)
    })
  } else if (cloning_exp==F) {
    data_all <- lapply(data_all, function(x) {
      x <- mutate(x,meth=meth_avg)
      return(x)
    })
  }
    
  # Retrieve the list of all CpG positions
  data_cpg <- lapply(data_all, function(x) {
    x <- x[,c(1:4)] })
  data_cpg <- suppressMessages(purrr::reduce(data_cpg, full_join))
  data_cpg <- data_cpg[order(data_cpg$CG_nb),]
  
  # Add rows for missing positions with NA and add sample information
  data_all <- lapply(data_all, function(x) {
    
    # Add rows for missing CpG positions
    x <- suppressMessages(full_join(data_cpg,x))
    
    # Complete sample info in the new rows
    # check if there is already a value or NA for collection, rep and clone information
    if(length(unique(x$collection))!=1) {x$collection <- rep(unique(x$collection)[which(!is.na(unique(x$collection)))],length(x$collection))}
    x$group <- rep(unique(x$group)[which(!is.na(unique(x$group)))],length(x$group))
    if(length(unique(x$rep))!=1) {x$rep <- rep(unique(x$rep)[which(!is.na(unique(x$rep)))],length(x$rep))}
    if(length(unique(x$rep))!=1) {x$clone <- rep(unique(x$clone)[which(!is.na(unique(x$clone)))],length(x$clone))}
    x$alg_coord_start <- rep(unique(x$alg_coord_start)[which(!is.na(unique(x$alg_coord_start)))],length(x$alg_coord_start))
    x$alg_coord_end<- rep(unique(x$alg_coord_end)[which(!is.na(unique(x$alg_coord_end)))],length(x$alg_coord_end))
    
    return(x)
  })
  
  ####### COORDINATES MIN & MAX
  # Retrieve the minimal and maximal coordinates covered in all experiments
  min_coord_start <- min(unlist(lapply(data_all, function(x) unique(x[,"alg_coord_start"]))))
  max_coord_end <- max(unlist(lapply(data_all, function(x) unique(x[,"alg_coord_end"]))))
  coord <- list(start=min_coord_start,end=max_coord_end)
 
  
  ####### MAIN DATA TABLE
  all <- purrr::reduce(data_all, rbind)
  
  all$collection[is.na(all$collection)] <- ""
  all$rep[is.na(all$rep)] <- ""
  all$clone[is.na(all$clone)] <- ""
  
  # Add columns
  # (clone_ID = ID of clone for cloning)
  # (rep_ID = ID of replicate for direct)
  # (group_ID = ID unique to each condition, corresponds to mean of clones or replicates per group)
  all <- all %>% 
    dplyr::mutate(group_ID = gsub("^ ","",paste(collection,group))) %>% 
    dplyr::mutate(rep_ID = gsub("^ ","",gsub(" $","",paste(collection,group,rep)))) %>% 
    dplyr::mutate(clone_ID = gsub("^ ","",gsub(" $","",gsub("  "," ",paste(collection,group,rep,"clone",clone)))))
    
  # Change CG_nb to 1:n
  all$CG_nb <- as.factor(all$CG_nb)
  levels(all$CG_nb) <- c(1:length(unique(all$CG_nb)))
  
  
  # Duplicate start column as new column named position
  all[,"position"] <- all$start
  all[,"CG_coord"] <- paste0(all$chromosome,":",all$start,"-",all$end)
  all <- all[,c("CG_nb", "CG_coord", "chromosome", "start", "end", "position", 
                "meth_F", "meth_R", "meth_avg", "meth_sd", "meth",
                "collection", "group", "rep", "clone",
                "group_ID", "rep_ID", "clone_ID",
                "alg_coord_start", "alg_coord_end")]
  
  
  
  
  # Set group as ordered variable using GroupOrder levels
  if (length(setdiff(unique(all$group), GroupOrder))!=0 ) { # if a group is missing in GroupOrder, it goes at the end but it is not missing
    levels_group <- c(GroupOrder,setdiff(unique(all$group), GroupOrder))
    all$group <- ordered(all$group, levels=levels_group)
  } else {
    all$group <- ordered(all$group, levels=GroupOrder)
  }
  all <- all[order(all$group),]
  
  # Set collection as ordered by alphabetic order 
  all$collection <- ordered(all$collection, levels=sort(unique(all$collection)))
  all <- all[order(all$collection),]
  
  data = list()
  data[["all"]] <- all
  
  
  ####### TABLE
  # remove unnecessary information
  table <- all[,c("CG_nb","CG_coord","chromosome","start","end","position","meth","collection","group","rep","clone","group_ID","rep_ID","clone_ID")]
  
  
  # For cloning :
  if(cloning_exp==T) {
    
    table <- select(table,-c("rep","rep_ID"))
    
    ####### TABLE METHYLATION DATA OF CLONES
    tableClones <- table
    tableClones[,"samples"] <- tableClones$clone_ID
    table_clones <- list()
    for(i in unique(tableClones$group_ID)) {
      table_clones[[i]] <- tableClones[which(tableClones$group_ID==i),]
    }
    data[["table_clones"]] <- table_clones
    
    ####### TABLE METHYLATION MEANS PER GROUP
    table_group <- table %>% 
      group_by(position,collection,group) %>%
      dplyr::mutate(meth_sd = round(sd(meth,na.rm=T),2)) %>%
      dplyr::mutate(n_clones = sum(!is.na(meth))) %>% 
      dplyr::mutate(meth = round(mean(meth,na.rm=T),2)) %>% 
      select(-c("clone_ID","clone")) %>% 
      unique()
    table_group[,"samples"] <- table_group$group_ID
    data[["table_group"]] <- table_group
  }
  
  # For direct : 
  if(cloning_exp==F){
    
    table <- select(table,-c("clone","clone_ID"))
    
    ####### TABLE METHYLATION DATA OF EACH REPLICATES
    table_rep <- table 
    table_rep[,"samples"] <- table_rep$rep_ID
    data[["table_rep"]] <- table_rep
    
    ####### TABLE METHYLATION MEANS PER GROUP
    table_group <- table %>% 
      group_by(position,collection,group) %>%
      dplyr::mutate(meth_sd = round(sd(meth,na.rm=T),2)) %>%
      dplyr::mutate(n_rep = sum(!is.na(meth))) %>% 
      dplyr::mutate(meth = round(mean(meth,na.rm=T),2)) %>% 
      select(-c("rep_ID","rep")) %>% 
      unique()
    table_group[,"samples"] <- table_group$group_ID
    data[["table_group"]] <- table_group
    
  }
  
  data[["coord"]] <- coord
  return(data)
}





# TRANSFORM COMPLETE DATA TABLE TO METHYLATION DATA TABLE -----------------

methdata <- function(table){
  list_df <- list()
  for (i in unique(as.ordered(table$samples))) {
    meth_df_samples <- table[which(table$samples==i),c("chromosome","start","end","meth")]
    colnames(meth_df_samples)[4] <- i
    list_df <- list.append(list_df,meth_df_samples)
  }
  meth_df <- suppressMessages(purrr::reduce(list_df, full_join)) 
  meth_df <- as.data.frame(meth_df)
  return(meth_df)
}




# DENDROGRAM CLUSTERING FUNCTION ------------------------------------------

dendro_plot <- function(table, filename, plot_title) {
  
  # First, remove samples with only NA values (removed clones for example)
  table_complete <- table[complete.cases(table),]
  
  # Transform methylation data from table (df) to matrix to compute distance
  list_df <- list()
  for (i in unique(table_complete$samples)) {
    meth_df_samples <- table_complete[which(table_complete$samples==i),c("chromosome","start","end","meth")]
    colnames(meth_df_samples)[4] <- i
    list_df <- list.append(list_df,meth_df_samples)
  }
  meth_df_complete <- suppressMessages(purrr::reduce(list_df, full_join))
  meth_mat_complete <- as.matrix(t(meth_df_complete[,-c(1,2,3)]))
  
  # Get distances as matrix
  dist <- as.matrix(dist(meth_mat_complete))
  # Remove NA in matrix : cluster cannot be computed with NA in distance matrix
  dist <- as.dist(dist[rowSums(is.na(dist)) == 0, colSums(is.na(dist)) == 0, drop = FALSE]) # remove NA
  # Get hierarchical clustering
  hc <- hclust(dist)
  # Get dendrogram data
  hcdata <- dendro_data(hc)
  
  # Dendrogramm plot
  dendroplot <- ggdendrogram(hcdata, rotate = TRUE) + 
    labs(title = plot_title) +
    theme(
      line=element_line(size=1),
      plot.title.position = "plot",
      plot.title = element_text(hjust=0.5, size = 14, face = "bold", color = "black"),
      axis.text.x.bottom = element_text(size = 13, angle = 90, color = "black"),
      axis.text.y.left = element_text(size = 13, color = "black"),
      plot.margin = unit(c(1,1,1,1),"line"),
      aspect.ratio = length(unique(table$samples))/4)
  suppressMessages(
    ggsave(plot = dendroplot, filename = filename, 
           height = length(unique(table$samples))*70 + 400, 
           width = max( max(nchar(unique(table$samples)))*28 + 380 + 240 , nchar(plot_title)*30), dpi=300, units = c("px"))
  )
}





# LOLLIPOP METHYLATION PLOT -----------------------------------------------

lollipop_plot <- function(table, filename, coord, plot_title, SampleOrder, GroupOrder, plotType, MethLevels=c(2,4,5), pos_labels) {
  
  
  #************************************************************************************************
  # METHYLATIONS LEVELS
  if (MethLevels==2) { # for clones
    table <- table %>% 
      dplyr::mutate(meth_level = case_when(
        meth == 0 ~ "Unmethylated",
        meth == 100 ~ "Methylated",
        TRUE ~ as.character(meth)))
    table$meth_level <- factor(table$meth_level, levels=c("Unmethylated","Methylated"))
    grey_colors <- grey.colors(2, start=1, end=0)
    legend_name <- "Methylation status"
    legend_width <- 1200 #px
  }
  if (MethLevels==4) { # for direct replicates
    table <- table %>% 
      dplyr::mutate(meth_level = case_when(
        meth >= 0 & meth < 25 ~ "0-25",
        meth >= 25 & meth < 50 ~ "25-50",
        meth >= 50 & meth < 75 ~ "50-75",
        meth >= 75 & meth <= 100 ~ "75-100",
        TRUE ~ as.character(meth)))
    table$meth_level <- ordered(table$meth_level, levels=c("0-25","25-50","50-75","75-100"))
    grey_colors <- grey.colors(4, start=1, end=0)
    legend_name <- "Methylation levels"
    legend_width <- 1000 #px
  }
  if (MethLevels==5) { # for clones replicates 
    table <- table %>% 
      dplyr::mutate(meth_level = case_when(
        meth >= 0 & meth < 20 ~ "0-20",
        meth >= 20 & meth < 40 ~ "20-40",
        meth >= 40 & meth < 60 ~ "40-60",
        meth >= 60 & meth < 80 ~ "60-80",
        meth >= 80 & meth <= 100 ~ "80-100",
        TRUE ~ as.character(meth)))
    table$meth_level <- ordered(table$meth_level, levels=c("0-20","20-40","40-60","60-80","80-100"))
    grey_colors <- grey.colors(5, start=1, end=0)
    legend_name <- "Methylation levels"
    legend_width <- 1000 #px
  }
  
  #************************************************************************************************
  # Change position label depending on input request
  if(pos_labels=="coordinates") {
    pos_lab <-unique(table$CG_coord)
    angle <- 90
    vjust = 0.5
    height_pos <- max(nchar(table$CG_coord))*70
  } 
  if(pos_labels=="numbers") {
    pos_lab <- seq(from=1, to=length(unique(table$position)), by=1) 
    angle <- 0
    vjust = 0
    height_pos <- 150
  }
  if(pos_labels=="none") {
    pos_lab <- rep("",length(unique(table$position)))
    angle <- 0
    vjust = 0
    height_pos <- 150
  }
  
  
  #************************************************************************************************
  # Sort samples
  
  if(SampleOrder == "as-is"){
    # Reverse order of samples on y axis, to read from top to bottom
    table$samples <- ordered(table$samples, levels=rev(unique(table$samples)))
  }
  if(SampleOrder == "by-group"){
    # Order data frame by GroupOrder
    table <- table[order(ordered(table$group, levels = GroupOrder)),]
    # Reverse order of samples on y axis, to read from top to bottom
    table$samples <- ordered(table$samples, levels=rev(unique(table$samples)))
  }
  
  if(SampleOrder == "by-meth"){
    # Order data frame by methylation mean on sequence (all positions)
    table <- table %>% 
      group_by(samples) %>%
      dplyr::mutate(meth_mean = round(mean(meth,na.rm=T),2))
    table <- table[order(table$meth_mean),]
    # Reverse order of samples on y axis, to read from top to bottom
    table$samples <- ordered(table$samples, levels=rev(unique(table$samples)))
  }
  
  if(SampleOrder == "clust"){
    # First, remove samples with only NA values (removed clones for example)
    table <- table %>% 
      group_by(samples) %>%
      dplyr::mutate(meth_mean = round(mean(meth,na.rm=T),2))
    table_complete <- table[complete.cases(table),]
    
    # Transform methylation data from table (df) to matrix to compute distance
    list_df <- list()
    for (i in unique(table_complete$samples)) {
      meth_df_samples <- table_complete[which(table_complete$samples==i),c("chromosome","start","end","meth")]
      colnames(meth_df_samples)[4] <- i
      list_df <- list.append(list_df,meth_df_samples)
    }
    meth_df_complete <- suppressMessages(purrr::reduce(list_df, full_join)) 
    meth_mat_complete <- as.matrix(t(meth_df_complete[,-c(1,2,3)]))
      
    # Get distances as matrix
    dist <- as.matrix(dist(meth_mat_complete))
    # Remove NA in matrix : cluster cannot be computed with NA in distance matrix
    dist <- as.dist(dist[rowSums(is.na(dist)) == 0, colSums(is.na(dist)) == 0, drop = FALSE]) # remove NA
    # Get hierarchical clustering
    hc <- hclust(dist)
    # Get dendrogram data
    hcdata <- dendro_data(hc)
    
    # Order complete meth data by dendro order
    table$samples <- ordered(table$samples, levels=c(setdiff(table$samples, table_complete$samples),hcdata$labels$label))
  }
  
  
  #************************************************************************************************
  # Proportional
  if(plotType=="proportional"){
    
    # LENGHT AND RATIO
    width_x <- (coord$end - coord$start)*0.1
    height_y <- length(unique(table$samples))
    

    # PLOT
    ggplot(table)+
      geom_segment(aes(x=coord$start, xend=coord$end, y=samples, yend=samples), color="black") +
      geom_point(data=subset(table,!is.na(meth_level)), aes(x=position, y=samples, fill=meth_level), size=6, shape=21, col="black") +
      geom_point(data=subset(table,is.na(meth_level)), aes(x=position, y=samples), color="black", size=4, shape=4) + 
      scale_fill_manual(breaks = levels(table$meth_level), values = grey_colors, limits=levels(table$meth_level), name = legend_name) +
      scale_x_continuous(breaks=unique(table$position), labels = pos_lab, limits=c(coord$start,coord$end)) +
      labs(title = plot_title) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = "black"),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            axis.text.x = element_text(angle = angle,size = 13, color = "black", vjust = vjust),
            axis.text.y = element_text(size = 13, color = "black"),
            legend.title = element_text(size=11),
            legend.key.size = unit(2, 'mm'),
            legend.text = element_text(size=10),
            plot.margin = unit(c(1,1,1,1),"line"),
            aspect.ratio = height_y/width_x)+
      guides(fill = guide_legend(override.aes = list(size=3)))
    
    # Plot width and height on dependencies
    width_plot <- max( max(nchar(as.character(unique(table$samples))))*40 + (coord$end - coord$start)*18 + legend_width + 240 , nchar(plot_title)*60)
    height_plot <- height_pos + length(unique(table$samples))*160 + 175 + 240 
    
    # SAVE
    ggsave(filename = filename, width = width_plot, height = height_plot, dpi=600, units = c("px"))
  }
  
  #************************************************************************************************
  # Condensed
  if(plotType=="condensed"){
    
    table$CG_nb <- as.integer(table$CG_nb)
    
    # LENGHT AND RATIO
    width_x <- length(unique(table$position))
    height_y <- length(unique(table$samples))
    
    # PLOT
    ggplot(table)+
      geom_segment(aes(x=min(CG_nb), xend=max(CG_nb), y=samples, yend=samples), color="transparent") +
      geom_point(data=subset(table,!is.na(meth_level)), aes(x=CG_nb, y=samples, fill=meth_level), size=7, shape=21, col="black") +
      geom_point(data=subset(table,is.na(meth_level)), aes(x=CG_nb, y=samples), color="black", size=4, shape=4) +
      scale_fill_manual(breaks = levels(table$meth_level), values = grey_colors, limits = levels(table$meth_level), name = legend_name) +
      scale_x_discrete(breaks = unique(as.factor(table$CG_nb)), labels = pos_lab, limits = unique(as.factor(table$CG_nb))) +
      labs(title = plot_title) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = "black"),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            axis.text.x = element_text(angle = angle, size = 13, color = "black", vjust = vjust),
            axis.text.y = element_text(size = 13, color = "black"),
            legend.title = element_text(size=11),
            legend.text = element_text(size=10),
            legend.key.size = unit(2, 'mm'),
            plot.margin = unit(c(1,1,1,1),"line"),
            aspect.ratio = height_y/width_x)+
      guides(fill = guide_legend(override.aes = list(size=3)))
    
    # Plot width and height on dependencies
    width_plot <- max( max(nchar(as.character(unique(table$samples))))*40 + length(unique(table$position))*240 + legend_width + 240 , nchar(plot_title)*70)
    height_plot <- height_pos + length(unique(table$samples))*160 + 175 + 240 
    
    # SAVE
    ggsave(filename = filename, width = width_plot, height = height_plot, dpi=600, units = c("px"))
  }
  
}







# GENOMIC METHYLATION PLOT FUNCTION ---------------------------------------

genomic_plot <- function(table, filename, coord, plot_title, SampleOrder, GroupOrder, MethLevels=c(2,4,5), genome, bases_colors) {
  
  
  #************************************************************************************************
  # Sort samples
  
  if(SampleOrder == "as-is"){
    # Reverse order of samples on y axis, to read from top to bottom
    table$samples <- ordered(table$samples, levels=rev(unique(table$samples)))
  }
  if(SampleOrder == "by-group"){
    # Order data frame by GroupOrder
    table <- table[order(ordered(table$group, levels = GroupOrder)),]
    # Reverse order of samples on y axis, to read from top to bottom
    table$samples <- ordered(table$samples, levels=rev(unique(table$samples)))
  }
  
  if(SampleOrder == "by-meth"){
    # Order data frame by methylation mean on sequence (all positions)
    table <- table %>% 
      group_by(samples) %>%
      dplyr::mutate(meth_mean = round(mean(meth,na.rm=T),2))
    table <- table[order(table$meth_mean),]
    # Reverse order of samples on y axis, to read from top to bottom
    table$samples <- ordered(table$samples, levels=rev(unique(table$samples)))
  }
  
  if(SampleOrder == "clust"){
    # First, remove samples with only NA values (removed clones for example)
    table <- table %>% 
      group_by(samples) %>%
      dplyr::mutate(meth_mean = round(mean(meth,na.rm=T),2))
    table_complete <- table[complete.cases(table),]
    
    # Transform methylation data from table (df) to matrix to compute distance
    list_df <- list()
    for (i in unique(table_complete$samples)) {
      meth_df_samples <- table_complete[which(table_complete$samples==i),c("chromosome","start","end","meth")]
      colnames(meth_df_samples)[4] <- i
      list_df <- list.append(list_df,meth_df_samples)
    }
    meth_df_complete <- suppressMessages(purrr::reduce(list_df, full_join)) 
    meth_mat_complete <- as.matrix(t(meth_df_complete[,-c(1,2,3)]))
    
    # Get distances as matrix
    dist <- as.matrix(dist(meth_mat_complete))
    # Remove NA in matrix : cluster cannot be computed with NA in distance matrix
    dist <- as.dist(dist[rowSums(is.na(dist)) == 0, colSums(is.na(dist)) == 0, drop = FALSE]) # remove NA
    # Get hierarchical clustering
    hc <- hclust(dist)
    # Get dendrogram data
    hcdata <- dendro_data(hc)
    
    # Order complete meth data by dendro order
    table$samples <- ordered(table$samples, levels=c(setdiff(table$samples, table_complete$samples),hcdata$labels$label))
  }
  
  
  #************************************************************************************************
  # Get Methylation as data frame with samples as ordered
  list_df <- list()
  for (i in unique(as.character(table$samples))) {
    meth_df_samples <- table[which(table$samples==i),c("chromosome","start","end","meth")]
    colnames(meth_df_samples)[4] <- i
    list_df <- list.append(list_df,meth_df_samples)
  }
  meth_df <- suppressMessages(purrr::reduce(list_df, full_join)) 
  
  
  #************************************************************************************************
  # Data track
  # change column names to conserve spaces and - in data track title
  colnames(meth_df) <- gsub(" ","__",colnames(meth_df))
  colnames(meth_df) <- gsub("-","_",colnames(meth_df))
  genomeID <- metadata(genome)$genome
  track_meth <- DataTrack(range = meth_df, genome = genomeID, chromosome = unique(meth_df$chromosome), 
                          lwd.title=5, showColorBar=F, cex.axis=3, lineheight=0,
                          type="heatmap", showSampleNames=TRUE, cex.sampleNames=3, col.sampleNames="black", 
                          col.border.title="transparent", col="black", col.axis="transparent", col.frame="black", ylim=c(0,100),
                          lwd.border=10, lwd=10, min.height=10, min.width=10, lwd.grid=10, 
                          background.title = "transparent", name="Methylation", fontcolor.title="transparent",
                          ncolor=MethLevels, gradient=c("grey90","black"))
  rownames(track_meth@data) <- gsub("__"," ",rownames(track_meth@data))
  rownames(track_meth@data) <- gsub("_","-",rownames(track_meth@data))
  
  
  # CG track
  track_CG <- AnnotationTrack(range = meth_df, name = "CpG", chromosome = unique(meth_df$chromosome),
                              fontcolor.title="black", col=NULL, fill = "navy", 
                              background.title = "transparent", col.border.title = "transparent", rotate.title = F)
  
  # Sequence track
  sTrack <- SequenceTrack(sequence = genome, chromosome = unique(meth_df$chromosome), genome = genomeID,
                          add53=F , noLetters=T, fontcolor=bases_colors, name="sequence", col.border.title="transparent",cex=0.5)
  
  
  
  # Genomic axis track
  gtrack <- GenomeAxisTrack(cex=2.5,col="black",fontcolor="black", distFromAxis = 5, lwd=2, col.border.title="transparent")
  if (coord$end-coord$start<=100) { by <- 20}
  if (coord$end-coord$start >100 & coord$end-coord$start<250) { by <- 50} 
  if (coord$end-coord$start>=250) { by <- 100}
  plot_ticks <- seq(from=round(coord$start,-2), to=round(coord$end+100,-2), by=by)
  
  
  #************************************************************************************************
  track_meth_width = max(nchar(colnames(meth_df)[-c(1,2,3)]))*0.003+1
  
  width_plot <- (220 + max(nchar(colnames(meth_df[-c(1,2,3)])))*22) + ((coord$end-coord$start)*5+1000) + 50
  heigh_plot <- (ncol(meth_df[,-c(1,2,3)]))*56 + 400
  
  #************************************************************************************************
  png(filename = filename, width = width_plot, height = heigh_plot)
  plotTracks(list(track_meth,track_CG,sTrack,gtrack),
             sizes = c(ncol(meth_df[,-c(1,2,3)]),1,1,2.5),
             main=plot_title, cex.main = 3,
             cex.title = 3, title.width = track_meth_width,
             from=coord$start, to=coord$end, add53=TRUE,labelPos="below",ticksAt=plot_ticks)
  dev.off()
}






# BOXPLOTS FOR EACH POSITION ----------------------------------------------

boxplot_pos <- function(table, filename, plot_title, cloning_exp, p_label=c("pval","psign"), plot_colors) {
  
  if (cloning_exp==T) { table[,"samples"] <- table$clone }
  if (cloning_exp==F) { table[,"samples"] <- table$rep }
  
  # Parameters
  nb <- list(
    CG = length(unique(table$CG_coord)),
    collection = length(levels(table$collection)),
    group = length(levels(table$group)))
  
  list_comb <- as.list(data.frame(t(combinations(levels(table$group), k=2))))
  group_colors <- as.vector(plot_colors[1:nb$group])
  
  # Boxplot or points for clones
  if (cloning_exp==T & length(unique(table$samples))>=3) {
    label_y <- seq(length.out = length(list_comb), from=110, by=12)
    
    boxplot <- ggplot(table, aes(x = group, y = meth)) + 
      facet_grid(collection ~ CG_coord) + # grid
      geom_dotplot(aes(col=group),fill="white", stroke = 2, binaxis = 'y', stackdir='center', binwidth = 1, dotsize = 5, na.rm = T) +
      stat_summary(aes(col=group),fun = "mean", size = 1.5, shape = 18) +
      scale_y_continuous(limits = c(-10, max(label_y)+10), breaks = seq(0,100,50),  minor_breaks = seq(0, 100, 10))
  }
  if (cloning_exp==F & length(unique(table$samples))>=3) {
    label_y <- seq(length.out = length(list_comb), from=100, by=12)
    
    boxplot <- ggplot(table, aes(x = group, y = meth)) + 
      facet_grid(collection ~ CG_coord) + # grid
      geom_boxplot(aes(col=group),lwd=1, fill="transparent", na.rm=T) + 
      stat_boxplot(aes(col=group), geom = "errorbar", width=0.3, size=1, na.rm=T) +
      stat_summary(aes(col=group),fun = "mean", size = 1.5, shape = 18) +
      scale_y_continuous(limits = c(0, max(label_y)+10), breaks = seq(0,100,20), minor_breaks = seq(0, 100, 10))
  } 
  else if (length(unique(table$samples))<3) {
    label_y <- 100
    
    boxplot <- ggplot(table, aes(x = group, y = meth)) + 
      facet_grid(collection ~ CG_coord) + # grid
      geom_dotplot(aes(col=group),fill="white", stroke = 2, binaxis = 'y', stackdir='center', binwidth = 1, dotsize = 5, na.rm = T) +
      scale_y_continuous(limits = c(0, 100), breaks = seq(0,100,20), minor_breaks = seq(0, 100, 10))
  }
  
  # Boxplot 
  boxplot <- boxplot +
    theme_bw() + # theme
    labs(
      title = plot_title,
      x="Groups",
      y="Methylation percentages") +
    scale_color_manual(values=group_colors) +
    theme(
      line=element_line(size=1),
      plot.title.position = "plot",
      plot.title = element_text(hjust=0.5, size = 19, face = "bold", color = "black"), # center title
      axis.title.x.bottom = element_text(size = 16, face = "bold", vjust = -2),
      axis.title.y.left = element_text(size = 16, face = "bold", vjust = 2),
      axis.text.x.bottom = element_text(size = 15, color = "black"),
      axis.text.y.left = element_text(size = 14, color = "black"),
      strip.text.x = element_text(size = 14, color = "black"),
      strip.text.y = element_text(size = 16, color = "black"),
      panel.grid.major.x = element_blank(),
      legend.position = "none",
      panel.spacing = unit(1, "line"),
      plot.margin = unit(c(1,1,1,1),"line"),
      aspect.ratio = 3/nb$group) # ratio height/width of plots
  
  # p value format : p value
  if (p_label=="pval") { 
    boxplot <- boxplot + 
      stat_compare_means(method = "t.test", paired = FALSE, 
                         label.y = label_y, aes(label = "..p.format.."), 
                         comparisons = list_comb, na.rm = TRUE, size = 4, bracket.size=0.8) 
    }

  # p value format : significativity symbol
  if (p_label=="psign") {
    boxplot <- boxplot + 
      stat_compare_means(method = "t.test", paired = FALSE, 
                         label.y = label_y, aes(label = "..p.signif.."), 
                         comparisons = list_comb, na.rm = TRUE, size = 4, bracket.size=0.8) 
    }
 
  # save ggplot
  width_plot <- max(((450*nb$group)*nb$CG)+400, nchar(plot_title)*35)
  height_plot <- (1400*nb$collection)+400
  ggsave(filename = filename, width = width_plot, height = height_plot, dpi=300, units = c("px"),limitsize = FALSE)
}




# BOXPLOT FOR MEANS OF ALL POSITIONS --------------------------------------

boxplot_mean <- function(table, filename, plot_title, p_label=c("pval","psign"), plot_colors) {
  
  # Parameters
  nb <- list(
    collection = length(unique(table$collection)),
    group = length(levels(table$group)))
  
  list_comb <- as.list(data.frame(t(combinations(levels(table$group), k=2))))
  label_y <- seq(length.out = length(list_comb), from=100, by=12)
  group_colors <- as.vector(plot_colors[1:nb$group])
  
  
  # Boxplot 
  boxplot <- ggplot(table, aes(x = group, y = meth)) + 
    geom_boxplot(aes(col=group),lwd=1, fill="transparent", na.rm=T) + 
    stat_boxplot(aes(col=group), geom = "errorbar", width=0.3, size=1, na.rm=T) +
    stat_summary(aes(col=group),fun = "mean", size = 1.5, shape = 18) +
    theme_bw() + # theme
    scale_y_continuous(limits = c(0, max(label_y)+10), breaks = seq(0,100,20), minor_breaks = seq(0, 100, 10)) + 
    labs(
      title = plot_title,
      x="Groups",
      y="Methylation percentages") +
    scale_color_manual(values=group_colors) +
    theme(
      line=element_line(size=1),
      plot.title.position = "plot",
      plot.title = element_text(hjust=0.5,size = 19, face = "bold", color = "black"), # center title
      axis.title.x.bottom = element_text(size = 16, face = "bold", vjust = -2),
      axis.title.y.left = element_text(size = 16, face = "bold", vjust = 2),
      axis.text.x.bottom = element_text(size = 15, color = "black"),
      axis.text.y.left = element_text(size = 14, color = "black"),
      strip.text.x = element_text(size = 16, color = "black"),
      panel.grid.major.x = element_blank(),
      legend.position = "none",
      panel.spacing = unit(1, "line"),
      plot.margin = unit(c(1,1,1,1),"line"),
      aspect.ratio = 3/nb$group) # ratio height/width of plots
  
  
  # one panel per collection if more than one
  if (nb$collection>1) { 
    boxplot <- boxplot + facet_wrap(~collection) 
    } 
  
  # p value format : p value
  if (p_label=="pval") { 
    boxplot <- boxplot + 
      stat_compare_means(method = "t.test", paired = FALSE, 
                         label.y = label_y, aes(label = "..p.format.."), 
                         comparisons = list_comb, na.rm = TRUE, size = 4, bracket.size=0.8) 
    }
  
  # p value format : significativity symbol
  if (p_label=="psign") {
      boxplot <- boxplot + 
        stat_compare_means(method = "t.test", paired = FALSE, 
                           label.y = label_y, aes(label = "..p.signif.."), 
                           comparisons = list_comb, na.rm = TRUE, size = 4, bracket.size=0.8) 
      }
  
  # save ggplot
  width_plot <- max(((450*nb$group)*nb$collection)+800, nchar(plot_title)*40)
  height_plot <- 1800
  ggsave(filename = filename, width = width_plot, height = height_plot, dpi=300, units = c("px"),limitsize = FALSE)
}





# METHYLATION PROFILE PLOT FUNCTION ---------------------------------------

profile_plot <- function(table, filename, plot_title, plotType=c("proportional","condensed"), p_label=c("pval","psign"), pos_labels, plot_colors, plot_shapes){
  
  #************************************************************************************************
  # Parameters
  nb <- list(
    CG = length(unique(table$CG_nb)),
    group = length(levels(table$group)))
  group_colors <- as.vector(plot_colors[1:nb$group])
  group_shapes <- as.vector(plot_shapes[1:nb$group])
  
  
  #************************************************************************************************
  # Plot type
  if(plotType=="proportional") {
    table[,"CG"] <- table$position
    ratio <- 1/2
    width_x <- 5000
  }
  if(plotType=="condensed") {
    table[,"CG"] <- as.integer(table$CG_nb)
    ratio <- 1/(nb$CG*0.2)
    width_x <- nb$CG*500
  }
  
  #************************************************************************************************
  # Compute KW test
  KW_df <- table %>% 
    group_by(CG,group) %>%
    dplyr::mutate(meth_mean = round(mean(meth,na.rm=T),2)) # Get mean values for each group at each position
  
  KW_no_data <- KW_df[which(is.na(KW_df$meth_mean)),] # mean of meth NA reveal that no data is available for one entire group for one position
  KW_no_data <- KW_no_data %>% 
    select(CG,group) %>% 
    unique()
  KW_no_data <- KW_no_data %>% 
    group_by(CG) %>% 
    mutate("nb_group"=length(unique(group))) # Get the number of group concerned with no data

  cg_no_data <- KW_no_data[which(length(unique(table$group))-KW_no_data$nb_group ==1),"CG"]$CG # if the difference between number total of groups and number of groups with no data is 1, KW test not possible
  KW_df <- KW_df[which(KW_df$CG %in% setdiff(KW_df$CG,cg_no_data)),] # Remove CG positions for which 1 grop only have data
  
  
  KW_df <- KW_df %>% 
    group_by(CG) %>% 
    dplyr::mutate(p=kruskal.test(meth ~ group)$p.value) %>% # KW test on all values
    add_significance("p") %>% 
    select(CG,group,meth_mean,p,p.signif) %>% 
    unique()
  
  KW_df <- KW_df %>% 
    group_by(CG) %>% 
    dplyr::mutate(meth_max = max(meth_mean)) %>%  # Get max of mean values, for y position of label
    select(CG,meth_max,p,p.signif) %>% 
    unique()
  KW_df$p <- p_format(KW_df$p)
  KW_df$p.signif[which(KW_df$p.signif=="ns")] <- ""
  
  #************************************************************************************************
  # Get means of groups
  table <- table %>% 
    group_by(CG,group_ID) %>%
    dplyr::mutate(meth = round(mean(meth,na.rm=T),2)) %>% 
    select(c(CG_nb,CG_coord,position,CG,meth,group,group_ID)) %>% 
    unique()
  
  
  #************************************************************************************************
  # Change position label depending on input request
  if(pos_labels=="coordinates") {
    pos_lab <-unique(table$CG_coord)
    angle <- 90
    vjust = 0.5
    height_pos <- max(nchar(table$CG_coord))*70
  } 
  if(pos_labels=="numbers") {
    pos_lab <- seq(from=1, to=nb$CG, by=1) 
    angle <- 0
    vjust = 0
    height_pos <- 150
  }
  if(pos_labels=="none") {
    pos_lab <- rep("",nb$CG)
    angle <- 0
    vjust = 0
    height_pos <- 150
  }
  
  
  #************************************************************************************************
  # Plot
  
  plot <- ggplot(table, aes(x=CG, y=meth)) +
    geom_line(aes(col=group)) +
    geom_point(aes(col=group, shape=group), fill="white", size=2) + 
    scale_x_continuous(breaks = unique(table$CG), labels = pos_lab) +
    scale_y_continuous(limits = c(0,110), breaks = seq(0,100,20), minor_breaks = seq(0,100,10)) +
    scale_fill_manual(values=group_colors) +
    scale_color_manual(values=group_colors) +
    scale_shape_manual(values=group_shapes) +
    labs(title = plot_title,
         y = "Methylation percentages") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 13, color = "black"),
          axis.text.x = element_text(angle = angle,size = 12, color = "black", vjust = vjust),
          axis.text.y = element_text(size = 12, color = "black"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size=11),
          plot.margin = unit(c(1,1,1,1),"line"),
          aspect.ratio = ratio)
  
  
  
  
  #************************************************************************************************  
  # Kruskal Wallis significance
  if (p_label=="pval") {
    plot <- plot +
      geom_text(data = KW_df, aes(x = CG, y = meth_max, label = p), size = 3.5, vjust=-2)
  }
  
  if (p_label=="psign") {
    plot <- plot +
      geom_text(data = KW_df, aes(x = CG, y = meth_max, label = p.signif), size = 3.5, vjust=-2)
  }
  
  #************************************************************************************************
  # save ggplot
  width_plot <- width_x + 240 # margins 120x2
  height_plot <- height_pos + 2500 + 650 
  
  ggsave(filename = filename, width = width_plot, height = height_plot, dpi=600, units = c("px"))
  
}








