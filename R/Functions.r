## Main functions
####################################################
glm.predict.                   <- function(x,slope,intercept){
  y <- 1/(1+exp( -( slope*x + intercept ) ))
  return(y)
}
LASTparser_SplitForwardReverse <- function(inFile, ref_label){
  #In LAST file, detect lines with wildtype sequence idicating if read was on forward or reverse strand
  indR <- as.numeric(system(paste0("awk \'/^s/ && $5==\"-\" && $2!~/",ref_label,"/ {print NR}\' ",inFile), intern = T))
  indF <- as.numeric(system(paste0("awk \'/^s/ && $5==\"+\" && $2!~/",ref_label,"/ {print NR}\' ",inFile), intern = T))

  #Add lines before and after of the same read
  indR <- sort(c(indR-2,indR-1,indR, indR+1,indR+2))
  indF <- sort(c(indF-2,indF-1,indF, indF+1,indF+2))

  #Produce separate datasets for F or R
  data <- readLines(inFile)
  data_f <- data[indF]
  data_r <- data[indR]

  #Save files
  outFileF <- sub(pattern=".maf", replacement = "_forward.maf", inFile)
  outFileR <- sub(pattern=".maf", replacement = "_reverse.maf", inFile)

  write(x = data_f, file =  outFileF)
  write(x = data_r, file = outFileR)

  return(0)
}
LASTparser_MafToTxt            <- function(inFile_prefix, ref_label){
  inFile <- paste0(inFile_prefix, ".maf")
  outFile <- inFile_prefix

  #Parse
  wt_sequence        <- system(paste0("awk \'$2~/",ref_label,"/ {print $7}\' ",inFile, "> ", outFile, "_wtSeqs.txt" ))
  wt_initialposition <- system(paste0("awk \'$2~/",ref_label,"/ {print $3}\' ",inFile, "> ", outFile, "_pos.txt"  ) )
  read_sequence      <- system(paste0("awk \'/^s/ && $2!~/",ref_label,"/ {print $7}\' ",inFile, "> ", outFile, "_readSeqs.txt"  ))
  read_qscore        <- system(paste0("awk \'/^q/ {print $3}\' ",inFile, "> ", outFile , "_readQscore.txt"))
  read_name          <- system(paste0("awk \'/^s/ && $2!~/",ref_label,"/ {print $2}\' ",inFile, "> ", outFile, "_readName.txt" ))

  return(0)
}
LASTparser_TxtToFasta          <- function(inFile_prefix, length_sequencing, ORF_ini, ORF_end){
  ## Filenames to load
  file.wt   <- paste0(inFile_prefix,"_wtSeqs.txt")
  file.seqs <- paste0(inFile_prefix,"_readSeqs.txt")
  file.ini  <- paste0(inFile_prefix,"_pos.txt")
  file.qual <- paste0(inFile_prefix,"_readQscore.txt")
  file.name <- paste0(inFile_prefix,"_readName.txt")

  #Create files to store results
  OutFileSequence <- paste0(inFile_prefix,"_FakeAli.fasta")
  if(file.exists(OutFileSequence)){file.remove(OutFileSequence)}
  OutFileQuality <- paste0(inFile_prefix,"_FakeAliQ.fasta")
  if(file.exists(OutFileQuality)){file.remove(OutFileQuality)}

  ## Load files
  wt <- scan(file.wt, what=character(), quiet = T)
  if(length(wt)==0){next()}
  wt <- toupper(wt);  wts <- do.call(strsplit,list(wt,split=""))
  seqs    <- scan(file.seqs, what=character(), quiet = T)
  seqs    <- toupper(seqs);   seqss <-  do.call(strsplit,list(seqs,split=""))
  iniPos  <- scan(file.ini, what=numeric(),quiet = T)
  quality <- scan(file.qual, what=character(), quiet = T, sep="\n")
  quality <- do.call(strsplit, list(quality, split=""))
  names   <- scan(file.name, what=character(), quiet = T)

  length_seq <- ORF_end-ORF_ini+1
  for(i in seq_along(wt)){

    p.in <- which(wts[[i]]!="-")
    seq.aux <- seqss[[i]][p.in]
    quality.aux <- quality[[i]][p.in]

    #Fill initial positions (if not read) with gaps
    if(iniPos[i]>0){
      beg <- rep("-",iniPos[i])
      beg_quality <- rep("~", iniPos[i])
    }else{
      beg <- NULL; beg_quality <- NULL
    }

    #Fill final positions (if not read) with gaps
    largo <-  iniPos[i] + length(seq.aux)
    if( largo <length_sequencing){
      end <- rep("-",length_sequencing - largo )
      end_quality <- rep("~", length_sequencing-largo)
    }else{
      end <- NULL; end_quality <- NULL
    }

    #Join beggining - sequence - end
    seq.out     <- c(beg, seq.aux, end)
    quality.out <- c(beg_quality,quality.aux,end_quality)

    #Remove uninteresting positions (out of the ORF)
    seq.out     <- seq.out[ORF_ini:ORF_end]
    quality.out <- quality.out[ORF_ini:ORF_end]


    #* Qscore value for gaps given by nearest nucleotides
    quality.out.aux <- quality.out
    ind <- which(seq.out=="-")
    for(j in ind){
      #skip if there are only gaps before or after
      start <- seq.out[1:j]
      end   <- seq.out[j:length_seq]

      #Replace ~ by the mean of the next non gaps qualities:
      #-if all positions before are gaps, replace by the first Qscore no gap after
      if(all(start == "-")){
        nend   <- which(end[-1]!="-")[1]
        quality.out[j] <- quality.out[j+nend]
        next()
      }
      #-if all positions after are gaps, replace by the first Qscore no gap before
      if(all(end == "-")){
        nstart <- which(rev(start)[-1]!="-")[1]
        quality.out[j] <- quality.out[j-nstart]
        next()
      }
      #-if it has nucleaotides before and after, replace by the mean of the Qscores of nearest neighbours
      nstart <- which(rev(start)[-1]!="-")[1]
      nend   <- which(end[-1]       !="-")[1]

      #Mean of probabilities
      aux.qual       <- quality.out.aux[c(j-nstart, j+nend)]
      aux.qual.q     <- mapvalues(x=aux.qual, from=ascii$Symbol, to=ascii$Q, warn_missing = F)
      aux.qual.q     <- as.numeric(aux.qual.q)
      mean.qscore    <- round(mean(aux.qual.q))

      aux.symbol     <- mapvalues(x=mean.qscore, from=ascii$Q, to=(ascii$Symbol), warn_missing = F)
      quality.out[j] <- levels(ascii$Symbol)[aux.symbol]
    }

    #Write in output files
    cat(">", names[i],"\n", sep="", append = T, file=OutFileSequence)
    cat(seq.out,"\n", sep="", append = T, file=OutFileSequence)

    cat(">",names[i],i,"\n", sep="", append = T, file=OutFileQuality)
    cat(quality.out,"\n", sep="", append = T, file=OutFileQuality)

  }
  return(0)
}
parse_maf_to_alignment         <- function(filename, LAST_reference_label, length_sequencing, ORF_ini,ORF_end, auxiliary_files=NULL){


  #file names for training set
  filename_wt_sufix <-  sub(pattern = ".maf", replacement = "",x=filename)
  filename_wt_F <-  paste0(filename_wt_sufix, "_forward")
  filename_wt_R <-  paste0(filename_wt_sufix, "_reverse")

  #Split reads in forward and reverse
  LASTparser_SplitForwardReverse(inFile = filename,
                                 ref_label = LAST_reference_label)
  #Parse files to tables
  LASTparser_MafToTxt(inFile_prefix = filename_wt_F,ref_label = LAST_reference_label)
  LASTparser_MafToTxt(inFile_prefix = filename_wt_R,ref_label = LAST_reference_label)

  #Parse  prior_muta tables to FASTAs | construct 'fake' MSA
  LASTparser_TxtToFasta(filename_wt_F, length_sequencing = length_sequencing,ORF_ini =ORF_ini ,ORF_end = ORF_end )
  LASTparser_TxtToFasta(filename_wt_R, length_sequencing = length_sequencing,ORF_ini =ORF_ini ,ORF_end = ORF_end )

  if(!is.null(auxiliary_files)){
    if(!dir.exists(auxiliary_files)){dir.create(auxiliary_files)}
    system(paste0("mv ", filename_wt_F,".maf ", filename_wt_R,".maf ",auxiliary_files))
  }else{
    system(paste0("rm ", filename_wt_F,".maf ", filename_wt_R,".maf "))
  }
  return(0)

}

calculate_prior_errors         <- function(inFile_prefix, save=F){
  prior_error <- read.table(paste0(inFile_prefix, "_countsPNQ.txt"), header = T)
  prior_error <- prior_error%>% group_by(position, nucleotide)%>%
    dplyr::summarise(counts=sum(counts))%>%
    mutate(isWT=0) %>%ungroup()

  prior_error <- prior_error%>% group_by(position)%>%
    dplyr::mutate(prior.error = counts/sum(counts))

  for(i in seq_along(reference_sequence)){ prior_error$isWT[which(prior_error$position==i & prior_error$nucleotide==reference_sequence[i])] <- 1}
  prior_error <- prior_error[,c("position", "nucleotide","prior.error")] %>% ungroup()
  if(save){write.table(prior_error, file = paste0(inFile_prefix, "_prior_errors.txt"), row.names = F)}
  return(prior_error)
}
calculate_prior_mutations      <- function(rates.matrix, mean.n.mut, wildtype_sequence, save=F, filename="tablePriorMutations.txt"){
  composition_wt <- c(table(wildtype_sequence), length(wildtype_sequence))
  names(composition_wt)[5] <- "-"

  mutations_rate              <- apply(rates.matrix, 2, sum, na.rm=T)
  expected_mutations_perbase  <- composition_wt*mutations_rate
  normalization_factor        <- sum(expected_mutations_perbase)/mean.n.mut

  expected_mutation_rate <- rates.matrix / normalization_factor
  expected_mutation_rate <- melt(expected_mutation_rate, varnames = c("from","to"))
  expected_mutation_rate <- expected_mutation_rate[!is.na(expected_mutation_rate$value),]

  colnames(expected_mutation_rate) <-  c("wt.base","nucleotide","p_mutation")
  if(save){write.table(expected_mutation_rate, file = filename, row.names = F)}
  return(expected_mutation_rate)
}

count_nucleotides_per_qscore   <- function(inFile_prefix){
  #Load reads from fasta
  sequences <- readLines(paste0(inFile_prefix, "_FakeAli.fasta"))
  quality   <- readLines(paste0(inFile_prefix,"_FakeAliQ.fasta"))

  #keep relevant lines of files
  lines_data <- seq(from = 2,to = length(sequences), by=2)
  sequences  <- sequences[lines_data]
  quality    <- quality[lines_data]

  #convert reads into a list of vectors
  sequences.list <- strsplit(sequences, split="")
  quality.list   <- strsplit(quality, split="", fixed = T)

  #convert list of vectors into a matrix
  sequences.matrix <- do.call(rbind,sequences.list)
  quality.matrix   <- do.call(rbind,quality.list)

  #For relevant positions in the alignment, count occurences of nucleotide per position and Qscore
  for(i in 1:ncol(sequences.matrix)){
    df <- data.frame(nucleotide=sequences.matrix[,i],quality=quality.matrix[,i])
    df.counts <- df%>% group_by(nucleotide,quality)%>% tally(name = "counts")%>%mutate(position=i)
    write.table(x=df.counts[,c("position", "nucleotide","quality","counts")],
                file=paste0(inFile_prefix,"_countsPNQ.txt"),
                append = i!=1, col.names = i==1,row.names = F)
  }
}
fit_logistic_regression        <- function(inFile_prefix, reference_sequence, prior_error, prior_mutation){
  data <- read.table(paste0(inFile_prefix, "_countsPNQ.txt"), header=T)
  #Pre-editing data
  data <- dplyr::as_tibble(data)%>%
    #dplyr::filter(nucleotide!="-")%>%                                             # remove gaps
    dplyr::mutate(qualityQ=mapvalues(quality, from = ascii$Symbol, to=ascii$Q, warn_missing = F))%>% # map Qscores symbols to numbers
    dplyr::mutate(qualityQ= as.numeric(levels(qualityQ)[qualityQ]))%>%             # Qscore numbers to numeric
    select(-quality)%>%                                                            # remove Qscore symbols column
    dplyr::mutate(wt.base = reference_sequence[position])                                    # Add wildtype base
  data$nucleotide <- droplevels(data$nucleotide)                                   # remove - from levels of nucleotide

  ## Wildtype matrix: keep rows with wildtype reads
  data_wt <- data %>%
    filter(nucleotide==wt.base)%>%
    select(-nucleotide)
  colnames(data_wt)[colnames(data_wt)=="counts"] <- "counts.wt"

  ## Data with mutations (errors)
  data_mut <- data %>%                                             #start from data
    filter(nucleotide!=wt.base)%>%                                 #keep only mutations
    full_join(expand(.,position, nucleotide, qualityQ), by = c("position", "nucleotide", "qualityQ"))%>%         #complete all combinations position - nucleotide - qualityQ
    mutate(wt.base = reference_sequence[position]) %>%             #fill wildtype base for missing values (new rows)
    full_join(data_wt,   by=c("position", "qualityQ","wt.base"))%>%  # add wildtype info in new columns
    left_join(prior_mutation, by=c("wt.base", "nucleotide"))%>%      # add prior of being mutated
    left_join(prior_error,    by=c("position", "nucleotide"))%>%     # add prior of being an error
    filter(nucleotide!=wt.base)%>%                                   # remove wldtype rows
    mutate(counts=replace_na(counts,0))%>%                           # fill NA with 0
    mutate(counts.wt=replace_na(counts.wt,0))%>%                     # fill NA with 0
    mutate(prior.error=replace_na(prior.error,0))                    # fill NA with 0

  ## Count number of errors/good reads by position and nucleotide
  total_counts <- data_mut %>%
    group_by(position, nucleotide)%>%
    dplyr::summarise(total.counts.mut=sum(counts,na.rm=T),
                     total.counts.wt=sum(counts.wt,na.rm=T))%>%
    mutate(total.counts=total.counts.mut+total.counts.wt)
  data_mut <- full_join(data_mut, total_counts,    by=c("position", "nucleotide"))

  ## reweight counts by prior probabilities
  data_mut <- data_mut%>%
    mutate(pc = p_mutation / (p_mutation+prior.error),                            #prob of being correct
           pi = prior.error / (p_mutation+prior.error))%>%                        #prob of being an error
    mutate(counts.scaled    = (counts / total.counts.mut * pi * total.counts ),   #counts errors re-weighted
           counts.wt.scaled = (counts.wt / total.counts.wt * pc * total.counts )) #counts wildtype re-weighted

  ## Fit data:
  data_fits <- data_mut %>%
    select(position, nucleotide)%>%
    distinct(position, nucleotide)%>%
    arrange(position, nucleotide)
  data_fits$regular <- rep(NA, nrow(data_fits))
  data_fits$priors  <- rep(NA, nrow(data_fits))
  n_i <- nrow(data_fits)
  cat("\n Fitting \n")
  pb = txtProgressBar(min = 0, max = nrow(data_fits), initial = 0, width=10)
  for (i in 1:nrow(data_fits)){
    setTxtProgressBar(pb,i)
    #Keep data from this position & nucleotide
    aux_df <- data_mut %>%
      filter(position==data_fits$position[i] & nucleotide ==data_fits$nucleotide[i])%>%
      select(qualityQ,counts, counts.wt,counts.scaled,counts.wt.scaled) %>%
      mutate(tot=counts+counts.wt, proportion.wt=counts.wt/tot)%>%
      mutate(tot.scaled=counts.scaled+counts.wt.scaled, proportion.wt.scaled=counts.wt.scaled/tot.scaled)

    qval  <- aux_df$qualityQ
    yval  <- aux_df$proportion.wt
    yvals <- aux_df$proportion.wt.scaled

    #Naive fit (no reweight)
    data_fits$regular[i] <-
      list( glm(yval[!is.na(yval)]~ qval[!is.na(yval)],family = "binomial"))

    #Corrected fit by prior data
    if(! (all(is.na(aux_df$counts.wt.scaled)) | all(is.na(aux_df$counts.scaled))) ){
      data_fits$priors[i]  <-
        list(glm(yvals[!is.na(yvals)]~ qval[!is.na(yvals)],family = "binomial"))
    }
  }

  save(data_fits, file=paste0(inFile_prefix, "_apriori_fits.Rdata"))
  save(data_mut, file=paste0(inFile_prefix, "_apriori_data.Rdata"))
}
evaluate_fits                  <- function(inFile_prefix, data_fits, reference_sequence, ascii){
  frequencies <- read.table(paste0(inFile_prefix,"_countsPNQ.txt"), header=T)               # Load frequency data

  frequencies <- frequencies%>%
    mutate(qualityQ=mapvalues(quality, from = ascii$Symbol, to=ascii$Q, warn_missing = FALSE))%>%     # Convert Qscores to numerical representation
    mutate(qualityQ= as.numeric(levels(qualityQ)[qualityQ]))%>%                                                                                                  # de Qscore a probabilidades (p_minION)
    group_by(position)%>%
    mutate(p_error_minION=10^(-qualityQ/10))%>%               # Convert Qscores to p_error
    ungroup()%>%
    mutate(wt.base=reference_sequence[position])%>%           # Column with reference nucleotide
    mutate(isWT = nucleotide==wt.base) %>%                    # Logical column oindicating if read is reference or different nucleotide
    select(-quality)                                          # Remove Qscore column

  frequencies$p_right_priors_model <- rep(NA)
  frequencies$p_right_null_model  <- rep(NA)
  cat("\n Evaluating for", inFile_prefix,"\n")
  pb = txtProgressBar(min = 0, max = nrow(frequencies), initial = 0, width=10)
  for(i in 1:nrow(frequencies)){
    setTxtProgressBar(pb,i)
    # skip if it is a wildtype or a gaps:
    if(frequencies$isWT[i]==1){next()}
    #if(frequencies$nucleotide[i]=="-"){next()}
    # match read with fit to be used
    index  <- which(data_fits$position== frequencies$position[i] & data_fits$nucleotide==frequencies$nucleotide[i])
    aux_df <- data.frame(qualityQ=frequencies$qualityQ[i])
    if(is.na(data_fits$priors[index][[1]])){next()}

    #fit models
    f_prior <- data_fits$priors[index][[1]]
    f_regular <- data_fits$regular[index][[1]]
    frequencies$p_right_priors_model [i] <- glm.predict.(x=aux_df[[1]],slope = f_prior$coefficients[2], intercept = f_prior$coefficients[1])
    frequencies$p_right_null_model [i]   <- glm.predict.(x=aux_df[[1]],slope = f_regular$coefficients[2], intercept = f_regular$coefficients[1])
  }

  #Save results
  write.table(frequencies, file = paste0(inFile_prefix, "_FrequenciesModelPriors.txt"), row.names = F)
}
SplitSequencesInFiles          <- function(inFile_prefix, reference_sequence, columns=NULL){
  # Load data
  #data of qscores
  qtable <- read.table(paste0(inFile_prefix,"_FrequenciesModelPriors.txt"), header = T)
  qtable <- qtable%>% select(position,nucleotide,qualityQ,p_right_priors_model,p_right_null_model)
  sequences <- readLines(paste0(inFile_prefix, "_FakeAli.fasta"))
  quality   <- readLines(paste0(inFile_prefix,"_FakeAliQ.fasta"))

  # Generate folder to store
  bc <- strsplit(inFile_prefix, split = "_")[[1]][1]
  if(!dir.exists(bc)){dir.create(bc)}

  #Parse into matrix
  lines_data <- seq(from = 2,to = length(sequences), by=2)
  names_sequences <-  sub(pattern = ">", replacement = "", x=sequences[lines_data-1])
  sequences <- sequences[lines_data]
  quality <- quality[lines_data]

  sequences.list <- strsplit(sequences, split="")
  quality.list <- strsplit(quality, split="", fixed = T)

  sequences.matrix <- do.call(rbind,sequences.list)
  quality.matrix <- do.call(rbind,quality.list)

  if(!is.null(columns)){
    sequences.matrix <- sequences.matrix[,columns]
    quality.matrix <- quality.matrix[,columns]

  }
  # Make a table with relevant information to store
  nc <- ncol(sequences.matrix)
  for(i in 1:length(sequences)){
    df <- data.frame( nucleotide= sequences.matrix[i,], q = quality.matrix[i,], position=1:nc) %>%
      mutate(qualityQ=plyr::mapvalues(q, from = ascii$Symbol, to=ascii$Q, warn_missing = F)) %>%
      mutate(qualityQ= as.numeric(levels(qualityQ)[qualityQ]), isWT = 0) %>%
      mutate(p_error_minION=10^(-qualityQ/10))%>%select(-q)%>%
      mutate(isWT = as.numeric(nucleotide==reference_sequence[position]))%>%
      left_join(qtable, by = c("nucleotide","position","qualityQ"))
    df=left_join(df,qtable,by=c("nucleotide", "position", "qualityQ", "p_right_priors_model", "p_right_null_model"))

    write.table(df, file=paste0(bc,"/",bc,"_", names_sequences[i],".txt"), row.names = F)
  }
}



LRminion_train_model    <- function(filename_wildtype, LAST_reference_label, length_sequencing, ORF_ini,ORF_end, auxiliary_files=NULL, verbose=TRUE){

  if(verbose){cat("Processing", filename_wildtype, "\n")}
  #file names for training set
  filename_wt_sufix <-  sub(pattern = ".maf", replacement = "",x=filename_wildtype)
  filename_wt_F <-  paste0(filename_wt_sufix, "_forward")
  filename_wt_R <-  paste0(filename_wt_sufix, "_reverse")

  if(verbose){cat("\t Parsing maf file \n")}
  parse_maf_to_alignment(filename_wildtype, LAST_reference_label, length_sequencing, ORF_ini,ORF_end, auxiliary_files=NULL)

  if(verbose){cat("\t Counting nucleotides \n")}
  #Count occurrences per position/nucleotide
  count_nucleotides_per_qscore(filename_wt_F)
  count_nucleotides_per_qscore(filename_wt_R)

  if(verbose){cat("\t Calculating p_prior-error \n")}
  # Calculate priors errors
  priors_errors_F <- calculate_prior_errors(inFile_prefix = filename_wt_F, save = T)
  priors_errors_R <- calculate_prior_errors(inFile_prefix = filename_wt_R, save = T)

  if(verbose){cat("\t Calculating p_prior-right\n")}
  ## Load/calculate information of priors probabilities
  prior_mutation <- calculate_prior_mutations(rates.matrix = mutationratePCR, mean.n.mut = mean.n.mut, wildtype_sequence = reference_sequence, save = T)

  #Fit wt
  fit_logistic_regression(filename_wt_F, reference_sequence=reference_sequence, prior_error = priors_errors_F, prior_mutation = prior_mutation)
  fit_logistic_regression(filename_wt_R, reference_sequence=reference_sequence, prior_error = priors_errors_R, prior_mutation = prior_mutation)

  return(0)

}
LRminion_evaluate_model <- function(filename_to_evaluate,filename_wildtype, LAST_reference_label, length_sequencing, ORF_ini,ORF_end,reference_sequence, auxiliary_files=NULL, verbose=T){

  if(verbose){cat("\t Loading fits data \n")}
  filename_wt_sufix <-  sub(pattern = ".maf", replacement = "",x=filename_wildtype)
  filename_wt_F <-  paste0(filename_wt_sufix, "_forward")
  filename_wt_R <-  paste0(filename_wt_sufix, "_reverse")

  #Load fits
  load(paste0(filename_wt_F, "_apriori_fits.Rdata"))
  fitsF <- data_fits
  rm(data_fits)

  load(paste0(filename_wt_R, "_apriori_fits.Rdata"))
  fitsR <- data_fits
  rm(data_fits)

  #Prepare files
  for(i in seq_along(filename_to_evaluate)){
    if(verbose){cat("Processing", filename_to_evaluate[i], "\n")}

    #file names for training set
    filename <- filename_to_evaluate[i]
    filename_sufix <-  sub(pattern = ".maf", replacement = "",x=filename)
    filename_F <-  paste0(filename_sufix, "_forward")
    filename_R <-  paste0(filename_sufix, "_reverse")

    if(verbose){cat("\t Parsing maf file \n")}
    #Parse files
    parse_maf_to_alignment(filename, LAST_reference_label, length_sequencing, ORF_ini, ORF_end, auxiliary_files)

    if(verbose){cat("\t Counting nucleotides \n")}
    #Count occurrences per position/nucleotide
    count_nucleotides_per_qscore(filename_F)
    count_nucleotides_per_qscore(filename_R)

    if(verbose){cat("\t Evaluating fits \n")}
    evaluate_fits(inFile_prefix = filename_F, data_fits = fitsF, reference_sequence = reference_sequence,ascii=ascii)
    evaluate_fits(inFile_prefix = filename_R, data_fits = fitsR, reference_sequence = reference_sequence,ascii=ascii)

    if(verbose){cat("\t Storing data into folders \n")}
    SplitSequencesInFiles(filename_F,reference_sequence)
    SplitSequencesInFiles(filename_R,reference_sequence)


  }

  return(0)
}


