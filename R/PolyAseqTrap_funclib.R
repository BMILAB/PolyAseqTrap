#' @importFrom Rsamtools BamFile ScanBamParam scanBam
#' @importFrom stringr str_extract_all str_sub str_count str_extract
#' @importFrom plyr rbind.fill
#' @importFrom Biostrings DNAStringSet reverseComplement getSeq
#' @importFrom dplyr %>% group_by  summarise bind_rows  rename mutate
#' @importFrom GenomicRanges findOverlaps reduce GRanges makeGRangesFromDataFrame
#' @importFrom IRanges extractList IRanges
#' @importFrom BiocGenerics start
#' @importFrom tibble tibble
#' @importFrom outliers scores
#' @importFrom pbmcapply pbmcmapply
#' @importFrom limma strsplit2
#' @importFrom methods as
NULL



options(stringsAsFactors=F)

#' A FindPTA function
#'
#' This function allows you to identify and quantify polyA sites from different
#' 3'seq data.
#'
#' @name FindPTA
#' @param bam The file name of the 'BAM' file to be processed.
#' @param yieldSize Number of records to yield each time the file is read from
#'                   with Rsamtools::scanBam.
#' @param reverse If TRUE, reverse complement the reads from the BAM file.
#'                BAM files store reads mapping to the minus strand as though they
#'                are on the plus strand. Default is FALSE.
#' @param bsgenome A BSgenome object for reference genome. e.g. 'BSgenome.Hsapiens.UCSC.hg38'.
#' @param threeUTRregion A GRanges object for 3'UTR region, will use level-8 reads in identification of polyA sites.
#'                       Recommended to use if input reads lack polyA tails.
#' @param cutoffCount Minimum count threshold for PACs identified at the V8 level
#'                    to be used for polyA determination.
#' @param ext3UTRlen An integer denoting the extended length of 3'UTR region. Default is 1000 nt.
#' @param isDRS Set to TRUE if input is Direct RNA Sequencing (DRS) data.
#'              Default is FALSE.
#' @param adjust.chr  Set to TRUE to add 'chr' prefix if input chromosome names
#'                    in BAM file lack 'chr' prefix but the genome sequence includes it.
#'                    Default is FALSE.
#' @param run.quantify  If TRUE (default), quantify polyA sites.
#' @param d A distance to group nearby polyA sites into a polyA site cluster (PAC). Default is 24 nt.
#' @param poly Specifies the type of polyA tail at the 3' end of the sequence.
#'             Must be either 'T' (for a poly(T) stretch, e.g., TTTTTT-like tails)
#'             or 'A' (for a poly(A) stretch, e.g., AAAAAA-like tails).
#'             Default is 'T'.
#' @return A list containing two elements about alignment results and PACs.
#' \itemize{
#'   \item \code{pa.table}: A \code{data.frame} with the following columns:
#'     \describe{
#'       \item{\code{\bold{readName}}}{Name of the read.}
#'       \item{\code{\bold{flag}}}{Bitwise FLAG.}
#'       \item{\code{\bold{cigar}}}{CIGAR string.}
#'       \item{\code{\bold{chr}}}{The name of the chromosome.}
#'       \item{\code{\bold{start}}}{The starting position of the feature in the chromosome.}
#'       \item{\code{\bold{end}}}{The ending position of the feature in the chromosome.}
#'       \item{\code{\bold{strand}}}{The strand direction of the feature in the chromosome.}
#'       \item{\code{\bold{coord}}}{The coordinate of the potential polyA site.}
#'       \item{\code{\bold{seq}}}{The sequence of the input trimmed read.}
#'       \item{\code{\bold{softClipFragment}}}{The fragment of the soft-clipping.}
#'       \item{\code{\bold{trimmed_seq}}}{The trimmed polyA tail.}
#'       \item{\code{\bold{unmapped_seq}}}{The sequence in the original read that cannot be aligned to the reference genome.}
#'       \item{\code{\bold{reference_seq}}}{The reference sequence around the potential polyA site.}
#'       \item{\code{\bold{ref_match_trimmed}}}{The part of the reference sequence corresponding to \code{trimmed_seq}.}
#'       \item{\code{\bold{ref_match_unmapped}}}{The part of the reference sequence corresponding to \code{unmapped_seq}.}
#'       \item{\code{\bold{is_Arich}}}{Whether it is an A-rich polyA site.}
#'       \item{\code{\bold{A_number}}}{The total number of adenines (A) within the reference sequence, spanning 10 bp upstream and downstream of the polyA site.}
#'       \item{\code{\bold{diff_a}}}{The total number of non-genomic As at the 3’ end.}
#'       \item{\code{\bold{level}}}{The subclass or confidence level of the aligned read.}
#'       \item{\code{\bold{class}}}{The category of the aligned read.}
#'       \item{\code{\bold{use.as.count}}}{Whether the read is used in the quantification of PACs.}
#'       \item{\code{\bold{PAid}}}{The identity of PACs.}
#'     }
#'   \item \code{pa.coord}: A \code{data.frame} with the following columns:
#'     \describe{
#'       \item{\code{\bold{seqnames}}}{The name of the chromosome of PACs.}
#'       \item{\code{\bold{start}}}{The starting position of PACs in the chromosome.}
#'       \item{\code{\bold{end}}}{The ending position of PACs in the chromosome.}
#'       \item{\code{\bold{width}}}{The width of PACs.}
#'       \item{\code{\bold{coord}}}{The representative coordinate of PACs.}
#'       \item{\code{\bold{coord_readName}}}{The name of the read supporting the representative polyA site.}
#'       \item{\code{\bold{coord_unmappedA}}}{The total number of non-genomic As at the 3’ end of the representative read.}
#'       \item{\code{\bold{coord_level}}}{The confidence level of the representative polyA site.}
#'       \item{\code{\bold{UPA_coord}}}{The coordinates of aligned reads in corresponding PACs.}
#'       \item{\code{\bold{UPA_readName}}}{The name of aligned reads in corresponding PACs.}
#'       \item{\code{\bold{UPA_level}}}{The confidence level label of aligned reads in corresponding PACs.}
#'       \item{\code{\bold{PAid}}}{The identity of PACs.}
#'       \item{\code{\bold{total.count}}}{The total counts of PACs based on all categories of aligned reads.}
#'       \item{\code{\bold{count_with.polyA}}}{The total counts of PACs based on C1-C3 of aligned reads.}
#'     }
#' }
#' @export
#'





FindPTA <- function(bam=NULL,yieldSize=26000000,reverse=FALSE,
                          bsgenome=NULL,d=24,poly="T",adjust.chr=FALSE,
                          threeUTRregion=NULL,
                          cutoffCount=5,ext3UTRlen=1000,isDRS=FALSE,
                          run.quantify=TRUE){
  ## convert bam into data.frame
  if((poly!="T" & poly!="A")){
    stop("Please make sure input of poly is `T` or `A`")
  }
  print("The program is reading the BAM files.")
  bamfile <- Rsamtools::BamFile(bam,yieldSize = yieldSize)
  param <- Rsamtools::ScanBamParam(reverseComplement = reverse,
                 what = c('qname','rname','flag','pos','cigar','seq'))
  aln.raw <- Rsamtools::scanBam(bamfile,param = param)

  aln <- as.data.frame(do.call(cbind,lapply(aln.raw[[1]], as.character)))
  rm(aln.raw)
  gc()


  ## Process bam files to obtain trimmed reads,
  # soft-clipped fragment, and reference sequences
  aln <- process_data(aln=aln,bsgenome=bsgenome,adjust.chr=adjust.chr,poly=poly)
  gc()


  aln.tail <- subset(aln,unmapreads_width>0)
  aln.notail <- subset(aln,unmapreads_width<=0)
  rm(aln)


  print("Detecting reads with polyA tail and adjust polyA coordinate")
  aln.tail$statue.unmapped.reference <- (aln.tail$unmapreads == aln.tail$ref_match_unmapped)
  index.gene.a <- which(aln.tail$statue.unmapped.reference =="FALSE")
  index.pure.a <- grep("^A+$",aln.tail$unmapreads)
  index.pure.a.c <- intersect(index.gene.a,  index.pure.a )
  aln.tail.pure <- aln.tail[index.pure.a.c,]
  aln.tail.nopure <- aln.tail[-index.pure.a.c,]

  print("C1: processing reads with 100% An tails")
  aln.result <- process_tail_pure(aln.tail.pure=aln.tail.pure)
  rm(aln.tail.pure)



  print("C2: processing reads with !100% An tails")
  aln.tail.nopure$priority <- ""
  aln.tail.nopure$ref_a <- apply(nchar(str_extract_all(aln.tail.nopure$ref_match_unmapped,pattern = 'A+$',simplify = T)),MARGIN = 1,max)+
    apply(nchar(str_extract_all(aln.tail.nopure$reference_rest,pattern = '^A+',simplify = T)),MARGIN = 1,max)

  index.in <- which(is.infinite(aln.tail.nopure$ref_a))
  aln.tail.nopure$ref_a[index.in]<-0
  aln.tail.nopure$unmppreads_a <- apply(nchar(str_extract_all(aln.tail.nopure$unmapreads,pattern = 'A+$',simplify = T)),MARGIN = 1,max)
  aln.tail.nopure$diff_a <- aln.tail.nopure$unmppreads_a - aln.tail.nopure$ref_a
  index.pure.a <- grep("^A+$",aln.tail.nopure$trimseq)
  index.nopure.a <- setdiff(c(1:nrow(aln.tail.nopure)),index.pure.a)
  index.with.soft <- which(aln.tail.nopure$softClipLen>0)
  index.without.soft <- which(aln.tail.nopure$softClipLen==0)
  index.soft.pure.a <- intersect( index.pure.a,index.with.soft)
  aln.tail.nopure.v1 <- aln.tail.nopure[index.soft.pure.a,]
  index.soft.nopure.a  <- intersect( index.nopure.a,index.with.soft)
  aln.tail.nopure.v2 <- aln.tail.nopure[index.soft.nopure.a,]

  index.nosoft.nopure.a  <- intersect( index.nopure.a,index.without.soft)
  aln.tail.nopure.v3 <- aln.tail.nopure[ index.nosoft.nopure.a ,]

  #=========[1.2.1] aln.tail.nopure: all soft mapping from reference genome, is IP
  #----IP---------------------------------------
  index.nosoft.pure.a <- intersect( index.pure.a,index.without.soft)
  if(length(index.nosoft.pure.a )>0){
    print("Excluding reads with perfectly matched polyA tails but polyA tail actually from reference genome")
    aln.tail.nopure.v4 <- aln.tail.nopure[index.nosoft.pure.a,]
    aln.tail.nopure.v4 <- process_tail_nopure_v4(aln.tail.nopure.v4=aln.tail.nopure.v4)
    aln.result <- plyr::rbind.fill(aln.result,aln.tail.nopure.v4)
  }


  #=========[1.2.4] aln.tail.nopure.v3: with  no-soft and no-pure A tial
  if(nrow(aln.tail.nopure.v3)>0){
    aln.result <- process_tail_nopure_v3(aln.tail.nopure.v3=aln.tail.nopure.v3,
                                         aln.result=aln.result )
  }



  #=========[1.2.2] aln.tail.nopure.v1: with soft and pure A tial
  if(nrow(aln.tail.nopure.v1)>0){
    print("Detecting V5 polyA sites")
    aln.result  <- process_tail_nopure_v1(aln.tail.nopure.v1=aln.tail.nopure.v1,
                                          aln.result=aln.result )
  }
  #=========[1.2.3] aln.tail.nopure.v2: with soft and no-pure A tial
  if(nrow(aln.tail.nopure.v2)>0){
    print("Detecting V6 polyA sites")
    aln.result  <- process_tail_nopure_v2(aln.tail.nopure.v2=aln.tail.nopure.v2,
                                          aln.result=aln.result )
  }


  #=================For notail Part============================================
  print("C3: reads without polyA tails")
  if(nrow(aln.notail)>0){
    #print("Condition 3: M3_Notail_IP")
    aln.notail$ref_a <- apply(nchar(str_extract_all(aln.notail$reference_sequence,pattern = '^A+',simplify = T)),MARGIN = 1,max)
    index.in <- which(is.infinite(aln.notail$ref_a))
    aln.notail$ref_a[ index.in]<-0
    aln.notail$priority <- "M3_NoTail"
    aln.notail$priority[which(aln.notail$ref_a>=6)] <- "M3_IP"
    aln.result <- plyr::rbind.fill(aln.result,aln.notail)
    rm(aln.notail)
  }

  #aln.result$level <- aln.result$priority
  aln.result$level <- factor(aln.result$priority,
                             levels=c("M1_V1", "M2_C3_V1","M2_C3_V3","M2_C3_V4","M2_C1_V1","M2_C2_V3","M2_C3_NoTail",
                                      "M1_V2", "M2_C4_IP","M2_C3_V2" ,"M2_C3_V5"  ,"M2_C1_V3" ,
                                      "M2_C1_Other" , "M2_C1_V2","M2_C2_NoTail", "M2_C2_IP","M2_C2_V6"   ,
                                      "M2_C2_V7", "M2_C2_Other","M2_C2_V4" ,
                                      "M2_C2_V5" , "M3_NoTail"  , "M3_IP" , "M2_C3_IP" ),
                             labels=c("V1","V2","V3","V4","V5","V6","V7",
                                      rep("Count",17)))

  ### Adjust SNP
  print("Identifying potential polyA sites affected by SNPs")
  aln.result <- correctSNP(reads=aln.result)

  if(isDRS){
    index <- which(aln.result$strand=="-")
    if(length(index)>0){
      aln.result$strand[index] <- "+"
      aln.result$strand[-index] <- "-"
      aln.result$coord[index] <- aln.result$end[index]
      aln.result$coord[-index] <- aln.result$start[-index]
      aln.result$flag[index] <- "16"
      aln.result$flag[-index] <- "0"
    }else{
      aln.result$strand <- "-"
      aln.result$coord <- aln.result$start
      aln.result$flag <- "0"
    }

  }




  ### Add A rich information
  print("Searching for A-rich fragments surrounding polyA sites")
  aln.result <- Arich(reads=aln.result,bsgenome =bsgenome )
  if(length(threeUTRregion)>0){
    print("C3: detecting V8 polyA sites")
    aln.result <- recycle_v8(aln.result ,cutoffCount=cutoffCount,threeUTRregion=threeUTRregion)
    aln.result$level <- factor(aln.result$level,
                               levels=c("V1","V2","V3","V4","V5","V6","V7","V8","Count"))
    aln.result$class <- factor(aln.result$level,
                               levels=c("V1","V2","V3","V4","V5","V6","V7","V8","Count"),
                               labels= c("C1","C1","C2","C2","C1","C2","C2","C3","Count"))
  }else{
    aln.result$level <- factor(aln.result$level,
                               levels=c("V1","V2","V3","V4","V5","V6","V7","Count"))
    aln.result$class <- factor(aln.result$level,
                               levels=c("V1","V2","V3","V4","V5","V6","V7","Count"),
                               labels= c("C1","C1","C2","C2","C1","C2","C2","Count"))
  }



  aln.result<-aln.result[,c( "readName","flag","cigar" , "chr","start","end","coord" ,"strand",
     "seq","softClipFragment","trimseq" , "unmapreads","reference_sequence" ,
     "ref_match_trimmed" ,"ref_match_unmapped","is_Arich", "A_number","diff_a",
      "level", "class")]
  colnames(aln.result) <-c( "readName","flag","cigar" , "chr","start","end","coord" ,"strand",
                              "seq","softClipFragment","trimmed_seq" , "unmapped_seq","reference_seq" ,
                              "ref_match_trimmed" ,"ref_match_unmapped"  ,"is_Arich", "A_number","diff_a",
                              "level", "class")

  if(run.quantify){
    print("Completed: Identification and quantification of polyA sites")
    result.list <- resut.PA(aln.result=aln.result,d=d)
    return(result.list)
  }else{
    print("Completed: Identification of polyA sites")
    result.list <- list(pa.table=aln.result)
  }
  return(result.list)

}




#' A process_data function
#'
#' This function processes a BAM file to extract alignment information and
#' polyA tail details.
#'
#' @name process_data
#' @param aln A alignment result from bam file.
#' @param bsgenome A BSgenome object for reference genome. e.g. 'BSgenome.Hsapiens.UCSC.hg38'.
#' @param adjust.chr  Set to TRUE to add 'chr' prefix if input chromosome names
#'                    in BAM file lack 'chr' prefix but the genome sequence includes it.
#'                    Default is FALSE.
#' @param poly Specifies the type of polyA tail at the 3' end of the sequence.
#'             Must be either 'T' (for a poly(T) stretch, e.g., TTTTTT-like tails)
#'             or 'A' (for a poly(A) stretch, e.g., AAAAAA-like tails).
#'             Default is 'T'.
#' @return  A data frame containing alignment results, including chromosome, start, end,
#'         strand, trimmed read sequence, polyA tail sequence, and reference sequence
#'         around potential polyA sites.
#' @export

process_data <- function(aln=NULL,bsgenome=NULL,adjust.chr=FALSE,poly=NULL){
  trimmed.flag <- grep("\\_",  aln$qname)
  if(length(  trimmed.flag)==0){
    aln$qname<- paste0(aln$qname,"_")
  }

  index.temp <- which(aln$flag == '256')
  aln$flag[index.temp] <- '0'
  index.temp <- which(aln$flag == '272')
  aln$flag[index.temp] <- '16'

  colnames(aln)[c(1,3,4)] <- c('readName','chr','start')

  if(adjust.chr){
    aln$chr <- paste0("chr",aln$chr)
  }

  aln <- aln[which(aln$chr %in% as.character(seqnames(bsgenome))),]

  aln$seq_width <-as.numeric(width(aln$seq))

  aln$start <- as.numeric(aln$start)

  if(poly=="T"){
    aln$strand <- ifelse(aln$flag == '0','-','+')
    aln$trimseq <- as.character(reverseComplement(DNAStringSet(gsub("\\S+\\_","",aln$readName))))
  }
  if(poly=="A"){
    aln$flag.raw <- aln$flag
    aln$flag  <- ifelse(aln$flag.raw == '0','16','0')
    aln$strand <- ifelse(aln$flag == '0','-','+')
    aln$trimseq <- as.character((DNAStringSet(gsub("\\S+\\_","",aln$readName))))
  }

  aln$width <- rowSums(apply(str_extract_all(aln$cigar,'\\d+(?=([MND]))',simplify = T),2,as.numeric),na.rm = T)

  aln$end <- aln$start+aln$width-1
  aln <- aln %>% dplyr::mutate(coord=ifelse(flag=='0',start,end))

  aln$sofClipLen.total <- rowSums(apply(str_extract_all(aln$cigar,'\\d+(?=S)',simplify = T),2,as.numeric),na.rm = T)
  aln$softClipLen <- ifelse(aln$flag == '0',as.numeric(str_extract(aln$cigar,"^\\d+(?=S)")),as.numeric(str_extract(aln$cigar,"\\d+(?=S$)")))
  aln$softClipLen[is.na(aln$softClipLen)] <-  0
  aln$softClipLen[is.na(aln$sofClipLen.total)] <-  0


  aln$softClipFragment <- as.character(reverseComplement(DNAStringSet(str_sub(aln$seq, start = 1, end = aln$softClipLen))))
  index <- which(aln$flag=="16")
  if(length( index)>0){
    aln$softClipFragment[index] <-str_sub(aln$seq[index],start=  aln$seq_width[index]-aln$softClipLen[index]+1, end=aln$seq_width[index])
  }
  aln$unmapreads <- paste0(aln$softClipFragment,aln$trimseq,sep="")

  aln$readName <-  as.character(limma::strsplit2(aln$readName,"\\_")[,1])
  aln$unmapreads_width <-as.numeric(width(aln$unmapreads))

  ###Getting reference sequence information
  print("Extract reference sequence around polyA coordinate.")
  aln <- subset(aln, chr %in% as.character(seqnames(bsgenome)))
  aln <- get_referenceseq(aln=aln,bsgenome=bsgenome)
  aln$ref_match_unmapped <- str_sub(aln$reference_sequence,start = 1,end = aln$unmapreads_width)
  aln$ref_match_trimmed <- str_sub(aln$reference_sequence,start = aln$softClipLen +1 ,end = aln$unmapreads_width)
  aln$reference_rest <- str_sub(aln$reference_sequence,start = aln$unmapreads_width+1,end = aln$unmapreads_width+10)
  aln$ref_rest_a <- apply(nchar(str_extract_all(aln$reference_rest,pattern = '^A+',simplify = T)),MARGIN = 1,max)
  aln$coord_raw <- aln$coord
  rownames(aln) <- NULL
  return(aln)
}



#' A process_tail_pure function
#'
#' This function to identify polyA sites from reads with perfectly matched polyA tails.
#'
#' @name process_tail_pure
#' @param aln.tail.pure A alignment result about reads with perfectly matched polyA tails.
#' @return  A data frame of identified polyA sites.
#' @export

process_tail_pure <- function(aln.tail.pure=NULL){
  aln.tail.pure$ref_a <- apply(nchar(str_extract_all(aln.tail.pure$ref_match_unmapped,pattern = '^A+',simplify = T)),MARGIN = 1,max)
  index.temp <-which(is.infinite(aln.tail.pure$ref_a))
  aln.tail.pure$ref_a[index.temp] <- 0

  aln.tail.pure$diff_a <- aln.tail.pure$unmapreads_width -  aln.tail.pure$ref_a
  print("Detecting V1 polyA sites")
  aln.tail.pure$priority <- "M1_V1"
  aln.tail.pure$priority[which(aln.tail.pure$diff_a<=1)] <- "M1_V2"

  #-----------Adjust polyA coordiante-------------------------------------------
  index.plus <- which(aln.tail.pure$flag=="0")
  index.minus <- which(aln.tail.pure$flag=="16")
  aln.tail.pure$coord[index.plus] <-  aln.tail.pure$coord[index.plus] -  aln.tail.pure$ref_a[index.plus]
  aln.tail.pure$start[index.plus] <-  aln.tail.pure$coord[index.plus]
  aln.tail.pure$coord[index.minus] <-  aln.tail.pure$coord[index.minus] +  aln.tail.pure$ref_a[index.minus]
  aln.tail.pure$end[index.minus] <-  aln.tail.pure$coord[index.minus]
  return(aln.tail.pure)
}


#' process_tail_nopure_v4
#'
#' This function to exclude reads with perfectly matched polyA tails, but polyA tails
#' actually from reference genome
#'
#' @name process_tail_pure
#' @param aln.tail.nopure.v4  A alignment result about reads with perfectly matched polyA tails.
#' @return  A data frame of alignment result.
#' @export


process_tail_nopure_v4<-function(aln.tail.nopure.v4=NULL){
  index.plus <- which(aln.tail.nopure.v4$flag=="0")
  index.minus <- which(aln.tail.nopure.v4$flag=="16")
  aln.tail.nopure.v4$coord[index.plus] <-  aln.tail.nopure.v4$coord[index.plus] -  aln.tail.nopure.v4$ref_a[index.plus]
  aln.tail.nopure.v4$start[index.plus] <-  aln.tail.nopure.v4$start[index.plus] -  aln.tail.nopure.v4$ref_a[index.plus]
  aln.tail.nopure.v4$coord[index.minus] <-  aln.tail.nopure.v4$coord[index.minus] +  aln.tail.nopure.v4$ref_a[index.minus]
  aln.tail.nopure.v4$end[index.minus] <-  aln.tail.nopure.v4$end[index.minus] +  aln.tail.nopure.v4$ref_a[index.minus]
  aln.tail.nopure.v4$priority <- "M2_C4_IP"
  return(aln.tail.nopure.v4)
}




#' A process_tail_nopure_v3  function
#'
#' This function to detect polyA sites from reads with partially matched polyA tails.
#'
#' @name process_tail_nopure_v3
#' @param aln.result A alignment result from bam file.
#' @param aln.tail.nopure.v3 A alignment result of reads with partially matched polyA tails.
#' @return  A data frame containing more detailed alignment results.
#' @export

process_tail_nopure_v3 <- function( aln.tail.nopure.v3=NULL,aln.result=NULL){
  aln.tail.nopure.v3$trimseq_end_a <- apply(nchar(str_extract_all(aln.tail.nopure.v3$trimseq,pattern = 'A+$',simplify = T)),MARGIN = 1,max)
  aln.tail.nopure.v3$trimseq_upper <- str_sub(aln.tail.nopure.v3$trimseq,start = 1,end = as.numeric(width(aln.tail.nopure.v3$trimseq))-aln.tail.nopure.v3$trimseq_end_a)
  aln.tail.nopure.v3$ref_match_trimmed_upper <- str_sub(aln.tail.nopure.v3$ref_match_trimmed,start = 1,end = as.numeric(width(aln.tail.nopure.v3$trimseq))-aln.tail.nopure.v3$trimseq_end_a)
  aln.tail.nopure.v3$is.subset.genome <- (aln.tail.nopure.v3$trimseq_upper == aln.tail.nopure.v3$ref_match_trimmed_upper)
  aln.tail.nopure.v3$ref_match_trimmed_down <- str_sub( aln.tail.nopure.v3$ref_match_trimmed,
                                                        start = as.numeric(width( aln.tail.nopure.v3$trimseq))- aln.tail.nopure.v3$trimseq_end_a+1,
                                                        end =as.numeric(width( aln.tail.nopure.v3$trimseq)) )
  aln.tail.nopure.v3$refseq_a <- apply(nchar(str_extract_all( aln.tail.nopure.v3$ref_match_trimmed_down,pattern = '^A+',simplify = T)),MARGIN = 1,max)
  aln.tail.nopure.v3$diff_a <-  aln.tail.nopure.v3$trimseq_end_a -   aln.tail.nopure.v3$refseq_a


  index.keep <- which(aln.tail.nopure.v3$diff_a>0&aln.tail.nopure.v3$is.subset.genome=="TRUE")
  index.plus <- intersect(which(aln.tail.nopure.v3$flag=="0"),index.keep)
  index.minus <- intersect(which(aln.tail.nopure.v3$flag=="16"),index.keep)
  aln.tail.nopure.v3$coord[index.plus] <-  aln.tail.nopure.v3$coord[index.plus]-
    as.numeric(width(aln.tail.nopure.v3$ref_match_trimmed_upper[index.plus]))-aln.tail.nopure.v3$refseq_a[index.plus]
  aln.tail.nopure.v3$start[index.plus] <- aln.tail.nopure.v3$coord[index.plus]

  aln.tail.nopure.v3$coord[index.minus] <-  aln.tail.nopure.v3$coord[index.minus] +
    as.numeric(width(aln.tail.nopure.v3$ref_match_trimmed_upper[index.minus])) +aln.tail.nopure.v3$refseq_a[index.minus]

  aln.tail.nopure.v3$end[index.minus] <-   aln.tail.nopure.v3$coord[index.minus]
  aln.tail.nopure.v3.sub <- aln.tail.nopure.v3[index.keep,]
  aln.tail.nopure.v3.other <- aln.tail.nopure.v3[-index.keep,]


  if(nrow(aln.tail.nopure.v3.sub)>0){
    print("Processing reads without soft-clippingand and with partially matched polyA tails")
    print("Detecting V2 polyA sites")
    aln.tail.nopure.v3.sub$priority <- "M2_C3_V1"
    aln.tail.nopure.v3.sub$priority[which(aln.tail.nopure.v3.sub$diff_a==1)] <- "M2_C3_V2"
    aln.tail.nopure.v3.sub$priority[which(aln.tail.nopure.v3.sub$diff_a==0)] <- "M2_C3_V3"
    aln.result <- plyr::rbind.fill(aln.result,aln.tail.nopure.v3.sub)
    rm(aln.tail.nopure.v3.sub)
  }

  index.other<- which(aln.tail.nopure.v3.other$diff_a==0&aln.tail.nopure.v3.other$is.subset.genome=="TRUE")
  aln.tail.nopure.v3.sub <- aln.tail.nopure.v3.other[index.other,]
  if(length(index.other)==0){
    aln.tail.nopure.v3.other <- aln.tail.nopure.v3.other
  }else{
    aln.tail.nopure.v3.other <- aln.tail.nopure.v3.other[-index.other,]
  }

  if(nrow(aln.tail.nopure.v3.sub)>0){
    print("Detecting part of polyA tail from reference genome")
    index.plus <- which(aln.tail.nopure.v3.sub$flag=="0")
    index.minus <- which(aln.tail.nopure.v3.sub$flag=="16")
    aln.tail.nopure.v3.sub$coord[index.plus] <-  aln.tail.nopure.v3.sub$coord[index.plus]-aln.tail.nopure.v3.sub$unmapreads_width[index.plus]
    aln.tail.nopure.v3.sub$start[index.plus] <- aln.tail.nopure.v3.sub$coord[index.plus]

    aln.tail.nopure.v3.sub$coord[index.minus] <-  aln.tail.nopure.v3.sub$coord[index.minus] +aln.tail.nopure.v3.sub$unmapreads_width[index.minus]
    aln.tail.nopure.v3.sub$end[index.minus] <-   aln.tail.nopure.v3.sub$coord[index.minus]

    print("Detecting V7 polyA sites")
    aln.tail.nopure.v3.sub$priority <- "M2_C3_NoTail"
    index <- which(aln.tail.nopure.v3.sub$ref_rest_a>5)
    aln.tail.nopure.v3.sub$priority[index] <- "M2_C3_IP"
    aln.result <- plyr::rbind.fill(aln.result,aln.tail.nopure.v3.sub)
    rm(aln.tail.nopure.v3.sub)
  }

  if(nrow(aln.tail.nopure.v3.other)>0){
    print("Detecting V3 and V4 polyA sites")
    aln.tail.nopure.v3.other$trimseq_start_a <- apply(nchar(str_extract_all(aln.tail.nopure.v3.other$trimseq,pattern = '^A+',simplify = T)),MARGIN = 1,max)
    aln.tail.nopure.v3.other$red_start_a <- apply(nchar(str_extract_all(aln.tail.nopure.v3.other$reference_sequence,pattern = '^A+',simplify = T)),MARGIN = 1,max)

    aln.tail.nopure.v3.other$bais_a <-aln.tail.nopure.v3.other$trimseq_start_a - aln.tail.nopure.v3.other$red_start_a
    aln.tail.nopure.v3.other <-  aln.tail.nopure.v3.other%>%mutate(bais_aa=ifelse(bais_a>0,red_start_a,trimseq_start_a))

    index.plus <- which(aln.tail.nopure.v3.other$flag=="0")
    index.minus <- which(aln.tail.nopure.v3.other$flag=="16")
    aln.tail.nopure.v3.other$coord[index.plus] <-  aln.tail.nopure.v3.other$coord[index.plus]-aln.tail.nopure.v3.other$bais_aa[index.plus]
    aln.tail.nopure.v3.other$start[index.plus] <- aln.tail.nopure.v3.other$coord[index.plus]

    aln.tail.nopure.v3.other$coord[index.minus] <-  aln.tail.nopure.v3.other$coord[index.minus] +aln.tail.nopure.v3.other$bais_aa[index.minus]
    aln.tail.nopure.v3.other$end[index.minus] <-   aln.tail.nopure.v3.other$coord[index.minus]
    aln.tail.nopure.v3.other$priority <- "M2_C3_V5"
    aln.tail.nopure.v3.other$priority[which( aln.tail.nopure.v3.other$trimseq_end_a>1)] <- "M2_C3_V4"
    aln.tail.nopure.v3.other$priority[which( aln.tail.nopure.v3.other$trimseq_end_a>5)] <- "M2_C3_V3"

    aln.result <- plyr::rbind.fill(aln.result,aln.tail.nopure.v3.other)
    rm(aln.tail.nopure.v3.other)
  }
  return(aln.result)
}


#' A process_tail_nopure_v1  function
#'
#' This function to detect polyA sites from reads with perfectly matched polyA tails and soft-clipping.
#'
#' @name process_tail_nopure_v1
#' @param aln.result A alignment result from bam file.
#' @param aln.tail.nopure.v1 A alignment result of reads with perfectly matched polyA tails and soft-clipping.
#' @return  A data frame containing more detailed alignment results.
#' @export
#'

process_tail_nopure_v1 <- function(aln.tail.nopure.v1=NULL,aln.result=NULL){
  aln.tail.nopure.v1$trimseq_a <- apply(nchar(str_extract_all(aln.tail.nopure.v1$trimseq,pattern = 'A+$',simplify = T)),MARGIN = 1,max)
  aln.tail.nopure.v1$reference_sequence_1 <- str_sub(aln.tail.nopure.v1$reference_sequence,start=aln.tail.nopure.v1$softClipLen+1,end=aln.tail.nopure.v1$unmapreads_width+10)
  aln.tail.nopure.v1$ref_a <- apply(nchar(str_extract_all(aln.tail.nopure.v1$reference_sequence_1,pattern = '^A+',simplify = T)),MARGIN = 1,max)

  index.in <- which(is.infinite(aln.tail.nopure.v1$ref_a))
  aln.tail.nopure.v1$ref_a[ index.in]<-0

  aln.tail.nopure.v1$diff_a <- aln.tail.nopure.v1$trimseq_a  - aln.tail.nopure.v1$ref_a

  index.1 <- which(aln.tail.nopure.v1$diff_a>1)
  index.2 <- which(aln.tail.nopure.v1$diff_a==1)
  index.3 <- which(aln.tail.nopure.v1$softClipLen<=5)
  index.4 <- which(aln.tail.nopure.v1$softClipLen>5)

  index13 <- intersect(index.1,index.3)
  index23 <- intersect(index.2,index.3)
  index14 <- intersect(index.1,index.4)

  index1234 <- union(union(index13, index23),index14)

  aln.tail.nopure.v1$priority[index13] <- "M2_C1_V1"

  aln.tail.nopure.v1$priority[index23] <- "M2_C1_V2"
  aln.tail.nopure.v1$priority[index14] <- "M2_C1_V3"
  aln.tail.nopure.v1$priority[-index1234] <- "M2_C1_Other"
  index.plus <- which(aln.tail.nopure.v1$flag=="0")
  index.minus <- which(aln.tail.nopure.v1$flag=="16")
  aln.tail.nopure.v1$coord[index.plus] <-  aln.tail.nopure.v1$coord[index.plus] -  aln.tail.nopure.v1$ref_a[index.plus] - aln.tail.nopure.v1$softClipLen[index.plus]
  aln.tail.nopure.v1$start[index.plus] <-   aln.tail.nopure.v1$coord[index.plus]
  aln.tail.nopure.v1$coord[index.minus] <-  aln.tail.nopure.v1$coord[index.minus] +  aln.tail.nopure.v1$ref_a[index.minus]+aln.tail.nopure.v1$softClipLen[index.minus]
  aln.tail.nopure.v1$end[index.minus] <-   aln.tail.nopure.v1$coord[index.minus]


  aln.result <- plyr::rbind.fill(aln.result,aln.tail.nopure.v1)
  rm(aln.tail.nopure.v1)
  return(aln.result )
}



#' process_tail_nopure_v2  function
#'
#' This function to detect polyA sites from reads with perfectly matched polyA tails and soft-clipping.
#'
#' @name process_tail_nopure_v2
#' @param aln.result A alignment result from bam file.
#' @param aln.tail.nopure.v2 A alignment result of reads with perfectly matched polyA tails and soft-clipping.
#' @return  A data frame containing more detailed alignment results.
#' @export
#'

process_tail_nopure_v2 <- function(aln.tail.nopure.v2=NULL, aln.result=NULL  ){
  index <- which(aln.tail.nopure.v2$trimseq=="")
  aln.tail.nopure.v2$priority <- ""
  aln.tail.nopure.v2$priority[index] <-  "M2_C2_NoTail"
  aln.result <- plyr::rbind.fill(aln.result,aln.tail.nopure.v2[index,])
  aln.tail.nopure.v2 <- aln.tail.nopure.v2[-index,]
  if(nrow(aln.tail.nopure.v2)>0){
    aln.tail.nopure.v2$is.subset.genome <- apply(aln.tail.nopure.v2, 1, function(row) grepl(row["trimseq"], row["reference_sequence"]))
    aln.tail.nopure.v2.true <- subset(aln.tail.nopure.v2,is.subset.genome=="TRUE")
    aln.tail.nopure.v2.false <-subset(aln.tail.nopure.v2,is.subset.genome=="FALSE")

    #for true, trimmed sequence from reference genome
    #for true, trimmed sequence from reference genome
    if( nrow(aln.tail.nopure.v2.true)>0){
      aln.tail.nopure.v2.true$priority <- "M2_C2_IP"
      index.plus <- which(aln.tail.nopure.v2.true$flag=="0")
      index.minus <- which(aln.tail.nopure.v2.true$flag=="16")
      aln.tail.nopure.v2.true$coord[index.plus] <-  aln.tail.nopure.v2.true$coord[index.plus]  - aln.tail.nopure.v2.true$unmapreads_width[index.plus]
      aln.tail.nopure.v2.true$start[index.plus] <-  aln.tail.nopure.v2.true$start[index.plus] - aln.tail.nopure.v2.true$unmapreads_width[index.plus]
      aln.tail.nopure.v2.true$coord[index.minus] <-  aln.tail.nopure.v2.true$coord[index.minus] +aln.tail.nopure.v2.true$unmapreads_width[index.minus]
      aln.tail.nopure.v2.true$end[index.minus] <-  aln.tail.nopure.v2.true$end[index.minus] +  aln.tail.nopure.v2.true$unmapreads_width[index.minus]
      aln.result <- plyr::rbind.fill(aln.result,aln.tail.nopure.v2.true)
    }


    if(nrow(aln.tail.nopure.v2.false)>0){
      aln.tail.nopure.v2.false$trimseq_end_a <- apply(nchar(str_extract_all(aln.tail.nopure.v2.false$trimseq,pattern = 'A+$',simplify = T)),MARGIN = 1,max)
      aln.tail.nopure.v2.false$trimseq_upper <- str_sub(aln.tail.nopure.v2.false$trimseq,start = 1,end = as.numeric(width(aln.tail.nopure.v2.false$trimseq))-aln.tail.nopure.v2.false$trimseq_end_a)
      aln.tail.nopure.v2.false$ref_match_trimmed_upper <- str_sub(aln.tail.nopure.v2.false$ref_match_trimmed,start = 1,end = as.numeric(width(aln.tail.nopure.v2.false$trimseq))-aln.tail.nopure.v2.false$trimseq_end_a)
      aln.tail.nopure.v2.false$is.subset.genome<- (aln.tail.nopure.v2.false$trimseq_upper==aln.tail.nopure.v2.false$ref_match_trimmed_upper)
      aln.tail.nopure.v2.true <- subset(aln.tail.nopure.v2.false,is.subset.genome=="TRUE")

      if(!isEmpty(aln.tail.nopure.v2.true)){
        aln.tail.nopure.v2.true$ref_match_trimmed_down <- str_sub(aln.tail.nopure.v2.true$ref_match_trimmed,
                                                                  start = as.numeric(width(aln.tail.nopure.v2.true$trimseq))-aln.tail.nopure.v2.true$trimseq_end_a+1,
                                                                  end =as.numeric(width(aln.tail.nopure.v2.true$trimseq)) )
        aln.tail.nopure.v2.true$refseq_a <- apply(nchar(str_extract_all(aln.tail.nopure.v2.true$ref_match_trimmed_down,pattern = '^A+',simplify = T)),MARGIN = 1,max)
        index.in <- which(is.infinite(aln.tail.nopure.v2.true$refseq_a))
        aln.tail.nopure.v2.true$refseq_a[ index.in]<-0

        aln.tail.nopure.v2.true$diff_a <- aln.tail.nopure.v2.true$trimseq_end_a -  aln.tail.nopure.v2.true$refseq_a

        index.plus <- which(aln.tail.nopure.v2.true$flag=="0")
        index.minus <- which(aln.tail.nopure.v2.true$flag=="16")
        aln.tail.nopure.v2.true$coord[index.plus] <-  aln.tail.nopure.v2.true$coord[index.plus]-
          as.numeric(width(aln.tail.nopure.v2.true$ref_match_trimmed_upper[index.plus]))-
          aln.tail.nopure.v2.true$softClipLen[index.plus] - aln.tail.nopure.v2.true$refseq_a[index.plus]
        aln.tail.nopure.v2.true$start[index.plus] <- aln.tail.nopure.v2.true$coord[index.plus]

        aln.tail.nopure.v2.true$coord[index.minus] <-  aln.tail.nopure.v2.true$coord[index.minus] +
          as.numeric(width(aln.tail.nopure.v2.true$ref_match_trimmed_upper[index.minus])) +
          aln.tail.nopure.v2.true$softClipLen[index.minus]+aln.tail.nopure.v2.true$refseq_a[index.minus]

        aln.tail.nopure.v2.true$end[index.minus] <-   aln.tail.nopure.v2.true$coord[index.minus]

        aln.tail.nopure.v2.true$priority <- "M2_C2_Other"
        index <- which(aln.tail.nopure.v2.true$diff_a>1&aln.tail.nopure.v2.true$softClipLen<=5)
        aln.tail.nopure.v2.true$priority[index] <- "M2_C2_V6"
        index <- which(aln.tail.nopure.v2.true$diff_a>6&aln.tail.nopure.v2.true$softClipLen>5)
        aln.tail.nopure.v2.true$priority[index] <- "M2_C2_V7"
        aln.result <- plyr::rbind.fill(aln.result,aln.tail.nopure.v2.true)
        rm(aln.tail.nopure.v2.true)
      }


      aln.tail.nopure.v2.false <- subset(aln.tail.nopure.v2.false,is.subset.genome=="FALSE")

      if(!isEmpty(aln.tail.nopure.v2.false)){
        aln.tail.nopure.v2.false$priority <- "M2_C2_V5"
        index <- which(aln.tail.nopure.v2.false$trimseq_end_a>7 & aln.tail.nopure.v2.false$softClipLen<=5)

        aln.tail.nopure.v2.false$priority[index] <- "M2_C2_V3"
        index1 <- which( aln.tail.nopure.v2.false$softClipLen<=5)
        index2<-  setdiff(index1,index)
        aln.tail.nopure.v2.false$priority[index2] <- "M2_C2_V4"
        aln.result <- plyr::rbind.fill(aln.result,aln.tail.nopure.v2.false)
        rm(aln.tail.nopure.v2.false)
      }
    }



  }
  return(aln.result )
}


#' A correctSNP  function
#'
#' This function to detect potential polyA sites affected by SNPs.
#'
#' @name correctSNP
#' @param reads A data frame containing alignment results.
#' @return  A data frame containing corrected polyA sites affected by SNPs
#' @export
#'

correctSNP <- function(reads=NULL){

  reads$ID <- paste0("PAS",1:nrow(reads))
  reads.umapped <- subset(reads,softClipLen>1)
  reads.umapped$unsoft_cp <- gsub("A+$","",reads.umapped$softClipFragment)
  reads.umapped$unsoft_cp_width <- width(reads.umapped$unsoft_cp)
  reads.umapped <- subset(reads.umapped,unsoft_cp_width>1)
  reads.umapped$ref_softFragment <-  str_sub(reads.umapped$ref_match_unmapped,1,reads.umapped$unsoft_cp_width)
  reads.umapped$unsoft_cp_A <- apply(nchar(str_extract_all(reads.umapped$unsoft_cp,pattern = 'A',simplify = T)),MARGIN = 1,sum)
  reads.umapped$unsoft_cp_C <- apply(nchar(str_extract_all(reads.umapped$unsoft_cp,pattern = 'C',simplify = T)),MARGIN = 1,sum)
  reads.umapped$unsoft_cp_T <- apply(nchar(str_extract_all(reads.umapped$unsoft_cp,pattern = 'T',simplify = T)),MARGIN = 1,sum)
  reads.umapped$unsoft_cp_G <- apply(nchar(str_extract_all(reads.umapped$unsoft_cp,pattern = 'G',simplify = T)),MARGIN = 1,sum)


  reads.umapped$ref_cp_A <- apply(nchar(str_extract_all(reads.umapped$ref_softFragment,pattern = 'A',simplify = T)),MARGIN = 1,sum)
  reads.umapped$ref_cp_C <- apply(nchar(str_extract_all(reads.umapped$ref_softFragment,pattern = 'C',simplify = T)),MARGIN = 1,sum)
  reads.umapped$ref_cp_T <- apply(nchar(str_extract_all(reads.umapped$ref_softFragment,pattern = 'T',simplify = T)),MARGIN = 1,sum)
  reads.umapped$ref_cp_G <- apply(nchar(str_extract_all(reads.umapped$ref_softFragment,pattern = 'G',simplify = T)),MARGIN = 1,sum)


  reads.umapped$diff_ref_unsoft <- abs(reads.umapped$unsoft_cp_A-reads.umapped$ref_cp_A )+
    abs(reads.umapped$unsoft_cp_C-reads.umapped$ref_cp_C )+
    abs(reads.umapped$unsoft_cp_T-reads.umapped$ref_cp_T )+
    abs(reads.umapped$unsoft_cp_G-reads.umapped$ref_cp_G )


  reads.umapped.snp <- subset(reads.umapped, diff_ref_unsoft<=3)
  reads.umapped.snp <- reads.umapped.snp[reads.umapped.snp$softClipLen>reads.umapped.snp$diff_ref_unsoft,]



  reads.umapped.snp$soft_end_a <-  apply(nchar(str_extract_all(reads.umapped.snp$softClipFragment,pattern = 'A+$',simplify = T)),MARGIN = 1,max)

  if( length(unique(reads.umapped.snp$soft_end_a))==1&&is.infinite(unique(reads.umapped.snp$soft_end_a))){
    reads.umapped.snp$soft_end_a <-0
  }

  reads.umapped.snp$unmmaped_end_a_seq <-  str_sub(reads.umapped.snp$unmapreads,1+reads.umapped.snp$unsoft_cp_width,reads.umapped.snp$unmapreads_width)


  reads.umapped.snp$ref_unmmaped_end_a_seq <-  str_sub(reads.umapped.snp$ref_match_unmapped,1+reads.umapped.snp$unsoft_cp_width,reads.umapped.snp$unmapreads_width)
  reads.umapped.snp$ref_unmmaped_end_a_seq <- paste0(reads.umapped.snp$ref_unmmaped_end_a_seq,reads.umapped.snp$reference_rest)


  reads.umapped.snp$unmmaped_end_a_seq_top <- apply(nchar(str_extract_all(reads.umapped.snp$unmmaped_end_a_seq,pattern = '^A+',simplify = T)),MARGIN = 1,sum)

  if( length(unique(reads.umapped.snp$unmmaped_end_a_seq_top))==1&&is.infinite(unique(reads.umapped.snp$unmmaped_end_a_seq_top))){
    reads.umapped.snp$unmmaped_end_a_seq_top <-0
  }

  reads.umapped.snp$ref_unmmaped_end_a_seq_top <- apply(nchar(str_extract_all(reads.umapped.snp$ref_unmmaped_end_a_seq,pattern = '^A+',simplify = T)),MARGIN = 1,sum)

  if( length(unique(reads.umapped.snp$ref_unmmaped_end_a_seq_top))==1&&is.infinite(unique(reads.umapped.snp$ref_unmmaped_end_a_seq_top))){
    reads.umapped.snp$ref_unmmaped_end_a_seq_top <-0
  }

  reads.umapped.snp$unmmaped_end_a_seq_top_diff <- reads.umapped.snp$unmmaped_end_a_seq_top -reads.umapped.snp$ref_unmmaped_end_a_seq_top

  reads.umapped.snp.need <- subset(reads.umapped.snp,unmmaped_end_a_seq_top_diff>=2)



  reads.umapped.snp.need$temp_width <- width(reads.umapped.snp.need$unmmaped_end_a_seq)
  reads.umapped.snp.need$temp_width_diff <- reads.umapped.snp.need$temp_width-reads.umapped.snp.need$unmmaped_end_a_seq_top
  reads.umapped.snp.need <- subset(reads.umapped.snp.need,temp_width_diff==0)
  index.plus <- which(reads.umapped.snp.need$flag=="0")
  index.minus <- which(reads.umapped.snp.need$flag=="16")
  reads.umapped.snp.need$coord[index.plus] <-  reads.umapped.snp.need$coord_raw[index.plus] -  reads.umapped.snp.need$ref_unmmaped_end_a_seq_top[index.plus]-reads.umapped.snp.need$softClipLen[index.plus]
  reads.umapped.snp.need$coord[index.minus] <-  reads.umapped.snp.need$coord_raw[index.minus] +  reads.umapped.snp.need$ref_unmmaped_end_a_seq_top[index.minus]+reads.umapped.snp.need$softClipLen[index.minus]

  reads.umapped.snp.need$level_n <- as.character(reads.umapped.snp.need$level)
  reads.umapped.snp.need$level_n  <- paste0("recycle_",reads.umapped.snp.need$level_n)
  reads$level <- as.character(reads$level)
  reads$level_n <- reads$level

  index <- which(reads$ID %in% reads.umapped.snp.need$ID)

  reads.output <- plyr::rbind.fill(reads[-index,],reads.umapped.snp.need)

  reads.output$level_raw <- reads.output$level
  reads.output$level_correct <- reads.output$level_n
  reads.output$level_n <- NULL

  reads.output$level <- reads.output$level_correct


  reads.output$level <- factor(as.character(reads.output$level) ,
                               levels=c("V1","V2","V3","V4","recycle_V5","recycle_Count","V5","V6","V7","Count"),
                               labels=c("V1","V2","V3","V4","V4", "V4","V5","V6","V7","Count"))


  reads.output$level <- as.character(reads.output$level)
  return( reads.output )
}



#' A Arich function
#'
#' This function to search A-rich fragments surrounding polyA sites.
#'
#' @name Arich
#' @param reads A data frame containing alignment results and identified polyA sites.
#' @param bsgenome A BSgenome object for reference genome. e.g. 'BSgenome.Hsapiens.UCSC.hg38'.
#' @return  A data frame containing information about polyA sites near A-rich regions
#' @export
#'

Arich <- function(reads=NULL,bsgenome=NULL){
  reads$reference_sequence <- ""
  index.plus <- which(reads$strand=="-")
  index.minus <- which(reads$strand=="+")

  sub.info.plus <- data.frame(chrom=reads$chr[index.plus],
                              start=as.numeric(reads$coord[index.plus]-10),
                              end=as.numeric(reads$coord[index.plus])+10)
  reads$reference_sequence[index.plus] <- as.character(reverseComplement(getSeq(bsgenome, as(sub.info.plus, "GRanges"))) )


  sub.info.minus <- data.frame(chrom=reads$chr[index.minus],
                               start=as.numeric(reads$coord[index.minus])-10,
                               end=as.numeric(reads$coord[index.minus])+10)
  reads$reference_sequence[index.minus] <- as.character(getSeq(bsgenome, as(sub.info.minus, "GRanges")))


  index<- grep('AAAAAA',reads$reference_sequence)
  reads$is_Arich <- "No"
  reads$is_Arich[index] <- "Yes"
  reads$A_number <- stringr::str_count(reads$reference_sequence, pattern = "A")
  reads$A_number_p <-reads$A_number/nchar(reads$reference_sequence[1])
  index2<- which(reads$A_number_p>0.7)
  length(intersect(index,index2))
  length(setdiff(index2,index))
  reads$is_Arich[index2] <- "Yes"
  reads$unmapreads_A_count <- apply(nchar(str_extract_all(reads$unmapreads,pattern = 'A',simplify = T)),MARGIN = 1,sum)
  reads$unmapreads_nonA_count <- reads$unmapreads_width - reads$unmapreads_A_count
  return(reads)
}


#' A recycle_v8 function
#'
#' This function to detect V8 polyA sites.
#'
#' @name recycle_v8
#' @param aln.result A data frame contains alignment results and identified polyA sites.
#' @param threeUTRregion A GRanges object for 3'UTR region, will use level-8 reads in identification of polyA sites.
#'                       Recommended to use if input reads lack polyA tails.
#' @param cutoffCount Minimum count threshold for PACs identified at the V8 level
#'                    to be used for polyA determination.
#' @param ext3UTRlen An integer denoting the extended length of 3'UTR region. Default is 1000 nt.
#' @return Identifed V8 polyA sites.
#' @export
#'

recycle_v8<- function(aln.result=NULL,cutoffCount=NULL,
                      threeUTRregion=NULL,ext3UTRlen=1000){

  plantdpd.data.po<- subset(aln.result,level=="Count")
  plantdpd.data.level <- subset(aln.result,level!="Count")
  plantdpd.data.po.region <- GRanges(seqnames=plantdpd.data.po$chr,
                                     ranges=IRanges(start= plantdpd.data.po$coord,
                                                    end=plantdpd.data.po$coord),
                                     strand=plantdpd.data.po$strand)
  plantdpd.data.po.region.24 <- GenomicRanges::reduce(plantdpd.data.po.region ,
                                                      min.gapwidth=24)

  ov <- GenomicRanges::findOverlaps(  plantdpd.data.po.region.24,
                                      plantdpd.data.po.region)
  ov <- as.data.frame(ov)
  ov <- ov[order(ov$queryHits),]
  ov <- ov[!duplicated(ov$subjectHits), ]
  ov$count <- 1

  ov.combined <- ov  %>% dplyr::group_by(queryHits) %>%
    dplyr::summarise(total.count=sum(count))
  ov$total.count <- ov.combined$total.count[match(ov$queryHits,ov.combined$queryHits)]


  ov.3utr <- GenomicRanges::findOverlaps(plantdpd.data.po.region.24,
                                         threeUTRregion,
                                         maxgap=ext3UTRlen )

  ov$at_3utr <- "No"
  ov$at_3utr[which(ov$queryHits %in% unique(ov.3utr@from)) ] <- "Yes"
  ov.keep <- subset(ov,total.count>=cutoffCount&at_3utr=="Yes")

  plantdpd.data.po$level <- as.character(plantdpd.data.po$level)
  plantdpd.data.po$level[unique(ov.keep$subjectHits)] <- "V8"
  plantdpd.data.level$level <- as.character(plantdpd.data.level$level)

  aln.result <- plyr::rbind.fill(plantdpd.data.po,plantdpd.data.level)
  return(aln.result)
}




#' A resut.PA function
#'
#' TThis function to classify PACs and quantify PACs.
#'
#' @name resut.PA
#' @param aln.result A data frame contains alignment results and identified polyA sites.
#' @param d A distance to group nearby polyA sites into a polyA site cluster (PAC). Default is 24 nt.
#' @return A list containing two element. The 'pa.table' element
#'         is a mapped table of 3'reads, including sequence alignment results,
#'         base compositions of the polyA tail, polyA site position, and
#'         the confidence level of reads used to identify polyA sites (from V1 to V8).
#'         the 'pa.coord'  element contains identified and quantified PACs.
#' @export
#'

resut.PA <- function(aln.result=NULL,d=NULL){
  #print("Identification and quantification of polyA")
  poly.coord <- priority.pa(aln.result=aln.result,d=d)
  poly.coord$queryHits <- NULL
  poly.coord$PAid <- paste0("PA",1:nrow(poly.coord))
  poly.region <- GenomicRanges::makeGRangesFromDataFrame(poly.coord)
  all.pa.region <-GenomicRanges::GRanges(seqnames = as.character(aln.result$chr),
                                         ranges = IRanges(start = as.integer(aln.result$coord),
                                                          end = as.integer(aln.result$coord)),
                                         strand = as.character(aln.result$strand))
  ov <- GenomicRanges::findOverlaps(  poly.region , all.pa.region)
  ov <- as.data.frame(ov)
  ov <- ov[order(ov$queryHits),]
  ov <- ov[!duplicated(ov$subjectHits), ]

  aln.result$use.as.count <- 0
  aln.result$use.as.count[ov$subjectHits] <- 1
  aln.result$PAid <- NA
  aln.result$PAid[ov$subjectHits] <- paste0("PA",ov$queryHits)


  ov$level <- aln.result$level[ov$subjectHits]
  ov$count <- 1
  ov$count_pa <- 0
  ov$count_pa[grep("^V",as.character(ov$level))] <- 1

  ov <- ov  %>% dplyr::group_by(queryHits) %>%
    dplyr::summarise(total.count=sum(count),count.with.polyA=sum(count_pa))

  ov<- ov[order(ov$queryHits),]
  poly.coord$total.count <- 0
  poly.coord$count_with.polyA <- 0
  poly.coord$total.count[ov$queryHits] <- ov$total.count
  poly.coord$count_with.polyA[ov$queryHits] <- ov$count.with.polyA
  #print("Completed!!!")
  return(list(pa.table=aln.result,pa.coord=poly.coord))
}





priority.pa <- function(aln.result=NULL,d=24){
  #===========Identify polyA position: highest==================================
  primary.pa <- subset(aln.result,level %in% c("V1","V2"))
  if(nrow(primary.pa)==0){
    stop("Please make sure that the input reads contain sequence of poly A tail information: such as @SRR16867064.13677986_AAAAAAAAAAAA")
  }

  ov.primary <- Identify.pa(primary.pa ,d=d)
  primary.region <- GenomicRanges::makeGRangesFromDataFrame(ov.primary)



  #===========Identify polyA position: Second==================================
  second.pa <- subset(aln.result,level %in% c("V3","V4"))
  second.region <- GenomicRanges::GRanges(seqnames = as.character(second.pa$chr),
                                          ranges = IRanges(start = as.integer(second.pa$coord),
                                                           end = as.integer(second.pa$coord)),
                                          strand = as.character(second.pa$strand))
  ov <- GenomicRanges::findOverlaps( second.region , primary.region)
  second.pa <- second.pa[-ov@from,]
  if(nrow(second.pa )>0){
    ov.second <- Identify.pa(second.pa,d=d)

    ov.primary <- plyr::rbind.fill( ov.primary ,ov.second )
    primary.region <- GenomicRanges::makeGRangesFromDataFrame(ov.primary)
  }


  #===========Identify polyA position: Third==================================
  thrid.pa <- subset(aln.result,level %in% c("V5","V6"))
  thrid.region <- GenomicRanges::GRanges(seqnames = as.character(thrid.pa$chr),
                                         ranges = IRanges(start = as.integer(thrid.pa$coord),
                                                          end = as.integer(thrid.pa$coord)),
                                         strand = as.character(thrid.pa$strand))
  ov <- GenomicRanges::findOverlaps( thrid.region , primary.region)
  thrid.pa <- thrid.pa[-ov@from,]
  if(nrow(thrid.pa )>0){
    ov.thrid <- Identify.pa(thrid.pa,d=d)

    ov.primary <- plyr::rbind.fill( ov.primary ,ov.thrid )
    primary.region <- GenomicRanges::makeGRangesFromDataFrame(ov.primary)
  }


  #===========Identify polyA position: Third==================================
  fourth.pa <- subset(aln.result,level %in% c("V7","V8"))
  fourth.region <- GRanges(seqnames = as.character(fourth.pa$chr),
                           ranges = IRanges(start = as.integer(fourth.pa$coord),
                                            end = as.integer(fourth.pa$coord)),
                           strand = as.character(fourth.pa$strand))
  ov <- GenomicRanges::findOverlaps( fourth.region , primary.region)
  fourth.pa <- fourth.pa[-ov@from,]
  if(nrow(fourth.pa )>0){

    fourth.pa$label <- paste(fourth.pa$chr,fourth.pa$coord,fourth.pa$strand,sep="_")
    temp <- as.data.frame(table( fourth.pa$label))
    index <- match(fourth.pa$label, temp$Var1)
    fourth.pa$diff_a_raw <- fourth.pa$diff_a
    fourth.pa$diff_a <- temp$Freq[index]

    ov.fourth <- Identify.pa(fourth.pa,d=d)

    ov.primary <- plyr::rbind.fill( ov.primary ,ov.fourth )
    primary.region <- GenomicRanges::makeGRangesFromDataFrame(ov.primary)
  }
  return(ov.primary)
}


Identify.pa <- function(pa.data=NULL,d=24){
  pa.data$label <- paste(pa.data$chr,pa.data$coord,pa.data$strand,sep="_")
  pa.data <- pa.data[order(pa.data$label,pa.data$level,-pa.data$diff_a),]
  pa.data <- pa.data[!duplicated(pa.data$label),]


  primary.region <- GenomicRanges::GRanges(seqnames = as.character(pa.data$chr),
                                           ranges = IRanges(start = as.integer(pa.data$coord),
                                                            end = as.integer(pa.data$coord)),
                                           strand = as.character(pa.data$strand))
  primary.region.reduce <- GenomicRanges::reduce(primary.region,min.gapwidth=d)
  ov <- GenomicRanges::findOverlaps( primary.region.reduce,primary.region)
  ov <- as.data.frame(ov)

  ov <- cbind(ov,pa.data[ov$subjectHits,c("coord","readName","diff_a","level")])
  ov <- ov[order(ov$queryHits,ov$level,-ov$diff_a),]

  ov.uq <- ov[!duplicated(ov$queryHits),]
  colnames(ov.uq )[3:6] <- c("coord","coord_readName","coord_unmappedA","coord_level")
  ov <- ov  %>% dplyr::group_by(queryHits) %>%
    dplyr::summarise(UPA_coord=paste(coord,collapse = ","),
                     UPA_readName=paste(readName,collapse = ","),
                     UPA_unmappedA=paste(diff_a,collapse = ","),
                     UPA_level=paste(level,collapse = ","))
  index <- match(ov$queryHits,ov.uq$queryHits)
  ov <- cbind(ov.uq[index,c(1,3:6)],ov[,-1] )
  index <- match(1:length(primary.region.reduce),ov$queryHits)
  ov <- cbind(as.data.frame(primary.region.reduce),ov[index,])
  return(ov)
}




#' A simpleCluster function
#'
#' Cluster poly(A) sites into groups by a certain distance iterative.
#'
#' @name simpleCluster
#' @param points.gr A GRanges object constructed from poly(A) sites.
#' @param max.gapwidth A cutoff to limit the maximum distance between two adjacent sites in a PAC.
#' @return A reduced GRanges object.
#' @importFrom GenomicRanges reduce
#' @importFrom IRanges extractList
#' @importFrom BiocGenerics start
#' @export
#'
simpleCluster <- function(points.gr,max.gapwidth=24){

  # Cluster points by distance
  range.gr = reduce(points.gr,min.gapwidth=max.gapwidth,with.revmap=T,ignore.strand=FALSE)

  # Sum the score
  range.gr$score = sum(IRanges::extractList(points.gr$score,range.gr$revmap))

  # Select the center
  idx = BiocGenerics::which.max(IRanges::extractList(points.gr$score,range.gr$revmap))
  range.gr$center = start(points.gr)[idx+c(0,cumsum(lengths(range.gr$revmap))[1:(length(range.gr$revmap)-1)])]

  # Return result
  return(range.gr)
}



#' A findPeaks function
#'
#' Cluster the one-dimensional poly(A) sites into groups based on local density

#' @name findPeaks
#' @usage findPeaks(sub_pos,sub_wts,min_delta=24)
#' @param sub_pos The genomic coordinates of poly(A) sites.
#' @param sub_wts The numbers of reads support each poly(A) site in 'sub_pos'.
#' @param min_delta A cutoff to limit the minimum distance between two centers.
#' @return A tibble object with five columns 'group', 'start', 'end', 'sum.wts', 'center'.
#' @importFrom dplyr tibble %>% group_by summarise
#' @importFrom outliers scores
#' @export
#'

findPeaks <- function(sub_pos,sub_wts,min_delta=24){
  # Declare
  group = pos = wts = NULL
  ND <- length(sub_pos)

  if(ND<=2){
    res <- tibble::tibble(group=1,start=min(sub_pos),end=max(sub_pos),sum.wts=sum(sub_wts),center=sub_pos[which.max(sub_wts)])
    return(res)
  }

  dc <- 0.5
  #paste('Computing Rho with gaussian kernel of radius: ',dc,collapse = '')

  rho = rep(0,ND)
  #
  # Gaussian kernel
  #
  for(i in 1:(ND-1)){
    for(j in (i+1):ND){
      tmp_dist = abs(sub_pos[i]-sub_pos[j])
      tmp_wts = sub_wts[i]*sub_wts[j]
      rho[i]=rho[i]+exp(-(tmp_dist/dc)^2)*sub_wts[j];
      rho[j]=rho[j]+exp(-(tmp_dist/dc)^2)*sub_wts[i];
    }
  }
  for(i in 1:ND){
    rho[i]=rho[i]+sub_wts[i];
  }

  maxd <- abs(sub_pos[1] - sub_pos[ND])
  rho_sorted = sort(rho,decreasing = T,index.return=T)
  ordrho <- rho_sorted$ix
  rho_sorted <- rho_sorted$x

  delta <- rep(-1.,ND)
  nneigh <- rep(0,ND)

  for(ii in 2:ND){
    delta[ordrho[ii]] = maxd
    for(jj in 1:(ii-1)){
      tmp_dist = abs(sub_pos[ordrho[ii]]-sub_pos[ordrho[jj]])
      if(tmp_dist<=delta[ordrho[ii]]){
        delta[ordrho[ii]] = tmp_dist
        nneigh[ordrho[ii]] = ordrho[jj]
      }
    }
  }
  delta[ordrho[1]]=max(delta)

  decision.data <- data.frame(rho=rho,delta=delta,rho.delta=rho*delta)
  outlier <- rho>max(rho)/2 & delta>max(delta)/2 & delta>min_delta

  tmp <- decision.data$rho.delta
  tmp[which(outlier)] <- 0
  outlier <- (outlier | outliers::scores(tmp, type="chisq", prob=0.99))& delta>min_delta

  if(!sum(outlier)){
    res <- tibble::tibble(group=1,start=min(sub_pos),end=max(sub_pos),sum.wts=sum(sub_wts),center=sub_pos[which.max(sub_wts)])
    return(res)
  }

  NCLUST = 0
  cl = rep(-1,ND)
  icl = c()
  for(i in 1:ND){
    if(outlier[i]){
      NCLUST = NCLUST + 1
      cl[i] = NCLUST
      icl[NCLUST] = i
    }
  }
  #paste('NUMBER OF CLUSTERS: ',NCLUST,collapse = '')

  # Assignation
  for(i in 1:ND){
    if (cl[ordrho[i]]==-1){
      cl[ordrho[i]] = cl[nneigh[ordrho[i]]];
    }
  }

  center <- rep(1,ND)
  center[icl] <- 2
  res <- data.frame(pos= sub_pos,wts = sub_wts,group = cl,center=center)
  res <- res %>% dplyr::group_by(group) %>% dplyr::summarise(start=min(pos),end=max(pos),sum.wts = sum(wts),center = pos[center==2],.groups = 'drop')
}

#' A split_pac function
#'
#' Cluster poly(A) sites into groups by a certain distance iterative.
#'
#' @name split_pac
#' @param pa.data A list of polyA information, the output of the FindPTA function
#' @param d A cutoff limiting the max distance between two adjacent poly(A) sites in a poly(A) site cluster (PAC), default value is 24 nt.
#' @param mc.cores An integer indicating the number of processors/cores to perform the clustering, default value is 4.
#' @return A list contains re-defined PACs based on a weighted density peak clustering method
#' @export
#'

split_pac <- function(pa.data=NULL, d=24, mc.cores=4){
  pa.data$pa.coord$split.test <- "No"
  pa.data$pa.coord$split.test[which(pa.data$pa.coord$width>d)]  <- "Yes"

  pa.data.coord.part1 <- subset(pa.data$pa.coord,split.test=="No")
  pa.data.coord.part2 <- subset(pa.data$pa.coord,split.test=="Yes")
  if(nrow(pa.data.coord.part2)>0){
    print('Re-clustering by weighted density peak clustering algorithm!')
    print("Detail refer to R package-QuantifyPolyA")
    pa.data.coord.part2.region <- makeGRangesFromDataFrame(pa.data.coord.part2)

    count.data <-subset(pa.data$pa.table,PAid %in% pa.data.coord.part2$PAid)
    count.data$start <-   count.data$coord
    count.data$end <- count.data$coord

    count.region <- GenomicRanges::makeGRangesFromDataFrame(count.data)
    ov <- GenomicRanges::findOverlaps(count.region,pa.data.coord.part2.region)
    count.data <- count.data[unique(ov@from),]

    split.reads.data <- count.data[,c("readName","chr","coord","strand","level","PAid")]
    rm(count.data)
    split.reads.data$score <- 1
    split.reads.data.split <- split.reads.data %>%
      dplyr::group_by(chr,strand,coord) %>%
      dplyr::summarise(score = sum(score),.groups = 'keep')

    points.gr <- GRanges(seqnames = split.reads.data.split$chr,
                         ranges = IRanges(start =as.integer(split.reads.data.split$coord),  width = 1),
                         strand = as.character(split.reads.data.split$strand),
                         score = as.integer(split.reads.data.split$score))
    points.gr <- sort(points.gr)


    simple.clusters <- simpleCluster(points.gr, max.gapwidth = d)
    simple.clusters.df <- as.data.frame(simple.clusters)
    simple.clusters.df$split_label <- NA

    idx <- which(simple.clusters.df$width > d)
    pos <- IRanges::extractList(points.gr@ranges@start,simple.clusters$revmap[idx])
    wts <- IRanges::extractList(points.gr$score,simple.clusters$revmap[idx])
    split.clusters <- pbmcapply::pbmcmapply(findPeaks,pos,wts,
                                            SIMPLIFY = F, mc.cores = mc.cores )


    lens <- sapply(split.clusters,nrow)

    split.clusters.df <- dplyr::bind_rows(split.clusters[lens!=1], .id = "split_label")


    split.clusters.df$seqnames <- rep(simple.clusters.df$seqnames[idx[lens!=1]],times = lens[lens!=1])
    split.clusters.df$strand <- rep(simple.clusters.df$strand[idx[lens!=1]],times = lens[lens!=1])
    split.clusters.df$width <- split.clusters.df$end - split.clusters.df$start + 1
    split.clusters.df <- dplyr::rename(split.clusters.df,score=sum.wts)
    # Generate final clusters.
    polyA <- simple.clusters.df[-idx[lens!=1],c("seqnames","start","end","width","strand","score","center","split_label")]
    polyA <- rbind(polyA,split.clusters.df[,c("seqnames","start","end","width","strand","score","center","split_label")])
    colnames(polyA ) <- c("seqnames","start","end","width","strand","total.count","coord","split_label")

    split.reads.data <- split.reads.data[order(split.reads.data$chr,
                                               split.reads.data$coord,
                                               split.reads.data$level),]
    split.reads.data$label <- paste(split.reads.data$chr,
                                    split.reads.data$coord,
                                    split.reads.data$strand,sep = "_")

    polyA$label <- paste(polyA$seqnames,
                         polyA$coord,
                         polyA$strand,sep = "_")

    index <- match( polyA$label,split.reads.data$label)
    polyA$coord_level <- split.reads.data$level[index]
    polyA$PAid <- split.reads.data$PAid[index]

    polyA$label <- NULL
    pa.data.coord.part1.keep <- pa.data.coord.part1[,c("seqnames","start","end","width","strand","total.count","coord","coord_level","PAid")]
    pa.data.coord.part1.keep$split_label <- NA
    pa.data.coord.part1.keep$use.split.test <- "No"
    polyA$use.split.test  <- "Yes"
    polyA <- plyr::rbind.fill(pa.data.coord.part1.keep,
                              polyA)

    pa.data$split.clusters <- polyA
  }
  return(pa.data)


}


get_referenceseq<-function(aln=NULL,bsgenome=NULL){
  aln$reference_sequence <- ""
  index.plus <- which(aln$flag=="0")
  index.minus <- which(aln$flag=="16")
  sub.info.plus <- data.frame(chrom=aln$chr[index.plus],
                              start=as.numeric(aln$coord[index.plus]-10-aln$unmapreads_width[index.plus]),
                              end=as.numeric(aln$coord[index.plus])-1)


  aln$reference_sequence[index.plus] <- as.character(reverseComplement(getSeq(bsgenome, as(sub.info.plus, "GRanges"))) )


  sub.info.minus <- data.frame(chrom=aln$chr[index.minus],
                               start=as.numeric(aln$coord[index.minus])+1,
                               end=as.numeric(aln$coord[index.minus])+10+aln$unmapreads_width[index.minus])
  aln$reference_sequence[index.minus] <- as.character(getSeq(bsgenome, as(sub.info.minus, "GRanges")))
  return(aln)
}


#' A generateFASTA function
#'
#' Extract sequence from upstream 100 nt to downstream 100 nt region of polyA sites.
#'
#' @name generateFASTA
#' @param reads A data frame contains alignment result and identified A-rich polyA sites.
#' @param bsgenome A BSgenome object for reference genome. e.g. 'BSgenome.Hsapiens.UCSC.hg38'.
#' @param output.name The name of output fastq file.
#' @return  A fastq file
#' @export
#'

generateFASTA <- function(reads=NULL,bsgenome=NULL,output.name=NULL){

  reads <- subset(reads, level!="Count")
  reads <- subset(reads,is_Arich=="Yes")

  if(nrow(  reads)>0){
    fa_name_1 <- paste0(reads$chr,"_",
                        reads$coord,"_",
                        reads$strand,"_")
    reads <- reads[!duplicated(fa_name_1),]

    fa_name <- paste0(
      reads$chr,"_",
      reads$coord,"_",
      reads$strand,"_0")

    index.plus <- which(reads$flag=="0")
    index.minus <- which(reads$flag=="16")

    ###For on strand '-'
    sub.info.plus <- data.frame(chrom=as.character(reads$chr[index.plus ]),
                                start=as.integer(reads$coord)[index.plus]-99,
                                end=as.integer(reads$coord)[index.plus]+100)
    seq.plus <- reverseComplement(getSeq(bsgenome, as(sub.info.plus, "GRanges")))
    names(seq.plus) <- fa_name[index.plus]

    ###For on strand '+'
    sub.info.minus <- data.frame(chrom=as.character(reads$chr[index.minus]),
                                 start=as.integer(reads$coord)[index.minus]-100,
                                 end=as.integer(reads$coord)[index.minus]+99)
    seq.minus <- getSeq(bsgenome, as(sub.info.minus, "GRanges"))
    names(seq.minus) <- fa_name[index.minus]

    seq <- c(seq.plus,seq.minus)

    writeXStringSet(seq, output.name,format="fasta", width=200)
  }else{
    print("There are no polyA sites near the A-rich region.")
  }

}



#' A findTailAT function
#'
#' findTailAT finds and trims A or T tails from a fastq file.
#'
#' @name findTailAT
#' @param infile input fastq file.
#' @param poly should be one of "A", "T"
#' @param ml min length after trim. Default is 20.
#' @param mg margin from the start (poly=T) or to the end (poly=A) (default is 5)
#' @param mp min length of successive poly (default is 8)
#' @param mtail min length of trimmed tail (default is 8)
#' @param mper min percent of 'A\\/T' in trimmed tail (default is 0.75)
#' @param mm mismatch between TxxTTT (default is 2)
#' @param mr minT in reg (default is 3)
#' @param reg regexp string, default is 1 for loose searching, or 2 for stricter searching.
#'   1 = \code{'^.{0,mg}?(T{mr,}[^T]{0,mm})*T{mp,}([^T]{0,mm}T{mr,})*'}
#'   2 = \code{'^.{0,mg}?(T{mp,}([^T]{0,mm}T{mr,})*)'}.
#' @param bar the last bast position of three barcode, default is 0.
#' @param deep T(default) or FALSE, if T, then when no match was found for reg, will do deeper searching for tails like TTTTTCTTTTCTCTTTTTTTT..
#' @param review T (default) or F, if T then when no match was found for reg or deep, do deeper search for tails like AACCCCTTTTTTTTTTTTTTTTT...
#' @param debug T (default) or F, if T then output debugging files.
#' @param odir output path, if nULL then the same as infile dir.
#' @param suf suffix for output files, like <infile><.suf>.T/A.fq, default is NULL.
#' @return character. Output several files:
#'   \itemize{
#'     \item "*.T.fq" - Contains reads with names annotated to indicate polyT stretches at their 3' ends.
#'     \item "*.A.fq" - Contains reads with names annotated to indicate polyA stretches at their 3' ends.
#'   }
#' These files can serve as input for alignment tools.
#' @examples
#' file_T <- system.file("extdata", "SRR1168402_T.fastq", package = "PolyAseqTrap")
#' findTailAT(infile=file_T, odir=NULL,
#' poly='T', ml=20, mp=5, mg=10, mm=2, deep=FALSE,
#' mtail=6, mper=0.75, mr=3, review=TRUE, debug=TRUE,
#' bar=0, reg=1, suf=NULL)
#'
#' file_A <- system.file("extdata", "SRR11837378_A.fastq", package = "PolyAseqTrap")
#' findTailAT(infile=file_A, odir=NULL,
#' poly='A', ml=20, mp=5, mg=10, mm=2, deep=0, mtail=6, mper=0.75, mr=3,
#' review=TRUE, debug=TRUE, bar=0, reg=1, suf=NULL)
#' @export

findTailAT<-function(infile=NULL,
                     poly,
                     ml=20,
                     mg=5,
                     mp=8,
                     mtail=8,
                     mper=0.75,
                     mm=2,
                     mr=3,
                     reg=1,
                     bar=0,
                     deep=TRUE,
                     review=TRUE,
                     debug=TRUE,
                     odir=NULL,
                     suf=NULL
) {


  .isFq <- function(filePath) {
    conn <- file(filePath, "r")
    lines <- readLines(conn, n = 4)
    close(conn)
    if (length(lines) < 4) {
      return(FALSE)
    }
    if (!grepl("^@", lines[1])) {
      return(FALSE)
    }
    return(TRUE)
  }


  reg1='^.{0,mg}?(T{mr,}[^T]{0,mm})*T{mp,}([^T]{0,mm}T{mr,})*'
  reg2='^.{0,mg}?(T{mp,}([^T]{0,mm}T{mr,})*)'


  if (!(poly %in% c('A', 'T', 'A|T', 'A&T'))) stop("error poly(=A T A|T A&T)")

  findA=0; findT=0
  if ( length(grep('A', poly))>0) findA=1
  if ( length(grep('T', poly))>0) findT=1

  if (deep & findA) {
    stop("TODO: deep and findA is not ready yet")
  }


  if (!file.exists(infile)) {
    stop("file not exists! (", infile, ")\n")
  }

  if (!.isFq(infile)) stop("only support fastq file")

  inpath=paste0(dirname(infile), '/')

  infname=basename(infile)
  infname=gsub("\\..*$", "", infname)

  if (is.null(odir)) odir=inpath

  if (!is.null(suf)) suf=paste0('.', suf)


  cat(sprintf("poly=%s\ninfile=%s\noutput dir=%s\nmin length after trim(ml)=%d\nmin continue tail length(mp)=%d\nmax margin to the end(mg)=%d\nmax mismatches between xTTT(mm)=%d\n",
              poly, infile, odir, ml, mp, mg, mm))

  cat(sprintf("min reg T(mr)=%d\nmin tail length(mtail)=%d\nmin T/A percent(mper)=%.2f\n",
              mr, mtail, mper ))

  cat(sprintf("deep=%s\n", deep))

  if(debug) cat("DEBUG MODE, output to .test.file too \n")

  cntNotailT=0
  cntShortT=0
  cntTotal=0
  cntFinalT=0
  cntNotailA=0
  cntShortA=0
  cntFinalA=0
  cntMissA=0
  cntMissT=0
  cntBadT=0
  cntBadA=0


  regT=sprintf("^.{0,%d}?(T{%d,}[^T]{0,%d})*T{%d,}([^T]{0,%d}T{%d,})*",
               mg, mr, mm, mp, mm, mr)
  if (reg==2) {
    regT=sprintf("^.{0,%d}?T{%d,}([^T]{0,%d}T{%d,})*",
                 mg, mp, mm, mr)
  }
  cat("regT=", regT, "\n")

  regA=sprintf("(A{%d,}[^A]{0,%d})*A{%d,}[^A]{0,1}A+([^A]{0,%d}A{%d,})*.{0,%d}?$",
               mr, mm, mp, mm, mr, mg )
  if (reg==2) {
    regA=sprintf("A{%d,}([^A]{0,%d}A{%d,})*.{0,%d}?$",
                 mp, mm, mr, mg  )
  }
  cat("regA=", regA, "\n")

  if (deep & findT) {
    cat(sprintf("deep reg=^.{0,%d}?(T{1,}[^T]{0,2})*T{8,}\n", mg))
  }
  if(review & findT){
    cat(sprintf("review_reg=T{%d,}([^T]{0,%d}T{%d,}){0,}\n", mp, mm, mr))
  }

  md1=0
  md2=0
  md3=0


  OAfile <- paste0(odir, infname, suf, ".A.fq")
  TESTAfile = paste0(odir, infname, suf, ".testA.fq")

  OTfile <- paste0(odir, infname, suf, ".T.fq")
  SHORTTfile =paste0(odir, infname, suf, ".shortT.fq")


  cntTotal=0

  fh1=file(infile, "r")
  lines <- readLines(fh1, n = 4)

  if (findT) {
    OT=file(OTfile, "w")
    SHORTT=file(SHORTTfile, "w")
  }


  if (findA) {
    OA=file(OAfile, "w")
    TESTA=file(TESTAfile, "w")
  }

  while(length(lines) > 0) {

    tempID = lines[1]
    lines[1]=unlist(strsplit(tempID, split='\\s+', perl=TRUE))[1]
    cntTotal=cntTotal+1

    haveT=0
    haveA=0
    shortT=0
    notailT=0
    shortA=0
    notailA=0
    badTailT=0
    badTailA=0

    if (findT) {
      seq=lines[2]
      seqT=seq
      seqT=gsub(regT, '', seq, perl=TRUE)
      length=nchar(seqT);
      seqlen=nchar(seq)
      if (length<seqlen) {
        Ts=substr(seq, bar+1, seqlen-length)
        Tcnt=stringr::str_count(Ts, 'T')

        if(length<ml){
          shortT=1
          badTailT=0
        } else if ((seqlen-length-bar)<mtail | Tcnt/nchar(Ts)<mper) {
          badTailT=1
          shortT=0
        } else {
          md1=md1+1
          haveT=1
          trimS=Ts  #trimmed part sequence
          qualT=substr(lines[4], nchar(lines[4])-length+1, stop=1e6);
          writeLines(c( paste0(lines[1],'_',trimS ), seqT, lines[3], qualT), OT)
          cntFinalT=cntFinalT+1
          badTailT=0
          shortT=0
        }
      }

      if (haveT & deep) {
        seqT=lines[2]
        deepreg=sprintf("^.{0,%d}?(T{1,}[^T]{0,2})*T{%d,}", mg, mp)
        seqT=gsub(deepreg, '', seqT, perl=TRUE)
        seqlen=nchar(lines[2])
        length=nchar(seqT)
        if (length<seqlen){
          Ts=substr(seq, bar+1, seqlen-length)
          Tcnt=stringr::str_count(Ts, 'T')

          if (length<ml) {
            shortT=1;
            badTailT=0;
          } else if	((seqlen-length-bar)<mtail | Tcnt/nchar(Ts)<mper) {#tail???? or T%<mper
            badTailT=1;
            shortT=0;
          }else{
            md2=md2+1
            haveT=1;
            trimS=Ts; #trimmed part sequence
            qualT=substr(lines[4], nchar(lines[4])-length+1, stop=1e6)
            writeLines(c( paste0(lines[1],'_',trimS ), seqT, lines[3], qualT), OT)
            cntFinalT=cntFinalT+1
            badTailT=0
            shortT=0
          }
        }
      }

      if(haveT & review){
        seqT=lines[2]
        deepmod=sprintf("^.{0,%d}?T{$mp,}([^T]{0,%d}T{%d,}){0,}", mg, mm, mr)

        matches=gregexpr(deepmod, seqT, perl=TRUE)
        pos_start <- matches[[1]][1]
        match_length <- attr(matches[[1]], "match.length")[1]

        if(pos_start != -1){

          pos_start=pos_start-1
          pos_end = pos_start+match_length

          seqT =substr(seqT, pos_end+1, stop=1e6)

          length=nchar(seqT)

          Ts=substr(lines[2], bar+1, pos_end)
          Tcnt=stringr::str_count(Ts, 'T')

          if (length<ml) {
            shortT=1;
            badTailT=0;
          } else if ( nchar(Ts)<mtail | Tcnt/nchar(Ts)<mper) {#tail???? or T%<mper
            badTailT=1;
            shortT=0;
          } else{
            md3=md3+1
            haveT=1;
            trimS=Ts; #trimmed part sequence
            qualT=substr(lines[4], pos_end, stop=1e6)
            writeLines(c( paste0(lines[1],'_',trimS ), seqT, lines[3], qualT), OT)
            cntFinalT=cntFinalT+1
            badTailT=0;
            shortT=0;
          }
        }
      }


      if (!haveT & !shortT) {
        notailT=1;
        cntNotailT=cntNotailT+1
        writeLines(c( paste0(lines[1],'_', '' ), lines[2], lines[3], lines[4]), OT)
      }
      if(badTailT) {
        cntBadT=cntBadT+1
      }
      if(shortT) {
        cntShortT=cntShortT+1
        if (debug) writeLines(c(lines[1], lines[2]), SHORTT)
      }

    }#end find poly(T)


    if (findA) {

      seq=lines[2]
      seqA=gsub(regA, '', seq, perl=TRUE)
      length=nchar(seqA)
      seqlen=nchar(seq)

      if (length<seqlen) {
        As=substr(seq, length+1, seqlen)
        Acnt = stringr::str_count(As, 'A')

        if (length<ml) { #default ml=20 or 40
          shortA=1;
          badTailA=0;
        } else if (seqlen-length<mtail | Acnt/nchar(As)<mper) { #tail????
          badTailA=1;
        } else {
          haveA=1;
          shortA=0;
          trimS=As; #trimmed part sequence
          qualA=substr(lines[4],1,length);
          writeLines(c( paste0(lines[1],'_',trimS ), seqA, lines[3], qualA), OA)
          cntFinalA=cntFinalA+1
        }
      } else {
        notailA=1;
        badTailA=0;
        shortA=0;
      }

      if(!haveA & review){

        deepmod="A{10,}([^A]{0,2}A{3,})?.{0,50}?$"
        seq=lines[2]
        seqlen=nchar(lines[2])
        seqA=gsub(deepmod, '', seq, perl=TRUE)
        length=nchar(seqA)

        if (length<seqlen) {

          As=substr(seq, length+1, seqlen);
          Acnt = stringr::str_count(As, 'A')

          if (length<ml) { #default ml=20 or 40
            shortA=1;
            badTailA=0;
          } else if (seqlen-length<mtail) {
            badTailA=1;
            shortA=0;
          } else {
            haveA=1;
            trimS=As; #trimmed part sequence
            qualA=substr(lines[4], 1, length)
            writeLines(c( paste0(lines[1],'_',trimS ), seqA, lines[3], qualA), OA)
            cntFinalA=cntFinalA+1
          }
        }else{
          notailA=1;
        }

      }

      if (!haveA & !shortA) {
        notailA=1;
        cntNotailA=cntNotailA+1
        writeLines(c( paste0(lines[1],'_', '' ), lines[2], lines[3], lines[4]), OA)
      }

      if(!haveA & shortA){
        cntShortA=cntShortA+1
        writeLines(paste0("Short:*******", seq) , TESTA)
      }
      if(!haveA & badTailA & !shortA){
        cntBadA=cntBadA+1
        writeLines(paste0("BadA:*******", seq) , TESTA)
      }
      if(!haveA & notailA & !shortA){
        cntNotailA=cntNotailA+1
        writeLines(paste0("noTail:*******", seq) , TESTA)
      }

    } #~findA

    lines <- readLines(fh1, n = 4)

  } #end while


  close(fh1)

  if(findA){
    close(OA)
    close(TESTA)
  }


  if(findT){
    close(OT)
    if(debug) close(SHORTT)
  }



  cat(sprintf("\ntotal\t%d\n", cntTotal))

  if (findA) {
    cat(sprintf("[polyA]\nfinal\t%d\nnotail\t%d\ntooshort\t%d\nbadTail\t%d\nmissby(A|T)\t%d\n",
                cntFinalA, cntNotailA, cntShortA, cntBadA, cntMissA))
    if (cntTotal!=cntFinalA+cntNotailA+cntShortA+cntMissA+cntBadA) {
      cat("cntTotal!=cntFinalA+cntNotailA+cntShortA+cntMissA+cntBadA\n")
    }
  }


  if (findT) {
    cat(sprintf("[polyT]\nfinal\t%d\nnotail\t%d\ntooshort\t%d\nbadTail\t%d\nmissby(A|T)\t%d\n",
                cntFinalT, cntNotailT, cntShortT, cntBadT, cntMissT))
    cat(sprintf("region\t%d\ndeep\t%d\nreview\t%d\n", md1, md2, md3))
    if (cntTotal!=cntFinalT+cntNotailT+cntShortT+cntMissT+cntBadT) {
      cat("cntTotal!=cntFinalT+cntNotailT+cntShortT+cntMissT+cntBadT\n")
    }
  }

}






