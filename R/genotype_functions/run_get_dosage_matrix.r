generate_genotype_matrix=function(vcfs,
                                  gds_file="merge_test_seqArray.gds",
                                  outdir=".",
                                  parallel=FALSE,
                                  preprocess=FALSE,
                                  autosomalonly=TRUE,
                                  minmaf=0.05){



  if (preprocess) {
    SeqArray::seqVCF2GDS(vcfs, gds_file,parallel=parallel)
  }

  genofile <- SeqArray::seqOpen(gds_file)

  if (autosomalonly) {
    SeqArray::seqSetFilterChrom(genofile,1:22)
  }

  SeqArray::seqSetFilterCond(genofile,maf=minmaf)

  #get snp info
  annot<-SeqArray::seqGetData(genofile,"annotation/id")

  #get sample names
  sample.id<-SeqArray::seqGetData(genofile,"sample.id")

  #get genotype dosage matrix (0,1,2)
  geno_mat<-t(SeqArray::seqGetData(genofile,"$dosage"))

  #create SNP locations file
  position<-SeqArray::seqGetData(genofile,"position")
  chrom<-SeqArray::seqGetData(genofile,"chromosome")
  chrompos_mat<-data.frame(annot,chrom,position)
  chrompos_mat<-chrompos_mat[!duplicated(chrompos_mat$annot),]

  #clean / reformat
  colnames(geno_mat)<-sample.id
  geno_mat<-as.data.frame(geno_mat)
  geno_mat$snp<-annot
  geno_mat<-geno_mat[!duplicated(geno_mat$snp),]
  rownames(geno_mat)<-geno_mat$snp

  #
  chrompos_mat$chrom<-paste0("chr",chrompos_mat$chrom)
  snpnumber<-length(rownames(geno_mat))

  message("Checking snps for build and converting to rsids..")

  
  chrompos_mat=check_snps(chrompos_mat)
  geno_mat=geno_mat[match(chrompos_mat$old_snp,geno_mat$snp),]
  geno_mat$snp=chrompos_mat$annot
  rownames(geno_mat)=geno_mat$snp
  message(paste0("A total of ",snpnumber," snps were kept."))

  data.table::fwrite(geno_mat,"genotype_012mat.csv")
  

  #add alt. Allele Freq
  allele<-SeqArray::seqGetData(genofile,"allele")
  allele<-strsplit(allele,",")
  allele<-as.data.frame(do.call(rbind,allele))


    
  #check if MAF info is available. If not, calculate it
  maf_avail<-tryCatch({
        out<-SeqArray::seqGetData(genofile,"annotation/info/MAF")
    },
     error=function(e){
        out<-"error"
         return(out)
    }
  )

  if(length(maf_avail)>1){
      maf<-SeqArray::seqGetData(genofile,"annotation/info/MAF")
      }else if(maf_avail=="error"){
      maf<- SeqArray::seqAlleleFreq(genofile, minor=TRUE)
  }
  
  af_df<-data.frame(ref=allele$V1,
    alt=allele$V2,
    maf=maf,
    snp=annot)
    af_df$ref<-as.character(af_df$ref)
    af_df$alt<-as.character(af_df$alt)
    af_df$snp<-as.character(af_df$snp)


    af_df<-af_df[match(chrompos_mat$old_snp,af_df$snp),]
    af_df$snp<-chrompos_mat$annot



      #remove extra columns
  chrompos_mat=chrompos_mat[,c("annot","chrom","position")]
  data.table::fwrite(chrompos_mat,"snp_chromlocations.csv")



  data.table::fwrite(af_df,"MAF_mat.csv")

  SeqArray::seqClose(genofile)
}