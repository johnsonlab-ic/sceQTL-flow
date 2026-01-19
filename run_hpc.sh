indir=/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/projects/pipelines/TEST_DATA
export NXF_LOG_FILE="$EPHEMERAL/NEXTFLOW/nextflow.log"


nextflow run main.nf \
-profile imperial \
--outdir /rds/general/user/ah3918/ephemeral/eqtl_pipepline_outs/test_data \
--gds_file ${indir}/final_geno_440samples.gds \
--single_cell_file ${indir}/roche_ms_decontx.rds \
--celltype_column "CellType" \
--individual_column "Individual_ID" \
--counts_slot "counts" \
--counts_assay "decontXcounts" \
-with-report pipeline_report.html \
--report true



#YAZAR

git pull;nextflow run -c nextflow.config main.nf \
-profile imperial \
--gds_file ${indir}/YAZAR/merge_test_seqArray.gds \
--single_cell_file ${indir}YAZAR/Yazar_annotateddata.rds \
--celltype_column "eight_Cell_Types" \
--individual_column "donor_id" \
--counts_slot "data" \
--counts_assay "RNA" \
--optimize_pcs true \
--report true




#### epilepsy 

outdir=/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/projects/epilepsy_eqtls/sc_eqtls
export NXF_LOG_FILE="$EPHEMERAL/NEXTFLOW/nextflow.log"
cov_file=$outdir/epilepsy_metadata_matrixeqtl.csv

nextflow run -c nextflow.config main.nf \
-profile imperial \
--outdir $outdir \
--gds_file ${outdir}/BONN_post_imputation_QC.gds \
--single_cell_file ${outdir}/matched_seurat_demuxlet.rds \
--cov_file $cov_file \
--celltype_column "cell_type" \
--individual_column "Individual_ID" \
--counts_assay "originalexp" \
--counts_slot "counts" \
-with-report pipeline_report.html \
--report true \
--covariates_to_include "Sex,Pathology,AgeAtSurgery"



### yazar

outdir=/rds/general/user/ah3918/ephemeral/YAZAR_OUTS
indir=/rds/general/user/ah3918/projects/puklandmarkproject/ephemeral/Liv_Yazar/
export NXF_LOG_FILE="$EPHEMERAL/NEXTFLOW/nextflow.log"

nextflow run -c nextflow.config main.nf \
-profile imperial \
--outdir $outdir \
--gds_file ${indir}/merge_test_seqArray.gds \
--single_cell_file ${indir}/Yazar_annotateddata.rds \
--celltype_column "eight_Cell_Types" \
--individual_column "donor_id" \
--counts_slot "data" \
--counts_assay "RNA" \
-with-report pipeline_report.html \
--report true 

