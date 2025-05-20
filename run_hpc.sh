indir=/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines//TEST_DATA/
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