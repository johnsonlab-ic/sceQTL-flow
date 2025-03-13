workdir=/var/lib/docker/alex_tmp/NF_WORK/
datadir=/var/lib/docker/alex_tmp/data/
mount_path=/home/ah3918/rds/
indir=${mount_path}/live/Users/Alex/pipelines//TEST_DATA/
outdir=${workdir}/test_data

# mkdir -p ${datadir}/test_data/
# cp -r ${indir}/* ${datadir}/test_data/




nextflow run main.nf \
-w ${workdir} \
-profile "offline" \
--outdir ${outdir} \
--gds_file ${datadir}/test_data/final_geno_440samples.gds \
--single_cell_file ${datadir}//test_data/roche_ms_decontx.rds \
--celltype_column "CellType" \
--individual_column "Individual_ID" \
--counts_slot "counts" \
--counts_assay "decontXcounts" \
-with-report pipeline_report.html \
--report true \
--cis_distance 1000000 \
--optimize_pcs true \
--filter_chr "chr1,chr2,chr3,chr4,chr5"


#testing with false

nextflow run main.nf \
-w ${workdir} \
-profile "offline" \
--outdir ${outdir} \
--gds_file ${datadir}/test_data/test_geno.gds \
--single_cell_file ${datadir}//test_data/roche_ms_decontx.rds \
--celltype_column "CellType" \
--individual_column "Individual_ID" \
--counts_slot "counts" \
--counts_assay "decontXcounts" \
--optimize_pcs false \
-with-report pipeline_report.html \
--report true


#testing with YAZAR
workdir=/var/lib/docker/alex_tmp/NF_WORK/
datadir=/var/lib/docker/alex_tmp/data/test_data/YAZAR/
mount_path=/home/ah3918/rds/
indir=${mount_path}/live/Users/Alex/pipelines//TEST_DATA/
outdir=${workdir}/test_data


nextflow run main.nf \
-w ${workdir} \
-profile "offline" \
--outdir ${outdir} \
--gds_file ${datadir}/merge_test_seqArray.gds \
--single_cell_file ${datadir}/Yazar_annotateddata.rds  \
--celltype_column "eight_Cell_Types" \
--individual_column "donor_id" \
--counts_slot "data" \
--counts_assay "RNA" \
--optimize_pcs true \
-with-report pipeline_report.html \
--report true
