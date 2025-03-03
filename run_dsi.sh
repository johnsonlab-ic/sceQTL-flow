workdir=/var/lib/docker/alex_tmp/NF_WORK/
datadir=/var/lib/docker/alex_tmp/data/
mount_path=/home/ah3918/rds/
indir=${mount_path}/live/Users/Alex/pipelines//TEST_DATA/
outdir=${mount_path}/ephemeral/eqtl_pipepline_outs/test_data

# mkdir -p ${datadir}/test_data/
# cp -r ${indir}/* ${datadir}/test_data/

git pull;nextflow run -c nextflow.config main.nf \
-w ${workdir} \
-profile offline \
--outdir ${outdir} \
--gds_file ${datadir}/test_data/test_geno.gds \
--single_cell_file ${datadir}//test_data/roche_ms_decontx.rds \
--celltype_column "CellType" \
--individual_column "Individual_ID" \
--counts_slot "counts" \
--counts_assay "decontXcounts" \
--optimize_pcs true \
-N a.haglund19@imperial.ac.uk \
-with-report pipeline_report.html \
--report true
