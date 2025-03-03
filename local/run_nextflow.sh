
cd /Users/ah3918/Dropbox/LANDMARK/PROJECTS/eQTL_PIPELINE

nextflow run -c /Users/ah3918/Dropbox/LANDMARK/PROJECTS/eQTL_PIPELINE/NEXTFLOW/nextflow_local.config \
/Users/ah3918/Dropbox/LANDMARK/PROJECTS/eQTL_PIPELINE/NEXTFLOW/eqtl_pipeline.nf -profile offline \
--gds_file /Users/ah3918/Dropbox/LANDMARK/PROJECTS/eQTL_PIPELINE/data/inputs/test_geno.gds \
--single_cell_file /Users/ah3918/Dropbox/LANDMARK/PROJECTS/eQTL_PIPELINE/data/inputs/roche_ms_decontx.rds \
--outdir /Users/ah3918/Dropbox/LANDMARK/PROJECTS/eQTL_PIPELINE/data/outputs


### on the HPC 

cd /rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines/eQTL_PIPELINE

