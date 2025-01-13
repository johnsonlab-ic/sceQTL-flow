DOCKER_HUB_USERNAME="ah3918"
indir="/Users/ah3918/Dropbox/WORK/LANDMARK/PROJECTS/eQTL_PIPELINE/"
# Build the genotype_image
echo "Building genotype_image..."
docker buildx build --platform linux/amd64,linux/arm64 -t ${DOCKER_HUB_USERNAME}/genotype_image:latest -f ${indir}/Images/Dockerfile.genotype --push .

# Build the expression_image
echo "Building expression_image..."
docker buildx build --platform linux/amd64,linux/arm64 -t ${DOCKER_HUB_USERNAME}/expression_image:latest -f ${indir}/Images/Dockerfile.expression --push .

# Build the report_image
echo "Building report_image..."
docker buildx build --platform linux/amd64 -t ${DOCKER_HUB_USERNAME}/report_image:latest -f ${indir}/Images/Dockerfile.reports --push .

echo "All images have been built and pushed successfully."