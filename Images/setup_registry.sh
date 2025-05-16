# Add your GitHub Personal Access Token
# Create token at: GitHub → Settings → Developer Settings → Personal Access Tokens
# Token needs at least 'write:packages' permission
export CR_PAT=YOUR_TOKEN
# Login to GitHub Container Registry
echo $CR_PAT | docker login ghcr.io -u HaglundA --password-stdin

# Build and tag Docker images for GitHub Container Registry
# Format: ghcr.io/USERNAME/IMAGE_NAME:TAG

# Build Expression image
docker build -t ghcr.io/haglunda/eqtl-expression:latest -f Images/Dockerfile.expression .

# Build Genotype image
docker build -t ghcr.io/haglunda/eqtl-genotype:latest -f Images/Dockerfile.genotype .

# Build Python image
docker build -t ghcr.io/haglunda/eqtl-python:latest -f Images/Dockerfile.python .

# Build Reports image
docker build -t ghcr.io/haglunda/eqtl-reports:latest -f Images/Dockerfile.reports .

# Push images to GitHub Container Registry
docker push ghcr.io/haglunda/eqtl-expression:latest
docker push ghcr.io/haglunda/eqtl-genotype:latest
docker push ghcr.io/haglunda/eqtl-python:latest
docker push ghcr.io/haglunda/eqtl-reports:latest

echo "All images have been built and pushed to GitHub Container Registry!"