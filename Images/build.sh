#!/bin/bash
# Build and optionally push Docker image from the devcontainer
# Usage: ./build.sh [--push] [--no-cache] [--tag TAG]
#
# Examples:
#   ./build.sh                              # Build local image
#   ./build.sh --push                       # Build and push to registry
#   ./build.sh --push --no-cache            # Force rebuild and push
#   ./build.sh --tag v1.2.3 --push          # Build with custom tag and push

set -e

REGISTRY="ghcr.io/johnsonlab-ic"
IMAGE_NAME="sceqtl-flow"
DEFAULT_TAG="latest"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Parse arguments
PUSH=false
NO_CACHE=""
CUSTOM_TAG=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --push)
            PUSH=true
            shift
            ;;
        --no-cache)
            NO_CACHE="--no-cache"
            shift
            ;;
        --tag)
            CUSTOM_TAG="$2"
            shift 2
            ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            exit 1
            ;;
    esac
done

TAG="${CUSTOM_TAG:-$DEFAULT_TAG}"
FULL_IMAGE="${REGISTRY}/${IMAGE_NAME}:${TAG}"

# Verify Dockerfile exists
if [ ! -f "$SCRIPT_DIR/../.devcontainer/Dockerfile" ]; then
    echo -e "${RED}Error: Dockerfile not found at $SCRIPT_DIR/../.devcontainer/Dockerfile${NC}"
    exit 1
fi

echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${GREEN}Unified sceQTL-flow Docker Image Builder${NC}"
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""
echo -e "${BLUE}Build Configuration:${NC}"
echo -e "  ${GREEN}Registry:${NC} ${REGISTRY}"
echo -e "  ${GREEN}Image:${NC} ${IMAGE_NAME}"
echo -e "  ${GREEN}Tag:${NC} ${TAG}"
echo -e "  ${GREEN}Full Image:${NC} ${FULL_IMAGE}"
echo -e "  ${GREEN}Dockerfile:${NC} .devcontainer/Dockerfile"
echo -e "  ${GREEN}Push to Registry:${NC} $([ "$PUSH" = true ] && echo 'YES' || echo 'NO')"
if [ -n "$NO_CACHE" ]; then
    echo -e "  ${YELLOW}Force Rebuild:${NC} YES (--no-cache)"
fi
echo ""

# Build the image
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${GREEN}Building Docker image: ${FULL_IMAGE}${NC}"
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""

docker build \
    --network=host \
    $NO_CACHE \
    -t "$FULL_IMAGE" \
    -f "$SCRIPT_DIR/../.devcontainer/Dockerfile" \
    "$SCRIPT_DIR/.."

if [ $? -eq 0 ]; then
    echo ""
    echo -e "${GREEN}✓ Successfully built: ${FULL_IMAGE}${NC}"
else
    echo ""
    echo -e "${RED}✗ Build failed!${NC}"
    exit 1
fi

# Push if requested
if [ "$PUSH" = true ]; then
    echo ""
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${GREEN}Pushing image to registry: ${FULL_IMAGE}${NC}"
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo ""
    
    docker push "$FULL_IMAGE"
    
    if [ $? -eq 0 ]; then
        echo ""
        echo -e "${GREEN}✓ Successfully pushed: ${FULL_IMAGE}${NC}"
    else
        echo ""
        echo -e "${RED}✗ Push failed!${NC}"
        exit 1
    fi
fi

echo ""
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${GREEN}Build complete!${NC}"
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""
echo -e "${YELLOW}Next steps:${NC}"
echo "  • Update nextflow.config to use the new unified image:"
echo "      All processes → ${FULL_IMAGE}"
echo "  • Test the pipeline with the new image"
echo "  • (Optional) Tag with version numbers:"
echo "      ./build.sh --tag v1.0.0 --push"
echo ""
