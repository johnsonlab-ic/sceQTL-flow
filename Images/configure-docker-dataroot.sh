#!/bin/bash
# Script to configure Docker to use /data partition instead of /var/lib
# Run this script with sudo

set -e

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${GREEN}Docker Data Root Configuration${NC}"
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""

# Check if running as root
if [ "$EUID" -ne 0 ]; then
    echo -e "${RED}Please run as root (use sudo)${NC}"
    exit 1
fi

# Backup existing daemon.json if it exists
if [ -f /etc/docker/daemon.json ]; then
    echo -e "${YELLOW}Backing up existing /etc/docker/daemon.json${NC}"
    cp /etc/docker/daemon.json /etc/docker/daemon.json.backup.$(date +%Y%m%d_%H%M%S)
fi

# Create new data-root directory
NEW_DATA_ROOT="/data/docker"
echo -e "${BLUE}Creating new Docker data directory: ${NEW_DATA_ROOT}${NC}"
mkdir -p "$NEW_DATA_ROOT"

# Stop Docker service
echo -e "${YELLOW}Stopping Docker service...${NC}"
systemctl stop docker

# Copy Docker daemon configuration
echo -e "${BLUE}Configuring Docker to use ${NEW_DATA_ROOT}${NC}"
cat > /etc/docker/daemon.json <<EOF
{
  "data-root": "${NEW_DATA_ROOT}",
  "storage-driver": "overlay2"
}
EOF

# Optional: Copy existing Docker data to new location (commented out by default)
# Uncomment if you want to preserve existing images/containers
# echo -e "${YELLOW}Copying existing Docker data (this may take a while)...${NC}"
# rsync -aP /var/lib/docker/ "$NEW_DATA_ROOT/"

# Start Docker service
echo -e "${YELLOW}Starting Docker service...${NC}"
systemctl start docker

# Verify new data root
echo ""
echo -e "${GREEN}✓ Docker reconfigured successfully!${NC}"
echo ""
echo -e "${BLUE}Verification:${NC}"
docker info | grep -i "Docker Root Dir"
echo ""
echo -e "${YELLOW}Note: Existing images were NOT migrated. You'll need to rebuild or pull images.${NC}"
echo -e "${YELLOW}Old Docker data is still in /var/lib/docker (you can remove it manually if needed)${NC}"
echo ""
