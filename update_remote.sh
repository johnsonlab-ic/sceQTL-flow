#!/bin/bash
# Script to update the Git remote URL after repository rename

# Display current remotes
echo "Current Git remotes:"
git remote -v

# Update the organization remote to point to the new repository name using HTTPS
git remote set-url organization https://github.com/johnsonlab-ic/sceQTL-flow.git

# Check if pushall alias exists
if git config --get alias.pushall > /dev/null; then
  echo "Git alias 'pushall' already exists:"
  git config --get alias.pushall
else
  # Create pushall alias to push to both remotes
  git config alias.pushall '!git push origin $(git rev-parse --abbrev-ref HEAD) && git push organization $(git rev-parse --abbrev-ref HEAD)'
  echo "Created Git alias 'pushall' to push to both remotes"
fi

# Display updated remotes
echo -e "\nUpdated Git remotes:"
git remote -v

echo -e "\nRemote URL for 'organization' has been updated to point to 'sceQTL-flow' using SSH"
echo -e "\nIMPORTANT: To push to the organization repository:"
echo "1. Make sure your SSH key is added to your GitHub account"
echo "2. Ensure your SSH key is added to the ssh-agent (run: ssh-add ~/.ssh/id_rsa)"
echo "3. Verify you have the necessary permissions in the johnsonlab-ic organization"
echo "4. You can test your SSH connection with: ssh -T git@github.com"
echo -e "\nAfter confirming SSH access works, you can use 'git pushall' to push to both repositories"
