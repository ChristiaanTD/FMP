#!/bin/bash

# Update package list
sudo apt-get update

# Install Java (required for Nextflow)
sudo apt-get install -y openjdk-11-jdk

# Install Nextflow
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

# Verify Nextflow installation
nextflow -v


# Install Docker
sudo apt-get install -y apt-transport-https ca-certificates curl software-properties-common
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg
echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update
sudo apt-get install -y docker-ce docker-ce-cli containerd.io


# Verify Docker installation
sudo systemctl status docker

# Add user to group
sudo usermod -aG docker $USER 

# Print message for user to logout and login back
echo "Installation complete. Please log out and log back in for the group changes to take effect."
