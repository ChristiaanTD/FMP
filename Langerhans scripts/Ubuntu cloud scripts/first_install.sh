#!/bin/bash

#-------------This script can be used if necessary as a first install on your cloud server to intstall Nextflow and Docker-----------#

# first update package list
sudo apt-get update

# next we install Java (if you don't have it installed already)
sudo apt-get install -y openjdk-11-jdk

# install Nextflow
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

# Verify Nextflow installation
nextflow -v


# install Docker 
sudo apt-get install -y apt-transport-https ca-certificates curl software-properties-common
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg
echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update
sudo apt-get install -y docker-ce docker-ce-cli containerd.io


# verify Docker installation
sudo systemctl status docker

# add user to group
sudo usermod -aG docker $USER 

# print message so you don't forget to log out and back on for group changes to take effect!
echo "Installation complete. Please log out and log back in for the group changes to take effect."
