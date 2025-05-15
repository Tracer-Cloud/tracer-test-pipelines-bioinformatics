#!/bin/bash
set -e

echo "[*] Installing Nextflow"
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/


curl -sSL https://install.tracer.cloud/installation-script-development.sh | bash && source ~/.bashrc
