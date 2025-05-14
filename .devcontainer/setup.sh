#!/bin/bash
set -e

echo "[*] Installing Nextflow"
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/