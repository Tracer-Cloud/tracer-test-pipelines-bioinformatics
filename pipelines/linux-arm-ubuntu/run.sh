#!/bin/bash

echo "[INFO] Preparing environment: installing Java, Pixi and Nextflow (if needed)..."

# Step 0: Ensure Java is available
if ! command -v java &> /dev/null; then
    echo "[INFO] Java not found. Installing OpenJDK 17..."
    sudo apt update && sudo apt install -y openjdk-17-jdk
fi

# Verify Java installation
if command -v java &> /dev/null; then
    echo "[INFO] Java installed:"
    java -version
else
    echo "[ERROR] Java installation failed. Exiting."
    exit 1
fi

# Step 1: Ensure Pixi is installed
if ! command -v pixi &> /dev/null; then
    echo "[INFO] Pixi not found. Installing..."
    curl -sSf https://pixi.sh/install.sh | bash
    export PATH="$HOME/.pixi/bin:$PATH"
    echo 'export PATH="$HOME/.pixi/bin:$PATH"' >> ~/.bashrc
fi

# Step 2: Skip automatic tracer init
echo "[INFO] Skipping automatic tracer init. You can run it manually later with:"
echo "       tracer init"

# Step 3: Install Nextflow if not available
if ! command -v nextflow &> /dev/null; then
    echo "[INFO] Installing Nextflow..."
    curl -s https://get.nextflow.io | bash
    chmod +x nextflow
    sudo mv nextflow /usr/local/bin/
fi

# Step 4: Confirm readiness
if command -v nextflow &> /dev/null; then
    echo "[✅] Setup complete. You can now run manually:"
    echo
    echo "    tracer init"
    echo "    nextflow run main.nf -c nextflow.config"
else
    echo "[❌] Nextflow install failed. Please check the above output."
    exit 1
fi
