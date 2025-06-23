#!/bin/bash

echo "[INFO] Preparing environment: installing Pixi and Nextflow (if needed)..."

# 1. Ensure Pixi is installed
if ! command -v pixi &> /dev/null; then
    echo "[INFO] Pixi not found. Installing..."
    curl -sSf https://pixi.sh/install.sh | bash
    export PATH="$HOME/.pixi/bin:$PATH"
    echo 'export PATH="$HOME/.pixi/bin:$PATH"' >> ~/.bashrc
fi

# 2. Optionally initialize Pixi environment (safe to skip if not needed)
# You will manually run `tracer init` later
echo "[INFO] Skipping automatic tracer init. You can run it manually later with:"
echo "       tracer init"

# 3. Install Nextflow if not available
if ! command -v nextflow &> /dev/null; then
    echo "[INFO] Installing Nextflow..."
    curl -s https://get.nextflow.io | bash
    chmod +x nextflow
    sudo mv nextflow /usr/local/bin/
fi

# 4. Confirm everything is ready
if command -v nextflow &> /dev/null; then
    echo "[✅] Setup complete. You can now run manually:"
    echo
    echo "    tracer init"
    echo "    nextflow run main.nf -c nextflow.config"
else
    echo "[❌] Nextflow install failed. Please check for issues."
    exit 1
fi
