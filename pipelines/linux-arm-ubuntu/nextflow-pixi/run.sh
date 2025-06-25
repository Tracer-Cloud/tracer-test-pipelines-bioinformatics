#!/usr/bin/env bash
set -e

if ! command -v pixi &> /dev/null; then
    echo "Pixi not found. Installing Pixi..."
    curl -fsSL https://pixi.sh/install.sh | bash
    export PATH="$HOME/.pixi/bin:$PATH"
fi

pixi install
pixi run pipeline
