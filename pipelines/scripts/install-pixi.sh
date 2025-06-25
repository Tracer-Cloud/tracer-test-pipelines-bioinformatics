#!/usr/bin/env bash
set -e

# Install Pixi if not present
if ! command -v pixi &> /dev/null; then
    echo "Pixi not found. Installing Pixi..."
    curl -fsSL https://pixi.sh/install.sh | bash
    export PATH="$HOME/.pixi/bin:$PATH"
    # Source shell profile for current session
    if [ -f "$HOME/.zshrc" ]; then
        source "$HOME/.zshrc"
    fi
    if [ -f "$HOME/.bashrc" ]; then
        source "$HOME/.bashrc"
    fi
fi

echo "Pixi is installed and ready." 