#!/bin/bash

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to install Miniconda
install_miniconda() {
    print_status "Installing Miniconda..."
    
    # Create temporary directory
    TEMP_DIR=$(mktemp -d)
    cd "$TEMP_DIR"
    
    # Detect OS and architecture for appropriate installer
    OS=$(uname -s)
    ARCH=$(uname -m)
    
    # Use a specific stable version of Miniconda instead of latest
    MINICONDA_VERSION="Miniconda3-py39_4.12.0"
    
    if [ "$OS" = "Darwin" ]; then
        # macOS
        if [ "$ARCH" = "arm64" ]; then
            # Apple Silicon
            curl -L -o miniconda.sh "https://repo.anaconda.com/miniconda/${MINICONDA_VERSION}-MacOSX-arm64.sh"
        else
            # Intel Mac
            curl -L -o miniconda.sh "https://repo.anaconda.com/miniconda/${MINICONDA_VERSION}-MacOSX-x86_64.sh"
        fi
    elif [ "$OS" = "Linux" ]; then
        # Linux
        if [ "$ARCH" = "x86_64" ]; then
            curl -L -o miniconda.sh "https://repo.anaconda.com/miniconda/${MINICONDA_VERSION}-Linux-x86_64.sh"
        elif [ "$ARCH" = "aarch64" ]; then
            curl -L -o miniconda.sh "https://repo.anaconda.com/miniconda/${MINICONDA_VERSION}-Linux-aarch64.sh"
        else
            print_error "Unsupported architecture: $ARCH"
            exit 1
        fi
    fi
    
    # Check if download was successful
    if [ ! -f miniconda.sh ]; then
        print_error "Failed to download Miniconda installer"
        cd - > /dev/null
        rm -rf "$TEMP_DIR"
        return 1
    fi
    
    # Install Miniconda with -u option to update if exists
    print_status "Installing/Updating Miniconda (this may take a few minutes)..."
    if ! bash miniconda.sh -b -u -p "$HOME/miniconda3"; then
        print_error "Failed to install Miniconda"
        cd - > /dev/null
        rm -rf "$TEMP_DIR"
        return 1
    fi
    
    # Initialize conda
    if [ -f "$HOME/miniconda3/bin/conda" ]; then
        "$HOME/miniconda3/bin/conda" init bash
        "$HOME/miniconda3/bin/conda" init zsh
        
        # Add conda to current session
        export PATH="$HOME/miniconda3/bin:$PATH"
        
        print_success "Miniconda installed/updated successfully"
        print_status "Please restart your terminal or run 'source ~/.bashrc' (or ~/.zshrc) to activate conda"
    else
        print_error "Miniconda installation appears to have failed"
        cd - > /dev/null
        rm -rf "$TEMP_DIR"
        return 1
    fi
    
    # Clean up
    cd - > /dev/null
    rm -rf "$TEMP_DIR"
}


# Function to safely append to .zshrc
append_to_zshrc() {
    local line="$1"
    local zshrc_path="$HOME/.zshrc"
    
    # Create .zshrc if it doesn't exist
    if [ ! -f "$zshrc_path" ]; then
        print_status "Creating .zshrc file..."
        sudo touch "$zshrc_path"
        sudo chown "$USER:$(id -gn)" "$zshrc_path"
        sudo chmod 644 "$zshrc_path"
    fi
    
    # Check if we have write permissions
    if [ ! -w "$zshrc_path" ]; then
        print_status "Fixing permissions for .zshrc..."
        sudo chown "$USER:$(id -gn)" "$zshrc_path"
        sudo chmod 644 "$zshrc_path"
    fi
    
    # Check if the line already exists
    if ! grep -qF "$line" "$zshrc_path"; then
        print_status "Adding configuration to .zshrc..."
        echo "$line" >> "$zshrc_path"
    fi
}

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to check Java version
check_java() {
    if command_exists java; then
        JAVA_VERSION=$(java -version 2>&1 | head -n1 | cut -d'"' -f2 | cut -d'.' -f1)
        if [ "$JAVA_VERSION" -ge 11 ]; then
            print_success "Java $JAVA_VERSION is already installed and meets requirements (Java 11+)"
            return 0
        else
            print_warning "Java $JAVA_VERSION is installed but Nextflow requires Java 11+"
            return 1
        fi
    else
        print_status "Java is not installed"
        return 1
    fi
}

# Function to check Python version
check_python() {
    if command_exists python3; then
        PYTHON_VERSION=$(python3 --version 2>&1 | cut -d' ' -f2 | cut -d'.' -f1,2)
        print_success "Python $PYTHON_VERSION is already installed"
        return 0
    else
        print_status "Python3 is not installed"
        return 1
    fi
}

# Function to install Python on macOS
install_python_macos() {
    print_status "Installing Python3 on macOS..."
    
    if command_exists brew; then
        print_status "Using Homebrew to install Python3..."
        brew install python@3.11
    else
        print_status "Installing Python3 using official installer..."
        
        # Create temporary directory
        TEMP_DIR=$(mktemp -d)
        cd "$TEMP_DIR"
        
        # Download Python installer
        if [[ $(uname -m) == "arm64" ]]; then
            # Apple Silicon
            curl -L -o python-installer.pkg "https://www.python.org/ftp/python/3.11.9/python-3.11.9-macos11.pkg"
        else
            # Intel Mac
            curl -L -o python-installer.pkg "https://www.python.org/ftp/python/3.11.9/python-3.11.9-macos11.pkg"
        fi
        
        # Install Python
        print_status "Installing Python3 (may require administrator password)..."
        sudo installer -pkg python-installer.pkg -target /
        
        # Add Python to PATH
        append_to_zshrc 'export PATH="/usr/local/bin:$PATH"'
        append_to_zshrc 'export PATH="/Library/Frameworks/Python.framework/Versions/3.11/bin:$PATH"'
        export PATH="/Library/Frameworks/Python.framework/Versions/3.11/bin:$PATH"
        
        # Clean up
        cd - > /dev/null
        rm -rf "$TEMP_DIR"
        
        print_success "Python3 installed successfully"
    fi
}

# Function to install Java on macOS
install_java_macos() {
    print_status "Installing Java on macOS..."
    
    if command_exists brew; then
        print_status "Using Homebrew to install OpenJDK..."
        brew install openjdk@17
        
        # Detect if using Apple Silicon or Intel Mac for correct path
        if [[ $(uname -m) == "arm64" ]]; then
            # Apple Silicon Mac
            append_to_zshrc 'export PATH="/opt/homebrew/bin:$PATH"'
            append_to_zshrc 'export JAVA_HOME="/opt/homebrew/opt/openjdk@17"'
            export PATH="/opt/homebrew/bin:$PATH"
            export JAVA_HOME="/opt/homebrew/opt/openjdk@17"
        else
            # Intel Mac
            append_to_zshrc 'export PATH="/usr/local/bin:$PATH"'
            append_to_zshrc 'export JAVA_HOME="/usr/local/opt/openjdk@17"'
            export PATH="/usr/local/bin:$PATH"
            export JAVA_HOME="/usr/local/opt/openjdk@17"
        fi
        
    else
        print_warning "Homebrew not found. Installing Java using alternative method..."
        
        # Check if we're on macOS 10.15+ (supports direct download method)
        MACOS_VERSION=$(sw_vers -productVersion | cut -d. -f1,2)
        
        print_status "Downloading OpenJDK 17 from Adoptium..."
        
        # Create temporary directory
        TEMP_DIR=$(mktemp -d)
        cd "$TEMP_DIR"
        
        # Detect architecture and download appropriate JDK
        if [[ $(uname -m) == "arm64" ]]; then
            # Apple Silicon
            curl -L -o openjdk-17.tar.gz "https://github.com/adoptium/temurin17-binaries/releases/download/jdk-17.0.9%2B9/OpenJDK17U-jdk_aarch64_mac_hotspot_17.0.9_9.tar.gz"
        else
            # Intel Mac
            curl -L -o openjdk-17.tar.gz "https://github.com/adoptium/temurin17-binaries/releases/download/jdk-17.0.9%2B9/OpenJDK17U-jdk_x64_mac_hotspot_17.0.9_9.tar.gz"
        fi
        
        # Extract JDK
        tar -xzf openjdk-17.tar.gz
        
        # Create Java directory and move JDK
        sudo mkdir -p /Library/Java/JavaVirtualMachines/
        sudo mv jdk-17.0.9+9/Contents /Library/Java/JavaVirtualMachines/adoptopenjdk-17.jdk/
        
        # Set up environment variables
        append_to_zshrc 'export JAVA_HOME="/Library/Java/JavaVirtualMachines/adoptopenjdk-17.jdk/Contents/Home"'
        append_to_zshrc 'export PATH="$JAVA_HOME/bin:$PATH"'
        export JAVA_HOME="/Library/Java/JavaVirtualMachines/adoptopenjdk-17.jdk/Contents/Home"
        export PATH="$JAVA_HOME/bin:$PATH"
        
        # Clean up
        cd - > /dev/null
        rm -rf "$TEMP_DIR"
        
        print_success "Java installed successfully without Homebrew"
    fi
}

# Check and install Python if needed
if ! check_python; then
    install_python_macos
fi

# Check and install Java if needed
if ! check_java; then
    install_java_macos
fi

# Check if conda/mamba is installed
if ! command -v conda &> /dev/null && ! command -v mamba &> /dev/null; then
    echo "❌ Conda/Mamba not found. Installing miniconda ..."
    echo ""
    install_miniconda
fi

# Check if Nextflow is installed
if ! command -v nextflow &> /dev/null; then
    echo "📦 Installing Nextflow..."
    curl -s https://get.nextflow.io | bash
    sudo mv nextflow /usr/local/bin/
    echo "✅ Nextflow installed"
else
    echo "✅ Nextflow already installed: $(nextflow -version | head -1)"
fi

# Create test data
echo "📄 Creating test FASTA files..."
mkdir -p test_data

# Create sample FASTA file 1
cat > test_data/sample1.fasta << 'EOF'
>sequence1
ATCGATCGATCGATCG
>sequence2
GCTAGCTAGCTAGCTA
>sequence3
TTTTAAAACCCCGGGG
EOF

# Create sample FASTA file 2
cat > test_data/sample2.fasta << 'EOF'
>seq_A
AAAAAAAAAA
>seq_B
TTTTTTTTTT
>seq_C
CCCCCCCCCC
>seq_D
GGGGGGGGGG
>seq_E
ATCGATCGAT
EOF

echo "✅ Test data created in test_data/"

# Create results directory
mkdir -p results

echo ""
echo "🚀 Setup complete! Now you can run the pipeline:"
echo ""
echo "   # Run minimal example"
echo "   ./run.sh"