#!/bin/bash

# Nextflow Installation Script for macOS and Linux
# This script installs Java (if needed) and Nextflow

set -e  # Exit on any error

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

# Function to check Miniconda installation
check_miniconda() {
    if [ -f "$HOME/miniconda3/bin/conda" ]; then
        print_success "Miniconda is already installed"
        return 0
    else
        print_status "Miniconda is not installed"
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
        echo 'export PATH="/usr/local/bin:$PATH"' >> ~/.zshrc
        echo 'export PATH="/Library/Frameworks/Python.framework/Versions/3.11/bin:$PATH"' >> ~/.zshrc
        export PATH="/Library/Frameworks/Python.framework/Versions/3.11/bin:$PATH"
        
        # Clean up
        cd - > /dev/null
        rm -rf "$TEMP_DIR"
        
        print_success "Python3 installed successfully"
    fi
}

# Function to install Python on Linux
install_python_linux() {
    print_status "Installing Python3 on Linux..."
    
    # Detect Linux distribution
    if [ -f /etc/debian_version ]; then
        # Debian/Ubuntu
        print_status "Detected Debian/Ubuntu system"
        sudo apt update
        sudo apt install -y python3 python3-pip python3-venv
        
    elif [ -f /etc/redhat-release ]; then
        # RedHat/CentOS/Fedora
        print_status "Detected RedHat/CentOS/Fedora system"
        if command_exists dnf; then
            sudo dnf install -y python3 python3-pip
        elif command_exists yum; then
            sudo yum install -y python3 python3-pip
        fi
        
    elif [ -f /etc/arch-release ]; then
        # Arch Linux
        print_status "Detected Arch Linux system"
        sudo pacman -S --noconfirm python python-pip
        
    else
        print_error "Unsupported Linux distribution for automatic Python installation"
        print_error "Please install Python3 manually and run this script again"
        exit 1
    fi
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

# Function to install Java on macOS
install_java_macos() {
    print_status "Installing Java on macOS..."
    
    if command_exists brew; then
        print_status "Using Homebrew to install OpenJDK..."
        brew install openjdk@17
        
        # Detect if using Apple Silicon or Intel Mac for correct path
        if [[ $(uname -m) == "arm64" ]]; then
            # Apple Silicon Mac
            echo 'export PATH="/opt/homebrew/bin:$PATH"' >> ~/.zshrc
            echo 'export JAVA_HOME="/opt/homebrew/opt/openjdk@17"' >> ~/.zshrc
            export PATH="/opt/homebrew/bin:$PATH"
            export JAVA_HOME="/opt/homebrew/opt/openjdk@17"
        else
            # Intel Mac
            echo 'export PATH="/usr/local/bin:$PATH"' >> ~/.zshrc
            echo 'export JAVA_HOME="/usr/local/opt/openjdk@17"' >> ~/.zshrc
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
        echo 'export JAVA_HOME="/Library/Java/JavaVirtualMachines/adoptopenjdk-17.jdk/Contents/Home"' >> ~/.zshrc
        echo 'export PATH="$JAVA_HOME/bin:$PATH"' >> ~/.zshrc
        export JAVA_HOME="/Library/Java/JavaVirtualMachines/adoptopenjdk-17.jdk/Contents/Home"
        export PATH="$JAVA_HOME/bin:$PATH"
        
        # Clean up
        cd - > /dev/null
        rm -rf "$TEMP_DIR"
        
        print_success "Java installed successfully without Homebrew"
    fi
}

# Function to install Java on Linux
install_java_linux() {
    print_status "Installing Java on Linux..."
    
    # Detect Linux distribution
    if [ -f /etc/debian_version ]; then
        # Debian/Ubuntu
        print_status "Detected Debian/Ubuntu system"
        sudo apt update
        sudo apt install -y openjdk-17-jdk
        
    elif [ -f /etc/redhat-release ] || [ -f /etc/system-release ]; then
        # RedHat/CentOS/Fedora/Amazon Linux
        print_status "Detected RedHat/CentOS/Fedora/Amazon Linux system"
        if command_exists dnf; then
            sudo dnf install -y java
        elif command_exists yum; then
            # For Amazon Linux and other yum-based systems
            sudo yum update -y
            sudo yum install -y java
        fi
        
    elif [ -f /etc/arch-release ]; then
        # Arch Linux
        print_status "Detected Arch Linux system"
        sudo pacman -S --noconfirm jdk17-openjdk
        
    else
        print_error "Unsupported Linux distribution"
        print_error "Please install Java 17 manually and run this script again"
        exit 1
    fi
    
    # Verify Java installation
    if command_exists java; then
        JAVA_VERSION=$(java -version 2>&1 | head -n1 | cut -d'"' -f2 | cut -d'.' -f1)
        if [ "$JAVA_VERSION" -ge 11 ]; then
            print_success "Java $JAVA_VERSION installed successfully"
        else
            print_error "Java installation failed - version $JAVA_VERSION is less than required version 11"
            exit 1
        fi
    else
        print_error "Java installation failed - java command not found"
        exit 1
    fi
}

# Function to install Nextflow
install_nextflow() {
    print_status "Installing Nextflow..."
    
    # Create bin directory if it doesn't exist
    mkdir -p ~/bin
    
    # Download Nextflow
    curl -s https://get.nextflow.io | bash
    
    # Move to bin directory
    mv nextflow ~/bin/
    
    # Make executable
    chmod +x ~/bin/nextflow
    
    # Add to PATH if not already there
    if [[ ":$PATH:" != *":$HOME/bin:"* ]]; then
        echo 'export PATH="$HOME/bin:$PATH"' >> ~/.bashrc
        echo 'export PATH="$HOME/bin:$PATH"' >> ~/.zshrc
        export PATH="$HOME/bin:$PATH"
    fi
}


## INstall make on amazon linux

# Function to check Docker installation
check_docker() {
    if command_exists docker; then
        print_success "Docker is already installed"
        return 0
    else
        print_status "Docker is not installed"
        return 1
    fi
}

# Function to install Docker on Linux
install_docker_linux() {
    print_status "Installing Docker on Linux..."
    
    # Detect Linux distribution
    if [ -f /etc/debian_version ]; then
        # Debian/Ubuntu
        print_status "Detected Debian/Ubuntu system"
        
        # Update package index
        sudo apt update
        
        # Install prerequisites
        sudo apt install -y apt-transport-https ca-certificates curl software-properties-common
        
        # Add Docker's official GPG key
        curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
        
        # Add Docker repository
        sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable"
        
        # Update package index again
        sudo apt update
        
        # Install Docker
        sudo apt install -y docker-ce docker-ce-cli containerd.io
        
    elif [ -f /etc/redhat-release ] || [ -f /etc/system-release ]; then
        # RedHat/CentOS/Fedora/Amazon Linux
        print_status "Detected RedHat/CentOS/Fedora/Amazon Linux system"
        
        if command_exists dnf; then
            # Fedora
            sudo dnf -y install dnf-plugins-core
            sudo dnf config-manager --add-repo https://download.docker.com/linux/fedora/docker-ce.repo
            sudo dnf install -y docker-ce docker-ce-cli containerd.io
        elif command_exists yum; then
            # Amazon Linux/CentOS/RHEL
            sudo yum update -y
            sudo yum install -y docker
            sudo systemctl start docker
            sudo systemctl enable docker
        fi
        
    elif [ -f /etc/arch-release ]; then
        # Arch Linux
        print_status "Detected Arch Linux system"
        sudo pacman -S --noconfirm docker
        
    else
        print_error "Unsupported Linux distribution for automatic Docker installation"
        print_error "Please install Docker manually and run this script again"
        exit 1
    fi
    
    # Start and enable Docker service
    sudo systemctl start docker
    sudo systemctl enable docker
    
    # Add current user to docker group
    sudo usermod -aG docker $USER
    
    # Verify Docker installation
    if command_exists docker; then
        DOCKER_VERSION=$(docker --version)
        print_success "Docker installed successfully: $DOCKER_VERSION"
    else
        print_error "Docker installation failed - docker command not found"
        exit 1
    fi
    
    print_warning "You may need to log out and log back in for the docker group changes to take effect"
}

# Function to check if make is installed
check_make() {
    if command_exists make; then
        print_success "Make is already installed"
        return 0
    else
        print_status "Make is not installed"
        return 1
    fi
}

# Function to install make on Linux
install_make_linux() {
    print_status "Installing Make on Linux..."
    
    # Detect Linux distribution
    if [ -f /etc/debian_version ]; then
        # Debian/Ubuntu
        print_status "Detected Debian/Ubuntu system"
        sudo apt update
        sudo apt install -y make
        
    elif [ -f /etc/redhat-release ] || [ -f /etc/system-release ]; then
        # RedHat/CentOS/Fedora/Amazon Linux
        print_status "Detected RedHat/CentOS/Fedora/Amazon Linux system"
        if command_exists dnf; then
            sudo dnf install -y make
        elif command_exists yum; then
            sudo yum install -y make
        fi
        
    elif [ -f /etc/arch-release ]; then
        # Arch Linux
        print_status "Detected Arch Linux system"
        sudo pacman -S --noconfirm make
        
    else
        print_error "Unsupported Linux distribution for automatic make installation"
        print_error "Please install make manually and run this script again"
        exit 1
    fi
    
    print_success "Make installed successfully"
}

# Function to verify installation
verify_installation() {
    print_status "Verifying installation..."
    
    # Source the shell configuration
    if [ -f ~/.zshrc ]; then
        source ~/.zshrc
    elif [ -f ~/.bashrc ]; then
        source ~/.bashrc
    fi
    
    # Check Java
    if check_java; then
        print_success "Java verification passed"
    else
        print_error "Java verification failed"
        return 1
    fi
    
    # Check Python
    if check_python; then
        print_success "Python3 verification passed"
    else
        print_warning "Python3 verification failed - may need terminal restart"
    fi
    
    # Check Miniconda
    if check_miniconda; then
        print_success "Miniconda verification passed"
    else
        print_warning "Miniconda verification failed - may need terminal restart"
    fi
    
    # Check Docker on Linux
    if [ "$(uname -s)" = "Linux" ]; then
        if check_docker; then
            print_success "Docker verification passed"
        else
            print_warning "Docker verification failed - may need terminal restart"
        fi
    fi
    
    # Check Nextflow
    if command_exists nextflow; then
        NF_VERSION=$(nextflow -version 2>&1 | head -n1)
        print_success "Nextflow installed: $NF_VERSION"
    else
        print_error "Nextflow installation failed"
        return 1
    fi
}

# Main installation function
main() {
    print_status "Starting installation of Nextflow, Java, Python3, Miniconda, and Docker..."
    print_status "Detected OS: $(uname -s)"
    
    # Detect operating system
    OS=$(uname -s)
    
    case $OS in
        "Darwin")
            print_status "macOS detected"
            
            # Check and install Java if needed
            if ! check_java; then
                install_java_macos
            fi
            
            # Check and install Python if needed
            if ! check_python; then
                install_python_macos
            fi
            ;;
            
        "Linux")
            print_status "Linux detected"
            
            # Check and install make if needed
            if ! check_make; then
                install_make_linux
            fi
            
            # Check and install Java if needed
            if ! check_java; then
                install_java_linux
            fi
            
            # Check and install Python if needed
            if ! check_python; then
                install_python_linux
            fi
            
            # Check and install Docker if needed
            if ! check_docker; then
                install_docker_linux
            fi
            ;;
            
        *)
            print_error "Unsupported operating system: $OS"
            print_error "This script supports macOS and Linux only"
            exit 1
            ;;
    esac
    
    # Check and install Miniconda if needed
    if ! check_miniconda; then
        install_miniconda
    fi
    
    # Install Nextflow
    if command_exists nextflow; then
        print_success "Nextflow is already installed"
        nextflow -version
    else
        install_nextflow
    fi
    
    # Verify installation
    verify_installation
    
    print_success "Installation completed successfully!"
    print_status "Installed components:"
    print_status "  - Java (OpenJDK 17)"
    print_status "  - Python3"
    print_status "  - Miniconda"
    if [ "$(uname -s)" = "Linux" ]; then
        print_status "  - Docker"
    fi
    print_status "  - Nextflow"
    print_status ""
    print_status "You may need to restart your terminal or run 'source ~/.bashrc' (or ~/.zshrc) to use all tools"
    print_status "Test your installations with:"
    print_status "  - nextflow run hello"
    print_status "  - python3 --version"
    print_status "  - conda --version"
    if [ "$(uname -s)" = "Linux" ]; then
        print_status "  - docker --version"
    fi
}

# Run main function
main "$@"
