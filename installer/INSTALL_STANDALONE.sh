#!/usr/bin/env bash
set -Eeuo pipefail

###############################################################################
# INSTALL_STANDALONE.sh - Combustion Toolbox Standalone Installer
#
# Automatically detects the OS, downloads the latest standalone installer
# from GitHub, extracts it, and runs the appropriate installation.
#
# Supported platforms: macOS, Linux, Windows
#   
# Usage:
#   sh INSTALL_STANDALONE.sh
#
# @author: Alberto Cuadra Lara
#          Postdoctoral researcher - Group Fluid Mechanics
#          Universidad Carlos III de Madrid
###############################################################################

#-----------------------------#
#   GLOBAL VARIABLES & LOGS   #
#-----------------------------#

# Colors for logging
GREEN='\033[0;32m'
RED='\033[0;31m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Temporary folders/files
INSTALLER_FILE=""
EXTRACT_DIR="/tmp/combustion_toolbox_install"

# Logging helpers using printf
log_info()    { printf "%b[INFO]%b %s\n"    "$CYAN"  "$NC" "$*"; }
log_success() { printf "%b[SUCCESS]%b %s\n" "$GREEN" "$NC" "$*"; }
log_error()   { printf "%b[ERROR]%b %s\n"   "$RED"   "$NC" "$*"; }

# Cleanup function called on exit (if needed)
cleanup() {
  # Optionally remove extracted folder or leftover files
  :
}
trap cleanup EXIT

#-----------------------------#
#      STEP 1: DETECT OS      #
#-----------------------------#
detect_os() {
  log_info "Detecting operating system..."

  if [[ "$OSTYPE" == "darwin"* ]]; then
      PLATFORM="macos"
  elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
      PLATFORM="linux"
  elif [[ "$OSTYPE" == "msys" ]] || [[ "$OSTYPE" == "cygwin" ]]; then
      PLATFORM="windows"
  else
      log_error "Unsupported operating system."
      exit 1
  fi

  log_success "Detected OS: $PLATFORM"
}

#-----------------------------#
#  STEP 2: FETCH LATEST INFO  #
#-----------------------------#
fetch_latest_release() {
  local repo="CombustionToolbox/combustion_toolbox"
  local api_url="https://api.github.com/repos/$repo/releases/latest"

  log_info "Fetching the latest release from GitHub..."
  RELEASE_DATA=$(curl -s "$api_url" || true)

  if [[ -z "$RELEASE_DATA" ]]; then
    log_error "Could not fetch release data. Check your internet connection or GitHub status."
    exit 1
  fi
}

#-----------------------------#
#  STEP 3: FIND INSTALLER URL #
#-----------------------------#
find_installer_url() {
  # Look for a .zip asset that matches our PLATFORM
  INSTALLER_URL=$(echo "$RELEASE_DATA" \
    | sed -n "s/.*\"browser_download_url\": \"\\([^\"]*${PLATFORM}[^\"]*\\.zip\\)\".*/\\1/p")

  if [[ -z "$INSTALLER_URL" ]]; then
    log_error "No suitable .zip installer found for your platform ($PLATFORM)."
    exit 1
  fi

  log_info "Found installer URL: $INSTALLER_URL"
}

#-----------------------------#
#  STEP 4: DOWNLOAD INSTALLER #
#-----------------------------#
download_installer() {
  local installer_name
  installer_name=$(basename "$INSTALLER_URL")
  INSTALLER_FILE="/tmp/$installer_name"

  log_info "Downloading installer to $INSTALLER_FILE ..."
  curl -L -o "$INSTALLER_FILE" "$INSTALLER_URL"

  if [[ ! -f "$INSTALLER_FILE" ]]; then
    log_error "Download failed. Installer file not found."
    exit 1
  fi

  log_success "Download complete: $INSTALLER_FILE"
}

#-----------------------------#
#  STEP 5: EXTRACT INSTALLER  #
#-----------------------------#
extract_installer() {
  log_info "Extracting installer into $EXTRACT_DIR ..."
  mkdir -p "$EXTRACT_DIR"
  unzip -o "$INSTALLER_FILE" -d "$EXTRACT_DIR" >/dev/null

  if [[ $? -ne 0 ]]; then
    log_error "Extraction failed."
    exit 1
  fi

  # Remove any __MACOSX metadata folder
  if [[ -d "$EXTRACT_DIR/__MACOSX" ]]; then
    rm -rf "$EXTRACT_DIR/__MACOSX"
  fi

  log_success "Extraction complete."
}

#-----------------------------#
#  STEP 6: RUN INSTALLER      #
#-----------------------------#
run_installer() {
  log_info "Starting the installation process..."

  if [[ "$PLATFORM" == "windows" ]]; then
    local exe
    exe=$(find "$EXTRACT_DIR" -iname "*.exe" | head -n 1 || true)
    if [[ -z "$exe" ]]; then
      log_error "No .exe file found after extraction."
      exit 1
    fi
    log_info "Running $exe ..."
    start "" "$exe"

  elif [[ "$PLATFORM" == "macos" ]]; then
    local app
    app=$(find "$EXTRACT_DIR" -maxdepth 2 -iname "*.app" | head -n 1 || true)
    if [[ -z "$app" ]]; then
      log_error "No .app found after extraction."
      exit 1
    fi

    log_info "Applying executable permissions to $app/Contents/MacOS ..."
    chmod -R +x "$app/Contents/MacOS"

    log_info "Opening application..."
    open "$app"

  elif [[ "$PLATFORM" == "linux" ]]; then
    local install_file
    install_file=$(find "$EXTRACT_DIR" -iname "*.install" | head -n 1 || true)
    if [[ -z "$install_file" ]]; then
      log_error "No .install file found after extraction."
      exit 1
    fi

    chmod +x "$install_file"
    log_info "Running $install_file ..."
    "$install_file" &

  else
    log_error "Unsupported platform."
    exit 1
  fi

  log_success "Installation process started. Follow the on-screen instructions."
}

#-----------------------------#
#           MAIN              #
#-----------------------------#
main() {
  detect_os
  fetch_latest_release
  find_installer_url
  download_installer
  extract_installer
  run_installer
}

main "$@"
