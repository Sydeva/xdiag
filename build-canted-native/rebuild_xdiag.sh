#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
BUILD_NATIVE="$PROJECT_ROOT/build-native"
INSTALL_NATIVE="$PROJECT_ROOT/install-native"
BUILD_LOG="$SCRIPT_DIR/rebuild_xdiag.build.log"
INSTALL_LOG="$SCRIPT_DIR/rebuild_xdiag.install.log"
STATUS_FILE="$SCRIPT_DIR/rebuild_xdiag.status"
TIMESTAMP_FILE="$SCRIPT_DIR/rebuild_xdiag.artifacts"

rm -f "$STATUS_FILE"

if [[ ! -d "$BUILD_NATIVE" ]]; then
  echo "Missing build directory: $BUILD_NATIVE" | tee "$STATUS_FILE"
  exit 1
fi

if [[ ! -f "$BUILD_NATIVE/CMakeCache.txt" ]]; then
  echo "Missing CMake cache: $BUILD_NATIVE/CMakeCache.txt" | tee "$STATUS_FILE"
  exit 1
fi

cmake --build "$BUILD_NATIVE" --target xdiag --clean-first --parallel "$(nproc)" \
  >"$BUILD_LOG" 2>&1

cmake --install "$BUILD_NATIVE" >"$INSTALL_LOG" 2>&1

stat -c '%Y %n' \
  "$PROJECT_ROOT/build-native/xdiag/libxdiag.a" \
  "$INSTALL_NATIVE/lib/libxdiag.a" \
  > "$TIMESTAMP_FILE"

echo "OK" > "$STATUS_FILE"
echo "Rebuilt xdiag from $BUILD_NATIVE and installed into $INSTALL_NATIVE"
echo "Build log: $BUILD_LOG"
echo "Install log: $INSTALL_LOG"
echo "Artifacts: $TIMESTAMP_FILE"

