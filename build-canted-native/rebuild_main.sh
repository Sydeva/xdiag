#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
BUILD_DIR="$PROJECT_ROOT/build-canted-native"
SOURCE_DIR="$PROJECT_ROOT/cantedAFM"
BUILD_LOG="$SCRIPT_DIR/rebuild_main.build.log"
CONFIGURE_LOG="$SCRIPT_DIR/rebuild_main.configure.log"
VERIFY_LOG="$SCRIPT_DIR/rebuild_main.verify.log"
STATUS_FILE="$SCRIPT_DIR/rebuild_main.status"
ARTIFACT_FILE="$SCRIPT_DIR/rebuild_main.artifacts"
EXECUTABLE="$BUILD_DIR/main"
XDIAG_LIB="$PROJECT_ROOT/install-native/lib/libxdiag.a"
LINK_FILE="$BUILD_DIR/CMakeFiles/main.dir/link.txt"

APP_SELECTOR="${1:-GSwPG}"
GSWPG_ENTRY_SELECTOR="${2:-main}"
case "$APP_SELECTOR" in
  GSwPG|gswpg)
    APP_SELECTOR="GSwPG"
    ;;
  strucfac|STRUCFAC)
    APP_SELECTOR="strucfac"
    ;;
  *)
    echo "Usage: $0 [GSwPG|strucfac] [main|mainani]" | tee "$STATUS_FILE"
    exit 1
    ;;
esac

case "$GSWPG_ENTRY_SELECTOR" in
  main|mainani)
    ;;
  *)
    echo "Usage: $0 [GSwPG|strucfac] [main|mainani]" | tee "$STATUS_FILE"
    exit 1
    ;;
esac

if [[ "$APP_SELECTOR" != "GSwPG" && "$GSWPG_ENTRY_SELECTOR" != "main" ]]; then
  echo "Entry selector '$GSWPG_ENTRY_SELECTOR' is only valid with GSwPG" | tee "$STATUS_FILE"
  exit 1
fi

OBJECT_BASENAME="main.cpp"
if [[ "$APP_SELECTOR" == "GSwPG" ]]; then
  OBJECT_BASENAME="${GSWPG_ENTRY_SELECTOR}.cpp"
fi

OBJECT_FILE="$BUILD_DIR/CMakeFiles/main.dir/${APP_SELECTOR}/${OBJECT_BASENAME}.o"

rm -f "$STATUS_FILE"

if [[ ! -d "$BUILD_DIR" ]]; then
  echo "Missing build directory: $BUILD_DIR" | tee "$STATUS_FILE"
  exit 1
fi

if [[ ! -f "$BUILD_DIR/CMakeCache.txt" ]]; then
  echo "Missing CMake cache: $BUILD_DIR/CMakeCache.txt" | tee "$STATUS_FILE"
  exit 1
fi

cmake -S "$SOURCE_DIR" -B "$BUILD_DIR" -DCANTED_AFM_APP="$APP_SELECTOR" \
  -DCANTED_AFM_GSWPG_ENTRY="$GSWPG_ENTRY_SELECTOR" \
  >"$CONFIGURE_LOG" 2>&1

cmake --build "$BUILD_DIR" --target main --clean-first --parallel "$(nproc)" \
  >"$BUILD_LOG" 2>&1

if [[ ! -x "$EXECUTABLE" ]]; then
  echo "Missing executable after build: $EXECUTABLE" | tee "$STATUS_FILE"
  exit 1
fi

if [[ ! -f "$OBJECT_FILE" ]]; then
  echo "Missing object file after build: $OBJECT_FILE" | tee "$STATUS_FILE"
  exit 1
fi

if [[ ! -f "$XDIAG_LIB" ]]; then
  echo "Missing installed xdiag library: $XDIAG_LIB" | tee "$STATUS_FILE"
  exit 1
fi

{
  echo "Selected app: $APP_SELECTOR"
  echo "GSwPG entry: $GSWPG_ENTRY_SELECTOR"
  echo "Executable: $EXECUTABLE"
  echo "Object: $OBJECT_FILE"
  echo "xdiag library: $XDIAG_LIB"
  echo "Configure log: $CONFIGURE_LOG"
  echo "Link line:"
  cat "$LINK_FILE"
  echo
  echo "Timestamps:"
  stat -c '%Y %n' "$EXECUTABLE" "$OBJECT_FILE" "$XDIAG_LIB"
} >"$VERIFY_LOG"

stat -c '%Y %n' "$EXECUTABLE" "$OBJECT_FILE" "$XDIAG_LIB" >"$ARTIFACT_FILE"

echo "OK" > "$STATUS_FILE"
echo "Rebuilt main in $BUILD_DIR using $APP_SELECTOR/$OBJECT_BASENAME"
echo "Configure log: $CONFIGURE_LOG"
echo "Build log: $BUILD_LOG"
echo "Verify log: $VERIFY_LOG"
echo "Artifacts: $ARTIFACT_FILE"

