# build-canted-native

This directory holds the native cantedAFM build tree plus helper scripts/logs for rebuilding `xdiag` and the `cantedAFM` `main` executable.

## What to select

- `CANTED_AFM_APP`: `GSwPG` or `strucfac`
- `CANTED_AFM_GSWPG_ENTRY`: `main` or `mainani` (only valid when `CANTED_AFM_APP=GSwPG`)

Valid combinations:

- `GSwPG` + `main` -> builds `cantedAFM/GSwPG/main.cpp`
- `GSwPG` + `mainani` -> builds `cantedAFM/GSwPG/mainani.cpp`
- `strucfac` + `main` -> builds `cantedAFM/strucfac/main.cpp`

## Recommended rebuild flow

1. Rebuild/install native `xdiag` (updates `install-native/lib/libxdiag.a`).
2. Rebuild `cantedAFM/main` with the app/entry selectors.

```bash
bash /home/t30/all/ge45hub/CLionProjects/xdiagAFM/build-canted-native/rebuild_xdiag.sh
bash /home/t30/all/ge45hub/CLionProjects/xdiagAFM/build-canted-native/rebuild_main.sh GSwPG mainani
```

## `rebuild_main.sh` usage

```bash
bash /home/t30/all/ge45hub/CLionProjects/xdiagAFM/build-canted-native/rebuild_main.sh [GSwPG|strucfac] [main|mainani]
```

Examples:

```bash
bash /home/t30/all/ge45hub/CLionProjects/xdiagAFM/build-canted-native/rebuild_main.sh GSwPG
bash /home/t30/all/ge45hub/CLionProjects/xdiagAFM/build-canted-native/rebuild_main.sh GSwPG mainani
bash /home/t30/all/ge45hub/CLionProjects/xdiagAFM/build-canted-native/rebuild_main.sh strucfac main
```

`strucfac mainani` is rejected by design.

## Direct CMake equivalent

```bash
cmake -S /home/t30/all/ge45hub/CLionProjects/xdiagAFM/cantedAFM -B /home/t30/all/ge45hub/CLionProjects/xdiagAFM/build-canted-native -DCANTED_AFM_APP=GSwPG -DCANTED_AFM_GSWPG_ENTRY=mainani
cmake --build /home/t30/all/ge45hub/CLionProjects/xdiagAFM/build-canted-native --target main --clean-first --parallel "$(nproc)"
```

## Logs and status files

- `rebuild_main.status`: last status (`OK` or error message)
- `rebuild_main.configure.log`: CMake configure output for `cantedAFM`
- `rebuild_main.build.log`: build output for `main`
- `rebuild_main.verify.log`: selected options, link line, and timestamps
- `rebuild_main.artifacts`: timestamp snapshot for executable/object/lib
- `rebuild_xdiag.status`: last status for native library rebuild
- `rebuild_xdiag.build.log`: build output for `xdiag`
- `rebuild_xdiag.install.log`: install output for native prefix
- `rebuild_xdiag.artifacts`: timestamp snapshot for built/installed `libxdiag.a`

