# Binary Executables Directory

This directory contains external binary executables required for the pipeline.

## Required Binaries

### HDock (Optional - for docking)

**Required for**: Mode 1 (separate antigen + framework with docking)

**Files needed**:
- `hdock` - Main docking executable
- `createpl` - Post-processing tool

**Download from**: http://huanglab.phys.hust.edu.cn/software/hdocklite/

**Installation**:
```bash
# 1. Download HDock from website
# 2. Extract the archive
# 3. Copy binaries to this directory
cp /path/to/hdock bin/hdock
cp /path/to/createpl bin/createpl

# 4. Make executable
chmod +x bin/hdock
chmod +x bin/createpl
```

**Usage in config**:
```yaml
docking:
  enabled: true
  hdock_bin: ./bin/hdock
  createpl_bin: ./bin/createpl
```

## Notes

- These binaries are **not** included in the repository
- HDock is only needed if using docking (Mode 1)
- For pre-docked complexes (Mode 2) or de novo design (Mode 3), HDock is not required
- Make sure binaries are executable (`chmod +x`)

## Platform Compatibility

HDock provides binaries for:
- Linux (x86_64)
- macOS (may need Rosetta 2 on Apple Silicon)

If binaries don't work on your system, you may need to compile from source or use a compatible platform.

## Verification

Test that HDock is working:
```bash
./bin/hdock -h
./bin/createpl -h
```

You should see help messages if installed correctly.