#!/usr/bin/env python3
"""
Validate pipeline installation and dependencies

Run this after installing to check that everything is set up correctly.
"""

import sys
import os
import subprocess
from pathlib import Path


def check_python_version():
    """Check Python version"""
    print("=" * 60)
    print("Checking Python version...")
    version = sys.version_info
    if version.major == 3 and version.minor >= 8:
        print(f"✓ Python {version.major}.{version.minor}.{version.micro}")
        return True
    else:
        print(f"✗ Python {version.major}.{version.minor}.{version.micro} (need ≥3.8)")
        return False


def check_package(package_name, import_name=None):
    """Check if Python package is installed"""
    if import_name is None:
        import_name = package_name

    try:
        __import__(import_name)
        print(f"✓ {package_name}")
        return True
    except ImportError:
        print(f"✗ {package_name} (install with: pip install {package_name})")
        return False


def check_python_packages():
    """Check required Python packages"""
    print("\n" + "=" * 60)
    print("Checking Python packages...")

    required = {
        'torch': 'torch',
        'numpy': 'numpy',
        'biopython': 'Bio',
        'abnumber': 'abnumber',
        'iglm': 'iglm',
        'transformers': 'transformers',
        'pyyaml': 'yaml',
    }

    results = []
    for pkg, import_name in required.items():
        results.append(check_package(pkg, import_name))

    return all(results)


def check_optional_packages():
    """Check optional Python packages"""
    print("\n" + "=" * 60)
    print("Checking optional packages...")

    optional = {
        'chai_lab': 'chai_lab',
        'scipy': 'scipy',
        'pandas': 'pandas',
    }

    for pkg, import_name in optional.items():
        check_package(pkg, import_name)


def check_antibody_pipeline():
    """Check if antibody_pipeline is installed"""
    print("\n" + "=" * 60)
    print("Checking antibody_pipeline installation...")

    try:
        from antibody_pipeline import AntibodyDesignPipeline, PipelineConfig
        print("✓ antibody_pipeline installed")
        return True
    except ImportError:
        print("✗ antibody_pipeline not installed")
        print("  Install with: pip install -e .")
        return False


def check_rfdiffusion():
    """Check RFdiffusion installation"""
    print("\n" + "=" * 60)
    print("Checking RFdiffusion...")

    rfdiffusion_path = Path("./RFdiffusion")
    run_inference = rfdiffusion_path / "scripts" / "run_inference.py"

    if run_inference.exists():
        print(f"✓ RFdiffusion found at {rfdiffusion_path}")
        return True
    else:
        print(f"✗ RFdiffusion not found at {rfdiffusion_path}")
        print("  Clone with: git clone https://github.com/RosettaCommons/RFdiffusion.git")
        return False


def check_hdock():
    """Check HDock installation"""
    print("\n" + "=" * 60)
    print("Checking HDock (optional)...")

    hdock_path = Path("./bin/hdock")
    createpl_path = Path("./bin/createpl")

    if hdock_path.exists() and createpl_path.exists():
        print(f"✓ HDock found in bin/")

        # Check if executable
        if os.access(hdock_path, os.X_OK):
            print("✓ hdock is executable")
        else:
            print("✗ hdock is not executable (run: chmod +x bin/hdock)")
            return False

        if os.access(createpl_path, os.X_OK):
            print("✓ createpl is executable")
        else:
            print("✗ createpl is not executable (run: chmod +x bin/createpl)")
            return False

        return True
    else:
        print("✗ HDock not found in bin/")
        print("  Download from: http://huanglab.phys.hust.edu.cn/software/hdocklite/")
        print("  (Optional - only needed for Mode 1: separate antigen + framework)")
        return False


def check_directories():
    """Check directory structure"""
    print("\n" + "=" * 60)
    print("Checking directory structure...")

    required_dirs = [
        "antibody_pipeline",
        "examples",
        "docs",
        "data/input",
        "data/output",
        "data/models",
        "bin",
    ]

    results = []
    for dir_path in required_dirs:
        path = Path(dir_path)
        if path.exists():
            print(f"✓ {dir_path}/")
            results.append(True)
        else:
            print(f"✗ {dir_path}/ (creating...)")
            path.mkdir(parents=True, exist_ok=True)
            results.append(True)

    return all(results)


def check_cuda():
    """Check CUDA availability"""
    print("\n" + "=" * 60)
    print("Checking CUDA...")

    try:
        import torch
        if torch.cuda.is_available():
            print(f"✓ CUDA available (device: {torch.cuda.get_device_name(0)})")
            return True
        else:
            print("⚠ CUDA not available (will use CPU - slower)")
            return False
    except ImportError:
        print("✗ PyTorch not installed")
        return False


def main():
    """Run all checks"""
    print("=" * 60)
    print("ANTIBODY DESIGN PIPELINE - INSTALLATION VALIDATION")
    print("=" * 60)

    checks = []

    # Critical checks
    checks.append(("Python version", check_python_version()))
    checks.append(("Python packages", check_python_packages()))
    checks.append(("Pipeline package", check_antibody_pipeline()))
    checks.append(("RFdiffusion", check_rfdiffusion()))
    checks.append(("Directory structure", check_directories()))

    # Optional checks
    check_optional_packages()
    check_hdock()
    check_cuda()

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    critical_passed = all(result for name, result in checks)

    for name, result in checks:
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{status}: {name}")

    print("\n" + "=" * 60)
    if critical_passed:
        print("✓ Installation is complete and ready to use!")
        print("\nNext steps:")
        print("1. Review examples: ls examples/")
        print("2. Run example: python examples/run_example.py")
        print("3. Or create config: python -m antibody_pipeline --create-config config.yaml")
        return 0
    else:
        print("✗ Installation incomplete. Please fix the errors above.")
        print("\nFor help, see:")
        print("- README.md")
        print("- docs/QUICK_START.md")
        return 1


if __name__ == "__main__":
    sys.exit(main())