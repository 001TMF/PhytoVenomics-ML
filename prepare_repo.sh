#!/bin/bash

echo "Preparing Phytovenomics repository for GitHub..."

# Make sure we're in the project root
if [ ! -f "setup.py" ] || [ ! -d "phytovenomics" ]; then
  echo "Error: Run this script from the project root directory"
  exit 1
fi

# Clean pycache files
echo "Cleaning Python cache files..."
find . -type d -name "__pycache__" -exec rm -rf {} +
find . -name "*.pyc" -delete
find . -name "*.pyo" -delete
find . -name "*.pyd" -delete
find . -name ".pytest_cache" -exec rm -rf {} +
find . -name ".coverage" -delete
find . -name "htmlcov" -exec rm -rf {} +

# Check for large files
echo "Checking for large files (>10MB)..."
large_files=$(find . -type f -size +10M | grep -v ".git" | sort -h)
if [ ! -z "$large_files" ]; then
  echo "Warning: Large files found that might not be suitable for GitHub:"
  echo "$large_files"
  echo ""
  echo "Consider adding these to .gitignore if they shouldn't be tracked."
fi

# Check for sensitive data patterns
echo "Checking for potential API keys or credentials..."
sensitive_files=$(grep -r -l "api[_-]key\|password\|secret\|credential" --include="*.py" --include="*.ipynb" --include="*.json" --include="*.yaml" --include="*.yml" .)
if [ ! -z "$sensitive_files" ]; then
  echo "Warning: Files containing potential sensitive information found:"
  echo "$sensitive_files"
  echo ""
  echo "Review these files to ensure no credentials are being committed."
fi

# Verify .gitignore exists
if [ ! -f ".gitignore" ]; then
  echo "Warning: No .gitignore file found. Creating a basic one..."
  cat > .gitignore << 'GITIGNORE'
# Python
__pycache__/
*.py[cod]
*$py.class
*.so
.Python
build/
develop-eggs/
dist/
downloads/
eggs/
.eggs/
lib/
lib64/
parts/
sdist/
var/
wheels/
*.egg-info/
.installed.cfg
*.egg

# Testing
.coverage
htmlcov/
.pytest_cache/
.tox/

# Environments
.env
.venv
env/
venv/
ENV/
env.bak/
venv.bak/

# IDE
.idea/
.vscode/
*.swp

# OS specific
.DS_Store

# Large data files
*.h5
*.hdf5
*.pkl
*.npy
*.npz
GITIGNORE
fi

# Make script executable
chmod +x scripts/*.py
chmod +x install.sh
chmod +x setup_igfold.sh

echo "Repository preparation complete!"
echo "Next steps:"
echo "  1. Review the warnings above (if any)"
echo "  2. Run 'git init' if this is a new repository"
echo "  3. Run 'git add .' to stage files"
echo "  4. Run 'git commit -m "Initial commit"' to commit files"
echo "  5. Add your remote with 'git remote add origin <repo-url>'"
echo "  6. Push with 'git push -u origin main'"
