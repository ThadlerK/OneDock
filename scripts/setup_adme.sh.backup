#!/bin/bash

# Setup script for ADME Integration
# This script installs the necessary dependencies for the ADME screening module

echo "========================================"
echo "OneDock ADME Integration - Setup"
echo "========================================"

# Check if pip is available
if ! command -v pip &> /dev/null; then
    echo "❌ Error: pip not found. Please install Python 3 and pip first."
    exit 1
fi

echo ""
echo "📦 Installing Python dependencies..."
echo ""

# Install core dependencies needed for ADME
pip install --no-cache-dir \
    requests>=2.31.0 \
    beautifulsoup4>=4.12.0 \
    lxml>=4.9.0 \
    pandas>=2.0.0 \
    plotly>=5.14.0

echo ""
echo "✅ Dependencies installed!"
echo ""

# Test the installation
echo "🧪 Running validation tests..."
echo ""

python tests/test_adme_integration.py

if [ $? -eq 0 ]; then
    echo ""
    echo "========================================"
    echo "✅ Setup Complete!"
    echo "========================================"
    echo ""
    echo "You can now use the ADME screening feature."
    echo ""
    echo "To start the application:"
    echo "  streamlit run app/Home.py --server.address=0.0.0.0 --server.port=8501"
    echo ""
    echo "For documentation, see:"
    echo "  docs/ADME_INTEGRATION.md"
    echo ""
else
    echo ""
    echo "⚠️  Some tests failed. Please check the output above."
    echo ""
fi
