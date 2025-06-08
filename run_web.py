"""
Web interface runner for VaxGenAI
"""

import os
import sys
import logging
from pathlib import Path

# Add the project root to the path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.web_interface import WebInterface
from src.utils import setup_logging

def main():
    """Run the web interface"""
    # Setup logging
    logger = setup_logging()
    logger.info("Starting VaxGenAI web interface")
    
    # Create web interface
    web_interface = WebInterface()
    
    # Run web interface
    web_interface.run(host="0.0.0.0", port=8000)

if __name__ == "__main__":
    main()

