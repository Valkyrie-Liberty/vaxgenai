"""
Utility functions for the VaxGenAI platform
"""

import os
import logging
from pathlib import Path
from datetime import datetime

def setup_logging(log_level=logging.INFO):
    """Set up logging configuration"""
    log_dir = Path("logs")
    log_dir.mkdir(exist_ok=True)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = log_dir / f"vaxgenai_{timestamp}.log"
    
    # Configure logging
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    
    logger = logging.getLogger("vaxgenai")
    logger.info(f"Logging initialized. Log file: {log_file}")
    return logger

def get_project_root():
    """Return the project root directory"""
    return Path(__file__).parent.parent

def ensure_directory(directory_path):
    """Ensure a directory exists, create it if it doesn't"""
    Path(directory_path).mkdir(parents=True, exist_ok=True)
    return directory_path

def get_data_path():
    """Return the path to the data directory"""
    data_path = get_project_root() / "data"
    ensure_directory(data_path)
    return data_path

def get_results_path():
    """Return the path to the results directory"""
    results_path = get_project_root() / "results"
    ensure_directory(results_path)
    return results_path

def get_models_path():
    """Return the path to the models directory"""
    models_path = get_project_root() / "models"
    ensure_directory(models_path)
    return models_path

