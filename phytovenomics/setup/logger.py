#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Logger Class

This class provides logging functionality for the Phytovenomics setup process,
with support for different log levels and colored output.
"""

import os
import sys
import logging
from datetime import datetime
from pathlib import Path


class Logger:
    """
    Provides logging functionality for the setup process.
    
    Attributes:
        log_path (str): Path to log file
        log_level (str): Log level (INFO, WARNING, ERROR, DEBUG)
        logger (logging.Logger): Logger instance
    """
    
    # ANSI color codes for colored terminal output
    COLORS = {
        'RESET': '\033[0m',
        'INFO': '\033[94m',  # Blue
        'SUCCESS': '\033[92m',  # Green
        'WARNING': '\033[93m',  # Yellow
        'ERROR': '\033[91m',  # Red
        'BOLD': '\033[1m'
    }
    
    def __init__(self, log_path=None, log_level='INFO'):
        """
        Initialize the Logger.
        
        Args:
            log_path (str): Path to log file, or None for no file logging
            log_level (str): Log level (INFO, WARNING, ERROR, DEBUG)
        """
        self.log_path = log_path
        self.log_level = log_level
        
        # Create logger
        self.logger = logging.getLogger('phytovenomics_setup')
        self.logger.setLevel(logging.DEBUG)  # Capture all logs
        
        # Clear existing handlers
        self.logger.handlers = []
        
        # Add console handler
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(getattr(logging, log_level))
        console_handler.setFormatter(logging.Formatter('%(message)s'))
        self.logger.addHandler(console_handler)
        
        # Add file handler if log path is provided
        if log_path:
            # Create log directory if it doesn't exist
            log_dir = os.path.dirname(log_path)
            if log_dir:
                Path(log_dir).mkdir(exist_ok=True, parents=True)
            
            file_handler = logging.FileHandler(log_path)
            file_handler.setLevel(logging.DEBUG)  # Log everything to file
            file_formatter = logging.Formatter(
                '%(asctime)s - %(levelname)s - %(message)s',
                datefmt='%Y-%m-%d %H:%M:%S'
            )
            file_handler.setFormatter(file_formatter)
            self.logger.addHandler(file_handler)
    
    def log_info(self, message):
        """
        Log an informational message.
        
        Args:
            message (str): Message to log
        """
        colored_message = f"{self.COLORS['INFO']}[INFO] {message}{self.COLORS['RESET']}"
        self.logger.info(message)
        if self.log_level != 'INFO':
            print(colored_message)
    
    def log_success(self, message):
        """
        Log a success message.
        
        Args:
            message (str): Message to log
        """
        colored_message = f"{self.COLORS['SUCCESS']}[SUCCESS] {message}{self.COLORS['RESET']}"
        self.logger.info(f"[SUCCESS] {message}")
        print(colored_message)
    
    def log_warning(self, message):
        """
        Log a warning message.
        
        Args:
            message (str): Message to log
        """
        colored_message = f"{self.COLORS['WARNING']}[WARNING] {message}{self.COLORS['RESET']}"
        self.logger.warning(message)
        print(colored_message)
    
    def log_error(self, message):
        """
        Log an error message.
        
        Args:
            message (str): Message to log
        """
        colored_message = f"{self.COLORS['ERROR']}[ERROR] {message}{self.COLORS['RESET']}"
        self.logger.error(message)
        print(colored_message)
    
    def log_debug(self, message):
        """
        Log a debug message.
        
        Args:
            message (str): Message to log
        """
        self.logger.debug(message)
        if self.log_level == 'DEBUG':
            print(f"[DEBUG] {message}")
