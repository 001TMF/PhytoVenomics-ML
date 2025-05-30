#!/usr/bin/env python
# antibody_design/cdr_processor.py

import re
import logging
from typing import Dict, List, Optional, Tuple, Union

logger = logging.getLogger("Phytovenomics.CDRProcessor")

class CDRProcessor:
    """
    Module for identifying and manipulating CDR regions in antibody sequences.
    
    This class provides functionality to:
    1. Extract CDR sequences from antibody sequences
    2. Mask CDRs with placeholder tokens
    3. Replace CDRs with new sequences
    4. Identify framework regions
    """
    
    def __init__(self):
        """Initialize the CDR processor with pattern definitions."""
        # Define CDR patterns based on Kabat numbering scheme
        # These patterns are simplified and work for most common antibodies
        self.cdr_patterns = {
            # Heavy chain CDRs
            "CDRH1": r"(?:C[A-Z]{1,2}[KRQAG][A-Z]{5,7}G[YF][A-Z]{1,2})([A-Z]{10,12})(?:W[VIM][KRQAG][A-Z][A-Z])",  # ~10-12 residues
            "CDRH2": r"(?:[WY][VIM][KRQAG][A-Z])([A-Z]{16,19})(?:[KRQAG][LIVFTA][TSAGN][A-Z])",  # ~16-19 residues
            "CDRH3": r"(?:C[A-Z]{1,2}[KRQAG])([A-Z]{3,25})(?:W[GY][A-Z][A-Z])",  # Variable length (3-25 residues)
            
            # Light chain CDRs
            "CDRL1": r"(?:C[A-Z]{2}[A-Z][A-Z][A-Z][A-Z]{1,2})([A-Z]{10,17})(?:[A-Z][A-Z][LIVMFTA][A-Z])",  # ~10-17 residues
            "CDRL2": r"(?:[A-Z][A-Z][LIVMFTA][A-Z])([A-Z]{7})(?:[TSAGN][A-Z][A-Z])",  # ~7 residues
            "CDRL3": r"(?:C[A-Z]{2}[A-Z][A-Z][A-Z])([A-Z]{7,11})(?:[FW]G[A-Z]G[A-Z][A-Z][A-Z][A-Z])"  # ~7-11 residues
        }
        
        # Common structural markers for antibodies (helps identify VH/VL)
        self.vh_markers = ["WVRQAPGKGLEWV", "WVRQAPGQGLEWM", "WVRQPPGKGLEWIG"]
        self.vl_markers = ["WYQQKPGKAPKLL", "WYQQKPGQAPRLL", "WYQQKPGKAPKLLI"]
        
        logger.info("CDR processor initialized with standard Kabat-based patterns")
    
    def extract_cdr_sequences(self, antibody_sequence: str) -> Dict[str, str]:
        """
        Extract CDR sequences from an antibody sequence.
        
        Args:
            antibody_sequence: Amino acid sequence of the antibody
            
        Returns:
            Dictionary mapping CDR names to their sequences
        """
        # Clean the sequence (remove any non-amino acid characters)
        sequence = re.sub(r"[^A-Za-z]", "", antibody_sequence).upper()
        
        # Identify if we have heavy chain, light chain, or both
        is_vh = any(marker in sequence for marker in self.vh_markers)
        is_vl = any(marker in sequence for marker in self.vl_markers)
        
        cdrs = {}
        
        # Search for all CDRs
        for cdr_name, pattern in self.cdr_patterns.items():
            # Skip light chain CDRs if we only have heavy chain
            if cdr_name.startswith("CDRL") and not is_vl and is_vh:
                continue
                
            # Skip heavy chain CDRs if we only have light chain
            if cdr_name.startswith("CDRH") and not is_vh and is_vl:
                continue
                
            match = re.search(pattern, sequence)
            if match:
                cdr_sequence = match.group(1)
                cdrs[cdr_name] = cdr_sequence
                logger.debug(f"Found {cdr_name}: {cdr_sequence}")
            else:
                logger.debug(f"{cdr_name} not found in sequence")
        
        if not cdrs:
            # Try alternative patterns for non-standard antibodies
            cdrs = self._extract_cdrs_alternative(sequence)
        
        # Log the results
        if cdrs:
            logger.info(f"Extracted {len(cdrs)} CDRs from antibody sequence")
        else:
            logger.warning("No CDRs found in the antibody sequence")
        
        return cdrs
    
    def _extract_cdrs_alternative(self, sequence: str) -> Dict[str, str]:
        """
        Use alternative patterns for non-standard antibodies.
        
        Args:
            sequence: Cleaned antibody amino acid sequence
            
        Returns:
            Dictionary mapping CDR names to their sequences
        """
        # Alternative patterns (more relaxed, but may have false positives)
        alt_patterns = {
            # Alternative heavy chain CDRs
            "CDRH1": r"(?:C[A-Z]{2}[A-Z]{5,7}G[YF])([A-Z]{10,12})(?:W[A-Z]{3})",
            "CDRH2": r"(?:W[A-Z]{3})([A-Z]{16,19})(?:[A-Z]{4})",
            "CDRH3": r"(?:C[A-Z]{2})([A-Z]{3,25})(?:W[GY])",
            
            # Alternative light chain CDRs
            "CDRL1": r"(?:C[A-Z]{6})([A-Z]{10,17})(?:[A-Z]{4})",
            "CDRL2": r"(?:[A-Z]{4})([A-Z]{7})(?:[A-Z]{3})",
            "CDRL3": r"(?:C[A-Z]{5})([A-Z]{7,11})(?:F[GA])"
        }
        
        cdrs = {}
        for cdr_name, pattern in alt_patterns.items():
            match = re.search(pattern, sequence)
            if match:
                cdr_sequence = match.group(1)
                cdrs[cdr_name] = cdr_sequence
                logger.debug(f"Found {cdr_name} using alternative pattern: {cdr_sequence}")
        
        return cdrs
    
    def get_cdr_positions(self, antibody_sequence: str) -> Dict[str, Tuple[int, int]]:
        """
        Get the start and end positions of CDRs within the antibody sequence.
        
        Args:
            antibody_sequence: Amino acid sequence of the antibody
            
        Returns:
            Dictionary mapping CDR names to (start, end) position tuples
        """
        # Clean the sequence
        sequence = re.sub(r"[^A-Za-z]", "", antibody_sequence).upper()
        
        cdr_positions = {}
        
        # Find positions for each CDR
        for cdr_name, pattern in self.cdr_patterns.items():
            match = re.search(pattern, sequence)
            if match:
                start_pos = match.start(1)  # Start position of the CDR
                end_pos = match.end(1)      # End position of the CDR
                cdr_positions[cdr_name] = (start_pos, end_pos)
                logger.debug(f"{cdr_name} located at positions {start_pos}-{end_pos}")
        
        if not cdr_positions:
            # Try alternative patterns
            for cdr_name, pattern in self._extract_cdrs_alternative(sequence):
                match = re.search(pattern, sequence)
                if match:
                    start_pos = match.start(1)
                    end_pos = match.end(1)
                    cdr_positions[cdr_name] = (start_pos, end_pos)
        
        return cdr_positions
    
    def mask_cdr_in_sequence(self, antibody_sequence: str, cdr_name: str) -> Optional[str]:
        """
        Mask a specific CDR in the antibody sequence with <mask> tokens.
        
        Args:
            antibody_sequence: Amino acid sequence of the antibody
            cdr_name: Name of the CDR to mask (e.g., "CDRH3")
            
        Returns:
            Sequence with the specified CDR masked, or None if CDR not found
        """
        # Clean the sequence
        sequence = re.sub(r"[^A-Za-z]", "", antibody_sequence).upper()
        
        # Find the pattern for the specified CDR
        if cdr_name not in self.cdr_patterns:
            logger.error(f"Invalid CDR name: {cdr_name}")
            return None
        
        pattern = self.cdr_patterns[cdr_name]
        
        # Attempt to find and mask the CDR
        match = re.search(pattern, sequence)
        if not match:
            logger.warning(f"{cdr_name} not found in sequence")
            return None
        
        # Get the CDR sequence and its position
        cdr_sequence = match.group(1)
        start_pos = match.start(1)
        end_pos = match.end(1)
        
        # Create the masked sequence
        masked_sequence = sequence[:start_pos] + "<mask>" + sequence[end_pos:]
        
        logger.info(f"Masked {cdr_name} ({cdr_sequence}) in antibody sequence")
        return masked_sequence
    
    def replace_cdr_in_sequence(self, antibody_sequence: str, cdr_name: str, new_cdr_sequence: str) -> Optional[str]:
        """
        Replace a specific CDR in the antibody sequence with a new sequence.
        
        Args:
            antibody_sequence: Amino acid sequence of the antibody
            cdr_name: Name of the CDR to replace (e.g., "CDRH3")
            new_cdr_sequence: New amino acid sequence for the CDR
            
        Returns:
            Sequence with the specified CDR replaced, or None if CDR not found
        """
        # Clean the sequence and new CDR sequence
        sequence = re.sub(r"[^A-Za-z]", "", antibody_sequence).upper()
        new_cdr_sequence = re.sub(r"[^A-Za-z]", "", new_cdr_sequence).upper()
        
        # Find the pattern for the specified CDR
        if cdr_name not in self.cdr_patterns:
            logger.error(f"Invalid CDR name: {cdr_name}")
            return None
        
        pattern = self.cdr_patterns[cdr_name]
        
        # Attempt to find the CDR
        match = re.search(pattern, sequence)
        if not match:
            logger.warning(f"{cdr_name} not found in sequence")
            return None
        
        # Get the CDR position
        original_cdr = match.group(1)
        start_pos = match.start(1)
        end_pos = match.end(1)
        
        # Validate the length of the new CDR
        # For CDR-H3, allow more flexibility in length
        if cdr_name == "CDRH3":
            if len(new_cdr_sequence) < 3 or len(new_cdr_sequence) > 25:
                logger.warning(f"New {cdr_name} sequence length ({len(new_cdr_sequence)}) is outside typical range (3-25)")
                # Continue anyway, as CDRH3 has high variability
        else:
            # For other CDRs, prefer similar length to original
            length_diff = abs(len(new_cdr_sequence) - len(original_cdr))
            if length_diff > 3:
                logger.warning(f"New {cdr_name} sequence length ({len(new_cdr_sequence)}) differs significantly "
                              f"from original ({len(original_cdr)})")
                # Continue anyway, but log the warning
        
        # Create the modified sequence
        modified_sequence = sequence[:start_pos] + new_cdr_sequence + sequence[end_pos:]
        
        logger.info(f"Replaced {cdr_name} in antibody sequence "
                   f"(original: {original_cdr}, new: {new_cdr_sequence})")
        return modified_sequence
    
    def get_framework_regions(self, antibody_sequence: str) -> Dict[str, str]:
        """
        Get the framework regions of the antibody sequence.
        
        Args:
            antibody_sequence: Amino acid sequence of the antibody
            
        Returns:
            Dictionary mapping framework region names to their sequences
        """
        # Clean the sequence
        sequence = re.sub(r"[^A-Za-z]", "", antibody_sequence).upper()
        
        # Extract CDR positions
        cdr_positions = self.get_cdr_positions(sequence)
        
        if not cdr_positions:
            logger.warning("No CDRs found, cannot identify framework regions")
            return {}
        
        # Sort CDRs by position
        sorted_cdrs = sorted(cdr_positions.items(), key=lambda x: x[1][0])
        
        # Determine if we have heavy chain, light chain, or both
        has_vh = any(cdr.startswith("CDRH") for cdr in cdr_positions.keys())
        has_vl = any(cdr.startswith("CDRL") for cdr in cdr_positions.keys())
        
        framework_regions = {}
        
        # Extract framework regions
        if has_vh:
            vh_cdrs = [(name, pos) for name, pos in sorted_cdrs if name.startswith("CDRH")]
            vh_cdrs.sort(key=lambda x: x[1][0])  # Sort by start position
            
            if vh_cdrs:
                # FWR1-H: Start to CDRH1
                if "CDRH1" in dict(vh_cdrs):
                    h1_start = dict(vh_cdrs)["CDRH1"][0]
                    framework_regions["FWR1-H"] = sequence[:h1_start]
                
                # Extract other framework regions
                for i in range(len(vh_cdrs) - 1):
                    current_cdr_name, (_, current_end) = vh_cdrs[i]
                    next_cdr_name, (next_start, _) = vh_cdrs[i + 1]
                    fw_name = f"FWR{i+2}-H"
                    framework_regions[fw_name] = sequence[current_end:next_start]
                
                # FWR4-H: After last CDR-H
                _, (_, last_end) = vh_cdrs[-1]
                # Find the end of the VH region (usually around WGQG sequence)
                vh_end_match = re.search(r"W[GY][QA]G[TG][TV][VL][TS][VS]S", sequence[last_end:])
                if vh_end_match:
                    vh_end = last_end + vh_end_match.end()
                    framework_regions["FWR4-H"] = sequence[last_end:vh_end]
                else:
                    # If we can't find the end marker, take a reasonable portion
                    framework_regions["FWR4-H"] = sequence[last_end:last_end+20]
        
        if has_vl:
            vl_cdrs = [(name, pos) for name, pos in sorted_cdrs if name.startswith("CDRL")]
            vl_cdrs.sort(key=lambda x: x[1][0])  # Sort by start position
            
            if vl_cdrs:
                # FWR1-L: Start to CDRL1
                if "CDRL1" in dict(vl_cdrs):
                    l1_start = dict(vl_cdrs)["CDRL1"][0]
                    framework_regions["FWR1-L"] = sequence[:l1_start]
                
                # Extract other framework regions
                for i in range(len(vl_cdrs) - 1):
                    current_cdr_name, (_, current_end) = vl_cdrs[i]
                    next_cdr_name, (next_start, _) = vl_cdrs[i + 1]
                    fw_name = f"FWR{i+2}-L"
                    framework_regions[fw_name] = sequence[current_end:next_start]
                
                # FWR4-L: After last CDR-L
                _, (_, last_end) = vl_cdrs[-1]
                # Find the end of the VL region
                vl_end_match = re.search(r"F[GA][A-Z]G[A-Z]{2}[GA][A-Z]{2}", sequence[last_end:])
                if vl_end_match:
                    vl_end = last_end + vl_end_match.end()
                    framework_regions["FWR4-L"] = sequence[last_end:vl_end]
                else:
                    # If we can't find the end marker, take a reasonable portion
                    framework_regions["FWR4-L"] = sequence[last_end:last_end+15]
        
        if framework_regions:
            logger.info(f"Extracted {len(framework_regions)} framework regions from antibody sequence")
        else:
            logger.warning("Failed to extract framework regions")
        
        return framework_regions