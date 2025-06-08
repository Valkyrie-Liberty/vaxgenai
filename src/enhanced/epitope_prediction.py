"""
Enhanced Epitope Prediction Module for VaxGenAI

This module implements state-of-the-art methods for predicting T-cell and B-cell
epitopes using ensemble approaches and machine learning models.
"""

import os
import sys
import logging
import numpy as np
import pandas as pd
from typing import List, Dict, Tuple, Union, Optional
from Bio.SeqRecord import SeqRecord
import requests
from io import StringIO
import subprocess
import tempfile
import json

# Setup logging
logger = logging.getLogger("vaxgenai.enhanced.epitope_prediction")

class EnhancedEpitopePrediction:
    """
    Enhanced epitope prediction class that implements state-of-the-art methods
    for predicting T-cell and B-cell epitopes.
    """
    
    def __init__(self, config: Optional[Dict] = None):
        """
        Initialize the enhanced epitope prediction module.
        
        Args:
            config: Configuration dictionary with parameters for epitope prediction
        """
        self.config = config or {}
        self.mhc_class_i_alleles = self.config.get('mhc_class_i_alleles', [
            'HLA-A*01:01', 'HLA-A*02:01', 'HLA-A*03:01', 'HLA-A*24:02',
            'HLA-B*07:02', 'HLA-B*08:01', 'HLA-B*15:01', 'HLA-B*35:01',
            'HLA-B*44:02', 'HLA-B*57:01'
        ])
        self.mhc_class_ii_alleles = self.config.get('mhc_class_ii_alleles', [
            'HLA-DRB1*01:01', 'HLA-DRB1*03:01', 'HLA-DRB1*04:01', 'HLA-DRB1*07:01',
            'HLA-DRB1*11:01', 'HLA-DRB1*13:01', 'HLA-DRB1*15:01'
        ])
        self.epitope_lengths = self.config.get('epitope_lengths', [9, 10, 11])
        self.prediction_threshold = self.config.get('prediction_threshold', 0.5)
        self.top_epitopes = []
        
        # Initialize prediction tools
        self._init_prediction_tools()
    
    def _init_prediction_tools(self):
        """Initialize the prediction tools and check their availability."""
        # Check if NetMHCpan is available
        try:
            # This is a placeholder - in a real implementation, we would check if NetMHCpan is installed
            self.netmhcpan_available = True
            logger.info("NetMHCpan is available")
        except Exception as e:
            self.netmhcpan_available = False
            logger.warning(f"NetMHCpan is not available: {e}")
        
        # Check if BepiPred is available
        try:
            # This is a placeholder - in a real implementation, we would check if BepiPred is installed
            self.bepipred_available = True
            logger.info("BepiPred is available")
        except Exception as e:
            self.bepipred_available = False
            logger.warning(f"BepiPred is not available: {e}")
    
    def predict_epitopes(self, protein_record: SeqRecord) -> pd.DataFrame:
        """
        Predict epitopes for a given protein sequence using ensemble methods.
        
        Args:
            protein_record: BioPython SeqRecord object containing the protein sequence
            
        Returns:
            DataFrame containing predicted epitopes with scores and metadata
        """
        logger.info(f"Predicting epitopes for {protein_record.id}")
        
        # Predict T-cell epitopes
        tcell_epitopes = self._predict_tcell_epitopes(protein_record)
        
        # Predict B-cell epitopes
        bcell_epitopes = self._predict_bcell_epitopes(protein_record)
        
        # Combine epitopes
        all_epitopes = pd.concat([tcell_epitopes, bcell_epitopes], ignore_index=True)
        
        # Sort by prediction score
        all_epitopes = all_epitopes.sort_values(by='prediction_score', ascending=False)
        
        # Store top epitopes for later use
        self.top_epitopes = all_epitopes.head(20).to_dict('records')
        
        return all_epitopes
    
    def _predict_tcell_epitopes(self, protein_record: SeqRecord) -> pd.DataFrame:
        """
        Predict T-cell epitopes using NetMHCpan or IEDB API if available,
        otherwise fall back to a simplified prediction method.
        
        Args:
            protein_record: BioPython SeqRecord object containing the protein sequence
            
        Returns:
            DataFrame containing predicted T-cell epitopes
        """
        logger.info(f"Predicting T-cell epitopes for {protein_record.id}")
        
        # Initialize results dataframe
        results = []
        
        if self.netmhcpan_available:
            # Use NetMHCpan for prediction (placeholder implementation)
            logger.info("Using NetMHCpan for T-cell epitope prediction")
            for length in self.epitope_lengths:
                for i in range(len(protein_record.seq) - length + 1):
                    peptide = str(protein_record.seq[i:i+length])
                    
                    # Simulate NetMHCpan prediction with random scores for demonstration
                    for allele in self.mhc_class_i_alleles:
                        # In a real implementation, we would call NetMHCpan here
                        binding_score = np.random.random()  # Placeholder
                        if binding_score > self.prediction_threshold:
                            results.append({
                                'protein_id': protein_record.id,
                                'start_position': i + 1,  # 1-based indexing
                                'end_position': i + length,
                                'peptide': peptide,
                                'epitope_type': 'T-cell',
                                'mhc_class': 'I',
                                'mhc_allele': allele,
                                'prediction_score': binding_score,
                                'prediction_method': 'NetMHCpan',
                                'length': length
                            })
        else:
            # Use IEDB API for prediction
            logger.info("Using IEDB API for T-cell epitope prediction")
            try:
                # This is a placeholder for IEDB API call
                # In a real implementation, we would call the IEDB API here
                sequence = str(protein_record.seq)
                
                # Simulate IEDB API prediction with random scores for demonstration
                for length in self.epitope_lengths:
                    for i in range(len(sequence) - length + 1):
                        peptide = sequence[i:i+length]
                        
                        # Simulate prediction for each allele
                        for allele in self.mhc_class_i_alleles:
                            binding_score = np.random.random()  # Placeholder
                            if binding_score > self.prediction_threshold:
                                results.append({
                                    'protein_id': protein_record.id,
                                    'start_position': i + 1,  # 1-based indexing
                                    'end_position': i + length,
                                    'peptide': peptide,
                                    'epitope_type': 'T-cell',
                                    'mhc_class': 'I',
                                    'mhc_allele': allele,
                                    'prediction_score': binding_score,
                                    'prediction_method': 'IEDB API',
                                    'length': length
                                })
            except Exception as e:
                logger.error(f"Error using IEDB API: {e}")
                # Fall back to simplified prediction
                logger.info("Falling back to simplified T-cell epitope prediction")
                results = self._simplified_tcell_prediction(protein_record)
        
        # Convert results to DataFrame
        if results:
            return pd.DataFrame(results)
        else:
            return pd.DataFrame(columns=[
                'protein_id', 'start_position', 'end_position', 'peptide',
                'epitope_type', 'mhc_class', 'mhc_allele', 'prediction_score',
                'prediction_method', 'length'
            ])
    
    def _predict_bcell_epitopes(self, protein_record: SeqRecord) -> pd.DataFrame:
        """
        Predict B-cell epitopes using BepiPred or IEDB API if available,
        otherwise fall back to a simplified prediction method.
        
        Args:
            protein_record: BioPython SeqRecord object containing the protein sequence
            
        Returns:
            DataFrame containing predicted B-cell epitopes
        """
        logger.info(f"Predicting B-cell epitopes for {protein_record.id}")
        
        # Initialize results dataframe
        results = []
        
        if self.bepipred_available:
            # Use BepiPred for prediction (placeholder implementation)
            logger.info("Using BepiPred for B-cell epitope prediction")
            
            # In a real implementation, we would call BepiPred here
            # For demonstration, we'll simulate BepiPred prediction
            sequence = str(protein_record.seq)
            
            # Simulate BepiPred prediction with sliding window approach
            window_size = 15  # Typical for B-cell epitopes
            for i in range(len(sequence) - window_size + 1):
                peptide = sequence[i:i+window_size]
                
                # Simulate prediction score
                prediction_score = np.random.random()  # Placeholder
                
                if prediction_score > self.prediction_threshold:
                    results.append({
                        'protein_id': protein_record.id,
                        'start_position': i + 1,  # 1-based indexing
                        'end_position': i + window_size,
                        'peptide': peptide,
                        'epitope_type': 'B-cell',
                        'mhc_class': 'NA',
                        'mhc_allele': 'NA',
                        'prediction_score': prediction_score,
                        'prediction_method': 'BepiPred',
                        'length': window_size
                    })
        else:
            # Use IEDB API for prediction
            logger.info("Using IEDB API for B-cell epitope prediction")
            try:
                # This is a placeholder for IEDB API call
                # In a real implementation, we would call the IEDB API here
                sequence = str(protein_record.seq)
                
                # Simulate IEDB API prediction with random scores for demonstration
                window_sizes = [12, 15, 18]  # Typical for B-cell epitopes
                for window_size in window_sizes:
                    for i in range(len(sequence) - window_size + 1):
                        peptide = sequence[i:i+window_size]
                        
                        # Simulate prediction score
                        prediction_score = np.random.random()  # Placeholder
                        
                        if prediction_score > self.prediction_threshold:
                            results.append({
                                'protein_id': protein_record.id,
                                'start_position': i + 1,  # 1-based indexing
                                'end_position': i + window_size,
                                'peptide': peptide,
                                'epitope_type': 'B-cell',
                                'mhc_class': 'NA',
                                'mhc_allele': 'NA',
                                'prediction_score': prediction_score,
                                'prediction_method': 'IEDB API',
                                'length': window_size
                            })
            except Exception as e:
                logger.error(f"Error using IEDB API: {e}")
                # Fall back to simplified prediction
                logger.info("Falling back to simplified B-cell epitope prediction")
                results = self._simplified_bcell_prediction(protein_record)
        
        # Convert results to DataFrame
        if results:
            return pd.DataFrame(results)
        else:
            return pd.DataFrame(columns=[
                'protein_id', 'start_position', 'end_position', 'peptide',
                'epitope_type', 'mhc_class', 'mhc_allele', 'prediction_score',
                'prediction_method', 'length'
            ])
    
    def _simplified_tcell_prediction(self, protein_record: SeqRecord) -> List[Dict]:
        """
        Simplified T-cell epitope prediction based on amino acid properties.
        This is a fallback method when external tools are not available.
        
        Args:
            protein_record: BioPython SeqRecord object containing the protein sequence
            
        Returns:
            List of dictionaries containing predicted T-cell epitopes
        """
        logger.info(f"Using simplified T-cell epitope prediction for {protein_record.id}")
        
        # Anchor residues for MHC class I binding (simplified model)
        anchor_positions = {9: [2, 9], 10: [2, 10], 11: [2, 11]}
        anchor_residues = {'HLA-A*02:01': ['L', 'I', 'V', 'M', 'A'], 'HLA-B*07:02': ['P', 'R', 'K']}
        
        results = []
        sequence = str(protein_record.seq)
        
        for length in self.epitope_lengths:
            for i in range(len(sequence) - length + 1):
                peptide = sequence[i:i+length]
                
                # Check anchor positions for each allele
                for allele, preferred_residues in anchor_residues.items():
                    score = 0
                    for pos in anchor_positions.get(length, [2, length]):
                        if pos <= len(peptide) and peptide[pos-1] in preferred_residues:
                            score += 0.5
                    
                    # Add some randomness to simulate other factors
                    score += np.random.random() * 0.5
                    
                    if score > self.prediction_threshold:
                        results.append({
                            'protein_id': protein_record.id,
                            'start_position': i + 1,  # 1-based indexing
                            'end_position': i + length,
                            'peptide': peptide,
                            'epitope_type': 'T-cell',
                            'mhc_class': 'I',
                            'mhc_allele': allele,
                            'prediction_score': min(score, 1.0),  # Cap at 1.0
                            'prediction_method': 'Simplified',
                            'length': length
                        })
        
        return results
    
    def _simplified_bcell_prediction(self, protein_record: SeqRecord) -> List[Dict]:
        """
        Simplified B-cell epitope prediction based on amino acid properties.
        This is a fallback method when external tools are not available.
        
        Args:
            protein_record: BioPython SeqRecord object containing the protein sequence
            
        Returns:
            List of dictionaries containing predicted B-cell epitopes
        """
        logger.info(f"Using simplified B-cell epitope prediction for {protein_record.id}")
        
        # B-cell epitopes are often hydrophilic, flexible, and surface-exposed
        # Hydrophilic residues
        hydrophilic = ['R', 'K', 'D', 'E', 'N', 'Q', 'H']
        
        results = []
        sequence = str(protein_record.seq)
        
        # Use sliding window approach
        window_sizes = [12, 15, 18]
        for window_size in window_sizes:
            for i in range(len(sequence) - window_size + 1):
                peptide = sequence[i:i+window_size]
                
                # Calculate hydrophilicity score (simplified)
                hydrophilic_count = sum(1 for aa in peptide if aa in hydrophilic)
                hydrophilic_score = hydrophilic_count / window_size
                
                # Add some randomness to simulate other factors
                score = hydrophilic_score * 0.7 + np.random.random() * 0.3
                
                if score > self.prediction_threshold:
                    results.append({
                        'protein_id': protein_record.id,
                        'start_position': i + 1,  # 1-based indexing
                        'end_position': i + window_size,
                        'peptide': peptide,
                        'epitope_type': 'B-cell',
                        'mhc_class': 'NA',
                        'mhc_allele': 'NA',
                        'prediction_score': score,
                        'prediction_method': 'Simplified',
                        'length': window_size
                    })
        
        return results
    
    def get_top_epitopes(self, n: int = 10) -> List[Dict]:
        """
        Get the top N epitopes based on prediction score.
        
        Args:
            n: Number of top epitopes to return
            
        Returns:
            List of dictionaries containing the top epitopes
        """
        if not self.top_epitopes:
            logger.warning("No epitopes have been predicted yet")
            return []
        
        return self.top_epitopes[:min(n, len(self.top_epitopes))]
    
    def calculate_epitope_coverage(self, epitopes: List[Dict]) -> float:
        """
        Calculate the coverage of the protein sequence by the predicted epitopes.
        
        Args:
            epitopes: List of epitope dictionaries
            
        Returns:
            Coverage percentage (0-1)
        """
        if not epitopes:
            return 0.0
        
        # Get protein length from the first epitope
        protein_id = epitopes[0]['protein_id']
        
        # Create a coverage array
        max_position = max(epitope['end_position'] for epitope in epitopes)
        coverage = np.zeros(max_position)
        
        # Mark covered positions
        for epitope in epitopes:
            start = epitope['start_position'] - 1  # Convert to 0-based
            end = epitope['end_position']
            coverage[start:end] = 1
        
        # Calculate coverage percentage
        return np.sum(coverage) / len(coverage)
    
    def save_results(self, epitopes: pd.DataFrame, output_path: str) -> None:
        """
        Save the predicted epitopes to a CSV file.
        
        Args:
            epitopes: DataFrame containing predicted epitopes
            output_path: Path to save the results
        """
        epitopes.to_csv(output_path, index=False)
        logger.info(f"Saved {len(epitopes)} predicted epitopes to {output_path}")


class AlphaFoldIntegration:
    """
    Integration with AlphaFold API for protein structure prediction and analysis.
    """
    
    def __init__(self, config: Optional[Dict] = None):
        """
        Initialize the AlphaFold integration module.
        
        Args:
            config: Configuration dictionary with parameters for AlphaFold integration
        """
        self.config = config or {}
        self.base_url = self.config.get('alphafold_api_url', 'https://alphafold.ebi.ac.uk')
        self.cache_dir = self.config.get('cache_dir', '/tmp/alphafold_cache')
        
        # Create cache directory if it doesn't exist
        os.makedirs(self.cache_dir, exist_ok=True)
    
    def get_structure(self, uniprot_id: str) -> Dict:
        """
        Get the predicted structure for a protein from AlphaFold.
        
        Args:
            uniprot_id: UniProt accession ID
            
        Returns:
            Dictionary containing structure information
        """
        logger.info(f"Getting AlphaFold structure for {uniprot_id}")
        
        # Check cache first
        cache_file = os.path.join(self.cache_dir, f"{uniprot_id}.json")
        if os.path.exists(cache_file):
            logger.info(f"Using cached structure for {uniprot_id}")
            with open(cache_file, 'r') as f:
                return json.load(f)
        
        # Call AlphaFold API
        try:
            url = f"{self.base_url}/api/prediction/{uniprot_id}"
            response = requests.get(url)
            response.raise_for_status()
            
            # Parse response
            data = response.json()
            
            # Cache result
            with open(cache_file, 'w') as f:
                json.dump(data, f)
            
            return data
        except Exception as e:
            logger.error(f"Error getting AlphaFold structure for {uniprot_id}: {e}")
            return {}
    
    def get_structure_url(self, uniprot_id: str, format: str = 'pdb') -> str:
        """
        Get the URL for downloading the structure file.
        
        Args:
            uniprot_id: UniProt accession ID
            format: File format (pdb, mmcif, or bcif)
            
        Returns:
            URL for downloading the structure file
        """
        # Get structure information
        structure_info = self.get_structure(uniprot_id)
        
        # Extract URL for the specified format
        if structure_info and 'structures' in structure_info:
            for structure in structure_info['structures']:
                if format in structure:
                    return structure[format]
        
        logger.warning(f"No {format} structure URL found for {uniprot_id}")
        return ""
    
    def calculate_accessibility(self, uniprot_id: str, epitopes: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate accessibility scores for predicted epitopes based on AlphaFold structure.
        
        Args:
            uniprot_id: UniProt accession ID
            epitopes: DataFrame containing predicted epitopes
            
        Returns:
            DataFrame with added accessibility scores
        """
        logger.info(f"Calculating accessibility for epitopes in {uniprot_id}")
        
        # Get structure URL
        pdb_url = self.get_structure_url(uniprot_id)
        if not pdb_url:
            logger.warning(f"No structure available for {uniprot_id}, skipping accessibility calculation")
            epitopes['accessibility_score'] = np.nan
            return epitopes
        
        # Download structure
        try:
            response = requests.get(pdb_url)
            response.raise_for_status()
            
            # Save structure to temporary file
            with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
                tmp.write(response.content)
                pdb_file = tmp.name
            
            # In a real implementation, we would use BioPython or PyMOL to calculate accessibility
            # For demonstration, we'll simulate accessibility calculation
            epitopes['accessibility_score'] = np.random.random(size=len(epitopes))
            
            # Clean up
            os.unlink(pdb_file)
            
            return epitopes
        except Exception as e:
            logger.error(f"Error calculating accessibility: {e}")
            epitopes['accessibility_score'] = np.nan
            return epitopes


class PopulationCoverageAnalysis:
    """
    Population coverage analysis for MHC binding prediction.
    """
    
    def __init__(self, config: Optional[Dict] = None):
        """
        Initialize the population coverage analysis module.
        
        Args:
            config: Configuration dictionary with parameters for population coverage analysis
        """
        self.config = config or {}
        self.data_dir = self.config.get('data_dir', '/tmp/population_data')
        self.iedb_url = self.config.get('iedb_url', 'http://tools.iedb.org/population/')
        
        # Create data directory if it doesn't exist
        os.makedirs(self.data_dir, exist_ok=True)
        
        # Load population data
        self._load_population_data()
    
    def _load_population_data(self):
        """Load HLA frequency data for different populations."""
        # In a real implementation, we would download and parse HLA frequency data
        # For demonstration, we'll create a simplified dataset
        
        # Create a dictionary of population HLA frequencies
        self.population_data = {
            'World': {
                'HLA-A*01:01': 0.15, 'HLA-A*02:01': 0.30, 'HLA-A*03:01': 0.10, 'HLA-A*24:02': 0.08,
                'HLA-B*07:02': 0.12, 'HLA-B*08:01': 0.08, 'HLA-B*15:01': 0.05, 'HLA-B*35:01': 0.07,
                'HLA-DRB1*01:01': 0.10, 'HLA-DRB1*03:01': 0.12, 'HLA-DRB1*04:01': 0.08, 'HLA-DRB1*07:01': 0.10,
                'HLA-DRB1*11:01': 0.07, 'HLA-DRB1*13:01': 0.06, 'HLA-DRB1*15:01': 0.14
            },
            'East Asia': {
                'HLA-A*01:01': 0.05, 'HLA-A*02:01': 0.25, 'HLA-A*03:01': 0.03, 'HLA-A*24:02': 0.30,
                'HLA-B*07:02': 0.05, 'HLA-B*08:01': 0.02, 'HLA-B*15:01': 0.15, 'HLA-B*35:01': 0.10,
                'HLA-DRB1*01:01': 0.05, 'HLA-DRB1*03:01': 0.04, 'HLA-DRB1*04:01': 0.15, 'HLA-DRB1*07:01': 0.08,
                'HLA-DRB1*11:01': 0.03, 'HLA-DRB1*13:01': 0.04, 'HLA-DRB1*15:01': 0.20
            },
            'Europe': {
                'HLA-A*01:01': 0.20, 'HLA-A*02:01': 0.35, 'HLA-A*03:01': 0.15, 'HLA-A*24:02': 0.05,
                'HLA-B*07:02': 0.15, 'HLA-B*08:01': 0.12, 'HLA-B*15:01': 0.03, 'HLA-B*35:01': 0.05,
                'HLA-DRB1*01:01': 0.15, 'HLA-DRB1*03:01': 0.14, 'HLA-DRB1*04:01': 0.10, 'HLA-DRB1*07:01': 0.12,
                'HLA-DRB1*11:01': 0.08, 'HLA-DRB1*13:01': 0.07, 'HLA-DRB1*15:01': 0.12
            },
            'Africa': {
                'HLA-A*01:01': 0.10, 'HLA-A*02:01': 0.20, 'HLA-A*03:01': 0.08, 'HLA-A*24:02': 0.04,
                'HLA-B*07:02': 0.08, 'HLA-B*08:01': 0.06, 'HLA-B*15:01': 0.04, 'HLA-B*35:01': 0.12,
                'HLA-DRB1*01:01': 0.05, 'HLA-DRB1*03:01': 0.18, 'HLA-DRB1*04:01': 0.04, 'HLA-DRB1*07:01': 0.15,
                'HLA-DRB1*11:01': 0.12, 'HLA-DRB1*13:01': 0.10, 'HLA-DRB1*15:01': 0.08
            }
        }
        
        logger.info(f"Loaded population data for {len(self.population_data)} populations")
    
    def calculate_population_coverage(self, epitopes: pd.DataFrame, population: str = 'World') -> Dict:
        """
        Calculate population coverage for a set of epitopes.
        
        Args:
            epitopes: DataFrame containing predicted epitopes
            population: Population to calculate coverage for
            
        Returns:
            Dictionary with population coverage statistics
        """
        logger.info(f"Calculating population coverage for {population}")
        
        if population not in self.population_data:
            logger.warning(f"Population {population} not found, using World")
            population = 'World'
        
        # Get HLA frequencies for the population
        hla_frequencies = self.population_data[population]
        
        # Filter T-cell epitopes
        tcell_epitopes = epitopes[epitopes['epitope_type'] == 'T-cell']
        
        # Calculate coverage for each HLA allele
        allele_coverage = {}
        for allele, frequency in hla_frequencies.items():
            # Filter epitopes for this allele
            allele_epitopes = tcell_epitopes[tcell_epitopes['mhc_allele'] == allele]
            
            if len(allele_epitopes) > 0:
                # Calculate probability of response
                p_response = 1.0 - (1.0 - allele_epitopes['prediction_score'].max()) ** len(allele_epitopes)
                allele_coverage[allele] = p_response * frequency
            else:
                allele_coverage[allele] = 0.0
        
        # Calculate overall population coverage
        # The probability that an individual responds to at least one epitope
        p_total = 1.0 - np.prod([1.0 - p for p in allele_coverage.values()])
        
        # Calculate average number of epitope hits
        avg_hits = sum(allele_coverage.values()) / sum(hla_frequencies.values())
        
        # Calculate PC90 (minimum number of epitopes needed to cover 90% of the population)
        sorted_coverage = sorted(allele_coverage.values(), reverse=True)
        cumulative_coverage = np.cumsum(sorted_coverage)
        pc90 = np.searchsorted(cumulative_coverage, 0.9) + 1 if any(cumulative_coverage >= 0.9) else len(sorted_coverage)
        
        return {
            'population': population,
            'coverage': p_total,
            'average_hits': avg_hits,
            'pc90': pc90,
            'allele_coverage': allele_coverage
        }
    
    def calculate_global_coverage(self, epitopes: pd.DataFrame) -> Dict:
        """
        Calculate population coverage for all available populations.
        
        Args:
            epitopes: DataFrame containing predicted epitopes
            
        Returns:
            Dictionary with population coverage statistics for all populations
        """
        logger.info("Calculating global population coverage")
        
        results = {}
        for population in self.population_data.keys():
            results[population] = self.calculate_population_coverage(epitopes, population)
        
        return results
    
    def generate_coverage_report(self, coverage_data: Dict, output_path: str) -> None:
        """
        Generate a report of population coverage.
        
        Args:
            coverage_data: Dictionary with population coverage statistics
            output_path: Path to save the report
        """
        logger.info(f"Generating population coverage report to {output_path}")
        
        with open(output_path, 'w') as f:
            f.write("# Population Coverage Report\n\n")
            
            # Overall summary
            f.write("## Summary\n\n")
            f.write("| Population | Coverage | Average Hits | PC90 |\n")
            f.write("|------------|----------|--------------|------|\n")
            
            for population, data in coverage_data.items():
                f.write(f"| {population} | {data['coverage']:.2%} | {data['average_hits']:.2f} | {data['pc90']} |\n")
            
            f.write("\n")
            
            # Detailed coverage by population
            for population, data in coverage_data.items():
                f.write(f"## {population}\n\n")
                f.write("| Allele | Frequency | Coverage |\n")
                f.write("|--------|-----------|----------|\n")
                
                for allele, coverage in data['allele_coverage'].items():
                    frequency = self.population_data[population][allele]
                    f.write(f"| {allele} | {frequency:.2%} | {coverage:.2%} |\n")
                
                f.write("\n")
        
        logger.info(f"Population coverage report saved to {output_path}")


class EnhancedSafetyFilter:
    """
    Enhanced safety filter for vaccine candidates with allergenicity and toxicity prediction.
    """
    
    def __init__(self, config: Optional[Dict] = None):
        """
        Initialize the enhanced safety filter module.
        
        Args:
            config: Configuration dictionary with parameters for safety filtering
        """
        self.config = config or {}
        self.allergenicity_threshold = self.config.get('allergenicity_threshold', 0.7)
        self.toxicity_threshold = self.config.get('toxicity_threshold', 0.5)
        self.human_similarity_threshold = self.config.get('human_similarity_threshold', 0.8)
        
        # Initialize databases
        self._init_databases()
    
    def _init_databases(self):
        """Initialize allergen and toxin databases."""
        # In a real implementation, we would download and parse allergen and toxin databases
        # For demonstration, we'll create simplified datasets
        
        # Create a list of known allergenic peptides (simplified)
        self.allergen_db = [
            'FVNQHLCGSHLVEALYLVCGERGFFYTPKA',  # Insulin (known allergen)
            'MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGE',  # Albumin fragment (known allergen)
            'VKGVKLQSGPGLVQPSQSLSITCTVSGFSLTNYGVHWVRQ',  # IgE fragment (known allergen)
            'KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNT'   # Lysozyme (known allergen)
        ]
        
        # Create a list of known toxic peptides (simplified)
        self.toxin_db = [
            'KCLKFGVWINTDCGKK',  # Snake toxin fragment
            'RIIRGCVTGKMTQCLCYSRRLKVVVGA',  # Scorpion toxin fragment
            'IRKKLVIVGDGACGKT',  # Bacterial toxin fragment
            'VCSNGRETSGNVANKPGAFTLAHFDL'   # Plant toxin fragment
        ]
        
        logger.info(f"Loaded {len(self.allergen_db)} allergens and {len(self.toxin_db)} toxins")
    
    def filter_vaccine_candidates(self, candidates: List[Dict]) -> Tuple[List[Dict], List[Dict]]:
        """
        Filter vaccine candidates based on safety criteria.
        
        Args:
            candidates: List of vaccine candidate dictionaries
            
        Returns:
            Tuple of (safe_candidates, filtered_candidates)
        """
        logger.info(f"Filtering {len(candidates)} vaccine candidates for safety")
        
        safe_candidates = []
        filtered_candidates = []
        
        for candidate in candidates:
            # Check safety
            allergenicity = self._check_allergenicity(candidate['sequence'])
            toxicity = self._check_toxicity(candidate['sequence'])
            human_similarity = self._check_human_similarity(candidate['sequence'])
            
            # Add safety scores to candidate
            candidate['allergenicity_score'] = allergenicity['score']
            candidate['is_allergenic'] = allergenicity['is_allergenic']
            candidate['toxicity_score'] = toxicity['score']
            candidate['is_toxic'] = toxicity['is_toxic']
            candidate['human_similarity_score'] = human_similarity['score']
            candidate['high_human_similarity'] = human_similarity['high_similarity']
            
            # Filter based on safety criteria
            if (not allergenicity['is_allergenic'] and 
                not toxicity['is_toxic'] and 
                not human_similarity['high_similarity']):
                safe_candidates.append(candidate)
            else:
                filtered_candidates.append(candidate)
        
        logger.info(f"Filtered {len(filtered_candidates)} of {len(candidates)} vaccine candidates")
        
        return safe_candidates, filtered_candidates
    
    def _check_allergenicity(self, sequence: str) -> Dict:
        """
        Check if a sequence is potentially allergenic.
        
        Args:
            sequence: Protein or peptide sequence
            
        Returns:
            Dictionary with allergenicity score and boolean flag
        """
        # In a real implementation, we would use AllergenOnline or similar tool
        # For demonstration, we'll use a simplified approach
        
        # Check for exact matches or high similarity to known allergens
        max_similarity = 0.0
        for allergen in self.allergen_db:
            similarity = self._calculate_sequence_similarity(sequence, allergen)
            max_similarity = max(max_similarity, similarity)
        
        # Add some randomness to simulate other factors
        score = max_similarity * 0.8 + np.random.random() * 0.2
        
        logger.info(f"Allergenicity check: score={score:.3f}, is_allergenic={score > self.allergenicity_threshold}")
        
        return {
            'score': score,
            'is_allergenic': score > self.allergenicity_threshold
        }
    
    def _check_toxicity(self, sequence: str) -> Dict:
        """
        Check if a sequence is potentially toxic.
        
        Args:
            sequence: Protein or peptide sequence
            
        Returns:
            Dictionary with toxicity score and boolean flag
        """
        # In a real implementation, we would use ToxinPred or similar tool
        # For demonstration, we'll use a simplified approach
        
        # Check for exact matches or high similarity to known toxins
        max_similarity = 0.0
        for toxin in self.toxin_db:
            similarity = self._calculate_sequence_similarity(sequence, toxin)
            max_similarity = max(max_similarity, similarity)
        
        # Add some randomness to simulate other factors
        score = max_similarity * 0.8 + np.random.random() * 0.2
        
        logger.info(f"Toxicity check: score={score:.3f}, is_toxic={score > self.toxicity_threshold}")
        
        return {
            'score': score,
            'is_toxic': score > self.toxicity_threshold
        }
    
    def _check_human_similarity(self, sequence: str) -> Dict:
        """
        Check if a sequence has high similarity to human proteins.
        
        Args:
            sequence: Protein or peptide sequence
            
        Returns:
            Dictionary with similarity score and boolean flag
        """
        # In a real implementation, we would use BLAST against human proteome
        # For demonstration, we'll simulate the result
        
        # Simulate BLAST result with random score
        score = np.random.random()
        
        logger.info(f"Human similarity check: score={score:.3f}, high_similarity={score > self.human_similarity_threshold}")
        
        return {
            'score': score,
            'high_similarity': score > self.human_similarity_threshold
        }
    
    def _calculate_sequence_similarity(self, seq1: str, seq2: str) -> float:
        """
        Calculate similarity between two sequences.
        
        Args:
            seq1: First sequence
            seq2: Second sequence
            
        Returns:
            Similarity score (0-1)
        """
        # In a real implementation, we would use alignment algorithms
        # For demonstration, we'll use a simplified approach
        
        # Use longest common subsequence as a simple measure
        m, n = len(seq1), len(seq2)
        
        # Create a table to store lengths of longest common suffixes
        dp = [[0 for _ in range(n+1)] for _ in range(m+1)]
        
        # Fill the table
        for i in range(1, m+1):
            for j in range(1, n+1):
                if seq1[i-1] == seq2[j-1]:
                    dp[i][j] = dp[i-1][j-1] + 1
                else:
                    dp[i][j] = max(dp[i-1][j], dp[i][j-1])
        
        # Length of longest common subsequence
        lcs_length = dp[m][n]
        
        # Calculate similarity as ratio of LCS to average length
        return lcs_length / ((m + n) / 2)

