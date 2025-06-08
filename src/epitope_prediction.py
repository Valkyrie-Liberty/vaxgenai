"""
Epitope Prediction Module for VaxGenAI

This module predicts potential epitopes from protein sequences using various
prediction algorithms for B-cell and T-cell epitopes.
"""

import os
import logging
import numpy as np
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .utils import get_data_path, get_results_path

logger = logging.getLogger("vaxgenai.epitope_prediction")

class EpitopePrediction:
    """
    Predicts B-cell and T-cell epitopes from protein sequences
    """
    
    def __init__(self):
        """Initialize the epitope prediction module"""
        self.data_path = get_data_path()
        self.results_path = get_results_path()
        self.epitope_results = None
        logger.info("Epitope prediction module initialized")
    
    def predict_bcell_epitopes(self, protein_record, window_size=20, threshold=0.5):
        """
        Predict B-cell epitopes using a simple propensity scale method
        
        Args:
            protein_record: SeqRecord object with protein sequence
            window_size: Size of the sliding window
            threshold: Threshold for epitope prediction
            
        Returns:
            list: List of predicted epitopes with scores
        """
        # This is a simplified implementation using hydrophilicity as a proxy
        # In a real implementation, we would use more sophisticated methods
        
        # Kyte & Doolittle hydrophobicity scale (negative values are hydrophilic)
        hydrophobicity = {
            'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8,
            'G': -0.4, 'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8,
            'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5,
            'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
        }
        
        sequence = str(protein_record.seq).upper()
        scores = []
        
        # Calculate hydrophilicity (negative of hydrophobicity)
        for aa in sequence:
            if aa in hydrophobicity:
                scores.append(-hydrophobicity[aa])  # Convert to hydrophilicity
            else:
                scores.append(0)  # Default for unknown amino acids
        
        # Apply sliding window
        window_scores = []
        epitopes = []
        
        for i in range(len(sequence) - window_size + 1):
            window = sequence[i:i+window_size]
            avg_score = sum(scores[i:i+window_size]) / window_size
            window_scores.append(avg_score)
            
            # If score exceeds threshold, consider it a potential epitope
            if avg_score > threshold:
                epitopes.append({
                    'start': i + 1,  # 1-based indexing
                    'end': i + window_size,
                    'sequence': window,
                    'score': avg_score,
                    'type': 'B-cell'
                })
        
        logger.info(f"Predicted {len(epitopes)} B-cell epitopes")
        return epitopes
    
    def predict_tcell_epitopes(self, protein_record, peptide_length=9, threshold=0.5):
        """
        Predict T-cell epitopes using a simple scoring method
        
        Args:
            protein_record: SeqRecord object with protein sequence
            peptide_length: Length of peptides to consider
            threshold: Threshold for epitope prediction
            
        Returns:
            list: List of predicted epitopes with scores
        """
        # This is a simplified implementation
        # In a real implementation, we would use MHC binding prediction tools
        
        # Amino acid frequencies in known T-cell epitopes (fictional values for demonstration)
        aa_weights = {
            'A': 0.8, 'C': 0.3, 'D': 0.4, 'E': 0.5, 'F': 0.9,
            'G': 0.3, 'H': 0.6, 'I': 0.8, 'K': 0.6, 'L': 0.9,
            'M': 0.7, 'N': 0.5, 'P': 0.4, 'Q': 0.5, 'R': 0.6,
            'S': 0.5, 'T': 0.5, 'V': 0.8, 'W': 0.7, 'Y': 0.8
        }
        
        sequence = str(protein_record.seq).upper()
        epitopes = []
        
        for i in range(len(sequence) - peptide_length + 1):
            peptide = sequence[i:i+peptide_length]
            
            # Calculate score based on amino acid weights
            score = 0
            for aa in peptide:
                if aa in aa_weights:
                    score += aa_weights[aa]
                else:
                    score += 0.5  # Default for unknown amino acids
            
            score = score / peptide_length
            
            # If score exceeds threshold, consider it a potential epitope
            if score > threshold:
                epitopes.append({
                    'start': i + 1,  # 1-based indexing
                    'end': i + peptide_length,
                    'sequence': peptide,
                    'score': score,
                    'type': 'T-cell'
                })
        
        logger.info(f"Predicted {len(epitopes)} T-cell epitopes")
        return epitopes
    
    def predict_epitopes(self, protein_record):
        """
        Predict both B-cell and T-cell epitopes
        
        Args:
            protein_record: SeqRecord object with protein sequence
            
        Returns:
            DataFrame: Combined epitope predictions
        """
        bcell_epitopes = self.predict_bcell_epitopes(protein_record)
        tcell_epitopes = self.predict_tcell_epitopes(protein_record)
        
        # Combine epitopes
        all_epitopes = bcell_epitopes + tcell_epitopes
        
        # Convert to DataFrame
        epitope_df = pd.DataFrame(all_epitopes)
        
        # Sort by score (descending)
        epitope_df = epitope_df.sort_values('score', ascending=False)
        
        # Save results
        self.epitope_results = epitope_df
        
        # Save to CSV
        results_file = self.results_path / f"{protein_record.id}_epitopes.csv"
        epitope_df.to_csv(results_file, index=False)
        logger.info(f"Saved epitope predictions to {results_file}")
        
        return epitope_df
    
    def get_top_epitopes(self, n=10, epitope_type=None):
        """
        Get the top N epitopes from the prediction results
        
        Args:
            n: Number of top epitopes to return
            epitope_type: Type of epitopes to filter ('B-cell', 'T-cell', or None for all)
            
        Returns:
            DataFrame: Top N epitopes
        """
        if self.epitope_results is None:
            logger.error("No epitope prediction results available")
            return None
        
        df = self.epitope_results
        
        if epitope_type:
            df = df[df['type'] == epitope_type]
        
        return df.head(n)
    
    def save_epitopes_fasta(self, output_file=None):
        """
        Save predicted epitopes to a FASTA file
        
        Args:
            output_file: Path to output file (default: auto-generated)
            
        Returns:
            str: Path to the saved file
        """
        if self.epitope_results is None:
            logger.error("No epitope prediction results available")
            return None
        
        if output_file is None:
            output_file = self.results_path / "predicted_epitopes.fasta"
        
        records = []
        
        for i, row in self.epitope_results.iterrows():
            record = SeqRecord(
                Seq(row['sequence']),
                id=f"epitope_{i+1}",
                description=f"{row['type']} epitope, score={row['score']:.3f}, position={row['start']}-{row['end']}"
            )
            records.append(record)
        
        SeqIO.write(records, output_file, "fasta")
        logger.info(f"Saved epitopes to FASTA file: {output_file}")
        
        return output_file

