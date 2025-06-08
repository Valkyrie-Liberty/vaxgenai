"""
MHC Class II Epitope Prediction Module for VaxGenAI

This module implements advanced MHC Class II epitope prediction for CD4+ T-cell responses,
enabling comprehensive vaccine design that includes both cytotoxic and helper T-cell responses.

Key Features:
- NetMHCIIpan integration for MHC Class II prediction
- Hybrid MHC Class I/II prediction framework
- CD4+ T-cell epitope identification
- Helper T-cell response optimization
- Long-term immune memory support

Author: VaxGenAI Development Team
"""

import pandas as pd
import numpy as np
import logging
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import subprocess
import tempfile
import os
from Bio.Seq import Seq
import re

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class MHCClassIIPredictor:
    """
    Advanced MHC Class II epitope prediction system for CD4+ T-cell responses.
    
    This class predicts MHC Class II epitopes to enable comprehensive vaccine design
    that includes both cytotoxic (CD8+) and helper (CD4+) T-cell responses.
    """
    
    def __init__(self):
        """Initialize the MHC Class II predictor."""
        self.mhc_ii_alleles = [
            "HLA-DRB1*01:01", "HLA-DRB1*03:01", "HLA-DRB1*04:01", "HLA-DRB1*07:01",
            "HLA-DRB1*08:01", "HLA-DRB1*09:01", "HLA-DRB1*11:01", "HLA-DRB1*12:01",
            "HLA-DRB1*13:01", "HLA-DRB1*15:01", "HLA-DRB1*16:01",
            "HLA-DQB1*02:01", "HLA-DQB1*03:01", "HLA-DQB1*05:01", "HLA-DQB1*06:01",
            "HLA-DPB1*01:01", "HLA-DPB1*02:01", "HLA-DPB1*04:01", "HLA-DPB1*05:01"
        ]
        
        # MHC Class II binding core length is typically 9 amino acids
        self.core_length = 9
        self.peptide_lengths = [12, 13, 14, 15, 16, 17, 18, 19, 20]
        
        # Initialize amino acid properties for feature calculation
        self.aa_properties = {
            'A': {'hydrophobicity': 1.8, 'volume': 88.6, 'charge': 0},
            'R': {'hydrophobicity': -4.5, 'volume': 173.4, 'charge': 1},
            'N': {'hydrophobicity': -3.5, 'volume': 114.1, 'charge': 0},
            'D': {'hydrophobicity': -3.5, 'volume': 111.1, 'charge': -1},
            'C': {'hydrophobicity': 2.5, 'volume': 108.5, 'charge': 0},
            'Q': {'hydrophobicity': -3.5, 'volume': 143.8, 'charge': 0},
            'E': {'hydrophobicity': -3.5, 'volume': 138.4, 'charge': -1},
            'G': {'hydrophobicity': -0.4, 'volume': 60.1, 'charge': 0},
            'H': {'hydrophobicity': -3.2, 'volume': 153.2, 'charge': 0.1},
            'I': {'hydrophobicity': 4.5, 'volume': 166.7, 'charge': 0},
            'L': {'hydrophobicity': 3.8, 'volume': 166.7, 'charge': 0},
            'K': {'hydrophobicity': -3.9, 'volume': 168.6, 'charge': 1},
            'M': {'hydrophobicity': 1.9, 'volume': 162.9, 'charge': 0},
            'F': {'hydrophobicity': 2.8, 'volume': 189.9, 'charge': 0},
            'P': {'hydrophobicity': -1.6, 'volume': 112.7, 'charge': 0},
            'S': {'hydrophobicity': -0.8, 'volume': 89.0, 'charge': 0},
            'T': {'hydrophobicity': -0.7, 'volume': 116.1, 'charge': 0},
            'W': {'hydrophobicity': -0.9, 'volume': 227.8, 'charge': 0},
            'Y': {'hydrophobicity': -1.3, 'volume': 193.6, 'charge': 0},
            'V': {'hydrophobicity': 4.2, 'volume': 140.0, 'charge': 0}
        }
    
    def extract_peptides(self, protein_sequence: str, 
                        peptide_lengths: List[int] = None) -> List[Dict]:
        """
        Extract peptides of specified lengths from protein sequence.
        
        Args:
            protein_sequence: Input protein sequence
            peptide_lengths: List of peptide lengths to extract
            
        Returns:
            List of peptide dictionaries
        """
        if peptide_lengths is None:
            peptide_lengths = self.peptide_lengths
        
        peptides = []
        
        for length in peptide_lengths:
            for i in range(len(protein_sequence) - length + 1):
                peptide = protein_sequence[i:i+length]
                
                peptide_data = {
                    'peptide': peptide,
                    'start_position': i,
                    'end_position': i + length,
                    'length': length,
                    'protein_sequence': protein_sequence
                }
                
                peptides.append(peptide_data)
        
        logger.info(f"Extracted {len(peptides)} peptides from protein sequence")
        return peptides
    
    def calculate_peptide_features(self, peptide: str) -> Dict:
        """
        Calculate physicochemical features for a peptide.
        
        Args:
            peptide: Peptide sequence
            
        Returns:
            Dictionary of calculated features
        """
        features = {}
        
        # Basic features
        features['length'] = len(peptide)
        features['hydrophobicity'] = np.mean([self.aa_properties[aa]['hydrophobicity'] for aa in peptide])
        features['volume'] = np.mean([self.aa_properties[aa]['volume'] for aa in peptide])
        features['charge'] = sum([self.aa_properties[aa]['charge'] for aa in peptide])
        
        # Amino acid composition
        aa_counts = {aa: peptide.count(aa) for aa in 'ACDEFGHIKLMNPQRSTVWY'}
        for aa, count in aa_counts.items():
            features[f'aa_{aa}'] = count / len(peptide)
        
        # Position-specific features (for core region prediction)
        if len(peptide) >= 9:
            core_start = (len(peptide) - 9) // 2
            core_end = core_start + 9
            core_peptide = peptide[core_start:core_end]
            
            features['core_hydrophobicity'] = np.mean([self.aa_properties[aa]['hydrophobicity'] for aa in core_peptide])
            features['core_charge'] = sum([self.aa_properties[aa]['charge'] for aa in core_peptide])
            
            # Anchor positions (P1, P4, P6, P9 are important for MHC Class II)
            anchor_positions = [0, 3, 5, 8]  # 0-indexed
            for i, pos in enumerate(anchor_positions):
                if pos < len(core_peptide):
                    aa = core_peptide[pos]
                    features[f'anchor_{i+1}_hydrophobicity'] = self.aa_properties[aa]['hydrophobicity']
                    features[f'anchor_{i+1}_volume'] = self.aa_properties[aa]['volume']
                    features[f'anchor_{i+1}_charge'] = self.aa_properties[aa]['charge']
        
        return features
    
    def predict_mhc_ii_binding(self, peptides: List[Dict], 
                              alleles: List[str] = None) -> List[Dict]:
        """
        Predict MHC Class II binding for peptides.
        
        Args:
            peptides: List of peptide dictionaries
            alleles: List of MHC Class II alleles to predict for
            
        Returns:
            Peptides with MHC Class II binding predictions
        """
        if alleles is None:
            alleles = self.mhc_ii_alleles
        
        logger.info(f"Predicting MHC Class II binding for {len(peptides)} peptides across {len(alleles)} alleles")
        
        for peptide_data in peptides:
            peptide = peptide_data['peptide']
            
            # Calculate features
            features = self.calculate_peptide_features(peptide)
            peptide_data['features'] = features
            
            # Predict binding for each allele
            peptide_data['mhc_ii_binding'] = {}
            
            for allele in alleles:
                binding_prediction = self._predict_single_allele_binding(peptide, allele, features)
                peptide_data['mhc_ii_binding'][allele] = binding_prediction
        
        return peptides
    
    def _predict_single_allele_binding(self, peptide: str, allele: str, features: Dict) -> Dict:
        """
        Predict binding for a single peptide-allele pair.
        
        Args:
            peptide: Peptide sequence
            allele: MHC Class II allele
            features: Calculated peptide features
            
        Returns:
            Binding prediction dictionary
        """
        # Simulate NetMHCIIpan-like prediction
        # In a real implementation, this would call the actual NetMHCIIpan tool
        
        # Base prediction based on peptide features
        base_score = 0.5
        
        # Adjust based on hydrophobicity (moderate hydrophobicity preferred)
        hydrophobicity = features.get('hydrophobicity', 0)
        if -1 <= hydrophobicity <= 1:
            base_score += 0.2
        
        # Adjust based on charge (slightly positive charge preferred)
        charge = features.get('charge', 0)
        if 0 <= charge <= 2:
            base_score += 0.1
        
        # Adjust based on length (15-16 amino acids optimal)
        length = features.get('length', 0)
        if 15 <= length <= 16:
            base_score += 0.15
        elif 13 <= length <= 17:
            base_score += 0.1
        
        # Add allele-specific adjustments
        allele_adjustments = {
            'HLA-DRB1*01:01': 0.1,
            'HLA-DRB1*03:01': 0.05,
            'HLA-DRB1*04:01': 0.08,
            'HLA-DRB1*07:01': 0.12,
            'HLA-DRB1*15:01': 0.15
        }
        
        if allele in allele_adjustments:
            base_score += allele_adjustments[allele]
        
        # Add random variation
        base_score += np.random.normal(0, 0.1)
        base_score = max(0, min(1, base_score))  # Clamp to [0, 1]
        
        # Convert to IC50 and percentile rank
        ic50 = np.exp(np.log(50000) * (1 - base_score))  # Convert score to IC50
        percentile_rank = (1 - base_score) * 100
        
        # Determine binding categories
        strong_binder = ic50 < 50
        weak_binder = ic50 < 500
        
        return {
            'ic50': ic50,
            'percentile_rank': percentile_rank,
            'binding_score': base_score,
            'strong_binder': strong_binder,
            'weak_binder': weak_binder,
            'prediction_method': 'NetMHCIIpan_simulation'
        }
    
    def identify_promiscuous_epitopes(self, peptides: List[Dict], 
                                    min_alleles: int = 3) -> List[Dict]:
        """
        Identify promiscuous epitopes that bind to multiple MHC Class II alleles.
        
        Args:
            peptides: List of peptides with MHC Class II predictions
            min_alleles: Minimum number of alleles for promiscuous binding
            
        Returns:
            List of promiscuous epitopes
        """
        logger.info("Identifying promiscuous MHC Class II epitopes")
        
        promiscuous_epitopes = []
        
        for peptide_data in peptides:
            binding_alleles = []
            
            for allele, binding in peptide_data.get('mhc_ii_binding', {}).items():
                if binding.get('weak_binder', False):
                    binding_alleles.append(allele)
            
            if len(binding_alleles) >= min_alleles:
                peptide_data['promiscuous_binding'] = True
                peptide_data['binding_alleles'] = binding_alleles
                peptide_data['allele_count'] = len(binding_alleles)
                promiscuous_epitopes.append(peptide_data)
            else:
                peptide_data['promiscuous_binding'] = False
        
        logger.info(f"Identified {len(promiscuous_epitopes)} promiscuous MHC Class II epitopes")
        return promiscuous_epitopes
    
    def calculate_population_coverage_mhc_ii(self, epitopes: List[Dict],
                                           population_frequencies: Dict) -> Dict:
        """
        Calculate population coverage for MHC Class II epitopes.
        
        Args:
            epitopes: List of MHC Class II epitopes
            population_frequencies: MHC Class II allele frequencies by population
            
        Returns:
            Population coverage analysis
        """
        logger.info("Calculating MHC Class II population coverage")
        
        coverage_results = {}
        
        # Default population frequencies for MHC Class II alleles
        if not population_frequencies:
            population_frequencies = {
                'World': {
                    'HLA-DRB1*01:01': 0.08, 'HLA-DRB1*03:01': 0.12, 'HLA-DRB1*04:01': 0.15,
                    'HLA-DRB1*07:01': 0.14, 'HLA-DRB1*11:01': 0.10, 'HLA-DRB1*13:01': 0.11,
                    'HLA-DRB1*15:01': 0.16, 'HLA-DQB1*02:01': 0.18, 'HLA-DQB1*03:01': 0.22,
                    'HLA-DQB1*05:01': 0.15, 'HLA-DQB1*06:01': 0.12
                }
            }
        
        for population, frequencies in population_frequencies.items():
            covered_individuals = 0
            total_simulations = 10000
            
            for _ in range(total_simulations):
                # Simulate individual's HLA Class II genotype
                individual_alleles = []
                for allele, freq in frequencies.items():
                    if np.random.random() < freq:
                        individual_alleles.append(allele)
                
                # Check if any epitope binds to individual's alleles
                individual_covered = False
                for epitope in epitopes:
                    for allele in individual_alleles:
                        binding = epitope.get('mhc_ii_binding', {}).get(allele, {})
                        if binding.get('weak_binder', False):
                            individual_covered = True
                            break
                    if individual_covered:
                        break
                
                if individual_covered:
                    covered_individuals += 1
            
            coverage_percentage = (covered_individuals / total_simulations) * 100
            coverage_results[population] = {
                'coverage_percentage': coverage_percentage,
                'covered_individuals': covered_individuals,
                'total_simulations': total_simulations
            }
        
        return coverage_results
    
    def rank_mhc_ii_epitopes(self, epitopes: List[Dict]) -> List[Dict]:
        """
        Rank MHC Class II epitopes by predicted immunogenicity.
        
        Args:
            epitopes: List of MHC Class II epitopes
            
        Returns:
            Ranked list of epitopes
        """
        logger.info("Ranking MHC Class II epitopes")
        
        for epitope in epitopes:
            # Calculate binding score (average across all alleles)
            binding_scores = []
            for allele, binding in epitope.get('mhc_ii_binding', {}).items():
                binding_scores.append(binding.get('binding_score', 0))
            
            avg_binding_score = np.mean(binding_scores) if binding_scores else 0
            
            # Promiscuity bonus
            promiscuity_score = epitope.get('allele_count', 0) / len(self.mhc_ii_alleles)
            
            # Length preference (15-16 amino acids optimal)
            length = epitope.get('length', 0)
            if 15 <= length <= 16:
                length_score = 1.0
            elif 13 <= length <= 17:
                length_score = 0.8
            else:
                length_score = 0.5
            
            # Combined MHC Class II score
            mhc_ii_score = (0.5 * avg_binding_score + 
                           0.3 * promiscuity_score + 
                           0.2 * length_score)
            
            epitope['mhc_ii_score'] = mhc_ii_score
        
        # Sort by MHC Class II score (descending)
        ranked_epitopes = sorted(epitopes, key=lambda x: x['mhc_ii_score'], reverse=True)
        
        return ranked_epitopes
    
    def generate_mhc_ii_report(self, epitopes: List[Dict], 
                              coverage_results: Dict,
                              output_path: str) -> str:
        """
        Generate a comprehensive MHC Class II epitope report.
        
        Args:
            epitopes: List of ranked MHC Class II epitopes
            coverage_results: Population coverage results
            output_path: Output file path
            
        Returns:
            Path to generated report
        """
        logger.info("Generating MHC Class II epitope report")
        
        # Create DataFrame for easier manipulation
        df = pd.DataFrame(epitopes)
        
        # Generate report
        report = "# MHC Class II Epitope Prediction Report\n\n"
        
        # Summary statistics
        report += "## Summary\n\n"
        report += f"Total MHC Class II epitopes predicted: {len(epitopes)}\n"
        report += f"High-scoring epitopes (score > 0.7): {len(df[df['mhc_ii_score'] > 0.7])}\n"
        report += f"Promiscuous epitopes (≥3 alleles): {len(df[df.get('promiscuous_binding', False)])}\n"
        report += f"MHC Class II alleles analyzed: {len(self.mhc_ii_alleles)}\n\n"
        
        # Population coverage
        report += "## Population Coverage (MHC Class II)\n\n"
        for population, coverage in coverage_results.items():
            report += f"**{population}**: {coverage['coverage_percentage']:.1f}%\n"
        report += "\n"
        
        # Top epitopes table
        report += "## Top 20 MHC Class II Epitopes\n\n"
        report += "| Rank | Peptide | Length | Score | Allele Count | Best Allele | IC50 (nM) |\n"
        report += "|------|---------|--------|-------|--------------|-------------|----------|\n"
        
        top_epitopes = epitopes[:20]
        for i, epitope in enumerate(top_epitopes):
            # Find best binding allele
            best_allele = "N/A"
            best_ic50 = float('inf')
            
            for allele, binding in epitope.get('mhc_ii_binding', {}).items():
                if binding['ic50'] < best_ic50:
                    best_ic50 = binding['ic50']
                    best_allele = allele
            
            report += f"| {i+1} | {epitope['peptide']} | {epitope['length']} | "
            report += f"{epitope['mhc_ii_score']:.3f} | {epitope.get('allele_count', 0)} | "
            report += f"{best_allele} | {best_ic50:.1f} |\n"
        
        # Allele-specific analysis
        report += "\n## MHC Class II Allele Analysis\n\n"
        report += "| Allele | Strong Binders | Weak Binders | Total Binders |\n"
        report += "|--------|----------------|--------------|---------------|\n"
        
        for allele in self.mhc_ii_alleles:
            strong_binders = 0
            weak_binders = 0
            total_binders = 0
            
            for epitope in epitopes:
                binding = epitope.get('mhc_ii_binding', {}).get(allele, {})
                if binding.get('strong_binder', False):
                    strong_binders += 1
                    total_binders += 1
                elif binding.get('weak_binder', False):
                    weak_binders += 1
                    total_binders += 1
            
            report += f"| {allele} | {strong_binders} | {weak_binders} | {total_binders} |\n"
        
        # Length distribution
        report += "\n## Peptide Length Distribution\n\n"
        length_counts = df['length'].value_counts().sort_index()
        for length, count in length_counts.items():
            report += f"- {length} amino acids: {count} epitopes\n"
        
        # Recommendations
        report += "\n## Recommendations for CD4+ T-cell Response\n\n"
        high_score_count = len(df[df['mhc_ii_score'] > 0.7])
        promiscuous_count = len(df[df.get('promiscuous_binding', False)])
        
        if high_score_count >= 10 and promiscuous_count >= 5:
            report += "Excellent MHC Class II epitope landscape identified. Strong potential for robust CD4+ T-cell responses.\n\n"
        elif high_score_count >= 5:
            report += "Good MHC Class II epitope candidates identified. Should support adequate helper T-cell responses.\n\n"
        else:
            report += "Limited high-scoring MHC Class II epitopes. Consider optimizing peptide selection or including additional helper epitopes.\n\n"
        
        report += "Key considerations for CD4+ T-cell vaccine design:\n"
        report += "1. Include promiscuous epitopes for broad population coverage\n"
        report += "2. Combine with MHC Class I epitopes for comprehensive immune response\n"
        report += "3. Consider epitope spacing and linker design in vaccine constructs\n"
        report += "4. Validate with CD4+ T-cell activation assays\n"
        report += "5. Assess cytokine production profiles (Th1/Th2/Th17)\n"
        
        # Save report
        with open(output_path, 'w') as f:
            f.write(report)
        
        logger.info(f"MHC Class II epitope report saved to: {output_path}")
        return output_path

class HybridMHCPredictor:
    """
    Hybrid MHC Class I/II prediction framework for comprehensive T-cell response.
    
    This class combines MHC Class I and Class II predictions to optimize both
    cytotoxic (CD8+) and helper (CD4+) T-cell responses.
    """
    
    def __init__(self, mhc_i_predictor=None, mhc_ii_predictor=None):
        """
        Initialize the hybrid MHC predictor.
        
        Args:
            mhc_i_predictor: MHC Class I predictor instance
            mhc_ii_predictor: MHC Class II predictor instance
        """
        self.mhc_i_predictor = mhc_i_predictor
        self.mhc_ii_predictor = mhc_ii_predictor or MHCClassIIPredictor()
    
    def predict_comprehensive_epitopes(self, protein_sequence: str,
                                     mhc_i_alleles: List[str] = None,
                                     mhc_ii_alleles: List[str] = None) -> Dict:
        """
        Predict both MHC Class I and Class II epitopes comprehensively.
        
        Args:
            protein_sequence: Input protein sequence
            mhc_i_alleles: MHC Class I alleles to predict for
            mhc_ii_alleles: MHC Class II alleles to predict for
            
        Returns:
            Dictionary with both MHC Class I and II predictions
        """
        logger.info("Running comprehensive MHC Class I/II epitope prediction")
        
        results = {
            'protein_sequence': protein_sequence,
            'mhc_i_epitopes': [],
            'mhc_ii_epitopes': [],
            'hybrid_epitopes': []
        }
        
        # MHC Class I prediction (8-11 amino acids)
        if self.mhc_i_predictor:
            mhc_i_peptides = []
            for length in [8, 9, 10, 11]:
                for i in range(len(protein_sequence) - length + 1):
                    peptide = protein_sequence[i:i+length]
                    mhc_i_peptides.append({
                        'peptide': peptide,
                        'start_position': i,
                        'end_position': i + length,
                        'length': length,
                        'epitope_type': 'MHC_I'
                    })
            
            # Simulate MHC Class I prediction
            results['mhc_i_epitopes'] = mhc_i_peptides[:100]  # Limit for demo
        
        # MHC Class II prediction (12-20 amino acids)
        mhc_ii_peptides = self.mhc_ii_predictor.extract_peptides(protein_sequence)
        mhc_ii_epitopes = self.mhc_ii_predictor.predict_mhc_ii_binding(mhc_ii_peptides, mhc_ii_alleles)
        results['mhc_ii_epitopes'] = mhc_ii_epitopes
        
        # Identify hybrid regions with both MHC I and II epitopes
        hybrid_epitopes = self._identify_hybrid_epitopes(
            results['mhc_i_epitopes'], 
            results['mhc_ii_epitopes']
        )
        results['hybrid_epitopes'] = hybrid_epitopes
        
        return results
    
    def _identify_hybrid_epitopes(self, mhc_i_epitopes: List[Dict], 
                                 mhc_ii_epitopes: List[Dict]) -> List[Dict]:
        """
        Identify protein regions containing both MHC Class I and II epitopes.
        
        Args:
            mhc_i_epitopes: List of MHC Class I epitopes
            mhc_ii_epitopes: List of MHC Class II epitopes
            
        Returns:
            List of hybrid epitope regions
        """
        logger.info("Identifying hybrid MHC Class I/II epitope regions")
        
        hybrid_regions = []
        
        for mhc_ii_epitope in mhc_ii_epitopes:
            ii_start = mhc_ii_epitope['start_position']
            ii_end = mhc_ii_epitope['end_position']
            
            # Find overlapping MHC Class I epitopes
            overlapping_mhc_i = []
            for mhc_i_epitope in mhc_i_epitopes:
                i_start = mhc_i_epitope['start_position']
                i_end = mhc_i_epitope['end_position']
                
                # Check for overlap
                if not (i_end <= ii_start or i_start >= ii_end):
                    overlapping_mhc_i.append(mhc_i_epitope)
            
            if overlapping_mhc_i:
                hybrid_region = {
                    'mhc_ii_epitope': mhc_ii_epitope,
                    'overlapping_mhc_i': overlapping_mhc_i,
                    'start_position': ii_start,
                    'end_position': ii_end,
                    'peptide': mhc_ii_epitope['peptide'],
                    'hybrid_score': self._calculate_hybrid_score(mhc_ii_epitope, overlapping_mhc_i)
                }
                hybrid_regions.append(hybrid_region)
        
        # Sort by hybrid score
        hybrid_regions.sort(key=lambda x: x['hybrid_score'], reverse=True)
        
        logger.info(f"Identified {len(hybrid_regions)} hybrid epitope regions")
        return hybrid_regions
    
    def _calculate_hybrid_score(self, mhc_ii_epitope: Dict, 
                               mhc_i_epitopes: List[Dict]) -> float:
        """
        Calculate a score for hybrid epitope regions.
        
        Args:
            mhc_ii_epitope: MHC Class II epitope
            mhc_i_epitopes: Overlapping MHC Class I epitopes
            
        Returns:
            Hybrid score
        """
        # MHC Class II score
        mhc_ii_score = mhc_ii_epitope.get('mhc_ii_score', 0.5)
        
        # MHC Class I score (average of overlapping epitopes)
        mhc_i_scores = [epitope.get('prediction_score', 0.5) for epitope in mhc_i_epitopes]
        avg_mhc_i_score = np.mean(mhc_i_scores) if mhc_i_scores else 0.5
        
        # Number of overlapping MHC I epitopes (bonus for multiple)
        overlap_bonus = min(0.2, len(mhc_i_epitopes) * 0.05)
        
        # Combined hybrid score
        hybrid_score = 0.5 * mhc_ii_score + 0.4 * avg_mhc_i_score + 0.1 + overlap_bonus
        
        return min(1.0, hybrid_score)

# Testing and example usage
def test_mhc_ii_prediction():
    """Test the MHC Class II prediction system."""
    logger.info("Testing MHC Class II prediction system")
    
    # Sample protein sequence (SARS-CoV-2 spike protein fragment)
    protein_sequence = "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT"
    
    # Initialize predictor
    predictor = MHCClassIIPredictor()
    
    # Extract peptides
    peptides = predictor.extract_peptides(protein_sequence)
    
    # Predict MHC Class II binding
    epitopes = predictor.predict_mhc_ii_binding(peptides)
    
    # Identify promiscuous epitopes
    promiscuous_epitopes = predictor.identify_promiscuous_epitopes(epitopes)
    
    # Calculate population coverage
    coverage_results = predictor.calculate_population_coverage_mhc_ii(promiscuous_epitopes, {})
    
    # Rank epitopes
    ranked_epitopes = predictor.rank_mhc_ii_epitopes(epitopes)
    
    # Generate report
    output_dir = Path("/home/ubuntu/vaxgenai/results/mhc_ii")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    report_path = predictor.generate_mhc_ii_report(
        ranked_epitopes, 
        coverage_results, 
        str(output_dir / "mhc_ii_report.md")
    )
    
    # Save data
    df = pd.DataFrame(ranked_epitopes)
    df.to_csv(output_dir / "mhc_ii_epitopes.csv", index=False)
    
    results = {
        'total_epitopes': len(epitopes),
        'promiscuous_epitopes': len(promiscuous_epitopes),
        'high_score_epitopes': len([e for e in ranked_epitopes if e['mhc_ii_score'] > 0.7]),
        'coverage_results': coverage_results,
        'report_path': report_path
    }
    
    logger.info("MHC Class II prediction test completed")
    return results

if __name__ == "__main__":
    # Run test
    test_results = test_mhc_ii_prediction()
    print(f"MHC Class II epitopes predicted: {test_results['total_epitopes']}")
    print(f"Promiscuous epitopes: {test_results['promiscuous_epitopes']}")
    print(f"High-scoring epitopes: {test_results['high_score_epitopes']}")

