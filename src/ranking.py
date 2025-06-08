"""
Ranking Module for VaxGenAI

This module ranks vaccine candidates based on predicted efficacy, population coverage,
and manufacturability.
"""

import os
import logging
import numpy as np
import pandas as pd
from pathlib import Path

from .utils import get_results_path

logger = logging.getLogger("vaxgenai.ranking")

class RankingEngine:
    """
    Ranks vaccine candidates based on multiple criteria
    """
    
    def __init__(self):
        """Initialize the ranking engine"""
        self.results_path = get_results_path()
        logger.info("Ranking engine initialized")
    
    def calculate_efficacy_score(self, candidate):
        """
        Calculate efficacy score for a vaccine candidate
        
        Args:
            candidate: Vaccine candidate
            
        Returns:
            float: Efficacy score (0-1)
        """
        # This is a simplified implementation
        # In a real implementation, we would use more sophisticated methods
        
        if candidate['type'] == 'subunit':
            # For subunit vaccines, consider the number and quality of epitopes
            num_epitopes = candidate['num_epitopes']
            avg_epitope_score = np.mean([e['score'] for e in candidate['epitopes']])
            
            # Calculate efficacy score
            efficacy_score = 0.5 * (num_epitopes / 20) + 0.5 * avg_epitope_score
            
        elif candidate['type'] == 'mRNA':
            # For mRNA vaccines, consider sequence properties
            gc_content = candidate['gc_content']
            length = candidate['length']
            
            # Optimal GC content is around 0.5-0.6
            gc_score = 1.0 - abs(gc_content - 0.55) * 2
            
            # Length score (prefer shorter sequences for manufacturing)
            length_score = 1.0 - min(1.0, length / 5000)
            
            # Calculate efficacy score
            efficacy_score = 0.7 * gc_score + 0.3 * length_score
            
        elif candidate['type'] == 'peptide':
            # For peptide vaccines, consider the average epitope score
            avg_epitope_score = np.mean([e['score'] for e in candidate['epitopes']])
            
            # Calculate efficacy score
            efficacy_score = avg_epitope_score
            
        else:
            logger.warning(f"Unknown vaccine type: {candidate['type']}")
            efficacy_score = 0.5
        
        # Ensure score is between 0 and 1
        efficacy_score = max(0.0, min(1.0, efficacy_score))
        
        logger.info(f"Calculated efficacy score for {candidate['type']} vaccine: {efficacy_score:.3f}")
        return efficacy_score
    
    def calculate_population_coverage(self, candidate):
        """
        Calculate population coverage for a vaccine candidate
        
        Args:
            candidate: Vaccine candidate
            
        Returns:
            float: Population coverage score (0-1)
        """
        # This is a simplified implementation
        # In a real implementation, we would use MHC binding prediction across populations
        
        # For demonstration, we'll use a random score
        coverage_score = np.random.uniform(0.6, 0.9)
        
        logger.info(f"Calculated population coverage for {candidate['type']} vaccine: {coverage_score:.3f}")
        return coverage_score
    
    def calculate_manufacturability(self, candidate):
        """
        Calculate manufacturability score for a vaccine candidate
        
        Args:
            candidate: Vaccine candidate
            
        Returns:
            float: Manufacturability score (0-1)
        """
        if candidate['type'] == 'subunit':
            # For subunit vaccines, consider sequence length and complexity
            length = candidate['length']
            
            # Length score (prefer moderate length)
            length_score = 1.0 - abs(length - 300) / 300
            
            # Calculate manufacturability score
            manufacturability_score = length_score
            
        elif candidate['type'] == 'mRNA':
            # For mRNA vaccines, consider sequence properties
            gc_content = candidate['gc_content']
            length = candidate['length']
            
            # GC content score (prefer moderate GC content)
            gc_score = 1.0 - abs(gc_content - 0.55) * 2
            
            # Length score (prefer shorter sequences)
            length_score = 1.0 - min(1.0, length / 5000)
            
            # Calculate manufacturability score
            manufacturability_score = 0.5 * gc_score + 0.5 * length_score
            
        elif candidate['type'] == 'peptide':
            # For peptide vaccines, consider the number of epitopes
            num_epitopes = candidate['num_epitopes']
            
            # Number of epitopes score (prefer fewer epitopes for manufacturing)
            epitope_score = 1.0 - min(1.0, num_epitopes / 20)
            
            # Calculate manufacturability score
            manufacturability_score = epitope_score
            
        else:
            logger.warning(f"Unknown vaccine type: {candidate['type']}")
            manufacturability_score = 0.5
        
        # Ensure score is between 0 and 1
        manufacturability_score = max(0.0, min(1.0, manufacturability_score))
        
        logger.info(f"Calculated manufacturability for {candidate['type']} vaccine: {manufacturability_score:.3f}")
        return manufacturability_score
    
    def rank_candidates(self, candidates, weights=None):
        """
        Rank vaccine candidates based on multiple criteria
        
        Args:
            candidates: List of vaccine candidates
            weights: Dictionary of weights for each criterion
            
        Returns:
            DataFrame: Ranked candidates with scores
        """
        if not candidates:
            logger.error("No candidates to rank")
            return None
        
        # Default weights
        if weights is None:
            weights = {
                'efficacy': 0.5,
                'population_coverage': 0.3,
                'manufacturability': 0.2
            }
        
        # Calculate scores for each candidate
        ranked_candidates = []
        
        for candidate in candidates:
            # Calculate scores
            efficacy_score = self.calculate_efficacy_score(candidate)
            population_coverage = self.calculate_population_coverage(candidate)
            manufacturability = self.calculate_manufacturability(candidate)
            
            # Calculate overall score
            overall_score = (
                weights['efficacy'] * efficacy_score +
                weights['population_coverage'] * population_coverage +
                weights['manufacturability'] * manufacturability
            )
            
            # Create ranked candidate
            ranked_candidate = {
                'type': candidate['type'],
                'efficacy_score': efficacy_score,
                'population_coverage': population_coverage,
                'manufacturability': manufacturability,
                'overall_score': overall_score
            }
            
            # Add additional information based on vaccine type
            if candidate['type'] == 'subunit':
                ranked_candidate['length'] = candidate['length']
                ranked_candidate['num_epitopes'] = candidate['num_epitopes']
            elif candidate['type'] == 'mRNA':
                ranked_candidate['length'] = candidate['length']
                ranked_candidate['gc_content'] = candidate['gc_content']
            elif candidate['type'] == 'peptide':
                ranked_candidate['num_epitopes'] = candidate['num_epitopes']
            
            # Add safety information if available
            if 'safety_results' in candidate:
                if candidate['type'] == 'peptide':
                    ranked_candidate['is_safe'] = candidate['safety_results'].get('all_epitopes_safe', False)
                else:
                    ranked_candidate['is_safe'] = candidate['safety_results'].get('is_safe', False)
            
            ranked_candidates.append(ranked_candidate)
        
        # Convert to DataFrame
        df = pd.DataFrame(ranked_candidates)
        
        # Sort by overall score (descending)
        df = df.sort_values('overall_score', ascending=False)
        
        # Save to CSV
        output_file = self.results_path / "ranked_candidates.csv"
        df.to_csv(output_file, index=False)
        logger.info(f"Saved ranked candidates to {output_file}")
        
        return df
    
    def generate_report(self, ranked_candidates):
        """
        Generate a report for the ranked vaccine candidates
        
        Args:
            ranked_candidates: DataFrame of ranked candidates
            
        Returns:
            str: Path to the report file
        """
        # Create report
        report = []
        
        # Add header
        report.append("# Vaccine Candidate Ranking Report")
        report.append("")
        report.append("## Summary")
        report.append("")
        report.append(f"Total candidates evaluated: {len(ranked_candidates)}")
        report.append(f"Top candidate type: {ranked_candidates.iloc[0]['type']}")
        report.append(f"Top candidate overall score: {ranked_candidates.iloc[0]['overall_score']:.3f}")
        report.append("")
        
        # Add table of top candidates
        report.append("## Top Candidates")
        report.append("")
        report.append("| Rank | Type | Efficacy | Population Coverage | Manufacturability | Overall Score |")
        report.append("|------|------|----------|---------------------|-------------------|---------------|")
        
        for i, (_, candidate) in enumerate(ranked_candidates.head(5).iterrows()):
            report.append(f"| {i+1} | {candidate['type']} | {candidate['efficacy_score']:.3f} | {candidate['population_coverage']:.3f} | {candidate['manufacturability']:.3f} | {candidate['overall_score']:.3f} |")
        
        report.append("")
        
        # Add details for each vaccine type
        report.append("## Details by Vaccine Type")
        report.append("")
        
        for vaccine_type in ranked_candidates['type'].unique():
            report.append(f"### {vaccine_type.capitalize()} Vaccines")
            report.append("")
            
            type_candidates = ranked_candidates[ranked_candidates['type'] == vaccine_type]
            
            report.append(f"Number of candidates: {len(type_candidates)}")
            report.append(f"Average overall score: {type_candidates['overall_score'].mean():.3f}")
            report.append("")
            
            # Add specific details based on vaccine type
            if vaccine_type == 'subunit':
                report.append(f"Average length: {type_candidates['length'].mean():.1f} amino acids")
                report.append(f"Average number of epitopes: {type_candidates['num_epitopes'].mean():.1f}")
            elif vaccine_type == 'mRNA':
                report.append(f"Average length: {type_candidates['length'].mean():.1f} nucleotides")
                report.append(f"Average GC content: {type_candidates['gc_content'].mean():.3f}")
            elif vaccine_type == 'peptide':
                report.append(f"Average number of epitopes: {type_candidates['num_epitopes'].mean():.1f}")
            
            report.append("")
        
        # Add recommendations
        report.append("## Recommendations")
        report.append("")
        
        top_candidate = ranked_candidates.iloc[0]
        
        report.append(f"The top-ranked candidate is a {top_candidate['type']} vaccine with an overall score of {top_candidate['overall_score']:.3f}.")
        report.append("This candidate offers the best balance of efficacy, population coverage, and manufacturability.")
        report.append("")
        
        if 'is_safe' in top_candidate and top_candidate['is_safe']:
            report.append("This candidate has passed all safety checks and is ready for further evaluation.")
        else:
            report.append("Further safety evaluation is recommended before proceeding with this candidate.")
        
        report.append("")
        
        # Save report
        report_file = self.results_path / "vaccine_ranking_report.md"
        with open(report_file, 'w') as f:
            f.write('\n'.join(report))
        
        logger.info(f"Generated ranking report: {report_file}")
        
        return str(report_file)

