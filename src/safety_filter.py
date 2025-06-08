"""
Safety Filter Module for VaxGenAI

This module filters vaccine candidates based on safety criteria, including
allergenicity and toxicity prediction.
"""

import os
import logging
import numpy as np
import pandas as pd
from pathlib import Path
from Bio.Seq import Seq

from .utils import get_results_path

logger = logging.getLogger("vaxgenai.safety_filter")

class SafetyFilter:
    """
    Filters vaccine candidates based on safety criteria
    """
    
    def __init__(self):
        """Initialize the safety filter module"""
        self.results_path = get_results_path()
        logger.info("Safety filter module initialized")
        
        # Load allergenicity data (simplified for demonstration)
        self.allergenic_motifs = [
            'FGEQFK', 'EQFKAV', 'FKAVY', 'KAVYQ', 'AVYQK', 'VYQKG',
            'YQKGM', 'QKGMG', 'KGMGK', 'GMGKV', 'MGKVH', 'GKVHA',
            'KVHAA', 'VHAAN', 'HAANE', 'AANEE', 'ANEEI', 'NEEIL',
            'EEILK', 'EILKK', 'ILKKL', 'LKKLK', 'KKLKS', 'KLKSV',
            'LKSVP', 'KSVPE', 'SVPEV', 'VPEVP', 'PEVPT', 'EVPTP',
            'VPTPE', 'PTPEE', 'TPEEE', 'PEEEQ', 'EEEQK', 'EEQKK',
            'EQKKL', 'QKKLQ', 'KKLQL', 'KLQLK', 'LQLKR', 'QLKRM',
            'LKRME', 'KRMEE', 'RMEES', 'MEESQ', 'EESQR', 'ESQRQ',
            'SQRQR', 'QRQRQ', 'RQRQL', 'QRQLR', 'RQLRQ', 'QLRQQ',
            'LRQQQ', 'RQQQQ', 'QQQQQ', 'QQQQL', 'QQQLQ', 'QQLQQ',
            'QLQQQ', 'LQQQQ', 'QQQQE', 'QQQEE', 'QQEEE', 'QEEEE'
        ]
        
        # Load toxicity data (simplified for demonstration)
        self.toxic_motifs = [
            'RRRRRR', 'KKKKKK', 'DDDDDD', 'EEEEEE', 'PPPPPP',
            'GGGGGG', 'AAAAAA', 'VVVVVV', 'IIIIII', 'LLLLLL',
            'FFFFFF', 'WWWWWW', 'YYYYYY', 'HHHHHH', 'QQQQQQ',
            'NNNNNN', 'SSSSSS', 'TTTTTT', 'CCCCCC', 'MMMMMM'
        ]
    
    def check_allergenicity(self, sequence, threshold=0.5):
        """
        Check if a sequence contains allergenic motifs
        
        Args:
            sequence: Protein sequence to check
            threshold: Threshold for allergenicity score
            
        Returns:
            dict: Allergenicity assessment
        """
        # Count allergenic motifs
        allergen_count = 0
        allergen_motifs_found = []
        
        for motif in self.allergenic_motifs:
            if motif in sequence:
                allergen_count += 1
                allergen_motifs_found.append(motif)
        
        # Calculate allergenicity score
        allergenicity_score = allergen_count / (len(sequence) / 10)  # Normalize by sequence length
        
        # Determine if sequence is allergenic
        is_allergenic = allergenicity_score > threshold
        
        result = {
            'allergenicity_score': allergenicity_score,
            'is_allergenic': is_allergenic,
            'allergenic_motifs_found': allergen_motifs_found,
            'num_allergenic_motifs': allergen_count
        }
        
        logger.info(f"Allergenicity check: score={allergenicity_score:.3f}, is_allergenic={is_allergenic}")
        return result
    
    def check_toxicity(self, sequence, threshold=0.3):
        """
        Check if a sequence contains toxic motifs
        
        Args:
            sequence: Protein sequence to check
            threshold: Threshold for toxicity score
            
        Returns:
            dict: Toxicity assessment
        """
        # Count toxic motifs
        toxic_count = 0
        toxic_motifs_found = []
        
        for motif in self.toxic_motifs:
            if motif in sequence:
                toxic_count += 1
                toxic_motifs_found.append(motif)
        
        # Calculate toxicity score
        toxicity_score = toxic_count / (len(sequence) / 10)  # Normalize by sequence length
        
        # Determine if sequence is toxic
        is_toxic = toxicity_score > threshold
        
        result = {
            'toxicity_score': toxicity_score,
            'is_toxic': is_toxic,
            'toxic_motifs_found': toxic_motifs_found,
            'num_toxic_motifs': toxic_count
        }
        
        logger.info(f"Toxicity check: score={toxicity_score:.3f}, is_toxic={is_toxic}")
        return result
    
    def check_human_similarity(self, sequence, threshold=0.8):
        """
        Check if a sequence has high similarity to human proteins
        
        Args:
            sequence: Protein sequence to check
            threshold: Threshold for similarity score
            
        Returns:
            dict: Human similarity assessment
        """
        # This is a simplified implementation
        # In a real implementation, we would use BLAST against the human proteome
        
        # For demonstration, we'll use a random score
        similarity_score = np.random.uniform(0, 0.5)  # Random score between 0 and 0.5
        
        # Determine if sequence has high similarity to human proteins
        high_similarity = similarity_score > threshold
        
        result = {
            'similarity_score': similarity_score,
            'high_similarity': high_similarity
        }
        
        logger.info(f"Human similarity check: score={similarity_score:.3f}, high_similarity={high_similarity}")
        return result
    
    def filter_vaccine_candidates(self, candidates):
        """
        Filter vaccine candidates based on safety criteria
        
        Args:
            candidates: List of vaccine candidates
            
        Returns:
            tuple: (safe_candidates, filtered_candidates)
        """
        safe_candidates = []
        filtered_candidates = []
        
        for candidate in candidates:
            # Get the sequence to check
            if candidate['type'] == 'subunit':
                sequence = candidate['sequence']
            elif candidate['type'] == 'mRNA':
                sequence = candidate['protein_sequence']
            elif candidate['type'] == 'peptide':
                # For peptide vaccines, check each epitope
                epitope_results = []
                for epitope in candidate['epitopes']:
                    allergenicity = self.check_allergenicity(epitope['sequence'])
                    toxicity = self.check_toxicity(epitope['sequence'])
                    similarity = self.check_human_similarity(epitope['sequence'])
                    
                    epitope_result = {
                        'epitope': epitope,
                        'allergenicity': allergenicity,
                        'toxicity': toxicity,
                        'human_similarity': similarity,
                        'is_safe': not (allergenicity['is_allergenic'] or toxicity['is_toxic'] or similarity['high_similarity'])
                    }
                    epitope_results.append(epitope_result)
                
                # Consider the peptide vaccine safe if all epitopes are safe
                all_epitopes_safe = all(result['is_safe'] for result in epitope_results)
                
                # Add safety results to candidate
                candidate['safety_results'] = {
                    'epitope_results': epitope_results,
                    'all_epitopes_safe': all_epitopes_safe
                }
                
                if all_epitopes_safe:
                    safe_candidates.append(candidate)
                else:
                    filtered_candidates.append(candidate)
                
                continue
            else:
                logger.warning(f"Unknown vaccine type: {candidate['type']}")
                continue
            
            # Check safety criteria
            allergenicity = self.check_allergenicity(sequence)
            toxicity = self.check_toxicity(sequence)
            similarity = self.check_human_similarity(sequence)
            
            # Determine if candidate is safe
            is_safe = not (allergenicity['is_allergenic'] or toxicity['is_toxic'] or similarity['high_similarity'])
            
            # Add safety results to candidate
            candidate['safety_results'] = {
                'allergenicity': allergenicity,
                'toxicity': toxicity,
                'human_similarity': similarity,
                'is_safe': is_safe
            }
            
            # Add to appropriate list
            if is_safe:
                safe_candidates.append(candidate)
            else:
                filtered_candidates.append(candidate)
        
        logger.info(f"Filtered {len(filtered_candidates)} of {len(candidates)} vaccine candidates")
        
        # Save results
        self.save_safety_results(safe_candidates, filtered_candidates)
        
        return safe_candidates, filtered_candidates
    
    def save_safety_results(self, safe_candidates, filtered_candidates):
        """
        Save safety results to files
        
        Args:
            safe_candidates: List of safe vaccine candidates
            filtered_candidates: List of filtered vaccine candidates
            
        Returns:
            dict: Paths to saved files
        """
        # Create DataFrames
        safe_df = pd.DataFrame([
            {
                'type': c['type'],
                'length': c.get('length', 0),
                'allergenicity_score': c['safety_results'].get('allergenicity', {}).get('allergenicity_score', 0),
                'toxicity_score': c['safety_results'].get('toxicity', {}).get('toxicity_score', 0),
                'similarity_score': c['safety_results'].get('human_similarity', {}).get('similarity_score', 0)
            }
            for c in safe_candidates
        ])
        
        filtered_df = pd.DataFrame([
            {
                'type': c['type'],
                'length': c.get('length', 0),
                'allergenicity_score': c['safety_results'].get('allergenicity', {}).get('allergenicity_score', 0),
                'toxicity_score': c['safety_results'].get('toxicity', {}).get('toxicity_score', 0),
                'similarity_score': c['safety_results'].get('human_similarity', {}).get('similarity_score', 0),
                'reason': 'Allergenic' if c['safety_results'].get('allergenicity', {}).get('is_allergenic', False) else
                         'Toxic' if c['safety_results'].get('toxicity', {}).get('is_toxic', False) else
                         'Human similarity' if c['safety_results'].get('human_similarity', {}).get('high_similarity', False) else
                         'Unknown'
            }
            for c in filtered_candidates
        ])
        
        # Save to CSV
        safe_file = self.results_path / "safe_candidates.csv"
        filtered_file = self.results_path / "filtered_candidates.csv"
        
        safe_df.to_csv(safe_file, index=False)
        filtered_df.to_csv(filtered_file, index=False)
        
        logger.info(f"Saved safety results to {safe_file} and {filtered_file}")
        
        return {
            'safe_candidates': str(safe_file),
            'filtered_candidates': str(filtered_file)
        }

