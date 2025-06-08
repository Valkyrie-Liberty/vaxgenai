"""
Enhanced safety filter for VaxGenAI

This module provides an enhanced safety filter for vaccine candidates
that checks for allergenicity, toxicity, and human similarity.
"""

import logging
import numpy as np
from typing import List, Dict, Tuple, Optional

# Setup logging
logger = logging.getLogger("vaxgenai.enhanced.safety_filter")

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
            # Get the sequence to check
            sequence = self._get_sequence_from_candidate(candidate)
            
            if not sequence:
                logger.warning(f"Could not find sequence in candidate: {candidate.keys()}")
                filtered_candidates.append(candidate)
                continue
            
            # Check safety
            allergenicity = self._check_allergenicity(sequence)
            toxicity = self._check_toxicity(sequence)
            human_similarity = self._check_human_similarity(sequence)
            
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
    
    def _get_sequence_from_candidate(self, candidate: Dict) -> str:
        """
        Get the sequence from a vaccine candidate, handling different field names.
        
        Args:
            candidate: Vaccine candidate dictionary
            
        Returns:
            Sequence string or empty string if not found
        """
        # Try different field names that might contain the sequence
        sequence_fields = ['sequence', 'vaccine_sequence', 'protein_sequence']
        
        for field in sequence_fields:
            if field in candidate and candidate[field]:
                return candidate[field]
        
        # If no sequence field is found, return empty string
        return ""
    
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

