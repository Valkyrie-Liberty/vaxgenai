"""
Immune Evasion Modeling Module for VaxGenAI

This module implements advanced immune evasion and suppression modeling to address
challenges in vaccine development for diseases with immune evasion mechanisms.

Key Features:
- Immune evasion mechanism prediction
- Immunosuppressive microenvironment modeling
- Checkpoint inhibitor compatibility assessment
- Conserved region identification for resistant epitopes
- Tumor microenvironment integration

Author: VaxGenAI Development Team
"""

import pandas as pd
import numpy as np
import logging
from pathlib import Path
from typing import List, Dict, Tuple, Optional
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class ImmuneEvasionModeler:
    """
    Advanced immune evasion modeling system for vaccine design.
    
    This class models immune evasion mechanisms and immunosuppressive
    microenvironments to design vaccines that can overcome these challenges.
    """
    
    def __init__(self):
        """Initialize the immune evasion modeler."""
        self.evasion_mechanisms = {
            'mhc_downregulation': {
                'description': 'MHC Class I downregulation',
                'diseases': ['cancer', 'hiv', 'cmv'],
                'resistance_strategies': ['mhc_ii_targeting', 'nk_cell_activation']
            },
            'antigen_loss': {
                'description': 'Antigen loss variants',
                'diseases': ['cancer', 'hiv'],
                'resistance_strategies': ['conserved_epitopes', 'multi_epitope_targeting']
            },
            'immunosuppression': {
                'description': 'Immunosuppressive microenvironment',
                'diseases': ['cancer'],
                'resistance_strategies': ['checkpoint_inhibitors', 'adjuvants']
            },
            'mutation_escape': {
                'description': 'Immune escape mutations',
                'diseases': ['hiv', 'influenza', 'cancer'],
                'resistance_strategies': ['conserved_regions', 'broad_targeting']
            },
            'regulatory_t_cells': {
                'description': 'Regulatory T-cell suppression',
                'diseases': ['cancer'],
                'resistance_strategies': ['treg_depletion', 'th1_polarization']
            },
            'pd1_pdl1_pathway': {
                'description': 'PD-1/PD-L1 checkpoint inhibition',
                'diseases': ['cancer'],
                'resistance_strategies': ['pd1_blockade', 'pdl1_blockade']
            }
        }
        
        self.disease_profiles = {
            'leukemia': {
                'primary_evasion': ['mhc_downregulation', 'immunosuppression', 'regulatory_t_cells'],
                'microenvironment': 'immunosuppressive',
                'checkpoint_pathways': ['pd1_pdl1', 'ctla4'],
                'resistance_factors': ['minimal_residual_disease', 'blast_heterogeneity']
            },
            'pancreatic_cancer': {
                'primary_evasion': ['immunosuppression', 'regulatory_t_cells', 'mhc_downregulation'],
                'microenvironment': 'highly_immunosuppressive',
                'checkpoint_pathways': ['pd1_pdl1', 'ctla4', 'tim3', 'lag3'],
                'resistance_factors': ['stromal_barrier', 'myeloid_suppression', 'hypoxia']
            },
            'hiv': {
                'primary_evasion': ['mutation_escape', 'antigen_loss', 'mhc_downregulation'],
                'microenvironment': 'chronic_inflammation',
                'checkpoint_pathways': ['pd1_pdl1', 'ctla4'],
                'resistance_factors': ['viral_latency', 'high_mutation_rate', 'immune_exhaustion']
            }
        }
        
        # Initialize models
        self._initialize_evasion_models()
    
    def _initialize_evasion_models(self):
        """Initialize machine learning models for immune evasion prediction."""
        logger.info("Initializing immune evasion prediction models")
        
        # Evasion susceptibility classifier
        self.evasion_classifier = RandomForestClassifier(
            n_estimators=100,
            max_depth=10,
            random_state=42
        )
        
        # Train with simulated data
        self._train_evasion_models()
    
    def _train_evasion_models(self):
        """Train models with simulated immune evasion data."""
        logger.info("Training immune evasion models with simulated data")
        
        # Generate simulated training data
        n_samples = 3000
        features = []
        evasion_susceptibility = []
        
        for _ in range(n_samples):
            # Simulate epitope features
            feature_vector = {
                'conservation_score': np.random.beta(2, 2),  # Conservation level
                'mutation_frequency': np.random.exponential(0.1),  # Mutation frequency
                'mhc_binding_strength': np.random.beta(3, 2),  # MHC binding
                'expression_level': np.random.lognormal(1, 0.5),  # Expression
                'tissue_specificity': np.random.beta(2, 3),  # Tissue specificity
                'immunodominance': np.random.beta(2, 2),  # Immunodominance
                'structural_accessibility': np.random.beta(3, 2)  # Structural accessibility
            }
            
            features.append(list(feature_vector.values()))
            
            # Simulate evasion susceptibility based on features
            evasion_risk = (
                0.3 * (1 - feature_vector['conservation_score']) +  # Low conservation = high risk
                0.2 * min(1, feature_vector['mutation_frequency']) +  # High mutation = high risk
                0.2 * (1 - feature_vector['mhc_binding_strength']) +  # Weak binding = high risk
                0.1 * (1 - feature_vector['expression_level'] / 5) +  # Low expression = high risk
                0.1 * feature_vector['immunodominance'] +  # High immunodominance = high risk
                0.1 * np.random.random()  # Random component
            )
            
            evasion_susceptibility.append(1 if evasion_risk > 0.5 else 0)
        
        # Train classifier
        X = np.array(features)
        y = np.array(evasion_susceptibility)
        self.evasion_classifier.fit(X, y)
        
        logger.info("Immune evasion model training completed")
    
    def assess_evasion_susceptibility(self, epitopes: List[Dict], 
                                    disease_type: str = 'cancer') -> List[Dict]:
        """
        Assess immune evasion susceptibility for epitopes.
        
        Args:
            epitopes: List of epitope dictionaries
            disease_type: Type of disease (cancer, hiv, etc.)
            
        Returns:
            Epitopes with evasion susceptibility assessments
        """
        logger.info(f"Assessing immune evasion susceptibility for {len(epitopes)} epitopes")
        
        for epitope in epitopes:
            # Calculate evasion-related features
            evasion_features = self._calculate_evasion_features(epitope, disease_type)
            epitope['evasion_features'] = evasion_features
            
            # Predict evasion susceptibility
            feature_vector = [
                evasion_features['conservation_score'],
                evasion_features['mutation_frequency'],
                evasion_features['mhc_binding_strength'],
                evasion_features['expression_level'],
                evasion_features['tissue_specificity'],
                evasion_features['immunodominance'],
                evasion_features['structural_accessibility']
            ]
            
            evasion_probability = self.evasion_classifier.predict_proba([feature_vector])[0][1]
            epitope['evasion_susceptibility'] = evasion_probability
            
            # Assess specific evasion mechanisms
            epitope['evasion_mechanisms'] = self._assess_specific_mechanisms(epitope, disease_type)
            
            # Calculate resistance score
            epitope['evasion_resistance_score'] = self._calculate_resistance_score(epitope, disease_type)
        
        return epitopes
    
    def _calculate_evasion_features(self, epitope: Dict, disease_type: str) -> Dict:
        """
        Calculate features related to immune evasion.
        
        Args:
            epitope: Epitope dictionary
            disease_type: Type of disease
            
        Returns:
            Dictionary of evasion-related features
        """
        peptide = epitope.get('peptide', '')
        
        # Conservation score (simulate based on peptide properties)
        # In reality, this would be calculated from sequence alignments
        conservation_score = self._estimate_conservation_score(peptide)
        
        # Mutation frequency (simulate based on position and disease)
        mutation_frequency = self._estimate_mutation_frequency(epitope, disease_type)
        
        # MHC binding strength (from existing predictions)
        mhc_binding_strength = self._extract_mhc_binding_strength(epitope)
        
        # Expression level (from existing data or estimates)
        expression_level = epitope.get('expression_tpm', 5.0) / 10.0  # Normalize
        
        # Tissue specificity (simulate)
        tissue_specificity = np.random.beta(2, 3)
        
        # Immunodominance (simulate based on binding and expression)
        immunodominance = min(1, mhc_binding_strength * expression_level)
        
        # Structural accessibility (simulate)
        structural_accessibility = np.random.beta(3, 2)
        
        return {
            'conservation_score': conservation_score,
            'mutation_frequency': mutation_frequency,
            'mhc_binding_strength': mhc_binding_strength,
            'expression_level': min(1, expression_level),
            'tissue_specificity': tissue_specificity,
            'immunodominance': immunodominance,
            'structural_accessibility': structural_accessibility
        }
    
    def _estimate_conservation_score(self, peptide: str) -> float:
        """Estimate conservation score based on peptide properties."""
        # Simulate conservation based on amino acid composition
        # Hydrophobic and aromatic residues tend to be more conserved
        hydrophobic_aa = 'AILMFWYV'
        aromatic_aa = 'FWY'
        
        hydrophobic_fraction = sum(1 for aa in peptide if aa in hydrophobic_aa) / len(peptide)
        aromatic_fraction = sum(1 for aa in peptide if aa in aromatic_aa) / len(peptide)
        
        # Simulate conservation score
        conservation_score = 0.3 + 0.4 * hydrophobic_fraction + 0.3 * aromatic_fraction
        conservation_score += np.random.normal(0, 0.1)  # Add noise
        
        return max(0, min(1, conservation_score))
    
    def _estimate_mutation_frequency(self, epitope: Dict, disease_type: str) -> float:
        """Estimate mutation frequency based on disease type and epitope properties."""
        base_mutation_rates = {
            'hiv': 0.3,  # High mutation rate
            'cancer': 0.1,  # Moderate mutation rate
            'leukemia': 0.08,  # Lower mutation rate
            'pancreatic_cancer': 0.12  # Moderate-high mutation rate
        }
        
        base_rate = base_mutation_rates.get(disease_type, 0.1)
        
        # Adjust based on epitope properties
        peptide = epitope.get('peptide', '')
        
        # CpG dinucleotides are mutation hotspots
        cpg_sites = peptide.count('CG')
        cpg_adjustment = cpg_sites * 0.05
        
        # Repetitive sequences have higher mutation rates
        unique_aa = len(set(peptide))
        repetitive_adjustment = (1 - unique_aa / len(peptide)) * 0.1
        
        mutation_frequency = base_rate + cpg_adjustment + repetitive_adjustment
        return min(1, mutation_frequency)
    
    def _extract_mhc_binding_strength(self, epitope: Dict) -> float:
        """Extract MHC binding strength from epitope data."""
        # Look for existing MHC binding data
        mhc_binding = epitope.get('mhc_binding', {})
        if not mhc_binding:
            mhc_binding = epitope.get('mhc_ii_binding', {})
        
        if mhc_binding:
            # Find strongest binding
            best_binding = 0
            for allele, binding_data in mhc_binding.items():
                if isinstance(binding_data, dict):
                    binding_score = binding_data.get('binding_score', 0)
                    if binding_data.get('strong_binder', False):
                        binding_score = max(binding_score, 0.8)
                    elif binding_data.get('weak_binder', False):
                        binding_score = max(binding_score, 0.5)
                    best_binding = max(best_binding, binding_score)
            return best_binding
        
        # Default moderate binding strength
        return 0.5
    
    def _assess_specific_mechanisms(self, epitope: Dict, disease_type: str) -> Dict:
        """
        Assess susceptibility to specific immune evasion mechanisms.
        
        Args:
            epitope: Epitope dictionary
            disease_type: Type of disease
            
        Returns:
            Dictionary of mechanism-specific assessments
        """
        mechanisms = {}
        
        disease_profile = self.disease_profiles.get(disease_type, {})
        primary_evasion = disease_profile.get('primary_evasion', [])
        
        for mechanism in primary_evasion:
            if mechanism == 'mhc_downregulation':
                # Assess susceptibility to MHC downregulation
                mhc_dependence = epitope['evasion_features']['mhc_binding_strength']
                mechanisms[mechanism] = {
                    'susceptibility': mhc_dependence,
                    'resistance_strategies': ['mhc_ii_targeting', 'nk_cell_activation']
                }
            
            elif mechanism == 'antigen_loss':
                # Assess susceptibility to antigen loss
                expression_dependence = epitope['evasion_features']['expression_level']
                conservation = epitope['evasion_features']['conservation_score']
                susceptibility = expression_dependence * (1 - conservation)
                mechanisms[mechanism] = {
                    'susceptibility': susceptibility,
                    'resistance_strategies': ['conserved_epitopes', 'multi_epitope_targeting']
                }
            
            elif mechanism == 'immunosuppression':
                # Assess susceptibility to immunosuppressive environment
                immunodominance = epitope['evasion_features']['immunodominance']
                mechanisms[mechanism] = {
                    'susceptibility': 1 - immunodominance,  # Low immunodominance = high susceptibility
                    'resistance_strategies': ['checkpoint_inhibitors', 'adjuvants']
                }
            
            elif mechanism == 'mutation_escape':
                # Assess susceptibility to mutation escape
                mutation_freq = epitope['evasion_features']['mutation_frequency']
                conservation = epitope['evasion_features']['conservation_score']
                susceptibility = mutation_freq * (1 - conservation)
                mechanisms[mechanism] = {
                    'susceptibility': susceptibility,
                    'resistance_strategies': ['conserved_regions', 'broad_targeting']
                }
            
            elif mechanism == 'regulatory_t_cells':
                # Assess susceptibility to Treg suppression
                tissue_specificity = epitope['evasion_features']['tissue_specificity']
                mechanisms[mechanism] = {
                    'susceptibility': tissue_specificity,  # Tissue-specific antigens more susceptible
                    'resistance_strategies': ['treg_depletion', 'th1_polarization']
                }
        
        return mechanisms
    
    def _calculate_resistance_score(self, epitope: Dict, disease_type: str) -> float:
        """
        Calculate overall resistance score to immune evasion.
        
        Args:
            epitope: Epitope dictionary
            disease_type: Type of disease
            
        Returns:
            Resistance score (0-1, higher is better)
        """
        # Base resistance from conservation and low mutation frequency
        conservation = epitope['evasion_features']['conservation_score']
        mutation_resistance = 1 - epitope['evasion_features']['mutation_frequency']
        
        # MHC binding stability
        mhc_stability = epitope['evasion_features']['mhc_binding_strength']
        
        # Expression level (moderate levels preferred)
        expression = epitope['evasion_features']['expression_level']
        expression_score = 1 - abs(expression - 0.5) * 2  # Optimal around 0.5
        
        # Structural accessibility
        accessibility = epitope['evasion_features']['structural_accessibility']
        
        # Disease-specific adjustments
        disease_weight = 1.0
        if disease_type in ['pancreatic_cancer']:
            disease_weight = 0.8  # More challenging environment
        elif disease_type in ['hiv']:
            disease_weight = 0.7  # High mutation rate
        
        # Combined resistance score
        resistance_score = (
            0.3 * conservation +
            0.25 * mutation_resistance +
            0.2 * mhc_stability +
            0.15 * expression_score +
            0.1 * accessibility
        ) * disease_weight
        
        return max(0, min(1, resistance_score))
    
    def design_evasion_resistant_vaccines(self, epitopes: List[Dict], 
                                        disease_type: str) -> Dict:
        """
        Design vaccine strategies resistant to immune evasion.
        
        Args:
            epitopes: List of epitopes with evasion assessments
            disease_type: Type of disease
            
        Returns:
            Dictionary with vaccine design recommendations
        """
        logger.info(f"Designing evasion-resistant vaccine strategies for {disease_type}")
        
        # Filter epitopes by resistance score
        resistant_epitopes = [e for e in epitopes if e['evasion_resistance_score'] > 0.6]
        
        # Group epitopes by resistance mechanisms
        mechanism_groups = {}
        for epitope in resistant_epitopes:
            for mechanism, data in epitope['evasion_mechanisms'].items():
                if mechanism not in mechanism_groups:
                    mechanism_groups[mechanism] = []
                if data['susceptibility'] < 0.5:  # Low susceptibility = good resistance
                    mechanism_groups[mechanism].append(epitope)
        
        # Design multi-layered vaccine strategy
        vaccine_strategy = {
            'primary_epitopes': resistant_epitopes[:10],  # Top 10 resistant epitopes
            'backup_epitopes': resistant_epitopes[10:20],  # Backup epitopes
            'resistance_mechanisms': mechanism_groups,
            'adjuvant_recommendations': self._recommend_adjuvants(disease_type),
            'combination_therapy': self._recommend_combination_therapy(disease_type),
            'delivery_optimization': self._optimize_delivery_for_evasion(disease_type)
        }
        
        return vaccine_strategy
    
    def _recommend_adjuvants(self, disease_type: str) -> List[Dict]:
        """Recommend adjuvants based on disease-specific evasion mechanisms."""
        adjuvant_recommendations = []
        
        disease_profile = self.disease_profiles.get(disease_type, {})
        
        if 'immunosuppression' in disease_profile.get('primary_evasion', []):
            adjuvant_recommendations.append({
                'type': 'TLR_agonist',
                'examples': ['CpG_ODN', 'Poly_IC', 'LPS'],
                'mechanism': 'Overcome immunosuppression via innate immune activation'
            })
        
        if 'regulatory_t_cells' in disease_profile.get('primary_evasion', []):
            adjuvant_recommendations.append({
                'type': 'Th1_polarizing',
                'examples': ['IL-12', 'IFN-gamma', 'CFA'],
                'mechanism': 'Promote Th1 responses and reduce Treg activity'
            })
        
        if disease_type in ['cancer', 'leukemia', 'pancreatic_cancer']:
            adjuvant_recommendations.append({
                'type': 'checkpoint_inhibitor',
                'examples': ['anti-PD1', 'anti-PDL1', 'anti-CTLA4'],
                'mechanism': 'Block inhibitory checkpoint pathways'
            })
        
        return adjuvant_recommendations
    
    def _recommend_combination_therapy(self, disease_type: str) -> List[Dict]:
        """Recommend combination therapies to overcome immune evasion."""
        combinations = []
        
        if disease_type in ['cancer', 'leukemia', 'pancreatic_cancer']:
            combinations.append({
                'strategy': 'vaccine_plus_checkpoint_inhibition',
                'components': ['vaccine', 'anti-PD1', 'anti-CTLA4'],
                'rationale': 'Enhance T-cell activation and overcome exhaustion'
            })
            
            combinations.append({
                'strategy': 'vaccine_plus_adoptive_transfer',
                'components': ['vaccine', 'CAR-T', 'TIL'],
                'rationale': 'Combine active and passive immunotherapy'
            })
        
        if disease_type == 'hiv':
            combinations.append({
                'strategy': 'vaccine_plus_latency_reversal',
                'components': ['vaccine', 'HDAC_inhibitors', 'PKC_agonists'],
                'rationale': 'Activate latent virus for immune targeting'
            })
        
        return combinations
    
    def _optimize_delivery_for_evasion(self, disease_type: str) -> Dict:
        """Optimize vaccine delivery to overcome evasion mechanisms."""
        delivery_optimization = {
            'route': 'intramuscular',  # Default
            'formulation': 'standard',
            'timing': 'prime_boost'
        }
        
        disease_profile = self.disease_profiles.get(disease_type, {})
        
        if disease_profile.get('microenvironment') == 'highly_immunosuppressive':
            delivery_optimization.update({
                'route': 'intratumoral_plus_systemic',
                'formulation': 'nanoparticle_encapsulated',
                'timing': 'multiple_boost_extended'
            })
        
        if 'mhc_downregulation' in disease_profile.get('primary_evasion', []):
            delivery_optimization.update({
                'formulation': 'mhc_ii_targeting',
                'adjuvant': 'TLR_agonist'
            })
        
        return delivery_optimization
    
    def generate_evasion_resistance_report(self, epitopes: List[Dict], 
                                         vaccine_strategy: Dict,
                                         disease_type: str,
                                         output_path: str) -> str:
        """
        Generate a comprehensive immune evasion resistance report.
        
        Args:
            epitopes: List of epitopes with evasion assessments
            vaccine_strategy: Vaccine strategy recommendations
            disease_type: Type of disease
            output_path: Output file path
            
        Returns:
            Path to generated report
        """
        logger.info("Generating immune evasion resistance report")
        
        # Create DataFrame for easier manipulation
        df = pd.DataFrame(epitopes)
        
        # Generate report
        report = f"# Immune Evasion Resistance Report - {disease_type.title()}\n\n"
        
        # Summary statistics
        report += "## Summary\n\n"
        report += f"Total epitopes analyzed: {len(epitopes)}\n"
        report += f"Evasion-resistant epitopes (score > 0.6): {len(df[df['evasion_resistance_score'] > 0.6])}\n"
        report += f"High-risk epitopes (susceptibility > 0.7): {len(df[df['evasion_susceptibility'] > 0.7])}\n"
        report += f"Average resistance score: {df['evasion_resistance_score'].mean():.3f}\n"
        report += f"Average evasion susceptibility: {df['evasion_susceptibility'].mean():.3f}\n\n"
        
        # Disease-specific evasion profile
        disease_profile = self.disease_profiles.get(disease_type, {})
        report += f"## {disease_type.title()} Immune Evasion Profile\n\n"
        report += f"**Primary Evasion Mechanisms:**\n"
        for mechanism in disease_profile.get('primary_evasion', []):
            description = self.evasion_mechanisms.get(mechanism, {}).get('description', mechanism)
            report += f"- {description}\n"
        
        report += f"\n**Microenvironment:** {disease_profile.get('microenvironment', 'Unknown')}\n"
        report += f"**Checkpoint Pathways:** {', '.join(disease_profile.get('checkpoint_pathways', []))}\n\n"
        
        # Top resistant epitopes
        resistant_epitopes = df[df['evasion_resistance_score'] > 0.6].sort_values('evasion_resistance_score', ascending=False)
        
        report += "## Top 15 Evasion-Resistant Epitopes\n\n"
        report += "| Rank | Peptide | Resistance Score | Susceptibility | Conservation | Mutation Freq |\n"
        report += "|------|---------|------------------|----------------|--------------|---------------|\n"
        
        for i, (_, epitope) in enumerate(resistant_epitopes.head(15).iterrows()):
            conservation = epitope['evasion_features']['conservation_score']
            mutation_freq = epitope['evasion_features']['mutation_frequency']
            
            report += f"| {i+1} | {epitope['peptide']} | "
            report += f"{epitope['evasion_resistance_score']:.3f} | "
            report += f"{epitope['evasion_susceptibility']:.3f} | "
            report += f"{conservation:.3f} | {mutation_freq:.3f} |\n"
        
        # Mechanism-specific analysis
        report += "\n## Evasion Mechanism Analysis\n\n"
        
        for mechanism in disease_profile.get('primary_evasion', []):
            mechanism_data = self.evasion_mechanisms.get(mechanism, {})
            report += f"### {mechanism_data.get('description', mechanism)}\n\n"
            
            # Count epitopes by susceptibility to this mechanism
            susceptible_count = 0
            resistant_count = 0
            
            for epitope in epitopes:
                mech_data = epitope.get('evasion_mechanisms', {}).get(mechanism, {})
                if mech_data:
                    if mech_data['susceptibility'] > 0.5:
                        susceptible_count += 1
                    else:
                        resistant_count += 1
            
            report += f"- Susceptible epitopes: {susceptible_count}\n"
            report += f"- Resistant epitopes: {resistant_count}\n"
            
            # Resistance strategies
            strategies = mechanism_data.get('resistance_strategies', [])
            if strategies:
                report += f"- Recommended resistance strategies: {', '.join(strategies)}\n"
            
            report += "\n"
        
        # Vaccine strategy recommendations
        report += "## Vaccine Strategy Recommendations\n\n"
        
        # Primary epitopes
        primary_epitopes = vaccine_strategy.get('primary_epitopes', [])
        report += f"### Primary Epitopes ({len(primary_epitopes)} selected)\n\n"
        for i, epitope in enumerate(primary_epitopes):
            report += f"{i+1}. {epitope['peptide']} (Resistance: {epitope['evasion_resistance_score']:.3f})\n"
        
        # Adjuvant recommendations
        adjuvants = vaccine_strategy.get('adjuvant_recommendations', [])
        if adjuvants:
            report += "\n### Recommended Adjuvants\n\n"
            for adj in adjuvants:
                report += f"**{adj['type']}**\n"
                report += f"- Examples: {', '.join(adj['examples'])}\n"
                report += f"- Mechanism: {adj['mechanism']}\n\n"
        
        # Combination therapy
        combinations = vaccine_strategy.get('combination_therapy', [])
        if combinations:
            report += "### Combination Therapy Recommendations\n\n"
            for combo in combinations:
                report += f"**{combo['strategy']}**\n"
                report += f"- Components: {', '.join(combo['components'])}\n"
                report += f"- Rationale: {combo['rationale']}\n\n"
        
        # Delivery optimization
        delivery = vaccine_strategy.get('delivery_optimization', {})
        if delivery:
            report += "### Delivery Optimization\n\n"
            report += f"- Route: {delivery.get('route', 'Standard')}\n"
            report += f"- Formulation: {delivery.get('formulation', 'Standard')}\n"
            report += f"- Timing: {delivery.get('timing', 'Standard')}\n\n"
        
        # Recommendations
        report += "## Key Recommendations\n\n"
        
        resistant_count = len(df[df['evasion_resistance_score'] > 0.6])
        if resistant_count >= 10:
            report += "Excellent pool of evasion-resistant epitopes identified. Proceed with multi-epitope vaccine design.\n\n"
        elif resistant_count >= 5:
            report += "Good selection of evasion-resistant epitopes. Consider combining with additional strategies.\n\n"
        else:
            report += "Limited evasion-resistant epitopes. Strongly recommend combination with checkpoint inhibitors and adjuvants.\n\n"
        
        report += "Implementation priorities:\n"
        report += "1. Focus on epitopes with resistance score > 0.6\n"
        report += "2. Implement recommended combination therapies\n"
        report += "3. Use appropriate adjuvants for disease-specific evasion mechanisms\n"
        report += "4. Consider optimized delivery strategies\n"
        report += "5. Plan for monitoring of immune escape variants\n"
        report += "6. Validate resistance predictions in relevant disease models\n"
        
        # Save report
        with open(output_path, 'w') as f:
            f.write(report)
        
        logger.info(f"Immune evasion resistance report saved to: {output_path}")
        return output_path

# Testing and example usage
def test_immune_evasion_modeling():
    """Test the immune evasion modeling system."""
    logger.info("Testing immune evasion modeling system")
    
    # Sample epitopes with basic data
    test_epitopes = [
        {
            'peptide': 'YLQPRTFLL',
            'mhc_binding': {'HLA-A*02:01': {'binding_score': 0.8, 'strong_binder': True}},
            'expression_tpm': 8.5
        },
        {
            'peptide': 'KIADYNYKL',
            'mhc_binding': {'HLA-A*02:01': {'binding_score': 0.7, 'weak_binder': True}},
            'expression_tpm': 12.3
        },
        {
            'peptide': 'SLYNTVATL',
            'mhc_binding': {'HLA-A*02:01': {'binding_score': 0.9, 'strong_binder': True}},
            'expression_tpm': 15.2
        },
        {
            'peptide': 'GILGFVFTL',
            'mhc_binding': {'HLA-A*02:01': {'binding_score': 0.6, 'weak_binder': True}},
            'expression_tpm': 6.8
        }
    ]
    
    # Test different disease types
    disease_types = ['leukemia', 'pancreatic_cancer', 'hiv']
    
    results = {}
    
    for disease_type in disease_types:
        logger.info(f"Testing immune evasion modeling for {disease_type}")
        
        # Initialize modeler
        modeler = ImmuneEvasionModeler()
        
        # Assess evasion susceptibility
        assessed_epitopes = modeler.assess_evasion_susceptibility(test_epitopes.copy(), disease_type)
        
        # Design evasion-resistant vaccine
        vaccine_strategy = modeler.design_evasion_resistant_vaccines(assessed_epitopes, disease_type)
        
        # Generate report
        output_dir = Path(f"/home/ubuntu/vaxgenai/results/immune_evasion/{disease_type}")
        output_dir.mkdir(parents=True, exist_ok=True)
        
        report_path = modeler.generate_evasion_resistance_report(
            assessed_epitopes,
            vaccine_strategy,
            disease_type,
            str(output_dir / "evasion_resistance_report.md")
        )
        
        # Save data
        df = pd.DataFrame(assessed_epitopes)
        df.to_csv(output_dir / "evasion_analysis.csv", index=False)
        
        results[disease_type] = {
            'total_epitopes': len(assessed_epitopes),
            'resistant_epitopes': len([e for e in assessed_epitopes if e['evasion_resistance_score'] > 0.6]),
            'average_resistance': np.mean([e['evasion_resistance_score'] for e in assessed_epitopes]),
            'report_path': report_path
        }
    
    logger.info("Immune evasion modeling test completed")
    return results

if __name__ == "__main__":
    # Run test
    test_results = test_immune_evasion_modeling()
    for disease, results in test_results.items():
        print(f"\n{disease.title()} Results:")
        print(f"  Resistant epitopes: {results['resistant_epitopes']}/{results['total_epitopes']}")
        print(f"  Average resistance score: {results['average_resistance']:.3f}")

