"""
Vaccine Delivery System Integration Module for VaxGenAI

This module integrates epitope predictions with various vaccine delivery platforms
to optimize vaccine design for manufacturing and clinical efficacy.

Key Features:
- mRNA vaccine design and optimization
- Viral vector platform integration
- Peptide vaccine formulation
- Adjuvant selection and optimization
- Manufacturing compatibility assessment
- Platform-specific design recommendations

Author: VaxGenAI Development Team
"""

import pandas as pd
import numpy as np
import logging
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import json
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class VaccineDeliverySystemIntegrator:
    """
    Vaccine delivery system integration for optimized vaccine design.
    
    This class integrates epitope predictions with various delivery platforms
    to create manufacturable and effective vaccines.
    """
    
    def __init__(self):
        """Initialize the vaccine delivery system integrator."""
        
        # Platform-specific parameters
        self.platform_parameters = {
            'mRNA': {
                'max_sequence_length': 4000,  # nucleotides
                'optimal_length': 1500,
                'gc_content_range': (40, 60),  # percentage
                'codon_optimization': True,
                'stability_requirements': {
                    'min_half_life': 24,  # hours
                    'storage_temp': -70,  # Celsius
                    'formulation': 'lipid_nanoparticle'
                },
                'manufacturing_complexity': 'medium',
                'cost_per_dose': 15.0,  # USD
                'scalability': 'high'
            },
            'viral_vector': {
                'max_insert_size': 8000,  # base pairs
                'optimal_insert_size': 3000,
                'vector_types': ['adenovirus', 'lentivirus', 'AAV'],
                'immunogenicity': 'high',
                'manufacturing_complexity': 'high',
                'cost_per_dose': 45.0,
                'scalability': 'medium',
                'safety_considerations': ['vector_immunity', 'integration_risk']
            },
            'peptide': {
                'max_peptides': 20,
                'optimal_peptides': 8,
                'peptide_length_range': (8, 25),
                'solubility_requirements': True,
                'adjuvant_required': True,
                'manufacturing_complexity': 'low',
                'cost_per_dose': 25.0,
                'scalability': 'high',
                'stability': 'high'
            },
            'protein_subunit': {
                'max_protein_size': 500,  # amino acids
                'optimal_size': 200,
                'expression_systems': ['E.coli', 'yeast', 'mammalian'],
                'purification_required': True,
                'adjuvant_required': True,
                'manufacturing_complexity': 'medium',
                'cost_per_dose': 20.0,
                'scalability': 'high'
            },
            'DNA': {
                'max_sequence_length': 6000,  # base pairs
                'optimal_length': 2000,
                'delivery_methods': ['electroporation', 'gene_gun', 'lipofection'],
                'immunogenicity': 'moderate',
                'manufacturing_complexity': 'low',
                'cost_per_dose': 8.0,
                'scalability': 'very_high',
                'regulatory_complexity': 'high'
            }
        }
        
        # Adjuvant database
        self.adjuvants = {
            'alum': {
                'type': 'mineral_salt',
                'mechanism': 'depot_formation',
                'immune_response': 'Th2',
                'compatibility': ['peptide', 'protein_subunit'],
                'safety_profile': 'excellent',
                'cost': 'low'
            },
            'MF59': {
                'type': 'oil_in_water_emulsion',
                'mechanism': 'immune_activation',
                'immune_response': 'Th1/Th2_balanced',
                'compatibility': ['peptide', 'protein_subunit'],
                'safety_profile': 'good',
                'cost': 'medium'
            },
            'AS01': {
                'type': 'liposome_based',
                'mechanism': 'TLR4_activation',
                'immune_response': 'Th1_dominant',
                'compatibility': ['peptide', 'protein_subunit'],
                'safety_profile': 'good',
                'cost': 'high'
            },
            'CpG_ODN': {
                'type': 'TLR9_agonist',
                'mechanism': 'innate_immune_activation',
                'immune_response': 'Th1_dominant',
                'compatibility': ['peptide', 'protein_subunit', 'DNA'],
                'safety_profile': 'good',
                'cost': 'medium'
            },
            'poly_IC': {
                'type': 'TLR3_agonist',
                'mechanism': 'type_I_interferon',
                'immune_response': 'Th1_dominant',
                'compatibility': ['peptide', 'protein_subunit'],
                'safety_profile': 'moderate',
                'cost': 'medium'
            }
        }
        
        # Manufacturing considerations
        self.manufacturing_factors = {
            'scale_requirements': {
                'clinical_trial': 1000,      # doses
                'regional_deployment': 100000,
                'global_deployment': 1000000000
            },
            'quality_standards': ['GMP', 'FDA', 'EMA', 'WHO'],
            'cold_chain_requirements': {
                'mRNA': -70,
                'viral_vector': -20,
                'peptide': 4,
                'protein_subunit': 4,
                'DNA': 4
            }
        }
    
    def design_platform_specific_vaccines(self, epitopes: List[Dict],
                                        target_platforms: List[str] = None,
                                        disease_context: str = 'cancer') -> Dict:
        """
        Design vaccines for specific delivery platforms.
        
        Args:
            epitopes: List of predicted epitopes
            target_platforms: Platforms to design for
            disease_context: Disease context for optimization
            
        Returns:
            Platform-specific vaccine designs
        """
        logger.info(f"Designing vaccines for {len(epitopes)} epitopes across platforms")
        
        if target_platforms is None:
            target_platforms = ['mRNA', 'viral_vector', 'peptide', 'protein_subunit']
        
        vaccine_designs = {}
        
        for platform in target_platforms:
            if platform in self.platform_parameters:
                design = self._design_for_platform(epitopes, platform, disease_context)
                vaccine_designs[platform] = design
        
        # Compare platforms and provide recommendations
        platform_comparison = self._compare_platforms(vaccine_designs, disease_context)
        
        result = {
            'vaccine_designs': vaccine_designs,
            'platform_comparison': platform_comparison,
            'recommendations': self._generate_platform_recommendations(
                vaccine_designs, platform_comparison, disease_context
            ),
            'manufacturing_assessment': self._assess_manufacturing_feasibility(vaccine_designs)
        }
        
        return result
    
    def _design_for_platform(self, epitopes: List[Dict], platform: str, 
                           disease_context: str) -> Dict:
        """Design vaccine for specific platform."""
        
        platform_params = self.platform_parameters[platform]
        
        if platform == 'mRNA':
            return self._design_mrna_vaccine(epitopes, platform_params, disease_context)
        elif platform == 'viral_vector':
            return self._design_viral_vector_vaccine(epitopes, platform_params, disease_context)
        elif platform == 'peptide':
            return self._design_peptide_vaccine(epitopes, platform_params, disease_context)
        elif platform == 'protein_subunit':
            return self._design_protein_subunit_vaccine(epitopes, platform_params, disease_context)
        elif platform == 'DNA':
            return self._design_dna_vaccine(epitopes, platform_params, disease_context)
        else:
            return {'error': f'Unknown platform: {platform}'}
    
    def _design_mrna_vaccine(self, epitopes: List[Dict], params: Dict, 
                           disease_context: str) -> Dict:
        """Design mRNA vaccine."""
        
        # Select epitopes for mRNA design
        selected_epitopes = self._select_epitopes_for_mrna(epitopes, params)
        
        # Design mRNA sequence
        mrna_sequence = self._construct_mrna_sequence(selected_epitopes, disease_context)
        
        # Optimize codon usage
        optimized_sequence = self._optimize_codons(mrna_sequence)
        
        # Calculate stability metrics
        stability_metrics = self._calculate_mrna_stability(optimized_sequence)
        
        # Design lipid nanoparticle formulation
        lnp_formulation = self._design_lnp_formulation(optimized_sequence)
        
        # Calculate manufacturing parameters
        manufacturing_params = self._calculate_mrna_manufacturing(optimized_sequence, params)
        
        design = {
            'platform': 'mRNA',
            'selected_epitopes': selected_epitopes,
            'mrna_sequence': {
                'length': len(optimized_sequence),
                'gc_content': self._calculate_gc_content(optimized_sequence),
                'codon_optimization_score': 0.85,  # Simulated
                'sequence': optimized_sequence[:100] + '...'  # Truncated for display
            },
            'stability_metrics': stability_metrics,
            'formulation': lnp_formulation,
            'manufacturing': manufacturing_params,
            'predicted_efficacy': self._predict_mrna_efficacy(selected_epitopes, stability_metrics),
            'safety_profile': self._assess_mrna_safety(optimized_sequence, lnp_formulation),
            'regulatory_considerations': self._assess_mrna_regulatory(optimized_sequence)
        }
        
        return design
    
    def _select_epitopes_for_mrna(self, epitopes: List[Dict], params: Dict) -> List[Dict]:
        """Select optimal epitopes for mRNA vaccine design."""
        
        # Sort epitopes by prediction score
        sorted_epitopes = sorted(epitopes, key=lambda x: x.get('prediction_score', 0), reverse=True)
        
        selected = []
        total_length = 0
        max_aa_length = params['optimal_length'] // 3  # Convert to amino acids
        
        for epitope in sorted_epitopes:
            epitope_length = len(epitope.get('sequence', epitope.get('peptide', '')))
            
            if total_length + epitope_length <= max_aa_length:
                selected.append(epitope)
                total_length += epitope_length
                
                # Add linker length
                if len(selected) > 1:
                    total_length += 4  # GGSG linker
        
        return selected
    
    def _construct_mrna_sequence(self, epitopes: List[Dict], disease_context: str) -> str:
        """Construct mRNA sequence from epitopes."""
        
        # Start with signal peptide for secretion
        signal_peptide = "ATGAAGTGGGTGACCTTCATCAGCCTGCTGTTCCTGTTTCTGAGCGCCGCCAGAAGC"
        
        # Add epitopes with linkers
        epitope_sequences = []
        for i, epitope in enumerate(epitopes):
            aa_seq = epitope.get('sequence', epitope.get('peptide', ''))
            # Convert amino acid to DNA (simplified - use optimal codons)
            dna_seq = self._aa_to_optimal_dna(aa_seq)
            epitope_sequences.append(dna_seq)
        
        # Join with GGSG linker
        linker_dna = "GGCGGCAGCGGC"  # GGSG
        epitope_region = linker_dna.join(epitope_sequences)
        
        # Add stop codon
        stop_codon = "TGA"
        
        # Construct full sequence
        full_sequence = signal_peptide + epitope_region + stop_codon
        
        return full_sequence
    
    def _aa_to_optimal_dna(self, aa_sequence: str) -> str:
        """Convert amino acid sequence to optimal DNA codons."""
        
        # Optimal codon table (simplified)
        optimal_codons = {
            'A': 'GCC', 'R': 'CGC', 'N': 'AAC', 'D': 'GAC', 'C': 'TGC',
            'Q': 'CAG', 'E': 'GAG', 'G': 'GGC', 'H': 'CAC', 'I': 'ATC',
            'L': 'CTG', 'K': 'AAG', 'M': 'ATG', 'F': 'TTC', 'P': 'CCC',
            'S': 'AGC', 'T': 'ACC', 'W': 'TGG', 'Y': 'TAC', 'V': 'GTC'
        }
        
        dna_sequence = ""
        for aa in aa_sequence.upper():
            if aa in optimal_codons:
                dna_sequence += optimal_codons[aa]
            else:
                dna_sequence += "NNN"  # Unknown amino acid
        
        return dna_sequence
    
    def _optimize_codons(self, dna_sequence: str) -> str:
        """Optimize codon usage for expression."""
        # In real implementation, use sophisticated codon optimization
        # For now, return the input sequence
        return dna_sequence
    
    def _calculate_gc_content(self, dna_sequence: str) -> float:
        """Calculate GC content of DNA sequence."""
        gc_count = dna_sequence.count('G') + dna_sequence.count('C')
        return (gc_count / len(dna_sequence)) * 100 if dna_sequence else 0
    
    def _calculate_mrna_stability(self, mrna_sequence: str) -> Dict:
        """Calculate mRNA stability metrics."""
        
        gc_content = self._calculate_gc_content(mrna_sequence)
        
        # Simplified stability calculation
        base_stability = 0.5
        
        # GC content contribution
        optimal_gc = 50
        gc_penalty = abs(gc_content - optimal_gc) / 100
        gc_contribution = max(0, 0.3 - gc_penalty)
        
        # Length contribution
        length_score = min(1.0, len(mrna_sequence) / 1500) * 0.2
        
        overall_stability = base_stability + gc_contribution + length_score
        
        return {
            'overall_stability_score': min(1.0, overall_stability),
            'gc_content': gc_content,
            'predicted_half_life_hours': 12 + (overall_stability * 24),
            'degradation_resistance': 'moderate' if overall_stability > 0.7 else 'low',
            'storage_stability': 'good' if overall_stability > 0.6 else 'moderate'
        }
    
    def _design_lnp_formulation(self, mrna_sequence: str) -> Dict:
        """Design lipid nanoparticle formulation."""
        
        # Calculate mRNA mass (simplified)
        mrna_mass_da = len(mrna_sequence) * 330  # Average nucleotide mass
        
        # LNP composition (molar ratios)
        lnp_composition = {
            'ionizable_lipid': 'ALC-0315',  # Pfizer-like
            'helper_lipid': 'DSPC',
            'cholesterol': 'cholesterol',
            'peg_lipid': 'ALC-0159',
            'molar_ratios': [50, 10, 38.5, 1.5],  # Typical ratios
            'n_p_ratio': 6.0  # Nitrogen to phosphate ratio
        }
        
        # Formulation parameters
        formulation_params = {
            'particle_size_nm': 80,
            'pdi': 0.1,  # Polydispersity index
            'zeta_potential_mv': -5,
            'encapsulation_efficiency': 0.85,
            'mrna_concentration_mg_ml': 1.0,
            'buffer': 'PBS_pH_7.4',
            'osmolality_mosm_kg': 300
        }
        
        return {
            'composition': lnp_composition,
            'parameters': formulation_params,
            'stability_profile': {
                'storage_temp_c': -70,
                'shelf_life_months': 6,
                'freeze_thaw_cycles': 1
            },
            'manufacturing_complexity': 'medium',
            'cost_impact': 'moderate'
        }
    
    def _calculate_mrna_manufacturing(self, mrna_sequence: str, params: Dict) -> Dict:
        """Calculate mRNA manufacturing parameters."""
        
        sequence_length = len(mrna_sequence)
        
        # Manufacturing metrics
        manufacturing = {
            'synthesis_method': 'in_vitro_transcription',
            'template_dna_required': True,
            'purification_steps': ['DNase_treatment', 'column_chromatography', 'ultrafiltration'],
            'yield_percentage': 65,
            'purity_percentage': 95,
            'endotoxin_level_eu_ml': 0.1,
            'production_time_days': 7,
            'batch_size_doses': 10000,
            'cost_per_dose_usd': params['cost_per_dose'],
            'scalability': params['scalability'],
            'quality_control_tests': [
                'identity_sequencing',
                'purity_hplc',
                'integrity_gel_electrophoresis',
                'endotoxin_lal',
                'sterility'
            ]
        }
        
        return manufacturing
    
    def _predict_mrna_efficacy(self, epitopes: List[Dict], stability: Dict) -> Dict:
        """Predict mRNA vaccine efficacy."""
        
        # Base efficacy from epitope quality
        epitope_scores = [e.get('prediction_score', 0.5) for e in epitopes]
        base_efficacy = np.mean(epitope_scores) if epitope_scores else 0.5
        
        # Stability contribution
        stability_bonus = stability['overall_stability_score'] * 0.2
        
        # Platform-specific factors
        mrna_bonus = 0.15  # mRNA platforms typically show good efficacy
        
        overall_efficacy = min(1.0, base_efficacy + stability_bonus + mrna_bonus)
        
        return {
            'predicted_efficacy': overall_efficacy,
            'confidence_interval': [overall_efficacy - 0.1, overall_efficacy + 0.1],
            'factors': {
                'epitope_quality': base_efficacy,
                'stability_contribution': stability_bonus,
                'platform_bonus': mrna_bonus
            },
            'expected_immune_response': {
                'antibody_response': 'strong',
                'cd8_response': 'strong',
                'cd4_response': 'moderate',
                'memory_formation': 'good'
            }
        }
    
    def _assess_mrna_safety(self, mrna_sequence: str, lnp_formulation: Dict) -> Dict:
        """Assess mRNA vaccine safety profile."""
        
        safety_assessment = {
            'overall_safety_score': 0.85,  # mRNA generally safe
            'known_risks': [
                'injection_site_reactions',
                'mild_systemic_reactions',
                'rare_allergic_reactions'
            ],
            'contraindications': [
                'severe_immunocompromise',
                'known_peg_allergy'
            ],
            'special_populations': {
                'pregnancy': 'use_with_caution',
                'pediatric': 'age_dependent',
                'elderly': 'generally_safe',
                'immunocompromised': 'reduced_efficacy_expected'
            },
            'long_term_safety': 'good_profile_expected',
            'regulatory_precedent': 'established_covid_vaccines'
        }
        
        return safety_assessment
    
    def _assess_mrna_regulatory(self, mrna_sequence: str) -> Dict:
        """Assess regulatory considerations for mRNA vaccine."""
        
        regulatory = {
            'regulatory_pathway': 'biologics_license_application',
            'precedent': 'covid_mrna_vaccines',
            'key_requirements': [
                'sequence_characterization',
                'manufacturing_controls',
                'nonclinical_studies',
                'clinical_trials_phase_1_3'
            ],
            'timeline_estimate_months': 36,
            'major_regulatory_bodies': ['FDA', 'EMA', 'PMDA', 'Health_Canada'],
            'special_considerations': [
                'novel_antigen_characterization',
                'immunogenicity_assessment',
                'duration_of_immunity'
            ]
        }
        
        return regulatory
    
    def _design_viral_vector_vaccine(self, epitopes: List[Dict], params: Dict,
                                   disease_context: str) -> Dict:
        """Design viral vector vaccine."""
        
        # Select epitopes for viral vector
        selected_epitopes = self._select_epitopes_for_vector(epitopes, params)
        
        # Choose optimal vector
        vector_choice = self._select_optimal_vector(selected_epitopes, disease_context)
        
        # Design insert sequence
        insert_design = self._design_vector_insert(selected_epitopes, vector_choice)
        
        # Assess vector-specific considerations
        vector_assessment = self._assess_vector_considerations(vector_choice, insert_design)
        
        design = {
            'platform': 'viral_vector',
            'selected_epitopes': selected_epitopes,
            'vector_system': vector_choice,
            'insert_design': insert_design,
            'vector_assessment': vector_assessment,
            'manufacturing': self._calculate_vector_manufacturing(vector_choice, params),
            'predicted_efficacy': self._predict_vector_efficacy(selected_epitopes, vector_choice),
            'safety_profile': self._assess_vector_safety(vector_choice),
            'regulatory_considerations': self._assess_vector_regulatory(vector_choice)
        }
        
        return design
    
    def _select_epitopes_for_vector(self, epitopes: List[Dict], params: Dict) -> List[Dict]:
        """Select epitopes for viral vector vaccine."""
        
        # Sort by prediction score
        sorted_epitopes = sorted(epitopes, key=lambda x: x.get('prediction_score', 0), reverse=True)
        
        selected = []
        total_length = 0
        max_aa_length = params['optimal_insert_size'] // 3
        
        for epitope in sorted_epitopes:
            epitope_length = len(epitope.get('sequence', epitope.get('peptide', '')))
            
            if total_length + epitope_length <= max_aa_length:
                selected.append(epitope)
                total_length += epitope_length
                
                if len(selected) > 1:
                    total_length += 4  # Linker
        
        return selected
    
    def _select_optimal_vector(self, epitopes: List[Dict], disease_context: str) -> Dict:
        """Select optimal viral vector system."""
        
        # Vector scoring based on context
        vector_scores = {
            'adenovirus_5': {
                'immunogenicity': 0.9,
                'safety': 0.8,
                'manufacturing': 0.7,
                'regulatory': 0.9,
                'pre_existing_immunity': 0.3  # High pre-existing immunity
            },
            'adenovirus_26': {
                'immunogenicity': 0.85,
                'safety': 0.85,
                'manufacturing': 0.7,
                'regulatory': 0.8,
                'pre_existing_immunity': 0.7  # Lower pre-existing immunity
            },
            'modified_vaccinia_ankara': {
                'immunogenicity': 0.8,
                'safety': 0.9,
                'manufacturing': 0.6,
                'regulatory': 0.8,
                'pre_existing_immunity': 0.8
            }
        }
        
        # Select based on disease context
        if disease_context == 'cancer':
            # For cancer, prioritize immunogenicity
            selected_vector = 'adenovirus_5'
        else:
            # For infectious diseases, balance factors
            selected_vector = 'adenovirus_26'
        
        vector_info = {
            'vector_type': selected_vector,
            'properties': vector_scores[selected_vector],
            'capacity_bp': 8000,
            'tropism': 'broad',
            'replication': 'replication_deficient'
        }
        
        return vector_info
    
    def _design_vector_insert(self, epitopes: List[Dict], vector_choice: Dict) -> Dict:
        """Design insert sequence for viral vector."""
        
        # Construct insert with promoter, epitopes, and regulatory elements
        insert_elements = {
            'promoter': 'CMV_immediate_early',
            'signal_peptide': 'tissue_plasminogen_activator',
            'epitope_cassette': self._construct_epitope_cassette(epitopes),
            'kozak_sequence': 'GCCACC',
            'polya_signal': 'SV40_polya',
            'total_length_bp': 0
        }
        
        # Calculate total length
        estimated_length = 600 + len(insert_elements['epitope_cassette']) * 3  # Rough estimate
        insert_elements['total_length_bp'] = estimated_length
        
        return insert_elements
    
    def _construct_epitope_cassette(self, epitopes: List[Dict]) -> str:
        """Construct epitope cassette for vector insert."""
        
        epitope_sequences = []
        for epitope in epitopes:
            aa_seq = epitope.get('sequence', epitope.get('peptide', ''))
            epitope_sequences.append(aa_seq)
        
        # Join with GGSG linker
        cassette = 'GGSG'.join(epitope_sequences)
        
        return cassette
    
    def _assess_vector_considerations(self, vector_choice: Dict, insert_design: Dict) -> Dict:
        """Assess vector-specific considerations."""
        
        considerations = {
            'pre_existing_immunity': {
                'prevalence': vector_choice['properties']['pre_existing_immunity'],
                'impact_on_efficacy': 'moderate' if vector_choice['properties']['pre_existing_immunity'] < 0.5 else 'high',
                'mitigation_strategies': ['prime_boost', 'alternative_vectors']
            },
            'insert_stability': {
                'genetic_stability': 'good',
                'expression_level': 'high',
                'duration': 'transient'
            },
            'biodistribution': {
                'target_tissues': ['muscle', 'draining_lymph_nodes'],
                'systemic_exposure': 'limited',
                'clearance_mechanism': 'immune_mediated'
            }
        }
        
        return considerations
    
    def _calculate_vector_manufacturing(self, vector_choice: Dict, params: Dict) -> Dict:
        """Calculate viral vector manufacturing parameters."""
        
        manufacturing = {
            'production_system': 'HEK293_cells',
            'purification_method': 'cesium_chloride_gradient',
            'yield_particles_per_ml': 1e12,
            'purity_percentage': 90,
            'production_time_weeks': 4,
            'batch_size_doses': 50000,
            'cost_per_dose_usd': params['cost_per_dose'],
            'scalability': params['scalability'],
            'quality_control': [
                'vector_genome_titer',
                'infectious_titer',
                'purity_sds_page',
                'sterility',
                'endotoxin'
            ]
        }
        
        return manufacturing
    
    def _predict_vector_efficacy(self, epitopes: List[Dict], vector_choice: Dict) -> Dict:
        """Predict viral vector vaccine efficacy."""
        
        epitope_scores = [e.get('prediction_score', 0.5) for e in epitopes]
        base_efficacy = np.mean(epitope_scores) if epitope_scores else 0.5
        
        # Vector-specific bonus
        vector_bonus = vector_choice['properties']['immunogenicity'] * 0.2
        
        # Pre-existing immunity penalty
        immunity_penalty = (1 - vector_choice['properties']['pre_existing_immunity']) * 0.1
        
        overall_efficacy = min(1.0, base_efficacy + vector_bonus - immunity_penalty)
        
        return {
            'predicted_efficacy': overall_efficacy,
            'confidence_interval': [overall_efficacy - 0.15, overall_efficacy + 0.1],
            'expected_immune_response': {
                'antibody_response': 'strong',
                'cd8_response': 'very_strong',
                'cd4_response': 'strong',
                'memory_formation': 'excellent'
            }
        }
    
    def _assess_vector_safety(self, vector_choice: Dict) -> Dict:
        """Assess viral vector safety profile."""
        
        safety = {
            'overall_safety_score': vector_choice['properties']['safety'],
            'known_risks': [
                'injection_site_reactions',
                'flu_like_symptoms',
                'rare_thrombotic_events'
            ],
            'contraindications': [
                'severe_immunodeficiency',
                'pregnancy_caution'
            ],
            'monitoring_requirements': [
                'platelet_count',
                'coagulation_parameters'
            ]
        }
        
        return safety
    
    def _assess_vector_regulatory(self, vector_choice: Dict) -> Dict:
        """Assess viral vector regulatory considerations."""
        
        regulatory = {
            'regulatory_pathway': 'biologics_license_application',
            'precedent': 'covid_vector_vaccines',
            'timeline_estimate_months': 42,
            'special_requirements': [
                'vector_characterization',
                'biodistribution_studies',
                'shedding_studies'
            ]
        }
        
        return regulatory
    
    def _design_peptide_vaccine(self, epitopes: List[Dict], params: Dict,
                              disease_context: str) -> Dict:
        """Design peptide vaccine."""
        
        # Select optimal peptides
        selected_peptides = self._select_peptides_for_vaccine(epitopes, params)
        
        # Optimize peptide properties
        optimized_peptides = self._optimize_peptide_properties(selected_peptides)
        
        # Select adjuvant
        adjuvant_selection = self._select_optimal_adjuvant(optimized_peptides, disease_context)
        
        # Design formulation
        formulation_design = self._design_peptide_formulation(optimized_peptides, adjuvant_selection)
        
        design = {
            'platform': 'peptide',
            'selected_peptides': optimized_peptides,
            'adjuvant_system': adjuvant_selection,
            'formulation': formulation_design,
            'manufacturing': self._calculate_peptide_manufacturing(optimized_peptides, params),
            'predicted_efficacy': self._predict_peptide_efficacy(optimized_peptides, adjuvant_selection),
            'safety_profile': self._assess_peptide_safety(optimized_peptides, adjuvant_selection),
            'regulatory_considerations': self._assess_peptide_regulatory()
        }
        
        return design
    
    def _select_peptides_for_vaccine(self, epitopes: List[Dict], params: Dict) -> List[Dict]:
        """Select optimal peptides for vaccine."""
        
        # Sort by prediction score
        sorted_epitopes = sorted(epitopes, key=lambda x: x.get('prediction_score', 0), reverse=True)
        
        selected = []
        for epitope in sorted_epitopes[:params['optimal_peptides']]:
            peptide_length = len(epitope.get('sequence', epitope.get('peptide', '')))
            
            # Check length constraints
            if params['peptide_length_range'][0] <= peptide_length <= params['peptide_length_range'][1]:
                selected.append(epitope)
        
        return selected
    
    def _optimize_peptide_properties(self, peptides: List[Dict]) -> List[Dict]:
        """Optimize peptide properties for vaccine use."""
        
        optimized = []
        for peptide in peptides:
            sequence = peptide.get('sequence', peptide.get('peptide', ''))
            
            # Calculate properties
            properties = self._calculate_peptide_properties(sequence)
            
            # Add optimization recommendations
            optimizations = self._recommend_peptide_optimizations(sequence, properties)
            
            optimized_peptide = peptide.copy()
            optimized_peptide.update({
                'properties': properties,
                'optimizations': optimizations,
                'final_sequence': optimizations.get('optimized_sequence', sequence)
            })
            
            optimized.append(optimized_peptide)
        
        return optimized
    
    def _calculate_peptide_properties(self, sequence: str) -> Dict:
        """Calculate peptide physicochemical properties."""
        
        # Amino acid properties
        aa_properties = {
            'A': {'hydrophobic': 1.8, 'charge': 0, 'polar': False},
            'R': {'hydrophobic': -4.5, 'charge': 1, 'polar': True},
            'N': {'hydrophobic': -3.5, 'charge': 0, 'polar': True},
            'D': {'hydrophobic': -3.5, 'charge': -1, 'polar': True},
            'C': {'hydrophobic': 2.5, 'charge': 0, 'polar': False},
            'Q': {'hydrophobic': -3.5, 'charge': 0, 'polar': True},
            'E': {'hydrophobic': -3.5, 'charge': -1, 'polar': True},
            'G': {'hydrophobic': -0.4, 'charge': 0, 'polar': False},
            'H': {'hydrophobic': -3.2, 'charge': 0.5, 'polar': True},
            'I': {'hydrophobic': 4.5, 'charge': 0, 'polar': False},
            'L': {'hydrophobic': 3.8, 'charge': 0, 'polar': False},
            'K': {'hydrophobic': -3.9, 'charge': 1, 'polar': True},
            'M': {'hydrophobic': 1.9, 'charge': 0, 'polar': False},
            'F': {'hydrophobic': 2.8, 'charge': 0, 'polar': False},
            'P': {'hydrophobic': -1.6, 'charge': 0, 'polar': False},
            'S': {'hydrophobic': -0.8, 'charge': 0, 'polar': True},
            'T': {'hydrophobic': -0.7, 'charge': 0, 'polar': True},
            'W': {'hydrophobic': -0.9, 'charge': 0, 'polar': False},
            'Y': {'hydrophobic': -1.3, 'charge': 0, 'polar': True},
            'V': {'hydrophobic': 4.2, 'charge': 0, 'polar': False}
        }
        
        # Calculate properties
        hydrophobicity = sum(aa_properties.get(aa, {'hydrophobic': 0})['hydrophobic'] for aa in sequence) / len(sequence)
        net_charge = sum(aa_properties.get(aa, {'charge': 0})['charge'] for aa in sequence)
        polar_residues = sum(1 for aa in sequence if aa_properties.get(aa, {'polar': False})['polar'])
        
        properties = {
            'length': len(sequence),
            'molecular_weight_da': len(sequence) * 110,  # Approximate
            'hydrophobicity': hydrophobicity,
            'net_charge': net_charge,
            'polar_residue_fraction': polar_residues / len(sequence),
            'solubility_prediction': 'good' if hydrophobicity < 0 else 'moderate' if hydrophobicity < 2 else 'poor',
            'stability_prediction': 'stable' if abs(net_charge) < 3 else 'moderate',
            'aggregation_propensity': 'low' if hydrophobicity < 1 else 'moderate' if hydrophobicity < 3 else 'high'
        }
        
        return properties
    
    def _recommend_peptide_optimizations(self, sequence: str, properties: Dict) -> Dict:
        """Recommend peptide optimizations."""
        
        optimizations = {
            'optimized_sequence': sequence,
            'modifications': [],
            'rationale': []
        }
        
        # Solubility optimization
        if properties['solubility_prediction'] == 'poor':
            optimizations['modifications'].append('add_solubilizing_tag')
            optimizations['rationale'].append('Improve solubility with hydrophilic tag')
        
        # Stability optimization
        if properties['stability_prediction'] == 'moderate':
            optimizations['modifications'].append('cyclization')
            optimizations['rationale'].append('Cyclization to improve stability')
        
        # Aggregation prevention
        if properties['aggregation_propensity'] == 'high':
            optimizations['modifications'].append('spacer_residues')
            optimizations['rationale'].append('Add spacer residues to prevent aggregation')
        
        return optimizations
    
    def _select_optimal_adjuvant(self, peptides: List[Dict], disease_context: str) -> Dict:
        """Select optimal adjuvant for peptide vaccine."""
        
        # Score adjuvants based on context
        adjuvant_scores = {}
        
        for adj_name, adj_data in self.adjuvants.items():
            if 'peptide' in adj_data['compatibility']:
                score = 0.5  # Base score
                
                # Disease context bonus
                if disease_context == 'cancer' and adj_data['immune_response'] in ['Th1_dominant', 'Th1/Th2_balanced']:
                    score += 0.3
                elif disease_context == 'infectious' and adj_data['immune_response'] == 'Th1/Th2_balanced':
                    score += 0.2
                
                # Safety bonus
                if adj_data['safety_profile'] == 'excellent':
                    score += 0.2
                elif adj_data['safety_profile'] == 'good':
                    score += 0.1
                
                # Cost consideration
                if adj_data['cost'] == 'low':
                    score += 0.1
                elif adj_data['cost'] == 'medium':
                    score += 0.05
                
                adjuvant_scores[adj_name] = score
        
        # Select best adjuvant
        best_adjuvant = max(adjuvant_scores.items(), key=lambda x: x[1])
        
        selected_adjuvant = {
            'adjuvant_name': best_adjuvant[0],
            'selection_score': best_adjuvant[1],
            'properties': self.adjuvants[best_adjuvant[0]],
            'dosage_recommendation': self._recommend_adjuvant_dosage(best_adjuvant[0]),
            'administration_route': 'intramuscular'
        }
        
        return selected_adjuvant
    
    def _recommend_adjuvant_dosage(self, adjuvant_name: str) -> Dict:
        """Recommend adjuvant dosage."""
        
        dosage_recommendations = {
            'alum': {'amount_mg': 0.5, 'volume_ml': 0.5},
            'MF59': {'amount_ml': 0.5, 'volume_ml': 0.5},
            'AS01': {'amount_ml': 0.5, 'volume_ml': 0.5},
            'CpG_ODN': {'amount_ug': 1000, 'volume_ml': 0.5},
            'poly_IC': {'amount_ug': 50, 'volume_ml': 0.5}
        }
        
        return dosage_recommendations.get(adjuvant_name, {'amount': 'TBD', 'volume_ml': 0.5})
    
    def _design_peptide_formulation(self, peptides: List[Dict], adjuvant: Dict) -> Dict:
        """Design peptide vaccine formulation."""
        
        formulation = {
            'peptide_concentration_mg_ml': 1.0,
            'total_peptide_dose_mg': 0.1,
            'adjuvant_system': adjuvant,
            'buffer_system': 'phosphate_buffered_saline',
            'ph': 7.4,
            'osmolality_mosm_kg': 300,
            'preservative': 'none',
            'stabilizers': ['sucrose', 'polysorbate_80'],
            'storage_conditions': {
                'temperature_c': 4,
                'light_protection': True,
                'shelf_life_months': 24
            },
            'administration': {
                'route': 'intramuscular',
                'volume_ml': 0.5,
                'injection_site': 'deltoid_muscle'
            }
        }
        
        return formulation
    
    def _calculate_peptide_manufacturing(self, peptides: List[Dict], params: Dict) -> Dict:
        """Calculate peptide manufacturing parameters."""
        
        manufacturing = {
            'synthesis_method': 'solid_phase_peptide_synthesis',
            'purification_method': 'reverse_phase_hplc',
            'yield_percentage': 75,
            'purity_percentage': 98,
            'production_time_weeks': 2,
            'batch_size_doses': 100000,
            'cost_per_dose_usd': params['cost_per_dose'],
            'scalability': params['scalability'],
            'quality_control': [
                'mass_spectrometry',
                'amino_acid_analysis',
                'hplc_purity',
                'endotoxin_testing',
                'sterility'
            ]
        }
        
        return manufacturing
    
    def _predict_peptide_efficacy(self, peptides: List[Dict], adjuvant: Dict) -> Dict:
        """Predict peptide vaccine efficacy."""
        
        peptide_scores = [p.get('prediction_score', 0.5) for p in peptides]
        base_efficacy = np.mean(peptide_scores) if peptide_scores else 0.5
        
        # Adjuvant bonus
        adjuvant_bonus = adjuvant['selection_score'] * 0.15
        
        # Multiple peptide bonus
        multi_peptide_bonus = min(0.1, len(peptides) * 0.02)
        
        overall_efficacy = min(1.0, base_efficacy + adjuvant_bonus + multi_peptide_bonus)
        
        return {
            'predicted_efficacy': overall_efficacy,
            'confidence_interval': [overall_efficacy - 0.12, overall_efficacy + 0.08],
            'expected_immune_response': {
                'antibody_response': 'moderate',
                'cd8_response': 'strong',
                'cd4_response': 'strong',
                'memory_formation': 'good'
            }
        }
    
    def _assess_peptide_safety(self, peptides: List[Dict], adjuvant: Dict) -> Dict:
        """Assess peptide vaccine safety."""
        
        safety = {
            'overall_safety_score': 0.9,  # Peptides generally very safe
            'known_risks': [
                'injection_site_reactions',
                'adjuvant_related_reactions'
            ],
            'contraindications': [
                'known_peptide_allergies',
                'adjuvant_specific_contraindications'
            ],
            'special_populations': {
                'pregnancy': 'generally_safe',
                'pediatric': 'safe_with_dose_adjustment',
                'elderly': 'safe',
                'immunocompromised': 'safe_but_reduced_efficacy'
            }
        }
        
        return safety
    
    def _assess_peptide_regulatory(self) -> Dict:
        """Assess peptide vaccine regulatory considerations."""
        
        regulatory = {
            'regulatory_pathway': 'biologics_license_application',
            'precedent': 'established_peptide_vaccines',
            'timeline_estimate_months': 30,
            'key_requirements': [
                'peptide_characterization',
                'adjuvant_qualification',
                'nonclinical_studies',
                'clinical_trials'
            ]
        }
        
        return regulatory
    
    def _design_protein_subunit_vaccine(self, epitopes: List[Dict], params: Dict,
                                      disease_context: str) -> Dict:
        """Design protein subunit vaccine."""
        
        # Design protein construct
        protein_design = self._design_protein_construct(epitopes, params)
        
        # Select expression system
        expression_system = self._select_expression_system(protein_design, disease_context)
        
        # Design purification strategy
        purification_strategy = self._design_purification_strategy(protein_design, expression_system)
        
        # Select adjuvant
        adjuvant_selection = self._select_optimal_adjuvant(epitopes, disease_context)
        
        design = {
            'platform': 'protein_subunit',
            'protein_design': protein_design,
            'expression_system': expression_system,
            'purification_strategy': purification_strategy,
            'adjuvant_system': adjuvant_selection,
            'manufacturing': self._calculate_protein_manufacturing(protein_design, params),
            'predicted_efficacy': self._predict_protein_efficacy(protein_design, adjuvant_selection),
            'safety_profile': self._assess_protein_safety(protein_design, adjuvant_selection),
            'regulatory_considerations': self._assess_protein_regulatory()
        }
        
        return design
    
    def _design_protein_construct(self, epitopes: List[Dict], params: Dict) -> Dict:
        """Design protein construct for subunit vaccine."""
        
        # Select epitopes for protein design
        selected_epitopes = epitopes[:params['optimal_size'] // 15]  # Rough estimate
        
        # Design protein sequence
        protein_sequence = self._construct_protein_sequence(selected_epitopes)
        
        # Add tags and signals
        construct_design = {
            'signal_peptide': 'native_or_optimized',
            'epitope_regions': selected_epitopes,
            'linker_sequences': ['GGSG', 'EAAAK'],
            'purification_tag': 'his6_tag',
            'total_length_aa': len(protein_sequence),
            'molecular_weight_kda': len(protein_sequence) * 0.11,
            'predicted_structure': 'mixed_alpha_beta',
            'solubility_prediction': 'good'
        }
        
        return construct_design
    
    def _construct_protein_sequence(self, epitopes: List[Dict]) -> str:
        """Construct protein sequence from epitopes."""
        
        epitope_sequences = []
        for epitope in epitopes:
            seq = epitope.get('sequence', epitope.get('peptide', ''))
            epitope_sequences.append(seq)
        
        # Join with flexible linker
        protein_sequence = 'GGSG'.join(epitope_sequences)
        
        return protein_sequence
    
    def _select_expression_system(self, protein_design: Dict, disease_context: str) -> Dict:
        """Select optimal expression system."""
        
        expression_systems = {
            'E.coli': {
                'cost': 'low',
                'speed': 'fast',
                'yield': 'high',
                'post_translational_modifications': False,
                'endotoxin_risk': True,
                'scalability': 'excellent'
            },
            'yeast': {
                'cost': 'medium',
                'speed': 'medium',
                'yield': 'medium',
                'post_translational_modifications': True,
                'endotoxin_risk': False,
                'scalability': 'good'
            },
            'mammalian': {
                'cost': 'high',
                'speed': 'slow',
                'yield': 'low',
                'post_translational_modifications': True,
                'endotoxin_risk': False,
                'scalability': 'limited'
            }
        }
        
        # Select based on protein requirements
        if protein_design['molecular_weight_kda'] < 30:
            selected_system = 'E.coli'
        else:
            selected_system = 'yeast'
        
        return {
            'system': selected_system,
            'properties': expression_systems[selected_system],
            'strain_recommendation': 'BL21_DE3' if selected_system == 'E.coli' else 'Pichia_pastoris'
        }
    
    def _design_purification_strategy(self, protein_design: Dict, expression_system: Dict) -> Dict:
        """Design protein purification strategy."""
        
        purification_steps = []
        
        # Initial capture
        if 'his6_tag' in protein_design.get('purification_tag', ''):
            purification_steps.append({
                'step': 'immobilized_metal_affinity_chromatography',
                'resin': 'ni_nta',
                'expected_purity': 85,
                'yield': 80
            })
        
        # Polishing steps
        purification_steps.extend([
            {
                'step': 'size_exclusion_chromatography',
                'purpose': 'aggregate_removal',
                'expected_purity': 95,
                'yield': 90
            },
            {
                'step': 'endotoxin_removal',
                'method': 'polymyxin_b_chromatography',
                'target_level': '<1_EU_mg',
                'yield': 95
            }
        ])
        
        strategy = {
            'purification_steps': purification_steps,
            'overall_yield': 68,  # Product of individual yields
            'final_purity': 95,
            'endotoxin_level': '<1_EU_mg',
            'process_time_days': 3
        }
        
        return strategy
    
    def _calculate_protein_manufacturing(self, protein_design: Dict, params: Dict) -> Dict:
        """Calculate protein subunit manufacturing parameters."""
        
        manufacturing = {
            'fermentation_time_hours': 48,
            'harvest_yield_g_l': 2.0,
            'purification_yield_percentage': 68,
            'final_concentration_mg_ml': 10,
            'batch_size_doses': 75000,
            'production_time_weeks': 3,
            'cost_per_dose_usd': params['cost_per_dose'],
            'scalability': params['scalability'],
            'quality_control': [
                'sds_page_purity',
                'mass_spectrometry',
                'endotoxin_lal',
                'bioactivity_assay',
                'sterility'
            ]
        }
        
        return manufacturing
    
    def _predict_protein_efficacy(self, protein_design: Dict, adjuvant: Dict) -> Dict:
        """Predict protein subunit vaccine efficacy."""
        
        # Base efficacy from protein design
        base_efficacy = 0.6  # Typical for protein subunits
        
        # Adjuvant contribution
        adjuvant_bonus = adjuvant['selection_score'] * 0.2
        
        # Protein size bonus (larger proteins often more immunogenic)
        size_bonus = min(0.1, protein_design['molecular_weight_kda'] / 500)
        
        overall_efficacy = min(1.0, base_efficacy + adjuvant_bonus + size_bonus)
        
        return {
            'predicted_efficacy': overall_efficacy,
            'confidence_interval': [overall_efficacy - 0.1, overall_efficacy + 0.1],
            'expected_immune_response': {
                'antibody_response': 'strong',
                'cd8_response': 'moderate',
                'cd4_response': 'strong',
                'memory_formation': 'excellent'
            }
        }
    
    def _assess_protein_safety(self, protein_design: Dict, adjuvant: Dict) -> Dict:
        """Assess protein subunit vaccine safety."""
        
        safety = {
            'overall_safety_score': 0.88,
            'known_risks': [
                'injection_site_reactions',
                'adjuvant_related_reactions',
                'rare_allergic_reactions'
            ],
            'contraindications': [
                'known_protein_allergies',
                'adjuvant_contraindications'
            ],
            'special_populations': {
                'pregnancy': 'generally_safe',
                'pediatric': 'safe',
                'elderly': 'safe',
                'immunocompromised': 'safe_but_reduced_response'
            }
        }
        
        return safety
    
    def _assess_protein_regulatory(self) -> Dict:
        """Assess protein subunit vaccine regulatory considerations."""
        
        regulatory = {
            'regulatory_pathway': 'biologics_license_application',
            'precedent': 'established_protein_vaccines',
            'timeline_estimate_months': 33,
            'key_requirements': [
                'protein_characterization',
                'manufacturing_validation',
                'nonclinical_studies',
                'clinical_trials'
            ]
        }
        
        return regulatory
    
    def _design_dna_vaccine(self, epitopes: List[Dict], params: Dict,
                          disease_context: str) -> Dict:
        """Design DNA vaccine."""
        
        # Design DNA construct
        dna_construct = self._design_dna_construct(epitopes, params)
        
        # Select delivery method
        delivery_method = self._select_dna_delivery_method(disease_context)
        
        # Design formulation
        formulation = self._design_dna_formulation(dna_construct, delivery_method)
        
        design = {
            'platform': 'DNA',
            'dna_construct': dna_construct,
            'delivery_method': delivery_method,
            'formulation': formulation,
            'manufacturing': self._calculate_dna_manufacturing(dna_construct, params),
            'predicted_efficacy': self._predict_dna_efficacy(dna_construct, delivery_method),
            'safety_profile': self._assess_dna_safety(dna_construct),
            'regulatory_considerations': self._assess_dna_regulatory()
        }
        
        return design
    
    def _design_dna_construct(self, epitopes: List[Dict], params: Dict) -> Dict:
        """Design DNA vaccine construct."""
        
        # Select epitopes
        selected_epitopes = epitopes[:params['optimal_length'] // 100]  # Rough estimate
        
        # Design construct elements
        construct = {
            'promoter': 'CMV_immediate_early',
            'kozak_sequence': 'GCCACC',
            'signal_peptide': 'tissue_plasminogen_activator',
            'epitope_cassette': self._construct_epitope_cassette(selected_epitopes),
            'polya_signal': 'bovine_growth_hormone_polya',
            'backbone': 'pVAX1_derived',
            'antibiotic_resistance': 'kanamycin',
            'total_size_bp': 4500,  # Estimated
            'copy_number': 'high'
        }
        
        return construct
    
    def _select_dna_delivery_method(self, disease_context: str) -> Dict:
        """Select DNA vaccine delivery method."""
        
        delivery_methods = {
            'electroporation': {
                'efficiency': 'high',
                'invasiveness': 'moderate',
                'equipment_required': True,
                'cost': 'medium'
            },
            'gene_gun': {
                'efficiency': 'moderate',
                'invasiveness': 'low',
                'equipment_required': True,
                'cost': 'high'
            },
            'lipofection': {
                'efficiency': 'moderate',
                'invasiveness': 'low',
                'equipment_required': False,
                'cost': 'low'
            }
        }
        
        # Select based on context
        if disease_context == 'cancer':
            selected_method = 'electroporation'  # Higher efficiency needed
        else:
            selected_method = 'lipofection'  # Simpler administration
        
        return {
            'method': selected_method,
            'properties': delivery_methods[selected_method],
            'administration_details': self._get_delivery_details(selected_method)
        }
    
    def _get_delivery_details(self, method: str) -> Dict:
        """Get delivery method details."""
        
        details = {
            'electroporation': {
                'voltage': '100V',
                'pulse_duration': '5ms',
                'number_of_pulses': 3,
                'electrode_type': 'needle'
            },
            'gene_gun': {
                'particle_size': '1-3_microns',
                'pressure': '400_psi',
                'target_depth': 'epidermis'
            },
            'lipofection': {
                'lipid_formulation': 'cationic_liposomes',
                'dna_lipid_ratio': '1:3',
                'incubation_time': '30_minutes'
            }
        }
        
        return details.get(method, {})
    
    def _design_dna_formulation(self, dna_construct: Dict, delivery_method: Dict) -> Dict:
        """Design DNA vaccine formulation."""
        
        formulation = {
            'dna_concentration_mg_ml': 1.0,
            'dose_mg': 0.1,
            'buffer': 'tris_edta_ph_8.0',
            'stabilizers': ['sucrose', 'trehalose'],
            'preservative': 'none',
            'storage_conditions': {
                'temperature_c': 4,
                'light_protection': True,
                'shelf_life_months': 36
            }
        }
        
        # Add delivery-specific components
        if delivery_method['method'] == 'lipofection':
            formulation['lipid_components'] = {
                'cationic_lipid': 'DOTAP',
                'helper_lipid': 'DOPE',
                'ratio': '1:1'
            }
        
        return formulation
    
    def _calculate_dna_manufacturing(self, dna_construct: Dict, params: Dict) -> Dict:
        """Calculate DNA vaccine manufacturing parameters."""
        
        manufacturing = {
            'production_method': 'bacterial_fermentation',
            'host_strain': 'E.coli_DH5alpha',
            'plasmid_yield_mg_l': 50,
            'purification_method': 'alkaline_lysis_chromatography',
            'purity_percentage': 99,
            'endotoxin_level_eu_mg': 10,
            'production_time_weeks': 1,
            'batch_size_doses': 500000,
            'cost_per_dose_usd': params['cost_per_dose'],
            'scalability': params['scalability']
        }
        
        return manufacturing
    
    def _predict_dna_efficacy(self, dna_construct: Dict, delivery_method: Dict) -> Dict:
        """Predict DNA vaccine efficacy."""
        
        base_efficacy = 0.45  # DNA vaccines typically lower efficacy
        
        # Delivery method bonus
        delivery_bonus = 0.2 if delivery_method['method'] == 'electroporation' else 0.1
        
        # Construct design bonus
        construct_bonus = 0.1 if 'CMV' in dna_construct.get('promoter', '') else 0.05
        
        overall_efficacy = min(1.0, base_efficacy + delivery_bonus + construct_bonus)
        
        return {
            'predicted_efficacy': overall_efficacy,
            'confidence_interval': [overall_efficacy - 0.15, overall_efficacy + 0.1],
            'expected_immune_response': {
                'antibody_response': 'moderate',
                'cd8_response': 'strong',
                'cd4_response': 'moderate',
                'memory_formation': 'good'
            }
        }
    
    def _assess_dna_safety(self, dna_construct: Dict) -> Dict:
        """Assess DNA vaccine safety."""
        
        safety = {
            'overall_safety_score': 0.82,
            'known_risks': [
                'injection_site_reactions',
                'theoretical_integration_risk',
                'anti_dna_antibodies'
            ],
            'contraindications': [
                'pregnancy_caution',
                'immunodeficiency'
            ],
            'special_considerations': [
                'integration_monitoring',
                'long_term_follow_up'
            ]
        }
        
        return safety
    
    def _assess_dna_regulatory(self) -> Dict:
        """Assess DNA vaccine regulatory considerations."""
        
        regulatory = {
            'regulatory_pathway': 'investigational_new_drug',
            'precedent': 'limited_approved_dna_vaccines',
            'timeline_estimate_months': 48,
            'special_requirements': [
                'integration_studies',
                'biodistribution_studies',
                'long_term_safety_monitoring'
            ]
        }
        
        return regulatory
    
    def _compare_platforms(self, vaccine_designs: Dict, disease_context: str) -> Dict:
        """Compare vaccine platforms."""
        
        comparison_metrics = {}
        
        for platform, design in vaccine_designs.items():
            if 'error' not in design:
                metrics = {
                    'predicted_efficacy': design.get('predicted_efficacy', {}).get('predicted_efficacy', 0),
                    'safety_score': design.get('safety_profile', {}).get('overall_safety_score', 0),
                    'manufacturing_complexity': self._score_manufacturing_complexity(design.get('manufacturing', {})),
                    'cost_per_dose': design.get('manufacturing', {}).get('cost_per_dose_usd', 0),
                    'development_timeline': design.get('regulatory_considerations', {}).get('timeline_estimate_months', 36),
                    'scalability_score': self._score_scalability(design.get('manufacturing', {}))
                }
                
                # Calculate overall score
                weights = {
                    'predicted_efficacy': 0.3,
                    'safety_score': 0.25,
                    'manufacturing_complexity': 0.15,
                    'cost_per_dose': 0.1,
                    'development_timeline': 0.1,
                    'scalability_score': 0.1
                }
                
                # Normalize cost and timeline (lower is better)
                normalized_cost = max(0, 1 - (metrics['cost_per_dose'] / 50))
                normalized_timeline = max(0, 1 - (metrics['development_timeline'] / 60))
                
                overall_score = (
                    weights['predicted_efficacy'] * metrics['predicted_efficacy'] +
                    weights['safety_score'] * metrics['safety_score'] +
                    weights['manufacturing_complexity'] * metrics['manufacturing_complexity'] +
                    weights['cost_per_dose'] * normalized_cost +
                    weights['development_timeline'] * normalized_timeline +
                    weights['scalability_score'] * metrics['scalability_score']
                )
                
                metrics['overall_score'] = overall_score
                comparison_metrics[platform] = metrics
        
        # Rank platforms
        ranked_platforms = sorted(comparison_metrics.items(), 
                                key=lambda x: x[1]['overall_score'], reverse=True)
        
        return {
            'platform_metrics': comparison_metrics,
            'ranked_platforms': [p[0] for p in ranked_platforms],
            'top_recommendation': ranked_platforms[0][0] if ranked_platforms else None,
            'comparison_summary': self._generate_comparison_summary(comparison_metrics)
        }
    
    def _score_manufacturing_complexity(self, manufacturing: Dict) -> float:
        """Score manufacturing complexity (higher is simpler)."""
        
        complexity_map = {
            'low': 0.9,
            'medium': 0.6,
            'high': 0.3
        }
        
        complexity = manufacturing.get('manufacturing_complexity', 'medium')
        return complexity_map.get(complexity, 0.5)
    
    def _score_scalability(self, manufacturing: Dict) -> float:
        """Score scalability (higher is more scalable)."""
        
        scalability_map = {
            'very_high': 1.0,
            'high': 0.8,
            'medium': 0.6,
            'low': 0.4,
            'limited': 0.2
        }
        
        scalability = manufacturing.get('scalability', 'medium')
        return scalability_map.get(scalability, 0.5)
    
    def _generate_comparison_summary(self, metrics: Dict) -> Dict:
        """Generate platform comparison summary."""
        
        if not metrics:
            return {}
        
        # Find best in each category
        best_efficacy = max(metrics.items(), key=lambda x: x[1]['predicted_efficacy'])
        best_safety = max(metrics.items(), key=lambda x: x[1]['safety_score'])
        lowest_cost = min(metrics.items(), key=lambda x: x[1]['cost_per_dose'])
        fastest_development = min(metrics.items(), key=lambda x: x[1]['development_timeline'])
        
        summary = {
            'best_efficacy': {'platform': best_efficacy[0], 'score': best_efficacy[1]['predicted_efficacy']},
            'best_safety': {'platform': best_safety[0], 'score': best_safety[1]['safety_score']},
            'lowest_cost': {'platform': lowest_cost[0], 'cost': lowest_cost[1]['cost_per_dose']},
            'fastest_development': {'platform': fastest_development[0], 'months': fastest_development[1]['development_timeline']},
            'average_efficacy': np.mean([m['predicted_efficacy'] for m in metrics.values()]),
            'average_safety': np.mean([m['safety_score'] for m in metrics.values()]),
            'cost_range': [min(m['cost_per_dose'] for m in metrics.values()),
                          max(m['cost_per_dose'] for m in metrics.values())]
        }
        
        return summary
    
    def _generate_platform_recommendations(self, vaccine_designs: Dict,
                                         platform_comparison: Dict,
                                         disease_context: str) -> List[Dict]:
        """Generate platform-specific recommendations."""
        
        recommendations = []
        
        # Top platform recommendation
        if platform_comparison.get('top_recommendation'):
            top_platform = platform_comparison['top_recommendation']
            recommendations.append({
                'type': 'primary_recommendation',
                'platform': top_platform,
                'priority': 'high',
                'rationale': f'{top_platform} shows the best overall profile for {disease_context}',
                'next_steps': [
                    'Conduct detailed feasibility study',
                    'Initiate preclinical development',
                    'Engage with regulatory authorities'
                ]
            })
        
        # Alternative platform recommendations
        ranked_platforms = platform_comparison.get('ranked_platforms', [])
        if len(ranked_platforms) > 1:
            second_choice = ranked_platforms[1]
            recommendations.append({
                'type': 'alternative_recommendation',
                'platform': second_choice,
                'priority': 'medium',
                'rationale': f'{second_choice} as backup or parallel development option',
                'considerations': [
                    'Risk mitigation strategy',
                    'Different regulatory pathway',
                    'Complementary immune response profile'
                ]
            })
        
        # Context-specific recommendations
        if disease_context == 'cancer':
            recommendations.append({
                'type': 'context_specific',
                'recommendation': 'Consider combination approaches',
                'rationale': 'Cancer vaccines often benefit from multi-modal approaches',
                'suggestions': [
                    'Prime-boost strategies',
                    'Combination with checkpoint inhibitors',
                    'Personalized neoantigen approaches'
                ]
            })
        
        # Manufacturing recommendations
        recommendations.append({
            'type': 'manufacturing',
            'recommendation': 'Establish manufacturing partnerships early',
            'rationale': 'Manufacturing capabilities are critical for success',
            'actions': [
                'Identify qualified manufacturers',
                'Negotiate technology transfer agreements',
                'Plan for scale-up requirements'
            ]
        })
        
        return recommendations
    
    def _assess_manufacturing_feasibility(self, vaccine_designs: Dict) -> Dict:
        """Assess manufacturing feasibility across platforms."""
        
        feasibility_assessment = {
            'overall_feasibility': 'good',
            'platform_assessments': {},
            'critical_factors': [],
            'recommendations': []
        }
        
        for platform, design in vaccine_designs.items():
            if 'error' not in design:
                manufacturing = design.get('manufacturing', {})
                
                assessment = {
                    'complexity_score': self._score_manufacturing_complexity(manufacturing),
                    'scalability_score': self._score_scalability(manufacturing),
                    'cost_competitiveness': self._assess_cost_competitiveness(manufacturing),
                    'timeline_feasibility': self._assess_timeline_feasibility(manufacturing),
                    'risk_factors': self._identify_manufacturing_risks(platform, manufacturing)
                }
                
                feasibility_assessment['platform_assessments'][platform] = assessment
        
        # Identify critical factors
        feasibility_assessment['critical_factors'] = [
            'Supply chain reliability',
            'Quality control capabilities',
            'Regulatory compliance',
            'Scale-up challenges',
            'Cost optimization'
        ]
        
        # Generate recommendations
        feasibility_assessment['recommendations'] = [
            'Conduct detailed manufacturing feasibility studies',
            'Establish quality by design principles',
            'Plan for regulatory inspections',
            'Develop robust supply chain strategies',
            'Implement continuous improvement processes'
        ]
        
        return feasibility_assessment
    
    def _assess_cost_competitiveness(self, manufacturing: Dict) -> str:
        """Assess cost competitiveness."""
        
        cost_per_dose = manufacturing.get('cost_per_dose_usd', 25)
        
        if cost_per_dose < 15:
            return 'highly_competitive'
        elif cost_per_dose < 30:
            return 'competitive'
        elif cost_per_dose < 50:
            return 'moderate'
        else:
            return 'expensive'
    
    def _assess_timeline_feasibility(self, manufacturing: Dict) -> str:
        """Assess timeline feasibility."""
        
        production_time = manufacturing.get('production_time_weeks', 4)
        
        if production_time < 2:
            return 'very_fast'
        elif production_time < 4:
            return 'fast'
        elif production_time < 8:
            return 'moderate'
        else:
            return 'slow'
    
    def _identify_manufacturing_risks(self, platform: str, manufacturing: Dict) -> List[str]:
        """Identify platform-specific manufacturing risks."""
        
        risk_profiles = {
            'mRNA': [
                'Cold chain requirements',
                'Lipid nanoparticle stability',
                'RNA degradation',
                'Specialized equipment needs'
            ],
            'viral_vector': [
                'Cell culture contamination',
                'Vector stability',
                'Purification complexity',
                'Batch-to-batch variability'
            ],
            'peptide': [
                'Peptide aggregation',
                'Chemical stability',
                'Adjuvant compatibility',
                'Synthesis scale-up'
            ],
            'protein_subunit': [
                'Protein folding',
                'Endotoxin contamination',
                'Expression system limitations',
                'Purification yield'
            ],
            'DNA': [
                'Plasmid stability',
                'Endotoxin removal',
                'Supercoiling maintenance',
                'Delivery system integration'
            ]
        }
        
        return risk_profiles.get(platform, ['General manufacturing risks'])
    
    def generate_delivery_system_report(self, analysis_results: Dict,
                                      output_path: str) -> str:
        """
        Generate comprehensive delivery system integration report.
        
        Args:
            analysis_results: Results from delivery system analysis
            output_path: Output file path
            
        Returns:
            Path to generated report
        """
        logger.info("Generating vaccine delivery system integration report")
        
        report = "# Vaccine Delivery System Integration Report\n\n"
        
        # Executive Summary
        report += "## Executive Summary\n\n"
        
        if 'platform_comparison' in analysis_results:
            comparison = analysis_results['platform_comparison']
            top_platform = comparison.get('top_recommendation', 'Not determined')
            report += f"**Recommended Platform:** {top_platform.replace('_', ' ').title()}\n\n"
            
            if 'comparison_summary' in comparison:
                summary = comparison['comparison_summary']
                report += f"**Best Efficacy:** {summary.get('best_efficacy', {}).get('platform', 'N/A')} "
                report += f"({summary.get('best_efficacy', {}).get('score', 0):.1%})\n"
                report += f"**Best Safety:** {summary.get('best_safety', {}).get('platform', 'N/A')} "
                report += f"(Score: {summary.get('best_safety', {}).get('score', 0):.2f})\n"
                report += f"**Lowest Cost:** {summary.get('lowest_cost', {}).get('platform', 'N/A')} "
                report += f"(${summary.get('lowest_cost', {}).get('cost', 0):.2f}/dose)\n\n"
        
        # Platform Designs
        if 'vaccine_designs' in analysis_results:
            report += "## Platform-Specific Vaccine Designs\n\n"
            
            for platform, design in analysis_results['vaccine_designs'].items():
                if 'error' not in design:
                    report += f"### {platform.replace('_', ' ').title()} Vaccine\n\n"
                    
                    # Efficacy
                    if 'predicted_efficacy' in design:
                        efficacy = design['predicted_efficacy']
                        report += f"**Predicted Efficacy:** {efficacy.get('predicted_efficacy', 0):.1%}\n"
                        
                        immune_response = efficacy.get('expected_immune_response', {})
                        report += f"- Antibody Response: {immune_response.get('antibody_response', 'N/A').title()}\n"
                        report += f"- CD8+ T-cell Response: {immune_response.get('cd8_response', 'N/A').title()}\n"
                        report += f"- CD4+ T-cell Response: {immune_response.get('cd4_response', 'N/A').title()}\n"
                        report += f"- Memory Formation: {immune_response.get('memory_formation', 'N/A').title()}\n\n"
                    
                    # Safety
                    if 'safety_profile' in design:
                        safety = design['safety_profile']
                        report += f"**Safety Score:** {safety.get('overall_safety_score', 0):.2f}/1.0\n"
                        
                        known_risks = safety.get('known_risks', [])
                        if known_risks:
                            report += "**Known Risks:**\n"
                            for risk in known_risks:
                                report += f"- {risk.replace('_', ' ').title()}\n"
                        report += "\n"
                    
                    # Manufacturing
                    if 'manufacturing' in design:
                        manufacturing = design['manufacturing']
                        report += f"**Manufacturing:**\n"
                        report += f"- Cost per dose: ${manufacturing.get('cost_per_dose_usd', 0):.2f}\n"
                        report += f"- Production time: {manufacturing.get('production_time_weeks', 'N/A')} weeks\n"
                        report += f"- Batch size: {manufacturing.get('batch_size_doses', 'N/A'):,} doses\n"
                        report += f"- Scalability: {manufacturing.get('scalability', 'N/A').replace('_', ' ').title()}\n\n"
        
        # Platform Comparison
        if 'platform_comparison' in analysis_results:
            comparison = analysis_results['platform_comparison']
            
            if 'platform_metrics' in comparison:
                report += "## Platform Comparison\n\n"
                report += "| Platform | Efficacy | Safety | Cost/Dose | Timeline | Overall Score |\n"
                report += "|----------|----------|--------|-----------|----------|---------------|\n"
                
                for platform, metrics in comparison['platform_metrics'].items():
                    report += f"| {platform.replace('_', ' ').title()} | "
                    report += f"{metrics.get('predicted_efficacy', 0):.1%} | "
                    report += f"{metrics.get('safety_score', 0):.2f} | "
                    report += f"${metrics.get('cost_per_dose', 0):.2f} | "
                    report += f"{metrics.get('development_timeline', 0)} mo | "
                    report += f"{metrics.get('overall_score', 0):.3f} |\n"
                report += "\n"
        
        # Recommendations
        if 'recommendations' in analysis_results:
            report += "## Recommendations\n\n"
            
            for i, rec in enumerate(analysis_results['recommendations'], 1):
                report += f"### {i}. {rec.get('recommendation', 'Recommendation')}\n\n"
                report += f"**Platform:** {rec.get('platform', 'N/A').replace('_', ' ').title()}\n"
                report += f"**Priority:** {rec.get('priority', 'Medium').title()}\n"
                report += f"**Rationale:** {rec.get('rationale', 'Not specified')}\n"
                
                if 'next_steps' in rec:
                    report += "**Next Steps:**\n"
                    for step in rec['next_steps']:
                        report += f"- {step}\n"
                
                if 'considerations' in rec:
                    report += "**Considerations:**\n"
                    for consideration in rec['considerations']:
                        report += f"- {consideration}\n"
                
                if 'actions' in rec:
                    report += "**Actions:**\n"
                    for action in rec['actions']:
                        report += f"- {action}\n"
                
                report += "\n"
        
        # Manufacturing Feasibility
        if 'manufacturing_assessment' in analysis_results:
            assessment = analysis_results['manufacturing_assessment']
            report += "## Manufacturing Feasibility Assessment\n\n"
            
            report += f"**Overall Feasibility:** {assessment.get('overall_feasibility', 'Unknown').title()}\n\n"
            
            if 'critical_factors' in assessment:
                report += "**Critical Success Factors:**\n"
                for factor in assessment['critical_factors']:
                    report += f"- {factor}\n"
                report += "\n"
            
            if 'recommendations' in assessment:
                report += "**Manufacturing Recommendations:**\n"
                for rec in assessment['recommendations']:
                    report += f"- {rec}\n"
                report += "\n"
        
        # Implementation Timeline
        report += "## Implementation Timeline\n\n"
        report += "### Phase 1: Platform Selection and Optimization (Months 1-6)\n"
        report += "- Finalize platform selection based on recommendations\n"
        report += "- Optimize vaccine design for selected platform\n"
        report += "- Conduct detailed feasibility studies\n"
        report += "- Establish manufacturing partnerships\n\n"
        
        report += "### Phase 2: Preclinical Development (Months 7-18)\n"
        report += "- Manufacturing process development\n"
        report += "- Analytical method development\n"
        report += "- Preclinical safety and efficacy studies\n"
        report += "- Regulatory preparation\n\n"
        
        report += "### Phase 3: Clinical Development (Months 19-54)\n"
        report += "- Phase I safety studies\n"
        report += "- Phase II proof-of-concept studies\n"
        report += "- Phase III efficacy trials\n"
        report += "- Regulatory submission preparation\n\n"
        
        report += "### Phase 4: Regulatory and Commercial (Months 55+)\n"
        report += "- Regulatory review and approval\n"
        report += "- Commercial manufacturing scale-up\n"
        report += "- Market launch preparation\n"
        report += "- Post-market surveillance\n"
        
        # Save report
        with open(output_path, 'w') as f:
            f.write(report)
        
        logger.info(f"Delivery system integration report saved to: {output_path}")
        return output_path

# Testing and example usage
def test_vaccine_delivery_system_integration():
    """Test the vaccine delivery system integration."""
    logger.info("Testing vaccine delivery system integration")
    
    # Sample epitopes
    sample_epitopes = [
        {
            'id': 'epitope_1',
            'sequence': 'YLQPRTFLL',
            'epitope_type': 'T-cell',
            'prediction_score': 0.85
        },
        {
            'id': 'epitope_2',
            'sequence': 'KIADYNYKL',
            'epitope_type': 'T-cell',
            'prediction_score': 0.78
        },
        {
            'id': 'epitope_3',
            'sequence': 'SLYNTVATL',
            'epitope_type': 'T-cell',
            'prediction_score': 0.72
        },
        {
            'id': 'epitope_4',
            'sequence': 'GILGFVFTL',
            'epitope_type': 'T-cell',
            'prediction_score': 0.68
        },
        {
            'id': 'epitope_5',
            'sequence': 'FIAGLIAIV',
            'epitope_type': 'T-cell',
            'prediction_score': 0.65
        }
    ]
    
    # Initialize integrator
    integrator = VaccineDeliverySystemIntegrator()
    
    # Design platform-specific vaccines
    delivery_analysis = integrator.design_platform_specific_vaccines(
        sample_epitopes,
        target_platforms=['mRNA', 'viral_vector', 'peptide', 'protein_subunit'],
        disease_context='cancer'
    )
    
    # Generate report
    output_dir = Path("/home/ubuntu/vaxgenai/results/delivery_systems")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    report_path = integrator.generate_delivery_system_report(
        delivery_analysis,
        str(output_dir / "delivery_system_integration_report.md")
    )
    
    # Save analysis data
    with open(output_dir / "delivery_analysis.json", 'w') as f:
        json.dump(delivery_analysis, f, indent=2, default=str)
    
    results = {
        'platforms_analyzed': len(delivery_analysis['vaccine_designs']),
        'top_recommendation': delivery_analysis['platform_comparison'].get('top_recommendation'),
        'best_efficacy_platform': delivery_analysis['platform_comparison']['comparison_summary']['best_efficacy']['platform'],
        'lowest_cost_platform': delivery_analysis['platform_comparison']['comparison_summary']['lowest_cost']['platform'],
        'report_path': report_path
    }
    
    logger.info("Vaccine delivery system integration test completed")
    return results

if __name__ == "__main__":
    # Run test
    test_results = test_vaccine_delivery_system_integration()
    print(f"Delivery system integration completed")
    print(f"Platforms analyzed: {test_results['platforms_analyzed']}")
    print(f"Top recommendation: {test_results['top_recommendation']}")
    print(f"Best efficacy: {test_results['best_efficacy_platform']}")
    print(f"Lowest cost: {test_results['lowest_cost_platform']}")

