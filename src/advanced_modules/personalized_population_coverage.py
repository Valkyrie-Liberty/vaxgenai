"""
Personalized Population Coverage Module for VaxGenAI

This module implements personalized population coverage analysis to ensure
vaccine efficacy across diverse populations and address health equity concerns.

Key Features:
- Patient-specific HLA typing and coverage analysis
- Population-specific epitope optimization
- Health equity assessment and recommendations
- Personalized vaccine design strategies
- Global coverage optimization

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

class PersonalizedPopulationCoverageAnalyzer:
    """
    Personalized population coverage analysis system.
    
    This class analyzes epitope coverage for specific populations and individuals,
    optimizes vaccine designs for global coverage, and addresses health equity.
    """
    
    def __init__(self):
        """Initialize the personalized population coverage analyzer."""
        
        # Global HLA frequency data (simplified - in real implementation, use comprehensive databases)
        self.hla_frequencies = {
            'global': {
                'A*01:01': 0.142, 'A*02:01': 0.285, 'A*03:01': 0.128, 'A*11:01': 0.095,
                'A*24:02': 0.089, 'A*26:01': 0.067, 'A*68:01': 0.045, 'A*33:01': 0.034,
                'B*07:02': 0.089, 'B*08:01': 0.067, 'B*15:01': 0.045, 'B*35:01': 0.078,
                'B*40:01': 0.056, 'B*44:02': 0.067, 'B*51:01': 0.045, 'B*57:01': 0.034
            },
            'european': {
                'A*01:01': 0.168, 'A*02:01': 0.312, 'A*03:01': 0.145, 'A*11:01': 0.078,
                'A*24:02': 0.067, 'A*26:01': 0.089, 'A*68:01': 0.056, 'A*33:01': 0.023,
                'B*07:02': 0.112, 'B*08:01': 0.089, 'B*15:01': 0.056, 'B*35:01': 0.067,
                'B*40:01': 0.078, 'B*44:02': 0.089, 'B*51:01': 0.045, 'B*57:01': 0.034
            },
            'east_asian': {
                'A*01:01': 0.089, 'A*02:01': 0.234, 'A*03:01': 0.067, 'A*11:01': 0.156,
                'A*24:02': 0.178, 'A*26:01': 0.034, 'A*68:01': 0.023, 'A*33:01': 0.067,
                'B*07:02': 0.045, 'B*08:01': 0.023, 'B*15:01': 0.034, 'B*35:01': 0.089,
                'B*40:01': 0.034, 'B*44:02': 0.045, 'B*51:01': 0.067, 'B*57:01': 0.023
            },
            'african': {
                'A*01:01': 0.123, 'A*02:01': 0.198, 'A*03:01': 0.089, 'A*11:01': 0.134,
                'A*24:02': 0.067, 'A*26:01': 0.045, 'A*68:01': 0.078, 'A*33:01': 0.056,
                'B*07:02': 0.067, 'B*08:01': 0.045, 'B*15:01': 0.089, 'B*35:01': 0.123,
                'B*40:01': 0.045, 'B*44:02': 0.034, 'B*51:01': 0.078, 'B*57:01': 0.067
            },
            'south_asian': {
                'A*01:01': 0.134, 'A*02:01': 0.267, 'A*03:01': 0.112, 'A*11:01': 0.089,
                'A*24:02': 0.078, 'A*26:01': 0.056, 'A*68:01': 0.034, 'A*33:01': 0.045,
                'B*07:02': 0.078, 'B*08:01': 0.056, 'B*15:01': 0.067, 'B*35:01': 0.089,
                'B*40:01': 0.067, 'B*44:02': 0.056, 'B*51:01': 0.034, 'B*57:01': 0.045
            },
            'hispanic': {
                'A*01:01': 0.156, 'A*02:01': 0.289, 'A*03:01': 0.123, 'A*11:01': 0.067,
                'A*24:02': 0.089, 'A*26:01': 0.078, 'A*68:01': 0.045, 'A*33:01': 0.034,
                'B*07:02': 0.089, 'B*08:01': 0.067, 'B*15:01': 0.045, 'B*35:01': 0.078,
                'B*40:01': 0.056, 'B*44:02': 0.078, 'B*51:01': 0.045, 'B*57:01': 0.034
            }
        }
        
        # MHC Class II frequencies (simplified)
        self.hla_class_ii_frequencies = {
            'global': {
                'DRB1*01:01': 0.089, 'DRB1*03:01': 0.123, 'DRB1*04:01': 0.078,
                'DRB1*07:01': 0.134, 'DRB1*11:01': 0.067, 'DRB1*13:01': 0.089,
                'DRB1*15:01': 0.112, 'DQB1*02:01': 0.156, 'DQB1*03:01': 0.089,
                'DQB1*05:01': 0.078, 'DQB1*06:02': 0.067
            }
        }
        
        # Population demographics for health equity analysis
        self.population_demographics = {
            'global': {'population_size': 7800000000, 'healthcare_access': 0.65},
            'european': {'population_size': 750000000, 'healthcare_access': 0.85},
            'east_asian': {'population_size': 1650000000, 'healthcare_access': 0.72},
            'african': {'population_size': 1340000000, 'healthcare_access': 0.45},
            'south_asian': {'population_size': 1980000000, 'healthcare_access': 0.58},
            'hispanic': {'population_size': 650000000, 'healthcare_access': 0.62}
        }
        
        # Disease prevalence by population (simplified)
        self.disease_prevalence = {
            'leukemia': {
                'global': 0.000045, 'european': 0.000052, 'east_asian': 0.000038,
                'african': 0.000034, 'south_asian': 0.000041, 'hispanic': 0.000047
            },
            'pancreatic_cancer': {
                'global': 0.000089, 'european': 0.000098, 'east_asian': 0.000067,
                'african': 0.000078, 'south_asian': 0.000082, 'hispanic': 0.000091
            },
            'hiv': {
                'global': 0.0095, 'european': 0.0023, 'east_asian': 0.0012,
                'african': 0.0387, 'south_asian': 0.0034, 'hispanic': 0.0067
            }
        }
    
    def analyze_patient_specific_coverage(self, patient_hla_profile: Dict,
                                        epitopes: List[Dict],
                                        population: str = 'global') -> Dict:
        """
        Analyze epitope coverage for a specific patient's HLA profile.
        
        Args:
            patient_hla_profile: Patient's HLA typing results
            epitopes: List of predicted epitopes
            population: Population group for comparison
            
        Returns:
            Patient-specific coverage analysis
        """
        logger.info(f"Analyzing patient-specific coverage for {population} population")
        
        # Validate HLA profile
        validated_hla = self._validate_hla_profile(patient_hla_profile)
        
        # Calculate epitope binding for patient's HLA alleles
        patient_coverage = self._calculate_patient_epitope_binding(validated_hla, epitopes)
        
        # Compare with population average
        population_comparison = self._compare_with_population(
            patient_coverage, epitopes, population
        )
        
        # Identify coverage gaps
        coverage_gaps = self._identify_coverage_gaps(patient_coverage, epitopes)
        
        # Generate personalized recommendations
        recommendations = self._generate_personalized_recommendations(
            patient_coverage, coverage_gaps, validated_hla
        )
        
        analysis_result = {
            'patient_id': patient_hla_profile.get('patient_id', 'unknown'),
            'hla_profile': validated_hla,
            'population': population,
            'epitope_coverage': patient_coverage,
            'population_comparison': population_comparison,
            'coverage_gaps': coverage_gaps,
            'personalized_recommendations': recommendations,
            'coverage_summary': self._summarize_coverage(patient_coverage)
        }
        
        return analysis_result
    
    def _validate_hla_profile(self, hla_profile: Dict) -> Dict:
        """Validate and standardize HLA profile format."""
        validated = {
            'class_i': {},
            'class_ii': {},
            'typing_method': hla_profile.get('typing_method', 'unknown'),
            'resolution': hla_profile.get('resolution', '2-digit')
        }
        
        # Process Class I alleles
        for locus in ['A', 'B', 'C']:
            alleles = hla_profile.get(f'HLA-{locus}', [])
            if isinstance(alleles, str):
                alleles = [alleles]
            validated['class_i'][locus] = alleles[:2]  # Maximum 2 alleles per locus
        
        # Process Class II alleles
        for locus in ['DRB1', 'DQB1', 'DPB1']:
            alleles = hla_profile.get(f'HLA-{locus}', [])
            if isinstance(alleles, str):
                alleles = [alleles]
            validated['class_ii'][locus] = alleles[:2]
        
        return validated
    
    def _calculate_patient_epitope_binding(self, hla_profile: Dict, 
                                         epitopes: List[Dict]) -> Dict:
        """Calculate epitope binding for patient's specific HLA alleles."""
        coverage = {
            'total_epitopes': len(epitopes),
            'covered_epitopes': 0,
            'class_i_coverage': {},
            'class_ii_coverage': {},
            'epitope_details': []
        }
        
        for epitope in epitopes:
            epitope_covered = False
            binding_alleles = []
            
            # Check Class I binding
            if epitope.get('epitope_type') in ['T-cell', 'CD8']:
                for locus, alleles in hla_profile['class_i'].items():
                    for allele in alleles:
                        if allele and self._predicts_binding(epitope, allele, 'class_i'):
                            epitope_covered = True
                            binding_alleles.append(allele)
            
            # Check Class II binding
            elif epitope.get('epitope_type') in ['CD4', 'helper']:
                for locus, alleles in hla_profile['class_ii'].items():
                    for allele in alleles:
                        if allele and self._predicts_binding(epitope, allele, 'class_ii'):
                            epitope_covered = True
                            binding_alleles.append(allele)
            
            if epitope_covered:
                coverage['covered_epitopes'] += 1
            
            coverage['epitope_details'].append({
                'epitope_id': epitope.get('id', f"epitope_{len(coverage['epitope_details'])}"),
                'sequence': epitope.get('sequence', epitope.get('peptide', '')),
                'covered': epitope_covered,
                'binding_alleles': binding_alleles,
                'prediction_score': epitope.get('prediction_score', 0.5)
            })
        
        # Calculate coverage percentage
        coverage['coverage_percentage'] = (
            coverage['covered_epitopes'] / coverage['total_epitopes'] * 100
            if coverage['total_epitopes'] > 0 else 0
        )
        
        return coverage
    
    def _predicts_binding(self, epitope: Dict, hla_allele: str, hla_class: str) -> bool:
        """Predict if epitope binds to specific HLA allele."""
        # Simplified binding prediction - in real implementation, use actual prediction models
        
        epitope_seq = epitope.get('sequence', epitope.get('peptide', ''))
        if not epitope_seq:
            return False
        
        # Length constraints
        if hla_class == 'class_i' and len(epitope_seq) not in range(8, 12):
            return False
        elif hla_class == 'class_ii' and len(epitope_seq) not in range(12, 26):
            return False
        
        # Simplified binding prediction based on sequence properties
        binding_score = 0.0
        
        # Anchor residue preferences (simplified)
        if hla_class == 'class_i':
            if 'A*02:01' in hla_allele:
                # A*02:01 preferences: L, M at position 2; V, L, I at C-terminus
                if len(epitope_seq) >= 9:
                    if epitope_seq[1] in 'LM':
                        binding_score += 0.3
                    if epitope_seq[-1] in 'VLI':
                        binding_score += 0.3
            elif 'A*01:01' in hla_allele:
                # A*01:01 preferences: T, S at position 2; Y at C-terminus
                if len(epitope_seq) >= 9:
                    if epitope_seq[1] in 'TS':
                        binding_score += 0.3
                    if epitope_seq[-1] in 'Y':
                        binding_score += 0.3
        
        # Add base prediction score
        base_score = epitope.get('prediction_score', 0.5)
        binding_score += 0.4 * base_score
        
        # Add some randomness to simulate prediction uncertainty
        binding_score += np.random.normal(0, 0.1)
        
        return binding_score > 0.5
    
    def _compare_with_population(self, patient_coverage: Dict, 
                               epitopes: List[Dict], population: str) -> Dict:
        """Compare patient coverage with population average."""
        
        # Calculate population average coverage
        pop_frequencies = self.hla_frequencies.get(population, self.hla_frequencies['global'])
        
        # Simulate population coverage calculation
        population_coverage = 0.0
        for epitope in epitopes:
            epitope_pop_coverage = 0.0
            
            # Calculate coverage across population HLA frequencies
            for allele, frequency in pop_frequencies.items():
                if self._predicts_binding(epitope, allele, 'class_i'):
                    epitope_pop_coverage += frequency
            
            # Avoid double counting (simplified)
            epitope_pop_coverage = min(1.0, epitope_pop_coverage)
            population_coverage += epitope_pop_coverage
        
        avg_population_coverage = (
            population_coverage / len(epitopes) * 100 
            if epitopes else 0
        )
        
        comparison = {
            'patient_coverage': patient_coverage['coverage_percentage'],
            'population_average': avg_population_coverage,
            'relative_coverage': (
                patient_coverage['coverage_percentage'] / avg_population_coverage
                if avg_population_coverage > 0 else 1.0
            ),
            'coverage_category': self._categorize_coverage(
                patient_coverage['coverage_percentage'], avg_population_coverage
            )
        }
        
        return comparison
    
    def _categorize_coverage(self, patient_coverage: float, population_avg: float) -> str:
        """Categorize patient coverage relative to population."""
        relative_coverage = patient_coverage / population_avg if population_avg > 0 else 1.0
        
        if relative_coverage >= 1.2:
            return 'above_average'
        elif relative_coverage >= 0.8:
            return 'average'
        elif relative_coverage >= 0.6:
            return 'below_average'
        else:
            return 'poor'
    
    def _identify_coverage_gaps(self, patient_coverage: Dict, 
                              epitopes: List[Dict]) -> List[Dict]:
        """Identify epitopes not covered by patient's HLA profile."""
        gaps = []
        
        for detail in patient_coverage['epitope_details']:
            if not detail['covered']:
                gap_info = {
                    'epitope_id': detail['epitope_id'],
                    'sequence': detail['sequence'],
                    'prediction_score': detail['prediction_score'],
                    'gap_severity': self._assess_gap_severity(detail),
                    'alternative_alleles': self._suggest_alternative_alleles(detail)
                }
                gaps.append(gap_info)
        
        # Sort by gap severity
        gaps.sort(key=lambda x: x['gap_severity'], reverse=True)
        
        return gaps
    
    def _assess_gap_severity(self, epitope_detail: Dict) -> float:
        """Assess the severity of a coverage gap."""
        # Higher severity for high-scoring epitopes that are not covered
        base_severity = epitope_detail['prediction_score']
        
        # Increase severity for epitopes with no binding alleles
        if not epitope_detail['binding_alleles']:
            base_severity *= 1.5
        
        return min(1.0, base_severity)
    
    def _suggest_alternative_alleles(self, epitope_detail: Dict) -> List[str]:
        """Suggest HLA alleles that might bind this epitope."""
        # Simplified suggestion based on common alleles
        common_alleles = ['A*02:01', 'A*01:01', 'A*03:01', 'B*07:02', 'B*08:01']
        
        # In real implementation, use actual binding prediction
        suggested = []
        for allele in common_alleles:
            # Simulate binding prediction
            if np.random.random() > 0.7:  # 30% chance of binding
                suggested.append(allele)
        
        return suggested[:3]  # Return top 3 suggestions
    
    def _generate_personalized_recommendations(self, patient_coverage: Dict,
                                             coverage_gaps: List[Dict],
                                             hla_profile: Dict) -> List[Dict]:
        """Generate personalized vaccine recommendations."""
        recommendations = []
        
        # Coverage-based recommendations
        coverage_pct = patient_coverage['coverage_percentage']
        
        if coverage_pct >= 80:
            recommendations.append({
                'type': 'vaccine_strategy',
                'priority': 'high',
                'recommendation': 'Standard vaccine approach suitable',
                'rationale': f'Excellent epitope coverage ({coverage_pct:.1f}%)'
            })
        elif coverage_pct >= 60:
            recommendations.append({
                'type': 'vaccine_strategy',
                'priority': 'medium',
                'recommendation': 'Consider epitope optimization',
                'rationale': f'Good coverage ({coverage_pct:.1f}%) but room for improvement'
            })
        else:
            recommendations.append({
                'type': 'vaccine_strategy',
                'priority': 'high',
                'recommendation': 'Personalized vaccine design required',
                'rationale': f'Poor coverage ({coverage_pct:.1f}%) requires customization'
            })
        
        # Gap-specific recommendations
        if len(coverage_gaps) > 5:
            recommendations.append({
                'type': 'epitope_selection',
                'priority': 'high',
                'recommendation': 'Include additional epitopes targeting uncovered regions',
                'rationale': f'{len(coverage_gaps)} significant coverage gaps identified'
            })
        
        # HLA-specific recommendations
        rare_alleles = self._identify_rare_alleles(hla_profile)
        if rare_alleles:
            recommendations.append({
                'type': 'hla_consideration',
                'priority': 'medium',
                'recommendation': f'Special consideration for rare HLA alleles: {", ".join(rare_alleles)}',
                'rationale': 'Rare alleles may require specialized epitope selection'
            })
        
        # Delivery system recommendations
        recommendations.extend(self._recommend_delivery_systems(patient_coverage, hla_profile))
        
        return recommendations
    
    def _identify_rare_alleles(self, hla_profile: Dict) -> List[str]:
        """Identify rare HLA alleles in patient profile."""
        rare_alleles = []
        frequency_threshold = 0.05  # Consider alleles <5% frequency as rare
        
        for locus, alleles in hla_profile['class_i'].items():
            for allele in alleles:
                if allele:
                    # Check frequency in global population
                    global_freq = self.hla_frequencies['global'].get(allele, 0)
                    if global_freq < frequency_threshold:
                        rare_alleles.append(allele)
        
        return rare_alleles
    
    def _recommend_delivery_systems(self, patient_coverage: Dict, 
                                  hla_profile: Dict) -> List[Dict]:
        """Recommend optimal vaccine delivery systems."""
        recommendations = []
        
        coverage_pct = patient_coverage['coverage_percentage']
        
        # mRNA vaccine recommendations
        if coverage_pct >= 70:
            recommendations.append({
                'type': 'delivery_system',
                'priority': 'high',
                'recommendation': 'mRNA vaccine platform suitable',
                'rationale': 'Good epitope coverage supports mRNA approach',
                'delivery_details': {
                    'platform': 'mRNA',
                    'advantages': ['Rapid production', 'Strong immune response'],
                    'considerations': ['Cold chain requirements', 'Lipid nanoparticle optimization']
                }
            })
        
        # Peptide vaccine recommendations
        if len(patient_coverage['epitope_details']) <= 10:
            recommendations.append({
                'type': 'delivery_system',
                'priority': 'medium',
                'recommendation': 'Peptide vaccine with adjuvant',
                'rationale': 'Limited epitopes suitable for peptide approach',
                'delivery_details': {
                    'platform': 'peptide',
                    'advantages': ['Precise targeting', 'Safety profile'],
                    'considerations': ['Adjuvant selection', 'Peptide stability']
                }
            })
        
        # Viral vector recommendations
        if coverage_pct < 60:
            recommendations.append({
                'type': 'delivery_system',
                'priority': 'high',
                'recommendation': 'Viral vector for enhanced immunogenicity',
                'rationale': 'Poor coverage requires enhanced immune stimulation',
                'delivery_details': {
                    'platform': 'viral_vector',
                    'advantages': ['Strong immune activation', 'Broad epitope presentation'],
                    'considerations': ['Vector immunity', 'Manufacturing complexity']
                }
            })
        
        return recommendations
    
    def _summarize_coverage(self, patient_coverage: Dict) -> Dict:
        """Summarize patient coverage analysis."""
        return {
            'total_epitopes': patient_coverage['total_epitopes'],
            'covered_epitopes': patient_coverage['covered_epitopes'],
            'coverage_percentage': patient_coverage['coverage_percentage'],
            'coverage_quality': self._assess_coverage_quality(patient_coverage['coverage_percentage']),
            'key_metrics': {
                'high_confidence_covered': len([
                    e for e in patient_coverage['epitope_details'] 
                    if e['covered'] and e['prediction_score'] > 0.7
                ]),
                'total_binding_events': sum([
                    len(e['binding_alleles']) for e in patient_coverage['epitope_details']
                ])
            }
        }
    
    def _assess_coverage_quality(self, coverage_percentage: float) -> str:
        """Assess the quality of epitope coverage."""
        if coverage_percentage >= 80:
            return 'excellent'
        elif coverage_percentage >= 65:
            return 'good'
        elif coverage_percentage >= 50:
            return 'moderate'
        elif coverage_percentage >= 35:
            return 'poor'
        else:
            return 'very_poor'
    
    def optimize_global_coverage(self, epitopes: List[Dict],
                               target_populations: List[str] = None,
                               coverage_threshold: float = 0.7) -> Dict:
        """
        Optimize epitope selection for maximum global population coverage.
        
        Args:
            epitopes: List of candidate epitopes
            target_populations: Specific populations to optimize for
            coverage_threshold: Minimum coverage threshold
            
        Returns:
            Optimized epitope selection and coverage analysis
        """
        logger.info("Optimizing epitope selection for global population coverage")
        
        if target_populations is None:
            target_populations = list(self.hla_frequencies.keys())
        
        # Calculate coverage for each epitope across populations
        epitope_coverage_matrix = self._calculate_epitope_coverage_matrix(
            epitopes, target_populations
        )
        
        # Optimize epitope selection using greedy algorithm
        optimized_selection = self._greedy_epitope_selection(
            epitope_coverage_matrix, coverage_threshold
        )
        
        # Calculate final coverage statistics
        final_coverage = self._calculate_final_coverage(
            optimized_selection, epitope_coverage_matrix, target_populations
        )
        
        # Assess health equity implications
        equity_analysis = self._assess_health_equity(
            final_coverage, target_populations
        )
        
        optimization_result = {
            'selected_epitopes': optimized_selection,
            'population_coverage': final_coverage,
            'health_equity_analysis': equity_analysis,
            'optimization_summary': {
                'total_candidate_epitopes': len(epitopes),
                'selected_epitopes': len(optimized_selection),
                'selection_efficiency': len(optimized_selection) / len(epitopes),
                'average_coverage': np.mean(list(final_coverage.values())),
                'minimum_coverage': min(final_coverage.values()),
                'coverage_variance': np.var(list(final_coverage.values()))
            }
        }
        
        return optimization_result
    
    def _calculate_epitope_coverage_matrix(self, epitopes: List[Dict],
                                         populations: List[str]) -> Dict:
        """Calculate coverage matrix for epitopes across populations."""
        coverage_matrix = {}
        
        for i, epitope in enumerate(epitopes):
            epitope_id = epitope.get('id', f'epitope_{i}')
            coverage_matrix[epitope_id] = {}
            
            for population in populations:
                pop_frequencies = self.hla_frequencies.get(
                    population, self.hla_frequencies['global']
                )
                
                # Calculate population coverage for this epitope
                pop_coverage = 0.0
                covered_individuals = 0.0
                
                # Simplified coverage calculation
                for allele, frequency in pop_frequencies.items():
                    if self._predicts_binding(epitope, allele, 'class_i'):
                        # Assume individuals with this allele are covered
                        covered_individuals += frequency
                
                # Account for heterozygosity (simplified)
                pop_coverage = min(1.0, covered_individuals * 1.8)  # Rough heterozygosity adjustment
                
                coverage_matrix[epitope_id][population] = pop_coverage
        
        return coverage_matrix
    
    def _greedy_epitope_selection(self, coverage_matrix: Dict,
                                coverage_threshold: float) -> List[str]:
        """Select epitopes using greedy algorithm for maximum coverage."""
        selected_epitopes = []
        remaining_epitopes = list(coverage_matrix.keys())
        populations = list(next(iter(coverage_matrix.values())).keys())
        
        # Track current coverage for each population
        current_coverage = {pop: 0.0 for pop in populations}
        
        while remaining_epitopes and min(current_coverage.values()) < coverage_threshold:
            best_epitope = None
            best_improvement = 0.0
            
            # Find epitope that provides best coverage improvement
            for epitope_id in remaining_epitopes:
                # Calculate improvement if this epitope is added
                improvement = 0.0
                for population in populations:
                    epitope_coverage = coverage_matrix[epitope_id][population]
                    # Calculate marginal improvement (simplified)
                    new_coverage = min(1.0, current_coverage[population] + epitope_coverage * 0.5)
                    improvement += max(0, new_coverage - current_coverage[population])
                
                if improvement > best_improvement:
                    best_improvement = improvement
                    best_epitope = epitope_id
            
            if best_epitope:
                selected_epitopes.append(best_epitope)
                remaining_epitopes.remove(best_epitope)
                
                # Update current coverage
                for population in populations:
                    epitope_coverage = coverage_matrix[best_epitope][population]
                    current_coverage[population] = min(
                        1.0, current_coverage[population] + epitope_coverage * 0.5
                    )
            else:
                break  # No more improvement possible
        
        return selected_epitopes
    
    def _calculate_final_coverage(self, selected_epitopes: List[str],
                                coverage_matrix: Dict,
                                populations: List[str]) -> Dict:
        """Calculate final coverage for selected epitopes."""
        final_coverage = {}
        
        for population in populations:
            pop_coverage = 0.0
            
            for epitope_id in selected_epitopes:
                epitope_coverage = coverage_matrix[epitope_id][population]
                # Simplified coverage combination
                pop_coverage = min(1.0, pop_coverage + epitope_coverage * 0.4)
            
            final_coverage[population] = pop_coverage
        
        return final_coverage
    
    def _assess_health_equity(self, coverage_results: Dict,
                            populations: List[str]) -> Dict:
        """Assess health equity implications of coverage results."""
        
        equity_analysis = {
            'coverage_disparity': {},
            'equity_score': 0.0,
            'recommendations': [],
            'vulnerable_populations': []
        }
        
        # Calculate coverage disparities
        coverages = list(coverage_results.values())
        max_coverage = max(coverages)
        min_coverage = min(coverages)
        
        equity_analysis['coverage_disparity'] = {
            'max_coverage': max_coverage,
            'min_coverage': min_coverage,
            'disparity_ratio': max_coverage / min_coverage if min_coverage > 0 else float('inf'),
            'coverage_range': max_coverage - min_coverage
        }
        
        # Calculate equity score (higher is more equitable)
        coverage_variance = np.var(coverages)
        equity_analysis['equity_score'] = max(0, 1 - coverage_variance)
        
        # Identify vulnerable populations
        threshold = np.mean(coverages) * 0.8  # 80% of average coverage
        for population, coverage in coverage_results.items():
            if coverage < threshold:
                # Consider disease burden and healthcare access
                disease_burden = sum(self.disease_prevalence[disease].get(population, 0) 
                                   for disease in self.disease_prevalence.keys())
                healthcare_access = self.population_demographics.get(population, {}).get('healthcare_access', 0.5)
                
                vulnerability_score = disease_burden / healthcare_access if healthcare_access > 0 else float('inf')
                
                equity_analysis['vulnerable_populations'].append({
                    'population': population,
                    'coverage': coverage,
                    'vulnerability_score': vulnerability_score,
                    'disease_burden': disease_burden,
                    'healthcare_access': healthcare_access
                })
        
        # Generate equity recommendations
        if equity_analysis['coverage_disparity']['disparity_ratio'] > 1.5:
            equity_analysis['recommendations'].append({
                'type': 'disparity_reduction',
                'priority': 'high',
                'recommendation': 'Address significant coverage disparities between populations',
                'specific_actions': [
                    'Include population-specific epitopes',
                    'Consider alternative vaccine platforms',
                    'Develop targeted delivery strategies'
                ]
            })
        
        if equity_analysis['vulnerable_populations']:
            equity_analysis['recommendations'].append({
                'type': 'vulnerable_population_support',
                'priority': 'high',
                'recommendation': 'Provide additional support for vulnerable populations',
                'specific_actions': [
                    'Prioritize vaccine access',
                    'Develop culturally appropriate delivery methods',
                    'Consider subsidized or free vaccine programs'
                ]
            })
        
        return equity_analysis
    
    def generate_population_coverage_report(self, analysis_results: Dict,
                                          output_path: str) -> str:
        """
        Generate comprehensive population coverage report.
        
        Args:
            analysis_results: Results from coverage analysis
            output_path: Output file path
            
        Returns:
            Path to generated report
        """
        logger.info("Generating population coverage report")
        
        report = "# Personalized Population Coverage Analysis Report\n\n"
        
        # Executive Summary
        report += "## Executive Summary\n\n"
        
        if 'optimization_summary' in analysis_results:
            summary = analysis_results['optimization_summary']
            report += f"**Total candidate epitopes analyzed:** {summary['total_candidate_epitopes']}\n"
            report += f"**Epitopes selected for vaccine:** {summary['selected_epitopes']}\n"
            report += f"**Selection efficiency:** {summary['selection_efficiency']:.1%}\n"
            report += f"**Average population coverage:** {summary['average_coverage']:.1%}\n"
            report += f"**Minimum population coverage:** {summary['minimum_coverage']:.1%}\n\n"
        
        # Population Coverage Results
        if 'population_coverage' in analysis_results:
            report += "## Population Coverage Results\n\n"
            report += "| Population | Coverage | Quality Assessment |\n"
            report += "|------------|----------|-------------------|\n"
            
            for population, coverage in analysis_results['population_coverage'].items():
                quality = self._assess_coverage_quality(coverage * 100)
                report += f"| {population.replace('_', ' ').title()} | {coverage:.1%} | {quality.title()} |\n"
            report += "\n"
        
        # Health Equity Analysis
        if 'health_equity_analysis' in analysis_results:
            equity = analysis_results['health_equity_analysis']
            report += "## Health Equity Analysis\n\n"
            
            disparity = equity['coverage_disparity']
            report += f"**Coverage Disparity Ratio:** {disparity['disparity_ratio']:.2f}\n"
            report += f"**Coverage Range:** {disparity['coverage_range']:.1%}\n"
            report += f"**Equity Score:** {equity['equity_score']:.3f} (higher is more equitable)\n\n"
            
            # Vulnerable populations
            if equity['vulnerable_populations']:
                report += "### Vulnerable Populations\n\n"
                report += "| Population | Coverage | Vulnerability Score | Disease Burden | Healthcare Access |\n"
                report += "|------------|----------|-------------------|----------------|------------------|\n"
                
                for pop_data in equity['vulnerable_populations']:
                    report += f"| {pop_data['population'].replace('_', ' ').title()} | "
                    report += f"{pop_data['coverage']:.1%} | {pop_data['vulnerability_score']:.3f} | "
                    report += f"{pop_data['disease_burden']:.6f} | {pop_data['healthcare_access']:.1%} |\n"
                report += "\n"
        
        # Selected Epitopes
        if 'selected_epitopes' in analysis_results:
            report += "## Selected Epitopes for Global Coverage\n\n"
            selected = analysis_results['selected_epitopes']
            report += f"Total selected epitopes: {len(selected)}\n\n"
            
            for i, epitope_id in enumerate(selected[:15], 1):  # Show top 15
                report += f"{i}. {epitope_id}\n"
            
            if len(selected) > 15:
                report += f"... and {len(selected) - 15} more epitopes\n"
            report += "\n"
        
        # Recommendations
        if 'health_equity_analysis' in analysis_results and 'recommendations' in analysis_results['health_equity_analysis']:
            report += "## Recommendations\n\n"
            
            for i, rec in enumerate(analysis_results['health_equity_analysis']['recommendations'], 1):
                report += f"### {i}. {rec['recommendation']}\n\n"
                report += f"**Type:** {rec['type'].replace('_', ' ').title()}\n"
                report += f"**Priority:** {rec['priority'].title()}\n"
                
                if 'specific_actions' in rec:
                    report += "**Specific Actions:**\n"
                    for action in rec['specific_actions']:
                        report += f"- {action}\n"
                report += "\n"
        
        # Implementation Considerations
        report += "## Implementation Considerations\n\n"
        report += "### Vaccine Development Strategy\n"
        report += "1. **Multi-epitope Design**: Include epitopes optimized for global coverage\n"
        report += "2. **Platform Selection**: Consider mRNA, viral vector, or peptide platforms based on coverage needs\n"
        report += "3. **Adjuvant Strategy**: Select adjuvants that enhance responses across diverse populations\n"
        report += "4. **Manufacturing**: Plan for global distribution and cold chain requirements\n\n"
        
        report += "### Clinical Development\n"
        report += "1. **Phase I/II Trials**: Include diverse populations in early trials\n"
        report += "2. **Immunogenicity Assessment**: Monitor responses across different HLA backgrounds\n"
        report += "3. **Efficacy Endpoints**: Design endpoints that capture population-specific responses\n"
        report += "4. **Regulatory Strategy**: Engage with global regulatory authorities early\n\n"
        
        report += "### Access and Equity\n"
        report += "1. **Pricing Strategy**: Consider tiered pricing for different economic regions\n"
        report += "2. **Distribution**: Prioritize vulnerable populations and high-burden areas\n"
        report += "3. **Local Partnerships**: Collaborate with local healthcare systems\n"
        report += "4. **Technology Transfer**: Consider local manufacturing capabilities\n"
        
        # Save report
        with open(output_path, 'w') as f:
            f.write(report)
        
        logger.info(f"Population coverage report saved to: {output_path}")
        return output_path

# Testing and example usage
def test_personalized_population_coverage():
    """Test the personalized population coverage analysis system."""
    logger.info("Testing personalized population coverage analysis system")
    
    # Sample patient HLA profile
    patient_hla = {
        'patient_id': 'P001',
        'HLA-A': ['A*02:01', 'A*01:01'],
        'HLA-B': ['B*07:02', 'B*08:01'],
        'HLA-C': ['C*07:01', 'C*08:02'],
        'HLA-DRB1': ['DRB1*03:01', 'DRB1*07:01'],
        'HLA-DQB1': ['DQB1*02:01', 'DQB1*03:01'],
        'typing_method': 'NGS',
        'resolution': '4-digit'
    }
    
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
    
    # Initialize analyzer
    analyzer = PersonalizedPopulationCoverageAnalyzer()
    
    # Analyze patient-specific coverage
    patient_analysis = analyzer.analyze_patient_specific_coverage(
        patient_hla, sample_epitopes, 'european'
    )
    
    # Optimize for global coverage
    global_optimization = analyzer.optimize_global_coverage(
        sample_epitopes, 
        target_populations=['global', 'european', 'east_asian', 'african'],
        coverage_threshold=0.6
    )
    
    # Generate reports
    output_dir = Path("/home/ubuntu/vaxgenai/results/population_coverage")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Patient-specific report
    patient_report_path = analyzer.generate_population_coverage_report(
        {'patient_analysis': patient_analysis},
        str(output_dir / "patient_coverage_report.md")
    )
    
    # Global optimization report
    global_report_path = analyzer.generate_population_coverage_report(
        global_optimization,
        str(output_dir / "global_coverage_optimization_report.md")
    )
    
    # Save analysis data
    with open(output_dir / "patient_analysis.json", 'w') as f:
        json.dump(patient_analysis, f, indent=2, default=str)
    
    with open(output_dir / "global_optimization.json", 'w') as f:
        json.dump(global_optimization, f, indent=2, default=str)
    
    results = {
        'patient_coverage': patient_analysis['coverage_summary']['coverage_percentage'],
        'patient_coverage_quality': patient_analysis['coverage_summary']['coverage_quality'],
        'global_optimization_epitopes': len(global_optimization['selected_epitopes']),
        'average_global_coverage': global_optimization['optimization_summary']['average_coverage'],
        'equity_score': global_optimization['health_equity_analysis']['equity_score'],
        'patient_report_path': patient_report_path,
        'global_report_path': global_report_path
    }
    
    logger.info("Personalized population coverage analysis test completed")
    return results

if __name__ == "__main__":
    # Run test
    test_results = test_personalized_population_coverage()
    print(f"Personalized population coverage analysis completed")
    print(f"Patient coverage: {test_results['patient_coverage']:.1f}% ({test_results['patient_coverage_quality']})")
    print(f"Global optimization selected: {test_results['global_optimization_epitopes']} epitopes")
    print(f"Average global coverage: {test_results['average_global_coverage']:.1%}")
    print(f"Health equity score: {test_results['equity_score']:.3f}")

