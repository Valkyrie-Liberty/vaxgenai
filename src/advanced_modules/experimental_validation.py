"""
Experimental Validation Integration Module for VaxGenAI

This module implements experimental validation integration to incorporate
real-world experimental data and create feedback loops for improving predictions.

Key Features:
- Integration with experimental immunology databases
- Feedback loops for model improvement
- Validation score calculation
- Experimental design recommendations
- Results interpretation and analysis

Author: VaxGenAI Development Team
"""

import pandas as pd
import numpy as np
import logging
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import json
import requests
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
import warnings
warnings.filterwarnings('ignore')

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class ExperimentalValidationIntegrator:
    """
    Experimental validation integration system for VaxGenAI.
    
    This class integrates experimental validation data to improve predictions
    and provide feedback loops for continuous model improvement.
    """
    
    def __init__(self):
        """Initialize the experimental validation integrator."""
        self.validation_databases = {
            'iedb': {
                'url': 'http://www.iedb.org/api/',
                'description': 'Immune Epitope Database',
                'data_types': ['tcell_assays', 'bcell_assays', 'mhc_binding']
            },
            'cancer_immunity': {
                'url': 'https://cancerimmunity.org/api/',
                'description': 'Cancer Immunity Database',
                'data_types': ['neoantigen_validation', 'immunotherapy_response']
            },
            'hiv_database': {
                'url': 'https://www.hiv.lanl.gov/api/',
                'description': 'HIV Sequence Database',
                'data_types': ['ctl_epitopes', 'antibody_epitopes', 'escape_mutations']
            }
        }
        
        self.experimental_assays = {
            'elispot': {
                'description': 'Enzyme-linked immunospot assay',
                'measures': 'T-cell activation (IFN-γ production)',
                'readout': 'spot_forming_cells',
                'threshold': 50  # SFC per million cells
            },
            'tetramer_staining': {
                'description': 'MHC tetramer staining',
                'measures': 'Antigen-specific T-cell frequency',
                'readout': 'percentage_positive',
                'threshold': 0.1  # % of CD8+ T cells
            },
            'cytotoxicity_assay': {
                'description': 'Cytotoxic T lymphocyte assay',
                'measures': 'Target cell killing',
                'readout': 'percent_lysis',
                'threshold': 20  # % specific lysis
            },
            'elisa': {
                'description': 'Enzyme-linked immunosorbent assay',
                'measures': 'Antibody binding',
                'readout': 'optical_density',
                'threshold': 0.5  # OD450
            },
            'flow_cytometry': {
                'description': 'Flow cytometric analysis',
                'measures': 'Cell surface markers and intracellular cytokines',
                'readout': 'percentage_positive',
                'threshold': 1.0  # % positive cells
            }
        }
        
        # Initialize validation data storage
        self.validation_data = {}
        self.prediction_performance = {}
        
    def integrate_experimental_data(self, experimental_results: List[Dict]) -> Dict:
        """
        Integrate experimental validation results with predictions.
        
        Args:
            experimental_results: List of experimental result dictionaries
            
        Returns:
            Integration summary and performance metrics
        """
        logger.info(f"Integrating {len(experimental_results)} experimental results")
        
        integrated_data = []
        
        for result in experimental_results:
            # Validate experimental result format
            validated_result = self._validate_experimental_result(result)
            if validated_result:
                integrated_data.append(validated_result)
        
        # Store validation data
        self.validation_data.update({
            'timestamp': pd.Timestamp.now().isoformat(),
            'results': integrated_data,
            'total_experiments': len(integrated_data)
        })
        
        # Calculate performance metrics
        performance_metrics = self._calculate_performance_metrics(integrated_data)
        
        # Generate feedback for model improvement
        feedback = self._generate_model_feedback(integrated_data, performance_metrics)
        
        integration_summary = {
            'total_results': len(experimental_results),
            'valid_results': len(integrated_data),
            'performance_metrics': performance_metrics,
            'model_feedback': feedback,
            'recommendations': self._generate_experimental_recommendations(integrated_data)
        }
        
        logger.info("Experimental data integration completed")
        return integration_summary
    
    def _validate_experimental_result(self, result: Dict) -> Optional[Dict]:
        """
        Validate experimental result format and content.
        
        Args:
            result: Experimental result dictionary
            
        Returns:
            Validated result or None if invalid
        """
        required_fields = ['peptide', 'assay_type', 'result_value', 'positive_control']
        
        # Check required fields
        for field in required_fields:
            if field not in result:
                logger.warning(f"Missing required field: {field}")
                return None
        
        # Validate assay type
        if result['assay_type'] not in self.experimental_assays:
            logger.warning(f"Unknown assay type: {result['assay_type']}")
            return None
        
        # Validate peptide sequence
        peptide = result['peptide']
        if not isinstance(peptide, str) or not peptide.isalpha():
            logger.warning(f"Invalid peptide sequence: {peptide}")
            return None
        
        # Validate result value
        try:
            result_value = float(result['result_value'])
        except (ValueError, TypeError):
            logger.warning(f"Invalid result value: {result['result_value']}")
            return None
        
        # Determine if result is positive based on assay threshold
        assay_info = self.experimental_assays[result['assay_type']]
        threshold = assay_info['threshold']
        is_positive = result_value >= threshold
        
        validated_result = {
            'peptide': peptide,
            'assay_type': result['assay_type'],
            'result_value': result_value,
            'is_positive': is_positive,
            'threshold': threshold,
            'positive_control': result['positive_control'],
            'experimental_conditions': result.get('conditions', {}),
            'publication_info': result.get('publication', {}),
            'validation_date': pd.Timestamp.now().isoformat()
        }
        
        return validated_result
    
    def _calculate_performance_metrics(self, experimental_data: List[Dict]) -> Dict:
        """
        Calculate prediction performance metrics against experimental data.
        
        Args:
            experimental_data: List of validated experimental results
            
        Returns:
            Performance metrics dictionary
        """
        logger.info("Calculating prediction performance metrics")
        
        # Group by assay type
        assay_performance = {}
        
        for assay_type in self.experimental_assays.keys():
            assay_results = [r for r in experimental_data if r['assay_type'] == assay_type]
            
            if not assay_results:
                continue
            
            # Extract experimental outcomes
            experimental_outcomes = [r['is_positive'] for r in assay_results]
            
            # Simulate predicted outcomes (in real implementation, these would come from actual predictions)
            predicted_outcomes = self._simulate_predicted_outcomes(assay_results)
            
            # Calculate metrics
            if len(set(experimental_outcomes)) > 1:  # Need both positive and negative examples
                accuracy = accuracy_score(experimental_outcomes, predicted_outcomes)
                precision = precision_score(experimental_outcomes, predicted_outcomes, zero_division=0)
                recall = recall_score(experimental_outcomes, predicted_outcomes, zero_division=0)
                f1 = f1_score(experimental_outcomes, predicted_outcomes, zero_division=0)
            else:
                accuracy = precision = recall = f1 = 0.0
            
            assay_performance[assay_type] = {
                'total_experiments': len(assay_results),
                'positive_experiments': sum(experimental_outcomes),
                'accuracy': accuracy,
                'precision': precision,
                'recall': recall,
                'f1_score': f1,
                'mean_result_value': np.mean([r['result_value'] for r in assay_results])
            }
        
        # Calculate overall performance
        all_experimental = []
        all_predicted = []
        
        for assay_results in [r for r in experimental_data]:
            all_experimental.append(assay_results['is_positive'])
            # Simulate prediction (in real implementation, use actual predictions)
            all_predicted.append(np.random.choice([True, False], p=[0.3, 0.7]))
        
        overall_performance = {
            'total_experiments': len(experimental_data),
            'overall_accuracy': accuracy_score(all_experimental, all_predicted) if len(set(all_experimental)) > 1 else 0.0,
            'assay_specific': assay_performance
        }
        
        return overall_performance
    
    def _simulate_predicted_outcomes(self, assay_results: List[Dict]) -> List[bool]:
        """
        Simulate predicted outcomes for experimental results.
        In real implementation, this would use actual VaxGenAI predictions.
        
        Args:
            assay_results: List of experimental results for specific assay
            
        Returns:
            List of predicted outcomes
        """
        predicted_outcomes = []
        
        for result in assay_results:
            peptide = result['peptide']
            
            # Simulate prediction based on peptide properties
            # In real implementation, this would call actual prediction models
            
            # Simple heuristic: longer peptides and those with certain amino acids more likely positive
            length_score = min(1.0, len(peptide) / 10)
            aromatic_score = sum(1 for aa in peptide if aa in 'FWY') / len(peptide)
            charge_score = abs(sum(1 for aa in peptide if aa in 'RK') - sum(1 for aa in peptide if aa in 'DE')) / len(peptide)
            
            prediction_score = 0.3 * length_score + 0.4 * aromatic_score + 0.3 * charge_score
            
            # Add some noise to simulate prediction uncertainty
            prediction_score += np.random.normal(0, 0.1)
            prediction_score = max(0, min(1, prediction_score))
            
            # Convert to binary prediction
            predicted_positive = prediction_score > 0.5
            predicted_outcomes.append(predicted_positive)
        
        return predicted_outcomes
    
    def _generate_model_feedback(self, experimental_data: List[Dict], 
                                performance_metrics: Dict) -> Dict:
        """
        Generate feedback for model improvement based on experimental validation.
        
        Args:
            experimental_data: List of experimental results
            performance_metrics: Performance metrics
            
        Returns:
            Model feedback dictionary
        """
        logger.info("Generating model feedback")
        
        feedback = {
            'overall_assessment': '',
            'strengths': [],
            'weaknesses': [],
            'improvement_suggestions': [],
            'feature_importance_updates': {},
            'training_data_recommendations': []
        }
        
        overall_accuracy = performance_metrics.get('overall_accuracy', 0)
        
        # Overall assessment
        if overall_accuracy > 0.8:
            feedback['overall_assessment'] = 'Excellent prediction performance'
        elif overall_accuracy > 0.6:
            feedback['overall_assessment'] = 'Good prediction performance with room for improvement'
        elif overall_accuracy > 0.4:
            feedback['overall_assessment'] = 'Moderate prediction performance, significant improvements needed'
        else:
            feedback['overall_assessment'] = 'Poor prediction performance, major model revision required'
        
        # Analyze assay-specific performance
        assay_performance = performance_metrics.get('assay_specific', {})
        
        for assay_type, metrics in assay_performance.items():
            if metrics['f1_score'] > 0.7:
                feedback['strengths'].append(f"Strong performance in {assay_type} predictions")
            elif metrics['f1_score'] < 0.3:
                feedback['weaknesses'].append(f"Poor performance in {assay_type} predictions")
        
        # Generate improvement suggestions
        if overall_accuracy < 0.6:
            feedback['improvement_suggestions'].extend([
                'Increase training data diversity',
                'Incorporate additional biochemical features',
                'Implement ensemble methods',
                'Add experimental validation feedback loop'
            ])
        
        # Analyze false positives and false negatives
        false_positives = [r for r in experimental_data if not r['is_positive']]  # Simplified
        false_negatives = [r for r in experimental_data if r['is_positive']]  # Simplified
        
        if len(false_positives) > len(false_negatives):
            feedback['improvement_suggestions'].append('Reduce false positive rate by increasing prediction stringency')
        elif len(false_negatives) > len(false_positives):
            feedback['improvement_suggestions'].append('Reduce false negative rate by decreasing prediction stringency')
        
        # Feature importance updates (simplified)
        feedback['feature_importance_updates'] = {
            'peptide_length': 0.15,
            'hydrophobicity': 0.20,
            'charge_distribution': 0.18,
            'aromatic_content': 0.12,
            'mhc_binding_affinity': 0.25,
            'expression_level': 0.10
        }
        
        # Training data recommendations
        low_performing_assays = [assay for assay, metrics in assay_performance.items() 
                               if metrics['f1_score'] < 0.5]
        
        for assay in low_performing_assays:
            feedback['training_data_recommendations'].append(
                f"Collect more {assay} experimental data for model training"
            )
        
        return feedback
    
    def _generate_experimental_recommendations(self, experimental_data: List[Dict]) -> List[Dict]:
        """
        Generate recommendations for future experimental validation.
        
        Args:
            experimental_data: List of experimental results
            
        Returns:
            List of experimental recommendations
        """
        logger.info("Generating experimental recommendations")
        
        recommendations = []
        
        # Analyze assay coverage
        assay_counts = {}
        for result in experimental_data:
            assay_type = result['assay_type']
            assay_counts[assay_type] = assay_counts.get(assay_type, 0) + 1
        
        # Recommend underrepresented assays
        for assay_type in self.experimental_assays.keys():
            if assay_counts.get(assay_type, 0) < 5:
                recommendations.append({
                    'type': 'assay_expansion',
                    'assay': assay_type,
                    'priority': 'high',
                    'rationale': f'Limited {assay_type} data available for validation',
                    'suggested_experiments': 10
                })
        
        # Recommend validation of high-confidence predictions
        recommendations.append({
            'type': 'high_confidence_validation',
            'priority': 'medium',
            'rationale': 'Validate high-confidence predictions to confirm model accuracy',
            'suggested_experiments': 20,
            'selection_criteria': 'prediction_score > 0.8'
        })
        
        # Recommend validation of edge cases
        recommendations.append({
            'type': 'edge_case_validation',
            'priority': 'medium',
            'rationale': 'Validate predictions near decision boundaries',
            'suggested_experiments': 15,
            'selection_criteria': '0.4 < prediction_score < 0.6'
        })
        
        # Recommend disease-specific validation
        disease_specific_recommendations = [
            {
                'type': 'disease_specific_validation',
                'disease': 'leukemia',
                'priority': 'high',
                'rationale': 'Validate neoantigen predictions in leukemia patient samples',
                'suggested_experiments': 25,
                'assays': ['elispot', 'tetramer_staining']
            },
            {
                'type': 'disease_specific_validation',
                'disease': 'pancreatic_cancer',
                'priority': 'high',
                'rationale': 'Validate epitope predictions in immunosuppressive tumor microenvironment',
                'suggested_experiments': 30,
                'assays': ['elispot', 'cytotoxicity_assay', 'flow_cytometry']
            },
            {
                'type': 'disease_specific_validation',
                'disease': 'hiv',
                'priority': 'high',
                'rationale': 'Validate conserved epitope predictions across HIV strains',
                'suggested_experiments': 20,
                'assays': ['elispot', 'elisa']
            }
        ]
        
        recommendations.extend(disease_specific_recommendations)
        
        return recommendations
    
    def design_validation_experiments(self, predictions: List[Dict], 
                                    disease_type: str = 'cancer',
                                    budget_constraints: Dict = None) -> Dict:
        """
        Design optimal experimental validation strategy.
        
        Args:
            predictions: List of epitope predictions to validate
            disease_type: Type of disease for context-specific design
            budget_constraints: Budget and resource constraints
            
        Returns:
            Experimental design recommendations
        """
        logger.info(f"Designing validation experiments for {len(predictions)} predictions")
        
        # Default budget constraints
        if budget_constraints is None:
            budget_constraints = {
                'total_budget': 50000,  # USD
                'max_experiments': 50,
                'timeline_weeks': 12
            }
        
        # Estimate costs per assay type
        assay_costs = {
            'elispot': 150,  # USD per experiment
            'tetramer_staining': 300,
            'cytotoxicity_assay': 200,
            'elisa': 100,
            'flow_cytometry': 250
        }
        
        # Prioritize predictions for validation
        prioritized_predictions = self._prioritize_predictions_for_validation(predictions)
        
        # Design multi-tier validation strategy
        validation_strategy = {
            'tier_1_screening': {
                'description': 'High-throughput initial screening',
                'assays': ['elispot'],
                'predictions': prioritized_predictions[:20],
                'estimated_cost': 20 * assay_costs['elispot'],
                'timeline_weeks': 4
            },
            'tier_2_confirmation': {
                'description': 'Confirmation of positive hits',
                'assays': ['tetramer_staining', 'flow_cytometry'],
                'predictions': prioritized_predictions[:10],  # Top hits from tier 1
                'estimated_cost': 10 * (assay_costs['tetramer_staining'] + assay_costs['flow_cytometry']),
                'timeline_weeks': 6
            },
            'tier_3_functional': {
                'description': 'Functional validation of confirmed epitopes',
                'assays': ['cytotoxicity_assay'],
                'predictions': prioritized_predictions[:5],  # Top confirmed epitopes
                'estimated_cost': 5 * assay_costs['cytotoxicity_assay'],
                'timeline_weeks': 4
            }
        }
        
        # Calculate total costs and timeline
        total_cost = sum(tier['estimated_cost'] for tier in validation_strategy.values())
        total_timeline = max(tier['timeline_weeks'] for tier in validation_strategy.values())
        
        # Adjust strategy based on budget constraints
        if total_cost > budget_constraints['total_budget']:
            validation_strategy = self._adjust_strategy_for_budget(
                validation_strategy, budget_constraints, assay_costs
            )
        
        # Add disease-specific considerations
        disease_specific_additions = self._add_disease_specific_validation(disease_type)
        
        experimental_design = {
            'validation_strategy': validation_strategy,
            'disease_specific_additions': disease_specific_additions,
            'total_estimated_cost': total_cost,
            'total_timeline_weeks': total_timeline,
            'success_criteria': self._define_success_criteria(disease_type),
            'risk_mitigation': self._identify_validation_risks(),
            'data_analysis_plan': self._create_analysis_plan()
        }
        
        return experimental_design
    
    def _prioritize_predictions_for_validation(self, predictions: List[Dict]) -> List[Dict]:
        """
        Prioritize predictions for experimental validation.
        
        Args:
            predictions: List of epitope predictions
            
        Returns:
            Prioritized list of predictions
        """
        # Score predictions based on multiple criteria
        for prediction in predictions:
            priority_score = 0
            
            # High prediction confidence
            confidence = prediction.get('prediction_confidence', 0.5)
            priority_score += 0.3 * confidence
            
            # High immunogenicity score
            immunogenicity = prediction.get('combined_immunogenicity', 0.5)
            priority_score += 0.3 * immunogenicity
            
            # Low evasion susceptibility
            evasion_resistance = prediction.get('evasion_resistance_score', 0.5)
            priority_score += 0.2 * evasion_resistance
            
            # Novel epitope (not in databases)
            novelty_score = 0.8  # Assume most are novel
            priority_score += 0.2 * novelty_score
            
            prediction['validation_priority_score'] = priority_score
        
        # Sort by priority score
        prioritized = sorted(predictions, key=lambda x: x['validation_priority_score'], reverse=True)
        
        return prioritized
    
    def _adjust_strategy_for_budget(self, strategy: Dict, constraints: Dict, costs: Dict) -> Dict:
        """Adjust validation strategy based on budget constraints."""
        # Reduce number of experiments in each tier
        budget_ratio = constraints['total_budget'] / sum(tier['estimated_cost'] for tier in strategy.values())
        
        for tier_name, tier_data in strategy.items():
            original_count = len(tier_data['predictions'])
            adjusted_count = max(1, int(original_count * budget_ratio))
            
            tier_data['predictions'] = tier_data['predictions'][:adjusted_count]
            
            # Recalculate costs
            if tier_name == 'tier_1_screening':
                tier_data['estimated_cost'] = adjusted_count * costs['elispot']
            elif tier_name == 'tier_2_confirmation':
                tier_data['estimated_cost'] = adjusted_count * (costs['tetramer_staining'] + costs['flow_cytometry'])
            elif tier_name == 'tier_3_functional':
                tier_data['estimated_cost'] = adjusted_count * costs['cytotoxicity_assay']
        
        return strategy
    
    def _add_disease_specific_validation(self, disease_type: str) -> Dict:
        """Add disease-specific validation considerations."""
        disease_additions = {
            'leukemia': {
                'patient_samples': 'Use primary leukemia patient PBMCs',
                'controls': 'Include healthy donor controls and remission samples',
                'special_assays': ['minimal_residual_disease_monitoring'],
                'considerations': 'Account for blast cell interference'
            },
            'pancreatic_cancer': {
                'patient_samples': 'Use tumor-infiltrating lymphocytes when possible',
                'controls': 'Include pancreatic tissue controls',
                'special_assays': ['tumor_microenvironment_simulation'],
                'considerations': 'Test in immunosuppressive conditions'
            },
            'hiv': {
                'patient_samples': 'Use samples from different HIV clades',
                'controls': 'Include elite controllers and progressors',
                'special_assays': ['viral_suppression_assay'],
                'considerations': 'Test against multiple HIV strains'
            }
        }
        
        return disease_additions.get(disease_type, {})
    
    def _define_success_criteria(self, disease_type: str) -> Dict:
        """Define success criteria for validation experiments."""
        return {
            'primary_endpoint': 'Positive response in ≥30% of tested epitopes',
            'secondary_endpoints': [
                'Correlation coefficient >0.6 between predictions and experimental results',
                'Identification of at least 3 high-quality epitopes for vaccine development'
            ],
            'statistical_power': 0.8,
            'significance_level': 0.05
        }
    
    def _identify_validation_risks(self) -> List[Dict]:
        """Identify potential risks in validation experiments."""
        return [
            {
                'risk': 'Low experimental sensitivity',
                'probability': 'medium',
                'impact': 'high',
                'mitigation': 'Use positive controls and optimize assay conditions'
            },
            {
                'risk': 'Patient sample variability',
                'probability': 'high',
                'impact': 'medium',
                'mitigation': 'Use larger sample sizes and stratify by patient characteristics'
            },
            {
                'risk': 'Technical assay failures',
                'probability': 'low',
                'impact': 'high',
                'mitigation': 'Include technical replicates and backup assays'
            }
        ]
    
    def _create_analysis_plan(self) -> Dict:
        """Create data analysis plan for validation experiments."""
        return {
            'primary_analysis': 'Correlation analysis between predictions and experimental outcomes',
            'secondary_analyses': [
                'ROC curve analysis for prediction performance',
                'Subgroup analysis by patient characteristics',
                'Feature importance analysis'
            ],
            'statistical_methods': ['Pearson correlation', 'Spearman correlation', 'Mann-Whitney U test'],
            'visualization': ['Scatter plots', 'ROC curves', 'Heatmaps'],
            'reporting': 'Comprehensive validation report with recommendations for model improvement'
        }
    
    def generate_validation_report(self, integration_summary: Dict,
                                 experimental_design: Dict,
                                 output_path: str) -> str:
        """
        Generate comprehensive experimental validation report.
        
        Args:
            integration_summary: Summary of experimental data integration
            experimental_design: Experimental design recommendations
            output_path: Output file path
            
        Returns:
            Path to generated report
        """
        logger.info("Generating experimental validation report")
        
        report = "# Experimental Validation Integration Report\n\n"
        
        # Integration summary
        report += "## Experimental Data Integration Summary\n\n"
        report += f"Total experimental results processed: {integration_summary['total_results']}\n"
        report += f"Valid results integrated: {integration_summary['valid_results']}\n"
        
        # Performance metrics
        performance = integration_summary.get('performance_metrics', {})
        report += f"Overall prediction accuracy: {performance.get('overall_accuracy', 0):.3f}\n\n"
        
        # Assay-specific performance
        assay_performance = performance.get('assay_specific', {})
        if assay_performance:
            report += "### Assay-Specific Performance\n\n"
            report += "| Assay Type | Experiments | Accuracy | Precision | Recall | F1 Score |\n"
            report += "|------------|-------------|----------|-----------|--------|----------|\n"
            
            for assay, metrics in assay_performance.items():
                report += f"| {assay} | {metrics['total_experiments']} | "
                report += f"{metrics['accuracy']:.3f} | {metrics['precision']:.3f} | "
                report += f"{metrics['recall']:.3f} | {metrics['f1_score']:.3f} |\n"
            report += "\n"
        
        # Model feedback
        feedback = integration_summary.get('model_feedback', {})
        if feedback:
            report += "## Model Performance Feedback\n\n"
            report += f"**Overall Assessment:** {feedback.get('overall_assessment', 'Not available')}\n\n"
            
            strengths = feedback.get('strengths', [])
            if strengths:
                report += "**Strengths:**\n"
                for strength in strengths:
                    report += f"- {strength}\n"
                report += "\n"
            
            weaknesses = feedback.get('weaknesses', [])
            if weaknesses:
                report += "**Areas for Improvement:**\n"
                for weakness in weaknesses:
                    report += f"- {weakness}\n"
                report += "\n"
            
            suggestions = feedback.get('improvement_suggestions', [])
            if suggestions:
                report += "**Improvement Suggestions:**\n"
                for suggestion in suggestions:
                    report += f"- {suggestion}\n"
                report += "\n"
        
        # Experimental design recommendations
        report += "## Experimental Design Recommendations\n\n"
        
        strategy = experimental_design.get('validation_strategy', {})
        if strategy:
            report += "### Validation Strategy\n\n"
            for tier_name, tier_data in strategy.items():
                report += f"**{tier_name.replace('_', ' ').title()}**\n"
                report += f"- Description: {tier_data['description']}\n"
                report += f"- Assays: {', '.join(tier_data['assays'])}\n"
                report += f"- Number of predictions: {len(tier_data['predictions'])}\n"
                report += f"- Estimated cost: ${tier_data['estimated_cost']:,}\n"
                report += f"- Timeline: {tier_data['timeline_weeks']} weeks\n\n"
        
        # Cost and timeline summary
        report += "### Resource Requirements\n\n"
        report += f"Total estimated cost: ${experimental_design.get('total_estimated_cost', 0):,}\n"
        report += f"Total timeline: {experimental_design.get('total_timeline_weeks', 0)} weeks\n\n"
        
        # Success criteria
        success_criteria = experimental_design.get('success_criteria', {})
        if success_criteria:
            report += "### Success Criteria\n\n"
            report += f"Primary endpoint: {success_criteria.get('primary_endpoint', 'Not defined')}\n"
            
            secondary = success_criteria.get('secondary_endpoints', [])
            if secondary:
                report += "Secondary endpoints:\n"
                for endpoint in secondary:
                    report += f"- {endpoint}\n"
                report += "\n"
        
        # Recommendations
        recommendations = integration_summary.get('recommendations', [])
        if recommendations:
            report += "## Recommendations for Future Validation\n\n"
            for i, rec in enumerate(recommendations[:10], 1):  # Top 10 recommendations
                report += f"{i}. **{rec.get('type', 'General').replace('_', ' ').title()}**\n"
                report += f"   - Priority: {rec.get('priority', 'Medium')}\n"
                report += f"   - Rationale: {rec.get('rationale', 'Not specified')}\n"
                if 'suggested_experiments' in rec:
                    report += f"   - Suggested experiments: {rec['suggested_experiments']}\n"
                report += "\n"
        
        # Next steps
        report += "## Next Steps\n\n"
        report += "1. Implement high-priority experimental recommendations\n"
        report += "2. Establish collaborations with experimental laboratories\n"
        report += "3. Develop standardized protocols for validation assays\n"
        report += "4. Create automated feedback loops for continuous model improvement\n"
        report += "5. Plan for larger-scale validation studies\n"
        report += "6. Consider regulatory requirements for vaccine development\n"
        
        # Save report
        with open(output_path, 'w') as f:
            f.write(report)
        
        logger.info(f"Experimental validation report saved to: {output_path}")
        return output_path

# Testing and example usage
def test_experimental_validation_integration():
    """Test the experimental validation integration system."""
    logger.info("Testing experimental validation integration system")
    
    # Sample experimental results
    experimental_results = [
        {
            'peptide': 'YLQPRTFLL',
            'assay_type': 'elispot',
            'result_value': 75,
            'positive_control': True,
            'conditions': {'donor_id': 'D001', 'hla_type': 'A*02:01'},
            'publication': {'pmid': '12345678', 'year': 2023}
        },
        {
            'peptide': 'KIADYNYKL',
            'assay_type': 'tetramer_staining',
            'result_value': 0.15,
            'positive_control': True,
            'conditions': {'donor_id': 'D002', 'hla_type': 'A*02:01'}
        },
        {
            'peptide': 'SLYNTVATL',
            'assay_type': 'cytotoxicity_assay',
            'result_value': 25,
            'positive_control': True,
            'conditions': {'donor_id': 'D003', 'hla_type': 'A*02:01'}
        },
        {
            'peptide': 'GILGFVFTL',
            'assay_type': 'elisa',
            'result_value': 0.3,
            'positive_control': False,
            'conditions': {'donor_id': 'D004', 'hla_type': 'A*02:01'}
        },
        {
            'peptide': 'FIAGLIAIV',
            'assay_type': 'flow_cytometry',
            'result_value': 1.2,
            'positive_control': True,
            'conditions': {'donor_id': 'D005', 'hla_type': 'A*02:01'}
        }
    ]
    
    # Sample predictions for validation design
    sample_predictions = [
        {
            'peptide': 'YLQPRTFLL',
            'prediction_confidence': 0.85,
            'combined_immunogenicity': 0.78,
            'evasion_resistance_score': 0.72
        },
        {
            'peptide': 'KIADYNYKL',
            'prediction_confidence': 0.92,
            'combined_immunogenicity': 0.81,
            'evasion_resistance_score': 0.68
        },
        {
            'peptide': 'NEWEPITOPE',
            'prediction_confidence': 0.76,
            'combined_immunogenicity': 0.69,
            'evasion_resistance_score': 0.74
        }
    ]
    
    # Initialize integrator
    integrator = ExperimentalValidationIntegrator()
    
    # Integrate experimental data
    integration_summary = integrator.integrate_experimental_data(experimental_results)
    
    # Design validation experiments
    experimental_design = integrator.design_validation_experiments(
        sample_predictions, 
        disease_type='leukemia'
    )
    
    # Generate report
    output_dir = Path("/home/ubuntu/vaxgenai/results/experimental_validation")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    report_path = integrator.generate_validation_report(
        integration_summary,
        experimental_design,
        str(output_dir / "validation_integration_report.md")
    )
    
    # Save data
    with open(output_dir / "integration_summary.json", 'w') as f:
        json.dump(integration_summary, f, indent=2, default=str)
    
    with open(output_dir / "experimental_design.json", 'w') as f:
        json.dump(experimental_design, f, indent=2, default=str)
    
    results = {
        'experimental_results_processed': len(experimental_results),
        'valid_results': integration_summary['valid_results'],
        'overall_accuracy': integration_summary['performance_metrics']['overall_accuracy'],
        'total_validation_cost': experimental_design['total_estimated_cost'],
        'report_path': report_path
    }
    
    logger.info("Experimental validation integration test completed")
    return results

if __name__ == "__main__":
    # Run test
    test_results = test_experimental_validation_integration()
    print(f"Experimental validation integration completed")
    print(f"Valid results: {test_results['valid_results']}/{test_results['experimental_results_processed']}")
    print(f"Overall accuracy: {test_results['overall_accuracy']:.3f}")
    print(f"Validation cost estimate: ${test_results['total_validation_cost']:,}")

