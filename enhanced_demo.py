"""
Demo script for enhanced VaxGenAI modules

This script demonstrates the functionality of the enhanced VaxGenAI modules
using the SARS-CoV-2 spike protein sequence.
"""

import os
import sys
import logging
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO

# Add the project root to the path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.utils import setup_logging, get_data_path, get_results_path
from src.enhanced.epitope_prediction import (
    EnhancedEpitopePrediction,
    AlphaFoldIntegration,
    PopulationCoverageAnalysis
)
from src.enhanced.safety_filter import EnhancedSafetyFilter
from src.vaccine_design import VaccineDesign
from src.ranking import RankingEngine
from src.visualization import Visualization

def main():
    """Run the enhanced VaxGenAI demo"""
    # Setup logging
    logger = setup_logging()
    logger.info("Starting enhanced VaxGenAI demo")
    
    # Create results directory
    results_dir = get_results_path() / "enhanced"
    os.makedirs(results_dir, exist_ok=True)
    
    # Initialize modules
    epitope_predictor = EnhancedEpitopePrediction()
    alphafold = AlphaFoldIntegration()
    population_coverage = PopulationCoverageAnalysis()
    safety_filter = EnhancedSafetyFilter()
    vaccine_designer = VaccineDesign()
    ranking_engine = RankingEngine()
    visualizer = Visualization()
    
    # Load SARS-CoV-2 spike protein sequence
    spike_protein_path = get_data_path() / "spike_protein.fasta"
    logger.info(f"Loading sequence from {spike_protein_path}")
    
    try:
        protein_records = list(SeqIO.parse(spike_protein_path, "fasta"))
        if not protein_records:
            logger.error("Failed to load sequence")
            return
        
        protein_record = protein_records[0]
        logger.info(f"Loaded sequence: {protein_record.id}, length: {len(protein_record.seq)}")
        
        # Predict epitopes using enhanced prediction
        logger.info("Predicting epitopes using enhanced prediction")
        epitopes_df = epitope_predictor.predict_epitopes(protein_record)
        logger.info(f"Predicted {len(epitopes_df)} epitopes")
        
        # Save epitopes to CSV
        epitopes_csv = results_dir / "enhanced_epitopes.csv"
        epitope_predictor.save_results(epitopes_df, epitopes_csv)
        
        # Get top epitopes
        top_epitopes = epitope_predictor.get_top_epitopes(n=20)
        logger.info(f"Top 20 epitopes:\n{pd.DataFrame(top_epitopes).head()}")
        
        # Calculate epitope coverage
        coverage = epitope_predictor.calculate_epitope_coverage(top_epitopes)
        logger.info(f"Epitope coverage: {coverage:.2%}")
        
        # Integrate with AlphaFold for structure prediction
        logger.info("Integrating with AlphaFold for structure prediction")
        # Note: In a real implementation, we would use the actual UniProt ID for the spike protein
        # For demonstration, we'll use a placeholder
        uniprot_id = "P0DTC2"  # SARS-CoV-2 spike protein UniProt ID
        
        # Get structure information
        structure_info = alphafold.get_structure(uniprot_id)
        if structure_info:
            logger.info(f"Retrieved AlphaFold structure for {uniprot_id}")
            
            # Calculate accessibility for epitopes
            epitopes_with_accessibility = alphafold.calculate_accessibility(uniprot_id, epitopes_df)
            
            # Save updated epitopes to CSV
            accessibility_csv = results_dir / "epitopes_with_accessibility.csv"
            epitopes_with_accessibility.to_csv(accessibility_csv, index=False)
            logger.info(f"Saved epitopes with accessibility to {accessibility_csv}")
        else:
            logger.warning(f"Could not retrieve AlphaFold structure for {uniprot_id}")
            epitopes_with_accessibility = epitopes_df
        
        # Calculate population coverage
        logger.info("Calculating population coverage")
        coverage_data = population_coverage.calculate_global_coverage(epitopes_df)
        
        # Generate population coverage report
        coverage_report = results_dir / "population_coverage_report.md"
        population_coverage.generate_coverage_report(coverage_data, coverage_report)
        logger.info(f"Generated population coverage report: {coverage_report}")
        
        # Create compatibility layer for vaccine design and visualization
        logger.info("Creating compatibility layer for vaccine design and visualization")
        epitopes_for_vaccine = epitopes_with_accessibility.copy()
        epitopes_for_vaccine['score'] = epitopes_for_vaccine['prediction_score']
        epitopes_for_vaccine['sequence'] = epitopes_for_vaccine['peptide']
        epitopes_for_vaccine['type'] = epitopes_for_vaccine['epitope_type']
        
        # Design vaccine candidates
        logger.info("Designing vaccine candidates")
        vaccine_candidates = vaccine_designer.design_vaccines(epitopes_for_vaccine)
        logger.info(f"Designed {len(vaccine_candidates)} vaccine candidates")
        
        # Save vaccine candidates to CSV for inspection
        vaccine_csv = results_dir / "vaccine_candidates.csv"
        pd.DataFrame(vaccine_candidates).to_csv(vaccine_csv, index=False)
        logger.info(f"Saved vaccine candidates to {vaccine_csv}")
        
        # Inspect the structure of vaccine candidates
        logger.info("Inspecting vaccine candidates structure")
        vaccine_df = pd.DataFrame(vaccine_candidates)
        logger.info(f"Vaccine candidates columns: {vaccine_df.columns.tolist()}")
        
        # Filter vaccines for safety
        logger.info("Filtering vaccine candidates for safety")
        safe_candidates, filtered_candidates = safety_filter.filter_vaccine_candidates(vaccine_candidates)
        logger.info(f"Safe candidates: {len(safe_candidates)}, Filtered candidates: {len(filtered_candidates)}")
        
        # Save safe and filtered candidates
        safe_csv = results_dir / "safe_candidates.csv"
        filtered_csv = results_dir / "filtered_candidates.csv"
        pd.DataFrame(safe_candidates).to_csv(safe_csv, index=False)
        pd.DataFrame(filtered_candidates).to_csv(filtered_csv, index=False)
        
        # Rank vaccine candidates
        logger.info("Ranking vaccine candidates")
        if safe_candidates:
            ranked_candidates = ranking_engine.rank_candidates(safe_candidates)
            logger.info(f"Ranked {len(ranked_candidates)} vaccine candidates")
            
            # Generate report
            logger.info("Generating ranking report")
            report_path = ranking_engine.generate_report(ranked_candidates)
            logger.info(f"Generated ranking report: {report_path}")
            
            # Create visualizations
            logger.info("Creating visualizations")
            try:
                plots, html_report = visualizer.create_visualization_report(
                    epitopes_for_vaccine, ranked_candidates, len(protein_record.seq)
                )
                logger.info(f"Created {len(plots)} plots and HTML report: {html_report}")
            except Exception as e:
                logger.error(f"Error creating visualization report: {e}")
                plots = []
                html_report = "Error creating visualization report"
        else:
            logger.warning("No safe vaccine candidates to rank")
            ranked_candidates = pd.DataFrame()
            report_path = "No ranking report generated"
            html_report = "No visualization report generated"
        
        # Create population coverage visualization
        plt.figure(figsize=(10, 6))
        populations = list(coverage_data.keys())
        coverages = [data['coverage'] for data in coverage_data.values()]
        plt.bar(populations, coverages)
        plt.xlabel('Population')
        plt.ylabel('Coverage')
        plt.title('Population Coverage of Predicted Epitopes')
        plt.xticks(rotation=45)
        plt.tight_layout()
        
        # Save population coverage plot
        coverage_plot = results_dir / "population_coverage.png"
        plt.savefig(coverage_plot)
        logger.info(f"Saved population coverage plot to {coverage_plot}")
        
        # Create custom visualization for epitope distribution
        logger.info("Creating custom visualization for epitope distribution")
        plt.figure(figsize=(12, 6))
        
        # Filter to top 100 epitopes for better visualization
        top_100_epitopes = epitopes_df.sort_values('prediction_score', ascending=False).head(100)
        
        # Create scatter plot of epitope positions
        plt.scatter(
            top_100_epitopes['start_position'], 
            top_100_epitopes['prediction_score'],
            c=top_100_epitopes['epitope_type'].map({'T-cell': 'blue', 'B-cell': 'red'}),
            alpha=0.7
        )
        
        plt.xlabel('Position in Protein Sequence')
        plt.ylabel('Prediction Score')
        plt.title('Top 100 Predicted Epitopes by Position')
        plt.legend(['T-cell Epitopes', 'B-cell Epitopes'])
        plt.grid(True, linestyle='--', alpha=0.7)
        
        # Save epitope distribution plot
        epitope_plot = results_dir / "epitope_positions.png"
        plt.savefig(epitope_plot)
        logger.info(f"Saved epitope distribution plot to {epitope_plot}")
        
        # Print summary
        print("\n" + "="*80)
        print("Enhanced VaxGenAI Demo Summary")
        print("="*80)
        print(f"Sequence: {protein_record.id}")
        print(f"Sequence length: {len(protein_record.seq)} amino acids")
        print(f"Epitopes predicted: {len(epitopes_df)}")
        print(f"Epitope coverage: {coverage:.2%}")
        print(f"Population coverage (World): {coverage_data['World']['coverage']:.2%}")
        print(f"Vaccine candidates designed: {len(vaccine_candidates)}")
        print(f"Safe candidates: {len(safe_candidates)}")
        
        if not ranked_candidates.empty:
            print(f"Top ranked vaccine type: {ranked_candidates.iloc[0]['type']}")
            print(f"Top ranked vaccine score: {ranked_candidates.iloc[0]['overall_score']:.3f}")
        
        print("\nResults saved to:")
        print(f"- Enhanced epitopes: {epitopes_csv}")
        print(f"- Population coverage report: {coverage_report}")
        print(f"- Population coverage plot: {coverage_plot}")
        print(f"- Epitope distribution plot: {epitope_plot}")
        print(f"- Vaccine candidates: {vaccine_csv}")
        print(f"- Safe candidates: {safe_csv}")
        if not ranked_candidates.empty:
            print(f"- Ranking report: {report_path}")
            print(f"- Visualization report: {html_report}")
        print("="*80)
        
        logger.info("Enhanced VaxGenAI demo completed successfully")
    
    except Exception as e:
        logger.error(f"Error in enhanced VaxGenAI demo: {e}", exc_info=True)

if __name__ == "__main__":
    main()

