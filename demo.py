"""
Demo script for VaxGenAI

This script demonstrates the functionality of the VaxGenAI system using
the SARS-CoV-2 spike protein sequence.
"""

import os
import sys
import logging
from pathlib import Path

# Add the project root to the path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.input_module import InputModule
from src.epitope_prediction import EpitopePrediction
from src.vaccine_design import VaccineDesign
from src.safety_filter import SafetyFilter
from src.ranking import RankingEngine
from src.visualization import Visualization
from src.utils import setup_logging, get_data_path, get_results_path

def main():
    """Run the VaxGenAI demo"""
    # Setup logging
    logger = setup_logging()
    logger.info("Starting VaxGenAI demo")
    
    # Initialize modules
    input_module = InputModule()
    epitope_predictor = EpitopePrediction()
    vaccine_designer = VaccineDesign()
    safety_filter = SafetyFilter()
    ranking_engine = RankingEngine()
    visualizer = Visualization()
    
    # Load SARS-CoV-2 spike protein sequence
    spike_protein_path = get_data_path() / "spike_protein.fasta"
    logger.info(f"Loading sequence from {spike_protein_path}")
    
    protein_records = input_module.load_sequence(spike_protein_path)
    if not protein_records:
        logger.error("Failed to load sequence")
        return
    
    protein_record = protein_records[0]
    logger.info(f"Loaded sequence: {protein_record.id}, length: {len(protein_record.seq)}")
    
    # Predict epitopes
    logger.info("Predicting epitopes")
    epitopes_df = epitope_predictor.predict_epitopes(protein_record)
    logger.info(f"Predicted {len(epitopes_df)} epitopes")
    
    # Get top epitopes
    top_epitopes = epitope_predictor.get_top_epitopes(n=10)
    logger.info(f"Top 10 epitopes:\n{top_epitopes}")
    
    # Design vaccines
    logger.info("Designing vaccine candidates")
    vaccine_candidates = vaccine_designer.design_vaccines(epitopes_df)
    logger.info(f"Designed {len(vaccine_candidates)} vaccine candidates")
    
    # Filter vaccines for safety
    logger.info("Filtering vaccine candidates for safety")
    safe_candidates, filtered_candidates = safety_filter.filter_vaccine_candidates(vaccine_candidates)
    logger.info(f"Safe candidates: {len(safe_candidates)}, Filtered candidates: {len(filtered_candidates)}")
    
    # Rank vaccine candidates
    logger.info("Ranking vaccine candidates")
    ranked_candidates = ranking_engine.rank_candidates(safe_candidates)
    logger.info(f"Ranked {len(ranked_candidates)} vaccine candidates")
    
    # Generate report
    logger.info("Generating ranking report")
    report_path = ranking_engine.generate_report(ranked_candidates)
    logger.info(f"Generated ranking report: {report_path}")
    
    # Create visualizations
    logger.info("Creating visualizations")
    plots, html_report = visualizer.create_visualization_report(
        epitopes_df, ranked_candidates, len(protein_record.seq)
    )
    logger.info(f"Created {len(plots)} plots and HTML report: {html_report}")
    
    # Print summary
    print("\n" + "="*80)
    print("VaxGenAI Demo Summary")
    print("="*80)
    print(f"Sequence: {protein_record.id}")
    print(f"Sequence length: {len(protein_record.seq)} amino acids")
    print(f"Epitopes predicted: {len(epitopes_df)}")
    print(f"Vaccine candidates designed: {len(vaccine_candidates)}")
    print(f"Safe candidates: {len(safe_candidates)}")
    print(f"Top ranked vaccine type: {ranked_candidates.iloc[0]['type']}")
    print(f"Top ranked vaccine score: {ranked_candidates.iloc[0]['overall_score']:.3f}")
    print("\nResults saved to:")
    print(f"- Ranking report: {report_path}")
    print(f"- Visualization report: {html_report}")
    print("="*80)
    
    logger.info("VaxGenAI demo completed successfully")

if __name__ == "__main__":
    main()

