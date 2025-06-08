"""
Visualization Module for VaxGenAI

This module generates visualizations for epitope predictions, vaccine candidates,
and other results.
"""

import os
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .utils import get_results_path

logger = logging.getLogger("vaxgenai.visualization")

class Visualization:
    """
    Generates visualizations for VaxGenAI results
    """
    
    def __init__(self):
        """Initialize the visualization module"""
        self.results_path = get_results_path()
        
        # Create visualization directory
        self.viz_path = self.results_path / "visualizations"
        self.viz_path.mkdir(exist_ok=True)
        
        # Set up matplotlib
        plt.style.use('seaborn-v0_8-whitegrid')
        
        logger.info("Visualization module initialized")
    
    def plot_epitope_scores(self, epitopes_df, output_file=None):
        """
        Plot epitope scores
        
        Args:
            epitopes_df: DataFrame of predicted epitopes
            output_file: Path to output file
            
        Returns:
            str: Path to the saved plot
        """
        if output_file is None:
            output_file = self.viz_path / "epitope_scores.png"
        
        # Create figure
        plt.figure(figsize=(12, 6))
        
        # Plot epitope scores
        sns.barplot(x=epitopes_df.index, y='score', hue='type', data=epitopes_df)
        
        # Add labels and title
        plt.xlabel('Epitope Index')
        plt.ylabel('Score')
        plt.title('Predicted Epitope Scores')
        
        # Adjust layout
        plt.tight_layout()
        
        # Save figure
        plt.savefig(output_file, dpi=300)
        plt.close()
        
        logger.info(f"Saved epitope score plot to {output_file}")
        return str(output_file)
    
    def plot_epitope_positions(self, epitopes_df, protein_length, output_file=None):
        """
        Plot epitope positions on the protein sequence
        
        Args:
            epitopes_df: DataFrame of predicted epitopes
            protein_length: Length of the protein sequence
            output_file: Path to output file
            
        Returns:
            str: Path to the saved plot
        """
        if output_file is None:
            output_file = self.viz_path / "epitope_positions.png"
        
        # Create figure
        plt.figure(figsize=(12, 6))
        
        # Create position array
        positions = np.zeros(protein_length)
        
        # Mark epitope positions
        for _, epitope in epitopes_df.iterrows():
            start = epitope['start'] - 1  # Convert to 0-based indexing
            end = epitope['end']
            positions[start:end] += epitope['score']
        
        # Plot positions
        plt.plot(range(1, protein_length + 1), positions)
        
        # Add labels and title
        plt.xlabel('Position')
        plt.ylabel('Epitope Score')
        plt.title('Epitope Positions on Protein Sequence')
        
        # Adjust layout
        plt.tight_layout()
        
        # Save figure
        plt.savefig(output_file, dpi=300)
        plt.close()
        
        logger.info(f"Saved epitope position plot to {output_file}")
        return str(output_file)
    
    def plot_vaccine_scores(self, ranked_candidates, output_file=None):
        """
        Plot vaccine candidate scores
        
        Args:
            ranked_candidates: DataFrame of ranked vaccine candidates
            output_file: Path to output file
            
        Returns:
            str: Path to the saved plot
        """
        if output_file is None:
            output_file = self.viz_path / "vaccine_scores.png"
        
        # Create figure
        plt.figure(figsize=(12, 8))
        
        # Prepare data
        top_candidates = ranked_candidates.head(10)
        
        # Create index for plotting
        index = range(len(top_candidates))
        
        # Set bar width
        bar_width = 0.2
        
        # Plot bars
        plt.bar([i - bar_width for i in index], top_candidates['efficacy_score'], 
                width=bar_width, label='Efficacy')
        plt.bar([i for i in index], top_candidates['population_coverage'], 
                width=bar_width, label='Population Coverage')
        plt.bar([i + bar_width for i in index], top_candidates['manufacturability'], 
                width=bar_width, label='Manufacturability')
        plt.bar([i + 2*bar_width for i in index], top_candidates['overall_score'], 
                width=bar_width, label='Overall Score')
        
        # Add labels and title
        plt.xlabel('Candidate')
        plt.ylabel('Score')
        plt.title('Vaccine Candidate Scores')
        plt.xticks([i for i in index], [f"{i+1}: {t}" for i, t in enumerate(top_candidates['type'])])
        plt.legend()
        
        # Adjust layout
        plt.tight_layout()
        
        # Save figure
        plt.savefig(output_file, dpi=300)
        plt.close()
        
        logger.info(f"Saved vaccine score plot to {output_file}")
        return str(output_file)
    
    def plot_score_radar(self, ranked_candidates, output_file=None):
        """
        Plot radar chart of vaccine candidate scores
        
        Args:
            ranked_candidates: DataFrame of ranked vaccine candidates
            output_file: Path to output file
            
        Returns:
            str: Path to the saved plot
        """
        if output_file is None:
            output_file = self.viz_path / "vaccine_radar.png"
        
        # Create figure
        plt.figure(figsize=(10, 10))
        
        # Prepare data
        top_candidates = ranked_candidates.head(3)
        
        # Categories for radar chart
        categories = ['Efficacy', 'Population Coverage', 'Manufacturability']
        
        # Number of categories
        N = len(categories)
        
        # Create angles for radar chart
        angles = [n / float(N) * 2 * np.pi for n in range(N)]
        angles += angles[:1]  # Close the loop
        
        # Create subplot with polar projection
        ax = plt.subplot(111, polar=True)
        
        # Set category labels
        plt.xticks(angles[:-1], categories)
        
        # Set y-axis limits
        ax.set_ylim(0, 1)
        
        # Plot each candidate
        for i, (_, candidate) in enumerate(top_candidates.iterrows()):
            values = [candidate['efficacy_score'], candidate['population_coverage'], 
                     candidate['manufacturability']]
            values += values[:1]  # Close the loop
            
            # Plot values
            ax.plot(angles, values, linewidth=2, linestyle='solid', 
                   label=f"{i+1}: {candidate['type']}")
            ax.fill(angles, values, alpha=0.1)
        
        # Add legend
        plt.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1))
        
        # Add title
        plt.title('Vaccine Candidate Comparison')
        
        # Save figure
        plt.savefig(output_file, dpi=300)
        plt.close()
        
        logger.info(f"Saved vaccine radar plot to {output_file}")
        return str(output_file)
    
    def plot_epitope_types(self, epitopes_df, output_file=None):
        """
        Plot distribution of epitope types
        
        Args:
            epitopes_df: DataFrame of predicted epitopes
            output_file: Path to output file
            
        Returns:
            str: Path to the saved plot
        """
        if output_file is None:
            output_file = self.viz_path / "epitope_types.png"
        
        # Create figure
        plt.figure(figsize=(8, 6))
        
        # Count epitope types
        type_counts = epitopes_df['type'].value_counts()
        
        # Plot pie chart
        plt.pie(type_counts, labels=type_counts.index, autopct='%1.1f%%', 
               shadow=True, startangle=90)
        
        # Add title
        plt.title('Distribution of Epitope Types')
        
        # Equal aspect ratio ensures that pie is drawn as a circle
        plt.axis('equal')
        
        # Save figure
        plt.savefig(output_file, dpi=300)
        plt.close()
        
        logger.info(f"Saved epitope type plot to {output_file}")
        return str(output_file)
    
    def create_visualization_report(self, epitopes_df, ranked_candidates, protein_length):
        """
        Create a comprehensive visualization report
        
        Args:
            epitopes_df: DataFrame of predicted epitopes
            ranked_candidates: DataFrame of ranked vaccine candidates
            protein_length: Length of the protein sequence
            
        Returns:
            list: Paths to the generated plots
        """
        plots = []
        
        # Plot epitope scores
        epitope_scores_plot = self.plot_epitope_scores(epitopes_df)
        plots.append(epitope_scores_plot)
        
        # Plot epitope positions
        epitope_positions_plot = self.plot_epitope_positions(epitopes_df, protein_length)
        plots.append(epitope_positions_plot)
        
        # Plot epitope types
        epitope_types_plot = self.plot_epitope_types(epitopes_df)
        plots.append(epitope_types_plot)
        
        # Plot vaccine scores
        vaccine_scores_plot = self.plot_vaccine_scores(ranked_candidates)
        plots.append(vaccine_scores_plot)
        
        # Plot score radar
        vaccine_radar_plot = self.plot_score_radar(ranked_candidates)
        plots.append(vaccine_radar_plot)
        
        # Create HTML report
        html_report = self.create_html_report(plots)
        
        return plots, html_report
    
    def create_html_report(self, plot_paths):
        """
        Create an HTML report with all visualizations
        
        Args:
            plot_paths: List of paths to plot files
            
        Returns:
            str: Path to the HTML report
        """
        # Create HTML content
        html = []
        
        # Add header
        html.append("<!DOCTYPE html>")
        html.append("<html>")
        html.append("<head>")
        html.append("    <title>VaxGenAI Visualization Report</title>")
        html.append("    <style>")
        html.append("        body { font-family: Arial, sans-serif; margin: 20px; }")
        html.append("        h1 { color: #2c3e50; }")
        html.append("        h2 { color: #3498db; }")
        html.append("        .plot { margin: 20px 0; text-align: center; }")
        html.append("        .plot img { max-width: 100%; border: 1px solid #ddd; }")
        html.append("    </style>")
        html.append("</head>")
        html.append("<body>")
        
        # Add title
        html.append("    <h1>VaxGenAI Visualization Report</h1>")
        
        # Add plots
        for i, plot_path in enumerate(plot_paths):
            plot_name = Path(plot_path).stem.replace('_', ' ').title()
            
            html.append(f"    <div class='plot'>")
            html.append(f"        <h2>{plot_name}</h2>")
            html.append(f"        <img src='{plot_path}' alt='{plot_name}'>")
            html.append(f"    </div>")
        
        # Add footer
        html.append("</body>")
        html.append("</html>")
        
        # Save HTML report
        report_file = self.results_path / "visualization_report.html"
        with open(report_file, 'w') as f:
            f.write('\n'.join(html))
        
        logger.info(f"Created HTML visualization report: {report_file}")
        
        return str(report_file)

