"""
Vaccine Design Module for VaxGenAI

This module designs vaccine candidates based on predicted epitopes, including
mRNA sequences, protein subunits, and peptide vaccines.
"""

import os
import logging
import random
import numpy as np
import pandas as pd
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from .utils import get_results_path

logger = logging.getLogger("vaxgenai.vaccine_design")

class VaccineDesign:
    """
    Designs vaccine candidates based on predicted epitopes
    """
    
    def __init__(self):
        """Initialize the vaccine design module"""
        self.results_path = get_results_path()
        self.vaccine_candidates = []
        logger.info("Vaccine design module initialized")
    
    def design_subunit_vaccine(self, epitopes_df, max_length=200, min_epitopes=5):
        """
        Design a subunit vaccine by combining epitopes
        
        Args:
            epitopes_df: DataFrame of predicted epitopes
            max_length: Maximum length of the vaccine construct
            min_epitopes: Minimum number of epitopes to include
            
        Returns:
            dict: Vaccine candidate information
        """
        # Sort epitopes by score
        sorted_epitopes = epitopes_df.sort_values('score', ascending=False)
        
        # Initialize vaccine construct
        vaccine_seq = ""
        included_epitopes = []
        
        # Linker sequence between epitopes (commonly used GPGPG linker)
        linker = "GPGPG"
        
        # Add epitopes until we reach max_length or run out of epitopes
        for _, epitope in sorted_epitopes.iterrows():
            # Check if adding this epitope would exceed max_length
            if len(vaccine_seq) + len(epitope['sequence']) + len(linker) <= max_length:
                # Add linker if this is not the first epitope
                if vaccine_seq:
                    vaccine_seq += linker
                
                # Add epitope
                vaccine_seq += epitope['sequence']
                included_epitopes.append(epitope)
            
            # Check if we have enough epitopes
            if len(included_epitopes) >= min_epitopes and len(vaccine_seq) >= max_length / 2:
                break
        
        # Create vaccine candidate
        vaccine = {
            'type': 'subunit',
            'sequence': vaccine_seq,
            'length': len(vaccine_seq),
            'epitopes': included_epitopes,
            'num_epitopes': len(included_epitopes)
        }
        
        logger.info(f"Designed subunit vaccine with {len(included_epitopes)} epitopes, length: {len(vaccine_seq)}")
        self.vaccine_candidates.append(vaccine)
        
        return vaccine
    
    def design_mrna_vaccine(self, protein_sequence):
        """
        Design an mRNA vaccine for the given protein sequence
        
        Args:
            protein_sequence: Amino acid sequence for the vaccine
            
        Returns:
            dict: Vaccine candidate information
        """
        # Codon optimization table (human optimized codons)
        codon_table = {
            'A': ['GCC', 'GCT', 'GCA', 'GCG'],
            'C': ['TGC', 'TGT'],
            'D': ['GAC', 'GAT'],
            'E': ['GAG', 'GAA'],
            'F': ['TTC', 'TTT'],
            'G': ['GGC', 'GGT', 'GGA', 'GGG'],
            'H': ['CAC', 'CAT'],
            'I': ['ATC', 'ATT', 'ATA'],
            'K': ['AAG', 'AAA'],
            'L': ['CTG', 'CTC', 'TTG', 'CTT', 'TTA', 'CTA'],
            'M': ['ATG'],
            'N': ['AAC', 'AAT'],
            'P': ['CCC', 'CCT', 'CCA', 'CCG'],
            'Q': ['CAG', 'CAA'],
            'R': ['CGC', 'AGG', 'AGA', 'CGG', 'CGA', 'CGT'],
            'S': ['AGC', 'TCC', 'TCT', 'AGT', 'TCA', 'TCG'],
            'T': ['ACC', 'ACT', 'ACA', 'ACG'],
            'V': ['GTG', 'GTC', 'GTT', 'GTA'],
            'W': ['TGG'],
            'Y': ['TAC', 'TAT'],
            '*': ['TGA', 'TAA', 'TAG']
        }
        
        # Human optimized codon weights
        codon_weights = {
            'A': [0.4, 0.3, 0.2, 0.1],
            'C': [0.6, 0.4],
            'D': [0.6, 0.4],
            'E': [0.6, 0.4],
            'F': [0.6, 0.4],
            'G': [0.4, 0.3, 0.2, 0.1],
            'H': [0.6, 0.4],
            'I': [0.5, 0.4, 0.1],
            'K': [0.6, 0.4],
            'L': [0.4, 0.2, 0.2, 0.1, 0.05, 0.05],
            'M': [1.0],
            'N': [0.6, 0.4],
            'P': [0.4, 0.3, 0.2, 0.1],
            'Q': [0.6, 0.4],
            'R': [0.3, 0.2, 0.2, 0.1, 0.1, 0.1],
            'S': [0.3, 0.2, 0.2, 0.1, 0.1, 0.1],
            'T': [0.4, 0.3, 0.2, 0.1],
            'V': [0.4, 0.3, 0.2, 0.1],
            'W': [1.0],
            'Y': [0.6, 0.4],
            '*': [0.6, 0.3, 0.1]
        }
        
        # 5' UTR with Kozak sequence
        utr_5 = "GGGAAATAAGAGAGAAAAGAAGAGTAAGAAGAAATATAAGACCACC"
        
        # Start codon
        start_codon = "ATG"
        
        # 3' UTR
        utr_3 = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAA"
        
        # Optimize codons for the protein sequence
        mrna_seq = ""
        for aa in protein_sequence:
            if aa in codon_table:
                # Choose a codon based on weights
                codons = codon_table[aa]
                weights = codon_weights[aa]
                codon = random.choices(codons, weights=weights, k=1)[0]
                mrna_seq += codon
            else:
                # For unknown amino acids, use a stop codon
                mrna_seq += "TGA"
        
        # Add stop codon
        mrna_seq += "TGA"
        
        # Assemble the full mRNA sequence
        full_mrna = utr_5 + start_codon + mrna_seq + utr_3
        
        # Create vaccine candidate
        vaccine = {
            'type': 'mRNA',
            'protein_sequence': protein_sequence,
            'mrna_sequence': full_mrna,
            'length': len(full_mrna),
            'gc_content': (full_mrna.count('G') + full_mrna.count('C')) / len(full_mrna)
        }
        
        logger.info(f"Designed mRNA vaccine, length: {len(full_mrna)}, GC content: {vaccine['gc_content']:.2f}")
        self.vaccine_candidates.append(vaccine)
        
        return vaccine
    
    def design_peptide_vaccine(self, epitopes_df, num_epitopes=5):
        """
        Design a peptide vaccine by selecting top epitopes
        
        Args:
            epitopes_df: DataFrame of predicted epitopes
            num_epitopes: Number of epitopes to include
            
        Returns:
            dict: Vaccine candidate information
        """
        # Sort epitopes by score and take the top ones
        top_epitopes = epitopes_df.sort_values('score', ascending=False).head(num_epitopes)
        
        # Create vaccine candidate
        vaccine = {
            'type': 'peptide',
            'epitopes': top_epitopes.to_dict('records'),
            'num_epitopes': len(top_epitopes)
        }
        
        logger.info(f"Designed peptide vaccine with {len(top_epitopes)} epitopes")
        self.vaccine_candidates.append(vaccine)
        
        return vaccine
    
    def save_vaccine_candidates(self, output_prefix=None):
        """
        Save vaccine candidates to files
        
        Args:
            output_prefix: Prefix for output files
            
        Returns:
            dict: Paths to saved files
        """
        if not self.vaccine_candidates:
            logger.error("No vaccine candidates to save")
            return None
        
        if output_prefix is None:
            output_prefix = "vaccine_candidate"
        
        output_files = {}
        
        for i, vaccine in enumerate(self.vaccine_candidates):
            # Create base filename
            base_name = f"{output_prefix}_{vaccine['type']}_{i+1}"
            
            # Save as FASTA
            if vaccine['type'] == 'subunit' or vaccine['type'] == 'peptide':
                # For subunit and peptide vaccines, save protein sequences
                if vaccine['type'] == 'subunit':
                    seq = vaccine['sequence']
                    desc = f"Subunit vaccine with {vaccine['num_epitopes']} epitopes"
                else:  # peptide
                    # For peptide vaccines, save each epitope separately
                    records = []
                    for j, epitope in enumerate(vaccine['epitopes']):
                        record = SeqRecord(
                            Seq(epitope['sequence']),
                            id=f"{base_name}_epitope_{j+1}",
                            description=f"Peptide epitope, score={epitope['score']:.3f}"
                        )
                        records.append(record)
                    
                    fasta_file = self.results_path / f"{base_name}.fasta"
                    SeqIO.write(records, fasta_file, "fasta")
                    output_files[f"{vaccine['type']}_fasta"] = str(fasta_file)
                    continue
                
                record = SeqRecord(
                    Seq(seq),
                    id=base_name,
                    description=desc
                )
                
                fasta_file = self.results_path / f"{base_name}.fasta"
                SeqIO.write([record], fasta_file, "fasta")
                output_files[f"{vaccine['type']}_fasta"] = str(fasta_file)
            
            elif vaccine['type'] == 'mRNA':
                # For mRNA vaccines, save both protein and mRNA sequences
                protein_record = SeqRecord(
                    Seq(vaccine['protein_sequence']),
                    id=f"{base_name}_protein",
                    description="Protein sequence for mRNA vaccine"
                )
                
                mrna_record = SeqRecord(
                    Seq(vaccine['mrna_sequence']),
                    id=f"{base_name}_mrna",
                    description=f"mRNA sequence, GC content: {vaccine['gc_content']:.2f}"
                )
                
                protein_file = self.results_path / f"{base_name}_protein.fasta"
                mrna_file = self.results_path / f"{base_name}_mrna.fasta"
                
                SeqIO.write([protein_record], protein_file, "fasta")
                SeqIO.write([mrna_record], mrna_file, "fasta")
                
                output_files["mrna_protein_fasta"] = str(protein_file)
                output_files["mrna_sequence_fasta"] = str(mrna_file)
            
            # Save as JSON
            json_file = self.results_path / f"{base_name}.json"
            pd.Series(vaccine).to_json(json_file)
            output_files[f"{vaccine['type']}_json"] = str(json_file)
        
        logger.info(f"Saved {len(self.vaccine_candidates)} vaccine candidates")
        return output_files
    
    def design_vaccines(self, epitopes_df):
        """
        Design multiple types of vaccines based on epitopes
        
        Args:
            epitopes_df: DataFrame of predicted epitopes
            
        Returns:
            list: List of vaccine candidates
        """
        # Design subunit vaccine
        subunit_vaccine = self.design_subunit_vaccine(epitopes_df)
        
        # Design mRNA vaccine based on the subunit sequence
        mrna_vaccine = self.design_mrna_vaccine(subunit_vaccine['sequence'])
        
        # Design peptide vaccine
        peptide_vaccine = self.design_peptide_vaccine(epitopes_df)
        
        # Save all vaccine candidates
        output_files = self.save_vaccine_candidates()
        
        return self.vaccine_candidates

