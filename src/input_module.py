"""
Input Module for VaxGenAI

This module handles the input of pathogen sequences (DNA/RNA/protein) and performs
initial processing and validation.
"""

import os
import logging
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .utils import get_data_path

logger = logging.getLogger("vaxgenai.input")

class InputModule:
    """
    Handles the input and preprocessing of pathogen sequences
    """
    
    def __init__(self):
        """Initialize the input module"""
        self.data_path = get_data_path()
        logger.info(f"Input module initialized with data path: {self.data_path}")
    
    def validate_fasta(self, file_path):
        """
        Validate that a file is in FASTA format
        
        Args:
            file_path: Path to the FASTA file
            
        Returns:
            bool: True if valid, False otherwise
        """
        try:
            with open(file_path, "r") as handle:
                fasta = SeqIO.parse(handle, "fasta")
                return any(fasta)  # Returns True if there's at least one record
        except Exception as e:
            logger.error(f"Error validating FASTA file: {e}")
            return False
    
    def load_sequence(self, file_path):
        """
        Load a sequence from a FASTA file
        
        Args:
            file_path: Path to the FASTA file
            
        Returns:
            list: List of SeqRecord objects
        """
        if not self.validate_fasta(file_path):
            logger.error(f"Invalid FASTA file: {file_path}")
            return None
        
        try:
            records = list(SeqIO.parse(file_path, "fasta"))
            logger.info(f"Loaded {len(records)} sequences from {file_path}")
            return records
        except Exception as e:
            logger.error(f"Error loading sequence: {e}")
            return None
    
    def translate_to_protein(self, nucleotide_record):
        """
        Translate a nucleotide sequence to protein
        
        Args:
            nucleotide_record: SeqRecord object with nucleotide sequence
            
        Returns:
            SeqRecord: Protein sequence record
        """
        try:
            # Determine if the sequence is DNA or RNA
            seq_type = "DNA" if "T" in str(nucleotide_record.seq).upper() else "RNA"
            logger.info(f"Detected sequence type: {seq_type}")
            
            # Translate the sequence
            protein_seq = nucleotide_record.seq.translate()
            
            # Create a new SeqRecord for the protein sequence
            protein_record = SeqRecord(
                protein_seq,
                id=nucleotide_record.id + "_translated",
                name=nucleotide_record.name + "_protein",
                description=f"Translated from {seq_type} sequence"
            )
            
            logger.info(f"Translated sequence: {len(protein_seq)} amino acids")
            return protein_record
        except Exception as e:
            logger.error(f"Error translating sequence: {e}")
            return None
    
    def save_sequence(self, record, file_name=None):
        """
        Save a sequence to a FASTA file
        
        Args:
            record: SeqRecord object or list of SeqRecord objects
            file_name: Name of the output file (without path)
            
        Returns:
            str: Path to the saved file
        """
        if file_name is None:
            file_name = f"{record.id}.fasta" if hasattr(record, 'id') else "sequence.fasta"
        
        output_path = self.data_path / file_name
        
        try:
            SeqIO.write(record, output_path, "fasta")
            logger.info(f"Saved sequence to {output_path}")
            return output_path
        except Exception as e:
            logger.error(f"Error saving sequence: {e}")
            return None
    
    def fetch_from_database(self, accession, database="genbank"):
        """
        Fetch a sequence from a public database
        
        Args:
            accession: Accession number
            database: Database to fetch from (default: genbank)
            
        Returns:
            SeqRecord: Sequence record
        """
        try:
            from Bio import Entrez
            Entrez.email = "example@example.com"  # Should be a real email in production
            
            logger.info(f"Fetching {accession} from {database}")
            handle = Entrez.efetch(db=database, id=accession, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            handle.close()
            
            # Save the fetched sequence
            output_path = self.data_path / f"{accession}.fasta"
            SeqIO.write(record, output_path, "fasta")
            
            logger.info(f"Fetched and saved {accession} ({len(record.seq)} bp)")
            return record
        except Exception as e:
            logger.error(f"Error fetching sequence {accession}: {e}")
            return None

