# VaxGenAI API Reference

Complete API documentation for the VaxGenAI vaccine design platform.

## Core Modules

### InputModule

Handles input data processing and validation.

```python
from src.input_module import InputModule

input_module = InputModule()
```

#### Methods

##### `load_fasta(file_path: str) -> Dict`
Load protein sequences from FASTA file.

**Parameters:**
- `file_path` (str): Path to FASTA file

**Returns:**
- Dict containing protein ID and sequence

**Example:**
```python
protein_data = input_module.load_fasta("data/spike_protein.fasta")
print(protein_data['id'])      # 'spike_protein'
print(protein_data['sequence']) # 'MFVFLVLLPLVSSQ...'
```

##### `validate_sequence(sequence: str) -> bool`
Validate protein sequence format.

**Parameters:**
- `sequence` (str): Amino acid sequence

**Returns:**
- bool: True if valid, False otherwise

##### `load_fasta_chunks(file_path: str, chunk_size: int = 100) -> Iterator`
Load large FASTA files in chunks for memory efficiency.

**Parameters:**
- `file_path` (str): Path to FASTA file
- `chunk_size` (int): Number of sequences per chunk

**Returns:**
- Iterator yielding chunks of protein data

---

### EpitopePrediction

Core epitope prediction functionality.

```python
from src.epitope_prediction import EpitopePrediction

predictor = EpitopePrediction()
```

#### Methods

##### `predict_epitopes(protein_data: Dict, config: Dict = None) -> List[Dict]`
Predict T-cell and B-cell epitopes from protein sequence.

**Parameters:**
- `protein_data` (Dict): Protein data from InputModule
- `config` (Dict, optional): Configuration parameters

**Configuration Options:**
```python
config = {
    'mhc_alleles': ['HLA-A*02:01', 'HLA-B*07:02'],  # HLA alleles to test
    'epitope_length_range': (8, 11),                # Peptide length range
    'prediction_threshold': 0.5,                    # Minimum score threshold
    'include_b_cell': True,                         # Include B-cell epitopes
    'include_t_cell': True,                         # Include T-cell epitopes
    'n_jobs': 4,                                    # Parallel processing cores
    'batch_size': 1000                              # Batch size for processing
}
```

**Returns:**
- List of epitope dictionaries

**Epitope Dictionary Structure:**
```python
epitope = {
    'protein_id': str,           # Source protein identifier
    'start_position': int,       # Start position (0-based)
    'end_position': int,         # End position (exclusive)
    'peptide': str,              # Epitope sequence
    'epitope_type': str,         # 'T-cell' or 'B-cell'
    'prediction_score': float,   # Prediction confidence (0-1)
    'mhc_allele': str,          # HLA allele (T-cell only)
    'length': int                # Epitope length
}
```

##### `predict_t_cell_epitopes(sequence: str, alleles: List[str]) -> List[Dict]`
Predict T-cell epitopes for specific HLA alleles.

##### `predict_b_cell_epitopes(sequence: str) -> List[Dict]`
Predict B-cell epitopes using BepiPred algorithm.

---

### Enhanced Epitope Prediction

Advanced epitope prediction with improved algorithms.

```python
from src.enhanced.epitope_prediction import EnhancedEpitopePrediction

enhanced_predictor = EnhancedEpitopePrediction()
```

#### Methods

##### `predict_epitopes_enhanced(protein_data: Dict) -> List[Dict]`
Enhanced epitope prediction with higher accuracy.

**Features:**
- Improved machine learning models
- Better population coverage analysis
- Structural information integration

---

### VaccineDesign

Design vaccine candidates from predicted epitopes.

```python
from src.vaccine_design import VaccineDesign

designer = VaccineDesign()
```

#### Methods

##### `design_vaccines(epitopes: List[Dict], config: Dict = None) -> List[Dict]`
Design multiple types of vaccine candidates.

**Parameters:**
- `epitopes` (List[Dict]): Epitopes from prediction
- `config` (Dict, optional): Design configuration

**Configuration Options:**
```python
config = {
    'vaccine_types': ['mRNA', 'subunit', 'peptide'],  # Types to generate
    'max_epitopes_per_vaccine': 20,                   # Maximum epitopes
    'linker_type': 'flexible',                        # 'flexible', 'rigid', 'turn'
    'optimize_for': 'global_coverage',                # 'efficacy', 'coverage', 'manufacturability'
    'include_adjuvant': True,                         # Add adjuvant sequences
    'mrna_optimization': {
        'codon_optimize': True,
        'gc_content_target': 0.5,
        'avoid_motifs': ['AAAA', 'TTTT']
    }
}
```

**Returns:**
- List of vaccine candidate dictionaries

**Vaccine Candidate Structure:**
```python
candidate = {
    'type': str,                    # 'mRNA', 'subunit', or 'peptide'
    'sequence': str,                # Vaccine sequence
    'epitopes_included': int,       # Number of epitopes
    'epitope_list': List[str],      # List of included epitopes
    'design_rationale': str,        # Design explanation
    'predicted_properties': Dict    # Predicted characteristics
}
```

##### `design_mrna_vaccine(epitopes: List[Dict]) -> Dict`
Design mRNA vaccine with codon optimization.

##### `design_subunit_vaccine(epitopes: List[Dict]) -> Dict`
Design subunit vaccine with linker optimization.

##### `design_peptide_vaccine(epitopes: List[Dict]) -> Dict`
Design peptide vaccine for synthesis.

---

### SafetyFilter

Screen vaccine candidates for safety concerns.

```python
from src.safety_filter import SafetyFilter

safety_filter = SafetyFilter()
```

#### Methods

##### `filter_candidates(candidates: List[Dict]) -> List[Dict]`
Filter vaccine candidates for safety.

**Safety Checks:**
- Allergenicity prediction
- Toxicity screening
- Human similarity analysis
- Autoimmune risk assessment

**Returns:**
- List of safe vaccine candidates with safety scores

##### `assess_allergenicity(sequence: str) -> Dict`
Assess allergenicity risk of a sequence.

##### `assess_toxicity(sequence: str) -> Dict`
Assess toxicity risk of a sequence.

##### `assess_human_similarity(sequence: str) -> Dict`
Check similarity to human proteins.

---

### VaccineRanking

Rank and score vaccine candidates.

```python
from src.ranking import VaccineRanking

ranker = VaccineRanking()
```

#### Methods

##### `rank_candidates(candidates: List[Dict]) -> List[Dict]`
Rank vaccine candidates by multiple criteria.

**Ranking Criteria:**
- Efficacy score (epitope quality)
- Population coverage
- Manufacturability
- Safety profile

**Returns:**
- List of candidates sorted by overall score

##### `calculate_efficacy_score(candidate: Dict) -> float`
Calculate predicted vaccine efficacy.

##### `calculate_population_coverage(candidate: Dict) -> float`
Calculate population coverage percentage.

##### `calculate_manufacturability_score(candidate: Dict) -> float`
Assess manufacturing feasibility.

---

## Advanced Modules

### NeoantigenIdentification

Identify patient-specific neoantigens for cancer vaccines.

```python
from src.advanced_modules.neoantigen_identification import NeoantigenIdentification

neoantigen_id = NeoantigenIdentification()
```

#### Methods

##### `identify_neoantigens(vcf_file: str, hla_alleles: List[str], expression_data: str = None) -> List[Dict]`
Identify neoantigens from genomic data.

**Parameters:**
- `vcf_file` (str): VCF file with somatic mutations
- `hla_alleles` (List[str]): Patient HLA alleles
- `expression_data` (str, optional): RNA-seq expression data

**Returns:**
- List of neoantigen candidates

##### `process_mutations(vcf_file: str) -> List[Dict]`
Process VCF file to extract mutations.

##### `generate_mutant_peptides(mutations: List[Dict]) -> List[str]`
Generate mutant peptide sequences.

---

### MHCClassIIPrediction

Predict MHC Class II binding for CD4+ T-cell responses.

```python
from src.advanced_modules.mhc_class_ii_prediction import MHCClassIIPrediction

mhc_ii_predictor = MHCClassIIPrediction()
```

#### Methods

##### `predict_mhc_ii_binding(peptides: List[str], alleles: List[str]) -> List[Dict]`
Predict MHC Class II binding affinity.

**Parameters:**
- `peptides` (List[str]): Peptide sequences (12-25 amino acids)
- `alleles` (List[str]): MHC Class II alleles

**Returns:**
- List of binding predictions with IC50 values

---

### ImmunogenicityPrediction

Predict actual immunogenicity beyond MHC binding.

```python
from src.advanced_modules.immunogenicity_prediction import ImmunogenicityPrediction

immuno_predictor = ImmunogenicityPrediction()
```

#### Methods

##### `predict_immunogenicity(epitopes: List[Dict]) -> List[Dict]`
Predict T-cell activation and antibody response.

**Features:**
- TCR binding prediction
- T-cell activation modeling
- Antibody response prediction

---

### ImmuneEvasionModeling

Model immune evasion and design resistant vaccines.

```python
from src.advanced_modules.immune_evasion_modeling import ImmuneEvasionModeling

evasion_model = ImmuneEvasionModeling()
```

#### Methods

##### `identify_conserved_epitopes(sequences: List[str]) -> List[Dict]`
Identify epitopes in conserved regions.

##### `predict_escape_mutations(epitopes: List[Dict]) -> List[Dict]`
Predict potential escape mutations.

##### `design_escape_resistant_vaccine(epitopes: List[Dict]) -> Dict`
Design vaccine resistant to immune evasion.

---

### ConformationalBCellPrediction

Predict conformational B-cell epitopes from 3D structure.

```python
from src.advanced_modules.conformational_bcell_prediction import ConformationalBCellPrediction

conf_predictor = ConformationalBCellPrediction()
```

#### Methods

##### `predict_conformational_epitopes(protein_id: str, structure_file: str = None) -> List[Dict]`
Predict conformational epitopes using AlphaFold structures.

---

### ExperimentalValidation

Integrate experimental validation data.

```python
from src.advanced_modules.experimental_validation import ExperimentalValidation

exp_validator = ExperimentalValidation()
```

#### Methods

##### `integrate_assay_results(epitopes: List[Dict], assay_data: Dict) -> List[Dict]`
Integrate experimental assay results.

##### `update_models_with_feedback(validation_data: Dict) -> None`
Update prediction models with experimental feedback.

---

## Visualization

### VisualizationModule

Generate plots and visualizations.

```python
from src.visualization import VisualizationModule

viz = VisualizationModule()
```

#### Methods

##### `plot_epitope_distribution(epitopes: List[Dict], protein_length: int) -> str`
Plot epitope positions along protein sequence.

##### `plot_vaccine_comparison(candidates: List[Dict]) -> str`
Compare vaccine candidates across multiple metrics.

##### `plot_population_coverage(coverage_data: Dict) -> str`
Visualize population coverage analysis.

##### `generate_html_report(results: Dict) -> str`
Generate comprehensive HTML report.

---

## Utilities

### Utils

Common utility functions.

```python
from src.utils import *
```

#### Functions

##### `calculate_sequence_similarity(seq1: str, seq2: str) -> float`
Calculate sequence similarity using BLOSUM62.

##### `translate_dna_to_protein(dna_sequence: str) -> str`
Translate DNA sequence to protein.

##### `reverse_translate_protein(protein_sequence: str) -> str`
Reverse translate protein to optimized DNA.

##### `calculate_molecular_weight(sequence: str) -> float`
Calculate molecular weight of protein sequence.

##### `predict_secondary_structure(sequence: str) -> str`
Predict secondary structure (alpha helix, beta sheet, coil).

---

## Configuration

### Global Configuration

Set global configuration parameters.

```python
from src.config import Config

# Set global configuration
Config.set_global_config({
    'default_mhc_alleles': ['HLA-A*02:01', 'HLA-B*07:02'],
    'cache_enabled': True,
    'parallel_processing': True,
    'n_jobs': 4,
    'device': 'cuda',  # 'cuda' or 'cpu'
    'batch_size': 1000
})

# Get current configuration
config = Config.get_global_config()
```

---

## Error Handling

### Common Exceptions

```python
from src.exceptions import *

try:
    epitopes = predictor.predict_epitopes(protein_data)
except InvalidSequenceError as e:
    print(f"Invalid sequence: {e}")
except PredictionError as e:
    print(f"Prediction failed: {e}")
except InsufficientDataError as e:
    print(f"Insufficient data: {e}")
```

### Exception Types

- `InvalidSequenceError`: Invalid protein sequence format
- `PredictionError`: Epitope prediction failure
- `InsufficientDataError`: Not enough data for analysis
- `ConfigurationError`: Invalid configuration parameters
- `ModelNotFoundError`: Required model files missing

---

## Performance Optimization

### Parallel Processing

```python
# Enable parallel processing
epitopes = predictor.predict_epitopes(
    protein_data,
    config={'n_jobs': 8, 'batch_size': 2000}
)
```

### Caching

```python
# Enable result caching
predictor.enable_caching(cache_size=10000)

# Clear cache
predictor.clear_cache()
```

### GPU Acceleration

```python
# Use GPU for deep learning models
config = {'device': 'cuda', 'batch_size': 4000}
epitopes = enhanced_predictor.predict_epitopes_enhanced(
    protein_data, 
    config=config
)
```

---

## Data Formats

### Input Formats

**FASTA Format:**
```
>spike_protein
MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFS
NVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIV
```

**VCF Format (for neoantigens):**
```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
chr1    12345   .       A       T       60      PASS    .
chr2    67890   .       G       C       80      PASS    .
```

### Output Formats

**JSON Format:**
```json
{
  "epitopes": [
    {
      "protein_id": "spike_protein",
      "start_position": 145,
      "end_position": 153,
      "peptide": "KIADYNYKL",
      "epitope_type": "T-cell",
      "prediction_score": 0.892,
      "mhc_allele": "HLA-A*02:01"
    }
  ],
  "vaccine_candidates": [
    {
      "type": "mRNA",
      "sequence": "AUGAAGAUCGCUGACUACAACUACAAGCUG...",
      "overall_score": 0.812
    }
  ]
}
```

**CSV Format:**
```csv
protein_id,start_position,end_position,peptide,epitope_type,prediction_score,mhc_allele
spike_protein,145,153,KIADYNYKL,T-cell,0.892,HLA-A*02:01
spike_protein,267,275,YQTSNFRVQ,T-cell,0.845,HLA-A*02:01
```

---

For more examples and detailed usage, see the [examples/](../examples/) directory and [User Guide](user_guide.md).

