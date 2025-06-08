# VaxGenAI Quick Start Guide

Get up and running with VaxGenAI in 5 minutes! This guide will walk you through the basic setup and your first vaccine design.

## 🚀 Installation Options

### Option 1: Docker (Recommended)
```bash
# Clone and run with Docker
git clone https://github.com/Valkyrie-Liberty/vaxgenai.git
cd vaxgenai
docker build -t vaxgenai .
docker run -p 8000:8000 vaxgenai
```

### Option 2: Local Installation
```bash
# Clone the repository
git clone https://github.com/Valkyrie-Liberty/vaxgenai.git
cd vaxgenai

# Create virtual environment
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

## 🧪 Your First Vaccine Design

### 1. Basic Epitope Prediction
```python
from src.epitope_prediction import EpitopePrediction
from src.input_module import InputModule

# Load a protein sequence
input_module = InputModule()
protein_data = input_module.load_fasta("data/spike_protein.fasta")

# Predict epitopes
predictor = EpitopePrediction()
epitopes = predictor.predict_epitopes(protein_data)

print(f"Found {len(epitopes)} potential epitopes")
for epitope in epitopes[:5]:  # Show first 5
    print(f"  {epitope['peptide']} (score: {epitope['prediction_score']:.3f})")
```

### 2. Design Vaccine Candidates
```python
from src.vaccine_design import VaccineDesign

# Design vaccines from epitopes
designer = VaccineDesign()
candidates = designer.design_vaccines(epitopes)

print(f"Generated {len(candidates)} vaccine candidates:")
for candidate in candidates:
    print(f"  {candidate['type']}: {len(candidate['sequence'])} residues")
```

### 3. Rank and Select Best Candidate
```python
from src.ranking import VaccineRanking

# Rank candidates
ranker = VaccineRanking()
ranked = ranker.rank_candidates(candidates)

best = ranked[0]
print(f"Best candidate: {best['type']}")
print(f"Overall score: {best['overall_score']:.3f}")
print(f"Efficacy: {best['efficacy_score']:.3f}")
print(f"Population coverage: {best['population_coverage']:.3f}")
```

## 🌐 Web Interface

### Start the Web Interface
```bash
python run_web.py
```

Then open http://localhost:8000 in your browser.

### Using the Web Interface
1. **Upload FASTA File**: Click "Choose File" and select your protein sequence
2. **Configure Parameters**: Adjust prediction settings if needed
3. **Run Analysis**: Click "Analyze Sequence"
4. **View Results**: Explore epitopes, vaccine candidates, and visualizations
5. **Download Results**: Save reports and vaccine sequences

## 📊 Understanding Results

### Epitope Prediction Output
```python
epitope = {
    'protein_id': 'spike_protein',
    'start_position': 145,
    'end_position': 153,
    'peptide': 'KIADYNYKL',
    'epitope_type': 'T-cell',
    'prediction_score': 0.892,
    'mhc_allele': 'HLA-A*02:01'
}
```

### Vaccine Candidate Output
```python
candidate = {
    'type': 'mRNA',
    'sequence': 'AUGAAGAUCGCUGACUACAACUACAAGCUG...',
    'epitopes_included': 15,
    'efficacy_score': 0.856,
    'population_coverage': 0.696,
    'manufacturability': 0.874,
    'overall_score': 0.812
}
```

## 🎯 Common Use Cases

### Cancer Neoantigen Vaccine
```python
from src.advanced_modules.neoantigen_identification import NeoantigenIdentification

# Identify patient-specific neoantigens
neoantigen_id = NeoantigenIdentification()
neoantigens = neoantigen_id.identify_neoantigens(
    vcf_file="patient_mutations.vcf",
    hla_alleles=["HLA-A*02:01", "HLA-B*07:02"]
)

# Design personalized vaccine
candidates = designer.design_vaccines(neoantigens)
```

### HIV Vaccine Design
```python
# Load HIV protein sequences
hiv_proteins = input_module.load_fasta("data/hiv_proteins.fasta")

# Focus on conserved regions
from src.advanced_modules.immune_evasion_modeling import ImmuneEvasionModeling
evasion_model = ImmuneEvasionModeling()
conserved_epitopes = evasion_model.identify_conserved_epitopes(hiv_proteins)

# Design vaccine targeting conserved epitopes
hiv_vaccine = designer.design_vaccines(conserved_epitopes)
```

### Pandemic Response
```python
# Rapid analysis of emerging pathogen
emerging_pathogen = input_module.load_fasta("data/novel_coronavirus.fasta")

# Quick epitope prediction
epitopes = predictor.predict_epitopes(emerging_pathogen, mode="fast")

# Generate emergency vaccine candidates
emergency_vaccines = designer.design_vaccines(epitopes, priority="speed")
```

## 🔧 Configuration

### Basic Configuration
```python
# Configure epitope prediction
config = {
    'mhc_alleles': ['HLA-A*02:01', 'HLA-B*07:02', 'HLA-C*07:01'],
    'epitope_length_range': (8, 11),
    'prediction_threshold': 0.5,
    'population': 'global'
}

epitopes = predictor.predict_epitopes(protein_data, config=config)
```

### Advanced Configuration
```python
# Configure vaccine design
vaccine_config = {
    'max_epitopes_per_vaccine': 20,
    'linker_type': 'flexible',
    'optimize_for': 'global_coverage',
    'include_adjuvant': True,
    'mrna_optimization': {
        'codon_optimize': True,
        'gc_content_target': 0.5,
        'avoid_motifs': ['AAAA', 'TTTT']
    }
}

candidates = designer.design_vaccines(epitopes, config=vaccine_config)
```

## 📈 Performance Tips

### Speed Optimization
```python
# Use parallel processing
epitopes = predictor.predict_epitopes(
    protein_data, 
    n_jobs=8,  # Use 8 CPU cores
    batch_size=1000
)

# Cache results for repeated analysis
predictor.enable_caching(cache_size=10000)
```

### Memory Optimization
```python
# Process large datasets in chunks
for chunk in input_module.load_fasta_chunks("large_dataset.fasta", chunk_size=100):
    epitopes = predictor.predict_epitopes(chunk)
    # Process epitopes immediately
    candidates = designer.design_vaccines(epitopes)
```

## 🐛 Troubleshooting

### Common Issues

**Issue**: `ModuleNotFoundError: No module named 'src'`
```bash
# Solution: Make sure you're in the vaxgenai directory
cd vaxgenai
export PYTHONPATH="${PYTHONPATH}:$(pwd)"
```

**Issue**: `CUDA out of memory`
```python
# Solution: Reduce batch size or use CPU
config = {'device': 'cpu', 'batch_size': 256}
epitopes = predictor.predict_epitopes(protein_data, config=config)
```

**Issue**: Web interface not accessible
```bash
# Solution: Check if port 8000 is available
netstat -an | grep 8000
# If occupied, use different port
python run_web.py --port 8080
```

### Getting Help

1. **Check Documentation**: See [docs/](docs/) for detailed guides
2. **Search Issues**: Check [GitHub Issues](https://github.com/Valkyrie-Liberty/vaxgenai/issues)
3. **Ask Questions**: Create a new [GitHub Discussion](https://github.com/Valkyrie-Liberty/vaxgenai/discussions)
4. **Report Bugs**: Create a new [GitHub Issue](https://github.com/Valkyrie-Liberty/vaxgenai/issues/new)

## 🎓 Next Steps

1. **Explore Examples**: Check out [examples/](examples/) for real-world use cases
2. **Read API Docs**: See [docs/api_reference.md](docs/api_reference.md) for complete API
3. **Join Community**: Connect with other researchers using VaxGenAI
4. **Contribute**: Help improve VaxGenAI by contributing code or documentation

## 📚 Additional Resources

- **[User Guide](user_guide.md)**: Comprehensive usage instructions
- **[API Reference](api_reference.md)**: Complete API documentation
- **[Developer Guide](developer_guide.md)**: Contributing and extending VaxGenAI
- **[Examples](../examples/)**: Real-world usage examples
- **[FAQ](faq.md)**: Frequently asked questions

---

**Ready to design your first vaccine? Let's get started!** 🚀

