"""
Comprehensive example: Designing a personalized neoantigen vaccine for cancer treatment.

This example demonstrates how to use VaxGenAI to identify patient-specific neoantigens
and design a personalized cancer vaccine.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.advanced_modules.neoantigen_identification import NeoantigenIdentification
from src.advanced_modules.mhc_class_ii_prediction import MHCClassIIPrediction
from src.advanced_modules.immunogenicity_prediction import ImmunogenicityPrediction
from src.vaccine_design import VaccineDesign
from src.safety_filter import SafetyFilter
from src.ranking import VaccineRanking
from src.visualization import VisualizationModule
import json

def main():
    """
    Complete workflow for personalized cancer neoantigen vaccine design.
    """
    print("🧬 VaxGenAI: Personalized Cancer Neoantigen Vaccine Design")
    print("=" * 60)
    
    # Step 1: Patient Information
    print("\n📋 Step 1: Patient Information")
    patient_info = {
        'patient_id': 'PATIENT_001',
        'cancer_type': 'Pancreatic Adenocarcinoma',
        'hla_alleles': [
            'HLA-A*02:01', 'HLA-A*24:02',
            'HLA-B*07:02', 'HLA-B*15:01',
            'HLA-C*07:01', 'HLA-C*03:04',
            'HLA-DRB1*03:01', 'HLA-DRB1*15:01',
            'HLA-DQB1*02:01', 'HLA-DQB1*06:02'
        ],
        'tumor_stage': 'Stage III',
        'prior_treatments': ['Chemotherapy', 'Radiation']
    }
    
    print(f"Patient ID: {patient_info['patient_id']}")
    print(f"Cancer Type: {patient_info['cancer_type']}")
    print(f"HLA Alleles: {len(patient_info['hla_alleles'])} alleles identified")
    
    # Step 2: Neoantigen Identification
    print("\n🔍 Step 2: Neoantigen Identification from Genomic Data")
    
    # Initialize neoantigen identification
    neoantigen_id = NeoantigenIdentification()
    
    # Simulate patient genomic data (in real use, this would be from WES/RNA-seq)
    sample_mutations = [
        {
            'gene': 'KRAS',
            'position': 12,
            'wild_type': 'G',
            'mutant': 'D',
            'expression_level': 8.5,
            'mutation_type': 'missense'
        },
        {
            'gene': 'TP53',
            'position': 273,
            'wild_type': 'R',
            'mutant': 'H',
            'expression_level': 6.2,
            'mutation_type': 'missense'
        },
        {
            'gene': 'SMAD4',
            'position': 361,
            'wild_type': 'R',
            'mutant': 'H',
            'expression_level': 4.8,
            'mutation_type': 'missense'
        },
        {
            'gene': 'CDKN2A',
            'position': 58,
            'wild_type': 'G',
            'mutant': 'A',
            'expression_level': 3.1,
            'mutation_type': 'missense'
        }
    ]
    
    # Identify neoantigens
    neoantigens = neoantigen_id.identify_neoantigens_from_mutations(
        mutations=sample_mutations,
        hla_alleles=patient_info['hla_alleles'][:6]  # MHC Class I alleles
    )
    
    print(f"✅ Identified {len(neoantigens)} potential neoantigens")
    
    # Display top neoantigens
    print("\n🎯 Top Neoantigen Candidates:")
    for i, neoantigen in enumerate(neoantigens[:5], 1):
        print(f"  {i}. {neoantigen['peptide']} (Gene: {neoantigen['gene']}, Score: {neoantigen['prediction_score']:.3f})")
    
    # Step 3: MHC Class II Prediction for CD4+ T-cell Response
    print("\n🧪 Step 3: MHC Class II Prediction for Helper T-cell Response")
    
    mhc_ii_predictor = MHCClassIIPrediction()
    
    # Generate longer peptides for MHC Class II (15-mers)
    long_peptides = neoantigen_id.generate_long_peptides(
        mutations=sample_mutations,
        peptide_length=15
    )
    
    # Predict MHC Class II binding
    mhc_ii_alleles = [allele for allele in patient_info['hla_alleles'] if 'DR' in allele or 'DQ' in allele]
    mhc_ii_epitopes = mhc_ii_predictor.predict_mhc_ii_binding(
        peptides=[p['peptide'] for p in long_peptides],
        alleles=mhc_ii_alleles
    )
    
    print(f"✅ Identified {len(mhc_ii_epitopes)} MHC Class II epitopes")
    print(f"   Top epitope: {mhc_ii_epitopes[0]['peptide']} (IC50: {mhc_ii_epitopes[0]['ic50']:.1f} nM)")
    
    # Step 4: Immunogenicity Prediction
    print("\n⚡ Step 4: Immunogenicity Prediction")
    
    immuno_predictor = ImmunogenicityPrediction()
    
    # Combine MHC Class I and II epitopes
    all_epitopes = neoantigens + mhc_ii_epitopes
    
    # Predict actual immunogenicity
    immunogenic_epitopes = immuno_predictor.predict_immunogenicity(all_epitopes)
    
    print(f"✅ {len(immunogenic_epitopes)} epitopes predicted to be immunogenic")
    
    # Display immunogenicity scores
    print("\n🔥 Immunogenicity Predictions:")
    for epitope in immunogenic_epitopes[:5]:
        print(f"  {epitope['peptide']}: {epitope['immunogenicity_score']:.3f} "
              f"(T-cell activation: {epitope['tcr_binding_score']:.3f})")
    
    # Step 5: Vaccine Design
    print("\n💉 Step 5: Personalized Vaccine Design")
    
    designer = VaccineDesign()
    
    # Configure vaccine design for cancer
    cancer_config = {
        'vaccine_types': ['mRNA', 'peptide'],  # Focus on mRNA and peptide for cancer
        'max_epitopes_per_vaccine': 15,
        'prioritize_neoantigens': True,
        'include_adjuvant': True,
        'optimize_for': 'immunogenicity'
    }
    
    # Design vaccine candidates
    vaccine_candidates = designer.design_vaccines(immunogenic_epitopes, config=cancer_config)
    
    print(f"✅ Generated {len(vaccine_candidates)} vaccine candidates")
    
    for candidate in vaccine_candidates:
        print(f"  {candidate['type']}: {candidate['epitopes_included']} epitopes, "
              f"{len(candidate['sequence'])} residues/nucleotides")
    
    # Step 6: Safety Filtering
    print("\n🛡️ Step 6: Safety Assessment")
    
    safety_filter = SafetyFilter()
    safe_candidates = safety_filter.filter_candidates(vaccine_candidates)
    
    print(f"✅ {len(safe_candidates)} candidates passed safety screening")
    
    # Display safety scores
    for candidate in safe_candidates:
        safety = candidate['safety_assessment']
        print(f"  {candidate['type']}: Allergenicity {safety['allergenicity_score']:.3f}, "
              f"Toxicity {safety['toxicity_score']:.3f}")
    
    # Step 7: Ranking and Selection
    print("\n🏆 Step 7: Vaccine Candidate Ranking")
    
    ranker = VaccineRanking()
    ranked_candidates = ranker.rank_candidates(safe_candidates)
    
    print("📊 Final Rankings:")
    for i, candidate in enumerate(ranked_candidates, 1):
        print(f"  {i}. {candidate['type']} Vaccine")
        print(f"     Overall Score: {candidate['overall_score']:.3f}")
        print(f"     Efficacy: {candidate['efficacy_score']:.3f}")
        print(f"     Immunogenicity: {candidate.get('immunogenicity_score', 0):.3f}")
        print(f"     Safety: {candidate['safety_score']:.3f}")
        print(f"     Manufacturability: {candidate['manufacturability_score']:.3f}")
        print()
    
    # Step 8: Generate Report
    print("\n📄 Step 8: Generating Comprehensive Report")
    
    # Prepare results for visualization
    results = {
        'patient_info': patient_info,
        'neoantigens': neoantigens,
        'mhc_ii_epitopes': mhc_ii_epitopes,
        'immunogenic_epitopes': immunogenic_epitopes,
        'vaccine_candidates': ranked_candidates,
        'analysis_summary': {
            'total_mutations': len(sample_mutations),
            'neoantigens_identified': len(neoantigens),
            'immunogenic_epitopes': len(immunogenic_epitopes),
            'final_candidates': len(ranked_candidates)
        }
    }
    
    # Generate visualizations
    viz = VisualizationModule()
    
    # Plot neoantigen distribution
    neoantigen_plot = viz.plot_neoantigen_distribution(neoantigens)
    print(f"✅ Neoantigen distribution plot saved: {neoantigen_plot}")
    
    # Plot vaccine comparison
    comparison_plot = viz.plot_vaccine_comparison(ranked_candidates)
    print(f"✅ Vaccine comparison plot saved: {comparison_plot}")
    
    # Generate HTML report
    html_report = viz.generate_cancer_vaccine_report(results)
    print(f"✅ Comprehensive HTML report generated: {html_report}")
    
    # Save results to JSON
    results_file = "results/cancer_vaccine_results.json"
    os.makedirs(os.path.dirname(results_file), exist_ok=True)
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"✅ Results saved to: {results_file}")
    
    # Step 9: Clinical Recommendations
    print("\n🏥 Step 9: Clinical Recommendations")
    
    best_candidate = ranked_candidates[0]
    
    print("💊 Recommended Vaccine:")
    print(f"   Type: {best_candidate['type']}")
    print(f"   Epitopes: {best_candidate['epitopes_included']}")
    print(f"   Overall Score: {best_candidate['overall_score']:.3f}")
    
    print("\n📋 Clinical Protocol Recommendations:")
    print("   1. Prime-boost vaccination strategy")
    print("   2. Combine with checkpoint inhibitor therapy")
    print("   3. Monitor immune response via ELISpot assays")
    print("   4. Consider adjuvant to enhance immunogenicity")
    print("   5. Schedule follow-up at 2, 4, 8, and 12 weeks")
    
    print("\n🔬 Next Steps:")
    print("   1. Synthesize vaccine candidates")
    print("   2. In vitro immunogenicity testing")
    print("   3. Preclinical validation in mouse models")
    print("   4. IND filing for clinical trial")
    print("   5. Phase I safety and immunogenicity study")
    
    print("\n✨ Analysis Complete!")
    print(f"Personalized vaccine designed for {patient_info['patient_id']}")
    print(f"Ready for experimental validation and clinical development")

if __name__ == "__main__":
    main()

