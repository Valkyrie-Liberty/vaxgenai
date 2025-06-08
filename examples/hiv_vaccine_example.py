"""
Comprehensive example: Designing an HIV vaccine targeting conserved epitopes.

This example demonstrates how to use VaxGenAI to design an HIV vaccine that
targets conserved regions resistant to immune evasion.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.input_module import InputModule
from src.enhanced.epitope_prediction import EnhancedEpitopePrediction
from src.advanced_modules.immune_evasion_modeling import ImmuneEvasionModeling
from src.advanced_modules.conformational_bcell_prediction import ConformationalBCellPrediction
from src.advanced_modules.immunogenicity_prediction import ImmunogenicityPrediction
from src.vaccine_design import VaccineDesign
from src.safety_filter import SafetyFilter
from src.ranking import VaccineRanking
from src.visualization import VisualizationModule
import json

def main():
    """
    Complete workflow for HIV vaccine design targeting conserved epitopes.
    """
    print("🦠 VaxGenAI: HIV Vaccine Design - Conserved Epitope Targeting")
    print("=" * 65)
    
    # Step 1: Load HIV Protein Sequences
    print("\n📁 Step 1: Loading HIV Protein Sequences")
    
    input_module = InputModule()
    
    # Simulate HIV protein sequences (in real use, these would be from HIV databases)
    hiv_proteins = {
        'gag': {
            'id': 'HIV_gag',
            'sequence': 'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLYCVHQRIEIKDTKEALDKIEEEQNKSKKKAQQAAADTGHSNQVSQNYPIVQNIQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAAEWDRVHPVHAGPIAPGQMREPRGSDIAGTTSTLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPAATLEEMMTACQGVGGPGHKARVLAEAMSQVTNSATIMMQRGNFRNQRKIVKCFNCGKEGHTARNCRAPRKKGCWKCGKEGHQMKDCTERQANFLGKIWPSYKGRPGNFLQSRPEPTAPPEESFRSGVETTTPPQKQEPIDKELYPLTSLRSLFGNDPSSQ'
        },
        'pol': {
            'id': 'HIV_pol',
            'sequence': 'FFREDLAFLQGKAREFSSEQTRANSPTRRELQVWGRDNNSLSEAGADRQGTVSFSFPQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPISPIETVPVKLKPGMDGPKVKQWPLTEEKIKALVEICTEMEKEGKISKIGPENPYNTPVFAIKKKDSTKWRKLVDFRELNKRTQDFWEVQLGIPHPAGLKKKKSVTVLDVGDAYFSVPLDEDFRKYTAFTIPSINNETPGIRYQYNVLPQGWKGSPAIFQSSMTKILEPFRKQNPDIVIYQYMDDLYVGSDLEIGQHRTKIEELRQHLLRWGLTTPDKKHQKEPPFLWMGYELHPDKWTVQPIVLPEKDSWTVNDIQKLVGKLNWASQIYPGIKVRQLCKLLRGTKALTEVIPLTEEAELELAENREILKEPVHGVYYDPSKDLIAEIQKQGQGQWTYQIYQEPFKNLKTGKYARMRGAHTNDVKQLTEAVQKITTESIVIWGKTPKFKLPIQKETWETWWTEYWQATWIPEWEFVNTPPLVKLWYQLEKEPIVGAETFYVDGAANRETKLGKAGYVTNRGRQKVVTLTDTTNQKTELQAIHLALQDSGLEVNIVTDSQYALGIIQAQPDKSESELVSQIIEQLIKKEKVYLAWVPAHKGIGGNEQVDKLVSAGIRKVLFLDGIDKAQEEHEKYHSNWRAMASDFNLPPVVAKEIVASCDKCQLKGEAMHGQVDCSPGIWQLDCTHLEGKVILVAVHVASGYIEAEVIPAETGQETAYFLLKLAGRWPVKTIHTDNGSNFTSTTVKAACWWAGIKQEFGIPYNPQSQGVVESMNKELKKIIGQVRDQAEHLKTAVQMAVFIHNFKRKGGIGGYSAGERIVDIIATDIQTKELQKQITKIQNFRVYYRDSRNPLWKGPAKLLWKGEGAVVIQDNSDIKVVPRRKAKIIRDYGKQMAGDDCVASRQDED'
        },
        'env': {
            'id': 'HIV_env',
            'sequence': 'MRVKEKYQHLWRWGWRWGTMLLGMLMICSATEKLWVTVYYGVPVWKEATTTLFCASDAKAYDTEVHNVWATHACVPTDPNPQEVVLVNVTENFNMWKNDMVEQMHEDIISLWDQSLKPCVKLTPLCVSLKCTDLKNDTNTNSSSGRMIMEKGEIKNCSFNISTSIRGKVQKEYAFFYKLDIIPIDNDTTSYKLTSCNTSVITQACPKVSFEPIPIHYCAPAGFAILKCNNKTFNGTGPCTNVSTVQCTHGIRPVVSTQLLLNGSLAEEEVVIRSVNFTDNAKTIIVQLNTSVEINCTRPNNNTRKRIRIQRGPGRAFVTIGKIGNMRQAHCNISRAKWNNTLKQIASKLREQFGNNKTIIFKQSSGGDPEIVTHSFNCGGEFFYCNSTQLFNSTWFNSTWSTEGSNNTEGSDTITLPCRIKQIINMWQKVGKAMYAPPISGQIRCSSNITGLLLTRDGGNSNNESEIFRPGGGDMRDNWRSELYKYKVVKIEPLGVAPTKAKRRVVQREKRAVGIGALFLGFLGAAGSTMGAASMTLTVQARQLLSGIVQQQNNLLRAIEAQQHLLQLTVWGIKQLQARILAVERYLKDQQLLGIWGCSGKLICTTAVPWNASWSNKSLEQIWNHTTWMEWDREINNYTSLIHSLIEESQNQQEKNEQELLELDKWASLWNWFNITNWLWYIKLFIMIVGGLVGLRIVFAVLSIVNRVRQGYSPLSFQTHLPTPRGPDRPEGIEEEGGERDRDRSIRLVNGSLALIWDDLRSLCLFSYHRLRDLLLIVTRIVELLGRRGWEALKYWWNLLQYWSQELKNSAVSLLNATAIAVAEGTDRVIEVVQGACRAIRHIPRRIRQGLERILL'
        }
    }
    
    print(f"✅ Loaded {len(hiv_proteins)} HIV proteins:")
    for protein_id, data in hiv_proteins.items():
        print(f"   {protein_id.upper()}: {len(data['sequence'])} amino acids")
    
    # Step 2: Enhanced Epitope Prediction
    print("\n🎯 Step 2: Enhanced Epitope Prediction")
    
    enhanced_predictor = EnhancedEpitopePrediction()
    
    # Global HLA alleles for broad coverage
    global_hla_alleles = [
        'HLA-A*02:01', 'HLA-A*01:01', 'HLA-A*03:01', 'HLA-A*24:02', 'HLA-A*11:01',
        'HLA-B*07:02', 'HLA-B*08:01', 'HLA-B*15:01', 'HLA-B*35:01', 'HLA-B*44:02',
        'HLA-C*07:01', 'HLA-C*07:02', 'HLA-C*06:02', 'HLA-C*03:04', 'HLA-C*12:03'
    ]
    
    all_epitopes = []
    
    for protein_name, protein_data in hiv_proteins.items():
        print(f"\n   Analyzing {protein_name.upper()} protein...")
        
        epitopes = enhanced_predictor.predict_epitopes_enhanced(
            protein_data,
            config={
                'mhc_alleles': global_hla_alleles,
                'include_population_analysis': True,
                'structural_analysis': True
            }
        )
        
        # Add protein source information
        for epitope in epitopes:
            epitope['source_protein'] = protein_name
        
        all_epitopes.extend(epitopes)
        print(f"   ✅ Found {len(epitopes)} epitopes in {protein_name.upper()}")
    
    print(f"\n✅ Total epitopes identified: {len(all_epitopes)}")
    
    # Step 3: Immune Evasion Analysis
    print("\n🛡️ Step 3: Immune Evasion Analysis - Identifying Conserved Epitopes")
    
    evasion_model = ImmuneEvasionModeling()
    
    # Analyze sequence conservation across HIV strains
    conserved_epitopes = evasion_model.identify_conserved_epitopes(
        epitopes=all_epitopes,
        conservation_threshold=0.85  # 85% conservation across strains
    )
    
    print(f"✅ Identified {len(conserved_epitopes)} highly conserved epitopes")
    
    # Predict escape mutations
    escape_analysis = evasion_model.predict_escape_mutations(conserved_epitopes)
    
    print("🔍 Escape Mutation Analysis:")
    low_escape_epitopes = [e for e in escape_analysis if e['escape_probability'] < 0.1]
    print(f"   {len(low_escape_epitopes)} epitopes with <10% escape probability")
    
    # Display top conserved epitopes
    print("\n🎯 Top Conserved Epitopes (Escape-Resistant):")
    for i, epitope in enumerate(low_escape_epitopes[:5], 1):
        print(f"   {i}. {epitope['peptide']} ({epitope['source_protein'].upper()})")
        print(f"      Conservation: {epitope['conservation_score']:.1%}")
        print(f"      Escape Risk: {epitope['escape_probability']:.1%}")
        print(f"      HLA Coverage: {len(epitope['binding_alleles'])} alleles")
    
    # Step 4: Conformational B-cell Epitope Prediction
    print("\n🧬 Step 4: Conformational B-cell Epitope Prediction")
    
    conf_predictor = ConformationalBCellPrediction()
    
    # Predict conformational epitopes for broadly neutralizing antibodies
    conformational_epitopes = []
    
    for protein_name, protein_data in hiv_proteins.items():
        if protein_name == 'env':  # Focus on Env for antibody targets
            print(f"   Analyzing {protein_name.upper()} for conformational epitopes...")
            
            conf_epitopes = conf_predictor.predict_conformational_epitopes(
                protein_id=protein_data['id'],
                sequence=protein_data['sequence']
            )
            
            conformational_epitopes.extend(conf_epitopes)
            print(f"   ✅ Found {len(conf_epitopes)} conformational epitopes")
    
    print(f"✅ Total conformational epitopes: {len(conformational_epitopes)}")
    
    # Step 5: Immunogenicity Prediction
    print("\n⚡ Step 5: Immunogenicity Prediction")
    
    immuno_predictor = ImmunogenicityPrediction()
    
    # Combine T-cell and B-cell epitopes
    combined_epitopes = low_escape_epitopes + conformational_epitopes
    
    # Predict immunogenicity
    immunogenic_epitopes = immuno_predictor.predict_immunogenicity(combined_epitopes)
    
    print(f"✅ {len(immunogenic_epitopes)} epitopes predicted to be highly immunogenic")
    
    # Display top immunogenic epitopes
    print("\n🔥 Top Immunogenic Epitopes:")
    for epitope in immunogenic_epitopes[:5]:
        print(f"   {epitope['peptide']} ({epitope['source_protein'].upper()})")
        print(f"   Immunogenicity: {epitope['immunogenicity_score']:.3f}")
        print(f"   T-cell activation: {epitope.get('tcr_binding_score', 0):.3f}")
        print()
    
    # Step 6: HIV-Specific Vaccine Design
    print("\n💉 Step 6: HIV-Specific Vaccine Design")
    
    designer = VaccineDesign()
    
    # Configure vaccine design for HIV
    hiv_config = {
        'vaccine_types': ['mRNA', 'subunit'],  # Focus on mRNA and subunit for HIV
        'max_epitopes_per_vaccine': 25,  # More epitopes for HIV diversity
        'prioritize_conserved': True,
        'include_adjuvant': True,
        'optimize_for': 'broad_coverage',
        'hiv_specific_optimization': True
    }
    
    # Design vaccine candidates
    vaccine_candidates = designer.design_vaccines(immunogenic_epitopes, config=hiv_config)
    
    print(f"✅ Generated {len(vaccine_candidates)} HIV vaccine candidates")
    
    for candidate in vaccine_candidates:
        print(f"   {candidate['type']}: {candidate['epitopes_included']} epitopes")
        print(f"   Proteins covered: {len(set(e['source_protein'] for e in candidate['epitope_details']))}")
    
    # Step 7: Safety Assessment
    print("\n🛡️ Step 7: Safety Assessment")
    
    safety_filter = SafetyFilter()
    safe_candidates = safety_filter.filter_candidates(vaccine_candidates)
    
    print(f"✅ {len(safe_candidates)} candidates passed safety screening")
    
    # Step 8: Ranking with HIV-Specific Criteria
    print("\n🏆 Step 8: Ranking with HIV-Specific Criteria")
    
    ranker = VaccineRanking()
    
    # Custom ranking for HIV vaccines
    hiv_ranking_config = {
        'conservation_weight': 0.3,
        'escape_resistance_weight': 0.25,
        'immunogenicity_weight': 0.25,
        'population_coverage_weight': 0.2
    }
    
    ranked_candidates = ranker.rank_candidates(safe_candidates, config=hiv_ranking_config)
    
    print("📊 HIV Vaccine Rankings:")
    for i, candidate in enumerate(ranked_candidates, 1):
        print(f"   {i}. {candidate['type']} Vaccine")
        print(f"      Overall Score: {candidate['overall_score']:.3f}")
        print(f"      Conservation: {candidate.get('conservation_score', 0):.3f}")
        print(f"      Escape Resistance: {candidate.get('escape_resistance_score', 0):.3f}")
        print(f"      Immunogenicity: {candidate.get('immunogenicity_score', 0):.3f}")
        print(f"      Global Coverage: {candidate.get('population_coverage', 0):.3f}")
        print()
    
    # Step 9: HIV-Specific Analysis
    print("\n🔬 Step 9: HIV-Specific Analysis")
    
    best_candidate = ranked_candidates[0]
    
    # Analyze protein coverage
    protein_coverage = {}
    for epitope in best_candidate['epitope_details']:
        protein = epitope['source_protein']
        protein_coverage[protein] = protein_coverage.get(protein, 0) + 1
    
    print("🎯 Best Candidate Analysis:")
    print(f"   Type: {best_candidate['type']}")
    print(f"   Total Epitopes: {best_candidate['epitopes_included']}")
    print(f"   Protein Coverage:")
    for protein, count in protein_coverage.items():
        print(f"     {protein.upper()}: {count} epitopes")
    
    # Analyze HLA coverage
    hla_coverage = set()
    for epitope in best_candidate['epitope_details']:
        if 'binding_alleles' in epitope:
            hla_coverage.update(epitope['binding_alleles'])
    
    print(f"   HLA Allele Coverage: {len(hla_coverage)} alleles")
    print(f"   Estimated Global Coverage: {best_candidate.get('population_coverage', 0):.1%}")
    
    # Step 10: Generate Comprehensive Report
    print("\n📄 Step 10: Generating HIV Vaccine Report")
    
    # Prepare results
    results = {
        'analysis_type': 'HIV_vaccine_design',
        'proteins_analyzed': list(hiv_proteins.keys()),
        'total_epitopes': len(all_epitopes),
        'conserved_epitopes': len(conserved_epitopes),
        'escape_resistant_epitopes': len(low_escape_epitopes),
        'conformational_epitopes': len(conformational_epitopes),
        'immunogenic_epitopes': len(immunogenic_epitopes),
        'vaccine_candidates': ranked_candidates,
        'best_candidate': best_candidate,
        'protein_coverage': protein_coverage,
        'hla_coverage': list(hla_coverage)
    }
    
    # Generate visualizations
    viz = VisualizationModule()
    
    # Plot epitope conservation
    conservation_plot = viz.plot_epitope_conservation(conserved_epitopes)
    print(f"✅ Conservation analysis plot: {conservation_plot}")
    
    # Plot escape resistance
    escape_plot = viz.plot_escape_resistance(escape_analysis)
    print(f"✅ Escape resistance plot: {escape_plot}")
    
    # Plot vaccine comparison
    comparison_plot = viz.plot_vaccine_comparison(ranked_candidates)
    print(f"✅ Vaccine comparison plot: {comparison_plot}")
    
    # Generate HTML report
    html_report = viz.generate_hiv_vaccine_report(results)
    print(f"✅ Comprehensive HTML report: {html_report}")
    
    # Save results
    results_file = "results/hiv_vaccine_results.json"
    os.makedirs(os.path.dirname(results_file), exist_ok=True)
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"✅ Results saved to: {results_file}")
    
    # Step 11: Clinical Development Recommendations
    print("\n🏥 Step 11: Clinical Development Recommendations")
    
    print("💊 Recommended HIV Vaccine Strategy:")
    print(f"   Primary Candidate: {best_candidate['type']} vaccine")
    print(f"   Epitope Count: {best_candidate['epitopes_included']}")
    print(f"   Multi-protein targeting: {len(protein_coverage)} HIV proteins")
    
    print("\n📋 Clinical Protocol:")
    print("   1. Prime-boost vaccination regimen")
    print("   2. Combine with broadly neutralizing antibody therapy")
    print("   3. Monitor viral load and CD4+ T-cell counts")
    print("   4. Include diverse patient populations for HLA coverage")
    print("   5. Long-term follow-up for durability assessment")
    
    print("\n🔬 Development Priorities:")
    print("   1. In vitro validation with patient PBMCs")
    print("   2. Non-human primate challenge studies")
    print("   3. Phase I safety and immunogenicity trial")
    print("   4. Efficacy studies in high-risk populations")
    print("   5. Combination therapy optimization")
    
    print("\n🌍 Global Impact Potential:")
    print(f"   Estimated global population coverage: {best_candidate.get('population_coverage', 0):.1%}")
    print("   Targets conserved regions across HIV subtypes")
    print("   Designed for escape resistance")
    print("   Suitable for resource-limited settings")
    
    print("\n✨ HIV Vaccine Design Complete!")
    print("Ready for experimental validation and clinical development")

if __name__ == "__main__":
    main()

