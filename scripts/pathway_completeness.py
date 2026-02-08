#!/usr/bin/env python3
"""
KEGG Pathway Completeness Analyzer

Calculate completeness of key metabolic pathways using KEGG KO annotations.
"""

import sys
from pathlib import Path
from collections import defaultdict

sys.path.insert(0, str(Path(__file__).parent.parent))
from bennu.operators import Bennu

# Key metabolic pathways with their essential KOs
# Format: pathway_name -> list of (step_name, [required_KOs])
# A step is complete if ANY of its KOs are present

PATHWAYS = {
    "Wood-Ljungdahl (Acetyl-CoA pathway)": [
        ("Formate dehydrogenase", ["K00122", "K00123", "K00124", "K00125", "K00126", "K05299", "K22015", "K22516"]),
        ("Formyl-H4F synthetase", ["K01938"]),
        ("Methenyl-H4F cyclohydrolase", ["K01491"]),
        ("Methylene-H4F dehydrogenase", ["K01491", "K00297"]),
        ("Methylene-H4F reductase", ["K00297"]),
        ("Methyltransferase (to CFeSP)", ["K15023", "K14138"]),
        ("CO dehydrogenase/acetyl-CoA synthase", ["K00192", "K00193", "K00194", "K00197", "K00198", "K03520"]),
    ],

    "Methanogenesis (from CO2)": [
        ("Formylmethanofuran dehydrogenase", ["K00200", "K00201", "K00202", "K00203", "K00204", "K00205"]),
        ("Formylmethanofuran:H4MPT formyltransferase", ["K00672"]),
        ("Methenyl-H4MPT cyclohydrolase", ["K01499"]),
        ("Methylene-H4MPT dehydrogenase", ["K00319", "K00320"]),
        ("Methylene-H4MPT reductase", ["K00577"]),
        ("Methyl-H4MPT:CoM methyltransferase", ["K00577", "K00578", "K00579", "K00580", "K00581", "K00584"]),
        ("Methyl-CoM reductase (MCR)", ["K00399", "K00401", "K00402"]),  # THE key enzyme
        ("Heterodisulfide reductase", ["K03388", "K03389", "K03390", "K08264", "K08265"]),
    ],

    "Calvin Cycle (CO2 fixation)": [
        ("RuBisCO", ["K01601", "K01602"]),
        ("Phosphoglycerate kinase", ["K00927"]),
        ("GAPDH", ["K00134", "K00150"]),
        ("Triose-P isomerase", ["K01803"]),
        ("Aldolase", ["K01623", "K01624", "K11645"]),
        ("FBPase", ["K03841", "K02446", "K01086"]),
        ("Transketolase", ["K00615"]),
        ("SBPase", ["K01086", "K11532"]),
        ("Phosphoribulokinase (PRK)", ["K00855"]),  # THE key enzyme
        ("Ribose-5-P isomerase", ["K01807", "K01808"]),
        ("Ribulose-5-P epimerase", ["K01783"]),
    ],

    "TCA Cycle": [
        ("Citrate synthase", ["K01647", "K05942"]),
        ("Aconitase", ["K01681", "K01682"]),
        ("Isocitrate dehydrogenase", ["K00031", "K00030"]),
        ("2-Oxoglutarate dehydrogenase / 2-Oxoglutarate:ferredoxin oxidoreductase", ["K00174", "K00175", "K00164", "K00658", "K00382"]),
        ("Succinyl-CoA synthetase", ["K01899", "K01900", "K01902", "K01903"]),
        ("Succinate dehydrogenase", ["K00234", "K00235", "K00236", "K00237"]),
        ("Fumarase", ["K01676", "K01677", "K01678"]),
        ("Malate dehydrogenase", ["K00024", "K00025", "K00026"]),
    ],

    "Glycolysis/Gluconeogenesis": [
        ("Hexokinase/Glucokinase", ["K00844", "K00845", "K12407", "K00886"]),
        ("Glucose-6-P isomerase", ["K01810"]),
        ("Phosphofructokinase", ["K00850", "K16370", "K21071"]),
        ("Aldolase", ["K01623", "K01624", "K11645"]),
        ("Triose-P isomerase", ["K01803"]),
        ("GAPDH", ["K00134", "K00150"]),
        ("Phosphoglycerate kinase", ["K00927"]),
        ("Phosphoglycerate mutase", ["K01834", "K15633", "K15634", "K15635"]),
        ("Enolase", ["K01689"]),
        ("Pyruvate kinase", ["K00873"]),
    ],

    "Nitrogen Fixation": [
        ("Nitrogenase Fe protein (NifH)", ["K02588"]),
        ("Nitrogenase MoFe protein alpha (NifD)", ["K02586"]),
        ("Nitrogenase MoFe protein beta (NifK)", ["K02591"]),
    ],

    "Sulfate Reduction": [
        ("Sulfate adenylyltransferase", ["K00958", "K00957"]),
        ("Adenylylsulfate reductase", ["K00394", "K00395"]),
        ("Dissimilatory sulfite reductase", ["K11180", "K11181"]),
    ],

    "Hydrogenase (NiFe)": [
        ("NiFe-hydrogenase large subunit", ["K00437", "K06281", "K06282", "K18005", "K18006", "K18007", "K18008"]),
        ("NiFe-hydrogenase small subunit", ["K06281", "K03620", "K18017"]),
        ("Hydrogenase maturation (HypA-F)", ["K04651", "K04652", "K04653", "K04654", "K04655", "K04656"]),
    ],

    "Electron Transport (Archaeal)": [
        ("NADH dehydrogenase", ["K00330", "K00331", "K00332", "K00333", "K00334", "K00337", "K00338", "K00339", "K00340", "K00341", "K00342", "K00343"]),
        ("Ferredoxin:NADP oxidoreductase", ["K00384", "K03852"]),
        ("ATP synthase", ["K02117", "K02118", "K02119", "K02120", "K02121", "K02122", "K02123", "K02124", "K02125", "K02126"]),
    ],

    "CRISPR-Cas (Type I)": [
        ("Cas1", ["K15342"]),
        ("Cas2", ["K09951"]),
        ("Cas3", ["K07012"]),
        ("Cas5", ["K07013"]),
        ("Cas6", ["K19091", "K19122"]),
        ("Cas7", ["K19125", "K07014"]),
        ("Cas8", ["K07015", "K19126", "K19127"]),
    ],

    "CRISPR-Cas (Type III)": [
        ("Cas1", ["K15342"]),
        ("Cas2", ["K09951"]),
        ("Cas6", ["K19091"]),
        ("Cas10/Csm1/Cmr2", ["K19130", "K19137"]),
        ("Csm/Cmr proteins", ["K19131", "K19132", "K19133", "K19134", "K19138", "K19139", "K19140"]),
    ],
}


def get_ko_presence(b, genome=None):
    """Get set of KOs present in genome (or all genomes)."""
    if genome:
        query = f"""
            SELECT DISTINCT a.accession
            FROM annotations a
            JOIN proteins p ON a.protein_id = p.protein_id
            WHERE a.source = 'kegg' AND p.bin_id = '{genome}'
        """
    else:
        query = """
            SELECT DISTINCT accession
            FROM annotations
            WHERE source = 'kegg'
        """

    results = b.store.execute(query)
    return set(row[0] for row in results)


def calculate_pathway_completeness(ko_set, pathway_steps):
    """Calculate completeness of a pathway given KO presence."""
    completed_steps = []
    missing_steps = []

    for step_name, required_kos in pathway_steps:
        if any(ko in ko_set for ko in required_kos):
            completed_steps.append(step_name)
        else:
            missing_steps.append(step_name)

    completeness = len(completed_steps) / len(pathway_steps) * 100
    return completeness, completed_steps, missing_steps


def analyze_pathways(db_path, output_file=None):
    """Analyze pathway completeness across all genomes."""
    b = Bennu(db_path)

    # Get all genomes
    genomes = [row[0] for row in b.store.execute("SELECT DISTINCT bin_id FROM proteins")]

    print(f"Analyzing {len(genomes)} genomes for pathway completeness...")
    print()

    results = []

    # Dataset-wide analysis first
    all_kos = get_ko_presence(b)
    print("=" * 70)
    print("DATASET-WIDE PATHWAY COMPLETENESS")
    print("=" * 70)

    for pathway_name, steps in PATHWAYS.items():
        completeness, completed, missing = calculate_pathway_completeness(all_kos, steps)
        status = "COMPLETE" if completeness == 100 else "PARTIAL" if completeness > 50 else "INCOMPLETE"
        print(f"\n{pathway_name}: {completeness:.0f}% ({len(completed)}/{len(steps)} steps) [{status}]")

        if missing and completeness < 100:
            print(f"  Missing: {', '.join(missing[:3])}" + ("..." if len(missing) > 3 else ""))

        results.append({
            'genome': 'ALL',
            'pathway': pathway_name,
            'completeness': completeness,
            'completed_steps': len(completed),
            'total_steps': len(steps),
            'missing': missing
        })

    # Per-genome analysis
    print("\n" + "=" * 70)
    print("PER-GENOME PATHWAY PRESENCE (>50% complete)")
    print("=" * 70)

    pathway_genome_counts = defaultdict(int)

    for genome in genomes:
        genome_kos = get_ko_presence(b, genome)

        for pathway_name, steps in PATHWAYS.items():
            completeness, completed, missing = calculate_pathway_completeness(genome_kos, steps)

            if completeness > 50:
                pathway_genome_counts[pathway_name] += 1

            results.append({
                'genome': genome,
                'pathway': pathway_name,
                'completeness': completeness,
                'completed_steps': len(completed),
                'total_steps': len(steps),
                'missing': missing
            })

    print()
    for pathway_name in PATHWAYS:
        n_genomes = pathway_genome_counts[pathway_name]
        pct = n_genomes / len(genomes) * 100
        print(f"{pathway_name}: {n_genomes}/{len(genomes)} genomes ({pct:.0f}%)")

    # Save results if output file specified
    if output_file:
        import json
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\nResults saved to {output_file}")

    return results


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Analyze KEGG pathway completeness")
    parser.add_argument("db_path", help="Path to Bennu database")
    parser.add_argument("-o", "--output", help="Output JSON file")
    args = parser.parse_args()

    analyze_pathways(args.db_path, args.output)
