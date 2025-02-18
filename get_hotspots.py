#!/usr/bin/env python3
import argparse
import pandas as pd

def parse_protein_change(protein_change):
    """
    Given a string like 'E255K', return (WT_AA, ResidueIndex, Mut_AA).
    Example: 'E255K' -> ('E', 255, 'K')
    """
    pc = str(protein_change).strip()
    if len(pc) < 3:
        return None, None, None  # If format is off
    
    wt_aa = pc[0]
    mut_aa = pc[-1]

    # Attempt to extract the numeric portion
    # Typically the numeric portion is everything between the first and last character
    # e.g. E255K -> '255'
    residue_str = pc[1:-1]
    
    try:
        residue_index = int(residue_str)
    except ValueError:
        # In case the middle isn't purely digits
        residue_index = None

    return wt_aa, residue_index, mut_aa

def classify_tier(row):
    annotation = str(row['Annotation'])
    clinvar = str(row['ClinVar'])

    # Tier 4 criteria
    if ('level_1' in annotation or 'level_2' in annotation or
        'Pathogenic' in clinvar or 'FDA' in annotation):
        if ('CancerHotspot: yes' in annotation and
            'predictiveCount: ' in annotation):
            try:
                pred_count = int(annotation.split('predictiveCount: ')[1].split(',')[0])
                if pred_count >= 10:
                    return 4
            except:
                pass

    # Tier 3 criteria
    if ('level_R1' in annotation or 'level_R2' in annotation or
        'MyCancerGenome: present' in annotation):
        if 'predictiveCount: ' in annotation:
            try:
                pred_count = int(annotation.split('predictiveCount: ')[1].split(',')[0])
                if pred_count >= 5:
                    return 3
            except:
                pass

    # Tier 2 criteria
    if ('MutationAssessor: neutral' in annotation and 'SIFT: deleterious' in annotation or
        'predictiveCount: ' in annotation):
        try:
            pred_count = int(annotation.split('predictiveCount: ')[1].split(',')[0])
            if pred_count < 5:
                return 2
        except:
            pass

    # Tier 1 (benign)
    if 'Likely Neutral' in annotation or 'Benign' in clinvar:
        return 1

    # If none of the above
    return 0

def get_residue_groups(csv_path):
    # Read CSV
    df = pd.read_csv(csv_path)

    # Create new columns based on the "Protein Change" column
    df['WT AA'] = None
    df['Residue Index'] = None
    df['Mut AA'] = None

    # Parse out WT AA, Residue Index, and Mut AA
    for idx, row in df.iterrows():
        wt, residue_idx, mut = parse_protein_change(row['Protein Change'])
        df.at[idx, 'WT AA'] = wt
        df.at[idx, 'Residue Index'] = residue_idx
        df.at[idx, 'Mut AA'] = mut

    # Classify each row into tiers
    df['Tier'] = df.apply(classify_tier, axis=1)

    # Group by "Residue Index" and take the max tier
    # (Ensure residue index is numeric and drop rows where residue_index is None)
    df_filtered = df.dropna(subset=['Residue Index'])
    df_filtered['Residue Index'] = df_filtered['Residue Index'].astype(int)

    position_tiers = df_filtered.groupby('Residue Index')['Tier'].max()

    # Print residues by tier
    for tier in range(5):
        residues = position_tiers[position_tiers == tier].index.tolist()
        residues.sort()
        print(f"Tier {tier}: " + "+".join(str(r) for r in residues))

def main():
    parser = argparse.ArgumentParser(description="Parse and classify protein changes into tiers.")
    parser.add_argument("--csv", required=True, help="Path to the CSV file.")
    args = parser.parse_args()

    get_residue_groups(args.csv)

if __name__ == "__main__":
    main()

