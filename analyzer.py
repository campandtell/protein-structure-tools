#!/usr/bin/env python3
import requests
import pandas as pd
import os
from Bio import PDB
from Bio.PDB.Polypeptide import is_aa
import io
from concurrent.futures import ThreadPoolExecutor, as_completed

def analyze_genes(genes: list, output_dir: str = "results", max_workers: int = 40) -> dict:
    """
    Analyze structures for given gene symbols in parallel.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    results = {}
    
    def analyze_single_structure(pdb_id: str, gene: str, uniprot_id: str):
        """Worker function to analyze one structure"""
        try:
            info_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
            info = requests.get(info_url).json()
            
            method = info['rcsb_entry_info'].get('experimental_method')
            resolution = info['rcsb_entry_info'].get('resolution_combined', [None])[0]
            
            # Get gaps
            pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            response = requests.get(pdb_url)
            parser = PDB.PDBParser(QUIET=True)
            structure = parser.get_structure('protein', io.StringIO(response.text))
            
            gaps = {}
            gap_count = 0
            for chain in structure[0]:
                residues = [r.get_id()[1] for r in chain if is_aa(r) and 'CA' in r]
                if residues:
                    chain_gaps = []
                    for i in range(len(residues)-1):
                        if residues[i+1] - residues[i] > 1:
                            chain_gaps.append(f"{residues[i]}-{residues[i+1]}")
                            gap_count += 1
                    if chain_gaps:
                        gaps[chain.get_id()] = chain_gaps
            
            gap_str = "; ".join([f"Chain {chain}: {', '.join(gaps)}" for chain, gaps in gaps.items()])
            
            # Calculate quality score
            quality_score = 0
            if method == "X-RAY DIFFRACTION" and resolution:
                quality_score = 100 - (resolution * 20)
            elif method == "ELECTRON MICROSCOPY" and resolution:
                quality_score = 80 - (resolution * 15)
            else:
                quality_score = 30
                
            if gap_count > 0:
                quality_score -= (gap_count * 5)
            
            return {
                'gene': gene,
                'uniprot_id': uniprot_id,
                'pdb_id': pdb_id,
                'method': method,
                'resolution': resolution,
                'gap_count': gap_count,
                'gaps': gap_str if gaps else "No gaps",
                'quality_score': max(0, quality_score)
            }
        except Exception as e:
            print(f"Error analyzing {pdb_id}: {str(e)}")
            return None

    def analyze_one_gene(gene):
        """Process a single gene"""
        print(f"\nAnalyzing {gene}...")
        
        # Get UniProt ID
        uniprot_url = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{gene}+AND+organism_id:9606+reviewed:true&format=json"
        response = requests.get(uniprot_url)
        data = response.json()
        
        if not data.get('results'):
            print(f"No human UniProt entry found for {gene}")
            return gene, []
            
        uniprot_id = data['results'][0]['primaryAccession']
        print(f"Found human UniProt ID: {uniprot_id}")
        
        # Get PDB IDs
        entry_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}?format=json"
        response = requests.get(entry_url)
        data = response.json()
        
        pdb_ids = []
        for xref in data.get('uniProtKBCrossReferences', []):
            if xref.get('database') == 'PDB':
                pdb_ids.append(xref['id'])
        
        print(f"Found {len(pdb_ids)} PDB structures")
        
        # Analyze structures in parallel
        gene_results = []
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [
                executor.submit(analyze_single_structure, pdb_id, gene, uniprot_id)
                for pdb_id in pdb_ids
            ]
            
            for idx, future in enumerate(as_completed(futures), 1):
                result = future.result()
                if result:
                    gene_results.append(result)
                    print(f"Analyzed structure {idx}/{len(pdb_ids)} - Score: {result['quality_score']:.1f}")
        
        return gene, gene_results

    # Process genes in parallel at a higher level
    with ThreadPoolExecutor(max_workers=3) as executor:
        gene_futures = [executor.submit(analyze_one_gene, gene) for gene in genes]
        
        for future in as_completed(gene_futures):
            gene, gene_results = future.result()
            if gene_results:
                df = pd.DataFrame(gene_results)
                df = df.sort_values('quality_score', ascending=False)
                output_file = os.path.join(output_dir, f"{gene.lower()}_structures.csv")
                df.to_csv(output_file, index=False)
                print(f"\nResults saved to {output_file}")
                results[gene] = df
            
    return results

def download_structures(csv_files: list, min_quality: float = 50, output_dir: str = "pdbs", max_workers: int = 5):
    """
    Download high-quality structures in parallel.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    def download_single_structure(pdb_id: str, score: float):
        """Worker function to download one structure"""
        try:
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            response = requests.get(url)
            response.raise_for_status()
            
            output_file = os.path.join(output_dir, f"{pdb_id}.pdb")
            with open(output_file, 'wb') as f:
                f.write(response.content)
            print(f"Downloaded {pdb_id} (Score: {score:.1f})")
            return True
        except Exception as e:
            print(f"Error downloading {pdb_id}: {str(e)}")
            return False
    
    total_downloaded = 0
    for csv_file in csv_files:
        gene = os.path.splitext(os.path.basename(csv_file))[0].replace('_structures', '')
        print(f"\nProcessing structures for {gene}...")
        
        df = pd.read_csv(csv_file)
        good_structures = df[df['quality_score'] >= min_quality]
        
        if good_structures.empty:
            print(f"No structures found with quality score >= {min_quality}")
            continue
            
        print(f"Downloading {len(good_structures)} structures with quality score >= {min_quality}")
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [
                executor.submit(download_single_structure, row['pdb_id'], row['quality_score'])
                for _, row in good_structures.iterrows()
            ]
            
            for future in as_completed(futures):
                if future.result():
                    total_downloaded += 1
    
    print(f"\nDownloaded {total_downloaded} total structures to {output_dir}/")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Analyze and download protein structures")
    subparsers = parser.add_subparsers(dest='command')
    
    # Analyze command
    analyze_parser = subparsers.add_parser('analyze')
    analyze_parser.add_argument('genes', nargs='+', help='Gene symbols to analyze')
    analyze_parser.add_argument('--output-dir', default='results', 
                              help='Directory for CSV files')
    
    # Download command
    download_parser = subparsers.add_parser('download')
    download_parser.add_argument('csv_files', nargs='+', 
                               help='CSV files from analyze command')
    download_parser.add_argument('--min-quality', type=float, default=50, 
                               help='Minimum quality score (0-100)')
    download_parser.add_argument('--output-dir', default='pdbs',
                               help='Directory for PDB files')
    
    args = parser.parse_args()
    
    if args.command == 'analyze':
        analyze_genes(args.genes, args.output_dir)
    elif args.command == 'download':
        download_structures(args.csv_files, args.min_quality, args.output_dir)
    else:
        parser.print_help()
