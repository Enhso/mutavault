#!/usr/bin/env python3
"""
Biosecurity DNA Sequence Screening Tool
Analyzes a mystery DNA sequence through translation and BLAST searching
"""

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import time

# Mystery sequence to investigate
mystery_sequence = """ATGGTTCCTGGTCATATTGTTTCTGAAGCTAATGATCAACAACAAGAACAAGATGTTGGTGATGTTGTTTGTTTTGATATTTCTAGAGGTGCTGTTATTAGATGGAATGGTTCTCCTCCTGAAACTGCTGTTGGTATTGTTGCTTTGGATCATGGTCCTGCTCAAGATTTTGTTTCTGTTTATATTGGTGAACCTGTTAAATTGGATGAAGTTAAATGGGCTCAAGTTCCTGGTCCTGGTGCTTCTGAACCTGATTGTTTTTCTCATTTGAAAATTAGATTGGTT"""

def analyze_sequence(dna_sequence):
    """
    Perform biosecurity analysis on a DNA sequence
    """
    seq = Seq(dna_sequence.replace("\n", "").replace(" ", ""))
    
    print("="*70)
    print("BIOSECURITY SEQUENCE ANALYSIS")
    print("="*70)
    print(f"\nSequence length: {len(seq)} bp")
    print(f"GC content: {((seq.count('G') + seq.count('C')) / len(seq) * 100):.1f}%")
    
    # Step 1: Translate in all 6 reading frames
    print("\n" + "="*70)
    print("STEP 1: SIX-FRAME TRANSLATION ANALYSIS")
    print("="*70)
    
    frames_data = []
    
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            length = 3 * ((len(seq) - frame) // 3)
            trans = nuc[frame:frame + length].translate()
            
            # Find ORFs (stretches without stop codons)
            orfs = []
            current_orf = []
            for i, aa in enumerate(trans):
                if aa == '*':
                    if len(current_orf) > 20:  # Only report ORFs > 20 aa
                        orfs.append((len(current_orf), ''.join(current_orf)))
                    current_orf = []
                else:
                    current_orf.append(aa)
            
            if len(current_orf) > 20:
                orfs.append((len(current_orf), ''.join(current_orf)))
            
            frame_num = frame + 1 if strand == 1 else -(frame + 1)
            has_start = 'M' in str(trans[:50])
            stops = trans.count('*')
            
            frames_data.append({
                'frame': frame_num,
                'translation': str(trans),
                'has_start': has_start,
                'stops': stops,
                'orfs': orfs,
                'length': len(trans)
            })
            
            print(f"\nFrame {frame_num:+d}:")
            print(f"  Protein length: {len(trans)} aa")
            print(f"  Starts with M: {has_start}")
            print(f"  Stop codons: {stops}")
            if orfs:
                print(f"  ORFs found: {len(orfs)} (longest: {max(orfs, key=lambda x: x[0])[0]} aa)")
    
    # Step 2: Identify most likely coding frame
    print("\n" + "="*70)
    print("STEP 2: IDENTIFYING MOST LIKELY CODING FRAME")
    print("="*70)
    
    # Score frames: prefer those starting with M, few stops, long ORFs
    best_frame = None
    best_score = -1
    
    for frame_data in frames_data:
        score = 0
        if frame_data['has_start']:
            score += 100
        score -= frame_data['stops'] * 10
        if frame_data['orfs']:
            score += max(frame_data['orfs'], key=lambda x: x[0])[0]
        
        if score > best_score:
            best_score = score
            best_frame = frame_data
    
    print(f"\nBest candidate: Frame {best_frame['frame']:+d}")
    print(f"Score: {best_score}")
    print(f"\nTranslated protein sequence:")
    print(best_frame['translation'][:100] + "..." if len(best_frame['translation']) > 100 else best_frame['translation'])
    
    # Step 3: Analyze protein characteristics
    print("\n" + "="*70)
    print("STEP 3: PROTEIN CHARACTERISTICS")
    print("="*70)
    
    protein_seq = Seq(best_frame['translation'].replace('*', ''))
    
    print(f"\nProtein length: {len(protein_seq)} amino acids")
    print(f"Molecular weight: ~{len(protein_seq) * 110:.0f} Da (approximate)")
    
    # Amino acid composition analysis
    charged = protein_seq.count('R') + protein_seq.count('K') + protein_seq.count('D') + protein_seq.count('E')
    hydrophobic = protein_seq.count('A') + protein_seq.count('V') + protein_seq.count('L') + protein_seq.count('I')
    
    print(f"Charged residues: {charged} ({charged/len(protein_seq)*100:.1f}%)")
    print(f"Hydrophobic residues: {hydrophobic} ({hydrophobic/len(protein_seq)*100:.1f}%)")
    
    # Step 4: Perform BLAST searches
    print("\n" + "="*70)
    print("STEP 4: BLAST SEARCHES")
    print("="*70)
    
    # BLASTn search
    print("\n[4a] Nucleotide BLAST (BLASTn)")
    print("Searching NCBI nucleotide database (this may take 1-2 minutes)...")
    
    try:
        result_handle = NCBIWWW.qblast("blastn", "nt", str(seq), hitlist_size=5)
        blast_records = NCBIXML.parse(result_handle)
        
        print("\nTop BLASTn hits:")
        print("-" * 70)
        
        for i, blast_record in enumerate(blast_records):
            if blast_record.alignments:
                for j, alignment in enumerate(blast_record.alignments[:5]):
                    hsp = alignment.hsps[0]
                    identity = hsp.identities / hsp.align_length * 100
                    
                    print(f"\nHit {j+1}:")
                    print(f"  Sequence: {alignment.title[:80]}...")
                    print(f"  Identity: {identity:.1f}% ({hsp.identities}/{hsp.align_length})")
                    print(f"  E-value: {hsp.expect:.2e}")
                    print(f"  Score: {hsp.score}")
            else:
                print("\nNo significant matches found in nucleotide database")
            break
        
        result_handle.close()
        
    except Exception as e:
        print(f"\nBLASTn search failed: {e}")
        print("This could be due to network issues or NCBI server availability")
    
    # BLASTp search using translated protein
    print("\n" + "-" * 70)
    print("\n[4b] Protein BLAST (BLASTp)")
    print("Searching NCBI protein database with translated sequence...")
    print("(This may take 1-2 minutes)...")
    
    # Get clean protein sequence (remove stop codons)
    protein_query = best_frame['translation'].replace('*', '')
    
    try:
        result_handle = NCBIWWW.qblast("blastp", "nr", protein_query, hitlist_size=10)
        blast_records = NCBIXML.parse(result_handle)
        
        print("\nTop BLASTp hits:")
        print("-" * 70)
        
        hit_count = 0
        for i, blast_record in enumerate(blast_records):
            if blast_record.alignments:
                for j, alignment in enumerate(blast_record.alignments[:10]):
                    hsp = alignment.hsps[0]
                    identity = hsp.identities / hsp.align_length * 100
                    coverage = hsp.align_length / len(protein_query) * 100
                    
                    hit_count += 1
                    print(f"\nHit {hit_count}:")
                    print(f"  Protein: {alignment.title[:100]}")
                    if len(alignment.title) > 100:
                        print(f"           {alignment.title[100:200]}")
                    print(f"  Identity: {identity:.1f}% ({hsp.identities}/{hsp.align_length})")
                    print(f"  Coverage: {coverage:.1f}%")
                    print(f"  E-value: {hsp.expect:.2e}")
                    print(f"  Score: {hsp.score}")
                    
                    # Extract organism info if available
                    if '[' in alignment.title:
                        organism = alignment.title[alignment.title.find('[')+1:alignment.title.find(']')]
                        print(f"  Organism: {organism}")
            else:
                print("\nNo significant matches found in protein database")
            break
        
        result_handle.close()
        
    except Exception as e:
        print(f"\nBLASTp search failed: {e}")
        print("This could be due to network issues or NCBI server availability")
    
    # Step 5: Risk Assessment Summary
    print("\n" + "="*70)
    print("STEP 5: PRELIMINARY RISK ASSESSMENT")
    print("="*70)
    
    print("\nFactors to consider:")
    print("1. Sequence appears to encode a protein (has start codon, clear ORF)")
    print("2. Length suggests a small-medium protein domain")
    print("3. BLAST results would indicate organism of origin")
    print("4. Further analysis needed:")
    print("   - Cross-reference organism against select agent lists")
    print("   - Check protein function against known virulence factors")
    print("   - Analyze in context of full gene/genome if available")
    print("   - Consider intended use of synthesis")
    
    print("\n" + "="*70)
    print("RECOMMENDATION: Pending BLAST results and manual review")
    print("="*70)
    
    return best_frame

if __name__ == "__main__":
    analyze_sequence(mystery_sequence)