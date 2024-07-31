import argparse
from collections import defaultdict
import json
import os
import subprocess
from tempfile import NamedTemporaryFile
from ..utils import parseFasta, writeFasta
from Bio.Seq import Seq
from xxhash import xxh3_64
from multiprocessing import Pool


def dedupe(sequences):
    dupes = set()
    out_sequences = []
    for header, sequence in sequences:
        if sequence in dupes:
            continue
        dupes.add(sequence)
        out_sequences.append((header, sequence))
    return out_sequences

def worker(
    fasta_file,
    out_path,
    sequence_data,
    align_method,
):
    print(fasta_file)
    og_hashes = {}
    og_pubs = {}
    out_sequences = []
    sequences = []
    replacements = []
    insertions = []
    log = []
    for i, (header, seq) in enumerate(parseFasta(fasta_file, True)):
        og_hashes[xxh3_64(seq).hexdigest()] = header
        
        description, data = header.split(" ", 1)
        json_data = json.loads(data)
        pub_gene_id = json_data["pub_gene_id"]
        
        og_pubs[pub_gene_id] = header, seq, i
        
        out_sequences.append((header, seq))
    
    
    for json_data in sequence_data:
        pub_gene_id = json_data["pub_gene_id"]
        description = json_data["description"]
        
        if "isoform" in description:
            sequences.append(("isoform_///"+pub_gene_id+"///"+description, str(Seq(json_data["genome_sequence"]).translate()).strip("*")))
            # pub_to_transcripts[pub_gene_id].append((description, json_data["genome_sequence"]))
            
    with NamedTemporaryFile(mode="w") as raw_fa_file, NamedTemporaryFile(mode="w") as aln_file:
        
        short_to_full = {}
        for header, _ in out_sequences+sequences:
            short_to_full[header[:128]] = header
        
        writeFasta(raw_fa_file.name, out_sequences+sequences)
        if align_method == "clustal":
            subprocess.run(
                ["clustalo", "-i", raw_fa_file.name, "-o", aln_file.name, "--threads=1", "--force"]
            )
        else:
            subprocess.run(
                ["mafft", "--thread", "1", "--quiet", "--anysymbol", raw_fa_file.name],
                stdout=aln_file
            )
        
        aligned_sequences = [(short_to_full[header[:128]], sequence) for header, sequence in parseFasta(aln_file.name, True)]

    this_consensus = defaultdict(list)
    for header, sequence in aligned_sequences:
        if "isoform" not in header:
            for i, aa in enumerate(sequence):
                this_consensus[i].append(aa)

    pub_to_transcripts = defaultdict(list)
    for header, sequence in aligned_sequences:
        if "isoform_" in header:
            _, pub_gene_id, description = header.split("///")
            pub_to_transcripts[pub_gene_id].append((description, sequence))
  
    for pub_gene_id, aligned_sequences in pub_to_transcripts.items():
        # Dedupe
        aligned_sequences = dedupe(aligned_sequences)
                
        highest_scoring_index = None
        highest_score = -1
        
        for i, (_, sequence) in enumerate(aligned_sequences):
            score = sum(this_consensus[j].count(sequence[j]) for j in range(len(sequence)))
            if score > highest_score:
                highest_score = score
                highest_scoring_index = i
                
        best_head, best_seq = aligned_sequences[highest_scoring_index]
        best_seq = best_seq.replace("-", "")
        
        if pub_gene_id in og_pubs:
            if best_seq != og_pubs[pub_gene_id][1]:
                i = og_pubs[pub_gene_id][2]
                # insertions.append((i, (best_head, best_seq), pub_gene_id))
                replacements.append((i, (best_head, best_seq), pub_gene_id))
                log.append("Replaced: "+og_pubs[pub_gene_id][0]+" with "+best_head)
        else:
            insertions.append((0, (best_head, best_seq), pub_gene_id))
            log.append("Inserted: "+best_head)
    
    for i, (header, seq), pub_gene_id in replacements:
        out_sequences[i] = (og_pubs[pub_gene_id][0], seq)
    
    insertions.sort(key=lambda x: x[0])
    inserted = 0
    for i, (header, seq), pub_gene_id in insertions:
        out_sequences.insert(i + inserted, (header+" - "+og_pubs[pub_gene_id][0], seq))
        inserted += 1
        
    with NamedTemporaryFile(mode="w") as raw_fa_file, NamedTemporaryFile(mode="w") as aln_file:
        writeFasta(raw_fa_file.name, out_sequences)
        if align_method == "clustal":
            subprocess.run(
                ["clustalo", "-i", raw_fa_file.name, "-o", aln_file.name, "--threads=1", "--force"]
            )
        else:
            subprocess.run(
                ["mafft", "--thread", "1", "--quiet", "--anysymbol", raw_fa_file.name],
                stdout=aln_file
            )
            
        writeFasta(out_path, parseFasta(aln_file.name, True))

    return log

def main(args):
    print("Begin")
    folder = args.INPUT
    output = args.OUTPUT
    tsv_path = args.tsv
    
    if not os.path.exists(output):
        os.makedirs(output)
    
    fasta_files = os.listdir(folder)
    print("Grabbing TSV file")
    data = defaultdict(list)
    with open(tsv_path, "r") as f:
        for line in f:
            #pub_og_id,organism_name,pub_gene_id,description,assembly_accession,genome_sequence
            pub_og_id,organism_name,pub_gene_id,description,assembly_accession,genome_sequence = line.strip().split("\t")
            data[pub_og_id].append({
                "pub_og_id": pub_og_id,
                "organism_name": organism_name,
                "pub_gene_id": pub_gene_id,
                "description": description,
                "assembly_accession": assembly_accession,
                "genome_sequence": genome_sequence
            })
            
    arguments = []
    logs = []
            
    for file in fasta_files:
        fasta_file = os.path.join(folder, file)
        gene_name = file.split(".")[0]
        if args.processes <= 1:
            logs.append(worker(fasta_file, os.path.join(output, file), data[gene_name], args.align_method))
        else:
            arguments.append((fasta_file, os.path.join(output, file), data[gene_name], args.align_method))
            
    if args.processes > 1:
        with Pool(args.processes) as p:
            logs = p.starmap(worker, arguments)
            
    with open(os.path.join(output, "log.txt"), "w") as f:
        for log in logs:
            for line in log:
                f.write(line+"\n")
            
    print("Done!")