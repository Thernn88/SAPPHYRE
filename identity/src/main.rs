use std::collections::HashMap;
use std::env;
use std::fs;

fn read_fasta_file(filename: &str) -> HashMap<String, String> {
    let mut sequences = HashMap::new();
    let mut sequence_id = None;
    let mut sequence = String::new();
    let contents = fs::read_to_string(filename).expect("Failed to read file");

    for line in contents.lines() {
        if line.starts_with('>') {
            if let Some(id) = sequence_id.take() {
                sequences.insert(id, sequence);
                sequence = String::new();
            }
            sequence_id = Some(line[1..].to_string());
        } else {
            sequence.push_str(line.trim());
        }
    }

    if let Some(id) = sequence_id {
        sequences.insert(id, sequence);
    }

    sequences
}


fn pairwise_identity(seq1: &str, seq2: &str) -> f64 {
    if seq1.len() != seq2.len() {
        panic!("Sequences must be of equal length");
    }

    let mut identical = 0;
    let mut aligned_length = 0;

    for (a, b) in seq1.chars().zip(seq2.chars()) {
        if a == '-' || b == '-' {
            continue;
        }
        aligned_length += 1;
        if a == b {
            identical += 1;
        }
    }

    if aligned_length == 0 {
        return 0.0;
    }

    (identical as f64 / aligned_length as f64) * 100.0
}

fn has_overlap(seq1: &String, seq2: &String) -> bool {
    for (a,b) in seq1.chars().zip(seq2.chars()) {
        if a != '-' && b != '-' {
            return true;
        }
    }
    false
}

fn average_pairwise_identity(sequences: &HashMap<String, String>) -> f64 {
    let sequence_ids = sequences.keys().cloned().collect::<Vec<_>>();
    let mut num_pairs = 0;
    let mut total_identity = 0.0;

    for (i, seq_id1) in sequence_ids.iter().enumerate() {
        for (j, seq_id2) in sequence_ids.iter().enumerate() {
            if i >= j {
                continue;
            }
            if !(has_overlap(&sequences[seq_id1], &sequences[seq_id2])) { continue; }
            let identity = pairwise_identity(&sequences[seq_id1], &sequences[seq_id2]);
            total_identity += identity;
            num_pairs += 1;
        }
    }

    if num_pairs == 0 {
        return 0.0;
    }

    total_identity / num_pairs as f64
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        eprintln!("Usage: {} <fasta_file>", args[0]);
        std::process::exit(1);
    }

    let filename = &args[1];
    let sequences = read_fasta_file(filename);
    let avg_identity = average_pairwise_identity(&sequences);
    println!("Average pairwise identity: {:.2}%", avg_identity);
}
