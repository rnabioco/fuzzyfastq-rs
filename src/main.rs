use std::env;
use std::fs::{self, File};
use std::io::{BufReader, Read};
use std::path::Path;
use std::collections::HashMap;
use std::time::Instant;
use fastq::{Parser, Record};
use csv::Reader;
use flate2::read::MultiGzDecoder;

struct SequenceInfo {
    name: String,
    sequence: String,
    count: usize,
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 4 {
        eprintln!("Usage: {} <mode: 'seq' or 'csv'> <sequence or path_to_sequences.csv> <fastq_directory> [mismatch_percentage]", args[0]);
        return;
    }
    let mode = &args[1];
    let sequences_input = &args[2];
    let fastq_directory = &args[3];
    let mismatch_percentage: f64 = if args.len() > 4 {
        args[4].parse().expect("Invalid mismatch percentage")
    } else {
        0.0 // Default value if the argument is not provided
    };

    let mut sequences = HashMap::new();
    let mut sequence_order: Vec<String> = Vec::new();

    match mode.as_str() {
        "--seq" => {
            // Directly use the provided sequence
            let name = "Query".to_string();
            let sequence = sequences_input.to_uppercase();
            sequences.insert(name.clone(), SequenceInfo { name, sequence, count: 0 });
            sequence_order.push("Query".to_string());
        },
        "--csv" => {
            // Read sequences from the CSV file
            let mut rdr = Reader::from_path(sequences_input).expect("Failed to read CSV file");
            for result in rdr.records() {
                let record = result.expect("Error reading CSV record");
                let name = record[0].to_string();
                let sequence = record[1].to_uppercase();
                sequences.insert(name.clone(), SequenceInfo { name: name.clone(), sequence, count: 0 });
                sequence_order.push(name);
            }
        },
        _ => {
            eprintln!("Invalid mode. Use '--seq' for direct sequence input or '--csv' for CSV file input.");
            return;
        }
    }
    // Process each FASTQ file in the directory
    let paths = fs::read_dir(fastq_directory).expect("Failed to read directory");
    for path in paths {
        let path = path.expect("Error reading path").path();
        let file_ext = path.extension().and_then(std::ffi::OsStr::to_str).unwrap_or_default();

        if file_ext == "fastq" || file_ext == "fq" {
            process_fastq_file(&path, &mut sequences, mismatch_percentage, &sequence_order);
        } else if file_ext == "gz" {
            if let Some(stem) = path.file_stem().and_then(std::ffi::OsStr::to_str) {
                if stem.ends_with(".fastq") || stem.ends_with(".fq") {
                    process_fastq_file(&path, &mut sequences, mismatch_percentage, &sequence_order);
                }
            }
        }
    }
}

fn process_fastq_file(path: &Path, sequences: &mut HashMap<String, SequenceInfo>, mismatch_percentage: f64, sequence_order: &[String]) {
    let start = Instant::now();

    let file = File::open(path).expect("Failed to open FASTQ file");
    let file_ext = path.extension().and_then(std::ffi::OsStr::to_str).unwrap_or_default();
    
    let box_reader: Box<dyn Read> = if file_ext == "gz" {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let parser = Parser::new(box_reader);
    println!("Processing file: {:?}", path.file_name().unwrap());
    let mut total_reads = 0;

    let result = parser.each(|record| {
        total_reads += 1;
        let read = String::from_utf8_lossy(record.seq()).to_uppercase();
        println!("Sequence: {}", read);

        for seq_info in sequences.values_mut() {
            if is_match(&read, &seq_info.sequence, mismatch_percentage) {
                seq_info.count += 1;
            }
        }

        true
    });

    if let Err(e) = result {
        eprintln!("Error processing FASTQ file: {}", e);
        return;
    }
    println!("-----------------------------------\n");
    println!("Results for file: {:?}", path.file_name().unwrap());
    println!("Total reads: {}, Mismatch allowance: {:.2} \n", total_reads, mismatch_percentage);
    println!("Name, Sequence, Count, Percentage");
    for name in sequence_order.iter() {
        if let Some(seq_info) = sequences.get(name) {
            let percentage = (seq_info.count as f64 / total_reads as f64) * 100.0;
            println!("{}, {}, {}. {:.2}%", seq_info.name, seq_info.sequence, seq_info.count, percentage);
        }
    }

    let duration = start.elapsed();
    println!("\nTime taken for {:?}: {:?}\n", path.file_name().unwrap(), duration);

    // Reset the counts for sequences after processing each file
    for seq_info in sequences.values_mut() {
        seq_info.count = 0;
    }
}

fn is_match(read: &str, sequence: &str, mismatch_percentage: f64) -> bool {
    let mismatches_allowed = (sequence.len() as f64 * mismatch_percentage).round() as usize;
    let sequence_len = sequence.len();

    // Iterate over all possible alignments of the sequence within the read
    for start in 0..=read.len().saturating_sub(sequence_len) {
        let mut mismatches = 0;

        for (read_char, sequence_char) in read[start..].chars().zip(sequence.chars()) {
            if read_char != sequence_char {
                mismatches += 1;
                if mismatches > mismatches_allowed {
                    break;
                }
            }
        }

        // Return true as soon as a match is found
        if mismatches <= mismatches_allowed {
            return true;
        }
    }

    false
}