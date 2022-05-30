mod config;

use std::io::{Write, BufRead};
use std::collections::{HashMap, HashSet};

#[inline(always)]
fn make_features(
    bed_file_path: &str,
    genome_file_path: &str,
    whitelist_file: &str,
    out_file: &str
) {
    let filter_cells : HashSet<String> = {
        let file = std::fs::File::open(whitelist_file).unwrap();
        let reader = std::io::BufReader::new(file);
    
        let mut filter_cells = Vec::with_capacity(20_000);
        for line in reader.lines() {
            let line = line.unwrap();
            let cb = line.to_string();
    
            filter_cells.push(cb);
        }
        HashSet::from_iter(filter_cells)
    };

    let mut bed_reader = bio::io::bed::Reader::from_file(&bed_file_path).unwrap();
    let mut genome_reader = bio::io::fasta::IndexedReader::from_file(&genome_file_path).unwrap();

    let num_cells = filter_cells.len();
    let mut feats = HashMap::with_capacity(num_cells);
    for (line_num, record) in bed_reader.records().enumerate() {
        let bed_rec = record.expect("Error reading record.");
        let cb_str = bed_rec.name().unwrap().to_string();
        if ! filter_cells.contains(&cb_str) {
            continue;
        }

        genome_reader
            .fetch(bed_rec.chrom(), bed_rec.start(), bed_rec.end())
            .expect("Couldn't fetch interval");

        let mut seq = Vec::new();
        genome_reader.read(&mut seq).expect("Couldn't read the interval");

        if seq.len() < config::KMER {
            continue;
        }

        let total_kmers = seq.len() - config::KMER + 1;
        (0..total_kmers).for_each(|i| {
            let kmer = carina::barcode::cb_string_to_u64(&seq[i..i + config::KMER]).unwrap();

            let count = feats.entry(cb_str.to_owned())
                .or_insert(HashMap::new())
                .entry(kmer)
                .or_insert(0);
            *count += 1;
        });

        if line_num % 1_000_000 == 0 {
            print!("\r Done {}_M", line_num/1_000_000);
            std::io::stdout().flush().unwrap();
        }

        //if line_num > 10_000_000 { break; }
    }

    let mut row_names = Vec::with_capacity(num_cells);
    let mut mat: sprs::TriMat<usize> = sprs::TriMat::new((num_cells, 16384));
    for (row_index, (cb, cb_data)) in feats.into_iter().enumerate() {
        for (kmer, count) in cb_data {
            mat.add_triplet(row_index, kmer as usize, count);
        }
        row_names.push(cb);
    }

    let mat: sprs::CsMat<usize> = mat.to_csc();
    let mat_file_path = std::path::Path::new(out_file).join("mat.mtx");
    sprs::io::write_matrix_market(mat_file_path, &mat).unwrap();

    let file_path = std::path::Path::new(out_file).join("rows.txt");
    let mut file = std::io::BufWriter::new(std::fs::File::create(file_path).unwrap());
    for rname in row_names {
        writeln!(file, "{}", rname).unwrap();
    }
}


fn main() {
    let genome_file_path = "/brahms/srivastavaa/norma/data/genome/GRCh38.primary_assembly.genome.fa";
    let bed_file_path = "/brahms/srivastavaa/norma/data/multi/fragments.tsv";
    let whitelist_file = "/brahms/srivastavaa/norma/data/multi/cells.txt";
    let mat_file = "/brahms/srivastavaa/norma/data/multi/mat/";

    let kind = "tuple";
    match kind {
        "tuple" => {
            make_features(bed_file_path, genome_file_path, whitelist_file, mat_file);
        },
        _ => unreachable!()
    };

    println!("\nAll Done");
}
