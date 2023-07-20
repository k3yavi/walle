use crate::config;

use rust_htslib::tbx::Read;
use std::io::{Write, BufRead};
use std::collections::{HashMap};

pub fn get_whitelist_cells(whitelist_file: &str) -> HashMap<String, usize> {
    let file = std::fs::File::open(whitelist_file).unwrap();
    let reader = std::io::BufReader::new(file);
    
    let mut filter_cells = HashMap::new();
    for (row_id, line) in reader.lines().enumerate() {
        let line = line.unwrap();
        let cb = line.to_string();

        assert_eq!(filter_cells.contains_key(&cb), false);
        filter_cells.insert(cb, row_id);
    }
    
    filter_cells
}

pub fn get_genome_sequence(
    chr_name: &str,
    genome_reader: &mut bio::io::fasta::IndexedReader<std::fs::File>
) -> Vec<u8> {
    genome_reader.fetch_all(&chr_name)
        .expect("Couldn't fetch interval");

    let mut genome = Vec::new();
    genome_reader.read(&mut genome).expect("Couldn't read the interval");

    genome
}

pub fn get_bed_records(
    chr_len: u64,
    chr_name: &str,
    cells: &HashMap<String, usize>,
    bed_reader: &mut rust_htslib::tbx::Reader,
) -> Option<Vec<(usize, usize, usize)>> {
    // get bed entries
    match bed_reader.tid(&chr_name) {
        Ok(chr_id) => {
            bed_reader
                .fetch(chr_id, 0, chr_len)
                .expect("Could not seek to fetch region");

            let all_records: Vec<_> = bed_reader
                .records()
                .map(|x| {
                    let bed_line = String::from_utf8(x.unwrap())
                        .expect("UTF8 conversion error");
                    let toks: Vec<&str> = bed_line.split_whitespace()
                        .collect();
                
                    let cb = toks[3].to_string();
                    match cells.get(&cb) {
                        Some(&cb_index) => {
                            let start = toks[1].parse::<usize>().unwrap();
                            let end = toks[2].parse::<usize>().unwrap();
                            match end - start >= config::KMER {
                                true => Some((start, end, cb_index)),
                                false => None
                            }
                        },
                        None => None,
                    }
                })
                .flatten()
                .collect();
            Some(all_records)
        },
        _ => None
    }
}

fn _get_features(
    genome: Vec<u8>,
    num_cells: usize,
    records: Vec<(usize, usize, usize)>,
) -> HashMap<usize, HashMap<u64, usize>> {
    let pbar = indicatif::ProgressBar::new(records.len() as u64);
    pbar.set_style(
        indicatif::ProgressStyle::default_bar()
            .template(
                "{spinner:.red} [{elapsed_precise}] [{bar:40.cyan/blue}] {percent}/100% {msg}",
            )
            .progress_chars("╢▌▌░╟"),
    );

    let mut feats = HashMap::with_capacity(num_cells);
    for record in records {
        let cb_index = record.2;

        let total_kmers = record.1 - (config::KMER+record.0) + 2;
        (0..total_kmers).for_each(|i| {
            let start = record.0 + i;
            let end = start + config::KMER;

            let kmer = carina::barcode::cb_string_to_u64(&genome[start..end]).unwrap();
            let count = feats.entry(cb_index)
                .or_insert(HashMap::new())
                .entry(kmer)
                .or_insert(0);
            *count += 1;
        });
        pbar.inc(1);
    }
    feats
}

fn _write_features(
    out_file: &str,
    chr_name: &str,
    num_cells: usize,
    num_features: usize,
    feats: HashMap<usize, HashMap<u64, usize>>,
) {
    let mut mat: sprs::TriMat<usize> = sprs::TriMat::new((num_features, num_cells));
    for (cb_index, cb_data) in feats {
        for (kmer, count) in cb_data {
            mat.add_triplet(kmer as usize, cb_index, count);
        }
    }

    let mat: sprs::CsMat<usize> = mat.to_csc();
    let mat_file_path = std::path::Path::new(out_file).join(format!("mat_{}.mtx", &chr_name));
    sprs::io::write_matrix_market(mat_file_path, &mat).unwrap();
}

fn _make_features(
    bed_file_path: &str,
    genome_file_path: &str,
    whitelist_file: &str,
    out_file: &str
) {
    let filter_cells = get_whitelist_cells(whitelist_file);
    let num_cells = filter_cells.len();
    info!("Found total {} cells", num_cells);

    // genome iterator
    let mut genome_reader = bio::io::fasta::IndexedReader::from_file(&genome_file_path).unwrap();

    // bed iterator
    let mut bed_reader = rust_htslib::tbx::Reader::from_path(&bed_file_path)
        .unwrap_or_else(|_| panic!("Could not open {:?}", bed_file_path));
    bed_reader.set_threads(4).unwrap();

    {
        // col names
        let mut col_names = vec!["".to_string(); num_cells];
        for (cell, &index) in filter_cells.iter() {
            col_names[index] = cell.to_owned();
        }
    
        let file_path = std::path::Path::new(out_file).join(format!("cols.txt"));
        let mut file = std::io::BufWriter::new(std::fs::File::create(file_path).unwrap());
        for rname in col_names {
            writeln!(file, "{}", rname).unwrap();
        }    
    }

    let num_features: u64 = (4 as u64).pow(config::KMER as u32);
    {   // row names
        let file_path = std::path::Path::new(out_file).join(format!("rows.txt"));
        let mut file = std::io::BufWriter::new(std::fs::File::create(file_path).unwrap());
        for id in 0..num_features {
            writeln!(file, "{}", carina::barcode::u64_to_cb_string(id, config::KMER).unwrap()).unwrap();
        }
    }

    let num_chromosomes = 23;
    let mut mat = vec![vec![0; num_cells]; num_features as usize];
    //(1..num_chromosomes).rev().filter(|&x| x == 6).for_each(|chr_number| {
    (1..num_chromosomes).rev().for_each(|chr_number| {
        let chr_name = format!("chr{}", chr_number);
        let chr_len = config::CHR_LENS[chr_number-1];
        info!("Working on {} w/ length {}", chr_name, chr_len);

        let genome = get_genome_sequence(&chr_name, &mut genome_reader);
        info!("Sequence Extraction Complete w/ length {}", genome.len());

        if let Some(all_records) = get_bed_records(chr_len, &chr_name, &filter_cells, &mut bed_reader) {
            info!("BED record extraction complete w/ {} records", all_records.len());

            let feats = _get_features(genome, num_cells, all_records);
            info!("Feature extraction complete");

            for (cb_index, cb_data) in feats {
                for (kmer, count) in cb_data {
                    mat[kmer as usize][cb_index] += count;
                }
            }
        }
        //write_features(out_file, &chr_name, num_cells, num_features as usize, feats);
    });

    let mut tri_mat: sprs::TriMat<usize> = sprs::TriMat::new((num_features as usize, num_cells));
    for (row_index, row_data) in mat.into_iter().enumerate() {
        for (col_index, count) in row_data.into_iter().enumerate() {
            if count != 0 {
                tri_mat.add_triplet(row_index, col_index, count);
            }
        }
    }

    let mat: sprs::CsMat<usize> = tri_mat.to_csc();
    let mat_file_path = std::path::Path::new(out_file).join(format!("mat.mtx"));
    sprs::io::write_matrix_market(mat_file_path, &mat).unwrap();
}
