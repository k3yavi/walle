#[macro_use]
extern crate log;

mod old;
mod config;

use std::io::Write;
use std::collections::HashMap;

fn get_minimizers(kmers: Vec<u64>) -> Option<Vec<u64>> {
    let total_kmers = kmers.len();
    let group_size = config::WINDOW - config::KMER + 1;
    if total_kmers < group_size { return None }

    let total_groups = total_kmers - group_size + 1;
    let mut minimizers: Vec<u64> = (0..total_groups).map(|group_start_index|{
        let group_end_index = group_start_index + group_size;
        *kmers[group_start_index..group_end_index].iter().min().unwrap()
    }).collect();
    minimizers.dedup();

    Some(minimizers)
}

fn get_class_features(
    genome: Vec<u8>,
    records: Vec<(usize, usize, usize)>,
    feats: &mut HashMap<usize, HashMap<Vec<u64>, u32>>
) {
    let pbar = indicatif::ProgressBar::new(records.len() as u64);
    pbar.set_style(
        indicatif::ProgressStyle::default_bar()
            .template(
                "{spinner:.red} [{elapsed_precise}] [{bar:40.cyan/blue}] {percent}/100% {msg}",
            )
            .progress_chars("╢▌▌░╟"),
    );

    for record in records {
        let cb_index = record.2;

        let total_kmers = record.1 - (config::KMER+record.0) + 2;
        let all_kmers: Vec<u64> = (0..total_kmers).map(|i| {
            let start = record.0 + i;
            let end = start + config::KMER;

            carina::barcode::cb_string_to_u64(&genome[start..end]).unwrap()
        }).collect();
        
        let minimizers = get_minimizers(all_kmers);
        if let Some(minimizers) = minimizers {
            // let count = feats.entry(cb_index)
            //     .or_insert(HashMap::new())
            //     .entry(minimizers)
            //     .or_insert(0);
            // *count += 1;
            let num_minimizers = minimizers.len();
            if num_minimizers >= config::MINIMIZER_GSIZE {
                let num_minimizer_groups = num_minimizers - config::MINIMIZER_GSIZE + 1;
                (0..num_minimizer_groups).for_each(|minimizer_start_index| {
                    let minimizer_end_index = minimizer_start_index + config::MINIMIZER_GSIZE;
                    let count = feats.entry(cb_index)
                        .or_insert(HashMap::new())
                        .entry(minimizers[minimizer_start_index..minimizer_end_index].to_vec())
                        .or_insert(0);
                    *count += 1;    
                });
            }
        };

        pbar.inc(1);
    }
}

fn make_features(
    bed_file_path: &str,
    genome_file_path: &str,
    whitelist_file: &str,
    out_file: &str
) {
    let filter_cells = old::get_whitelist_cells(whitelist_file);
    let num_cells = filter_cells.len();
    info!("Found total {} cells", num_cells);

    // genome iterator
    let mut genome_reader = bio::io::fasta::IndexedReader::from_file(&genome_file_path).unwrap();

    // bed iterator
    let mut bed_reader = rust_htslib::tbx::Reader::from_path(&bed_file_path)
        .unwrap_or_else(|_| panic!("Could not open {:?}", bed_file_path));
    bed_reader.set_threads(4).unwrap();

    let num_chromosomes = 23;
    let mut feats: HashMap<usize, HashMap<Vec<u64>, u32>> = HashMap::new();
    //(1..num_chromosomes).rev().filter(|&x| x == 19).for_each(|chr_number| {
    (1..num_chromosomes).rev().for_each(|chr_number| {
        let chr_name = format!("chr{}", chr_number);
        let chr_len = config::CHR_LENS[chr_number-1];
        info!("Working on {} w/ length {}", chr_name, chr_len);

        let genome = old::get_genome_sequence(&chr_name, &mut genome_reader);
        info!("Sequence Extraction Complete w/ length {}", genome.len());

        if let Some(all_records) = old::get_bed_records(chr_len, &chr_name, &filter_cells, &mut bed_reader) {
            info!("BED record extraction complete w/ {} records", all_records.len());

            get_class_features(genome, all_records, &mut feats);
        }
    });

    let mut tuples = HashMap::new();
    for (cb, cb_data) in feats {
            for (kmers, freq) in cb_data { 
                //if kmers.len() < 20 { continue }
                assert_eq!(kmers.len(), config::MINIMIZER_GSIZE);
                let count = tuples.entry(kmers)
                    .or_insert(HashMap::new())
                    .entry(cb)
                    .or_insert(0);
                *count += freq;
            }
    }

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

    let num_features = tuples.len();
    {   // row names
        let row_file_path = std::path::Path::new(out_file).join(format!("rows.txt"));
        let mut row_file = std::io::BufWriter::new(std::fs::File::create(row_file_path).unwrap());
        for id in 0..num_features {
           writeln!(row_file, "feat_{}", id).unwrap();
        }
    }

    let mut tri_mat: sprs::TriMat<usize> = sprs::TriMat::new((num_features as usize, num_cells));
    for (row_index, (_kmers, row_data)) in tuples.into_iter().enumerate() {
        for (col_index, count) in row_data {
            if count != 0 {
                tri_mat.add_triplet(row_index, col_index, count as usize);
            }
        }
    }

    let mat: sprs::CsMat<usize> = tri_mat.to_csc();
    let mat_file_path = std::path::Path::new(out_file).join(format!("mat.mtx"));
    sprs::io::write_matrix_market(mat_file_path, &mat).unwrap();
    
    // let file_path = std::path::Path::new(out_file).join(format!("lens.txt"));
    // let mut file = std::io::BufWriter::new(std::fs::File::create(file_path).unwrap());
    // for (_, cb_data) in feats {
    //     for (kmers, freq) in cb_data { 
    //         let val = kmers.len();
    //         (0..freq).for_each(|_| {
    //             writeln!(file, "{}", val).unwrap();
    //         });
    //     }
    // }
}

fn main() {
    let genome_file_path = "/brahms/srivastavaa/norma/data/genome/GRCh38.primary_assembly.genome.fa";
    let bed_file_path = "/brahms/srivastavaa/norma/data/multi/single-cell/intersect.bed.gz";
    let whitelist_file = "/brahms/srivastavaa/norma/data/multi/cells.txt";
    //let mat_file = "/brahms/srivastavaa/norma/data/multi/single-cell/class.txt";
    let mat_file = "/brahms/srivastavaa/norma/data/multi/single-cell/mat_b";
    pretty_env_logger::init_timed();

    make_features(bed_file_path, genome_file_path, whitelist_file, mat_file);
    println!("\nAll Done");
}