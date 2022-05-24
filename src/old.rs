use std::collections::{HashMap, HashSet};
use std::io::{Write, BufRead};
use std::iter::FromIterator;

use crate::config;

pub fn make_tuples(
    bed_file_path: &str,
    genome_file_path: &str,
    out_file: &str
) {
    let mut bed_reader = bio::io::bed::Reader::from_file(&bed_file_path).unwrap();
    let mut reader = bio::io::fasta::IndexedReader::from_file(&genome_file_path).unwrap();

    let mut file = std::io::BufWriter::new(std::fs::File::create(out_file).unwrap());
    for (line_num, record) in bed_reader.records().enumerate() {
        let bed_rec = record.expect("Error reading record.");
        reader
            .fetch(bed_rec.chrom(), bed_rec.start(), bed_rec.end())
            .expect("Couldn't fetch interval");

        let cb_str = bed_rec.name().unwrap();
        let cb_int =
            carina::barcode::cb_string_to_u64(&cb_str.as_bytes()[0..cb_str.len() - 2]).unwrap();

        let mut seq = Vec::new();
        reader.read(&mut seq).expect("Couldn't read the interval");

        if seq.len() < config::KMER {
            continue;
        }

        let total_kmers = seq.len() - config::KMER + 1;
        (0..total_kmers).for_each(|i| {
            let kmer = carina::barcode::cb_string_to_u64(&seq[i..i + config::KMER]).unwrap();
            writeln!(file, "{}\t{}", kmer, cb_int).unwrap();
        });

        if line_num % 1_000_000 == 0 {
            print!("\r Done {}_M", line_num/1_000_000);
            std::io::stdout().flush().unwrap();
        }
    }
}

pub fn filter_cells(
    tuple_file: &str,
    whitelist_file: &str,
    cell_file: &str,
) {
    let file = std::fs::File::open(whitelist_file).unwrap();
    let reader = std::io::BufReader::new(file);

    let mut filter_cells = Vec::with_capacity(20_000);
    for line in reader.lines() {
        let line = line.unwrap();
        let cb = carina::barcode::cb_string_to_u64(&line.as_bytes()[0..line.len() - 2]).unwrap();

        filter_cells.push(cb);
    }
    let filter_cells : HashSet<u64> = HashSet::from_iter(filter_cells);

    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_path(tuple_file)
        .unwrap();

    let mut file = std::io::BufWriter::new(std::fs::File::create(cell_file).unwrap());
    for (line_num, line) in rdr.records().enumerate() {
        let record = line.unwrap();
        let values: Vec<u64> = record.into_iter().flat_map(str::parse::<u64>).collect();
        assert_eq!(values.len(), 2);


        if filter_cells.contains(&values[1]) {
            writeln!(file, "{}\t{}", values[0], values[1]).unwrap();
        }

        if line_num % 1_000_000 == 0 {
            print!("\r Done {}_M", line_num/1_000_000);
            std::io::stdout().flush().unwrap();
        }
    }
}

pub fn kmer_frequency(
    tuple_file: &str,
    freq_file: &str
) {
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_path(tuple_file)
        .unwrap();

    let mut freq = HashMap::new();
    for (line_num, line) in rdr.records().enumerate() {
        let record = line.unwrap();
        let values: Vec<u64> = record.into_iter().flat_map(str::parse::<u64>).collect();
        assert_eq!(values.len(), 2);

        let val = freq.entry(values[0])
            .or_insert(0);
        *val += 1;

        if line_num % 1_000_000 == 0 {
            print!("\r Done {}_M", line_num/1_000_000);
            std::io::stdout().flush().unwrap();
        }
    }

    let mut file = std::io::BufWriter::new(std::fs::File::create(freq_file).unwrap());
    for (kmer, count) in freq {
        writeln!(file, "{}\t{}", kmer, count).unwrap();
    }
}

pub fn filter_kmer(
    freq_file: &str,
    filter_file: &str
) {
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_path(freq_file)
        .unwrap();

    let mut file = std::io::BufWriter::new(std::fs::File::create(filter_file).unwrap());
    for (line_num, line) in rdr.records().enumerate() {
        let record = line.unwrap();
        let values: Vec<u64> = record.into_iter().flat_map(str::parse::<u64>).collect();
        assert_eq!(values.len(), 2);


        if values[1] > 1_00 {
            writeln!(file, "{}", values[0]).unwrap();
        }

        if line_num % 10_000_000 == 0 {
            print!("\r Done {}0_M", line_num/10_000_000);
            std::io::stdout().flush().unwrap();
        }
    }
}

pub fn filter_fragment(
    tuple_file: &str,
    filter_file: &str,
    fragment_file: &str,
) {

    let file = std::fs::File::open(filter_file).unwrap();
    let reader = std::io::BufReader::new(file);

    let mut filter_hash = Vec::with_capacity(1_000_000);
    for line in reader.lines() {
        let line = line.unwrap().parse::<u64>().unwrap();
        filter_hash.push(line);
    }
    let filter_hash : HashSet<u64> = HashSet::from_iter(filter_hash);

    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_path(tuple_file)
        .unwrap();

    let mut file = std::io::BufWriter::new(std::fs::File::create(fragment_file).unwrap());
    for (line_num, line) in rdr.records().enumerate() {
        let record = line.unwrap();
        let values: Vec<u64> = record.into_iter().flat_map(str::parse::<u64>).collect();
        assert_eq!(values.len(), 2);


        if filter_hash.contains(&values[0]) {
            writeln!(file, "{}\t{}", values[0], values[1]).unwrap();
        }

        if line_num % 1_000_000 == 0 {
            print!("\r Done {}_M", line_num/1_000_000);
            std::io::stdout().flush().unwrap();
        }
    }
}

pub fn get_stats(
    cell_file: &str,
    stats_file: &str,
) {
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_path(cell_file)
        .unwrap();

    let mut kmap = HashMap::new();
    for (line_num, line) in rdr.records().enumerate() {
        let record = line.unwrap();
        let values: Vec<u64> = record.into_iter().flat_map(str::parse::<u64>).collect();
        assert_eq!(values.len(), 2);

        let val = kmap.entry(values[0])
            .or_insert(HashMap::new())
            .entry(values[1])
            .or_insert(0);
        *val += 1;

        if line_num % 1_000_000 == 0 {
            print!("\r Done {}_M", line_num/1_000_000);
            std::io::stdout().flush().unwrap();
        }
    }

    let mut file = std::io::BufWriter::new(std::fs::File::create(stats_file).unwrap());
    for (kmer, kmer_data) in kmap {
        let values: Vec<f32> = kmer_data.into_iter().map(|(_,v)| v as f32).collect();
        
        let num_cb = values.len();
        let vec_sum: f32 = values.iter().sum();
        
        let vec_mean = vec_sum / 17335.0;
        let mut vec_var: f32 = values.into_iter().map(|val| f32::powf(val-vec_mean, 2.0) ).sum();
        vec_var += f32::powf(vec_mean, 2.0) * (17335 - num_cb) as f32;
        vec_var = vec_var / 17335.0;

        writeln!(file, "{}\t{}\t{}\t{}", kmer, num_cb, vec_mean, vec_var).unwrap();
    }
}

pub fn write_matrix(
    stats_file: &str,
    mat_file: &str,
) {
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_path(stats_file)
        .unwrap();
    
    let mut row_names = HashSet::new();
    let mut col_names = HashSet::new();

    let mut kmap = HashMap::new();
    for (line_num, line) in rdr.records().enumerate() {
        let record = line.unwrap();
        let values: Vec<u64> = record.into_iter().flat_map(str::parse::<u64>).collect();
        assert_eq!(values.len(), 2);

        let kmer = values[0];
        let cb = values[1];

        row_names.insert(kmer);
        col_names.insert(cb);

        let val = kmap.entry(kmer)
            .or_insert(HashMap::new())
            .entry(cb)
            .or_insert(0);
        *val += 1;

        if line_num % 1_000_000 == 0 {
            print!("\r Done {}_M", line_num/1_000_000);
            std::io::stdout().flush().unwrap();
        }
    }

    let mut row_map = HashMap::new();
    for (index, rname) in row_names.into_iter().enumerate() {
        row_map.insert(rname, index);
    }

    let mut col_map = HashMap::new();
    for (index, cname) in col_names.into_iter().enumerate() {
        col_map.insert(cname, index);
    }

    let mut mat: sprs::TriMat<usize> = sprs::TriMat::new((row_map.len(), col_map.len()));
    for (kmer, kmer_data) in kmap {
        let row_index = row_map[&kmer];
        for (cb, count) in kmer_data {
            let col_index = col_map[&cb];
            mat.add_triplet(row_index, col_index, count);
        }
    }

    let mat: sprs::CsMat<usize> = mat.to_csc();
    let mat_file_path = std::path::Path::new(mat_file).join("mat.mtx");
    sprs::io::write_matrix_market(mat_file_path, &mat).unwrap();

    let mut row_names = vec![0; row_map.len()];
    for (kmer, index) in row_map {
        row_names[index] = kmer;
    }

    let file_path = std::path::Path::new(mat_file).join("rows.txt");
    let mut file = std::io::BufWriter::new(std::fs::File::create(file_path).unwrap());
    for rname in row_names {
        writeln!(file, "{}", carina::barcode::u64_to_cb_string(rname, config::KMER).unwrap()).unwrap();
    }

    let mut col_names = vec![0; col_map.len()];
    for (cb, index) in col_map {
        col_names[index] = cb;
    }
    
    let file_path = std::path::Path::new(mat_file).join("cols.txt");
    let mut file = std::io::BufWriter::new(std::fs::File::create(file_path).unwrap());
    for cname in col_names {
        writeln!(file, "{}", carina::barcode::u64_to_cb_string(cname, 16).unwrap()).unwrap();
    }
}
