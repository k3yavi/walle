use std::io::Write;

mod config;

//fn second() {
//    let val = kmer_map.entry(kmer)
//    .or_insert(HashMap::new())
//    .entry(cb_int)
//    .or_insert(0);
//    *val += 1;
//}

fn first(
    bed_file_path: &str,
    genome_file_path: &str,
    out_file: &str
) {
    let mut bed_reader = bio::io::bed::Reader::from_file(&bed_file_path).unwrap();
    let mut reader = bio::io::fasta::IndexedReader::from_file(&genome_file_path).unwrap();

    let mut file = std::fs::File::create(out_file).unwrap();
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

        if line_num % 1_0_000 == 0 {
            println!("\r Done {}_0K", line_num/1_0_000);
            std::io::stdout().flush().unwrap();
        }
    }
}

fn main() {
    let genome_file_path = "/brahms/srivastavaa/norma/data/k4me1/GRCh38.primary_assembly.genome.fa";
    let bed_file_path = "/brahms/srivastavaa/norma/data/k4me1/H3K4me1_fragments.tsv";
    
    let out_file = "/brahms/srivastavaa/norma/data/k4me1/tuples.tsv";
    first(bed_file_path, genome_file_path, out_file);

    println!("\nAll Done");
}
