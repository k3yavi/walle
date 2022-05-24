mod config;
mod old;

fn main() {
    let genome_file_path = "/brahms/srivastavaa/norma/data/k4me3/full_data/GRCh38.primary_assembly.genome.fa";
    let bed_file_path = "/brahms/srivastavaa/norma/data/k4me3/full_data/H3K4me3_fragments.tsv";

    let tuple_file = "/brahms/srivastavaa/norma/data/k4me3/all.tuples.tsv";

    let whitelist_file = "/brahms/srivastavaa/norma/data/k4me3/cells.txt";
    let cell_file = "/brahms/srivastavaa/norma/data/k4me3/cell.filtered.tuples.tsv";

    let freq_file = "/brahms/srivastavaa/norma/data/k4me3/kmer.frequency.tsv";
    let filter_file = "/brahms/srivastavaa/norma/data/k4me3/kmer.selected.tsv";
    let fragment_file = "/brahms/srivastavaa/norma/data/k4me3/kmer.filtered.tuples.tsv";
    
    let stats_file = "/brahms/srivastavaa/norma/data/k4me3/kmer.stats.txt";
    let stats_in_file = "/brahms/srivastavaa/norma/data/k4me3/stats.selected.kmer.tsv";
    let stats_out_file = "/brahms/srivastavaa/norma/data/k4me3/stats.tuple.txt";

    let mat_file = "/brahms/srivastavaa/norma/data/k4me3/mat/";

    let kind = "last_matrix";
    match kind {
        "tuple" => {        
            old::make_tuples(bed_file_path, genome_file_path, tuple_file);
            println!("Done making tuples");
        },
        "cells" => {
            old::filter_cells(tuple_file, whitelist_file, cell_file);
            println!("Done filtering cells");
        },
        "freq" => {
            old::kmer_frequency(cell_file, freq_file);
            println!("Done generating frequency");
        },
        "filter" => {
            old::filter_kmer(freq_file, filter_file);
            println!("Done filtering kmers");
        },
        "fragment" => {
            old::filter_fragment(cell_file, filter_file, fragment_file);
            println!("Done filtering fragments");
        },
        "stats" => {
            old::get_stats(fragment_file, stats_file);
            println!("Done grouping cells");
        },
        "bulk" => {
            old::make_tuples(bed_file_path, genome_file_path, tuple_file);
            println!("Done making tuples");

            old::filter_cells(tuple_file, whitelist_file, cell_file);
            println!("Done filtering cells");

            old::kmer_frequency(cell_file, freq_file);
            println!("Done generating frequency");

            old::filter_kmer(freq_file, filter_file);
            println!("Done filtering kmers");

            old::filter_fragment(cell_file, filter_file, fragment_file);
            println!("Done filtering fragments");

            old::get_stats(fragment_file, stats_file);
            println!("Done grouping cells");
        },
        "last" => {
            old::filter_fragment(fragment_file, stats_in_file, stats_out_file);
            println!("Done last cells");
        },
        "matrix" => {
            old::write_matrix(stats_out_file, mat_file);
            println!("Done writing matrix");
        },
        "last_matrix" => {
            old::filter_fragment(fragment_file, stats_in_file, stats_out_file);
            println!("Done last cells");

            old::write_matrix(stats_out_file, mat_file);
            println!("Done writing matrix");
        },

        _ => unreachable!()
    };

    println!("\nAll Done");
}
