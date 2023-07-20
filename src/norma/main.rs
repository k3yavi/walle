pub mod em;

use std::path::{PathBuf, Path};
use std::collections::HashMap;

use clap::{App, Arg, SubCommand};

fn main() {
    let matches = App::new("norma")
        .version("0.1.0")
        .author("Avi Srivastava")
        .about("General Analyis")
        .subcommand(
            SubCommand::with_name("quantify")
                .about("Quantify")
                .arg(
                    Arg::with_name("input")
                        .short("i")
                        .takes_value(true)
                        .required(true)
                        .help("path to the input file"),
                )
                .arg(
                    Arg::with_name("output")
                        .short("o")
                        .takes_value(true)
                        .required(true)
                        .help("path to the output file"),
                ),
        )
        .get_matches();

    let ifile_path: PathBuf;
    let ofile_path: &Path;
    if let Some(sub_m) = matches.subcommand_matches("quantify") {
        ifile_path = carina::file::file_path_from_clap(sub_m, "input").unwrap();   
        let file_str = sub_m
            .value_of("output")
            .expect(&format!("can't find the flag: output"));
        ofile_path = Path::new(file_str);
    } else {
        panic!("Wrong subcommand");
    }

    //let file_path = Path::new("/home/srivastavaa/bin/norma/data/ctrl.tsv");

    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_path(ifile_path)
        .unwrap();

    let mut num_lines = 0;
    let mut triplets = HashMap::new();
    for line in rdr.records() {
        num_lines += 1;
        let record = line.unwrap();
        let values: Vec<String> = record.into_iter().flat_map(str::parse::<String>).collect();
        assert_eq!(values.len(), 3);

        let reference: String = values[0].clone();
        let query: String = values[1].clone();
        let _score: u8 = (values[2].parse::<f32>().unwrap() * 100.0) as u8;

        let val = triplets.entry(query)
            .or_insert(Vec::new());
        (*val).push(reference);
    }
    println!("Read {} lines", num_lines);

    let mut eqclass = HashMap::new();
    for (_, val) in triplets {
        let v = eqclass.entry(val)
            .or_insert(0);
        *v += 1
    }
    println!("Number of Eqclasses {:?}", eqclass.len()); 

    let (alphas, ref_names_hash) = em::optimize(eqclass);
    let mut ref_names = vec!["".to_string(); alphas.len()];
    for (name, idx) in ref_names_hash {
        ref_names[idx] = name;
    }

    //let out_file_path = Path::new("/home/srivastavaa/bin/norma/data/out.tsv");
    let mut wtr = csv::WriterBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_path(ofile_path)
        .unwrap();
    for i in 0..alphas.len() {
        wtr.write_record(&[ref_names[i].clone(), alphas[i].to_string()]).unwrap();
    }

    println!("All Done");
}
