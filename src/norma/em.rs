use std::collections::{HashMap, HashSet};

const MIN_OUTPUT_ALPHA: f32 = 0.01;
const ALPHA_CHECK_CUTOFF: f32 = 1e-2;

const MIN_ITER: u32 = 2;
const MAX_ITER: u32 = 100;
const REL_DIFF_TOLERANCE: f32 = 1e-2;

pub fn em_update(
    alphas_in: &[f32],
    alphas_out: &mut Vec<f32>,
    eqclasses: &HashMap<Vec<usize>, u32>,
) {
    // loop over all the eqclasses
    for (labels, count) in eqclasses {
        if labels.len() > 1 {
            let mut denominator: f32 = 0.0;
            for label in labels {
                denominator += alphas_in[*label as usize];
            }

            if denominator > 0.0 {
                let inv_denominator = *count as f32 / denominator;
                for label in labels {
                    let index = *label as usize;
                    let count = alphas_in[index] * inv_denominator;
                    alphas_out[index] += count;
                }
            }
        } else {
            let tidx = labels.get(0).expect("can't extract labels");
            alphas_out[*tidx as usize] += *count as f32;
        }
    }
}

pub fn optimize(
    eqclasses: HashMap<Vec<String>, u32>,
) -> (Vec<f32>, HashMap<String, usize>) {
    let ref_hash: HashSet<String> = eqclasses.keys()
        .flatten()
        .map(|x| x.to_owned())
        .collect();

    let mut ref_names = HashMap::new();
    for (cb_idx, cb_seq) in ref_hash.into_iter().enumerate() {
        ref_names.insert(cb_seq, cb_idx);
    }
    println!("Number of reference {}", ref_names.len());

    println!("Restructuring equivalence classes");
    let mut fast_eqclass = HashMap::new();
    let mut num_query = 0;
    for (labels, count) in eqclasses {
        let labels: Vec<usize> = labels.iter()
            .map(|x| *ref_names.get(x).unwrap())
            .collect();
        fast_eqclass.insert(labels, count);
        num_query += count;
    }

    println!("Number of query cells: {}", num_query);
    println!("Starting EM");
    let num_refs = ref_names.len();
    let mut alphas_in: Vec<f32> = vec![0.0; num_refs];
    let mut alphas_out: Vec<f32> = vec![0.0; num_refs];

    for (labels, count) in &fast_eqclass {
        if labels.len() == 1 {
            let idx = labels.get(0).expect("can't extract labels");
            alphas_in[*idx as usize] += *count as f32;
        }
    }

    let uni_prior = 1.0 / (num_refs as f32);
    alphas_in.iter_mut().for_each(|x| *x += uni_prior);


    let mut it_num: u32 = 0;
    let mut converged: bool = true;
    while it_num < MIN_ITER || (it_num < MAX_ITER && !converged) {
        // perform one round of em update
        em_update(&alphas_in, &mut alphas_out, &fast_eqclass);

        converged = true;
        let mut max_rel_diff = -f32::INFINITY;

        for index in 0..num_refs {
            if alphas_out[index] > ALPHA_CHECK_CUTOFF {
                let diff = alphas_in[index] - alphas_out[index];
                let rel_diff = diff.abs();

                max_rel_diff = match rel_diff > max_rel_diff {
                    true => rel_diff,
                    false => max_rel_diff,
                };

                if rel_diff > REL_DIFF_TOLERANCE {
                    converged = false;
                }
            } // end- in>out if

            alphas_in[index] = alphas_out[index];
            alphas_out[index] = 0.0_f32;
        } //end-for

        it_num += 1;
    }

    // update too small alphas
    alphas_in.iter_mut().for_each(|alpha| {
        if *alpha < MIN_OUTPUT_ALPHA {
            *alpha = 0.0_f32;
        }
    });

    (alphas_in, ref_names)
}
