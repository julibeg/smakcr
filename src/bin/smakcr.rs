use rayon::prelude::*;

use anyhow::Result;
use clap::{value_parser, Arg, ArgAction, Command};

use smakcr::KmerCounter;

fn main() -> Result<()> {
    let version = env!("CARGO_PKG_VERSION");
    let matches = Command::new("faSize")
        .version(version)
        .author("Andrea Talenti <andrea.talenti@ed.ac.uk>")
        .about("Print counts of K-mers of a given size across one or more FASTA files")
        .arg(
            Arg::new("FASTA")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Input FASTA file(s)"),
        )
        .arg(
            Arg::new("SIZE")
                .required(true)
                .short('k')
                .default_value("3")
                .value_parser(value_parser!(usize))
                .help("K-mer size"),
        )
        .arg(
            Arg::new("THREADS")
                .short('t')
                .default_value("1")
                .value_parser(value_parser!(usize))
                .help("Number of threads"),
        )
        .arg(
            Arg::new("OUT")
                .short('o')
                .long("output")
                .help("Output file"),
        )
        .arg(
            Arg::new("CANONICAL")
                .short('c')
                .long("canonical")
                .action(ArgAction::SetTrue)
                .help("Output canonical k-mers"),
        )
        .arg(
            Arg::new("ZERO")
                .short('z')
                .long("zero-counts")
                .action(ArgAction::SetTrue)
                .help("Also output k-mers with zero counts"),
        )
        .get_matches();

    let input_files: Vec<_> = matches.get_many::<String>("FASTA").unwrap().collect();
    let k: usize = *matches.get_one("SIZE").unwrap();
    let n_threads: usize = *matches.get_one("THREADS").unwrap();
    let output = matches.get_one::<String>("OUT");
    let canonical = matches.get_flag("CANONICAL");
    let write_zeros = matches.get_flag("ZERO");

    let counter = match (input_files.len(), n_threads) {
        (1, _) | (_, 1) => {
            // either just one file or just one thread --> iterate over the files and process with
            // `n_threads` threads at a time
            let mut counter = KmerCounter::new(k, n_threads);

            for &file in input_files.iter() {
                let reader = smakcr::FastxReader::from_file(file)?;
                counter.count(reader)?;
            }
            counter
        }
        (_, _) => {
            // many files and many threads --> use rayon (this is not very efficient when there are
            // more threads than files but we'll figure out a better way later)
            rayon::ThreadPoolBuilder::new()
                .num_threads(n_threads)
                .build_global()?;

            input_files
                .par_iter()
                .map(|&file| -> Result<KmerCounter> {
                    let reader = smakcr::FastxReader::from_file(file)?;
                    let mut counter = KmerCounter::new(k, 1);
                    counter.count(reader)?;
                    Ok(counter)
                })
                .try_reduce_with(|acc, counter| acc + counter)
                .ok_or_else(|| anyhow::anyhow!("No input files to process"))??
        }
    };

    counter.write(output, canonical, write_zeros)?;

    Ok(())
}
