use anyhow::Result;
use clap::{value_parser, Arg, ArgAction, Command};

use smakcr::KmerCounter;

fn main() -> Result<()> {
    let version = env!("CARGO_PKG_VERSION");
    let matches = Command::new("smakcr")
        .version(version)
        .author("Julian Libiseller-Egger <37619875+julibeg@users.noreply.github.com>")
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
            Arg::new("CHUNK_SIZE")
                .long("chunk-size")
                .default_value("1048576")
                .value_parser(value_parser!(usize))
                .help("Chunk size in bytes for parallel counting"),
        )
        .arg(
            Arg::new("QUEUE_MB")
                .long("queue-mb")
                .default_value("64")
                .value_parser(value_parser!(usize))
                .help("Queue memory limit in MB for parallel counting"),
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

    let input_files: Vec<String> = matches
        .get_many::<String>("FASTA")
        .unwrap()
        .cloned()
        .collect();
    let k: usize = *matches.get_one("SIZE").unwrap();
    let n_threads: usize = *matches.get_one("THREADS").unwrap();
    let output = matches.get_one::<String>("OUT");
    let canonical = matches.get_flag("CANONICAL");
    let write_zeros = matches.get_flag("ZERO");
    let chunk_size: usize = *matches.get_one("CHUNK_SIZE").unwrap();
    let queue_mb: usize = *matches.get_one("QUEUE_MB").unwrap();

    let mut counter = KmerCounter::new(k, n_threads);
    counter.count_paths(&input_files, chunk_size, queue_mb * 1024 * 1024)?;

    counter.write(output, canonical, write_zeros)?;

    Ok(())
}
