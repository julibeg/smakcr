use std::{
    fs::File,
    io::{stdout, BufWriter, Write},
    str, vec,
};

use anyhow::Result;
use bio::alphabets::dna;
use clap::{value_parser, Arg, ArgAction, Command};

use smakcr::read_fasta;

fn initialize_stuff(k: usize) -> ([usize; 256], usize, Vec<u32>) {
    measure_time::info_time!("initialize");
    // key map for converting 'ACGT' to 0, 1, 2, 3 (same for lowercase bases)
    let mut key_map = [usize::MAX; 256];
    for bases in [b"ACGT", b"acgt"].iter() {
        for (i, &b) in bases.iter().enumerate() {
            key_map[b as usize] = i;
        }
    }
    // define mask to set left-shifted bits to 0 (used below)
    let mask = 4usize.pow(k as u32) - 1;

    // initiate vector with counts
    let counts: Vec<u32> = vec![0; 4usize.pow(k as u32)];

    (key_map, mask, counts)
}

fn index_to_kmer(index: usize, k: usize) -> Vec<u8> {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut idx = index;
    let mut kmer = vec![b'N'; k];
    for item in kmer.iter_mut() {
        // extract the last two bits
        *item = bases[idx & 3];
        idx >>= 2;
    }
    // we still need to reverse to get the correct order
    kmer.reverse();
    kmer
}

fn kmer_to_index(kmer: &Vec<u8>, key_map: [usize; 256]) -> usize {
    let mut idx = 0;
    for b in kmer {
        idx = (idx << 2) + key_map[*b as usize];
    }
    idx
}

fn count_kmers_in_sequence(
    sequence: &[u8],
    ksize: usize,
    counts: &mut [u32],
    mask: usize,
    key_map: [usize; 256],
) {
    if sequence.len() < ksize {
        return;
    }

    // keep a skip counter to avoid counting kmers with unknown bases (and the first few nucleotides
    // that don't form a complete kmer)
    let mut skip = ksize - 1;
    let mut idx = 0;

    for b in sequence {
        idx = (idx << 2) & mask;
        let add = key_map[*b as usize];
        if add == usize::MAX {
            // unknown base -> skip this kmer and the next `ksize - 1` kmers
            skip = ksize - 1;
            continue;
        }
        idx += add;
        // only count the kmer if not skipping
        if skip == 0 {
            counts[idx] += 1;
        } else {
            skip -= 1;
        }
    }
}

fn count_kmers_in_all_files(
    input_files: &Vec<&String>,
    k: usize,
    counts: &mut [u32],
    mask: usize,
    key_map: [usize; 256],
) -> Result<()> {
    measure_time::info_time!("count");
    for fasta_fn in input_files {
        for record in read_fasta(fasta_fn)? {
            // count kmers in sequence and handle unknown base error
            count_kmers_in_sequence(record?.seq(), k, counts, mask, key_map);
        }
    }
    Ok(())
}

fn write_results(
    outfile: Option<&String>,
    counts: &[u32],
    k: usize,
    key_map: [usize; 256],
    canonical: bool,
    non_zero: bool,
) -> Result<()> {
    measure_time::info_time!("write results");
    // create output file or write to STDOUT
    let out: Box<dyn Write> = match outfile {
        Some(outfile) => Box::new(File::create(outfile)?),
        None => Box::new(stdout()),
    };

    let mut out_writer = BufWriter::new(out);

    if canonical {
        // look up the reverse complement of the kmer and print the smaller of the two (alongside
        // the sum of the counts)
        for (index, &count) in counts.iter().enumerate() {
            let kmer = index_to_kmer(index, k);
            // TODO: implement a way to get the index of the reverse complement directly
            let revcomp = dna::revcomp(&kmer);
            // only write if this kmer is canonical (i.e. smaller than the reverse complement)
            if kmer > revcomp {
                continue;
            }
            // for a palindromic kmer (identical to the reverse complement), write the counts;
            // otherwise, write the sum of both counts
            let combined_count = if kmer == revcomp {
                // kmer is palindrome
                count
            } else {
                // kmer is smaller than its reverse complement
                count + counts[kmer_to_index(&revcomp, key_map)]
            };
            if non_zero && combined_count == 0 {
                continue;
            }
            out_writer.write_all(
                format!("{}\t{}\n", str::from_utf8(&kmer)?, combined_count).as_bytes(),
            )?;
        }
    } else {
        // write all kmers and their counts
        for (index, &count) in counts.iter().enumerate() {
            if non_zero && count == 0 {
                continue;
            };
            let kmer = index_to_kmer(index, k);
            out_writer.write_all(format!("{}\t{}\n", str::from_utf8(&kmer)?, count).as_bytes())?;
        }
    }
    Ok(())
}

fn main() -> Result<()> {
    env_logger::init();
    measure_time::info_time!("total");
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
            Arg::new("NONZERO")
                .short('n')
                .long("non-zero")
                .action(ArgAction::SetTrue)
                .help("Output only non-zero counts"),
        )
        .get_matches();

    let input_files: Vec<_> = matches.get_many::<String>("FASTA").unwrap().collect();
    let k: usize = *matches.get_one("SIZE").unwrap();
    let output = matches.get_one::<String>("OUT");
    let canonical = matches.get_flag("CANONICAL");
    let non_zero = matches.get_flag("NONZERO");

    // key map for converting 'ACGT' to 0, 1, 2, 3
    let (key_map, mask, mut counts) = initialize_stuff(k);

    count_kmers_in_all_files(&input_files, k, &mut counts, mask, key_map)?;

    // write output (perform the reverse operation (going from index to kmer) for all indices and
    // print the non-zero counts)
    write_results(output, &counts, k, key_map, canonical, non_zero)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_count_kmers_in_sequence() {
        let k = 3;
        let (key_map, mask, mut counts) = initialize_stuff(k);
        count_kmers_in_sequence(b"AcGTaCgAaaA", k, &mut counts, mask, key_map);

        let acg_index = 0b00_01_10; // ACG; comes up twice
        let cgt_index = 0b01_10_11; // CGT
        let gta_index = 0b10_11_00; // GTA
        let tac_index = 0b11_00_01; // TAC
        let cga_index = 0b01_10_00; // CGA
        let gaa_index = 0b10_00_00; // GAA
        let aaa_index = 0b00_00_00; // AAA; comes up twice

        assert_eq!(counts[acg_index], 2);
        assert_eq!(counts[cgt_index], 1);
        assert_eq!(counts[gta_index], 1);
        assert_eq!(counts[tac_index], 1);
        assert_eq!(counts[cga_index], 1);
        assert_eq!(counts[gaa_index], 1);
        assert_eq!(counts[aaa_index], 2);

        // make sure all other counts are 0
        assert_eq!(counts.iter().sum::<u32>(), 9);
    }

    #[test]
    fn test_count_kmers_with_short_sequence() {
        let k = 5;
        let (key_map, mask, mut counts) = initialize_stuff(k);
        count_kmers_in_sequence(b"ACG", k, &mut counts, mask, key_map);
        assert!(counts.iter().all(|&count| count == 0));
    }

    #[test]
    fn test_index_to_kmer() {
        let k = 3;
        let kmer = index_to_kmer(0b00_01_10, k);
        assert_eq!(kmer, b"ACG");

        let acg_index = 0b00_01_10; // ACG
        let cgt_index = 0b01_10_11; // CGT
        let gta_index = 0b10_11_00; // GTA
        let tac_index = 0b11_00_01; // TAC
        let cga_index = 0b01_10_00; // CGA
        let gaa_index = 0b10_00_00; // GAA
        let aaa_index = 0b00_00_00; // AAA

        assert_eq!(index_to_kmer(acg_index, k), b"ACG");
        assert_eq!(index_to_kmer(cgt_index, k), b"CGT");
        assert_eq!(index_to_kmer(gta_index, k), b"GTA");
        assert_eq!(index_to_kmer(tac_index, k), b"TAC");
        assert_eq!(index_to_kmer(cga_index, k), b"CGA");
        assert_eq!(index_to_kmer(gaa_index, k), b"GAA");
        assert_eq!(index_to_kmer(aaa_index, k), b"AAA");
    }

    #[test]
    fn test_unknwn_bases() {
        let k = 3;
        let (key_map, mask, mut counts) = initialize_stuff(k);
        count_kmers_in_sequence(b"AaAxcGtx", k, &mut counts, mask, key_map);

        let aaa_index = 0b00_00_00; // AAA
        let cgt_index = 0b01_10_11; // CGT
        assert_eq!(counts[aaa_index], 1);
        assert_eq!(counts[cgt_index], 1);
        // make sure all other counts are 0
        assert_eq!(counts.iter().sum::<u32>(), 2);
    }

    // TODO: add tests for
    // - the whole program (i.e. against an existing file to make sure the output format etc is as
    //   expected)
}
