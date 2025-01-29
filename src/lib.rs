use anyhow::{Context, Result};
use bio::io::fasta;
use std::path::Path;
use std::{
    fs::File,
    io::{BufReader, Read},
    fmt::Display,
};

pub fn read_fasta<P: AsRef<Path> + Display>(
    path: P,
) -> Result<impl Iterator<Item = Result<bio::io::fasta::Record>>> {
    let sniffed_reader = niffler::from_path(&path);
    // if the file has fewer than 5 bytes, `niffler` can't sniff the compression format and will
    // return a `FileTooShort` error; this could be due to
    // * an empty file
    // * a file containing only a single FASTA record with the ID consisting only of a single
    //   character and the sequence being empty
    // * a file containing a single sequence with a one-character ID and one-character sequence and
    //   missing newline character at the end
    // we don't want to fail at this stage in these cases and thus handle the `FileTooShort` error
    // separately
    let reader = match sniffed_reader {
        Ok(rdr) => Ok(rdr.0),
        Err(niffler::Error::FileTooShort) => {
            Ok(Box::new(BufReader::new(File::open(&path)?)) as Box<dyn Read>)
        }
        Err(e) => Err(e).context(format!("Failed reading FASTA file \"{}\"", path)),
    }?;

    let fasta_reader = fasta::Reader::new(reader);

    // add a more specific error context to each record
    let records = fasta_reader
        .records()
        .map(|item| item.context("Could not parse FASTA record"));

    Ok(records)
}
// let record = record.context(format!("Failed parsing FASTA file \"{}\"", fasta_fn))?;