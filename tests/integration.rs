use std::collections::HashMap;
use std::path::Path;
use std::process::Command;

/// Parse tab-separated kmer counts (KMER\tCOUNT) into a HashMap.
fn parse_kmer_counts(output: &str) -> HashMap<String, usize> {
    output
        .lines()
        .filter(|line| !line.is_empty())
        .map(|line| {
            let mut parts = line.split('\t');
            let kmer = parts.next().unwrap().to_string();
            let count: usize = parts.next().unwrap().parse().unwrap();
            (kmer, count)
        })
        .collect()
}

/// Read a KMC reference fixture file into a HashMap.
fn load_kmc_reference(fixture_path: &str) -> HashMap<String, usize> {
    let content = std::fs::read_to_string(fixture_path)
        .unwrap_or_else(|e| panic!("Failed to read fixture {}: {}", fixture_path, e));
    parse_kmer_counts(&content)
}

/// Run smakcr with the given arguments and return stdout as a string.
fn run_smakcr(args: &[&str]) -> String {
    let binary = env!("CARGO_BIN_EXE_smakcr");
    let output = Command::new(binary)
        .args(args)
        .output()
        .expect("Failed to execute smakcr");

    assert!(
        output.status.success(),
        "smakcr failed with stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    String::from_utf8(output.stdout).expect("smakcr output is not valid UTF-8")
}

/// Compare smakcr output against a KMC reference fixture.
/// Only compares k-mers with count > 0 (since smakcr without -z skips zeros, matching KMC).
fn assert_matches_kmc(smakcr_output: &str, fixture_path: &str) {
    let smakcr_counts = parse_kmer_counts(smakcr_output);
    let kmc_counts = load_kmc_reference(fixture_path);

    // Check that every k-mer in KMC output is in smakcr output with the same count.
    for (kmer, &kmc_count) in &kmc_counts {
        let smakcr_count = smakcr_counts.get(kmer).copied().unwrap_or(0);
        assert_eq!(
            smakcr_count, kmc_count,
            "Mismatch for k-mer {}: smakcr={}, kmc={}",
            kmer, smakcr_count, kmc_count
        );
    }

    // Check that smakcr doesn't have extra k-mers that KMC doesn't.
    for (kmer, &smakcr_count) in &smakcr_counts {
        let kmc_count = kmc_counts.get(kmer).copied().unwrap_or(0);
        assert_eq!(
            smakcr_count, kmc_count,
            "Extra k-mer in smakcr output {}: smakcr={}, kmc={}",
            kmer, smakcr_count, kmc_count
        );
    }

    // Verify total counts match.
    let smakcr_total: usize = smakcr_counts.values().sum();
    let kmc_total: usize = kmc_counts.values().sum();
    assert_eq!(
        smakcr_total, kmc_total,
        "Total count mismatch: smakcr={}, kmc={}",
        smakcr_total, kmc_total
    );
}

fn fixtures_dir() -> &'static str {
    concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures")
}

fn test_fna_path() -> &'static str {
    concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/test.fna")
}

fn chrm_path() -> &'static str {
    concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/chrM.fa.gz")
}

fn chrun_combined_path() -> &'static str {
    concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/chrUn_combined.fa.gz")
}

fn chrun_per_scaffold_dir() -> &'static str {
    concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/chrUn-one-seq-per-file")
}

/// Collect all .fa.gz files in the chrUn per-scaffold directory as owned Strings.
fn chrun_per_scaffold_files() -> Vec<String> {
    let mut files: Vec<String> = std::fs::read_dir(chrun_per_scaffold_dir())
        .unwrap()
        .filter_map(|e| {
            let path = e.unwrap().path();
            if path.extension().and_then(|s| s.to_str()) == Some("gz") {
                Some(path.to_str().unwrap().to_string())
            } else {
                None
            }
        })
        .collect();
    files.sort();
    files
}

// =====================================================================
// test.fna — non-canonical
// =====================================================================

#[test]
fn test_fna_noncanonical_k3() {
    let output = run_smakcr(&[test_fna_path(), "-k", "3"]);
    assert_matches_kmc(&output, &format!("{}/test_fna_nc_k3_kmc.tsv", fixtures_dir()));
}

#[test]
fn test_fna_noncanonical_k5() {
    let output = run_smakcr(&[test_fna_path(), "-k", "5"]);
    assert_matches_kmc(&output, &format!("{}/test_fna_nc_k5_kmc.tsv", fixtures_dir()));
}

#[test]
fn test_fna_noncanonical_k7() {
    let output = run_smakcr(&[test_fna_path(), "-k", "7"]);
    assert_matches_kmc(&output, &format!("{}/test_fna_nc_k7_kmc.tsv", fixtures_dir()));
}

// =====================================================================
// test.fna — canonical
// =====================================================================

#[test]
fn test_fna_canonical_k3() {
    let output = run_smakcr(&[test_fna_path(), "-k", "3", "-c"]);
    assert_matches_kmc(&output, &format!("{}/test_fna_c_k3_kmc.tsv", fixtures_dir()));
}

#[test]
fn test_fna_canonical_k5() {
    let output = run_smakcr(&[test_fna_path(), "-k", "5", "-c"]);
    assert_matches_kmc(&output, &format!("{}/test_fna_c_k5_kmc.tsv", fixtures_dir()));
}

#[test]
fn test_fna_canonical_k7() {
    let output = run_smakcr(&[test_fna_path(), "-k", "7", "-c"]);
    assert_matches_kmc(&output, &format!("{}/test_fna_c_k7_kmc.tsv", fixtures_dir()));
}

// =====================================================================
// chrM — non-canonical
// =====================================================================

#[test]
fn chrm_noncanonical_k3() {
    let output = run_smakcr(&[chrm_path(), "-k", "3"]);
    assert_matches_kmc(&output, &format!("{}/chrM_nc_k3_kmc.tsv", fixtures_dir()));
}

#[test]
fn chrm_noncanonical_k5() {
    let output = run_smakcr(&[chrm_path(), "-k", "5"]);
    assert_matches_kmc(&output, &format!("{}/chrM_nc_k5_kmc.tsv", fixtures_dir()));
}

#[test]
fn chrm_noncanonical_k7() {
    let output = run_smakcr(&[chrm_path(), "-k", "7"]);
    assert_matches_kmc(&output, &format!("{}/chrM_nc_k7_kmc.tsv", fixtures_dir()));
}

// =====================================================================
// chrM — canonical
// =====================================================================

#[test]
fn chrm_canonical_k3() {
    let output = run_smakcr(&[chrm_path(), "-k", "3", "-c"]);
    assert_matches_kmc(&output, &format!("{}/chrM_c_k3_kmc.tsv", fixtures_dir()));
}

#[test]
fn chrm_canonical_k5() {
    let output = run_smakcr(&[chrm_path(), "-k", "5", "-c"]);
    assert_matches_kmc(&output, &format!("{}/chrM_c_k5_kmc.tsv", fixtures_dir()));
}

#[test]
fn chrm_canonical_k7() {
    let output = run_smakcr(&[chrm_path(), "-k", "7", "-c"]);
    assert_matches_kmc(&output, &format!("{}/chrM_c_k7_kmc.tsv", fixtures_dir()));
}

// =====================================================================
// Multi-threaded tests — single file with -t (seq_io::parallel path)
// Uses chrUn_all.fa.gz (39 sequences) to actually exercise parallelism.
// =====================================================================

#[test]
fn chrun_noncanonical_k3_multithreaded() {
    let output = run_smakcr(&[chrun_combined_path(), "-k", "3", "-t", "4"]);
    assert_matches_kmc(&output, &format!("{}/chrUn_nc_k3_kmc.tsv", fixtures_dir()));
}

#[test]
fn chrun_noncanonical_k5_multithreaded() {
    let output = run_smakcr(&[chrun_combined_path(), "-k", "5", "-t", "4"]);
    assert_matches_kmc(&output, &format!("{}/chrUn_nc_k5_kmc.tsv", fixtures_dir()));
}

#[test]
fn chrun_canonical_k5_multithreaded() {
    let output = run_smakcr(&[chrun_combined_path(), "-k", "5", "-t", "4", "-c"]);
    assert_matches_kmc(&output, &format!("{}/chrUn_c_k5_kmc.tsv", fixtures_dir()));
}

#[test]
fn chrun_noncanonical_k7_multithreaded() {
    let output = run_smakcr(&[chrun_combined_path(), "-k", "7", "-t", "2"]);
    assert_matches_kmc(&output, &format!("{}/chrUn_nc_k7_kmc.tsv", fixtures_dir()));
}

#[test]
fn chrun_canonical_k7_multithreaded() {
    let output = run_smakcr(&[chrun_combined_path(), "-k", "7", "-t", "4", "-c"]);
    assert_matches_kmc(&output, &format!("{}/chrUn_c_k7_kmc.tsv", fixtures_dir()));
}

// Also test test.fna multi-threaded (edge cases under threading)
#[test]
fn test_fna_noncanonical_k3_multithreaded() {
    let output = run_smakcr(&[test_fna_path(), "-k", "3", "-t", "4"]);
    assert_matches_kmc(&output, &format!("{}/test_fna_nc_k3_kmc.tsv", fixtures_dir()));
}

// =====================================================================
// Multi-file + multi-threaded tests (rayon path)
// Uses chrUn-one-seq-per-file/ (39 separate .fa.gz files) and compares
// against the same KMC references as the combined file.
// =====================================================================

#[test]
fn multi_file_noncanonical_k3_multithreaded() {
    let files = chrun_per_scaffold_files();
    let mut args: Vec<&str> = files.iter().map(|s| s.as_str()).collect();
    args.extend(["-k", "3", "-t", "4"]);

    let output = run_smakcr(&args);
    assert_matches_kmc(&output, &format!("{}/chrUn_nc_k3_kmc.tsv", fixtures_dir()));
}

#[test]
fn multi_file_noncanonical_k5_multithreaded() {
    let files = chrun_per_scaffold_files();
    let mut args: Vec<&str> = files.iter().map(|s| s.as_str()).collect();
    args.extend(["-k", "5", "-t", "4"]);

    let output = run_smakcr(&args);
    assert_matches_kmc(&output, &format!("{}/chrUn_nc_k5_kmc.tsv", fixtures_dir()));
}

#[test]
fn multi_file_canonical_k5_multithreaded() {
    let files = chrun_per_scaffold_files();
    let mut args: Vec<&str> = files.iter().map(|s| s.as_str()).collect();
    args.extend(["-k", "5", "-t", "4", "-c"]);

    let output = run_smakcr(&args);
    assert_matches_kmc(&output, &format!("{}/chrUn_c_k5_kmc.tsv", fixtures_dir()));
}

#[test]
fn multi_file_canonical_k7_multithreaded() {
    let files = chrun_per_scaffold_files();
    let mut args: Vec<&str> = files.iter().map(|s| s.as_str()).collect();
    args.extend(["-k", "7", "-t", "4", "-c"]);

    let output = run_smakcr(&args);
    assert_matches_kmc(&output, &format!("{}/chrUn_c_k7_kmc.tsv", fixtures_dir()));
}

// =====================================================================
// Output file tests
// =====================================================================

#[test]
fn test_output_to_file() {
    let tmp_dir = std::env::temp_dir();
    let outfile = tmp_dir.join("smakcr_test_output.tsv");
    let outfile_str = outfile.to_str().unwrap();

    run_smakcr(&[test_fna_path(), "-k", "3", "-o", outfile_str]);

    assert!(Path::new(outfile_str).exists(), "Output file was not created");
    let content = std::fs::read_to_string(outfile_str).unwrap();
    let counts = parse_kmer_counts(&content);
    assert!(!counts.is_empty(), "Output file is empty");

    // Compare against KMC reference
    let kmc_counts = load_kmc_reference(&format!("{}/test_fna_nc_k3_kmc.tsv", fixtures_dir()));
    assert_eq!(counts, kmc_counts);

    std::fs::remove_file(outfile_str).ok();
}

// =====================================================================
// Zero-count output tests
// =====================================================================

#[test]
fn test_zero_counts_flag() {
    let output = run_smakcr(&[test_fna_path(), "-k", "3", "-z"]);
    let counts = parse_kmer_counts(&output);
    // With -z, we should get all 4^3 = 64 k-mers
    assert_eq!(counts.len(), 64, "Expected 64 k-mers with -z flag, got {}", counts.len());
}

// =====================================================================
// Edge case tests
// =====================================================================

#[test]
fn test_empty_file() {
    let tmp_dir = std::env::temp_dir();
    let empty_file = tmp_dir.join("smakcr_test_empty.fna");
    std::fs::write(&empty_file, "").unwrap();

    let binary = env!("CARGO_BIN_EXE_smakcr");
    let output = Command::new(binary)
        .args([empty_file.to_str().unwrap(), "-k", "3"])
        .output()
        .expect("Failed to execute smakcr");

    // Should fail gracefully on empty file
    assert!(!output.status.success());

    std::fs::remove_file(&empty_file).ok();
}

#[test]
fn test_single_short_sequence() {
    let tmp_dir = std::env::temp_dir();
    let short_file = tmp_dir.join("smakcr_test_short.fna");
    std::fs::write(&short_file, ">seq1\nAC\n").unwrap();

    let output = run_smakcr(&[short_file.to_str().unwrap(), "-k", "3"]);
    let counts = parse_kmer_counts(&output);
    // Sequence "AC" is shorter than k=3, so no k-mers should be counted
    assert!(counts.is_empty(), "Expected no k-mers for sequence shorter than k");

    std::fs::remove_file(&short_file).ok();
}
