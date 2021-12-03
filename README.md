# code_utilities
Various scripts for common tasks, involving ones for DNA sequence handling.

### [abi_handling.py](https://github.com/liberjul/code_utilities/blob/main/abi_handling.py)

- Functions for use in cutting chromatograms by Ns (`gen_cut_fastas`) or by PHRED score(`gen_cut_fastas_phred`)
- Renaming using a metadata file (`rename_seqs`)
- Cutting BLAST record headers (`bin_nomen`)
- BLASTing sequences - `BLAST_it` and `batch_blast`

### [auto_sanger_seq_assem.py](https://github.com/liberjul/code_utilities/blob/main/auto_sanger_seq_assem.py)

- Functions for use in assembling reads from overlapping Sanger sequencing runs - `reverse_complement`, `seq_kmers`, `align_w_match`, `assem_seqs`

### [average_quality_fastq.py](https://github.com/liberjul/code_utilities/blob/main/average_quality_fastq.py)

- Script to determine the average PHRED score of bases in a FASTQ file, using PHRED+33

### [consensus_from_aln.py](https://github.com/liberjul/code_utilities/blob/main/consensus_from_aln.py)

- Script to determine the consensus sequence with degenerate bases for an alignment in FASTA format. Useful for primer design.

### [fasta_filename_to_header.py](https://github.com/liberjul/code_utilities/blob/main/fasta_filename_to_header.py)

- Script which inputs a FASTA file, and outputs a FASTA file with identical sequence but with its header changed to match the filename.

### [protein_fasta_from_kbase_prokka.py](https://github.com/liberjul/code_utilities/blob/main/protein_fasta_from_kbase_prokka.py)

- Script which inputs a genbank and GFF file from KBASE's prokka output, and creates a protein FASTA file as output

### [check_grnas_v_db.py](https://github.com/liberjul/code_utilities/blob/main/check_grnas_v_db.py)

- Script which inputs a FASTA file of gRNAS, a FASTA file of sequences potentially targeted by gRNAs (database), and CRISPR match parameters, then outputs a CSV file with of potential targets in the database.

### [check_grnas_v_db_par.py](https://github.com/liberjul/code_utilities/blob/main/check_grnas_v_db_par.py)

- Script which inputs a FASTA file of gRNAS, a FASTA file of sequences potentially targeted by gRNAs (database), and CRISPR match parameters, then outputs a CSV file with of potential targets in the database. Parallelized approach compared to [check_grnas_v_db.py](https://github.com/liberjul/code_utilities/blob/main/check_grnas_v_db.py), but speed increases are not yet substantial.

### [check_grnas_v_db_usearch.py](https://github.com/liberjul/code_utilities/blob/main/check_grnas_v_db_usearch.py)

- Script which inputs a FASTA file of gRNAS, a FASTA file of sequences potentially targeted by gRNAs (database), and CRISPR match parameters, then outputs a CSV file with of potential targets in the database. Uses USEARCH algorithm for finding matches. Very fast, preferred approach.

### [seq_from_filename.py](https://github.com/liberjul/code_utilities/blob/main/seq_from_filename.py)

- Function (`seq_from_filename`) which return a string sequence from a input FASTA filename.

### [seq_to_filename.py](https://github.com/liberjul/code_utilities/blob/main/seq_to_filename.py)

- Function (`seq_to_filename`) which takes a header string, a sequence string, and a filename then creates a corresponsing FASTA file.

### [update_citation_manager.py](https://github.com/liberjul/code_utilities/blob/main/update_citation_manager.py)

- Script to take BibTex files from two citation managers, identify the titles in the first but not second, and produce and importable BibTex file for updating the second database.

### [flatten_fasta.py](https://github.com/liberjul/code_utilities/blob/main/flatten_fasta.py)

- Script which takes a directory with `-d`/`--dir` and takes all FASTA files (`.fasta`) and converts their sequences to be one line per record. If file(s) is(are) large and order of records does not matter, `-i`/`--ignore_order` may improve speed. Does not alter headers or sequence, only formatting. Useful for copy-paste of FASTA contents.

### [process_chromatograms_eurofins.py](https://github.com/liberjul/code_utilities/blob/main/process_chromatograms_eurofins.py)

- Script which takes a directory with `-d`/`--dir`, finds all chromatogram files (`.ab1`) and a metadata file ending in `.csv` (see example [here](https://github.com/liberjul/code_utilities/blob/main/example_metadata.csv)) within that directory. Auto-trims chromatogram sequences, exports to FASTA files, and the chromatogram and FASTA files have names corresponding to the "Name" column in the CSV. "Barcode" and "Name" columns are required in the chromatogram.

Example usage:
```
python process_chromatograms_eurofins.py --dir ./2021_12_03_1954440-1/
```
