
Downloading, Processing, and Concatenating *S. cerevisiae* and *S. pombe* Genome Files
======================================================================================

**Supporting code and documentation for the *Bio-protocol* manuscript: "ChIP-seq Data Processing and Relative and Quantitative Signal Normalization for *Saccharomyces cerevisiae*."**

**Author:** *Kris Alavattam*

This notebook provides a step-by-step workflow for downloading, processing, and concatenating *Saccharomyces cerevisiae* and *Schizosaccharomyces pombe* FASTA and GFF3 genome files. The processed genome files are used for read alignment in the ChIP-seq data processing workflow described in the manuscript, while the GFF3 file supports annotation and interpretation of ChIP-seq signal visualizations.

*Note: Each code block is designed to be executed sequentially, as some variables and processed files generated in earlier steps are required for later steps.*

*Note: If using a high-performance computing cluster (HPCC), request an interactive node to ensure adequate resources for running code in the below chunks. The specific command (e.g., `grabnode` at Fred Hutch Cancer Center) will depend on the job scheduler setup. (This step is unnecessary if running the code on a local machine.)*

## Table of contents
<details>
<summary><i>Table of contents</i></summary>
<br />
<!-- MarkdownTOC -->

1. [1. Establish directory structure for genome files](#1-establish-directory-structure-for-genome-files)
1. [2. Download, extract, and organize *S. cerevisiae* genome files](#2-download-extract-and-organize-s-cerevisiae-genome-files)
1. [3. Download, extract, and organize *S. pombe* files](#3-download-extract-and-organize-s-pombe-files)
1. [4. Prepare *S. cerevisiae* FASTA file for concatenation](#4-prepare-s-cerevisiae-fasta-file-for-concatenation)
1. [5. Prepare *S. cerevisiae* GFF3 file for concatenation](#5-prepare-s-cerevisiae-gff3-file-for-concatenation)
1. [6. Prepare *S. pombe* FASTA file for concatenation](#6-prepare-s-pombe-fasta-file-for-concatenation)
1. [7. Prepare *S. pombe* GFF3 file for concatenation](#7-prepare-s-pombe-gff3-file-for-concatenation)
1. [8. Concatenate the processed FASTA and GFF3 files](#8-concatenate-the-processed-fasta-and-gff3-files)

<!-- /MarkdownTOC -->
</details>
<br />

<a id="1-establish-directory-structure-for-genome-files"></a>
## 1. Establish directory structure for genome files
<details>
<summary><i>Text: Establish directory structure for genome files.</i></summary>
<br />

This step sets up the directory structure for organizing *S. cerevisiae* and *S. pombe* genome files. Separate subdirectories are created for `raw` and processed (`proc`) FASTA and GFF3 files, along with a directory for concatenated genome files (`concat`). A temporary directory (`tmp`) is also included for intermediate processing steps.
</details>
<br />

<details>
<summary><i>Bash code: Establish directory structure for genome files.</i></summary>

```bash
#!/bin/bash

#  Optional: Request an interactive node
# grabnode  ## Uncomment to request 1 core, 20 GB memory, 1 day, no GPU ##

#  Define path variables
dir_bas="${HOME}/repos"  ## WARNING: Change as needed ##
dir_rep="${dir_bas}/protocol_chipseq_signal_norm"
dir_dat="${dir_rep}/data"
dir_gen="${dir_dat}/genomes"

# shellcheck disable=SC2086
{
    #  Create subdirectories for storing FASTA and GFF3 files
    mkdir -p ${dir_gen}/{cerevisiae,pombe}/{fasta,gff3}/{raw,proc}
    mkdir -p ${dir_gen}/concat/{fasta,gff3}/proc

    #  Create a temporary directory for intermediate files
    mkdir -p ${dir_gen}/cerevisiae/tmp
}
```
</details>
<br />

<a id="2-download-extract-and-organize-s-cerevisiae-genome-files"></a>
## 2. Download, extract, and organize *S. cerevisiae* genome files
<details>
<summary><i>Text: Download, extract, and organize </i>S. cerevisiae<i> genome files.</i></summary>
<br />

This step downloads the *S. cerevisiae* genome and annotation files from the [Saccharomyces Genome Database (SGD)](https://www.yeastgenome.org/). The genome files are provided as a compressed tarball&mdash;a compressed archive (`.tgz`) containing multiple related files (more information [here](https://en.wikipedia.org/wiki/Tar_(computing))). After downloading, the tarball is extracted in the temporary directory, `tmp`. The FASTA and GFF3 files are then moved to their respective `raw` data directories for further processing. To conserve space, `tmp` is removed after extraction.
</details>
<br />

<details>
<summary><i>Bash code: Download, extract, and organize </i>S. cerevisiae<i> genome files.</i></summary>

```bash
#!/bin/bash

#  Define S. cerevisiae URL, tarball, and file names
lnk_sc_fa_1="http://sgd-archive.yeastgenome.org/sequence"
lnk_sc_fa_2="S288C_reference/genome_releases"
tarball="S288C_reference_genome_R64-5-1_20240529.tgz"
fil_sc_fa="S288C_reference_sequence_R64-5-1_20240529"
unpack="${fil_sc_fa/sequence/genome}"
fil_sc_g3="saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz"

#  Download the tarball to the S. cerevisiae directory
curl \
    -o "${dir_gen}/cerevisiae/${tarball}" \
    "${lnk_sc_fa_1}/${lnk_sc_fa_2}/${tarball}"

#  Unpack the tarball in the temporary directory
tar \
    -xzf "${dir_gen}/cerevisiae/${tarball}" \
    -C "${dir_gen}/cerevisiae/tmp"

#  Move the FASTA and GFF3 files to corresponding raw directories
mv \
    "${dir_gen}/cerevisiae/tmp/${unpack}/${fil_sc_fa}.fsa.gz" \
    "${dir_gen}/cerevisiae/fasta/raw/"

mv \
    "${dir_gen}/cerevisiae/tmp/${unpack}/${fil_sc_g3}" \
    "${dir_gen}/cerevisiae/gff3/raw/"

#  Clean up the temporary directory (but retain the original tarball)
rm -rf "${dir_gen}/cerevisiae/tmp"
```
</details>
<br />

<a id="3-download-extract-and-organize-s-pombe-files"></a>
## 3. Download, extract, and organize *S. pombe* files
<details>
<summary><i>Text: Download, extract, and organize </i>S. pombe<i> files</i></summary>
<br />

This step downloads the *S. pombe* genome and annotation files from [Pombase](https://www.pombase.org/). Unlike *S. cerevisiae*, these files are available individually, so they are downloaded directly as compressed FASTA and GFF3 files. The downloaded files are stored in their respective `raw` data directories for further processing.
</details>
<br />

<details>
<summary><i>Bash code: Download, extract, and organize </i>S. pombe<i> files</i></summary>

```bash
#!/bin/bash

#  Define the S. pombe URL and file names
lnk_sp_1="https://www.pombase.org/data"
lnk_sp_2="releases/pombase-2024-11-01"
fil_sp_fa="Schizosaccharomyces_pombe_all_chromosomes.fa.gz"
fil_sp_g3="Schizosaccharomyces_pombe_all_chromosomes.gff3.gz"

#  Download and store the S. pombe FASTA file
curl \
    -o "${dir_gen}/pombe/fasta/raw/${fil_sp_fa}" \
    "${lnk_sp_1}/${lnk_sp_2}/fasta/chromosomes/${fil_sp_fa}"

#  Download and store the S. pombe GFF3 file
curl "${lnk_sp_1}/${lnk_sp_2}/gff/${fil_sp_g3%.gz}" \
    | gzip \
        > "${dir_gen}/pombe/gff3/raw/${fil_sp_g3}"
```
</details>
<br />

<a id="4-prepare-s-cerevisiae-fasta-file-for-concatenation"></a>
## 4. Prepare *S. cerevisiae* FASTA file for concatenation
<details>
<summary><i>Text: Prepare </i>S. cerevisiae<i> FASTA file for concatenation</i></summary>
<br />

This step processes the *S. cerevisiae* FASTA file to standardize chromosome names and simplify headers for downstream analysis. The script removes unnecessary prefixes (e.g., "chr"), renames the mitochondrial chromosome to "Mito," and retains only relevant sequence information. The processed FASTA file is then compressed and stored in the designated directory for further use in genome concatenation.
</details>
<br />

<details>
<summary><i>Bash code: Prepare </i>S. cerevisiae<i> FASTA file for concatenation</i></summary>

```bash
#!/bin/bash

#  Define directories and file names for raw and processed FASTA files
dir_sc_fa_un="${dir_gen}/cerevisiae/fasta/raw"
dir_sc_fa_pr="${dir_gen}/cerevisiae/fasta/proc"
fil_sc_fa_un="S288C_reference_sequence_R64-5-1_20240529.fsa.gz"
fil_sc_fa_pr="S288C_R64-5-1_proc.fasta.gz"

#  Process raw S. cerevisiae FASTA file
zcat "${dir_sc_fa_un}/${fil_sc_fa_un}" \
    | awk '
        #  Set input field separator (FS) to either "chromosome=" or
        #+ "location=", and no output field separator (OFS)
        BEGIN { FS="chromosome=|location="; OFS="" }

        #  Do main processing
        {
            #  Find header lines, which start with ">", and process them
            if ($0 ~ /^>/) {
                if ($2 ~ /mitochondrion/) {
                    #  Rename mitochondrial chromosome
                    print ">Mito"
                } else {
                    #  Extract chromosome names, printing them without "chr"
                    #+ prefixes
                    split($2, nam_chr, "]")
                    print ">" nam_chr[1]
                }
            } else {
                # Print sequence lines as they are
                print $0
            }
        }
    ' \
    | gzip \
        > "${dir_sc_fa_pr}/${fil_sc_fa_pr}"
```
</details>
<br />

<a id="5-prepare-s-cerevisiae-gff3-file-for-concatenation"></a>
## 5. Prepare *S. cerevisiae* GFF3 file for concatenation
<details>
<summary><i>Text: Prepare </i>S. cerevisiae<i> GFF3 file for concatenation</i></summary>
<br />

This step processes the *S. cerevisiae* GFF3 file to standardize chromosome names and convert systematic gene names to standard names to improve readability. The [AWK](https://en.wikipedia.org/wiki/AWK) scripting removes unnecessary prefixes (e.g., "chr"), renames "mt" to "Mito," and decodes URL-encoded characters and HTML entities. Additionally, it ensures that gene and ARS names follow standard conventions for increased legibility in downstream analysis and visualization. The processed GFF3 file is then compressed and stored for later use in genome concatenation.
</details>
<br />

<details>
<summary><i>Bash code: Prepare </i>S. cerevisiae<i> GFF3 file for concatenation</i></summary>

```bash
#!/bin/bash

#  Define directories and file names for raw and processed S. cerevisiae GFF3
#+ files
dir_sc_g3_un="${dir_gen}/cerevisiae/gff3/raw"
dir_sc_g3_pr="${dir_gen}/cerevisiae/gff3/proc"
fil_sc_g3_un="saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz"
fil_sc_g3_pr="S288C_R64-5-1_proc.gff3.gz"
fil_sc_g3_re="S288C_R64-5-1_proc_readable.gff3.gz"

#  Process raw S. cerevisiae GFF3 file
zcat "${dir_sc_g3_un}/${fil_sc_g3_un}" \
    | awk -F "\t" '
        #  Set in- and output field separators to tabs ("\t")
        BEGIN { OFS = FS; skip = 0 }

        #  Do main processing
        {
            #  Skip "###" line and all lines after it
            if ($0 ~ /^###/) { skip = 1 }
            if (skip == 1) { next }

            #  In column 1 ($1), remove "chr" prefix and rename "mt" to "Mito"
            gsub(/^chr/, "", $1)
            gsub(/^mt/, "Mito", $1)

            #  Convert URL-encoded characters in $9
            gsub(/%20/, " ", $9)
            gsub(/%2C/, ",", $9)
            gsub(/%3B/, ",", $9)  #TODO: Write a general note
            gsub(/%28/, "(", $9)
            gsub(/%29/, ")", $9)

            #  Convert HTML entities in $9
            gsub(/&#946/, "beta", $9)
            gsub(/&#8242/, " prime", $9)

            #  Convert literal percent symbol to word "percent"
            gsub(/%/, " percent", $9)

            #  In $7, convert "0" to "."
            if ($7 == "0") { $7 = "." }

            print $0
        }
    ' \
    | gzip \
        > "${dir_sc_g3_pr}/${fil_sc_g3_pr}"

#  Process GFF3 file for readability by modifying gene and ARS names
zcat "${dir_sc_g3_pr}/${fil_sc_g3_pr}" \
    | awk '
        BEGIN { FS = OFS = "\t" }

        #  Skip lines starting with a hash, as these are headers or comments
        $1 ~ /^#/ { print; next }

        #  Process lines where column 3 ($3) matches gene-related features
        $3 ~ /^(\
            blocked_reading_frame|gene|ncRNA_gene|\
            pseudogene|rRNA_gene|snRNA_gene|\
            snoRNA_gene|tRNA_gene|telomerase_RNA_gene\
        )$/ {
            #  If a "gene=" field exists in $9, extract the gene name
            if (match($9, /gene=[^;]+/)) {
                nam_gen = substr($9, RSTART + 5, RLENGTH - 5)
                
                #  Replace the "Name=" field in $9 with the extracted gene
                #+ name
                sub(/Name=[^;]+/, "Name=" nam_gen, $9)
            }
        }

        #  Process lines where $3 is "ARS"
        $3 == "ARS" {
            #  If an "Alias=" field exists in $9, extract the alias name
            if (match($9, /Alias=[^;]+/)) {
                nam_als = substr($9, RSTART + 6, RLENGTH - 6)

                #  Replace the "Name=" field in $9 with the extracted alias
                #+ name
                sub(/Name=[^;]+/, "Name=" nam_als, $9)
            }
        }

        #  Print all lines; modified lines are printed after processing
        { print }
    ' \
    | gzip \
        > "${dir_sc_g3_pr}/${fil_sc_g3_re}"
```
</details>
<br />

<a id="6-prepare-s-pombe-fasta-file-for-concatenation"></a>
## 6. Prepare *S. pombe* FASTA file for concatenation
<details>
<summary><i>Text: Prepare </i>S. pombe<i> FASTA file for concatenation</i></summary>
<br />

This step standardizes *S. pombe* chromosome names by prefixing them with "SP_" and simplifying longer names. These modifications ensure that *S. pombe* alignments are distinguishable from *S. cerevisiae* alignments after Bowtie2 mapping&mdash;which is essential for computing spike-in scaling factors.
</details>
<br />

<details>
<summary><i>Bash code: Prepare </i>S. pombe<i> FASTA file for concatenation</i></summary>

```bash
#!/bin/bash

#  Define directories and file names for raw and processed S. pombe FASTA files
dir_sp_fa_un="${dir_gen}/pombe/fasta/raw"
dir_sp_fa_pr="${dir_gen}/pombe/fasta/proc"
fil_sp_fa_un="Schizosaccharomyces_pombe_all_chromosomes.fa.gz"
fil_sp_fa_pr="972h-_2024-11-01_proc.fasta.gz"

#  Prepend substring "SP_" to chromosome names and simplify long chromosome/DNA
#+ names
zcat "${dir_sp_fa_un}/${fil_sp_fa_un}" \
    | sed -r '
        s/^>chr_II_telomeric_gap\ .*$/>SP_II_TG/g;
        s/^>I\ .*$/>SP_I/g;
        s/^>II\ .*$/>SP_II/g;
        s/^>III\ .*$/>SP_III/g;
        s/^>mating_type_region\ .*$/>SP_MTR/g;
        s/^>mitochondrial\ .*$/>SP_Mito/g
    ' \
    | gzip \
        > "${dir_sp_fa_pr}/${fil_sp_fa_pr}"
```
</details>
<br />

<a id="7-prepare-s-pombe-gff3-file-for-concatenation"></a>
## 7. Prepare *S. pombe* GFF3 file for concatenation
<details>
<summary><i>Text: Prepare </i>S. pombe<i> GFF3 file for concatenation</i></summary>
<br />

This step processes the *S. pombe* GFF3 file by standardizing chromosome names, prefixing them with "SP_" to distinguish them from *S. cerevisiae* chromosomes, and simplifying chromosome and DNA region names. These modifications ensure consistency with the processed FASTA file and facilitate downstream parsing and annotation.
</details>
<br />

<details>
<summary><i>Bash code: Prepare </i>S. pombe<i> GFF3 file for concatenation</i></summary>

```bash
#!/bin/bash

#  Define directories and file names for raw and processed S. pombe GFF3 files
dir_sp_g3_un="${dir_gen}/pombe/gff3/raw"
dir_sp_g3_pr="${dir_gen}/pombe/gff3/proc"
fil_sp_g3_un="Schizosaccharomyces_pombe_all_chromosomes.gff3.gz"
fil_sp_g3_pr="972h-_2024-11-01_proc.gff3.gz"

#  Prepend substring "SP_" to chromosome names and simplify long chromosome/DNA
#+ names
zcat "${dir_sp_g3_un}/${fil_sp_g3_un}" \
    | sed '
        s/^chr_II_telomeric_gap/SP_II_TG/g;
        s/^I/SP_I/g;
        s/^II/SP_II/g;
        s/^III/SP_III/g;
        s/^mating_type_region/SP_MTR/g;
        s/^mitochondrial/SP_Mito/g
    ' \
    | gzip \
        > "${dir_sp_g3_pr}/${fil_sp_g3_pr}"
```
</details>
<br />

<a id="8-concatenate-the-processed-fasta-and-gff3-files"></a>
## 8. Concatenate the processed FASTA and GFF3 files
<details>
<summary><i>Text: Concatenate the processed FASTA and GFF3 files</i></summary>
<br />

This final step concatenates the processed *S. cerevisiae* and *S. pombe* FASTA and GFF3 files into unified files for downstream data processing and analyses. The concatenated FASTA file supports Bowtie2 index generation, alignment, and IGV visualization, while the combined GFF3 file provides consistent, interpretable genome annotations for IGV.
</details>
<br />

<details>
<summary><i>Bash code: Concatenate the processed FASTA and GFF3 files</i></summary>

```bash
#!/bin/bash

#  Define variables for processed input and output FASTA and GFF3 files
sc="S288C_R64-5-1_proc"
sc_fa="${sc}.fasta.gz"
sc_g3="${sc}.gff3.gz"
sc_re="${sc}_readable.gff3.gz"

sp="972h-_2024-11-01_proc"
sp_fa="${sp}.fasta.gz"
sp_g3="${sp}.gff3.gz"
sp_re="${sp_g3}"

cc="sc_sp_proc"
cc_fa="${cc}.fasta.gz"
cc_g3="${cc}.gff3.gz"
cc_re="${cc}_readable.gff3.gz"

#  Concatenate processed FASTA files for S. cerevisiae and S. pombe
cat \
    "${dir_gen}/cerevisiae/fasta/proc/${sc_fa}" \
    "${dir_gen}/pombe/fasta/proc/${sp_fa}" \
        > "${dir_gen}/concat/fasta/proc/${cc_fa}"

#  Create uncompressed version of concatenated FASTA file for use in Bowtie2
#+ index generation
gunzip -c "${dir_gen}/concat/fasta/proc/${cc_fa}" \
    > "${dir_gen}/concat/fasta/proc/${cc_fa%.gz}"

#  Concatenate processed GFF3 files for S. cerevisiae and S. pombe
cat \
    "${dir_gen}/cerevisiae/gff3/proc/${sc_g3}" \
    "${dir_gen}/pombe/gff3/proc/${sp_g3}" \
        > "${dir_gen}/concat/gff3/proc/${cc_g3}"

#  Concatenate processed, easier-to-read GFF3 files
cat \
    "${dir_gen}/cerevisiae/gff3/proc/${sc_re}" \
    "${dir_gen}/pombe/gff3/proc/${sp_re}" \
        > "${dir_gen}/concat/gff3/proc/${cc_re}"
```
</details>
<br />
