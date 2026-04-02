#!/bin/bash

#  Function to output a string of FASTQ files such that sample-specific ones
#+ are separated by semicolons and, within semicolon-delimited substrings,
#+ sample-specific FASTQ pairs are separated by commas
function pair_fastqs() {
    awk '
        BEGIN { OFS = "" }
        {
            #  Check for lines ending with "_R1.atria" or simply "_R1" followed
            #+ by ".fastq.gz" or ".fq.gz"
            if ($0 ~ /_R1(\.atria)?\.(fastq|fq)(\.gz)?$/) {
                r1 = $0
                getline
                #  Ensure the next line is a matching "_R2" file
                if ($0 ~ /_R2(\.atria)?\.(fastq|fq)(\.gz)?$/) {
                    r2 = $0
                    print r1 ",", r2, ";"
                } else {
                    #  Print an error message and exit with a failure code if
                    #+ no matching "_R2" file is found
                    print "Error: Missing R2 file for " r1 > "/dev/stderr"
                    exit 1
                }
            } else {
                #  Handle unpaired files
                print $0 ";"
            }
        }
    ' 
}
