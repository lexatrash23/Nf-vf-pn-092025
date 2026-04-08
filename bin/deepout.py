#!/usr/bin/env python3

import sys

# Usage: python3 filter_3line.py input.3line output_fasta_all output_fasta_Oonly
if len(sys.argv) != 3:
    print("Usage: python3 deepoutput.py input.3line mature.fasta")
    sys.exit(1)

input_file = sys.argv[1]
mature_fasta = sys.argv[2]

with open(input_file, 'r') as infile, open(mature_fasta, 'w') as outfile:
    for line_number, line in enumerate(infile, start=1):
        if line_number % 3 == 1:
            header = line.strip()
            parts = header.split("|", 1)
            header_before = parts[0].strip()
            header_after = parts[1].strip() if len(parts) > 1 else ""
            print(header_after)
        elif line_number % 3 == 2:
            sequence = line.strip()
        else:
            label = line.strip()
            if "SP " in header_after:
                mature_seq = "".join(sequence[j] for j in range(len(sequence)) if label[j] == "O")
                if mature_seq:
                    outfile.write(header_before + "\n")
                    outfile.write(mature_seq + "\n")
            
