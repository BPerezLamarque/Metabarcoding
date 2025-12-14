#!/usr/bin/python3

import sys

if len(sys.argv) != 4:
    sys.stderr.write(
        "Usage: python3 chain_uc.py UNOISE_hits.uc clusters_OTU.uc output.uc\n"
    )
    sys.exit(1)

UNOISEFileName = sys.argv[1]
ClustersFileName = sys.argv[2]
OutputFileName = sys.argv[3]

# Helper function to strip ";size=XXXXX"
def strip_size(name):
    return name.split(";")[0] if name else name

# --------------------------------------------------
# Step 1: Read clusters_OTU.uc
# Build mapping: UNOISE_cluster -> OTU_cluster (size removed)
# --------------------------------------------------
UNOISE_to_OTU = {}

with open(ClustersFileName) as ClustersFile:
    LineNr = 0
    for Line in ClustersFile:
        LineNr += 1
        Line = Line.strip()
        if not Line or Line[0] == '#':
            continue

        Fields = Line.split("\t")
        if len(Fields) < 10:
            sys.stderr.write(
                f"Line {LineNr} in clusters_OTU.uc has < 10 fields\n"
            )
            sys.exit(1)

        cluster = strip_size(Fields[8])
        second_cluster = strip_size(Fields[9]) if Fields[9] != '*' else cluster

        UNOISE_to_OTU[cluster] = second_cluster

# --------------------------------------------------
# Step 2: Process UNOISE_hits.uc
# --------------------------------------------------
with open(UNOISEFileName) as UNOISEFile, open(OutputFileName, "w") as OutputFile:
    LineNr = 0
    for Line in UNOISEFile:
        LineNr += 1
        Line = Line.strip()
        if not Line or Line[0] == '#':
            OutputFile.write(Line + "\n")
            continue

        Fields = Line.split("\t")
        if len(Fields) < 10:
            sys.stderr.write(
                f"Line {LineNr} in UNOISE_hits.uc has < 10 fields\n"
            )
            sys.exit(1)

        # Strip size from read and first cluster
        read_name = strip_size(Fields[8])
        unoise_cluster = strip_size(Fields[9]) if Fields[9] != '*' else read_name

        # Lookup second cluster
        second_cluster = UNOISE_to_OTU.get(unoise_cluster, unoise_cluster)

        # If the second cluster is the read itself, write '*'
        Fields[8] = read_name
        Fields[9] = '*' if second_cluster == read_name else second_cluster

        OutputFile.write("\t".join(Fields) + "\n")


