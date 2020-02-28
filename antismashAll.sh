# This file is assuming a "run_antismash" script in "~/bin/"
# which will call docker run
# contant of that script is appendended in the end
# valid antismash options (help file) in the end

getAbsFilePath() {
  # $1 : relative filename
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

runAntismash() {
  # $1 : absolute filepath
  local abs=$(getAbsFilePath "$1")
  absout="${abs%.*}.antismash_out/"
  echo
  echo "==============================================================="
  echo
  # add "--genefinding-tool prodigal" in the end if input do not contain annotation
  # add "--taxon fungi" in the end if you use fungi sequence
  ~/bin/run_antismash $abs $absout --cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go --smcog-trees
}

# change this accroding to the downloaded genbank
readonly genbankFileExt=gbff

if [ $# -eq 1 ]; then
    for i in "${1}/"*.$genbankFileExt; do
        [ -f "$i" ] || break
        echo $i
        runAntismash "$i"
    done
else
    echo "Argument error, please give path to your folder containing gbk files as the only argument."
fi


########### run_antismash ###############
#!/bin/bash

# set -o errexit
# set -o nounset

# function realpath() {
#     echo $(readlink -f $1 2>/dev/null || python -c "import sys; import os; print(os.path.realpath(os.path.expanduser(sys.argv[1])))" $1)
# }

# # handle input file
# readonly INPUT_FILE=$(basename $1)
# readonly INPUT_DIR=$(dirname $(realpath $1))
# shift

# # handle output file
# readonly OUTPUT_DIR=$(realpath $1)
# shift

# # Links within the container
# readonly CONTAINER_SRC_DIR=/input
# readonly CONTAINER_DST_DIR=/output

# if [ ! -d ${OUTPUT_DIR} ]; then
#     mkdir ${OUTPUT_DIR}
# fi

# docker run \
#     --volume ${INPUT_DIR}:${CONTAINER_SRC_DIR}:ro \
#     --volume ${OUTPUT_DIR}:${CONTAINER_DST_DIR}:rw \
#     --detach=false \
#     --rm \
#     --user=$(id -u):$(id -g) \
#     antismash/standalone-nonfree:latest \
#     ${INPUT_FILE} \
#     $@
# # --rm for one time useage of the container
# # antismash/standalone-nonfree:latest -- change if your image is not this one



########### antiSMASH 5.1.1 #############

# usage: antismash [--taxon {bacteria,fungi}] [--output-dir OUTPUT_DIR]
#                  [--reuse-results PATH] [--limit LIMIT]
#                  [--minlength MINLENGTH] [--start START] [--end END]
#                  [--databases PATH] [--write-config-file PATH]
#                  [--without-fimo]
#                  [--executable-paths EXECUTABLE:PATH,EXECUTABLE2:PATH2,...]
#                  [-v] [-d] [--logfile PATH] [--list-plugins] [--check-prereqs]
#                  [--limit-to-record RECORD_ID] [-V] [--profiling]
#                  [--skip-sanitisation] [--skip-zip-file] [--minimal]
#                  [--enable-genefunctions] [--enable-tta]
#                  [--enable-lanthipeptides] [--enable-thiopeptides]
#                  [--enable-nrps-pks] [--enable-sactipeptides]
#                  [--enable-lassopeptides] [--enable-t2pks] [--enable-html]
#                  [--genefinding-tool {glimmerhmm,prodigal,prodigal-m,none,error}]
#                  [--genefinding-gff3 GFF3_FILE]
#                  [--hmmdetection-strictness {strict,relaxed,loose}]
#                  [--fullhmmer]
#                  [--fullhmmer-pfamdb-version FULLHMMER_PFAMDB_VERSION]
#                  [--cassis] [--cf-borders-only] [--cf-create-clusters]
#                  [--cf-min-cds INT] [--cf-mean-threshold FLOAT]
#                  [--cf-min-pfams INT] [--clusterhmmer]
#                  [--clusterhmmer-pfamdb-version CLUSTERHMMER_PFAMDB_VERSION]
#                  [--smcog-trees] [--tta-threshold TTA_THRESHOLD]
#                  [--cb-general] [--cb-subclusters] [--cb-knownclusters]
#                  [--cb-nclusters count] [--cb-min-homology-scale LIMIT]
#                  [--asf] [--pfam2go] [--html-title HTML_TITLE]
#                  [--html-description HTML_DESCRIPTION] [-h] [--help-showall]
#                  [-c CPUS]
#                  [SEQUENCE [SEQUENCE ...]]


# arguments:
#   SEQUENCE  GenBank/EMBL/FASTA file(s) containing DNA.

# --------
# Options
# --------
# -h, --help              Show this help text.
# --help-showall          Show full lists of arguments on this help text.
# -c CPUS, --cpus CPUS    How many CPUs to use in parallel. (default: 8)

# Basic analysis options:

#   --taxon {bacteria,fungi}
#                         Taxonomic classification of input sequence. (default:
#                         bacteria)

# Additional analysis:

#   --fullhmmer           Run a whole-genome HMMer analysis.
#   --cassis              Motif based prediction of SM gene cluster regions.
#   --cf-borders-only     Only annotate borders of existing clusters.
#   --cf-create-clusters  Find extra clusters.
#   --clusterhmmer        Run a cluster-limited HMMer analysis.
#   --smcog-trees         Generate phylogenetic trees of sec. met. cluster
#                         orthologous groups.
#   --tta-threshold TTA_THRESHOLD
#                         Lowest GC content to annotate TTA codons at (default:
#                         0.65).
#   --cb-general          Compare identified clusters against a database of
#                         antiSMASH-predicted clusters.
#   --cb-subclusters      Compare identified clusters against known subclusters
#                         responsible for synthesising precursors.
#   --cb-knownclusters    Compare identified clusters against known gene
#                         clusters from the MIBiG database.
#   --asf                 Run active site finder analysis.
#   --pfam2go             Run Pfam to Gene Ontology mapping module.

# Output options:

#   --output-dir OUTPUT_DIR
#                         Directory to write results to.
#   --html-title HTML_TITLE
#                         Custom title for the HTML output page (default is
#                         input filename).
#   --html-description HTML_DESCRIPTION
#                         Custom description to add to the output.

# Advanced options:

#   --reuse-results PATH  Use the previous results from the specified json
#                         datafile
#   --limit LIMIT         Only process the first <limit> records (default: -1).
#                         -1 to disable
#   --minlength MINLENGTH
#                         Only process sequences larger than <minlength>
#                         (default: 1000).
#   --start START         Start analysis at nucleotide specified.
#   --end END             End analysis at nucleotide specified
#   --databases PATH      Root directory of the databases (default:
#                         /usr/local/envs/antismash/lib/python3.7/site-
#                         packages/antismash/databases).
#   --write-config-file PATH
#                         Write a config file to the supplied path
#   --without-fimo        Run without FIMO (lowers accuracy of RiPP precursor
#                         predictions)
#   --executable-paths EXECUTABLE:PATH,EXECUTABLE2:PATH2,...
#                         A comma separated list of executable name->path pairs
#                         to override any on the system path.E.g.
#                         diamond=/alternate/path/to/diamond,hmmpfam2=hmm2pfam
#   --hmmdetection-strictness {strict,relaxed,loose}
#                         Defines which level of strictness to use for HMM-based
#                         cluster detection, (default: relaxed).

# Debugging & Logging options:

#   -v, --verbose         Print verbose status information to stderr.
#   -d, --debug           Print debugging information to stderr.
#   --logfile PATH        Also write logging output to a file.
#   --list-plugins        List all available sec. met. detection modules.
#   --check-prereqs       Just check if all prerequisites are met.
#   --limit-to-record RECORD_ID
#                         Limit analysis to the record with ID record_id
#   -V, --version         Display the version number and exit.
#   --profiling           Generate a profiling report, disables multiprocess
#                         python.
#   --skip-sanitisation   Skip input record sanitisation. Use with care.
#   --skip-zip-file       Do not create a zip of the output

# Debugging options for cluster-specific analyses:

#   --minimal             Only run core detection modules, no analysis modules
#                         unless explicitly enabled
#   --enable-genefunctions
#                         Enable Gene function annotations (default: enabled,
#                         unless --minimal is specified)
#   --enable-tta          Enable TTA detection (default: enabled, unless
#                         --minimal is specified)
#   --enable-lanthipeptides
#                         Enable Lanthipeptides (default: enabled, unless
#                         --minimal is specified)
#   --enable-thiopeptides
#                         Enable Thiopeptides (default: enabled, unless
#                         --minimal is specified)
#   --enable-nrps-pks     Enable NRPS/PKS analysis (default: enabled, unless
#                         --minimal is specified)
#   --enable-sactipeptides
#                         Enable sactipeptide detection (default: enabled,
#                         unless --minimal is specified)
#   --enable-lassopeptides
#                         Enable lassopeptide precursor prediction (default:
#                         enabled, unless --minimal is specified)
#   --enable-t2pks        Enable type II PKS analysis (default: enabled, unless
#                         --minimal is specified)
#   --enable-html         Enable HTML output (default: enabled, unless --minimal
#                         is specified)

# Gene finding options (ignored when ORFs are annotated):

#   --genefinding-tool {glimmerhmm,prodigal,prodigal-m,none,error}
#                         Specify algorithm used for gene finding: GlimmerHMM,
#                         Prodigal, Prodigal Metagenomic/Anonymous mode, or
#                         none. The 'error' option will raise an error if
#                         genefinding is attempted. The 'none' option will not
#                         run genefinding. (default: error).
#   --genefinding-gff3 GFF3_FILE
#                         Specify GFF3 file to extract features from.

# Full HMMer options:

#   --fullhmmer-pfamdb-version FULLHMMER_PFAMDB_VERSION
#                         PFAM database version number (e.g. 27.0) (default:
#                         latest).

# ClusterFinder options:

#   --cf-min-cds INT      Minimum size of a ClusterFinder cluster, in number of
#                         CDS features (default: 5).
#   --cf-mean-threshold FLOAT
#                         Minimum mean probability threshold (default: 0.6).
#   --cf-min-pfams INT    Minimum number of biosynthetic PFam domains in a
#                         ClusterFinder cluster (default: 5).

# Cluster HMMer options:

#   --clusterhmmer-pfamdb-version CLUSTERHMMER_PFAMDB_VERSION
#                         PFAM database version number (e.g. 27.0) (default:
#                         latest).

# NRPS/PKS options:

# ClusterBlast options:

#   --cb-nclusters count  Number of clusters from ClusterBlast to display,
#                         cannot be greater than 50. (default: 10)
#   --cb-min-homology-scale LIMIT
#                         A minimum scaling factor for the query BGC in
#                         ClusterBlast results. Valid range: 0.0 - 1.0. Warning:
#                         some homologous genes may no longer be visible!
#                         (default: 0.0)
