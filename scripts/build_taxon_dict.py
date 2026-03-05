#!/usr/bin/env python3

import argparse, logging, tarfile, subprocess, os, tempfile, pickle

parser = argparse.ArgumentParser(description="Downloads NCBI taxdump and saves as a pickle dictionary")
parser.add_argument('-o', '--outfile', help="Output pickle file", type=str, required=True)
parser.add_argument('-log', '--logfile', help="Logfile", type=str, required=False)
args = parser.parse_args()

if args.logfile:
    logging.basicConfig(filename=args.logfile, level=logging.INFO, 
                        format='%(levelname)s %(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
else:
    logging.basicConfig(level=logging.INFO, format='%(levelname)s %(asctime)s %(message)s')

logging.info("Downloading new_taxdump..")
fd, temp = tempfile.mkstemp()
os.close(fd) 

try:
    subprocess.run(
        ['wget', '-O', temp, 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz'], 
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True
    )
except (subprocess.CalledProcessError, FileNotFoundError):
    logging.info("wget failed or missing, trying curl...")
    subprocess.run(
        ['curl', '-o', temp, 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz'], 
        check=True
    )

logging.info("Creating taxon dictionary..")
taxon_dict = {}

with tarfile.open(temp, "r:gz") as tar:
    for line in tar.extractfile("rankedlineage.dmp"):
        line = line.decode().strip().replace('\t','').split('|')
        taxon_dict[line[0]] = line[1:-1]

with tarfile.open(temp, "r:gz") as tar:
    for line in tar.extractfile("merged.dmp"):
        line = line.decode().strip().replace('\t','').split('|')
        taxon_dict[line[0]] = taxon_dict[line[1]]

with tarfile.open(temp, "r:gz") as tar:
    for line in tar.extractfile("delnodes.dmp"):
        line = line.decode().strip().replace('\t','').split('|')
        taxon_dict[line[0]] = taxon_dict['1']

os.unlink(temp)

logging.info(f"Saving dictionary to {args.outfile}...")
# Save the dictionary as a binary file
with open(args.outfile, 'wb') as f:
    pickle.dump(taxon_dict, f, protocol=pickle.HIGHEST_PROTOCOL)

logging.info("Done!")
