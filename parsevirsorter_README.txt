To create the genomes list (of genomes done) use the following command:
find "data" -printf '%P\n' | awk -F 'VIRSorter' '$1{print substr($1,1,length($1)-1)}' | sort -u > genomes_list

To run virsorter you can use the bin/run_virsorter.bash script. This can also run in parallel.

There are dependencies for this.  They are hard coded to look for these in /home/ksenia/ (not me).
(I always had a hard time setting everything up. There are some things you have to alter in the code itself. IIRC there are also some hard coded path to a certain Franklin folder as that was the folder I by chance had the program in.)
It might be better to get a new version from the github. It seemed still maintained.
https://github.com/simroux/VirSorter

To move the correct files (and rename them properly) required for parsing, you can use the bin/move_virsorter.bash scripts.

For parsing you can make use of the bin/parsevirsorter.py script. There is a --help flag available. Output is to stdin and stderr.

$ find data/directory -type f -name '*.csv' | xargs python3 bin/parsevirsorter.py > phages.coords.tsv 2> parsevirsorter.log

Or parallelized

$ find data/directory -type f -name '*.csv' | parallel --xargs -P $NUMCORES python3 bin/parsevirsorter.py

The output is tab seperated with
genome_id, contig_id, phage_id, start_base, end_base, gene_start, gene_end

With phage_id as follows
1-3 phage with resp. confident, possible, doubtable
4-6 prophage with resp. confident, possible, doubtable

(This comes straight out the paper, so it is good if you read it)
https://peerj.com/articles/985/

The genes are numbered as found by the predict files.

If there are no numbers (NA), the whole contig is meant.

----

Virsorter does multiple analysis (also see paper) and these possibly indicate prophage regions. The parser takes all ranges of all these analysis.

The csv file only expresses these ranges as gene numbers. So the program has to match these gene numbers with the base coords in the predict file.

----

There are some errors in the output files. See the log for this. There are no fatal errors, but some files are broken, giving small loose of information. It looks like a flush problem, so maybe rerunning these can solve this.

There are also a few completely empty files, because of the program crashing. I don't know what boundary specific boundary condition this causes. Now is assumed there are no phages, but I am certain about that.

----

Contact:
Hielke Walinga
(hielkewalinga@gmail.com)