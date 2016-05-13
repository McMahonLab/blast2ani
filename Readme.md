## BLAST 2 ANI

This repo documents the workflow for the following manuscript:

Garcia, S. L., Stevens, S. L. R., Crary, B., Martinez-Garcia, M., Stepanauskas, R., Woyke, T., Tringe, S. G., Andersson,
S., Bertilsson, Malmstrom, R. R., McMahon, K. D. (prepared for submission). Genome level exploration of abundant and
uncultivated freshwater bacteria reveals distinct and dynamic populations.


### One Metagagenome

This section documents the workflow for blasting the reads from one metagenome, with no identity cutoff, against a set of single cell genomes.

Ran the following on each SAG individually.  Where:
	- blastn version = 2.2.31
	- SAG = SAG fna file (without rrna and contigs renamed for parsing ease) 
	- metagenome = the fna file of the metagenome reads
	- outname = name for resulting blastfile 

```
makeblastdb -in SAG -out SAG.db -dbtype nucl
blastn -task blastn -db SAG.db -query metagenome -out outname.blast -evalue 0.001 -outfmt 6
```

For each resulting blast file, filtered out hits shorter than 200bp (alignment length is the 4th column in the standard BLAST out format 6).
```
awk '($4 > 200)' outname.blast > outname.blast.len200
```

Wanted only the best hits (if a read hits that SAG/genome more than once).
Wrote blast_besthit.py to keep only the top hit (based on bit score).  If same bit score, keeps first one.
Ran the following command on each filtered blast result.

```
./blast_besthit.py -blast_in outname.blast.len200
```
