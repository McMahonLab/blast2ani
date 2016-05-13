## BLAST 2 ANI

This repo documents the workflow for the following manuscript:

Garcia, S. L., Stevens, S. L. R., Crary, B., Martinez-Garcia, M., Stepanauskas, R., Woyke, T., Tringe, S. G., Andersson,
S., Bertilsson, Malmstrom, R. R., McMahon, K. D. (prepared for submission). Genome level exploration of abundant and
uncultivated freshwater bacteria reveals distinct and dynamic populations.


### BLASTing metagenome closest to SAG collection against SAGs

This section documents the workflow for BLASTing the reads from one metagenome (collected closest to SAG collection), with no identity cutoff, against a set of single cell genomes.  Resulting in Figure 2 of the manuscript.

#### Running BLAST and filtering results

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
./blast\_besthit.py -blast_in outname.blast.len200
```
Resulting files have same name and end in .bbh

#### Non-competitive results

Parsed results into one file for each group (ME-acI, nonME-acI, ME-LD12, nonME-LD12, and ME-others)

Visualized results with vis_scripts/makeFig2\_plotseqdiscden.R

#### Competitive results

Also needed results that were competitive within acI and within LD12 groups.  Used blast_besthit.py to choose.

Concatenated all the within SAG bbh results for LD12 (PTXW.len150-vs-LD12\_norrna\_short.blast.len200.bbh) and acI (PTXW.len150-vs-acI\_norrna.fna\_short.blast.len200.bbh) into one file for each group.

Ran the following commands on them, -kb stand for keep all best, so it will count a read twice if the bit score is the same for 2 hits (since we can't tell which one it belongs to)
```
./blast\_besthit.py -bin LD12/PTXW.len150-vs-LD12\_norrna\_short.blast.len200.bbh -kb
./blast\_besthit.py -bin acI/PTXW.len150-vs-acI\_norrna.fna_short.blast.len200.bbh -kb
```

Parsed results into one file for each group (ME-acI, nonME-acI, ME-LD12, nonME-LD12)

Visualized results with vis\_scripts/competative\_plotseqdiscden.R

### BLASTing all metagenomes against all SAGs

#### Running BLAST 
Ran the following on each SAG and metagenome combo individually.  Where:
	- blastn version = 2.2.31
	- used already made blastdb's from blast single metagenome (info above)
	- SAG.db = SAG blast db made from fna file (without rrna and contigs renamed for parsing ease) 
	- metagenome = the fna file of the metagenome reads
	- outname = name for resulting blastfile 

```
blastn -task blastn -db SAG.db -query metagenome -out outname.blast -evalue 0.001 -outfmt 6 -perc_identity 95
```
I named all of the output files with this scheme: metagenome-vs-SAG.blast

#### Reformatting, filtering, and pooling resulting BLAST files

Reformatted the BLAST results to be more useful for me. Then filtered them only keeping hits longer than 200bp and > 97.5 identity.  Then I replaced the metagenome names with the month and year to pool the results by month.

##### Reformatting resulting BLAST files

All the raw blast results were in a folder called raw\_blast\_all
Ran the following commands to get the names of all the SAGs from the blast results and then use grep to combine and add the metagenome information in at the same time.
```
ls raw_blast_all/*.blast | cut -d - -f 3 | cut -d _ -f 1 | sort | uniq > SAGnames.txt
while read line; do echo $line; grep $line raw_blast_all/*$line*.blast > MEmeta-vs-$line.blast; done < SAGnames.txt
```

Removing path from metagenome name and then making that its own column.  All my metagenome names have '.len150' after the library name so I split on that.
```
for file in MEmeta-vs-A*; do sed -i.ibak 's/raw\_blast\_all\///g' $file; done
rm *.ibak
for file in MEmeta-vs-A*.blast; do sed -i.ibak 's/.len150/        len150/g' $file; done
rm *.ibak
```
Removed backups after checking that proccess was successful. 

##### Filtering resulting blast files

Ran the following for each reformatted blast result file.  Where:
	- blastfile = newly reformatted blast file

```
awk '($5 > 200)' blastfile | awk '($4 > 97.5)' > blastfile.len200.id975; done
```

##### Pooling Results by Month

Pooled data by month by replacing metagenome names with year and month by running poolBLASTS.py on each filtered and formated blast result file. Where:
	- blastfile.len200.id975 = the resulting file from reformatting and filtering.
	- sample_data.txt = a tab-separated file with the following columns in it (in order):
		sample, reads, bps, layer, date, year, month, day
		
```
./poolBLASTS.py blastfile.len200.id975 sample_data.txt
```


