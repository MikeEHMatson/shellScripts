#######################################################
########## Busco duplicate regions confirmation #######
#######################################################

# Names and locations of four bamfiles of interest to me
p618=/bamfiles/618-final_PacBio1.2.bam
p1306=/bamfiles/1306_MK_miSeq.pacBioV2.bam
p550=/bamfiles/550-final_pacBio1.2.bam
p6629=/bamfiles/6629-final_pacBio1.2.bam

# Location of reference genome fasta file
REF=reference_version1.fa

# Directory to write individual chromosome/contig fasta files to
chrDir=/blast/chrs

# Directory containing ncbi blast database files
dbDir=/blast/database
db=pacbio1.2_blastdb

# load ncbi blast module and samtools
module load ncbi-blast/2.6.0
module load samtools

# make blast database from reference sequence
makeblastdb -in $REF -out $dbDir/reference_v1 -dbtype nucl -parse_seqids

#### Output all contigs/chromosomes of reference genome as individual fasta files
# This file was generated from the program ALLMAPS and contains each contig in the reference in column 1 ($0) and its size in column 2 ($1)
awk '{print $1}' reverence_vsersion1.sizes > allChrList.txt
list='allChrList.txt'
allCtg=`cat $list`

for contig in $allCtg; 
do
samtools faidx $REF $contig > $chrDir/$contig.fa; 
done

### Example of a duplicated busco identified in my reference genome from a previously run busco analysis
##### EPrGT00050000001214 #####
##### this busco found on two contigs in reference_version1.fa, ctg32 and ctg64
##### Blast the smaller contig (ctg64) against entire genome and look for long alignment regions
contig=ctg64
## renamed here to match the reference genome that generated the bamfile. Ideally, this is unnecessary if the bamfile was generated drectly from reference_version1.fa
chr=chr000064F
# Calculate length of entire contig
wc -c < $chrDir/$contig.fa
blastn -query $chrDir/$contig.fa -out blastOut.$contig.txt -db $dbDir/reference_v1 -evalue 9e-150 -outfmt 7 -num_threads 8

#### Blast generates lots of hits, but one region was identified manually with high % identity and long alignment length.
#### Region happens to be on ctg32, which is expected
region=732980-754928
## Calculate length of alignment
declare -i regionLength
regionLength=($region)/-1

## Do the same with the longer contig and its region hit
contigPrime=ctg32
regionPrime=1540-23486
chrPrime=chr000032F

## Parese out blast hit region
samtools faidx $REF_1_2 $contig:$region > $contig.$region.fa
blastn -query $contig.$region.fa -out blastOut.$contig.$region.txt -db $dbDir/pacbio1.2_blastdb -evalue 9e-150 -outfmt 7 -num_threads 8;

## Calculate depth in regions for all four parents or sample bamfiles of interest.
## Generates single Test file and adds data to the end of each section.
parent=618
echo "Depth of $regionLength bp region: 
$contig:$region for parent $parent
From a total contig length of:
$(wc -c < $chrDir/$contig.fa)" > $contig:$region.depth.txt
samtools depth $p618 -r $chr:$region |  awk '{sum+=$3} END { print "Average = ",sum/NR}' >> $contig:$region.depth.txt
echo "

Depth of region: 
$contigPrime:$regionPrime for parent $parent
From a total contig length of:
$(wc -c < $chrDir/$contigPrime.fa)" >> $contig:$region.depth.txt
samtools depth $p618 -r $chrPrime:$regionPrime |  awk '{sum+=$3} END { print "Average = ",sum/NR}' >> $contig:$region.depth.txt
echo "

Depth of uniquely mapped reads
for $regionLength bp region 
$contig:$region for parent $parent" >> $contig:$region.depth.txt
samtools depth $p618 -r $chr:$region -Q 10|  awk '{sum+=$3} END { print "Average = ",sum/NR}' >> $contig:$region.depth.txt
echo "

Depth of entire 
$contig for $parent" >> $contig:$region.depth.txt
samtools depth $p618 -r $chr |  awk '{sum+=$3} END { print "Average = ",sum/NR}' >> $contig:$region.depth.txt
echo "

Depth of entire 
$contigPrime for $parent" >> $contig:$region.depth.txt
samtools depth $p618 -r $chrPrime |  awk '{sum+=$3} END { print "Average = ",sum/NR}' >> $contig:$region.depth.txt

parent=1306
echo "

Depth of $regionLength bp region: 
$contig:$region for parent $parent
From a total contig length of:
$(wc -c < $chrDir/$contig.fa)" >> $contig:$region.depth.txt
samtools depth $p1306 -r $chr:$region |  awk '{sum+=$3} END { print "Average = ",sum/NR}' >> $contig:$region.depth.txt
echo "

Depth of region: 
$contigPrime:$regionPrime for parent $parent
From a total contig length of:
$(wc -c < $chrDir/$contigPrime.fa)" >> $contig:$region.depth.txt
samtools depth $p1306 -r $chrPrime:$regionPrime |  awk '{sum+=$3} END { print "Average = ",sum/NR}' >> $contig:$region.depth.txt
echo "

Depth of uniquely mapped reads
for $regionLength bp region 
$contig:$region for parent $parent" >> $contig:$region.depth.txt
samtools depth $p1306 -r $chr:$region -Q 10|  awk '{sum+=$3} END { print "Average = ",sum/NR}' >> $contig:$region.depth.txt
echo "

Depth of entire 
$contig for $parent" >> $contig:$region.depth.txt
samtools depth $p1306 -r $chr |  awk '{sum+=$3} END { print "Average = ",sum/NR}' >> $contig:$region.depth.txt
echo "

Depth of entire 
$contigPrime for $parent" >> $contig:$region.depth.txt
samtools depth $p1306 -r $chrPrime |  awk '{sum+=$3} END { print "Average = ",sum/NR}' >> $contig:$region.depth.txt

parent=550
echo "

Depth of $regionLength bp region: 
$contig:$region for parent $parent
From a total contig length of:
$(wc -c < $chrDir/$contig.fa)" >> $contig:$region.depth.txt
samtools depth $p550 -r $chr:$region |  awk '{sum+=$3} END { print "Average = ",sum/NR}' >> $contig:$region.depth.txt
echo "

Depth of region: 
$contigPrime:$regionPrime for parent $parent
From a total contig length of:
$(wc -c < $chrDir/$contigPrime.fa)" >> $contig:$region.depth.txt
samtools depth $p550 -r $chrPrime:$regionPrime |  awk '{sum+=$3} END { print "Average = ",sum/NR}' >> $contig:$region.depth.txt
echo "

Depth of uniquely mapped reads
for $regionLength bp region 
$contig:$region for parent $parent" >> $contig:$region.depth.txt
samtools depth $p550 -r $chr:$region -Q 10|  awk '{sum+=$3} END { print "Average = ",sum/NR}' >> $contig:$region.depth.txt
echo "

Depth of entire 
$contig for $parent" >> $contig:$region.depth.txt
samtools depth $p550 -r $chr |  awk '{sum+=$3} END { print "Average = ",sum/NR}' >> $contig:$region.depth.txt
echo "

Depth of entire 
$contigPrime for $parent" >> $contig:$region.depth.txt
samtools depth $p550 -r $chrPrime |  awk '{sum+=$3} END { print "Average = ",sum/NR}' >> $contig:$region.depth.txt

parent=6629
echo "

Depth of $regionLength bp region: 
$contig:$region for parent $parent
From a total contig length of:
$(wc -c < $chrDir/$contig.fa)" >> $contig:$region.depth.txt
samtools depth $p6629 -r $chr:$region |  awk '{sum+=$3} END { print "Average = ",sum/NR}' >> $contig:$region.depth.txt
echo "

Depth of region: 
$contigPrime:$regionPrime for parent $parent
From a total contig length of:
$(wc -c < $chrDir/$contigPrime.fa)" >> $contig:$region.depth.txt
samtools depth $p6629 -r $chrPrime:$regionPrime |  awk '{sum+=$3} END { print "Average = ",sum/NR}' >> $contig:$region.depth.txt
echo "

Depth of uniquely mapped reads
for $regionLength bp region 
$contig:$region for parent $parent" >> $contig:$region.depth.txt
samtools depth $p6629 -r $chr:$region -Q 10|  awk '{sum+=$3} END { print "Average = ",sum/NR}' >> $contig:$region.depth.txt
echo "

Depth of entire 
$contig for $parent" >> $contig:$region.depth.txt
samtools depth $p6629 -r $chr |  awk '{sum+=$3} END { print "Average = ",sum/NR}' >> $contig:$region.depth.txt
echo "

Depth of entire 
$contigPrime for $parent" >> $contig:$region.depth.txt
samtools depth $p6629 -r $chrPrime |  awk '{sum+=$3} END { print "Average = ",sum/NR}' >> $contig:$region.depth.txt

