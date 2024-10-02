For Windows Users: Requires Ubuntu in order to run properly once installed For Linux and Mac Users: Should work normally running on the command line Required Python libraries for the tool: requests, pandas and biopython Make sure to have java installed as well.
```
pip install requests pandas biopython
```
First step once the repository is cloned into your enviroment is to index the human genome which can be downloaded here. https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/ Make sure to download the file Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa.gz specificially.

For indexing the genome, do this command in the command line this will be the longest process as it is a search for off target sites in the entire genome potentially taking 3 hours.
```
mkdir tmp
java -Xmx4g -jar FlashFry-assembly-1.15.jar \
 index \
 --tmpLocation ./tmp \
 --database GRCh38_cas9ngg_database \
 --reference Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa.gz \
 --enzyme spcas9ngg
```
Once the database is created simply run the code and enter your gene of interest and it will run properly.
```
python3 CRISPR-Webtool.py
```
