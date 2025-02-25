First step for first time users using command line:
For Windows Users: Requires Ubuntu in order to run properly once installed to install ubuntu you can follow this video here it also works for windows 11 and windows 10 https://www.youtube.com/watch?v=X-DHaQLrBi8&t=2s

For Linux and Mac Users: Will work normally running on the command line 

Second step:
Navigate to where you want to install One HDR Ready and install python which can be downloaded here https://www.python.org/downloads/
Required Python libraries for the tool: requests, pandas and biopython which can be installed with the command below
```
pip install requests pandas biopython
```
The tool uses crispor to generate sgRNA sequences and can be installed here https://github.com/maximilianh/crisporWebsite/blob/master/INSTALL.md

Once crispor is installed the program can be run and enter the gene of interest
```
python3 CRISPR-Webtool.py
```
