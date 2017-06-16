#Other Scripts for working with SNPbinner

##vcf\_to\_csv.pl
This script take a VCF file and splits it into per-chromosome CSV files that can be run through SNPbinner. The input VCF files should include both parents and progeny, with the parents appearing first in the VCF file. As an example, this could be performed using `vcf-merge` like so:
```vcf-merge positionsforparents.vcf.gz imputed.vcf.gz > parents_prog_snps.vcf```