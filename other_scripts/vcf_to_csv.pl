#!/usr/bin/perl
use strict;
use warnings;
use Array::Transpose;

#first use vcf-merge to merge parents and progeny
#vcf-merge output_raw_maf0.1_00.9_AN30_positionsforparents.vcf.gz imputed_AN60_gtgl_4.1.vcf.gz > parents_prog_snps.vcf

#Get input snp file
if (!$ARGV[0]){
    print "Provide an input multi vcf file.";
    exit 1;
}

my $infile = $ARGV[0];

#Open input file
open (INFILE, "<$infile")  || die "Cannot open input VCF file!";

#split into individual chromosome files
my $last_chromo = '';
my @header;

while (<INFILE>){
    chomp;

    if ($_ =~ /^##/){
	next;
    }

    elsif ($_ =~ /^#CHROM/){
	@header = split (/\t/, $_);
	splice @header, 0, 11;
    }

    elsif ($_ =~ /^SL2.50.*/){
	my @row = split (/\t/, $_);
	splice @row, 2, 7;
	my $chromo = $row[0];
	my $pos = $row[1];
	my $marker_name = $chromo . "_" . $pos;
	splice @row, 0, 2;
	my @parentb = split (/:/, $row[0]);
	my @parenta = split (/:/, $row[1]);
	splice @row, 0, 2;

	#disregard SNPs with missing data or heterozygous in one of the parents.
	if ($parenta[0] eq "0/1" || $parenta[0] eq "./." || $parentb[0] eq "0/1" || $parentb[0] eq "./." || $parenta[0] eq $parentb[0]){ next;}

	if ($chromo eq $last_chromo){
	    print OUTFILE $marker_name . "\t" . $pos;
	
	    foreach my $row (@row){
		my @geno = split (/:/, $row);

		if ($geno[0] eq $parenta[0]){
		    print OUTFILE "\ta";
		}

		elsif ($geno[0] eq "0/1"){
		    print OUTFILE "\th";
		}

		elsif ($geno[0] eq $parentb[0]){
		    print OUTFILE "\tb";
		}

		else {print OUTFILE "\t-";
		}
	    }
	    print OUTFILE "\n";
	}

	else {
	    open (OUTFILE, ">split$chromo.csv") || die "Cannot open output file!";
	    print OUTFILE "marker_name\tposition\t";
	    print OUTFILE join("\t", @header);
	    print OUTFILE "\n";
	    print OUTFILE $marker_name . "\t" . $pos;

            foreach my $row (@row){
                my @geno = split (/:/, $row);

                if ($geno[0] eq $parenta[0]){
                    print OUTFILE "\ta";
                }

                elsif ($geno[0] eq "0/1"){
                    print OUTFILE "\th";
                }

                elsif ($geno[0] eq $parentb[0]){
                    print OUTFILE "\tb";
                }

                else {print OUTFILE "\t-";
                }
	    }
	    print OUTFILE "\n";
	    $last_chromo = $chromo;
	}
    }
}
