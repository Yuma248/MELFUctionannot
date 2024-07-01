#!/usr/bin/perl -w
while (@ARGV){
$_=shift @ARGV;
if ($_=~ /^-i$/){$if=shift @ARGV;}
elsif ($_=~ /^-o$/){$out=shift @ARGV;}
elsif ($_=~ /^-t$/){$type=shift @ARGV;}
elsif ($_=~ /^-nc$/){$nc=shift @ARGV;}
}
if (not defined ($if && $out)){print "\nThis script creates a nexus file for several gene aligments to run IQTREE, each gene should be aligned and in its own phy format file.\n\nUsage:\nnexusIQ3.pl\n\t-i <path to folder with all the phy files>\n\t-o <path to nexus ouput>\nOPTIONAL\n\t-t <type DNA, protein, or codon, this is relevant just for CODON (DNA on exons), but IQtree recognize the other data types>\n\t-nc <not yet implemented but number of cores to run in parallel>\n\nFor example:\n\n.nexusIQ3.pl -i ./best15/ -o /mynexus.nex -t protein -nc 10 \n\n"; exit;}
if (not defined ($nc)) {$nc= 10;}
#if (not defined ($type)) {$type= "protein";}
##unless (-d $out){`mkdir $out`;}
my @files =`ls $if\/\*\.phy`;
`echo "\#nexus\nbegin sets\;" >> $out`;
my $nf = 0;
foreach $f (@files){
        $nf += 1;
        chomp $f;
        my $infl = `head -n 1 $f`;
        chomp $infl;
        my @inf =split (' ', $infl);
        if (not defined ($type)){
                our $line ="\tcharset part$nf = $f\: 1-$inf[1]\;";
        }elsif ($type == "CODON"){
                our $line ="\tcharset part$nf = $f\:CODON, 1-$inf[1]\;";
        }
        `echo "$line" >> $out`;
}
`echo "end\;" >> $out`;
