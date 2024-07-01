#!/usr/bin/perl -w
while (@ARGV){
$_=shift @ARGV;
if ($_=~ /^-if$/){$if=shift @ARGV;}
elsif ($_=~ /^-of$/){$of=shift @ARGV;}
elsif ($_=~ /^-gnc$/){$gnc=shift @ARGV;}
}
if (not defined ($if)){print "\nThis script calculates likelihood ratio tests for models of different ω ratios for branches. It uses paml results of model M0, M2a, M2a1 and M-2, the resutls should be in a folder with names format as <Gene_lineages_model.paml>.\n\nUsage:\nlrt.pl\n\t-if <path to a folder where paml resutls are>\n\t-of <path to output folder, deafult LRToutput>\n\t-gnc <number of genes to run in parallel, deaful 10. The gnc*lnc should not exceeds the toral numebr of core in system>\n\nFor example:\n\nlrt.pl -if /Tursiops/pamlretuslt/ -of /Tursiops/lrt/ -gnc 10\n\n"; exit;}
if (not defined ($of)) {$of="./LRToutput/";}
unless (-d $of){`mkdir $of`;}
if (not defined ($lnc)) {$lnc=1;}
if (not defined ($gnc)) {$gnc=10;}
our %lrt=(); 
our @genes=`ls  $if/ | grep "M2a." | cut -d "_" -f 1 | sort | uniq`;
our @lng =`ls  $if/ | grep "M2a." | cut -d "_" -f 2 | sort | uniq`;
my $if= quotemeta($if);
my $of= quotemeta($of);
chomp @genes;
chomp @lng;
use Parallel::Loops;
my $rawout=$of."/valuesT.txt";
`echo "Lineage\tGene\tW0\tWF\tWB\tLRTN\tSig\tLRT1\tSig" >> $rawout`;
my $plgen = Parallel::Loops->new($gnc);
$plgen->share( \%lrt);
$plgen->foreach (\@genes, sub{
my $gene = $_ ;
#our %glrt = ();
chomp $gene;
$gene=~ s/$if//;
$gene=~ s/\///;
foreach $lineage (@lng){
chomp $lineage;
$lineage=~ s/$if//;
$lineage=~ s/\///;

$lrt{$lineage}{$gene}{sig0}="NS";
$lrt{$lineage}{$gene}{sig1}="NS";
$lrt{$lineage}{$gene}{sig2n}="NS";
my $codemlr = $if."/".$gene."_".$lineage;
my $filename= $codemlr."_M0.paml";
if (-e $filename) {
$lrt{$lineage}{$gene}{M0}=`grep -m 1 "lnL" $codemlr\_M0.paml | awk '{print \$5}'`;
$lrt{$lineage}{$gene}{M2a1}=`grep "lnL" $codemlr\_M2a1.paml | awk '{print \$5}'`;
$lrt{$lineage}{$gene}{M2a}=`grep "lnL" $codemlr\_M2a.paml | awk '{print \$5}'`;
$lrt{$lineage}{$gene}{M2n}=`grep "lnL" $codemlr\_M2.paml | awk '{print \$5}'`;
$lrt{$lineage}{$gene}{M2s0}= -2 * log ( $lrt{$lineage}{$gene}{M2a} / $lrt{$lineage}{$gene}{M0});
$lrt{$lineage}{$gene}{M2s1}= -2 * log ( $lrt{$lineage}{$gene}{M2a} / $lrt{$lineage}{$gene}{M2a1});
$lrt{$lineage}{$gene}{M2s2n}= -2 * log ( $lrt{$lineage}{$gene}{M2a} / $lrt{$lineage}{$gene}{M2n});
$lrt{$lineage}{$gene}{W0}=`grep "omega (dN/dS) =" $codemlr\_M0.paml | awk '{print \$4}' | perl -pe 's/\n//'`;
$lrt{$lineage}{$gene}{WF}=`grep -A 5 "MLEs of dN/dS" $codemlr\_M2a.paml | tail -n 1 | awk '{print \$5}' | perl -pe 's/\n//'`;
$lrt{$lineage}{$gene}{WB}=`grep -A 5 "MLEs of dN/dS" $codemlr\_M2a1.paml | tail -n 2 | head -n 1 | awk '{print \$5}' | perl -pe 's/\n//'`;
if ($lrt{$lineage}{$gene}{M2s0} > 3.84){$lrt{$lineage}{$gene}{sig0}="\*";}
if ($lrt{$lineage}{$gene}{M2s1} > 3.84){$lrt{$lineage}{$gene}{sig1}="\*";}
if ($lrt{$lineage}{$gene}{M2s2n} > 6.63){$lrt{$lineage}{$gene}{sig2n}="\*";}
 
`echo "$lineage\t$gene\t$lrt{$lineage}{$gene}{W0}\t$lrt{$lineage}{$gene}{WF}\t$lrt{$lineage}{$gene}{WB}\t$lrt{$lineage}{$gene}{M2s0}\t$sig0\t$lrt{$lineage}{$gene}{M2s1}\t$sig1\t$lrt{$lineage}{$gene}{M2s2n}\t$sig2n" >> $rawout`;
}
}
});


