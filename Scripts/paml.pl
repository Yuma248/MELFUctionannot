
#!/usr/bin/perl -w
while (@ARGV){
$_=shift @ARGV;
if ($_=~ /^-if$/){$if=shift @ARGV;}
elsif ($_=~ /^-of$/){$of=shift @ARGV;}
elsif ($_=~ /^-tf$/){$tf=shift @ARGV;}
elsif ($_=~ /^-ll$/){$ll=shift @ARGV;}
elsif ($_=~ /^-lnc$/){$lnc=shift @ARGV;}
elsif ($_=~ /^-gnc$/){$gnc=shift @ARGV;}
elsif ($_=~ /^-mod$/){$mod=shift @ARGV;}
}
if (not defined ($if && $tf)){print "\nThis script use paml to calculate omega (ratio synonymous and nonsynonymous substitution) under both M0 (same omega all lineageas) and M2a (two omega ratios) models. It uses aligments of single copy completed genes from several linegeas.\n\nUsage:\npamlYuma.pl\n\t-if <path to a folder where aligmnts of all genes are>\n\t-ll <a file with the list of lineages name as in the aligments>\n\t-of <path to output folder, deafult output>\n\t-tf <path to folder with guide trees, label for each lineages to be test as \#1, see paml manual>\n\t-lnc <number lineageas to run in parallel, deaful 5>\n\t-gnc <number of genes to run in parallel, deaful 10. The gnc*lnc should not exceeds the toral numebr of core in system>\n\t-mod <what model of paml you want to run M2a(M2a) M2b(M2a1) M0(M0) or M1(M1),dafault all>\n\nFor example:\n\npamlYuma.pl -if /Tursiops/prankaligments/ -of /Tursiops/omega/ -ll /Tursiops/lineages.txt -tf /Tursiops/trees/ -lnc 5 -gnc 10 -mod M1,M2a1,M0,M2a\n\n"; exit;}
$of //= "./output/";
unless (-d $of){`mkdir $of`;}
my $ctlf=$of."/ctrdir/";
unless (-d $ctlf){`mkdir $ctlf`;}
$lnc //= 5;
$gnc //= 10;
$mod //= "M2a,M2a1,M0,M1";
our @MODS =split (/,/, $mod);
our %pmods = map { $_ => 1 } @MODS;
@genes=`ls $if\/*.phy`;
print "$genes[0]\n";
if ( $genes[0] =~ /best/){
	our $inext = ".best.phy";
	s/\.best.phy// for @genes;
}
elsif ( $genes[0] =~ /_dna.phy/){
	our $inext = "_dna.phy";
	s/\_dna.phy// for @genes;
}
elsif ( $genes[0] =~ /_dna_gb.phy/){
	our $inext = "_dna_gb.phy";
	s/\_dna_gb.phy// for @genes;
}
#my $if=quotemeta($if);
#my $of=quotemeta($of);
#my $tf=quotemeta($tf);
use Parallel::Loops;
chomp @genes;
my $plgen = Parallel::Loops->new($gnc);
$plgen->foreach (\@genes, sub{
my $gene = $_ ;
chomp $gene;
$gene=~ s/$if//;
$gene=~ s/\///;
my @rlng=`cat $ll`;
my @lng =grep ! /""/, @rlng;
print "$lng[0]\n";
my $pllin = Parallel::Loops->new($lnc);
$pllin->foreach (\@lng, sub{
my $lineage = $_ ;
chomp $lineage;
print "We are running the \"$lineage\" lineage\n"; 
next unless ($lineage ne "");
my $contfilM0=$ctlf."\/".$gene."_".$lineage."M0.ctl";
my $contfilM2a=$ctlf."\/".$gene."_".$lineage."M2a.ctl";
my $contfilM2a1=$ctlf."\/".$gene."_".$lineage."M2a1.ctl";
my $contfilM2=$ctlf."\/".$gene."_".$lineage."M-2.ctl";
my $contfilM1c=$ctlf."\/".$gene."_".$lineage."M1c.ctl";
my $contfilM2c=$ctlf."\/".$gene."_".$lineage."M2c.ctl";
my $input="$if\/$gene$inext";
my $outputM0="$of\/$gene\_$lineage\_M0.paml";
my $outputM2a="$of\/$gene\_$lineage\_M2a.paml";
my $outputM2a1="$of\/$gene\_$lineage\_M2a1.paml";
my $outputM2="$of\/$gene\_$lineage\_M2.paml";
my $outputM2c="$of\/$gene\_$lineage\_M2c.paml";
my $outputM1c="$of\/$gene\_$lineage\_M1c.paml";

my $treef= "$tf\/$lineage\.nex";
if (exists($pmods{"M2a1"})){
cnt_file($input, $outputM2a1, $treef, "M2a1", $contfilM2a1); 
`codeml $contfilM2a1`;
}
if (exists($pmods{"M2a"})){
cnt_file($input, $outputM2a, $treef, "M2a", $contfilM2a);
`codeml $contfilM2a`;
}
if (exists($pmods{"M0"})){
cnt_file($input, $outputM0, $treef, "M0", $contfilM0);
`codeml $contfilM0`;
}
if (exists($pmods{"M1"})){
cnt_file($input, $outputM2, $treef, "M1", $contfilM2);
`codeml $contfilM2`;
}
if (exists($pmods{"M1c"})){
cnt_file($input, $outputM1c, $treef, "M1c", $contfilM1c);
`codeml $contfilM1c`;
}
if (exists($pmods{"M2c"})){
cnt_file($input, $outputM2c, $treef, "M2c", $contfilM2c);
`codeml $contfilM2c`;
}

});

});

sub cnt_file {
	my ($seqf, $outf, $trf, $mod, $cntrp )  = @_ ;
	if ($mod eq "M2a1"){$modelD = 2; $NSistesD = 2; $fix_oD = 1; $omegaD = 1; $fix_rD = 1;}
	elsif ($mod eq "M2a"){$modelD = 2; $NSistesD = 2; $fix_oD = 0; $omegaD = 3; $fix_rD = 1;}
	elsif ($mod eq "M1"){$modelD = -2; $NSistesD = 0; $fix_oD = 0; $omegaD = 0.3; $fix_rD = 1;}
	elsif ($mod eq "M0"){$modelD = 0; $NSistesD = "0 1 2"; $fix_oD = 0; $omegaD = 0; $fix_rD = 0;}
#	elsif ($mod eq "M2c"){$modelD = 3; $NSistesD = 2; $fix_oD = 0; $omegaD = 0; $fix_rD = 0;}
	open(my $ocntr, '>', $cntrp) or die "Cannot open $cntrp: $!";
	print $ocntr "seqfile = $seqf\noutfile = $outf\ntreefile = $trf\nnoisy = 9\nverbose = 0\nrunmode = 0\nseqtype = 1\nCodonFreq = 2\nclock = 0\naaDist = 0\naaRatefile = wag.dat\nmodel = $modelD\nNSsites = $NSistesD\nicode = 0\nMgene = 0\nfix_kappa = 0\nkappa = 2\nfix_omega = $fix_oD\nomega = $omegaD\nfix_alpha = 1\nalpha = 0\nMalpha = 0\nncatG = 4\nfix_rho = $fix_rD\nrho = 0\ngetSE = 0\nRateAncestor = 0\nSmall_Diff = .5e-6\ncleandata = 0\nfix_blength = 0\nmethod = 0";
	close $ocntr;
    
}
