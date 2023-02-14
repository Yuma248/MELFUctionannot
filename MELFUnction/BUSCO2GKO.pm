package BUSCO2GKO;

my $LEVEL = 1;
sub getGKO{
my @arg = @_;

foreach $ar (@arg){
	if ($ar =~ /^-i/){our $input = (split(/ /,$ar))[1];}
	elsif($ar=~ /^-o/){our $output = (split(/ /, $ar))[1];}
	elsif ($ar=~ /^-snc/){$snc=(split(/ /, $ar))[1];}
	elsif ($ar=~ /^-tax/){$tax=(split(/ /, $ar))[1];}
}
if (not defined ($input)){print "\nThis script will extract the names of genes from busco codes and obtain the GO and KO terms for enrichment analyses. It requires a file with the list of codes one por row.\n\nUsage:\nBUSCO2GKO.pl\n\t-i <path to input list>\n\t-o <path to the output file, default GKO>\n Optional:\n\t-snc <number of runs in parallel, default 10>\n\t-tax <Taxon ID, sometimes the busco taxon ID is not recognized in other databases and you will have to change it, deafult the ID in the busco database you used>\n\t\n\nFor example:\nBUSCO2GKO -i /home/Yumafan/demultiplex/pamlgenes -o /home/Yumafan/nce -snc 10 -tax 9721\n\n"; exit;}
if (not defined ($output)){$output="./GKO";}
if (not defined ($snc)){$snc=10;}
my @busco = `cat $input | awk '{print \$1}'`; 
chomp @busco;
open(OUTFILE, "> $output") || die "could not open output file $output";
print OUTFILE "BUSCOID\tGENEID\tGO\tKO\n";
foreach $bg (@busco){
	if (not defined ($tax)){$tax= (split'at', $bg)[-1];}
	my $gID= `wget -q -O - \"https\:\/\/v10.orthodb.org\/tab\?id\=$bg\&species\=\" | head  -n 2 | tail -n 1 | awk -F '\t' '{print \$7}'`;
	$gID=~ s/ \n//g;
	chomp $gID;
	my $GO= `wget -qO- \"https://rest.uniprot.org/uniprotkb/search?query=gene\:$gID\+AND+taxonomy_id:$tax\&fields=go_id&format=tsv\" | grep -m 1 \"GO\:\"`;
	$GO=~ s/ //g;
	chomp $GO;
	my $KO= `wget -qO- \"https://rest.kegg.jp/find/genes/$gID\" | grep -w -m 1 \"$gID\" | awk \'{print \$1}\' | while read i; do wget -qO- \"https://rest.kegg.jp/link/ko/\${i}\" | cut -f 2,2; done`;
	$KO=~ s/ //g;
	chomp $KO;
	print OUTFILE "$bg\t$gID\t$GO\t$KO\n";
	
}
close OUTFILE;

}

1;
