#!/usr/bin/perl
while (@ARGV){
$_=shift @ARGV;
if ($_=~ /^-pdb$/){$pdb=shift @ARGV;}
elsif ($_=~ /^-ndb$/){$ndb=shift @ARGV;}
elsif ($_=~ /^-qry$/){$qry=shift @ARGV;}
elsif ($_=~ /^-nc$/){$nc=shift @ARGV;}
elsif ($_=~ /^-of$/){$out=shift @ARGV;}
}


if (not defined ($pdb && $qry && $out)){print "\nThis script annotate a fasta file to a protein local database (UniProt) and then create a table with results and GO terms.\n\nUsage:\nAGOUP.pl\n\t-pdb <path to protein database>\n\t-qry <path to fasta file of loci to be annotated>\n\t-of <path folder to save all outputs>\n\t-nc <number cores to run in parallel>\n\nFor example:\n\nAGOUP.pl -pdb /pygmy/uniprot.fasta -qry /pygmy/PPcandiloci.fasta -of /pygmy/blastout/ -nc 10 \n\n"; exit;}
if (not defined ($nc)) {$nc= 10;}
unless (-d $out){`mkdir $out`;} 

our $prefo = $qry;
$prefo =~ s/\.fasta//;
@prefc = split /\//, $prefo;
$pref = $prefc[-1];
our $qbase="$out/$pref";
##blast to protein database
`cat $qry | parallel -j $nc -k --block 3k --recstart '>' --pipe 'blastx -db $pdb -query - -evalue 1e-3 -outfmt 6' >> $qbase\.prot`;
#`awk '{if (\$3>=70 && \$11<0.001)print }' $qbase\.prot | sort -k 1,1 -k 11,11rn | awk '!seen[\$1]++' | awk '{print \$1, \$2}' | awk '!seen[\$2]++' > $qbase\_P.hit`;
`awk '{if (\$3>=70 && \$11<0.001)print }' $qbase\.prot | sort -k 1,1 -k 12,12rn | awk '!seen[\$1]++' > $qbase\_P.hit`;
`awk '{if (\$3>=70 && \$11<0.001)print }' $qbase\.prot | sort -k 1,1 -k 12,12rn | awk '!seen[\$1]++' | awk '{print \$1, \$2}' | awk '!seen[\$2]++' > $qbase\_PU.hit`;

##blast to nucleotide database
`cat $qry | parallel -j $nc -k --block 3k --recstart '>' --pipe 'blastn -db $ndb -query - -evalue 1e-3 -outfmt 6' >> $qbase\.gen`;
#`awk '{if (\$3>=70 && \$11<0.001)print }' $qbase\.gen | sort -k 1,1 -k 11,11rn | awk '!seen[\$1]++' | awk '{print \$1, \$2}' | awk '!seen[\$2]++' > $qbase\_P.hit`;
`awk '{if (\$3>=70 && \$11<0.001)print }' $qbase\.gen | sort -k 1,1 -k 12,12rn | awk '!seen[\$1]++' > $qbase\_G.hit`;
`awk '{if (\$3>=70 && \$11<0.001)print }' $qbase\.gen | sort -k 1,1 -k 12,12rn | awk '!seen[\$1]++' | awk '{print \$1, \$2}' | awk '!seen[\$2]++' > $qbase\_GU.hit`;
`awk 'FNR == NR { h[\$1]; next }; !(\$1 in h)' $qbase\_P.hit $qbase\_G.hit >> $qbase\_UG.hit`;



##Import the information of each gene
our $infgene="$out/genbank_info";
our $infgene2="$infgene/genbank_info";
unless (-d $infgene){`mkdir $infgene`;}
unless (-d $infgene2){`mkdir $infgene2`;}
#`cat  $qbase\_PU.hit | while read i; do feature=\$(echo \$i | cut -d " " -f 1 | perl -pe 's/\\|/\\\\|/g');  hit=\$(echo \$i | cut -d " " -f 2 | cut -d "|" -f 2); echo \"wget -q -O - http://www.uniprot.org/uniprot/\${hit}.txt\ > $infgene/\${feature}.info"; done > $out/gbcmd`;
`cat  $qbase\_PU.hit | while read i; do hit=\$(echo \$i | cut -d " " -f 2 | cut -d "|" -f 2); echo \"wget -q -O - http://www.uniprot.org/uniprot/\${hit}.txt\ >> $infgene/\${hit}.info"; done > $out/gbcmd`;
#`cat  $qbase\_UG.hit | while read i; do hit=\$(echo \$i | cut -d " " -f 2); echo \"wget -q -O - \\\"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=\${hit}&rettype=gb&retmode=text\\\" > $infgene2/\${hit}.info"; done > $out/gbcmd`;
#`cat  $qbase\_UG.hit | while read i; do hit=\$(echo \$i | cut -d " " -f 2); echo \"wget -q -O - \\\"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=\${hit}&rettype=gb&retmode=text\\\" | while read i; do if [[ \$i =~ "/product" ]]; then sname=$(echo $i | cut -d "=" -f 2 | sed 's/\"//g'); prot=$(grep -m 1 "$sname" ../db/Cetacea_UniProt.fasta | cut -d "|" -f 2); wget -q -O -http://www.uniprot.org/uniprot/${prot}.txt >> ./${prot}.info; break; fi; done

#`cat  $qbase\_UG.hit | while read i; do hit=\$(echo \$i | cut -d " " -f 2); echo \$(wget -q -O - \"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=\${hit}&rettype=gb&retmode=text\" | while read i; do if [[ \$i =~ "/product" ]]; then sname=$(echo $i | sed -e 's/\/product(.*\) \/protein_id/\1/'| cut -d "=" -f 2 | sed 's/\"//g'); prot=$(grep -m 1 "$sname" .../../prodb/EupercariaUniprot.fasta | cut -d "|" -f 2); wget -q -O -http://www.uniprot.org/uniprot/${prot}.txt >> ./${prot}.info; break; fi; done

`cat $out/gbcmd | parallel \"eval {}\"`;



##create annotation table
my $cmd="echo Locus\'\\t\'GeneAccession\'\\t\'MapPosition\'\\t\'GeneFullName\'\\t\'GeneShortName\'\\t\'GeneAltNames\'\\t\'GeneSubName\'\\t\'GeneGo\'\\t\'CellularComponent\'\\t\'MolecularFunction\'\\t\'BiologicalProcess >> $qbase\_annotation.csv";
system ($cmd);
my $cmd5="echo \\\"ID\\\"\$\'\\t\'\\\"Protein names\\\"\$\'\\t\'\\\"Gene names\\\" >> $qbase\_INF";
system ($cmd5);
our @hits=`cat  $qbase\_P.hit`;
foreach $hit (@hits){
	my $fname= "-";
	my $sname= "-";
	my $altname="-";
	my $subname="-";
	my $GOT="GO:";
	my $CC="CC:";
	my $MF="MF:";
	my $BP="BP:";
	my @hitinf=split /\t/, $hit;
	my @UP=split /\|/, $hitinf[1];
	my $UPa=$UP[1];
	our @geni=`cat $infgene/$UPa\.info`;
	foreach $line (@geni){
		chomp $line;
		if ($line =~ m/^DE   RecName: Full=/) {
			$fname = $line;
			$fname =~ s/DE   RecName: Full=//;
			$fname =~ s/\{.*//;
			$fname =~ s/\s+$//;
			$fname = quotemeta ($fname);
#			$fname =~ s/\>/\\>/g;
#			$fname =~ s/\'/\\'/g;
#			$fname =~ s/\,/\\,/g;
		}
		elsif ($line =~ m/^GN   Name=/) {
			$sname = $line;
			$sname =~ s/GN   Name=//;
			$sname =~ s/\{.*//;
			$sname =~ s/\s+$//;
			$sname = quotemeta ($sname);
#			$sname =~ s/\>/\\>/g;
#			$sname =~ s/\'/\\'/g;
#			$sname =~ s/\,/\\,/g;
		}
		elsif ($line =~ m/^DE   AltName: Full=/) {
			$altname = $line;
			$altname =~ s/DE   AltName: Full=//;
			$altname =~ s/\{.*//;
			$altname = quotemeta ($altname);
#			$altname =~ s/\>/\\>/g;
#			$altname =~ s/\'/\\'/g;
#			$altname =~ s/\,/\\,/g;
		}
		elsif ($line =~ m/^DE   SubName: Full=/) {
			$subname = $line;
			$subname =~ s/DE   SubName: Full=//;
			$subname =~ s/\{.*//;
			$subname = quotemeta ($subname);
#			$subname =~ s/\>/\\>/g;
#			$subname =~ s/\'/\\'/g;
#			$subname =~ s/\,/\\,/g;
		}
		elsif ($line =~ m/^DR   GO;/) {
			my $GOTi = $line;
			$GOTi =~ s/DR   GO;//;
			my @GOTin= split /;/, $GOTi;
			#$GOTin[0]=~ s/GO\://;
			$GOTin[0]=~ s/\s//g;
			$GOT= $GOT.$GOTin[0].",";
			if ($GOTin[1]=~ m/C:/){
			$GOTin[1]=~ s/C\://;
			$CC=$CC.$GOTin[1].",";
			
			}
			if ($GOTin[1]=~ m/F:/){
			$GOTin[1]=~ s/F\://;
			$MF=$MF.$GOTin[1].",";
			}
			if ($GOTin[1]=~ m/P:/){
			$GOTin[1]=~ s/P\://;
			$BP=$GOTin[1];
			}
		}
	}
$GOT =~ s/, $//;
$MF =~ s/\'/\\'/g;
$CC =~ s/\'/\\'/g;
$BP =~ s/\'/\\'/g;
my $cmd2="echo $hitinf[0]\'\\t\'$UP[-1]\'\\t\'$hitinf[8]-$hitinf[9]\'\\t\'$fname\'\\t\'$sname\'\\t\'$altname\'\\t\'$subname\'\\t\'$GOT\'\\t\'$CC\'\\t\'$MF\'\\t\'$BP >> $qbase\_annotation.csv";
#$cmd2 =~ s/\(/\\(/g;
#$cmd2 =~ s/\)/\\)/g;
$cmd2 =~ s/\;//g;
my $cmd3="echo $UP[-1]\'\\t\'$GOT >> $qbase\_GO.tab";
print "$cmd2\n";
system ($cmd2);
system ($cmd3);
my $cmd4="echo \\\"$UP[-1]\\\"\$\'\\t\'\\\"";
if ($fname ne "-"){
	$cmd4.="$fname\\\"\$\'\\t\'\\\"";
}elsif ($fname ne "-"){
        $cmd4.="$altname\\\"\$\'\\t\'\\\"";
}else{
	$cmd4.="$subname\\\"\$\'\\t\'\\\""; 
}
$cmd4.="$sname\\\" >> $qbase\_INF";
#$cmd4 =~ s/\(/\\(/g;
#$cmd4 =~ s/\)/\\)/g;
$cmd4 =~ s/\;//g;
}

##add not annotated genes table
our @hits2=`cat  $qbase\_UG.hit`;
foreach $hit2 (@hits2){
	my $fname= "-";
	my $sname= "-";
	my $altname="-";
	my $subname="-";
	my $GOT="GO:";
	my $CC="CC:";
	my $MF="MF:";
	my $BP="BP:";
	my @hitinf=split /\t/, $hit;
	# my @UP=split /\|/, $hitinf[1];
	my $UPa=$hiting[1];
	our @geni=`cat $infgene2/$UPa\.info`;
	foreach $line (@geni){
		chomp $line;
		if ($line =~ m/^DEFINITION  /) {
			$fname = $line;
			$fname =~ s/DEFINITION  //;
			$fname =~ s/\{.*//;
			$fname =~ s/\s+$//;
		}
		elsif ($line =~ m/\/gene=/) {
			$sname = $line;
			$sname =~ s/GN   Name=//;
			$sname =~ s/\{.*//;
			$sname =~ s/\s+$//;
		}
		elsif ($line =~ m/^DE   AltName: Full=/) {
			$altname = $line;
			$altname =~ s/DE   AltName: Full=//;
			$altname =~ s/\{.*//;
			$altname =~ s/\s+$//;
		}
		elsif ($line =~ m/^DE   SubName: Full=/) {
			$subname = $line;
			$subname =~ s/DE   SubName: Full=//;
			$subname =~ s/\{.*//;
			$subname =~ s/\s+$//;
		}
		elsif ($line =~ m/^DR   GO;/) {
			my $GOTi = $line;
			$GOTi =~ s/DR   GO;//;
			my @GOTin= split /;/, $GOTi;
			$GOTin[0]=~ s/GO\://;
			$GOT= $GOT.$GOTin[0].",";
			if ($GOTin[1]=~ m/C:/){
			$GOTin[1]=~ s/C\://;
			$CC=$CC.$GOTin[1].",";
			
			}
			if ($GOTin[1]=~ m/F:/){
			$GOTin[1]=~ s/F\://;
			$MF=$MF.$GOTin[1].",";
			}
			if ($GOTin[1]=~ m/P:/){
			$GOTin[1]=~ s/P\://;
			$BP=$GOTin[1];
			}
		}
	}
$MF =~ s/\'/\\'/g;
$CC =~ s/\'/\\'/g;
$BP =~ s/\'/\\'/g;
my $cmd3="echo $hitinf[0]\'\\t\'$UP[-1]\'\\t\'$hitinf[8]-$hitinf[9]\'\\t\'$fname\'\\t\'$sname\'\\t\'$altname\'\\t\'$subname\'\\t\'$GOT\'\\t\'$CC\'\\t\'$MF\'\\t\'$BP >> $qbase\_annotation.csv";
$cmd3 =~ s/\(/\\(/g;
$cmd3 =~ s/\)/\\)/g;d
$cmd3 =~ s/\;//g;
print "$cmd3\n";
system ($cmd3);
}
