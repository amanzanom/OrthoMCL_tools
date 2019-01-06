use Data::Dumper;

open (GENELST, "<$ARGV[0]");
while ($line= <GENELST>){
	chomp $line;
	if ($line =~ m/^([^\|]+)\|([^\|]+)/){
		$org= $1;
		$gene= $2;
		$orgGene{$org}{$gene}= 1;
		if (!defined($orgHash{$org})){
			$orgHash{$org}=1;
		}
	}
}
close (GENELST);

#print Dumper %orgGene;

open (OMCLGROUP, "<$ARGV[1]");
while ($line= <OMCLGROUP>){
	chomp $line;
	if ($line=~ m/^(\S+)\: (.+)$/){	
		%tempHash= %orgHash;
		$flag{'notCore'}= 0;
		$orthoClust= $1;
		$tempString= $2;
		@geneArray= ();
		while ($tempString=~ m/([^\|\s]+)\|([^\|\s]+)/g){
			$org= $1;
			$gene= $2;
			push (@geneArray, "$org|$gene");
			$tempHash{$org}++;
		}
#		print Dumper %tempHash;$pene=<STDIN>;
		foreach $org (keys %tempHash){
			# single-copy core
			if ($tempHash{$org} != 2){
				$flag{'notCore'}= 1;
			}
			
			#core
#			if ($tempHash{$org} < 2){
#				$flag{'notCore'}= 1;
#			}
		}
		if (!$flag{'notCore'}){
#			print $line . "\n";
			print $orthoClust . ": " . join(" ", (sort {lc($a) cmp lc($b)} @geneArray)) . "\n";
		}

	}
}
close (OMCLGROUP);


