use Data::Dumper;

@orgNames= ();
%temp=('string'=> undef
	);
%orgBin= ();

sub sortHashVal_byKey{
	# Usage: sortHashVal_byKey(<hashRef>, <sortString>, <boolean>, <boolean>);
	#        sortString is any of: -alpha (alpabetical)
	#                              -num (numerical)
	#                              -lex (lexical)
	#        first boolean refers to wheather or not sort case-insensitive (for alphabetical sorting)
	#        second boolean refers to wheather or not sort descending (or reverse)
	#        to sort case-insensitive and Descending, respectively
	# Example: printVersion('inFile.fasta', 'alpahbetical' );
	# Capture arguments
	my $hashRef= $_[0];
	my $sortString= $_[1];
	my $sortCaseIns= $_[2];
	my $sortDescend= $_[3];

	# Define variables
	my $key= undef;
	my @sortedKeysArray= ();
	my @valArray= ();
	
	if (!(scalar keys %$hashRef)){
		print STDERR "Empty hash!!!!, check usage\n";
		return(1);
	}
	
	if ($sortString eq 'alpha'){
		if ($sortCaseIns){
			@sortedKeysArray= sort {fc($a) cmp fc($b)} keys %$hashRef;
		}
		else {
			@sortedKeysArray= sort {$a cmp $b} keys %$hashRef;
		}
	}
	elsif ($sortString eq 'num'){
		@sortedKeysArray= sort {$a <=> $b} keys %$hashRef;
	}
	elsif ($sortString eq 'lex'){
		@sortedKeysArray= sort keys %$hashRef;
	}
	else {
		print STDERR "Unrecognized sort string!!!, check usage\n";
		return(1);
	}
	if ($sortDescend){
		@sortedKeysArray= reverse @sortedKeysArray;
	}
	
	for (my $i=0; $i<scalar(@sortedKeysArray); $i++){
		push(@valArray, $hashRef->{$sortedKeysArray[$i]})
	}
	
	return(\@valArray);
}

open ($orgLst_fH, '<', $ARGV[0]);
while ($line= <$orgLst_fH>){
	if ($line =~ m/^(\S+)/){
		push(@orgNames, $1);
	}
}
close ($orgLst_fH);

print STDOUT join("\t", @orgNames) . "\n";

open ($omclGroup_fH, '<', $ARGV[1]);
while ($line= <$omclGroup_fH>){
	chomp $line;
	if ($line=~ m/^(\S+)\: (.+)$/){
		$clustID= $1;
		$temp{'string'}= $2;
		print STDOUT $clustID . "\t";
		for (my $i=0; $i<scalar(@orgNames); $i++){
			$orgBin{$orgNames[$i]}= 0;
		}
		while ($temp{'string'}=~ m/([^\s\|]+)\|[^\s\|]+/g){
			$orgID= $1;
#			print STDOUT $orgID . "\n";
			if (defined($orgBin{$orgID}) && !$orgBin{$orgID}){
				$orgBin{$orgID}= 1;
			}
		}
#		print Dumper %orgBin;
		$tempLine= '';
		for (my $i=0; $i<scalar(@orgNames); $i++){
			$tempLine.= $orgBin{$orgNames[$i]} ."\t";
		}
		$tempLine=~ s/\t$/\n/;
		print STDOUT $tempLine;
#		print STDOUT join("\t", @{&sortHashVal_byKey(\%orgBin, 'alpha', 0, 0)}) . "\n";
#		print STDOUT $orgBin{$orgNames[$i]} ."\t";
	}
}
close ($omclGroup_fH);


