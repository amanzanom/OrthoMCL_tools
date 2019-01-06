#!/usr/bin/perl

use Getopt::Long;
use Parallel::ForkManager;
use File::Path;
use File::Copy qw(copy);

# Define subroutines
sub printVersion {
	# Usage: printVersion(<programName>, <version>, <typeglobRef>);
	# Example: printVersion('program_name', '3.2', \*STDERR);
	# Capture arguments
	my $software= $_[0];
	my $version= $_[1];
	my $fileHandle= $_[2];
	
	print $fileHandle "$software v$version\n";
	exit (0);
}

sub runSysCmd {
	# Usage: runSysCmd(<command>, <typeglobRef>);
	# Example: runSysCmd('blastp ...', \*STDERR);
	# Capture arguments
	my $cmd= $_[0];
	my $fileHandle= $_[1];
	
	print $fileHandle $cmd . "\n";
	
	return(system($cmd));
}

sub splitFASTA {
	# Usage: splitFASTA(<inFile.fasta>, <numSeq>, <numParts>);
	# Example: splitFASTA('inFile.fasta', 130000, 20);
	# Capture arguments
	my $inFASTA= $_[0];
	my $numSeq= $_[1];
	my $numParts= $_[2];

	# Define variables
	my $numSeqPerFile= 0;
	my @splittedFiles= ();
	my $fileCount= 0;
	my $numSeqRead= 0;
	my $splitFileName= '';
	my $inFileH= undef;
	my $splitFileH= undef;

	mkdir('temp_files');
	$numSeqPerFile= $numSeq/$numParts;
	if ($numSeqPerFile=~ m/\./){
		$numSeqPerFile=~ s/\..+$//;
		$numSeqPerFile++;
	}

	open($inFileH, '<', $inFASTA) || die "Unable to open file for reading: $inFASTA\n$!\n";
	while (my $line=<$inFileH>){
		if ($line=~ m/^>/){
			$numSeqRead++;
			if ($numSeqRead/($numSeqPerFile+1)==1){
				close($splitFileH);
				$numSeqRead= 1;
			}
			if ($numSeqRead==1){
				$fileCount++;
				$splitFileName= 'temp_files/tempFasta_' . $fileCount . '.fasta';
				push(@splittedFiles, $splitFileName);
				open($splitFileH, '>', $splitFileName) || die "Unable to open file for writting: $splitFileName\n$!\n";
			}
		}
		print $splitFileH $line;
		if (eof){
			close($splitFileH);
		}
	}
	close($inFileH);

	return(\@splittedFiles);
}

sub getGeneListFASTA {
	# Usage: getGeneListFASTA(<inFile.fasta>);
	# Example: getGeneListFASTA('inFile.fasta');
	# Capture arguments
	my $inFASTA= $_[0];

	# Define variables
	my $inFileH= undef;
	my $orgID= '';
	my $protID= '';
	my %geneList= ();

	open($inFileH, '<', $inFASTA) || die "Unable to open file for reading: $inFASTA\n$!\n";
	while (my $line=<$inFileH>){
		if ($line=~ m/^>(\S+)\|(\S+)/){
			$orgID= $1;
			$protID= $2;
			$geneList{$orgID}{$protID}= 1;
		}
	}
	close($inFileH);

	return(\%geneList);
}

sub addStrainSpecific {
	# Usage: addStrainSpecific(<groupFile.txt>, <hashRef>);
	# Example: addStrainSpecific('groups_1.5.txt', $geneList_ref);
	# Capture arguments
	my $groupFile= $_[0];
	my $geneList_ref= $_[1];

	# Define variables
	my $inFileH= undef;
	my $outFileH= undef;
	my $line= '';
	my $orthoClust= '';
	my $tempString= '';
	my $orgID= '';
	my $protID= '';
	my $orthoClust_prefix= '';
	my $orthoClust_suffix= '';
	$groupFile=~ m/(\S+)(\.txt)/;
	my $groupFile_ss= $1 . '_ss' . $2;
	
	open($inFileH, '<', $groupFile) || die "Unable to open file for reading: " . $groupFile . "\n$!\n";
		while (my $line=<$inFileH>){
			chomp $line;
			if ($line=~ m/^(\S+)\: (.+)$/){
				$orthoClust= $1;
				$tempString= $2;
				while ($tempString=~ m/([^\|\s]+)\|([^\|\s]+)/g){
					$orgID= $1;
					$protID= $2;
					delete($geneList_ref->{$orgID}{$protID});
				}
			}
		}
		close($inFileH);
		
		copy($groupFile, $groupFile_ss);
		
		$orthoClust=~ m/(\S+)_(\d+)/;
		$orthoClust_prefix= $1;
		$orthoClust_suffix= $2;

		open($outFileH, '>>', $groupFile_ss) || die "Unable to open file for writting: " . $groupFile_ss . "\n$!\n";
		foreach $orgID (sort {lc($a) cmp lc($b)} keys %{$geneList_ref}){
			foreach $protID (sort {lc($a) cmp lc($b)} keys %{$geneList_ref->{$orgID}}){
				print $outFileH $orthoClust_prefix . '_' . (++$orthoClust_suffix) . ': ' . $orgID . '|' . $protID . "\n";
			}
		}
}



GetOptions(\%opts, "config=s", "I=f@", "prefix=s", "sspecific!", "threads|t=i", "plength:i", "stop:i", "remove|r!", "help|h!");
        sub usage(){
                die "USAGE :: OrthoMCLPipe.pl -config filename -blast blast_command_line -I INT [-plength INT] [-stop INT] [-remove] [-help]\n\n
                config\t\tName of orthomcl config file [String]\n
                I\t\tInflation value(s) for mcl algorithm [Float >=1.2, <=5]\n
                prefix\t\tPrefix to use for cluster of orthologous groups [String]\n
                sspecific\t\tCreate file with strain specific genes added [Boolean]\n
                threads\t\tNumber of simultaneous blast threads to run (accomplished splitting input) [Integer >0]\n
                plength\t\tProtein minimum length (default 10) [Integer > 0]\n
                stop\t\tMaximum percent of stop codons in protein seuqences (default 20) [Integer >=0, <=100]\n
                remove|r\t\tWeather to remove or not files from database [Boolean]\n
                help|h\t\tPrint usage manual\n\n
                ";
        }
        
if ($opts{'help'} || !$opts{'config'} || !$opts{'prefix'} || !$opts{'I'}){
	&usage;
}

if ($opts{'sspecific'}){
	$addStrainSpecific= 1;
}
else {
	$addStrainSpecific= 0;
}

##### Define default values if needed#####
if (!$opts{'plength'}){
	$plength=10;
}
else {
	$plength=$opts{'plength'};
}
if (!$opts{'stop'}){
	$stop=20;
}
else {
	$stop=$opts{'stop'};
}
if ($opts{'remove'}){
	$r="all";
}
else {
	$r="no";
}
$ORTHOMCLBIN='/home/manzanoa/software/PHYLOGENY/ORTHOLOGUES/OrthoMCL/bin';

##### Capture DB username, password and name from config file#####
open (CONFIG, "$opts{'config'}") || die ("Unable to open: $opts{'config'}\n$!\n");
while (<CONFIG>){
	chomp $_;
	if ($_=~ /^#/){
		next;
	}
	elsif ($_=~ m/^percentMatchCutoff=(\d+)/){
		$perMatchCutoff= $1;
	}
	elsif ($_=~ /^dbConnectString=(\w+):(\w+):(\w+)/){
		$dB=$3;
	}
	elsif ($_=~ /^dbLogin=(\w+)/){
		$userName=$1;
	}
	elsif ($_=~ /^dbPassword=(\w+)/){
		$password=$1;
	}
}
close (CONFIG);

##### Create database in mySQL. If already exists, delete. #####
runSysCmd("mysql --user=".$userName." --password=".$password." -e 'DROP DATABASE IF EXISTS " . $dB . "; CREATE DATABASE " . $dB . "'", \*STDERR);

##### Begin Orthomcl Pipe described in http://orthomcl.org/common/downloads/software/v2.0/UserGuide.txt#####

###### Step 4 orthomclInstallSchema
runSysCmd("$ORTHOMCLBIN/orthomclInstallSchema ".$opts{'config'}." install_schema.log", \*STDERR);

####### Step 6 orthomclFilterFasta
runSysCmd("$ORTHOMCLBIN/orthomclFilterFasta compliantFasta/ ".$plength." ".$stop." ", \*STDERR);

###### Step 7 All-v-all BLAST
runSysCmd("makeblastdb -in goodProteins.fasta -input_type fasta -dbtype prot -out goodProteins_blastDB", \*STDERR);

$numProt=`grep ">" -c goodProteins.fasta`;
chomp $numProt;
print STDERR $cmd . "\n";

$splitFastaRef= splitFASTA('goodProteins.fasta', $numProt, $opts{'threads'});

$parallelBLAST = new Parallel::ForkManager($opts{'threads'});
$outFileBLASTCount= 0;
foreach $splitFastaFile (@$splitFastaRef){
	$outFileBLASTCount++;
	my $pid= $parallelBLAST-> start and next;
	my $cmd='blastp -task blastp -outfmt 6 -soft_masking true -seg yes -max_target_seqs 10000 -evalue 1e-5 -query ' . $splitFastaFile . ' -db goodProteins_blastDB -dbsize ' . $numProt . ' > temp_files/goodProteins_BLAST_part' . $outFileBLASTCount;
	print STDERR $cmd . "\n";
	system ($cmd);
	$parallelBLAST-> finish;
}
$parallelBLAST-> wait_all_children;

runSysCmd("cat temp_files/goodProteins_BLAST_part* > goodProteins_BLAST.out", \*STDERR);
unlink glob ("temp_files/* temp_files/*.*");
rmdir ("temp_files");

###### Step 8 orthomclBlastParser
if (-e 'similarSequences.txt'){
	unlink ('similarSequences.txt');
}
runSysCmd($ORTHOMCLBIN .'/orthomclBlastParser goodProteins_BLAST.out compliantFasta/ | awk \'$1!=$2 && $8>=' . $perMatchCutoff . '\' >> similarSequences.txt', \*STDERR);

###### Step 9 orthomclLoadBlast
runSysCmd("$ORTHOMCLBIN/orthomclLoadBlast ".$opts{'config'}." similarSequences.txt", \*STDERR);

###### Step 10 orthomclPairs
runSysCmd("$ORTHOMCLBIN/orthomclPairs ".$opts{'config'}." orthomcl_pairs.log cleanup=".$r, \*STDERR);

###### Step 11 orthomclDumpPairsFiles
if (-d 'pairs'){
	rmtree('pairs');
}
runSysCmd("$ORTHOMCLBIN/orthomclDumpPairsFiles ".$opts{'config'}, \*STDERR);

###### Step 11.5 (optional) add strain specific genes
if ($addStrainSpecific){
	$geneList_ref= getGeneListFASTA('goodProteins.fasta');
}

####### Step 12 & 13: Run all iterations with all Inflation values ######
foreach $I (@{$opts{'I'}}){
	
	###### Step 12 mcl
	runSysCmd("mcl mclInput --abc -I ".$I." -o mclOutput_".$I, \*STDERR);
	
	###### Step 13 orthomclMclToGroups
	runSysCmd("$ORTHOMCLBIN/orthomclMclToGroups " . $opts{'prefix'} . " 0001 < mclOutput_" . $I . " > groups_" . $I . ".txt", \*STDERR);
	
	if ($addStrainSpecific){
		addStrainSpecific(("groups_" . $I . ".txt"), $geneList_ref);
	}
}




