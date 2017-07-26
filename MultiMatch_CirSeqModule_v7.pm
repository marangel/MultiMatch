package MultiMatch_CirSeqModule_v7;
use strict;
use warnings;
use List::Util qw( min max sum);
use Exporter qw(import);
use Storable;
 
our @EXPORT_OK = qw(UnmmapedFilter CoordsCollector ReadRearranger MultiMatch ByQualityCounter AlternativeIndelSites);

sub UnmmapedFilter {

	my ($fastq,$sam,$outfile,$outfile2) = @_;
	if (-e $outfile) {
		die "UnmmapedReads.fastq file already exists; please delete it, rename it or change current directory.\n";
	}
	open (FASTQ,$fastq) or die "Unable to open $fastq file: $!";
	open (SAM,$sam) or die "Unable to open $sam file: $!";
	open (OUT,'>',$outfile) or die "Unable to write an output file: $!";
	open (OUT2,'>>',$outfile2) or die "Unable to count reads in fastq and sam input: $!";
	
	my %ids;
	my $temp;
	while (<SAM>) {
		unless (/^@/ || /^(\S+)\t4\t/) {
			$temp = [split(/\t/)];
			$ids{$temp -> [0]} = 1;
		}
	}
	close SAM;
	
	my ($i, $j) = 0;
	while (<FASTQ>) {
		my $seq = <FASTQ>;
		my $plus = <FASTQ>;
		my $qual = <FASTQ>;
		if (/^@(\S+)/ & $seq !~ /^@/) {
			
			if (!exists($ids{$1})) {
				print OUT $_,$seq,$plus,$qual;
				$j++; #Counter for filtered reads.
			}
			$i++; #Counter for input-fastq reads.
		}
	}
	#Print to the second output file the  number of reads in: the input fastq file, the input sam file and the output fastq file.
	if (-z $outfile2) {
		print OUT2 "ReadsInputfastq\tReadsInputSam\tReadsOutputfastq\n";
	}
	print OUT2 "$i\t"; print OUT2 scalar keys %ids; print OUT2 "\t$j\n";
	close OUT2;
	close FASTQ;
	close OUT;
	
}

sub ReadRearranger {

	my ($sam, $fastq, $cutoff, $firstIteration) = @_;
	my ($i, $id, $cigar, $seq, $qual, $rname, $leftS, $middle, $rightS,);
	my @dels;
	open (SAM, $sam) or die "Unable to open $sam file: $!";
	open (OUT, '>', $fastq) or die "Unable to write an output file: $!";

	while( <SAM> ) {
		
		unless(/^@/ || /^(\S+)\t4\t/) {
		
			if(/^(\S+)\t(?:.+?\t){4}(\S+)\t(?:.+?\t){3}([A|C|G|T|N]+)[\t|\s]+(\S+)[\t|\s]+.*$/) {
				($id, $cigar, $seq, $qual) = ($1, $2, $3, $4);
				if( $cigar =~ /^((?:\d+S+))*((?:\d+[M|D|I])+)((?:\d+S+))*$/ && length($seq) == length($qual) ) { ##
					($leftS, $middle, $rightS) =  ($1, $2, $3);##
					if( $cigar =~ /(?:(\d+)D)+/g ) {
						@dels = $1;
					}
					push @dels, 0;
					unless( defined($leftS) ){
						$leftS = 0;
					}
					unless( defined($rightS) ){
						$rightS = 0;
					}
					$leftS = sum(split(/\D/,$leftS));
					$middle = sum(split(/\D/,$middle)) - sum(@dels);
					$rightS = sum(split(/\D/,$rightS));
					if( $firstIteration ) {
						if( $leftS && !$rightS ) {
							$seq = substr($seq, $leftS).substr($seq, 0, $leftS);
							$qual = substr($qual, $leftS).substr($qual, 0, $leftS);
						}elsif( !$leftS && $rightS ) {
							$seq = substr($seq, $middle).substr($seq, 0, $middle);
							$qual = substr($qual, $middle).substr($qual, 0, $middle);
						}
					}else{
						$seq = substr( $seq, $leftS+$middle, $rightS ).substr( $seq, 0, $leftS );
						$qual = substr( $qual, $leftS+$middle, $rightS ).substr( $qual, 0, $leftS );
					}
					
					if ( length($seq) >= $cutoff || $firstIteration ) {
						print OUT "\@$id\n$seq\n+\n$qual\n";
					}
					undef($leftS); undef($middle); undef($rightS);
					undef @dels;
				}
			}

		}
		
	}
	
	close SAM;
	close OUT;

}

sub MultiMatch {

	my ($sam, $aln, $clipped, $mmcutoff) = @_;
	my ($id, $cigar, $seq, $qual, $xm);
	
	open (SAM, $sam) or die "Unable to open $sam file: $!";
	open (ALN, '>>', $aln) or die "Unable to write aligned reads: $!";
	open (CLIP, '>>', $clipped) or die "Unable to write clipped reads: $!";

	while (<SAM>) {
		
		if(/^(\S+)\t(?:.+?\t){4}(\S+)\t(?:.+?\t){3}([A|C|G|T|N]+)[\t|\s]+(\S+)[\t|\s]+(?:.+?\t)+XM\:i\:(\d+).*$/) {
			($id, $cigar, $seq, $qual, $xm) = ($1, $2, $3, $4, $5);
			if ($cigar =~ /^\d+M$/ && $xm <= $mmcutoff) { #full-length alignments, mismatches present (no DISH) if ($6 <= 3 && $2 =~ /^\d+M$/)
				print ALN $_;
			}elsif($cigar =~ /^(\d+)S\d+M$/) { #aligned reads with soft clipped bases	
				$seq = substr($seq,$1).substr($seq,0,$1);
				$qual = substr($qual,$1).substr($qual,0,$1);
				print CLIP "\@$id\n$seq\n+\n$qual\n";
			}elsif($cigar =~ /^(\d+)M\d+S$/) { #aligned reads with soft clipped bases
				$seq = substr($seq,$1).substr($seq,0,$1);
				$qual = substr($qual,$1).substr($qual,0,$1);
				print CLIP "\@$id\n$seq\n+\n$qual\n";
			}
		}
		
	}
	close SAM;
	close ALN;
	close CLIP;

}

sub ByQualityCounter {
	my ($ref, $key, $i, $quality, $pos, $cigar, $seq, $qual, $sep, $j, $cumulative, $leftS, $middle, $sam, $qthreshold, $output, $refseq, $storing, $storingfile, $indellist, $id, $alternativesThreshold, $lenoutput);
	my (@refseq, @retrieved, @seq, @qual, @length, @operation, @dels, @indels, @counter);
	my (%ntcode, %ids);
	if( scalar(@_) == 5 ){
		($sam, $qthreshold, $output, $refseq, $lenoutput) = @_;
		$indellist = 0;
		$storingfile = "EsteArchivoNoExiste";
		$storing = "dontStoreCounts";
	}elsif( scalar(@_) == 8 ){
		($sam, $qthreshold, $output, $refseq, $storing, $storingfile, $indellist, $alternativesThreshold) = @_;
	}

	open (REF, $refseq) or die "Cannot open refseq file: $!";
	while (<REF>) {
		$ref = $_;
	}
	close REF;
	@refseq = split(//, $ref);
	
	if ( -e $storingfile && -s $storingfile ){
		@retrieved = @{ retrieve ( $storingfile ) };
		%ntcode = %{$retrieved[0]};
	}else{ #if empty, create empty hash
		$ntcode{$_} = [] for qw(A C G T N); 
		for $key (keys %ntcode) {
			for $i (0..(scalar(@refseq)+1)) {
				push @{$ntcode{$key}},[(0) x 50];
			}
		}
	}
	
	if ( scalar(@_) == 8 ) {
	
		open (LIST, $indellist) or die "Cannot open indel list file: $!";
		while( <LIST> ){
			unless (/^readID/) {
				@indels = split(/\t/);
				if($alternativesThreshold eq "all"){
					$ids{$indels[0]} = 1;
				}elsif($indels[4] <= $alternativesThreshold){
					$ids{$indels[0]} = 1;
				}

				if ($storing eq "RemoveDouble" && $indels[3] eq "mm" && $indels[2] < 0) {
					for $i (0..(abs($indels[2])-1)) {
						$ntcode{$refseq[ $indels[1]-1+$i ]}[$indels[1]+$i][$qthreshold]--;
					}
				}
			}
		}
		close LIST;
		$indellist = 1;
	
	}
	

	open (SAM, $sam) or die "Unable to open $sam file: $!";
	
	while ( <SAM> ) {
		
		unless(/^@/ || /^(\S+)\t4\t/) { #no header, no unmapped reads
		
			if(/^(\S+)\t(?:.+?\t){2}(\d+)\t\d+\t(\S+)\t(?:.+?\t){3}([A|C|G|T|N]+)[\t|\s]+(\S+)[\t|\s]+.*$/) {
				($id, $pos, $cigar, $seq, $qual) = ($1, $2, $3, $4, $5);
				if( !$indellist && length($seq) == length($qual) && $cigar !~ /[I|D|S|H]/g ){ #just matches
				$counter[length($seq)]++;
				@seq = split(//, $seq);
				@qual = split(//, $qual);
					foreach $i (0..$#seq) {
						$quality = ord($qual[$i])-33;
						if( $quality >= $qthreshold ) {
							$ntcode{$seq[$i]}[$pos+$i][$quality]++;
						}
					}
				}elsif( $indellist && exists($ids{$id}) && length($seq) == length($qual) ){
					if( $cigar =~ /^((?:\d+S+))*((?:\d+[M|D|I])+)(?:\d+S+)*$/ ) { ##
						($leftS, $middle) =  ($1, $2);##
						if( $cigar =~ /(?:(\d+)D)+/g ) {
							@dels = $1;
						}
						push @dels, 0;
						unless( defined($leftS) ){
							$leftS = 0;
						}
						$cigar = $middle;
						$leftS = sum(split(/\D/,$leftS));
						$middle = sum(split(/\D/,$middle)) - sum(@dels);
						$seq = substr( $seq, $leftS, $middle );
						$qual = substr( $qual, $leftS, $middle );
						undef($leftS); undef($middle);
						undef @dels;
						@seq = split(//, $seq);
						@qual = split(//, $qual);
						@operation = split(/\d+/, $cigar);
						shift @operation;
						@length = split(/\D+/, $cigar);
						$cumulative = -1;
						foreach $j (0..$#length) {
								if ($operation[$j] eq 'I') {
									$pos -= $length[$j];
									$cumulative += $length[$j];
								}elsif ($operation[$j] eq 'D') {
									$pos += $length[$j];
								}elsif ($operation[$j] eq 'M') {
									foreach $i ($cumulative+1..$length[$j]+$cumulative) {
										$quality = ord($qual[$i])-33;
										if( $quality >= $qthreshold ) {
											$ntcode{$seq[$i]}[$pos+$i][$quality]++;
										}
									}
									$cumulative += $length[$j];
								}
						}
					}
				}
			}
		
		}
	
	}
	
	close SAM;
	$sep->{$_} = "\t" for qw(A C G);
	$sep->{'T'} = "\n";
	if ( $storing eq "storeCounts" ) {
		store [\%ntcode], $storingfile;
	}else{
		open (OUT, '>', $output) or die "Impossible to write an output file :$!";
		if (-z $output) {
			print OUT "Position\tReference\tA\tC\tG\tT\n";
		}
		foreach $i (1..$#refseq+1) {
			print OUT "$i\t",$refseq[$i-1],"\t";
			foreach $key (sort keys %ntcode) {
				unless( $key eq 'N' ) {
					if( scalar( grep $_, @{$ntcode{$key}[$i]} ) > 0 ) {
						print OUT sum(@{$ntcode{$key}[$i]});
					}else{
						print OUT '0';
					}
					print OUT $sep->{$key};
				}
			}
		}
		
		close OUT;
		
		if($storingfile eq "EsteArchivoNoExiste"){
			open(OUT, '>', $lenoutput) or die "Cannot output the read fragment length distribution :$!";
			print OUT "Length\tCounts\n";
			foreach $i (0..$#counter){
				if(defined($counter[$i])){
					print OUT "$i\t$counter[$i]\n";
				}
			}
		}
	}

}

sub CoordsCollector {

	my ($sam, $key, $nelem, $id, $cigstr, $start, $seq, $cumulative, $j, $i, $first, $second, $flag, $alternatives, $coordstring, $md, $xm, $mdcounts, $cumulative2, $cumulative3, $key3, $repeats, $fragment, $indelen, $flanklen, $q, $t, $flanking, $ref, $key2, $indelString, $flankingString);
	my (%allsam, %mdtags, %alternativeIDs);
	my (@pos, @cigar, @index, @operation, @length, @ins, @dels, @alternative);
	
	local *CoordsCollector::AlternativeIndelSites = sub {
		$alternatives = 1;
		$repeats = 0;
		if( $index[2] > 0 ) {
			$fragment = substr($ref, $index[0], abs($index[2]));
			$flanking = substr($ref, 0, $index[0]);
			unshift @alternative, $index[0] + 1;
		}else{
			if( exists($mdtags{$key2}) && ( (defined(@{$mdtags{$key2}[$t]}[0]) && max(@{$mdtags{$key2}[$t]}) >= $pos[$t+1]) || (defined(@{$mdtags{$key2}[$t+1]}[0]) && min(@{$mdtags{$key2}[$t+1]}) <= $index[0]) ) ) {
				return;
			}
			$fragment = substr($ref, $pos[$t+1] - 1, abs($index[2]));
			$flanking = substr($ref, 0, $pos[$t+1] - 1);
			unshift @alternative, $pos[$t+1];
		}
		FLANKINGONE:
		$indelen = length($fragment);
		$flanklen = length($flanking);
		$q = -1;
		while( substr($flanking, $q) eq substr($fragment, $q) && $indelen && $flanklen ) {
			$alternatives++;
			$key3 = "id".$index[2]."-".($alternative[0] - 1);#($alternative[0] + $q);
			unless( exists($alternativeIDs{$key3}) ) {
				unshift @alternative, $alternative[0] - 1;
				#push @alternative, $alternative[0] + $q;
				$alternativeIDs{$key3} = 0;
			}
			$q--;
			$indelen--;
			$flanklen--;
		}
		if( !$indelen && $flanklen ) {
			$flanking = substr($flanking, 0, -length($fragment));
			$repeats++;
			goto FLANKINGONE;
		}elsif( !$indelen ) {
			$repeats++;
		}
	
		if( $index[2] > 0 ) {
			$flanking = substr($ref, $pos[$t+1] - 1);	
		}else{
			$flanking = substr($ref, $index[0]);
		}
		unshift @alternative, $alternative[$#alternative];
		pop @alternative;
		FLANKINGTWO:
		$indelen = length($fragment);
		$flanklen = length($flanking);
		$q = 1;
		while( substr($flanking, 0, $q) eq substr($fragment, 0, $q) && $indelen && $flanklen ) {
			$alternatives++;
			$key3 = "id".$index[2]."-".($alternative[$#alternative] + 1);#($alternative[0] + $q);
			unless( exists($alternativeIDs{$key3}) ) {
				push @alternative, $alternative[$#alternative] + 1;
				#push @alternative, $alternative[0] + $q;
				$alternativeIDs{$key3} = 0;
			}
			$q++;
			$indelen--;
			$flanklen--;
		}
		if( !$indelen && $flanklen ) {
			$flanking = substr($flanking, length($fragment));
			$repeats++;
			goto FLANKINGTWO;
		}elsif( !$indelen ) {
			$repeats++;
		}
		
		shift @alternative;
		if( scalar(@alternative) > 1 ) {
			$key3 = "id".$index[2]."-".$alternative[0];
			$alternativeIDs{$key3} = 0;
			@alternative = sort{ $a <=> $b } @alternative;
			print OUT3 $index[2], "\t", join("\t", @alternative), "\n";
		}
		
		print OUT "$indelString\t$alternatives\t$repeats\t$flankingString";		
		
		if( $t > -1 ) {
			if( exists($mdtags{$key2}) && defined(@{$mdtags{$key2}[$t]}[0]) ) {
				print OUT scalar( @{$mdtags{$key2}[$t]} ), "\t";
			}else{
				print OUT "0\t";
			}
			if( exists($mdtags{$key2}) && defined(@{$mdtags{$key2}[$t+1]}[0]) ) {
				print OUT scalar( @{$mdtags{$key2}[$t+1]} ), "\n";
			}else{
				print OUT "0\n";
			}
		}

	};

	#Read reference sequence
	open (REF, $_[$#_-2]) or die "Cannot open refseq file: $!";
	while (<REF>) {
		$ref = $_;
	}
	close(REF);

	open (OUT, '>>', $_[$#_-1]) or die "Cannot append to an output file: $!";
	if (-z $_[$#_-1]) {
		print OUT "readID\tstart\tlength\tmethod\tpossible_aligning_options\trepeat_number\tflanking_left\tflanking_right\tmism_left\tmism_right\n";
	}
	open (OUT2, '>>', $_[$#_]) or die "Cannot append to an output file: $!";
	
	open (OUT3, '>>', "AlternativeIndelSites.txt") or die "Cannot append to an output file: $!";
	
	$t = -1; $flankingString = "";
	foreach $sam (0..$#_-3) {
		
		open( SAM, $_[$sam] ) or die "Can't open",$_[$sam],"sam file: $!";
		while( <SAM> ) {
			
			unless( /^@/ || /^\S+\t4\t/ ) {
				if( /^(\S+)\t(?:.+?\t){2}(\d+)\t\d+\t(\S+)\t(?:.+?\t){3}([A|C|G|T|N]+)[\t|\s]+(?:.+?\t)+XM\:i\:(\d+)\t(?:.+?\t)+MD\:Z\:(\S+).*$/ ) {
				#if( /^(\S+)\t(?:.+?\t){2}(\d+)\t\d+\t(\S+)\t.*$/ ) {
					($id, $start, $cigstr, $seq, $xm, $md) = ($1, $2, $3, $4, $5, $6);
					if( exists($allsam{$id}) && $allsam{$id}[$#{$allsam{$id}}] =~ /^(?:\d+S+)*((?:\d+[M|D|I])+)(?:\d+S+)*$/ ) {
						$allsam{$id}[$#{$allsam{$id}}] = $1;
						push @{$allsam{$id}}, $start, $cigstr;
					}else{
						push @{$allsam{$id}}, $start, $cigstr;
					}
					
					
					if($cigstr =~ /[I|D]/) {
						@operation = split(/\d+/, $cigstr);
						shift @operation;
						@length = split(/\D+/, $cigstr);
						$cumulative = $start; $cumulative2 = 0; $cumulative3 = $start;
						$flag = 0;
						undef $key2;
						foreach $j (0..$#length) {
								if ($operation[$j] eq 'I') {
									if(substr($seq,$cumulative2 - $length[$j], $length[$j]) eq substr($seq,$cumulative2, $length[$j]) || substr($seq,$cumulative2 + $length[$j], $length[$j]) eq substr($seq,$cumulative2, $length[$j]) ) {#make sure it is a duplication
										if($cumulative3 != $start){
											print OUT $cumulative - $cumulative3 - $flag, "\t", $xm, "\tNA\n";
										}
										$indelString = "$id\t$cumulative\t".(-$length[$j])."\tb2";
										@index = ($cumulative + $length[$j] - 1, 0, -1*$length[$j]);
										@pos = ($cumulative);
										$key2 = $id;
										push @alternative, $cumulative;
										CoordsCollector::AlternativeIndelSites();
										undef @alternative;
										print OUT $cumulative - $cumulative3 - $flag, "\t";
										$cumulative3 = $cumulative;
										$flag = 0;
									}
									#$cumulative += $length[$j];
								}elsif ($operation[$j] eq 'D') {
									if($cumulative3 != $start){
										print OUT $cumulative - $cumulative3 - $flag, "\t", $xm, "\tNA\n";
									}
									$indelString = "$id\t$cumulative\t$length[$j]\tb2";
									@index = ($cumulative-1, 0, $length[$j]);
									@pos = ($cumulative + $length[$j]);
									$key2 = $id;
									push @alternative, $cumulative;
									CoordsCollector::AlternativeIndelSites();
									undef @alternative;
									print OUT $cumulative - $cumulative3 - $flag, "\t";
									$cumulative3 = $cumulative;
									$flag = $length[$j];
									$cumulative += $length[$j];
									$cumulative2 -= $length[$j];
								}elsif ($operation[$j] eq 'M') {
									$cumulative += $length[$j];
								}
								$cumulative2 += $length[$j];
						}
						if( defined($key2) ){print OUT $cumulative - $cumulative3 - $flag, "\t", $xm, "\tNA\n";}

					}
					
					if( $xm ) {
						($mdcounts, $md) = $md =~ /^(\d+)(.*)$/;
						$mdcounts += $start;
						if( $cigstr =~ /(?:(\d+)D)+/g ) {
							@dels = $1;
						}
						while( length($md) ) {
							if( $md =~ s/^([A|C|G|T]{1})(.*)$/$2/ ) {
								push @{$mdtags{$id}[$sam]}, $mdcounts;
								$mdcounts++;
							}elsif( $md =~ /^\^/ ) {
								$md =~ s/\^[A|C|G|T]{$dels[0]}(.*)$/$1/;
								$mdcounts += shift(@dels);
							}elsif( $md =~ s/^(\d+)(.*)$/$2/ ) {
								$mdcounts += $1;
							}
						}
					}
					
				}
			}
			
		}

	}
	
	undef @pos; undef @index; undef $key;
	foreach $key (keys %allsam) {
		$key2 = $key;
		$nelem = scalar(@{$allsam{$key}});
		if( $nelem >= 4 ) { 
			foreach (1..$nelem/2) {
				push @pos, shift @{$allsam{$key}};
				push @cigar, shift @{$allsam{$key}};
			}
			@index = sort{ $pos[$a] <=> $pos[$b] } 0..$#pos;
			@pos = @pos[ @index ];
			@cigar = @cigar[ @index ];
			if( exists($mdtags{$key}) ) {
				@{$mdtags{$key}} = @{$mdtags{$key}}[ @index ];
			}
			undef @index;
			$nelem = $nelem/2 - 1;
			$flag = 1;
			foreach $i (1..$nelem) {
				$t = $i-1;
				# if( (($first) = $cigar[$i-1] =~ /^(?:\d+S)*((?:\d+[M|D|I])+)$/) && (($second) = $cigar[$i] =~ /^((?:\d+[M|D|I])+)(?:\d+S)*$/) ) {
				if( (($first) = $cigar[$i-1] =~ /^((?:\d+[M|D|I])+)$/) && (($second) = $cigar[$i] =~ /^((?:\d+[M|D|I])+)$/) ) {
					push @index, $pos[$i-1] + sum( split(/\D+/, $first) ) - 1;
					if( $cigar[$i-1] =~ /(?:(\d+)I)+/g ) {
						@ins = $1;
						$index[0] -= sum(@ins);
					}
					push @index, $pos[$i] + sum( split(/\D+/, $second) ) - 1;
					if( $cigar[$i] =~ /(?:(\d+)I)+/g ) {
						@ins = $1;
						$index[1] -= sum(@ins);
					}
					push @index, (max(@pos[$i-1..$i]) - min(@index[0..1]) - 1);

					if( $flag == 1 ) {
						$coordstring = "$key\t$pos[$i-1]\t$index[0]\t$pos[$i]\t$index[1]";
						$flag++;
					}else{
						$coordstring = $coordstring."\t$pos[$i]\t$index[1]";
					}


					if( $index[2] > 0 ) {
						push @alternative, (min(@index[0..1]) + 1);
						$indelString = "$key\t$alternative[0]\t$index[2]\tmm";
						$flankingString = ($index[0] - $pos[$i-1] + 1)."\t".($index[1] - $pos[$i] + 1)."\t";
						CoordsCollector::AlternativeIndelSites();
						$flag = 0;
					}elsif( $index[2] < 0 ) {
						push @alternative, max(@pos[$i-1..$i]);
						$indelString = "$key\t$alternative[0]\t$index[2]\tmm";
						$flankingString = ($index[0] - $pos[$i-1] + 1)."\t".($index[1] - $pos[$i] + 1)."\t";
						CoordsCollector::AlternativeIndelSites();
						$flag = 0;
					}

				}
				undef @index; undef @alternative;
			}
			if( $flag < 1 ){
				print OUT2 $coordstring, "\n";
			}
		}
		undef @pos; undef @cigar;

	}
	close OUT;
	close OUT2;
	close OUT3;	

}

1;
