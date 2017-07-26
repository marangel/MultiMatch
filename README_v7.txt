MultiMatch CirSeq Module
Mauricio Aguilar-Rangel
Frydman Lab, 2017

------------------------------------------- *** -------------------------------------------

MultiMatch is a set of scripts written in bash and perl that can be used for the detection of indels in CirSeq data sets (Ref to CirSeq paper),
as described in (Ref to our paper). 

------------------------------------------- *** -------------------------------------------

Before running MultiMatch:
	->MultiMatch runs under Perl 5, make sure you have it (you should, unless you have updated it to version 6)
	->make sure that bowtie2 is installed and added to your PATH
	->export to PERL5LIB the directory on which the MultiMatch perl module was saved, e.g. if the .pm file is located in /home/usr/Documents/Programs then:
		export PERL5LIB=/home/usr/Documents/Programs

------------------------------------------- *** -------------------------------------------

Running MultiMatch:

NOTE: when the arguments for the options that are explained below require to supply a "path to some file" you must use absolute paths
(e.g. /home/usr/Documents/Myfile.txt instead of ~/Documents/Myfile.txt)

options
	--mapconsensus [pathToConsensusSequences] [pathToBowtie2index] [MisMatchCutoff]
		this option can be used to align the CirSeq consensus sequences produced by ConsensusGeneration.py, a script that is part of the current CirSeq software
		distribution. It is faster than the aligning methods provided with the current CirSeq software distribution, and recovers more variant counts.
		Arguments:
		[pathToConsensusSequences] path to fastq consensus file (or fastq.gz)
		[pathToBowtie2index] path to bowtie2 index (including the prefix used in the  index)
		[MisMatchCutoff] an integer that will be used as the cutoff for the number of mismatches allowed in mapped reads.
		Example:
		--mapconsensus /home/usr/Documents/1_consensus.fastq.gz /home/usr/Documents/index/PolioVirus 2
		Output:
		It produces a file called AlignedConsensus.sam, located in the directory from where the command was called.
		
	--basecountconsensus [qualityCutoff] [pathToFastaReferenceFile]
		if --mapconsensus was used, this can be used to count the bases in AlignedConsensus.sam
		Arguments:
		[qualityCutoff] a phred quality cutoff (integer) to be used, bases with phredqual >= than this number will be counted
		[pathToFastaReferenceFile] the path to the reference fasta file
		Example:
		--basecountconsensus 20 /home/usr/Documents/PolioVirus.fasta
		Output:
		It produces a file called CountedBases_phred[qualityCutoff]_consensus.txt, where [qualityCutoff] is the supplied integer (e.g. CountedBases_phred20_consensus.txt)
		which contains the next columns:
			Position: nucleotide number (1 corresponds to the first nt in the provided reference fasta file)
			Reference: reference nucleotide
			A: mutations to A
			T: mutations to T
			C: mutations to C
			G: mutations to G
		It also outputs a file called ReadLengthDistribution.txt with the following columns:
		Length: read length
		Counts: number of times that read is present in the mapped consensus sequences
	
	--mapindel unmapped [pathToConsensusSequences or pathToCandidateReadsToLookForIndels] [pathToBowtie2index] [cutoffSoftClipped] [indelMappingMethod] [pathToFastaReferenceFile]
		is the core function of the module, whose purpose is to map deletions and duplications in a CirSeq data set. It can be invoked in stand-alone mode, or downstream of
		the --mapconsensus option, or --mapconsensus & --basecountconsensus.
		1) Invoked downstream of --mapconsensus:
		--mapindel [cutoffSoftClipped] [indelMappingMethod] [pathToFastaReferenceFile]
		2) Invoked downstream of --mapconsensus & --basecountconsensus:
		--mapindel [cutoffSoftClipped] [indelMappingMethod]
		3) Invoked in stand-alone mode, with consensus sequences and mapped consensus as input:
		--mapindel [pathToConsensusSequences] [pathToSamFile] [pathToBowtie2index] [cutoffSoftClipped] [indelMappingMethod] [pathToFastaReferenceFile]
		4) Invoked in stand-alone mode, with only consensus sequences as input:
		--mapindel unmapped [pathToConsensusSequences] [pathToBowtie2index] [cutoffSoftClipped] [indelMappingMethod] [pathToFastaReferenceFile]
			If you already have a fastq (or fastq.gz) file with the reads that didn't fully map to your reference sequence (i.e. candidates to have indels)
			you can use this option (recommended). As another option, this type of run can also take the consensus fastq (or fastq.gz) file as only input (no
			additional input sam file required as in number 3), and it will map the indels, but it will not retrieve a sam file with the mapped reads.
		
		Arguments:
		[cutoffSoftClipped] since MultiMatch is a split-read program, this integer indicates the minimum length of the fragments generated after read splitting
							(use a reasonable number, ~15)
		[indelMappingMethod] bowtie2 itself can align in a gapped fashion (mapping of short indels), however it relies in a gap penalty scoring, potentially increasing the
					false discovery rate of our method. If you set this option to "b2", multimatch will report multimatch identified indels as well as bowtie2 mapped indels.
					If set to "mm" indels will be mapped using multimatch only (this does not mean you'll necessary get less indels than if you use b2 method).
					The gap scoring parameters set for the "b2" mode are such that only very needed gaps will be reported (e.g. a long read with a single nt insertion in the middle).
					You are welcome to adjust those parameters in the indel_search function of the .sh file
		[pathToFastaReferenceFile] the path to the reference fasta file
		[pathToConsensusSequences] path to fastq consensus file (or fastq.gz)
		[pathToSamFile] path to the aligned consensus sequences in sam format (or sam.gz)
		[pathToBowtie2index] path to bowtie2 index (including the prefix used in the  index)
		Examples:
			1) --mapindel 15 mm /home/usr/Documents/PolioVirus.fasta
			2) --mapindel 15 mm
			3) --mapindel /home/usr/Documents/1_consensus.fastq.gz /home/usr/Documents/data.sam /home/usr/Documents/index/PolioVirus 15 mm /home/usr/Documents/PolioVirus.fasta
			4) --mapindel unmapped /home/usr/Documents/[1_consensus.fastq.gz or CandidateReads.fastq.gz] /home/usr/Documents/index/PolioVirus 15 mm /home/usr/Documents/PolioVirus.fasta
		Output:
		It produces three .txt files:
			Indels_List.txt contains the next columns:
				readID: ID of the read
				start: position at which the indel starts (inclusive)
				length: indel length, negative numbers for duplications, positive numbers for deletions
				method: how it was mapped (directly from bowtie2 or using multimatch)
				possible_aligning_options: if the sequences flanking the indel are full or partial repeats of it, then it's not possible to determine the exact
											position of the indel. Instead, one out of all the possible ones is reported (randomly chosen) and this colum
											tells the number of alternative aligning positions.
				repeat_number: if the previous column is non-zero, this one indicates  how many full repeats of the indel sequence were found
				flanking_left: number of nucleotides flanking the left end of the indel
				flanking_right: number of nucleotides flanking the right end of the indel
				mism_left: number of mismatches in the left-flanking read fragment (it is the number of mismatches in the entire read if the indel was identified with b2)
				mism_right: number of mismatches in the right-flanking read fragment (NA if the indel was identified using b2)
			Indels_Coordinates.txt contains the next columns:
				first column: read ID
				the rest of the columns: correspond to the starts and ends of the fragments onto which the original read was splitted. The first of those columns
				is the start of the first fragment, the second one is the end of the first fragment, the third the start of the second fragment, the fourth the end
				of the second fragment, and so on. (Not all the reads necesarily have the same number of start-end colums)
			Read_Stats.txt contains three columns:
				ReadsInputfastq: number of consensus reads in the input fastq file
				ReadsInputSam: number of consensus mapped reads in the provided sam file
				ReadsOutputfastq: number of consensus reads that will be used to look for indels
		When run in mode 1) or 2) it produces an additional file, called Unmmaped_Reads.fastq.gz, containing the reads that do not fully map to the reference sequence, i.e. the reads
		used to look for indels (candidate reads)
					
	--basecountindel [qualityCutoff]
		this can be used downstream of --mapindel to count the bases in the fragments flanking the identified indels (if any). Regions corresponding to insertions are counted
		only once.
		Example:
		--basecountindel 20
		Arguments:
		[qualityCutoff] a phred quality cutoff (integer) to be used, bases with phredqual >= than this number will be counted
		Output:
		It produces a file called CountedBases_phred[qualityCutoff]_indel.txt, where [qualityCutoff] is the supplied integer (e.g. CountedBases_phred20_indel.txt)
		which contains the next columns:
			Position: nucleotide number (1 corresponds to the first nt in the provided reference fasta file)
			Reference: reference nucleotide
			A: mutations to A
			T: mutations to T
			C: mutations to C
			G: mutations to G
		
------------------------------------------- *** -------------------------------------------

Example of a full run, i.e. the user has only the consensus sequences and wants to map mutations and indels, as well as base counting of the mapped reads:
MultiMatch_run_v5.sh --mapconsensus /home/usr/Documents/1_consensus.fastq.gz /home/usr/Documents/index/PolioVirus 2 --basecountconsensus 20 /home/usr/Documents/PolioVirus.fasta --mapindel 15 mm --basecountindel 20

------------------------------------------- *** -------------------------------------------
