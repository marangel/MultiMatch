#!/bin/bash

function indel_search
{
	#work in a temporary directory
	
	mkdir temp
	cd temp
	readcounter=$8
	clippedcutoff=$5
	outputindels=$6
	outputcoords=$7
	tr '\r' '\n' < ${10} > temprefseq.fasta
	sed -i '1d;:a;N;$!ba;s/\n//g' temprefseq.fasta
	refseq=`pwd`/temprefseq.fasta
	
	if [ "$4" == "b2" ];
	then
		bowtieMode=4
	else
		bowtieMode=500
	fi
		
	if [ "${11}" == "1" ];
	then
		qualitycutoff=${12}
		qualityOutput=${13}
		alternativeCoords="all" #this can be an integer, but don't change it unless you know what it is for
	fi
	export readcounter clippedcutoff intermediatefastq intermediatesam qualitycutoff refseq qualityOutput storecounts firstIter alternativeCoords outputindels outputcoords
	
	if [ "$1" != "unmapped" ];
	then
		#unzip consensus and alignment -if zipped-, or copy them to the temp directory
		gzip -t $1 2>/dev/null
		[[ $? -eq 0 ]]&&gunzip -c $1 > Consensus.fastq || cp $1 ./Consensus.fastq
		gzip -t $2 2>/dev/null
		[[ $? -eq 0 ]]&&gunzip -c $2 > Aligned.sam || cp $2 ./Aligned.sam

		unmaped=$9
		export readlen unmaped
		#obtain unmaped reads using UnmmapedFilter from MultiMatch_CirSeqModule perl module
		perl -e '
		use MultiMatch_CirSeqModule_v7 qw(UnmmapedFilter);
		UnmmapedFilter(`readlink -f Consensus.fastq`, `readlink -f Aligned.sam`, $ENV{unmaped}, $ENV{readcounter});
		'
	elif [ "$1" == "unmapped" ];
	then
		unmaped=$2
	fi
	
	i=1
	while [ "$i" != "0" ]
	do
		
		intermediatesam=AlignmentIntermediate_$i.sam		
		if [ "$i" == "1" ];
		then
			
			#primary alignment requesting long seed and with a high cost for gaps 
			bowtie2 --local --phred33 -L 30 --rdg 40,10 --rfg 40,10 --gbar $bowtieMode -x $3 -U $unmaped -S $intermediatesam &> /dev/null
			firstIter=1
		
		elif [ -s $intermediatefastq ];
		then

			#secondary alignment requesting short seed and with a high cost for gaps 
			bowtie2 --local --phred33 -L 20 --rdg 40,10 --rfg 40,10 --gbar $bowtieMode -x $3 -U $intermediatefastq -S $intermediatesam &> /dev/null
			firstIter=0
			
		else

			break
		
		fi
		
		#create a new fastq file containing reads made from the soft-clipped bases that were mapped in the first alignment. Min read length is indicated by $5
		intermediatefastq=TrimmedReads_$i.fastq
		perl -e '
			use MultiMatch_CirSeqModule_v7 qw(ReadRearranger);
			ReadRearranger($ENV{intermediatesam}, $ENV{intermediatefastq}, $ENV{clippedcutoff}, $ENV{firstIter} );
		'
		
		i=$(( i+1 ))
		
	done
	rm AlignmentIntermediate_1.sam
	
	#Calculate coordinates and indel lengths using CoordsCollector from CirSeq_Indels_plUtilities_v1 perl module
	perl -e'
	use MultiMatch_CirSeqModule_v7 qw(CoordsCollector);
	CoordsCollector(`ls AlignmentIntermediate* | sort -V`, $ENV{refseq}, $ENV{outputindels}, $ENV{outputcoords});
	'

	if [ "${11}" == "1" ];
	then

		storecounts="storeCounts"
		sams=( $(find ./ -type f -name "AlignmentIntermediate_*" | wc -l) )
		sams=$(( sams+1 ))
		for i in `seq 2 $sams`
		do

			intermediatesam=AlignmentIntermediate_$i.sam
			if [ "$i" == "$sams" ];
			then

				storecounts="RemoveDouble"
			fi

			perl -e '
			use MultiMatch_CirSeqModule_v7 qw(ByQualityCounter);
			ByQualityCounter($ENV{intermediatesam}, $ENV{qualitycutoff}, $ENV{qualityOutput}, $ENV{refseq}, $ENV{storecounts}, "indel_basecount", $ENV{outputindels}, $ENV{alternativeCoords});
			'

		done

	fi

	if [ "$1" != "unmapped" ];
	then
		#zip unmaped reads
		gzip $unmaped
	fi
	cd ..
	
	mv ./temp/AlternativeIndelSites.txt ./
	#remove temporary files
	rm -r temp

}


mapconsensus=0
basecountconsensus=0
mapindel=0
basecountindel=0
regex='^[0-9]+$'
while [[ $# -gt 1 ]]
do
	opt="$1"
	case $opt in
		--mapconsensus)
		mapconsensus=1
		pathFastq="$2"
		pathIndex="$3"
		cutoffMisMatch="$4"
		shift 3
		;;
		--basecountconsensus)
		basecountconsensus=1
		phredCutoff="$2"
		pathReference="$3"
		shift 2
		;;
		--mapindel)
		mapindel=1
		if [[ $2 == unmapped ]];
		then
			pathFastq="$2"
			pathSam="$3"
			pathIndex="$4"
			shift 3
		elif ! [[ $2 =~ $regex ]];
		then
			pathFastq="$2"
			pathSam="$3"
			pathIndex="$4"
			shift 3
		fi
		cutoffSclipped="$2"
		callBowtieIndels="$3"
		shift 2
		if [ -z ${pathReference+x} ];
		then
			pathReference="$2"
			shift
		fi
		;;
		--basecountindel)
		basecountindel=1
		phredCutoffIndel="$2"
		;;
		*)
		echo -e "Unknown option $opt provided, please see README file: aborting run\n"
		exit 1
		;;
	esac
shift
done




if [[ $mapconsensus == 1 ]];
then
	
	echo "Consensus Mapping: Started"

	pathSam=`pwd`/AlignedConsensus.sam
	export cutoffMisMatch pathSam

	bowtie2 --no-unal --no-hd --no-sq --local --phred33 -L 30 --rdg 40,40 --rfg 40,10 --gbar 500 -x $pathIndex -U $pathFastq -S Alignment_1.sam &> /dev/null

	perl -e'
	use MultiMatch_CirSeqModule_v7 qw(MultiMatch);
	MultiMatch("Alignment_1.sam", $ENV{pathSam}, "Rearranged.fastq", $ENV{cutoffMisMatch});
	'

	bowtie2 --no-unal --no-hd --no-sq --local --phred33 -L 30 --rdg 40,40 --rfg 40,10 --gbar 500 -x $pathIndex -U Rearranged.fastq -S Alignment_2.sam &> /dev/null
	#bowtie2 --no-unal --local --phred33 -L 20 -x $1 -U $2 -S Alignment_1.sam &> /dev/null

	rm Rearranged.fastq Alignment_1.sam

	perl -e'
	use MultiMatch_CirSeqModule_v7 qw(MultiMatch);
	MultiMatch("Alignment_2.sam", $ENV{pathSam}, "Rearranged.fastq", $ENV{cutoffMisMatch});
	'
	rm Rearranged.fastq Alignment_2.sam
	echo "Consensus Mapping: Finished"

fi

if [[ $basecountconsensus == 1 ]];
then

	echo "Base counting, Quality threshold "$phredCutoff": Started"

	tr '\r' '\n' < $pathReference > temprefseq.fasta
	sed -i '1d;:a;N;$!ba;s/\n//g' temprefseq.fasta
	refseq=`pwd`/temprefseq.fasta
	qOutput=`pwd`/CountedBases_phred$phredCutoff\_consensus.txt
	lenOutput=`pwd`/ReadLengthDistribution.txt
	export phredCutoff qOutput refseq lenOutput
	perl -e '
	use MultiMatch_CirSeqModule_v7 qw(ByQualityCounter);
	ByQualityCounter($ENV{pathSam}, $ENV{phredCutoff}, $ENV{qOutput}, $ENV{refseq}, $ENV{lenOutput});
	'
	rm temprefseq.fasta
	echo "Base counting: Finished"
		
fi

if [[ $mapindel == 1 ]];
then	
	echo "Indel Search: Started"
	indel_search $pathFastq $pathSam $pathIndex $callBowtieIndels $cutoffSclipped `pwd`/Indels_List.txt `pwd`/Indels_Coordinates.txt `pwd`/Read_Stats.txt `pwd`/Unmapped_Reads.fastq $pathReference $basecountindel $phredCutoffIndel `pwd`/CountedBases_phred$phredCutoffIndel\_indel.txt
	echo "Indel Search: Finished"
fi

if [[ $mapconsensus == 1 ]];
then
	gzip $pathSam
fi
