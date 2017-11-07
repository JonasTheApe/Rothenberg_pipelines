#!/bin/sh
set -e # Exit immediately if a command exits with a non-zero status.
#######################################################################
#######################################################################
##
##   This Script is designed to STAR map files on Sulu against a chosen reference genome and then process this with RSEM before taking all the redistributed reads (bam-file) creates a Homer tagdirectory that could be used for downstream analysis. Also count and fpkm result files are generated directly from rsem that could be used in downstream analysis. The parameters for STAR and RSEM are derived from the ENCODE pipeline https://github.com/ENCODE-DCC/long-rna-seq-pipeline/blob/master/DAC/STAR_RSEM.sh.

##
##   author: Jonas Ungerbäck
##   date:   June 2, 2017
##
##   Software used:  STAR 2.4.0, RSEM 1.2.25, Homer, Samtools, FastQC, picard-tools 2.9.2, Trimmomatic 0.33 with TruSeq3 adapter sequences, R (needs to be set in the PATH-variable).

###########The script also need a text-file defining the sample-folders. This should be made in advance.


# Public folder for the bigwig tracks
bwtrack_dir="/rothenberg/sulu/home/jonun/public_html/For_UCSC_genome_browser/bigwigs"


echo "#######################################################################"
echo "#######################################################################"
echo "##"
echo "##   This Script is designed to STAR map files on Sulu against a chosen reference genome and then process this with RSEM before taking all the redistributed reads (bam-file) creates a Homer tagdirectory that could be used for downstream analysis. Also count and fpkm result files are generated directly from rsem that could be used in downstream analysis. The parameters for STAR and RSEM are derived from the ENCODE pipeline https://github.com/ENCODE-DCC/long-rna-seq-pipeline/blob/master/DAC/STAR_RSEM.sh."

echo "##"
echo "##   author: Jonas Ungerbäck"
echo "##   date:   June 25, 2017"
echo "##"
echo "##   Software used:  STAR 2.4.0, RSEM 1.2.25, Homer, Samtools, FastQC, Trimmomatic 0.33 with TruSeq3 adapter sequences, picard-tools with CollectRnaSeqMetrics v. 1.105, R (needs to be set in the PATH-variable)."

echo ""
echo ""

#Variables to be changed by the user
threads=10
refdir_size_mm9="/rothenberg/sulu/data00/chromInfo/mm9.chr.sizes"
star_refdir_mm9="/rothenberg/sulu/data00/star-indexes/mm9_NCBIM37.66"
rsem_refdir_mm9="/rothenberg/sulu/data00/star-indexes/rsem-reference/rsem_mm9_NCBIM37.66/rsem_NCBIM37.66_gtf"
refFlat_mm9="/rothenberg/sulu/home/jonun/Rothenberg_NGS_pipeline/refFlat_for_STAR_rna_script/mm9_refFlat_from_gtf"
refdir_size_mm10="/rothenberg/sulu/data00/chromInfo/mm10.chr.sizes"
star_refdir_mm10="/rothenberg/sulu/data00/star-indexes/mm10_gencode_vM8_gtf"
rsem_refdir_mm10="/rothenberg/sulu/data00/star-indexes/rsem-reference/rsem_mm10_gencode_vM8_gtf/rsem_mm10_gencode_vM8_gtf"
refFlat_mm10="/rothenberg/sulu/home/jonun/Rothenberg_NGS_pipeline/refFlat_for_STAR_rna_script/mm10_refFlat_from_gtf"


# Fixed starting variables
PATH_TO_SAMPLE_FILE="0"
strandedness="0"
seqtype="0"
trimstatus="0"
directory="0"
refgenome="0"
process_check="0"
adapters="0"
rnaseq_qc="0"
refFlat="0"
exit_status="false"


##Write the path to the sample text-file. NOTE: In the long run this file is meant to be a tab-delimited text file containing all the information in the questions the script now provides. Each option should then be able to be lifted out with awk. This is howver a future prospect.
echo "Enter the path to your text-file with your sample (sample definition file)."
while [ ${exit_status} != 0 ];do
read -e PATH_TO_SAMPLE_FILE
    # check if the directory exists
    echo "Checking so all the paths and folders are there."
    [ -f ${PATH_TO_SAMPLE_FILE} ] && echo "and there is a file-path called '${PATH_TO_SAMPLE_FILE}'."
        for SAMPLE in $(cat ${PATH_TO_SAMPLE_FILE}); do
            # check if the directory exists
            echo "${SAMPLE}"
            [ -d ${SAMPLE} ] && echo "and there is a directory called '$(basename ${SAMPLE})' in '${SAMPLE}'."
            [ -d ${SAMPLE} ] || echo "${SAMPLE} does not exist. Start over and check your path and samples."
            [ -d $SAMPLE ] || exit
        done
    [ -f ${PATH_TO_SAMPLE_FILE} ] && break
    [ -f ${PATH_TO_SAMPLE_FILE} ] || echo "${PATH_TO_SAMPLE_FILE} does not exist. Try to enter the path again!"

done


#Enter if the data is stranded or unstranded and test if what is entered is a valid option.
echo ""
        while [ "${strandedness}" != stranded ] || [ "${strandedness}" != unstranded ]; do
        echo "Is the data 'stranded' or 'unstranded'? Enter 'stranded' or 'unstranded'. You must type something or the script while fail."
            read strandedness
                if [ ${strandedness} = stranded ] || [ ${strandedness} = unstranded ]; then break
                else
                    echo "!!!!${strandedness} is an invalid option!!!!"
                fi
        done

#Enter if the data is single or paired-end read and test if what is entered is a valid option.
echo ""
        while [ "${seqtype}" != single ] || [ "${seqtype}" != paired ]; do
        echo "Is the data single or paired-end sequenced? Enter 'single' or 'paired'. You must type something or the script will fail."
            read seqtype
                if [ ${seqtype} = single ] || [ ${seqtype} = paired ]; then break
                    else
                        echo "!!!!${seqtype} is an invalid option!!!!"
                    fi
        done

#Choose reference genome, mm9 or mm10.
echo ""
        while [ "${refgenome}" != mm9 ] || [ "${refgenome}" != mm10 ]; do
            echo "Which reference genome do you want to use? 'mm9' or 'mm10'? You must type something or the script will fail."
                read refgenome
                if [ ${refgenome} = mm9 ] || [ ${refgenome} = mm10 ]; then break
                    else
                        echo "!!!!${refgenome} is an invalid option!!!!"
                    fi
        done

        case "${refgenome}" in
        mm9)
        #OPTION: mm9-genome-chosen
            REF_DIR_SIZE=${refdir_size_mm9}
            STAR_REF_DIR=${star_refdir_mm9}
            RSEM_REF_DIR=${rsem_refdir_mm9}
        ;;
        #OPTION: paired end data data
        mm10)
            REF_DIR_SIZE=${refdir_size_mm10}
            STAR_REF_DIR=${star_refdir_mm10}
            RSEM_REF_DIR=${rsem_refdir_mm10}
        ;;
        esac

        echo "Your reference files are in ${REF_DIR_SIZE} (genome sizes), ${STAR_REF_DIR} (STAR reference directory) and ${RSEM_REF_DIR} (RSEM reference directory)."

##Decide if your data should be trimmed with Trimmomatic and check if what has been chosen is a valip option.
echo ""
while [ "${trimstatus}" != yes ] || [ "${trimstatus}" != no ]; do
echo "Do you want to trim you sequences with Trimmomatic? 'yes' or 'no'. You must type something or the script will fail."
    read trimstatus
            if [ ${trimstatus} = yes ]  || [ ${trimstatus} = no ]; then break
                else
                    echo "!!!!${trimstatus} is an invalid option!!!!"
                fi
done

    if [ ${trimstatus} = yes ]; then


    #Enter what kind of sequencing adapters you have

    while [ ${adapters} != nextera ] || [ ${adapters} != trueseq ]; do
        echo "What type of adapters were used to build your reference library. 'nextera' or 'trueseq'. You must type something or the script will fail. "
        IFS= read adapters

        case "${adapters}" in
            trueseq)
            #OPTION: rnaseq processing
                    echo "You are using Trueseq-adapters."

            ;;
            #OPTION: chipseq processing
            nextera)
                    echo "You are using Nextera-adapters."

            ;;
        esac

        if [ ${adapters} = trueseq ] || [ ${adapters} = nextera ]; then break
        else
            echo "!!!!${seqdata} is an invalid option!!!!"
        fi

    done

    else
        echo "You have chosen not to trim your data" && continue
    fi

### Setting the variables based on the selections.


dataType="${strandedness}|${seqtype}"

    #The combinations for dataType that exist are the following:
        #stranded|single
        #stranded|paired
        #unstranded|single
        #unstranded|paired

echo ""

## Decide whether or the genome alignemnt from STAR should be processed with Picard rnaseq metric

while [ "${rnaseq_qc}" != yes ] || [ "${rnaseq_qc}" != no ]; do


echo "Do you want to run Picard RnaSeqMetrics on your STAR output. Note: This takes place prior to RSEM processing of the data"
read rnaseq_qc
if [ ${rnaseq_qc} = no ] || [ ${rnaseq_qc} = yes ]; then break
else
echo "!!!!${rnaseq_qc} is an invalid option!!!!"
fi
done


echo "######################################################################################################"
echo "######################################################################################################"
echo ""
echo "You will begin to process *** $(cat ${PATH_TO_SAMPLE_FILE}) *****."
echo ""
echo "Your data is '${strandedness}' and '${seqtype}'-end sequenced. You have selected ${refgenome} as your reference genome."

if [ ${trimstatus} == yes ]; then
    echo "Your data will be trimmed with Trimmomatic with ${adapters}-adapters."
else
    echo "Your data will not be trimmed with Trimmomatic."
fi

echo ""

if [ ${rnaseq_qc} == yes ]; then
    echo "Picard RnaSeqMetrics will be run on the STAR output"
else
    echo "Picard RnaSeqMetrics will NOT be run on the STAR output"
fi


##A final sanity check so everything is alright. If the user has misstyped something, it is here give the option to quit and start over.
echo ""
while [ "${process_check}" != yes ] || [ "${process_check}" != no ]; do
    echo "######################################################################################################"

    echo "Are you sure you have set all you samples and conditions properly. If not type 'no', quit and start over. Otherwise type 'yes' and proceed accordingly. You must type something or the script will fail."
    read process_check
    if [ ${process_check} = no ]; then exit
    elif [ ${process_check} = yes ]; then break
    else
        echo "!!!!${process_check} is an invalid option!!!!"
    fi
done



echo "######################################################################################################"

for SAMPLE in $(cat ${PATH_TO_SAMPLE_FILE}); do

    j=$(basename ${SAMPLE})

        # go to the folder that contains fastq files
        cd ${SAMPLE}

       if [ ! -d ${j}${suffix} ]; then


        echo "Processing ${j}!"

        mkdir -p ${j}.rsem



        # Unzip the data (will depend on if the data is single and paired-end.

        cd ${j}.rsem

        # This file should contain the last argument/error so if the script stops it should be reported where here.
        RUNLOG=${j}.runlog.txt

        echo "############################################"
        echo "Run started by $(whoami) on $(date)" > $RUNLOG
        echo "############################################"
        echo "Unzipping and catenating ${j}.fastq.gz-files."
        echo "############################################"


_mydir="`pwd`"



echo "You are now in: $_mydir"

            mkdir -p ${j}.rsem_processed_files
            mkdir -p ${j}.for_genome_browser

        case "${seqtype}" in
        single)
            #OPTION: single-end data
            gunzip -c ${SAMPLE}/*.fastq.gz > ${SAMPLE}/${j}.rsem/${j}.fastq
            Trimmode=SE

            case "${adapters}" in
            trueseq)
                Trimadapter=/rothenberg/sulu/data00/Trimmomatic/adapters/TruSeq3-SE.fa
            ;;
            nextera)
               Trimadapter=/rothenberg/sulu/data00/Trimmomatic/adapters/NexteraPE-PE.fa
            ;;
            esac
        ;;
            #OPTION: paired end data data
        paired)
            gunzip -c ${SAMPLE}/*R1*.fastq.gz > ${SAMPLE}/${j}.rsem/${j}.R1.fastq
            gunzip -c ${SAMPLE}/*R2*.fastq.gz > ${SAMPLE}/${j}.rsem/${j}.R2.fastq
            Trimmode=PE

            case "${adapters}" in
            trueseq)
                Trimadapter=/rothenberg/sulu/data00/Trimmomatic/adapters/TruSeq3-PE.fa
            ;;
                nextera)
                Trimadapter=/rothenberg/sulu/data00/Trimmomatic/adapters/NexteraPE-PE.fa
            ;;
            esac
        ;;
        esac

        ## FastQC check of the untrimmed file.

        case "${seqtype}" in
        single)
            #OPTION: single-end data
            /rothenberg/sulu/data00/FastQC/./fastqc -t ${threads} ${j}.fastq
            TRIMfastqin=${j}.fastq
            TRIMfastqout=${j}.trimmed.fastq
        ;;
        paired)
            #OPTION: paired end data data
            /rothenberg/sulu/data00/FastQC/./fastqc -t ${threads} ${j}.R1.fastq
            /rothenberg/sulu/data00/FastQC/./fastqc -t ${threads} ${j}.R2.fastq
            TRIMfastqin="${j}.R1.fastq ${j}.R2.fastq"
            TRIMfastqout="${j}.R1.paired.trimmed.fastq ${j}.R1.unpaired.trimmed.fastq ${j}.R2.paired.trimmed.fastq ${j}.R2.unpaired.trimmed.fastq"

        ;;
        esac


        ###### How data is trimmed will depend on if the data

        if [ ${trimstatus} == yes ]; then

            ## Trimming of the fastq-file. Values taken from http://www.usadellab.org/cms/?page=trimmomatic.
            echo "Trimmomatic is processing ${TRIMfastqin}"
            echo " java -jar /rothenberg/sulu/data00/Trimmomatic/trimmomatic-0.33.jar ${Trimmode} -threads ${threads} -trimlog ${j}.trimlog.txt ${TRIMfastqin} ${TRIMfastqout} ILLUMINACLIP:${Trimadapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

            java -jar /rothenberg/sulu/data00/Trimmomatic/trimmomatic-0.33.jar ${Trimmode} -threads ${threads} -trimlog ${j}.trimlog.txt ${TRIMfastqin} ${TRIMfastqout} ILLUMINACLIP:${Trimadapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>> $RUNLOG
            STARfastq=${TRIMfastqout}
        else
            echo "Moving on without trimming" 
            STARfastq=${TRIMfastqin}
        fi


            # FastQC on the trimmed data if trimming has been chosen

        if [ ${trimstatus} == yes ]; then
            case "${seqtype}" in
            single)
                #OPTION: single-end data
                /rothenberg/sulu/data00/FastQC/./fastqc -t ${threads} ${TRIMfastqout}
            ;;
            paired)
                #OPTION: paired end data data
                /rothenberg/sulu/data00/FastQC/./fastqc -t ${threads} ${j}.R1.paired.trimmed.fastq
                /rothenberg/sulu/data00/FastQC/./fastqc -t ${threads} ${j}.R2.paired.trimmed.fastq
                /rothenberg/sulu/data00/FastQC/./fastqc -t ${threads} ${j}.R1.unpaired.trimmed.fastq
                /rothenberg/sulu/data00/FastQC/./fastqc -t ${threads} ${j}.R2.unpaired.trimmed.fastq
                TRIMfastqout="${j}.R1.paired.trimmed.fastq ${j}.R2.paired.trimmed.fastq"
                STARfastq=${TRIMfastqout}
            ;;
            esac
        else
            echo "No trimmed file to run FastQC on!"
        fi

    ##STAR mapping with  some parameters depending on the seqtype.

        case "$dataType" in
        "stranded|single"|"stranded|paired")
            #OPTION: stranded data
            STARparStrand=""
            STARparWig="--outWigStrand Stranded"
        ;;
            #OPTION: unstranded data
            "unstranded|single"|"unstranded|paired")
            STARparStrand="--outSAMstrandField intronMotif"
            STARparWig="--outWigStrand Unstranded"
        ;;
        esac

    echo "Starting to map ${STARfastq} with STAR."

    echo "STAR --genomeDir ${STAR_REF_DIR} --readFilesIn ${STARfastq} --runThreadN ${threads}  --quantMode TranscriptomeSAM --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 ${STARparStrand} ${STARparWig} --outFileNamePrefix ${j}.${refgenome}. --outReadsUnmapped Fastx"

        /rothenberg/sulu/data00/star/./STAR --genomeDir ${STAR_REF_DIR} --readFilesIn ${STARfastq} --runThreadN ${threads}  --quantMode TranscriptomeSAM --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 ${STARparStrand} ${STARparWig} --outFileNamePrefix ${j}.${refgenome}. --outReadsUnmapped Fastx 2>> $RUNLOG


        ########################################################################################

        #### Bigwig generation from the sorted STAR output sam-file
        echo "Bigwig generation from the STAR outputted sam-file"

        #Convert sam to bam
        echo "samtools sort -o ${j}.STAR.bam -T deleteme -@ ${threads} ${j}.${refgenome}.Aligned.out.sam"
        samtools view -bS ${j}.${refgenome}.Aligned.out.sam | samtools sort - deleteme -o > ${j}.STAR.bam 2>> $RUNLOG
        samtools index ${j}.STAR.bam
        # convert bam to bed
        echo "Creating and sorting a bed-file from ${j}.STAR.bam."
        # convert bam to bed
        bamToBed -i ${j}.STAR.bam  -split > ${j}.STAR.bed 2>> $RUNLOG
        # sort bed-file
        bedSort ${j}.STAR.bed ${j}.STAR.${refgenome}.bed 2>> $RUNLOG
        # get normaliztion factor for bigwig file
        scalevar_STAR=$( cat ${j}.STAR.${refgenome}.bed | wc -l ) 
        scalevar_STAR=$( echo "scale=4; 1000000/$scalevar_STAR" | bc ) 
        echo "**Inverse of Number of Reads (in Millions)**  $scalevar_STAR" 
        # compile reads
        echo "Generating ${j}.STAR.${refgenome}.bg"
        genomeCoverageBed -bg -i ${j}.STAR.${refgenome}.bed -g ${REF_DIR_SIZE} -scale $scalevar_STAR -split > ${j}.STAR.${refgenome}.bg 2>> $RUNLOG
        #convert bedgraph to bigwig file
        echo "Generating ${j}.STAR.${refgenome}.bw"
        bedGraphToBigWig ${j}.STAR.${refgenome}.bg ${REF_DIR_SIZE} ${j}.STAR.${refgenome}.bw 2>> $RUNLOG
        #Creating a symbolic link of the file to the
        ln -s "$directory/${j}.rsem/${j}.for_genome_browser/${j}.STAR.${refgenome}.bw" ${bwtrack_dir} 

#170602 Added a picard rna-seq metric step (in the long run maybe as a function?). Above a reflat choice must be made depending on genome. These refflat files may have to be updated.
# Also add that the bam file is saved insted of the ${j}.STAR.${refgenome}.bed file.
        #Cleaning up intermediate files

        if [ ${rnaseq_qc} == yes ]; then

            mkdir ${j}.${refgenome}.rnaseq_qc
            case "${refgenome}" in
                mm9)
                refFlat=${refFlat_mm9}

            ;;
                mm10)
                refFlat=${refFlat_mm10}
            ;;
            esac

             echo "picard-tools CollectRnaSeqMetrics I=${j}.STAR.bam O=${j}.${refgenome}.rnaseq_qc/${j}_rnaseq_metric.txt REF_FLAT=${refFlat} STRAND_SPECIFICITY=NONE ASSUME_SORTED=FALSE CHART_OUTPUT=${j}.${refgenome}.rnaseq_qc/${j}_rnaseq_metric_chart.pdf"

            picard-tools CollectRnaSeqMetrics I=${j}.STAR.bam O=${j}.${refgenome}.rnaseq_qc/${j}_rnaseq_metric.txt REF_FLAT=${refFlat} STRAND_SPECIFICITY=NONE ASSUME_SORTED=FALSE CHART_OUTPUT=${j}.${refgenome}.rnaseq_qc/${j}_rnaseq_metric_chart.pdf 2>> $RUNLOG  #This strand_specificity assumption is just true as long as the data is unstranded and should probably be updated in the future.

        else
            echo ""
            echo "No RnaSeqMetrics will be collected."  $RUNLOG
            echo ""
        fi


        mv ${j}.STAR.bam ${j}.STAR.${refgenome}.bw ${j}.STAR.${refgenome}.bg ${j}.for_genome_browser
        rm ${j}.STAR.bed ${j}.STAR.${refgenome}.bed

        ###########################################################################################



    #### prepare for RSEM: sort transcriptome BAM to ensure the order of the reads, to make RSEM output (not pme) deterministic. This is for single-end read data and the pipe is taken from the ENCODE script.

        echo "sorting ${j}.${refgenome}.Aligned.toTranscriptome.out.bam to ensure the order of the reads, to make RSEM output (not pme) deterministic."

        mv ${j}.${refgenome}.Aligned.toTranscriptome.out.bam Tr.bam


    ### Decides how the data is sorted before RSEM processing. This will depend if the data is single end or paired-end sequenced as well as stranded or unstranded.

    case "$dataType" in
        "stranded|single"|"unstranded|single")
        # single-end data
        cat <( samtools view -H Tr.bam ) <( samtools view -@ 10 Tr.bam | sort -T ./ ) | samtools view -@ 10 -bS - > ${j}.${refgenome}.Aligned.toTranscriptome.out.bam 2>> $RUNLOG
    ;;
        "stranded|paired"|"unstranded|paired")
        # paired-end data, merge mates into one line before sorting, and un-merge after sorting
        cat <( samtools view -H Tr.bam ) <( samtools view -@ 10 Tr.bam  | awk '{printf "%s", $0 " "; getline; print}' | sort -T ./ | tr ' ' '\n' ) | samtools view -@ 10 -bS - > ${j}.${refgenome}.Aligned.toTranscriptome.out.bam 2>> $RUNLOG
    ;;
    esac


    # Set some options for RSEM depending on datatype.


    case "$dataType" in
    "stranded|single")
        #OPTION: stranded single end
        RSEMparType="--forward-prob 0"
    ;;
    "stranded|paired")
        #OPTION: stranded paired end
        RSEMparType="--paired-end --forward-prob 0"
    ;;
    "unstranded|single")
        #OPTION: stranded single end
        RSEMparType=""
    ;;
    "unstranded|paired")
        #OPTION: stranded paired end
        RSEMparType="--paired-end"
    ;;
    esac


    #RSEM, rsem-calculate-expression on the Aligned.toTranscriptome.out.bam generated by STAR and ouputs output-genome-bam for Homer tagdir generations.


        echo "rsem-calculate expression on ${j}.${refgenome}.Aligned.toTranscriptome.out.bam!"

        echo "/rothenberg/sulu/data00/rsem/rsem-calculate-expression -p ${threads} --output-genome-bam --sampling-for-bam --estimate-rspd --seed 12345 ${RSEMparType} --calc-ci --append-names --bam ${j}.${refgenome}.Aligned.toTranscriptome.out.bam ${RSEM_REF_DIR} ${j}.Aligned.toTranscriptome.out"

## Sampling for bam added 160331. Testrun on a few samples over the weekend.

    /rothenberg/sulu/data00/rsem/rsem-calculate-expression -p ${threads} --output-genome-bam --sampling-for-bam --estimate-rspd --seed 12345 ${RSEMparType} --calc-ci --append-names --bam ${j}.${refgenome}.Aligned.toTranscriptome.out.bam ${RSEM_REF_DIR} ${j}.Aligned.toTranscriptome.out 2>> $RUNLOG


        echo "rsem-plot-model ${j}.Aligned.toTranscriptome.out ${j}.Aligned.toTranscriptome.out.pdf"
    /rothenberg/sulu/data00/rsem/rsem-plot-model ${j}.Aligned.toTranscriptome.out ${j}.Aligned.toTranscriptome.out.pdf 2>> $RUNLOG


        # create TagDirectory and keeping all the reads from the genome generated bam file
        echo ""
        echo "Creating a Homer tag directory keeping all reads from ${j}.Aligned.toTranscriptome.out.genome.sorted.bam!"

        suffix="_STAR_RSEM_${refgenome}_tagdir"
        echo "makeTagDirectory ${j}${suffix} ${j}.Aligned.toTranscriptome.out.genome.sorted.bam -keepAll -format sam"

        /rothenberg/sulu/data00/homer/bin/makeTagDirectory ${j}${suffix} ${j}.Aligned.toTranscriptome.out.genome.sorted.bam -keepAll -format sam 2>> $RUNLOG


        ### To get a proper scoring with Homer run analyzeRepeats.pl with rna and -count genes and -noadj options. Maybe -count exons would be more accurate?


    #### Bigwig generation from the rsem AlignedToTranscriptome.bam
        echo "Bigwig generation from the rsem AlignedToTranscriptome.bam."

        # convert bam to bed

        echo "Creating and sorting a bed-file from ${j}.Aligned.toTranscriptome.out.genome.sorted.bam."
        # convert bam to bed
        bamToBed -i ${j}.Aligned.toTranscriptome.out.genome.sorted.bam -split > ${j}.bed 2>> $RUNLOG
        # sort bed-file
        bedSort ${j}.bed ${j}.${refgenome}.s.bed 2>> $RUNLOG
        # get normaliztion factor for bigwig file
        scalevar=$( cat ${j}.${refgenome}.s.bed | wc -l ) 
        scalevar=$( echo "scale=4; 1000000/$scalevar" | bc )
        echo "**Inverse of Number of Reads (in Millions)**  $scalevar" 
        # compile reads
        echo "Generating ${j}.${refgenome}.bg"
        genomeCoverageBed -bg -i ${j}.${refgenome}.s.bed -g ${REF_DIR_SIZE} -scale $scalevar -split > ${j}.${refgenome}.bg 2>> $RUNLOG
        #convert bedgraph to bigwig file
        echo "Generating ${j}.${refgenome}.bw"
        bedGraphToBigWig ${j}.${refgenome}.bg ${REF_DIR_SIZE} ${j}.STAR.rsem.${refgenome}.bw 2>> $RUNLOG



        ## Cleaning up and moving files
        echo "Starting to clean up!"
        rm ${j}.${refgenome}.s.bed
        rm ${j}.Aligned.toTranscriptome.out.genome.bam
        rm ${j}.Aligned.toTranscriptome.out.transcript.bam
        rm ${j}.${refgenome}.Aligned.out.sam
        rm ${j}.bed
        rm Tr.bam


        echo "More cleaning up."
        mv  ${j}.${refgenome}.bg ${j}.for_genome_browser
        mv *.results ${j}.rsem_processed_files/
        mv ${j}.Aligned.toTranscriptome.out.pdf ${j}.rsem_processed_files/
        mv *.bam *.bai *.out *.mate1 *.tab ${j}.rsem_processed_files/
        mv *Aligned.toTranscriptome.out.stat ${j}.rsem_processed_files/
        rm *fastq


####Testing if the size of the bigwig is correct (not to small). If it is the script will assume something is wrong and exit. Most often an erronous file have the value of 384 and then you know that something is wrong and potentially could carry over to the next round

        bigwig_file_size=$(cat ${j}.STAR.rsem.${refgenome}.bw  | wc -c) 
        echo "The size of the bigwig is: ${bigwig_file_size}."
        minimum_size=1000000
        if [ ! -f ${j}.STAR.rsem.${refgenome}.bw ] || [ "${bigwig_file_size}" -le "${minimum_size}" ]; then
            echo "${j}.STAR.rsem.${refgenome}.bw  is to small. Something is wrong with the script and the rest of your samples will not be processed." && exit
        else
            echo "${j}.STAR.rsem.${refgenome}.bw seems ok. The script will continue."
        fi

        echo "Moving ${j}.STAR.rsem.${refgenome}.bw to ${bwtrack_dir}."
        mv ${j}.STAR.rsem.${refgenome}.bw ${j}.for_genome_browser
        ln -s "$directory/${j}.rsem/${j}.for_genome_browser/${j}.STAR.rsem.${refgenome}.bw" ${bwtrack_dir}
        else
         echo "${j}${suffix} is already there. skipping..."
        fi
done
