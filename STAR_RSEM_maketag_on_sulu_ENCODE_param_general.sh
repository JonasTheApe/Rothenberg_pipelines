#!/bin/sh
set -e # Exit immediately if a command exits with a non-zero status.
#######################################################################
#######################################################################
##
##   This Script is designed to STAR map files on Sulu against a chosen reference genome and then process this with RSEM before taking all the redistributed reads (bam-file) creates a Homer tagdirectory that could be used for downstream analysis. Also count and fpkm result files are generated directly from rsem that could be used in downstream analysis. The parameters for STAR and RSEM are derived from the ENCODE pipeline https://github.com/ENCODE-DCC/long-rna-seq-pipeline/blob/master/DAC/STAR_RSEM.sh.

##
##   author: Jonas Ungerbäck
##   date:   January 25, 2016
##
##   Software used:  STAR 2.4.0, RSEM 1.2.25, Homer, Samtools, FastQC, Trimmomatic 0.33 with TruSeq3 adapter sequences, R (needs to be set in the PATH-variable).


echo "#######################################################################"
echo "#######################################################################"
echo "##"
echo "##   This Script is designed to STAR map files on Sulu against a chosen reference genome and then process this with RSEM before taking all the redistributed reads (bam-file) creates a Homer tagdirectory that could be used for downstream analysis. Also count and fpkm result files are generated directly from rsem that could be used in downstream analysis. The parameters for STAR and RSEM are derived from the ENCODE pipeline https://github.com/ENCODE-DCC/long-rna-seq-pipeline/blob/master/DAC/STAR_RSEM.sh."

echo "##"
echo "##   author: Jonas Ungerbäck"
echo "##   date:   January 25, 2016"
echo "##"
echo "##   Software used:  STAR 2.4.0, RSEM 1.2.25, Homer, Samtools, FastQC, Trimmomatic 0.33 with TruSeq3 adapter sequences, R (needs to be set in the PATH-variable)."

echo ""
echo ""

#Variables to be changed by the user
threads=10
refdir_size_mm9="/rothenberg/sulu/data00/chromInfo/mm9.chr.sizes"
star_refdir_mm9="/rothenberg/sulu/data00/star-indexes/mm9_NCBIM37.66"
rsem_refdir_mm9="/rothenberg/sulu/data00/star-indexes/rsem-reference/rsem_mm9_NCBIM37.66/rsem_NCBIM37.66_gtf"
refdir_size_mm10="/rothenberg/sulu/data00/chromInfo/mm10.chr.sizes"
star_refdir_mm10="/rothenberg/sulu/data00/star-indexes/mm10_gencode_vM8_gtf"
rsem_refdir_mm10="/rothenberg/sulu/data00/star-indexes/rsem-reference/rsem_mm10_gencode_vM8_gtf/rsem_mm10_gencode_vM8_gtf"
# Public folder for the bigwig tracks
bwtrack_dir="/rothenberg/sulu/home/hiroyuki/public_html/bigwigs"


# Fixed starting variables
privrna="0"
strandedness="0"
seqtype="0"
trimstatus="0"
directory="0"
refgenome="0"
process_check="0"
adapters="0"


## Set the path to the folder where the sample folders are and press enter. Loop to test the path?

        while [ ! -d ${privrna} ] || [ ! -d ${directory} ]; do

            echo "Set the absolute path to your sample folder i.e '/rothenberg/sulu/data01/common/jonun/rnaseq/'"
           echo "You must type something or the script will fail!!!!"
                read -e privrna
                echo "These are the items in your chosen directory:"
                ls -lh ${privrna}

           echo "Enter the samples you wish to process with a space in between and press enter:"
               read samples

            # check if the directory exists
            echo "${privrna}"
            [ -d ${privrna} ] && echo "and there is a directory called '${privrna}'."
            [ -d ${privrna} ] || echo "${privrna} does not exist."
               for i in ${samples}; do
                   directory=${privrna}${i}

            # check if the samples exist in ${privrna}.
                    echo "$directory"
                    [ -d $directory ] && echo "and there is a directory with the name ${i} in ${privrna}"
                    [ -d $directory ] || echo "${i} does not exist in ${privrna}. Start over and check your path and samples."
                    [ -d $directory ] || continue
               done

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


echo "######################################################################################################"
echo ""
echo "You will begin to process your ${samples} samples ${samples} in ${privrna} and your data is '${strandedness}' and '${seqtype}'-end sequenced. You have selected ${refgenome} as your reference genome"

if [ ${trimstatus} == yes ]; then
    echo "Your data will be trimmed with Trimmomatic with ${adapters}-adapters."
else
    echo "Your data will not be trimmed with Trimmomatic."
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

for j in $samples; do
        directory=${privrna}${j}

        # go to the folder that contains fastq files
        cd $directory

       if [ ! -d ${j}${suffix} ]; then

        echo "Processing ${j}!"

        mkdir ${j}.rsem


        echo "Unzipping and catenating ${j}.fastq.gz-files."


        # Unzip the data (will depend on if the data is single and paired-end.

        cd ${j}.rsem

_mydir="`pwd`"

echo "You are now in: $_mydir"

            mkdir ${j}.rsem_processed_files
            mkdir ${j}.for_genome_browser

        case "${seqtype}" in
        single)
            #OPTION: single-end data
            gunzip -c ${privrna}${j}/*.fastq.gz > ${privrna}${j}/${j}.rsem/${j}.fastq
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
            gunzip -c ${privrna}${j}/*R1*.fastq.gz > ${privrna}${j}/${j}.rsem/${j}.R1.fastq
            gunzip -c ${privrna}${j}/*R2*.fastq.gz > ${privrna}${j}/${j}.rsem/${j}.R2.fastq
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

            java -jar /rothenberg/sulu/data00/Trimmomatic/trimmomatic-0.33.jar ${Trimmode} -threads ${threads} -trimlog ${j}.trimlog.txt ${TRIMfastqin} ${TRIMfastqout} ILLUMINACLIP:${Trimadapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
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

        /rothenberg/sulu/data00/star/./STAR --genomeDir ${STAR_REF_DIR} --readFilesIn ${STARfastq} --runThreadN ${threads}  --quantMode TranscriptomeSAM --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 ${STARparStrand} ${STARparWig} --outFileNamePrefix ${j}.${refgenome}. --outReadsUnmapped Fastx


        ########################################################################################

        #### Bigwig generation from the sorted STAR output sam-file
        echo "Bigwig generation from the STAR outputted sam-file"

        #Convert sam to bam
samtools view -bS -o ${j}.STAR.bam ${j}.${refgenome}.Aligned.out.sam # Here maybe I can pipe it via sort the the file comes out sorted already in this step.
        # convert bam to bed
        echo "Creating and sorting a bed-file from ${j}.STAR.bam."
        # convert bam to bed
        bamToBed -i ${j}.STAR.bam  -split > ${j}.STAR.bed
        # sort bed-file
        bedSort ${j}.STAR.bed ${j}.STAR.${refgenome}.bed
        # get normaliztion factor for bigwig file
        scalevar_STAR=$( cat ${j}.STAR.${refgenome}.bed | wc -l )
        scalevar_STAR=$( echo "scale=4; 1000000/$scalevar_STAR" | bc )
        echo "**Inverse of Number of Reads (in Millions)**  $scalevar_STAR"
        # compile reads
        echo "Generating ${j}.STAR.${refgenome}.bg"
        genomeCoverageBed -bg -i ${j}.STAR.${refgenome}.bed -g ${REF_DIR_SIZE} -scale $scalevar_STAR -split > ${j}.STAR.${refgenome}.bg
        #convert bedgraph to bigwig file
        echo "Generating ${j}.STAR.${refgenome}.bw"
        bedGraphToBigWig ${j}.STAR.${refgenome}.bg ${REF_DIR_SIZE} ${j}.STAR.${refgenome}.bw
        mv ${j}.STAR.${refgenome}.bw ${j}.STAR.${refgenome}.bg ${j}.STAR.${refgenome}.bed ${j}.for_genome_browser
        #Creating a symbolic link of the file to the
        ln -s "$directory/${j}.rsem/${j}.for_genome_browser/${j}.STAR.${refgenome}.bw" ${bwtrack_dir}

#170601 Add a picard rna-seq metric step (in the long run maybe as a function?). Above a reflat choice must be made depending on genome. Maybe one should add it conditionally with and if or a case.
# Also add that the bam file is saved insted of the ${j}.STAR.${refgenome}.bed file.
        #Cleaning up intermediate files
        rm ${j}.STAR.bed ${j}.STAR.bam

        ###########################################################################################



    #### prepare for RSEM: sort transcriptome BAM to ensure the order of the reads, to make RSEM output (not pme) deterministic. This is for single-end read data and the pipe is taken from the ENCODE script.

        echo "sorting ${j}.${refgenome}.Aligned.toTranscriptome.out.bam to ensure the order of the reads, to make RSEM output (not pme) deterministic."
        #mv ${j}.${refgenome}.Aligned.toTranscriptome.out.bam Tr.bam
        mv ${j}.${refgenome}.Aligned.toTranscriptome.out.bam Tr.bam


    ### Decides how the data is sorted before RSEM processing. This will depend if the data is single end or paired-end sequenced as well as stranded or unstranded.

    case "$dataType" in
        "stranded|single"|"unstranded|single")
        # single-end data
        cat <( samtools view -H Tr.bam ) <( samtools view -@ 10 Tr.bam | sort -T ./ ) | samtools view -@ 10 -bS - > ${j}.${refgenome}.Aligned.toTranscriptome.out.bam
    ;;
        "stranded|paired"|"unstranded|paired")
        # paired-end data, merge mates into one line before sorting, and un-merge after sorting
        cat <( samtools view -H Tr.bam ) <( samtools view -@ 10 Tr.bam  | awk '{printf "%s", $0 " "; getline; print}' | sort -T ./ | tr ' ' '\n' ) | samtools view -@ 10 -bS - > ${j}.${refgenome}.Aligned.toTranscriptome.out.bam
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

    /rothenberg/sulu/data00/rsem/rsem-calculate-expression -p ${threads} --output-genome-bam --sampling-for-bam --estimate-rspd --seed 12345 ${RSEMparType} --calc-ci --append-names --bam ${j}.${refgenome}.Aligned.toTranscriptome.out.bam ${RSEM_REF_DIR} ${j}.Aligned.toTranscriptome.out


        echo "rsem-plot-model ${j}.Aligned.toTranscriptome.out ${j}.Aligned.toTranscriptome.out.pdf"
    /rothenberg/sulu/data00/rsem/rsem-plot-model ${j}.Aligned.toTranscriptome.out ${j}.Aligned.toTranscriptome.out.pdf


        # create TagDirectory and keeping all the reads from the genome generated bam file
        echo ""
        echo "Creating a Homer tag directory keeping all reads from ${j}.Aligned.toTranscriptome.out.genome.sorted.bam!"

        suffix="_STAR_RSEM_${refgenome}_tagdir"
        echo "makeTagDirectory ${j}${suffix} ${j}.Aligned.toTranscriptome.out.genome.sorted.bam -keepAll -format sam"

        /rothenberg/sulu/data00/homer/bin/makeTagDirectory ${j}${suffix} ${j}.Aligned.toTranscriptome.out.genome.sorted.bam -keepAll -format sam


        ### To get a proper scoring with Homer run analyzeRepeats.pl with rna and -count genes and -noadj options. Maybe -count exons would be more accurate?


    #### Bigwig generation from the rsem AlignedToTranscriptome.bam
        echo "Bigwig generation from the rsem AlignedToTranscriptome.bam."

        # convert bam to bed

        echo "Creating and sorting a bed-file from ${j}.Aligned.toTranscriptome.out.genome.sorted.bam."
        # convert bam to bed
        bamToBed -i ${j}.Aligned.toTranscriptome.out.genome.sorted.bam -split > ${j}.bed
        # sort bed-file
        bedSort ${j}.bed ${j}.${refgenome}.s.bed
        # get normaliztion factor for bigwig file
        scalevar=$( cat ${j}.${refgenome}.s.bed | wc -l )
        scalevar=$( echo "scale=4; 1000000/$scalevar" | bc )
        echo "**Inverse of Number of Reads (in Millions)**  $scalevar"
        # compile reads
        echo "Generating ${j}.${refgenome}.bg"
        genomeCoverageBed -bg -i ${j}.${refgenome}.s.bed -g ${REF_DIR_SIZE} -scale $scalevar -split > ${j}.${refgenome}.bg
        #convert bedgraph to bigwig file
        echo "Generating ${j}.${refgenome}.bw"
        bedGraphToBigWig ${j}.${refgenome}.bg ${REF_DIR_SIZE} ${j}.STAR.rsem.${refgenome}.bw



        ## Cleaning up and moving files
        echo "Starting to clean up!"
        rm ${j}.${refgenome}.s.bed
        rm ${j}.Aligned.toTranscriptome.out.genome.bam
        rm ${j}.Aligned.toTranscriptome.out.transcript.bam
        rm ${j}.${refgenome}.Aligned.out.sam
        rm ${j}.bed
        rm Tr.bam


        echo "More cleaning up."
        mv ${j}.${refgenome}.bg ${j}.Aligned.toTranscriptome.out.genome.bam ${j}.for_genome_browser
        mv *.results ${j}.rsem_processed_files/
        mv ${j}.Aligned.toTranscriptome.out.pdf ${j}.rsem_processed_files/
        mv *.bai *.out *.mate1 *.tab ${j}.rsem_processed_files/
        mv *.bam ${j}.rsem_processed_files/
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
