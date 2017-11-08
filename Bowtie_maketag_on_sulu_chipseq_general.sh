#!/bin/sh
set -ueo # Exit immediately if a command exits with a non-zero status.
#Testing git. Tjena fan hur e leget va?
#######################################################################
#######################################################################
##
##   This Script is designed to Bowtie1 map files on Sulu against a chosen reference genome. It will also remove duplicates and/or repeater masked regions if the user chooses to do so.  It will also create a Homer tagdirectory that could be used for downstream analysis.

##
##   author: Jonas Ungerbäck
##   date:   January 26, 2016
##
##   Software used:  Bowtie 1.1.1, Bedtools, Homer, Samtools set in your path, FastQC, Trimmomatic 0.33 with TruSeq3 adapter sequences, R (needs to be set in the PATH-variable).


echo "#######################################################################"
echo "#######################################################################"
echo "##"
echo "##   This Script is designed to Bowtie1 map files on Sulu against a chosen reference genome. It will also remove duplicates and/or repeater masked regions if the user chooses to do so.  It will also create a Homer tagdirectory that could be used for downstream analysis."

echo "##"
echo "##   author: Jonas Ungerbäck"
echo "##   date:   January 26, 2016"
echo "##"
echo "##   Software used:  Bowtie 1.1.1, Bedtools, Homer, Samtools set in your path, FastQC, Trimmomatic 0.33 with TruSeq3 adapter sequences, R (needs to be set in the PATH-variable)."
echo ""
echo ""


#Variables to be changed by the user
threads=10
refdir_size_mm9="/rothenberg/sulu/data00/chromInfo/mm9.chr.sizes"
bowtie_refdir_mm9="/rothenberg/sulu/data00/bowtie-indexes/mm9"

refdir_size_mm10="/rothenberg/sulu/data00/chromInfo/mm10.chr.sizes"
bowtie_refdir_mm10="/rothenberg/sulu/data00/bowtie-indexes/mm10"

repeat_masker_mm9="/rothenberg/sulu/data01/common/jonun/chipseq/mm9_UCSC_rmsk.bed"
repeat_masker_mm10="/rothenberg/sulu/data01/common/jonun/chipseq/mm10_UCSC_rmsk.bed"
black_list_mm9="/rothenberg/sulu/data01/common/jonun/chipseq/mm9-blacklist.s.bed"
black_list_mm10="/rothenberg/sulu/data01/common/jonun/chipseq/mm10-blacklist.s.bed"

# Public folder for the bigwig tracks
bwtrack_dir="/rothenberg/sulu/home/jonun/public_html/For_UCSC_genome_browser/bigwigs"


# Fixed starting variables
privchip="0"
strandedness="0"
seqtype="0"
trimstatus="0"
directory="0"
refgenome="0"
rmdup="0"
rmsk="0"
process_check="0"
correctlength="0"

suffix="_${rmdup}_${rmsk}_${refgenome}_tagdir"


## Set the path to the folder where the sample folders are and press enter. Loop to test the path?

        while [ ! -d ${privchip} ] || [ ! -d ${directory} ]; do

            echo "Set the absolute path to your sample folder i.e '/rothenberg/sulu/data01/common/jonun/chipseq/'"
            echo "You must type something or the script will fail!!!!"
                read -e privchip
                echo "These are the items in your chosen directory:"
                ls -lh ${privchip}
               

           echo "Enter the samples you wish to process with a space in between and press enter:"
               read samples

                    # check if the directory exists
                    echo "${privchip}"
                    [ -d ${privchip} ] && echo "and there is a directory called '${privchip}'."
                    [ -d ${privchip} ] || echo "${privchip} does not exist."

                    for i in ${samples}; do
                            directory=${privchip}${i}
                            # check if the samples exist in ${privchip}.
                            echo "$directory"
                            [ -d $directory ] && echo "and there is a directory with the name ${i} in ${privchip}"
                            [ -d $directory ] || echo "${i} does not exist in ${privchip}. Start over and check your  path and samples."
                            [ -d $directory ] || continue

                    done

            done
echo ""
#Enter if the data is stranded or unstranded and test if what is entered is a valid option.

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

## This part will let you type in the distance between the mate pairs since this is important in Bowtie-mapping of paired-end sequencing data.

        while [ "${correctlength}" != yes ] || [ "${correctlength}" != no ]; do
            if [ ${seqtype} = paired ]; then
                echo "How far is it between your read-pairs (i.e) sequencing length+1? If you are unsure, type 250!"
                    read seqlength
                echo "Your read-length have been site to ${seqlength}. Is this correct? 'yes' or 'no'. You must type something or the script will fail."
                    read correctlength
                        if [ ${correctlength} = yes ] || [ ${seqtype} = no ]; then break
                        else
                            echo "!!!!${correctlength} is an invalid option!!!!"
                        fi
              else
                break
            fi
        done

#Choose reference genome, mm9 or mm10. Currently only full support for mm9
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
            BOWTIE_REF_DIR=${bowtie_refdir_mm9}
            REPEAT_MASKER_FILE=${repeat_masker_mm9}
            BLACK_LIST_FILE=${black_list_mm9}
        ;;
        #OPTION: paired end data data
        mm10)
            REF_DIR_SIZE=${refdir_size_mm10}
            BOWTIE_REF_DIR=${bowtie_refdir_mm10}
            REPEAT_MASKER_FILE=${repeat_masker_mm10}
            BLACK_LIST_FILE=${black_list_mm10}
        ;;
        esac

        echo "Your reference files are in ${REF_DIR_SIZE} (genome sizes), and ${BOWTIE_REF_DIR} (Bowtie reference directory)."

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


##############Ask if duplicates should be removed
echo ""
while [ "${rmdup}" != yes ] || [ "${rmdup}" != no ]; do
echo "Do you want to remove duplicate reads after mapping? 'yes' or 'no'. You must type something or the script will fail."
    read rmdup
        if [ ${rmdup} = yes ]  || [ ${rmdup} = no ]; then break
        else
                echo "!!!!${rmdup} is an invalid option!!!!"
        fi
done

    if [ ${rmdup} == yes ]; then
        echo "Samtools will be used to remove duplicates from your data."
    else
        echo "You have chosen not to remove duplicate reads."
    fi

##############Ask if masked and blacklisted regions should be removed
echo ""
while [ "${rmsk}" != yes ] || [ "${rmsk}" != no ]; do
echo "Do you want to subtract masked and blacklisted regions? Only choose if you have problems with unspecific binding to these sites. Enter 'yes' or 'no'. You must type something or the script will fail."
    read rmsk
        if [ ${rmsk} = yes ]  || [ ${rmsk} = no ]; then break
            else
                echo "!!!!${rmsk} is an invalid option!!!!"
        fi
done

    if [ ${rmsk} == yes ]; then
        echo "Bedtools will be used to remove regions defined by the UCSC '${refgenome}' RepeatMasker file as well as the ENCODE blacklist regions file from your data."
    else
        echo "You have chosen not subtract RepeatMasker regions."
    fi


### Setting the variables based on the selections.


dataType="${strandedness}|${seqtype}"

    #The combinations for dataType that exist are the following:
        #unstranded|single
        #unstranded|paired
        #stranded|single
        #stranded|paired



echo "######################################################################################################"
echo ""
echo "You will begin to process the '${samples}' samples in ${privchip} and your data is '${strandedness}' and '${seqtype}'-end sequenced. You have selected ${refgenome} as your reference genome".

    if [ ${trimstatus} == yes ]; then
        echo "Your data will be trimmed with Trimmomatic."
    else
        echo "Your data will not be trimmed with Trimmomatic."
    fi

rmdup_rmsk="${rmdup}|${rmsk}"

case ${rmdup_rmsk} in
    "yes|yes")
        echo "You will also remove duplicate reads with Samtools and masked regions defined by '${refgenome}' RepeatMasker file as well as the ENCODE blacklisted regions."
    ;;
    "yes|no")
        echo "You will also remove duplicate reads with Samtools."
    ;;
    "no|yes")
        echo "You will also remove masked regions defined by '${refgenome}' RepeatMasker file as well as the ENCODE blacklisted regions."
    ;;
esac


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
        directory=${privchip}${j}

        # go to the folder that contains fastq files
        cd $directory

       if [ ! -d ${j}${suffix} ]; then

        echo "Processing ${j}!"

        mkdir ${j}.bowtie
        echo "Unzipping and catenating ${j}.fastq.gz-files."


# Unzip the data (will depend on if the data is single and paired-end.

        cd ${j}.bowtie

        mkdir ${j}.processed_files
        mkdir ${j}.for_genome_browser

_mydir="`pwd`"

echo "You are now in: $_mydir"

 case "${seqtype}" in
       single)
            #OPTION: single-end data
            gunzip -c ${privchip}${j}/*.fastq.gz > ${privchip}${j}/${j}.bowtie/${j}.fastq
            Trimmode=SE
            Trimadapter=/rothenberg/sulu/data00/Trimmomatic/adapters/TruSeq3-SE.fa
       ;;
            #OPTION: paired-end data data
       paired)
            gunzip -c ${privchip}${j}/*R1*.fastq.gz > ${privchip}${j}/${j}.bowtie/${j}.R1.fastq
            gunzip -c ${privchip}${j}/*R2*.fastq.gz > ${privchip}${j}/${j}.bowtie/${j}.R2.fastq
            Trimmode=PE
            Trimadapter=/rothenberg/sulu/data00/Trimmomatic/adapters/TruSeq3-PE.fa
       ;;
       esac

        ## FastQC check of the untrimmed file.

        case "${seqtype}" in
        single)
            #OPTION: single-end data
            /rothenberg/sulu/data00/FastQC/./fastqc -t ${threads} ${j}.fastq
            TRIMfastqin=${j}.fastq
            TRIMfastqout=${j}.trimmed.fastq
            bowtiesettings="-v 3 -k 11 -m 10 -t --best --strata"
            unmapped_out="${j}.${refgenome}.unmapped.fq"
            rmdup_seqtype="-s"
        ;;
        paired)
            #OPTION: paired-end data data
            /rothenberg/sulu/data00/FastQC/./fastqc -t ${threads} ${j}.R1.fastq
            /rothenberg/sulu/data00/FastQC/./fastqc -t ${threads} ${j}.R2.fastq
            TRIMfastqin="${j}.R1.fastq ${j}.R2.fastq"
            TRIMfastqout="${j}.R1.paired.trimmed.fastq ${j}.R1.unpaired.trimmed.fastq ${j}.R2.paired.trimmed.fastq ${j}.R2.unpaired.trimmed.fastq"
            # These values are taken from Biostars Galaxy toolshed.
            bowtiesettings="-X ${seqlength} --fr --chunkmbs 1024 -t"
            unmapped_out="${j}.${refgenome}.unmapped"
            rmdup_seqtype=""

        ;;
        esac


        ###### Trimmomatic trimming of the data

        if [ ${trimstatus} == yes ]; then

            ## Trimming of the fastq-file. Values taken from http://www.usadellab.org/cms/?page=trimmomatic.
            echo "Trimmomatic is processing ${j}.fastq."
            echo " java -jar /rothenberg/sulu/data00/Trimmomatic/trimmomatic-0.33.jar ${Trimmode} -threads ${threads} -trimlog ${j}.trimlog.txt ${TRIMfastqin} ${TRIMfastqout} ILLUMINACLIP:${Trimadapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25"

            java -jar /rothenberg/sulu/data00/Trimmomatic/trimmomatic-0.33.jar ${Trimmode} -threads ${threads} -trimlog ${j}.trimlog.txt ${TRIMfastqin} ${TRIMfastqout} ILLUMINACLIP:${Trimadapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25

        else
            echo "Moving on without trimming"

            if [ ${seqtype} = single ]; then
                bowtiefastq=${TRIMfastqin}
            elif [ ${seqtype} = paired ]; then
                bowtiefastq="-1 ${j}.R1.fastq -2 ${j}.R2.fastq"
            else
                echo "What can of data type did you choose up there. Something's upt with the script!!!"
            fi
        fi


            # FastQC on the trimmed data if trimming has been chosen

        if [ ${trimstatus} == yes ]; then
            case "${seqtype}" in
            single)
                #OPTION: single-end data
                /rothenberg/sulu/data00/FastQC/./fastqc -t ${threads} ${TRIMfastqout}
                bowtiefastq=${TRIMfastqout}
            ;;
            paired)
                #OPTION: paired end data data
                /rothenberg/sulu/data00/FastQC/./fastqc -t ${threads} ${j}.R1.paired.trimmed.fastq
                /rothenberg/sulu/data00/FastQC/./fastqc -t ${threads} ${j}.R2.paired.trimmed.fastq
                /rothenberg/sulu/data00/FastQC/./fastqc -t ${threads} ${j}.R1.unpaired.trimmed.fastq
                /rothenberg/sulu/data00/FastQC/./fastqc -t ${threads} ${j}.R2.unpaired.trimmed.fastq
                bowtiefastq="-1 ${j}.R1.paired.trimmed.fastq -2 ${j}.R2.paired.trimmed.fastq"
            ;;
            esac
        else
            echo "No trimmed file to run FastQC on!"
        fi

    echo ""

    echo "Starting to map ${bowtiefastq} with Bowtie1."

    echo "/rothenberg/sulu/data00/bowtie/bowtie /rothenberg/sulu/data00/bowtie-indexes/${refgenome} -p ${threads} ${bowtiesettings} -q ${bowtiefastq} --un ${unmapped_out} --max ${j}.${refgenome}.repeat.fq -S ${j}.${refgenome}.sam"

  /rothenberg/sulu/data00/bowtie/bowtie /rothenberg/sulu/data00/bowtie-indexes/${refgenome} -p ${threads} ${bowtiesettings} -q ${bowtiefastq} --un ${unmapped_out} --max ${j}.${refgenome}.repeat.fq -S ${j}.${refgenome}.sam


#######################Remove duplicates and masked regions


        if [ ${rmdup} == yes ] || [ ${rmsk} == yes ]; then

            case ${rmdup_rmsk} in
            # duplicates, blacklisted and masked regions removed
            "yes|yes")


                suffix="_rmdup_rmsk_${refgenome}_repeat_addback_tagdir"
                bedname="rmdup.rmsk.repeat_addback"

                #SamToBam
                echo "Creating ${j}.bam"
                samtools view -bS -o ${j}.bam ${j}.${refgenome}.sam
                #Sort the bamfile
                samtools sort ${j}.bam ${j}.s
                #remove duplicate reads from the sorted bam-file
                echo "Removing duplicate reads from ${j}.s.bam"
                samtools rmdup ${rmdup_seqtype} ${j}.s.bam ${j}.rmdup.bam
                #Create the bed-file
                echo "bamToBed -i ${j}.rmdup.bam > ${j}.rmdup.bed"
                bamToBed -i ${j}.rmdup.bam > ${j}.rmdup.bed
                # sort bedfile
                echo "bedSort ${j}.rmdup.bed ${j}.rmdup.s.bed"
                bedSort ${j}.rmdup.bed ${j}.rmdup.s.bed

                #Subtract the bed information in the UCSC ReperatMask bed as well as the ENCODE blacklist regions from the rmdup bed-file.
                echo "Subtracting ${refgenome} masked and blacklisted regions from ${j}.rmdup.s.bed."

                echo "bedtools subtract -A -a ${j}.rmdup.s.bed -b ${BLACK_LIST_FILE} > ${j}.${bedname}.blacklist.bed"
                bedtools subtract -A -a ${j}.rmdup.s.bed -b ${BLACK_LIST_FILE} > ${j}.${bedname}.blacklist.bed
                echo "bedtools subtract -A -a ${j}.${bedname}.blacklist.bed -b ${REPEAT_MASKER_FILE} > ${j}.${bedname}.bed"
                bedtools subtract -A -a ${j}.${bedname}.blacklist.bed -b ${REPEAT_MASKER_FILE} > ${j}.${bedname}.bed


                #This part will add one "read" back over the repeat sequences. It will not add it back over the blasklisted regions since this messes up downstream processing.

                echo "This part will add one "read" back over the repeat sequences. It will not add it back over the blasklisted regions since this messes up downstream processing."
                echo "cat ${j}.${bedname}.bed ${REPEAT_MASKER_FILE} > ${j}.${bedname}.black.rmsk.bed"
                cat ${j}.${bedname}.bed ${REPEAT_MASKER_FILE} > ${j}.${bedname}.black.rmsk.bed
                bedSort ${j}.${bedname}.black.rmsk.bed ${j}.${bedname}.bed

                rm ${j}.${bedname}.black.rmsk.bed
                rm ${j}.${bedname}.blacklist.bed

            ;;
            #duplicates removed
            "yes|no")
                suffix="_rmdup_${refgenome}_tagdir"
                bedname="rmdup"

                #SamToBam
                echo "Creating ${j}.bam"
                samtools view -bS -o ${j}.bam ${j}.${refgenome}.sam
                #Sort the bamfile
                samtools sort ${j}.bam ${j}.s
                #remove duplicate reads from the sorted bam-file
                echo "Removing duplicate reads from ${j}.s.bam"
                samtools rmdup ${rmdup_seqtype} ${j}.s.bam ${j}.rmdup.bam
                #Create the bed-file
                echo "bamToBed -i ${j}.rmdup.bam > ${j}.rmdup.bed"
                bamToBed -i ${j}.rmdup.bam > ${j}.rmdup.bed
                # sort bedfile
                echo "bedSort ${j}.rmdup.bed ${j}.rmdup.s.bed"
                bedSort ${j}.rmdup.bed ${j}.${bedname}.bed

            ;;
            #masked and blacklisted regions removed
            "no|yes")
                suffix="_rmsk_${refgenome}_repeat_addback_tagdir"
                bedname="rmsk.repeat_addback"

                #SamToBam
                echo "Creating ${j}.bam"
                samtools view -bS -o ${j}.bam ${j}.${refgenome}.sam
                #Sort the bamfile
                samtools sort ${j}.bam ${j}.s
                #Create the bed-file
                echo "bamToBed -i ${j}.s.bam > ${j}.rmdup.bed"
                bamToBed -i ${j}.s.bam > ${j}.bed
                # sort bedfile
                echo "bedSort ${j}.rmdup.bed ${j}.rmdup.s.bed"
                bedSort ${j}.bed ${j}.s.bed

                #Subtract the bed information in the UCSC ReperatMask bed and blacklist files from the rmdup bed-file.
                echo "Subtracting ${refgenome} masked and blacklist regions from ${j}.rmdup.s.bed."

                echo "bedtools subtract -A -a ${j}.s.bed -b ${BLACK_LIST_FILE} > ${j}.${bedname}.blacklist.bed"
                bedtools subtract -A -a ${j}.s.bed  -b ${BLACK_LIST_FILE} > ${j}.${bedname}.blacklist.bed
                echo "bedtools subtract -A -a ${j}.rmdup.s.bed -b ${REPEAT_MASKER_FILE} > ${j}.${bedname}.bed"
                bedtools subtract -A -a ${j}.${bedname}.blacklist.bed -b ${REPEAT_MASKER_FILE} > ${j}.${bedname}.bed

                #This part will add one "read" back over the repeat sequences. It will not add it back over the blasklisted regions since this messes up downstream processing.

                echo "This part will add one "read" back over the repeat sequences. It will not add it back over the blasklisted regions since this messes up downstream processing."
                echo "cat ${j}.${bedname}.bed ${REPEAT_MASKER_FILE} > ${j}.${bedname}.black.rmsk.bed"
                cat ${j}.${bedname}.bed ${REPEAT_MASKER_FILE} > ${j}.${bedname}.black.rmsk.bed
                bedSort ${j}.${bedname}.black.rmsk.bed ${j}.${bedname}.bed

                rm ${j}.${bedname}.black.rmsk.bed
                rm ${j}.${bedname}.blacklist.bed


            ;;
            esac
        else
            suffix="_${refgenome}_tagdir"
            bedname=""

        fi


        # create TagDirectory and keeping  the reads from the genome generated bam file. The input here will depend on how the data has been processed above

        if  [ ${rmdup_rmsk} = "no|no" ]; then

            echo ""
            echo "Creating a Homer tag directory reads from ${j}.${refgenome}.sam"

            echo "makeTagDirectory ${j}${suffix} ${j}.${refgenome}.sam -format sam"

            /rothenberg/sulu/data00/homer/bin/makeTagDirectory ${j}${suffix} ${j}.${refgenome}.sam -format sam

            #### Bigwig generation

            # convert bam to bed
            echo "Creating and sorting a bed-file from ${j}.${refgenome}.sam."
            # convert sam to bam
            samtools view -bS -o ${j}.bam ${j}.${refgenome}.sam
            # convert bam to bed
            bamToBed -i ${j}.bam > ${j}.bed
            # sort bed-file
            bedSort ${j}.bed ${j}.${refgenome}.s.bed
            # get normaliztion factor for bigwig file
            scalevar=$( cat ${j}.${refgenome}.s.bed | wc -l )
            scalevar=$( echo "scale=4; 1000000/$scalevar" | bc )
            echo "**Inverse of Number of Reads (in Millions)**  $scalevar"
            # compile reads
            echo "Generating ${j}.${refgenome}.bg"
            genomeCoverageBed -bg -i ${j}.${refgenome}.s.bed -g ${REF_DIR_SIZE} -scale $scalevar > ${j}.${refgenome}.bg
            #convert bedgraph to bigwig file
            echo "Generating ${j}.${refgenome}.bw"
            bedGraphToBigWig ${j}.${refgenome}.bg ${REF_DIR_SIZE} ${j}.${refgenome}.bw
            mv ${j}.${refgenome}.bg ${j}.${refgenome}.s.bed ${j}.for_genome_browser
            scp ${j}.${refgenome}.bw ${j}.for_genome_browser ## The original file is used for a later check and then removed.
            ln -s  "$directory/${j}.bowtie/${j}.for_genome_browser/${j}.${refgenome}.bw" ${bwtrack_dir}

        else
            echo ""
            echo "Creating a Homer tag directory reads from ${j}.${bedname}.bed"

            echo "makeTagDirectory ${j}${suffix} ${j}.${bedname}.bed -format bed -forceBED"

            /rothenberg/sulu/data00/homer/bin/makeTagDirectory ${j}${suffix} ${j}.${bedname}.bed -forceBED

            # sort bed-file
            bedSort ${j}.${bedname}.bed ${j}.${refgenome}.s.bed
            # get normaliztion factor for bigwig file
            scalevar=$( cat ${j}.${refgenome}.s.bed | wc -l )
            scalevar=$( echo "scale=4; 1000000/$scalevar" | bc )
            echo "**Inverse of Number of Reads (in Millions)**  $scalevar"
            # compile reads
            echo "Generating ${j}.${refgenome}.bg"
            genomeCoverageBed -bg -i ${j}.${refgenome}.s.bed -g ${REF_DIR_SIZE} -scale $scalevar > ${j}.${refgenome}.bg
            #convert bedgraph to bigwig file
            echo "Generating ${j}.${refgenome}.bw"
            bedGraphToBigWig ${j}.${refgenome}.bg ${REF_DIR_SIZE} ${j}.${refgenome}.bw

            # Some renaming and moving of files for future understanding
            scp ${j}.${refgenome}.bw ${j}.${bedname}.${refgenome}.bw ## The original file is used for a later check and then removed.
            echo "Moving ${j}.${bedname}.${refgenome}.bw to ${j}.for_genome_browser"
            mv ${j}.${bedname}.${refgenome}.bw ${j}.for_genome_browser
            ln -s  "$directory/${j}.bowtie/${j}.for_genome_browser/${j}.${bedname}.${refgenome}.bw" ${bwtrack_dir}
            mv ${j}.${refgenome}.bg ${j}.${bedname}.${refgenome}.bg
            mv ${j}.${refgenome}.s.bed ${j}.${bedname}.${refgenome}.s.bed
            mv ${j}.${bedname}.${refgenome}.s.bed ${j}.${bedname}.${refgenome}.bg ${j}.for_genome_browser

        fi


        ## Cleaning up and moving files
        echo "Starting to clean up!"
        rm *bam
        rm *sam
        mv ${j}.${refgenome}.repeat.fq ${j}.${refgenome}.unmapped* ${j}.processed_files/
        rm *fastq
        rm *.bed

echo ""
####Testing if the size of the bigwig is correct (not to small). If it is the script will assume something is wrong and exit. Most often an erronous file have the value of 384 and then you know that something is wrong and potentially could carry over to the next round

        bigwig_file_size=$(cat ${j}.${refgenome}.bw | wc -c)
        echo "The size of the bigwig is: ${bigwig_file_size}."
        minimum_size=100000
        if [ ! -f ${j}.${refgenome}.bw ] || [ "${bigwig_file_size}" -le "${minimum_size}" ]; then
            echo "${j}.${refgenome}.bw is to small. Something is wrong with the script and the rest of your samples will not be processed." && exit
        else
            echo "${j}.${refgenome}.bw seems ok. The script will continue."
        fi

    rm ${j}.${refgenome}.bw

    else
      echo "${j}${suffix} is already there. skipping..."
    fi
done
