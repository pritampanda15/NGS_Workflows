
#!/bin/sh
####################################################################################################################################################################################################
#################		Panda-Exome-seq variant calling pieline accoring to GATK Best practices. (Somatic and Germline)							####################
####################################################################################################################################################################################################


####################################################################################################################################################################################################
# The script is made for a germline automated workflow that will call SNPs and INDELs
# In order to make this somatic one can simply stop post Step 2 and use Mutect2 or VarScan2 with their own custom scripts and filters on the processed bam files.
# For somatic Pipeline we use Mutect2 which you can find in the same pipeline.
# Call with following arguments
# sh master_script.sh  <output_basename> <fastq folder> <output_folder_loc> [cpus]
# you can run the script in above mentioned way through another processing script that will log the processing with time at each step
# Putting the bwa mem
# GATK bundle set : one can obtain these from gatk ftp (knonw as gatk bundle)
# ftp://ftp.broadinstitute.org/bundle/hg19/
# The bundles were downloaded from the GATK 2.8 while the GATK 4.0.2 & 3.8.0 was used for the analysis
# Please download all the required tools and install them in your Cluster Computing and modify the paths accordingly for effective running of the script						
#####################################################################################################################################################################################################




####################################################################################################################################################################################################
################## 			Initializing the variables required for calling the automated script									####################
####################################################################################################################################################################################################

bn=$1
floc=$2
outloc=$3
nproc=${4:-12}

####################################################################################################################################################################################################
################## 			Path to reference genome and Index files (This code is based on hg19)
##################				mention the location of the indexed hg19 file											####################
####################################################################################################################################################################################################

bwa_idx="/media/panda/hdd/Reference/ucsc.hg19.fasta"
#mention the location of the hg19 fasta file that will be used as the reference genome
ref="/media/panda/hdd/Reference/hg19/ucsc.hg19.fasta"
#Tools and environment required to run this script
#Path to STAR,java,gatk and picard tools
java_home="/home/panda/miniconda3/bin/java"
#location for bwa mem
bwamem="/home/panda/anaconda2/bin/bwa"
#location for the GATK toolkit jar file
gatk="/media/panda/hdd/Reference/GenomeAnalysisTK.jar"
#location of latest gatktool
gatk4="/media/panda/hdd/Reference/gatk-package-4.0.2.0-local.jar"
#location for picard for marking duplicates
picard="/media/panda/hdd/Reference/picard.jar"
#location for sambamba which is used for processing the bam files. can be found here http://lomereiter.github.io/sambamba/
sambamba="/home/panda/miniconda3/bin/sambamba"
#location for Annotation with SNPEff
snpEff="/media/panda/hdd/Reference/snpEff_latest_core/snpEff/snpEff.jar"
#location of hg19folder for annotation of SNPEff
hg19="/media/panda/hdd/Reference/snpEff_latest_core/snpEff/data/hg19"
vardict="/home/panda/bin/vardict"


####################################################################################################################################################################################################
################## 			Path to gatk bundle set files ( these files are used for the different base quality calibration steps that are useful for the variant
					#workflow. These bundles were downloaded from GATK resource bundle									####################
####################################################################################################################################################################################################

millsIndels="/media/panda/hdd/Reference/Mills_and_1000G_gold_standard.indels.hg19.vcf"
KGIndels="/media/panda/hdd/Reference/1000G_phase1.indels.hg19.vcf"
dbSNP138="/media/panda/hdd/Reference/dbsnp_138.hg19.vcf"
gnomad="/media/panda/hdd/Reference/af-only-gnomad.raw.sites.b37.vcf.gz"
#Exon target baits from Agilent SSv4 set (should be upgraded to new version if new agilent kit is used)
ss4exonbaits="/media/panda/hdd/Reference/Exon_SSV4_clean.list"
target="/media/panda/hdd/Reference/Exon_SSV4_clean.bed"

####################################################################################################################################################################################################
################## 			  Create an output directory														####################
####################################################################################################################################################################################################

mkdir -p "$outloc/${bn}_processed"
#fout="$outloc/${bn}_processed"
echo "output directory for fastq $outloc/${bn}"_processed" ..."

#fout=$outloc/$bn"_processed"
fout="$outloc/${bn}_processed"
echo "$fout ..."

####################################################################################################################################################################################################
################## 			  Performing assembly to create one fastq file for each read mates									####################
####################################################################################################################################################################################################


echo "performing assembly to create one fastq file for each read mates ..."
zcat $floc/*1.fastq.gz > $fout/${bn}_R1.fastq
zcat $floc/*2.fastq.gz > $fout/${bn}_R2.fastq

####################################################################################################################################################################################################
################## 			  Aligning the paired end read mates to the reference genome hg19 with bwa mem. This is for paired-end data				####################
####################################################################################################################################################################################################
 
echo -e "["$(date)"]\tAligning.."
$bwamem mem -t 8 -M -R "@RG\tID:${bn}\tLB:PairedEnd\tPL:Illumina\tPU:000000000-A442D\tSM:${bn}" $ref $fout/${bn}_R1.fastq $fout/${bn}_R2.fastq | samtools view -bS -F 0x04 - > $fout/${bn}"_Aligned.out.bam"


####################################################################################################################################################################################################
################## 			  Converting SAM to BAM	using PICARD: Its an optional command to convert the sam to bam files.						####################
####################################################################################################################################################################################################
#echo -e "["$(date)"]\tConverting SAM to BAM.."
#$java_home -Xmx20g -jar $picard SortSam INPUT=$fout/${bn}"_Aligned.sam" OUTPUT=$fout/${bn}"_Aligned.out.bam" SORT_ORDER=coordinate


####################################################################################################################################################################################################
################## 			  Sorting the bam files with sambamba													####################
####################################################################################################################################################################################################

echo -e "["$(date)"]\tSorting.."
$sambamba sort -o $fout/${bn}"_sorted.bam" -p -t 12 $fout/${bn}"_Aligned.out.bam"


####################################################################################################################################################################################################
################## 			  Remove the Aligned bam for memory management(Optional)										####################
####################################################################################################################################################################################################

#rm $fout/${bn}"_Aligned.out.bam"


####################################################################################################################################################################################################
################## 			  indexing the sorted aligned bam file													####################
####################################################################################################################################################################################################

echo -e "["$(date)"]\tIndexing.."
$sambamba index -p -t 12 $fout/${bn}"_sorted.bam"

####################################################################################################################################################################################################
################## 			Picard mark duplicates
#					Remove the reads due to PCR duplicates or one can mark them with a flag									####################
####################################################################################################################################################################################################


echo -e "["$(date)"]\tMarking duplicates.."
$java_home -Xmx32g -XX:-UseGCOverheadLimit -jar $picard MarkDuplicates I=$fout/${bn}"_sorted.bam" O=$fout/${bn}"_dupMarked.bam" M=$fout/${bn}"_dup.metrics" CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT 2>$fout/${bn}.MarkDuplicates.log

####################################################################################################################################################################################################
################## 			  Remove the Aligned bam,Sorted.bam for memory management(Optional)									####################
####################################################################################################################################################################################################

#rm $fout/${bn}"_Aligned.out.bam"
#rm $fout/${bn}"_sorted.bam"
#rm $fout/${bn}"_sorted.bam.bai"

####################################################################################################################################################################################################
################## 			  FOR RNASeq, SplitNCigarReads														####################
####################################################################################################################################################################################################

#echo -e "["$(date)"]\tSpliting reads.."
#$java_home -Xmx2g -jar $gatk -T SplitNCigarReads -R $ref -I $fout/${bn}"_dupMarked.bam" -o $fout/${bn}"_split.bam" -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS 2>$fout/${bn}.SplitNCigarReads.log


####################################################################################################################################################################################################
################## 			  #Now the GATK 3.5 process handles starts from below
					  #Step 1 Local realignment around the INDELS (please refer to http://gatkforums.broadinstitute.org/gatk/discussion/38/local-realignment-around-indels#latest)
					  #Create targets for indel realignment Step 1												####################
####################################################################################################################################################################################################

echo -e "["$(date)"]\tCreating targets for indel realignment.."
$java_home -Xmx20g -jar $gatk -T RealignerTargetCreator -R $ref -I $fout/${bn}"_dupMarked.bam" -o $fout/${bn}".intervals" -nt 12 -known $millsIndels -known $KGIndels -L $ss4exonbaits 2>$fout/${bn}.indel.log

####################################################################################################################################################################################################
################## 			  
					  #Perform indel realignment														####################
####################################################################################################################################################################################################


echo -e "["$(date)"]\tPerforming Indel Realignment.."
$java_home -Xmx20g -jar $gatk -T IndelRealigner -R $ref -I $fout/${bn}"_dupMarked.bam" -targetIntervals $fout/${bn}".intervals" -known $millsIndels -known $KGIndels -o $fout/${bn}"_realigned.bam" 2>$fout/${bn}.indel2.log


####################################################################################################################################################################################################
################## 			  Remove the Dedupmarked for memory management(Optional)										####################
####################################################################################################################################################################################################
#rm $fout/${bn}"_split.bam"
#rm $fout/${bn}"_split.bai"
#rm $fout/${bn}"_dupMarked.bam"
#rm $fout/${bn}"_dupMarked.bai"


####################################################################################################################################################################################################
################## 			#Step 2: Perform BQSR (base quality recalibration after realigment step. This is a standard test that is performed to Detect 
					#systematic errors in base quality scores.		  										####################
####################################################################################################################################################################################################


echo -e "["$(date)"]\tPerforming BQSR.."
$java_home -Xmx2g -jar $gatk -T BaseRecalibrator -I $fout/${bn}"_realigned.bam" -R $ref -knownSites $KGIndels -knownSites $millsIndels -knownSites $dbSNP138 -o $fout/${bn}"_recal.table" -L $ss4exonbaits 2>$fout/${bn}.BQSR.log

####################################################################################################################################################################################################
################## 			#Print recalibrated reads. This handle writes out sequence read data (for filtering, merging, subsetting etc)
						  																####################
####################################################################################################################################################################################################


echo -e "["$(date)"]\tPrinting recalibrated reads.."
$java_home -Xmx20g -jar $gatk -T PrintReads -R $ref -I $fout/${bn}"_realigned.bam" -nct 12 -BQSR $fout/${bn}"_recal.table" -o $fout/${bn}"_recal.bam" 2>$fout/${bn}.BQSR2.log

####################################################################################################################################################################################################
################## 			Step3:Call the variants post realignment and recalibration with Haplotype caller. It calls germline SNPs and indels via local re-assembly of haplotypes
						  																####################
####################################################################################################################################################################################################


echo -e "["$(date)"]\tRunning HaplotypeCaller.."
#$java_home -Xmx20g -jar $gatk -T HaplotypeCaller -R $ref -I $fout/${bn}"_recal.bam" -L $ss4exonbaits -dontUseSoftClippedBases -o $fout/${bn}".vcf" 2>$fout/${bn}.HaplotypeCaller.log
#$java_home -Xmx2g -jar $gatk -T HaplotypeCaller -R $ref -I $fout/${bn}"_recal.bam" -L $ss4exonbaits -dontUseSoftClippedBases -o $fout/${bn}".vcf" -gt_mode DISCOVERY --dbsnp $dbSNP138 -G Standard -G AS_Standard 2>$fout/${bn}.HaplotypeCaller.log

####################################################################################################################################################################################################
################## 			#Step 4: Call the Somatic pipeline using Mutect2 with latest GATK4.0 tools (Optional if you call somatic variants)
						  																####################
####################################################################################################################################################################################################

echo -e "["$(date)"]\tRunning Mutect2.."
#$java_home -Xmx20g -jar $gatk4 Mutect2 -R $ref -I $fout/${bn}"_recal.bam" -tumor ${bn} -L $ss4exonbaits -O $fout/${bn}"_Somatic.vcf"  2>$fout/${bn}.Mutect2.log
# starting from fastq files to vcf we have to specify the name of the $bn or the file name that you specified during running the command
#$java_home -Xmx20g -jar $gatk4 Mutect2 -R $ref -I $fout/${bn}"_recal.bam" -tumor ${bn} -L $ss4exonbaits -O $fout/${bn}"_Somatic.vcf" 2>$fout/${bn}.Mutect2.log


####################################################################################################################################################################################################
################## 			# Seperate out the SNPs from the from INDELs followed by filtration
					#Filter variants for SNPs and then filter for qaulity scores. Variant filtration methods
						  																####################
####################################################################################################################################################################################################


echo -e "["$(date)"]\tFiltering Point Variants.."
#$java_home -Xmx20g -jar $gatk -T SelectVariants -R $ref -V $fout/${bn}".vcf" -selectType SNP -o $fout/${bn}"_raw_snps.vcf"

####################################################################################################################################################################################################
################## 			# Filter over snps for scores that includes the combination of different Variant filter expression
					# Filter variant calls based on INFO and FORMAT annotations
						  																####################
####################################################################################################################################################################################################

#$java_home -Xmx20g -jar $gatk -T VariantFiltration -R $ref -V $fout/${bn}"_raw_snps.vcf" --filterExpression "DP < 5 || SB > -0.1 || QUAL < 50.0 || QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || DP > 500" --filterName "PANDA_snp_filters" -o $fout/${bn}"_filtered_snps.vcf" 2>$fout/${bn}.VariantFilter_snps.log

####################################################################################################################################################################################################
################## 			# Filtering INDEL Variants
					# Filter variant calls based on INFO and FORMAT annotations
						  																####################
####################################################################################################################################################################################################

echo -e "["$(date)"]\tFiltering INDEL Variants.."
###select indels with the SelectVariants handle and -selectType process to give raw indels
#$java_home -Xmx2g -jar $gatk -T SelectVariants -R $ref -V $fout/${bn}".vcf" -selectType INDEL -o $fout/${bn}"_raw_indels.vcf"

####################################################################################################################################################################################################
################## 			# Filter over raw indel indel files based on ased on INFO and FORMAT annotations
					# Filter variant calls based on INFO and FORMAT annotations
						  																####################
####################################################################################################################################################################################################

#$java_home -Xmx20g -jar $gatk -T VariantFiltration -R $ref -V $fout/${bn}"_raw_indels.vcf" --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" --filterName "DAS_indel_filters" -o $fout/${bn}"_filtered_indels.vcf" 2>$fout/${bn}.VariantFilter_indels.log

####################################################################################################################################################################################################
################## 			# Vardict for somatic clones below VAF<5% 												####################
####################################################################################################################################################################################################


$vardict -G $ref -b $fout/${bn}"_recal.bam" -z 1 -k 1 -c 1 -S 2 -E 3 -g 4 -f 0 $target > $fout/${bn}"_somatic_clones.txt"

#to run in a loop:
#for i in *.bam; do vardict -G '/media/panda/hdd/Reference/ucsc.hg19.fasta' -b $i -z 1 -k 1 -c 1 -S 2 -E 3 -g 4 28_panel.bed > $i.txt; done
echo -e "["$(date)"]\tDONE!"
