#!/bin/bash

# this bash script is to integrate the qubs and aux. perl scripts together 
# to streamlize the pipeline. 
# it will take the following form: 
#
# RNAseq_pipeline_{projname}.sh [steps_to_run] 
#
# where 
#   steps_to_run, is optional, if given, will only run steps passed in this param.
#      otherwise, will run from step 1 to 9
#      format shall be "4 5 6 7 8 9", "1 3", etc, separate by space
# for example, 
#     RNAseq_pipeline_APOE_Mye_Micro_trim_100.sh "4 5 6 7 8 9"

PROJNAME=Astrocyte_Ab
if [ $# -eq 1 ] 
then 
   STEPS_TO_RUN=($1)
else
   STEPS_TO_RUN=(1 2 3 4 5 6 7 8 9)
fi
PREV_JOB=""
CUR_JOB=""
# generate fastqc script: 
for i in ${STEPS_TO_RUN[@]}
do
    if [ $i -eq 1 ] 
    then
	CUR_JOB=RNASeq_${PROJNAME}_fastq_link
	CMD="qsub -N $CUR_JOB RNAseq_fastq_link_${PROJNAME}.qsub"
	echo "STEP$i: $CMD"
	eval $CMD
	PREV_JOB=$CUR_JOB
    fi
    if [ $i -eq 2 ] 
    then 
	CUR_JOB=aRNASeq_${PROJNAME}_fastqc
	if [ -z $PREV_JOB ]
	then
	    CMD="qsub -N $CUR_JOB RNAseq_fastqc_${PROJNAME}.qsub"
        else
	    CMD="qsub -hold_jid $PREV_JOB -N $CUR_JOB RNAseq_fastqc_${PROJNAME}.qsub"
        fi
	echo "STEP$i: $CMD"
	eval $CMD
	PREV_JOB=$CUR_JOB
    fi
    if [ $i -eq 3 ] 
    then 
	CUR_JOB=RNASeq_${PROJNAME}_multiqc
	if [ -z $PREV_JOB ]
	then
	    CMD="qsub -N $CUR_JOB  RNAseq_multiqc_${PROJNAME}.qsub"
        else
	    CMD="qsub -hold_jid $PREV_JOB -N $CUR_JOB RNAseq_multiqc_${PROJNAME}.qsub"
        fi
	echo "STEP$i: $CMD"
	eval $CMD
	PREV_JOB=$CUR_JOB
    fi
    if [ $i -eq 4 ] 
    then 
	CUR_JOB=aRNASeq_${PROJNAME}_align
	if [ -z $PREV_JOB ]
	then
	    CMD="qsub -N $CUR_JOB RNAseq_aln_${PROJNAME}.qsub"
        else
	    CMD="qsub -hold_jid $PREV_JOB -N $CUR_JOB RNAseq_aln_${PROJNAME}.qsub"
        fi

	echo "STEP$i: $CMD"
	eval $CMD
	PREV_JOB=$CUR_JOB
    fi
    if [ $i -eq 5 ] 
    then
	CUR_JOB=RNASeq_${PROJNAME}_mqc_aln
	if [ -z $PREV_JOB ]
	then
	    CMD="qsub -N $CUR_JOB RNAseq_multiqc_aln_${PROJNAME}.qsub"
        else
	    CMD="qsub -hold_jid $PREV_JOB -N $CUR_JOB RNAseq_multiqc_aln_${PROJNAME}.qsub"
        fi

	echo "STEP$i: $CMD"
	eval $CMD
	PREV_JOB=$CUR_JOB
    fi
    if [ $i -eq 6 ] 
    then
	CUR_JOB=aRNASeq_${PROJNAME}_idx_bam
	if [ -z $PREV_JOB ]
	then
	    CMD="qsub -N $CUR_JOB RNAseq_idx_bam_${PROJNAME}.qsub"
        else
	    CMD="qsub -hold_jid $PREV_JOB -N $CUR_JOB RNAseq_idx_bam_${PROJNAME}.qsub"
        fi
	echo "STEP$i: $CMD"
	eval $CMD
	PREV_JOB=$CUR_JOB
    fi
    if [ $i -eq 7 ] 
    then
	CUR_JOB=RNASeq_${PROJNAME}_fc
	if [ -z $PREV_JOB ]
	then
	    CMD="qsub -N $CUR_JOB RNAseq_fc_${PROJNAME}.qsub"
        else
	    CMD="qsub -hold_jid $PREV_JOB -N $CUR_JOB RNAseq_fc_${PROJNAME}.qsub"
        fi
	echo "STEP$i: $CMD"
	eval $CMD
	PREV_JOB=$CUR_JOB
    fi
    if [ $i -eq 8 ] 
    then
        CUR_JOB=RNASeq_${PROJNAME}_deseq2
	if [ -z $PREV_JOB ]
	then
	    CMD="qsub -N $CUR_JOB RNAseq_deseq2_${PROJNAME}.qsub"
        else
	    CMD="qsub -hold_jid $PREV_JOB -N $CUR_JOB RNAseq_deseq2_${PROJNAME}.qsub"
        fi
	echo "STEP$i: $CMD"
	eval $CMD
	PREV_JOB=$CUR_JOB
    fi
    if [ $i -eq 9 ] 
    then
	CUR_JOB=aRNASeq_${PROJNAME}_gsea
	if [ -z $PREV_JOB ]
	then
	    CMD="qsub -N $CUR_JOB RNAseq_gsea_${PROJNAME}.qsub"
        else
	    CMD="qsub -hold_jid $PREV_JOB -N $CUR_JOB RNAseq_gsea_${PROJNAME}.qsub"
        fi
	echo "STEP$i: $CMD"
	eval $CMD
	PREV_JOB=$CUR_JOB
    fi
    
done
