#!/bin/bash
sampleID="$1"
fastq="$2"
species="mmu"

# software
mixcrDir="/home/liuyi/02.software/mixcr_4.6"

mkdir -p /media/transport2/BGItrans/liyang13/backup/thymus_v1_v2/mixcr_for_v2/${sampleID}/

# mixcr analyze
${mixcrDir}/mixcr align \
	-p rna-seq \
	-Xmx50g \
	--threads 40 \
	--species ${species} \
	--keep-non-CDR3-alignments \
	-f \
	-OsaveOriginalReads=true \
	-OallowPartialAlignments=true \
	-OvParameters.geneFeatureToAlign="VGeneWithP" \
	-OdParameters.geneFeatureToAlign="DRegionWithP" \
	-OjParameters.geneFeatureToAlign="JRegionWithP" \
	--report /media/transport2/BGItrans/liyang13/backup/thymus_v1_v2/mixcr_for_v2/${sampleID}/${sampleID}.align.report.txt \
	--json-report /media/transport2/BGItrans/liyang13/backup/thymus_v1_v2/mixcr_for_v2/${sampleID}/${sampleID}.align.report.json \
	$fastq /media/transport2/BGItrans/liyang13/backup/thymus_v1_v2/mixcr_for_v2/${sampleID}/${sampleID}.vdjca

