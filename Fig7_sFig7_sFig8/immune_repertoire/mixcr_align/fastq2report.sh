###### input file #######
mask="/media/transport2/backup/crc_immune/fastq/Y01511D1/Y01511D1.barcodeToPos.h5"
fastq1="/media/transport2/backup/crc_immune/fastq/Y01511D1/Y01511D1_1.raw.gz"
fastq2="/media/transport2/backup/crc_immune/fastq/Y01511D1/Y01511D1_2.raw.gz"
name="Y01511D1"


#######  input file #######
outfastq="/storage/liuyi/tmp/crc_0205/Y01511D1_barcodemap.modify.fq.gz"
out="/storage/liuyi/tmp/crc_0205/Y01511D1_barcodemap.fq.gz"  # tmp fastq
#########################



mixcr="/home/liuyi/02.software/mixcr/mixcr"

# stmap is short for software ST_barcodemap  https://github.com/STOmics/ST_BarcodeMap
/usr/local/bin/stmap --in $mask\    
        --in1 $fastq1\
        --in2 $fastq2\
        --out $out\
        --mismatch 1\
        --thread 10\
        --umiLen -1

# process_fastq is custom script by C
/storage/liuyi/scirpt/process_v2_fastq/process_fastq $out $outfastq && rm  $out

[ -d "./${name}" ] || mkdir -p ./${name}
echo $name

$mixcr align \
    -t 20 \
    --species hsa \
    -f \
    -p rna-seq \
    -Xmx64g \
    -OsaveOriginalReads=true \
    --keep-non-CDR3-alignments \
    -OallowPartialAlignments=true \
    -OvParameters.geneFeatureToAlign="VGeneWithP" \
    -OdParameters.geneFeatureToAlign="DRegionWithP" \
    -OjParameters.geneFeatureToAlign="JRegionWithP" \
    --report ./tmp/${name}.align.report.txt \
    --json-report ./tmp/${name}.align.report.json \
    $outfastq \
    ./${name}/${name}.vdjca

vdjca="./${name}/${name}.vdjca"

mixcr_assamblePartial(){
    $mixcr assemblePartial \
      -f \
      -OminimalAssembleOverlap=15 \
      --report ./$name/$name.assemblePartial.report.txt \
      $1 \
      $2 
}

mixcr_assamblePartial $vdjca ./$name/$name.passembled.1.vdjca
mixcr_assamblePartial ./$name/$name.passembled.1.vdjca ./$name/$name.passembled.2.vdjca

$mixcr extend \
    -Xmx200g \
    -f \
    --report ./$name/$name.extend.report.txt \
    --json-report ./$name/$name.extend.report.json \
    ./$name/$name.passembled.2.vdjca \
    ./$name/$name.passembled.extended.vdjca


$mixcr assemble \
    -Xmx200g \
    -f \
    -OassemblingFeatures='CDR3' \
    -OseparateByC=true\
    -a \
    --report ./${name}/${name}.assemble.report.txt \
    --json-report ./${name}/${name}.assemble.report.json \
    ./$name/$name.passembled.extended.vdjca \
    ./$name/$name.clna


$mixcr exportClones\
    -f\
    --dont-split-files \
    --prepend-columns \
    -topChains \
    -isotype primary\
    ./$name/$name.clna ./$name/$name.contigs.tsv

$mixcr exportAlignments \
    --chains IG\
    -readIds \
    -descrsR1\
    -cloneId\
    -f \
    ./$name/$name.clna ./$name/$name.align.tsv

$mixcr exportAlignments \
    --drop-default-fields\
    -vHitsWithScore\
    -dHitsWithScore\
    -jHitsWithScore\
    -cHitsWithScore\
    -readIds \
    -descrsR1\
    -cloneId\
    -f \
    ./$name/$name.clna ./$name/$name.align.all.tsv
