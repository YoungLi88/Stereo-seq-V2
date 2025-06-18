mixcr="/home/liuyi/02.software/mixcr/mixcr"
name=$1
vdjca=$2

[ -d "./${name}" ] || mkdir -p ./${name}


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
