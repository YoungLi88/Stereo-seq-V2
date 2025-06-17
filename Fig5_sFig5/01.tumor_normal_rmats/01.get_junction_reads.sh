rmatsresult="/storage/liuyi/08.stereo_v2/HER2_breast/bam/tumor_normal_dup_rmats/out"
sample="A5"
bam="/storage/liuyi/08.stereo_v2/HER2_breast/bam/C03637A5.dup.total.bam"
outdir="./as"
cpu=10
binsize=100
gem="/home/lee/project/v2/tnbc/gem/C03637A5.lasso.bin1.Label1.gem.gz"
bin_type="normal"


mkdir -p ${outdir}/02.extract ${outdir}/03.intron
/home/lee/miniconda3/envs/rapids/bin/python /storage/liuyi/scirpt/spatial_AS/script/extract_intron.py $rmatsresult ${outdir}/02.extract $sample || exit 1
/home/lee/miniconda3/envs/rapids/bin/python /storage/liuyi/08.stereo_v2/HER2_breast/bam/intron_count.py $bam ${outdir}/02.extract/${sample}.Intron.set ${outdir}/03.intron/${sample}.intron.gz ${outdir}/02.extract/${sample}.RI.set ${outdir}/03.intron/${sample}.ri.gz ${cpu}
/home/lee/miniconda3/envs/rapids/bin/python /storage/liuyi/scirpt/spatial_AS/script/bulid_matrix.py ${outdir}/03.intron/${sample}.intron.gz ${outdir}/03.intron/count_bin${binsize} ${binsize} ${gem} ${bin_type}
/home/lee/miniconda3/envs/rapids/bin/python /storage/liuyi/scirpt/spatial_AS/script/bulid_matrix.py ${outdir}/03.intron/${sample}.ri.gz ${outdir}/03.intron/ri_bin${binsize} ${binsize} ${gem} ${bin_type}
