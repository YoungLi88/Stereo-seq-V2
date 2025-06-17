sample=""
input_bam=""
out_bam=""
matrics=""
java -jar /home/liuyi/02.software/picard.jar MarkDuplicates\
	-I $input_bam \
	-O $out_bam \
	-M $matrics\
	--BARCODE_TAG CB\
	--REMOVE_DUPLICATES true
