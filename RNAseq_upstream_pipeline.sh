###RNA_seq analysis and software install
##conda create -n rnaseq
##conda activate rnaseq
##conda install fastp hisat2 samtools stringtie
###start analysis
###Quality control
conda activate rnaseq
cd /datapool/hechaoying/anzhenghong/data/third/Vph
for i in 0 1 2 4 6 8;do
for j in 11 12 13;do
  {
  fastp -i raw/Vph10T16${i}dF${j}_1.fq.gz -I raw/Vph10T16${i}dF${j}_2.fq.gz -o clean/Vph10t16/Vph10T16${i}dF${j}_1_fastp.fq.gz -O clean/Vph10t16/Vph10T16${i}dF${j}_2_fastp.fq.gz --thread 4 -q 20 --html fastp_report/Vph10t16/Vph10T16${i}dF${j}.html --json fastp_report/Vph10t16/Vph10T16${i}dF${j}.json  -l 150
  } 
done
done
###mapping index
cd /datapool/hechaoying/anzhenghong/data/third/Vph/transcriptome
hisat2-build -p 4  Vph.fa Vph_transcriptomedb
####mapping
cd /datapool/hechaoying/anzhenghong/data/third/Vph/mapping_stage
for i in 0 1 2 4 6 8 ;do
for j in 11 12 13;do
  {
hisat2  -p 4  -x  ../transcriptome/Vph_transcriptomedb -1 ../clean/Vph10t16/Vph10T16${i}dF${j}_1_fastp.fq.gz  -2 ../clean/Vph10t16/Vph10T16${i}dF${j}_2_fastp.fq.gz -S Vph10T16${i}dF${j}.sam 2>Vph10T16${i}dF${j}.log 
samtools view -S Vph10T16${i}dF${j}.sam -b > Vph10T16${i}dF${j}.bam
samtools sort Vph10T16${i}dF${j}.bam -o Vph10T16${i}dF${j}_sorted.bam
samtools index Vph10T16${i}dF${j}_sorted.bam
}
done
done
##quantity
cd /datapool/hechaoying/anzhenghong/data/third/Vph/stringtie
for i in 0 1 2 4 6 8 ;do
for j in 11 12 13;do
 {
stringtie  -p 8 -e -G ../transcriptome/Vph.transdecoder.gff3 -o Vph10T16${i}dF${j}.gtf -A Vph10T16${i}dF${j}.tsv ../mapping_stage/Vph10T16${i}dF${j}_sorted.bam
}
done
done
###prepDE.py was downloaded from https://github.com/gpertea/stringtie/tree/master
##prepDE.py requires a sample_list(Vph10T16_samplelist.txt), where the first column contains sample names and the second column contains .gtf file paths.
python3 ./prepDE.py -i Vph10T16_samplelist.txt
