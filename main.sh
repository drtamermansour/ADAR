mkdir ADAR && cd ADAR
work_dir=$(pwd)
## adar_probes.fa: The probe sequences, copied from email (word doc), headers edited manually, and sequences converted to capital letters  
## adars_seq.fa: The ADAR genes sequences, copied from email, headers edited manually, and sequences converted to capital letters
## D_pealeii_ORFs.fa: The sequid transcriptome. copied from the dropbox 
grep -i "ADAR" D_pealeii_ORFs.fa | awk '{print $1}' | grep -A1 -Fwf - D_pealeii_ORFs.fa | grep -v "\-\-" > adar.D_pealeii_ORFs.fasta
## unwrap the bait sequences 
#awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < adars_seq.fa | tail -n +2 > adars_unwrapped.fa

## download RNAseq data
mkdir data && cd data
module load SRAToolkit/2.8.2
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR1522987 
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR1522988
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR1725213
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR1725235
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR1725236	
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR1725164
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR1725163
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR1725171
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR1725169
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR1725167
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR1725172

#target=ADAR1 #ADAR2 #ADAR3
module swap GNU GNU/4.9
module load trinity/2.4.0
> matching_baits.txt
> gene_exp.report
> gene_exp.read.report
for target in ADAR1 ADAR2 ADAR3;do
  mkdir -p $target
  grep -A1 $target adars_seq.fa > $target/${target}_seq.fa

  ## convert the bait to kmers
  bait=$work_dir/$target/${target}_seq.fa
  k=25
  while read header; read line;do
    step=1;
    end=$((${#line} +1 - $k));
    for i in $(seq 1 $step $end);do
      #echo $step $end $i;
      kmer=`expr substr $line $i $k`;
      echo $kmer;done
      #echo -e "$header.position.$i\n$kmer";done
  done < $bait | tr "atcg" "ATCG" | sort | uniq > $bait.kmers.txt
  ## find reverse complement kmers
  cat $bait.kmers.txt | rev | tr "ATGC" "TACG" > $bait.rckmers.txt
  ## make one list of all possible kmers
  cat $bait.kmers.txt $bait.rckmers.txt | sort | uniq > $bait.Allkmers.txt ## 7434, 4668, 3756

  ## fadars_seq.faind matching reads
  grep -B1 -F -f $bait.Allkmers.txt D_pealeii_ORFs.fa | grep -v "^\-\-" > $bait.D_pealeii_ORFs.fa ## each bait find single match
  grep -B1 -F -f $bait.Allkmers.txt adar_probes.fa | grep -v "^\-\-" > $bait.adar_probes.fa       ## each bait find single match
  grep -B1 -F -f $bait.Allkmers.txt adars_seq.fa | grep -v "^\-\-" > $bait.adars_seq.fa           ## each bait find single match

  > $bait.trans_trials.txt
  > $bait.probe_trials.txt
  > $bait.ORF_trials.txt
  while read kmer;do
    grep $kmer $bait.D_pealeii_ORFs.fa >> $bait.trans_trials.txt 
    grep $kmer $bait.adar_probes.fa >> $bait.probe_trials.txt
    grep $kmer $bait.adars_seq.fa >> $bait.ORF_trials.txt
  done < $bait.Allkmers.txt
  wc -l $bait.*_trials.txt >> matching_baits.txt

  #sample=SRR1522987 #SRR1522988 SRR1725213 SRR1725235 SRR1725236 SRR1725164 SRR1725163 SRR1725171 SRR1725169 SRR1725167 SRR1725172
  for sample in SRR1522987 SRR1522988 SRR1725213;do
    prey=$work_dir/$target/$sample/$target
    mkdir -p $work_dir/$target/$sample
    grep -B1 -A2 -F -f $bait.Allkmers.txt data/${sample}_1.fastq | grep -v "^\-\-$" > $prey.${sample}_1.fq # 816, 1168, 4280 #   , 16152, 
    echo $(( $(cat $prey.${sample}_1.fq | wc -l) / 4 )) $target/$sample/R1 >> gene_exp.report
    grep -B1 -A2 -F -f $bait.Allkmers.txt data/${sample}_2.fastq | grep -v "^\-\-$" > $prey.${sample}_2.fq # 812, 1156, 4240 #   , 15960,
    echo $(( $(cat $prey.${sample}_2.fq | wc -l) / 4 )) $target/$sample/R2 >> gene_exp.report

  ## Each reference ADAR gene was converted into all possible 25 kmers. ADAR1, ADAR2, ADAR3 produced 3717, 2334, and 1878 kmers respectively.
  ## I found no kmer overlap between the 3 genes! And all kmers are unique through the whole transcriptome

  ## D_pealeii transcriptome has 3 isoforms that share sequences with an ADAR gene: 
  ## comp134400_c0_seq1, comp135521_c1_seq2, comp136967_c2_seq2 respectively.
  #   * ADAR1 has 1398 (out of 3717) matches with comp134400_c0_seq1 and 976 matches with ADAR1_probe (perfect match)
  #   * ADAR2 has 1231 (out of 2334) matches with comp135521_c1_seq2 and 958 matches with ADAR2_probe (perfect match)
  #   * ADAR3 has 1878 (out of 1878) matches with comp136967_c2_seq2 and 996 matches with ADAR3_probe (perfect match)

  ## Using a similar approach, I identified all the reads that share any kmer sequence with ADAR genes: 
  ## For example: (SRR1522987: 204, 291, 1065) & (SRR1522988: 176, 4014, 1815) reads. 
  ## These numbers should be corrected for the gene length (3741, 2358, 1902 bp respectively) to be: 
  ## (SRR1522987: 102, 233, 1065) & (SRR1522988: 88, 3211, 1815) reads. 
  ## Notes that these no are no normalized for the read count so you can not compare the same gene across experiments. 
  ## We are only trying to find the "relative" expression of the 3 genes in the 3 experiments.  

    ## secondary baiting
    k=25
    cat $prey.${sample}_[12].fq | paste - - - - | awk '{print $2}' | while read line;do
      step=1;
      end=$((${#line} +1 - $k));
      for i in $(seq 1 $step $end);do
        #echo $step $end $i;
        kmer=`expr substr $line $i $k`;
        echo $kmer;done
        #echo -e "$header.position.$i\n$kmer";done
    done | sort | uniq > $prey.read.kmers.txt
    cat $prey.read.kmers.txt | rev | tr "ATGC" "TACG" > $prey.read.rckmers.txt
    cat $prey.read.kmers.txt $prey.read.rckmers.txt | sort | uniq > $prey.read.Allkmers.txt ## 11798, 11070,19626

    grep -B1 -A2 -F -f $prey.read.Allkmers.txt data/${sample}_1.fastq | grep -v "^\-\-$" > $prey.read.${sample}_1.fq # 55664 1162212 214588
    echo $(( $(cat $prey.read.${sample}_1.fq | wc -l) / 4 )) $target/$sample/R1 >> gene_exp.read.report
    grep -B1 -A2 -F -f $prey.read.Allkmers.txt data/${sample}_2.fastq | grep -v "^\-\-$" > $prey.read.${sample}_2.fq # 56824  202740 954028
    echo $(( $(cat $prey.read.${sample}_2.fq | wc -l) / 4 )) $target/$sample/R2 >> gene_exp.read.report

    ## De novo assembly of all reads
    cat $prey.read.${sample}_[12].fq > $prey.read.${sample}.fq
    Trinity --seqType fq --max_memory 24G --CPU 6 --output $prey.trinity --single $prey.read.${sample}.fq &> $prey.trinity.log

    Trinity --seqType fq --max_memory 24G --CPU 6 --output $prey.trinity.noNorm  --no_normalize_reads --single $prey.read.${sample}.fq &> $prey.trinity.noNorm.log

    Trinity --seqType fq --max_memory 24G --CPU 6 --output $prey.trinity.noMerge  --no_normalize_reads --no_path_merging --single $prey.read.${sample}.fq &> $prey.trinity.noMerge.log

    Trinity --seqType fq --max_memory 24G --CPU 6 --output $prey.trinity.relax  --no_normalize_reads --min_contig_length 120 --min_glue 1 --path_reinforcement_distance 1 --no_path_merging --single $prey.read.${sample}.fq &> $prey.trinity.relax.log
    cat $prey.trinity.relax/Trinity.fasta | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | tail -n +2 > $prey.trinity.relax/Trinity_unwrapped.fasta
    cat $bait.adar_probes.fa > ~/temp/$target.$sample.probe_Trinity.fasta
    cat $prey.trinity.relax/Trinity.fasta | awk '{print $1}' >> ~/temp/$target.$sample.probe_Trinity.fasta
    ## blast aganist the target gene sequence then download the fasta with successful queries 
  done

done
 
## ADAR2
## pos 43
bait="GGTAATTCTTCCGAGGATCCTAAAATTAAGATTGA"
revBait=$(echo $bait | rev | tr "ATGC" "TACG")
for asm in ADAR2/SRR*/ADAR2.trinity.relax/Trinity_unwrapped.fasta;do
echo -e "\n$asm";
cat $asm | grep -B1 $bait | grep -v "^\-\-$" | awk '{print $1}'
cat $asm | grep -B1 $revBait | grep -v "^\-\-$" | awk '{print $1}' | paste - - | while read header seq;do echo $header; echo $seq | rev | tr "ATGC" "TACG";done 
done > ADAR2.3samples.pos43

grep GGTAATTCTTCCGAGGATCCTAAAATTAAGATTGA ADAR2/SRR*/ADAR2.trinity.relax/Trinity_unwrapped.fasta >> ADAR2.3samples.pos43
grep TCAATCTTAATTTTAGGATCCTCGGAAGAATTACC ADAR2/SRR*/ADAR2.trinity.relax/Trinity_unwrapped.fasta >> ADAR2.3samples.pos43


