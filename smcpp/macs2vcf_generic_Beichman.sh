#! /bin/bash
#$ -cwd
#$ -l h_rt=24:00:00,h_data=2G,highp
#$ -N m2vSMCpp
#$ -o output.txt
#$ -e error.txt
#$ -m abe
#$ -M ab08028
#$ -t 1-200

# February 14, 2017
# author: Annabel Beichman
# contact: annabel.beichman@gmail.com
# Purpose: convert MaCS output to vcf format 

j=$SGE_TASK_ID
model= # name of model you're running (e.g. gutenkunst)
rundate=# date you ran simulations

flag= # flag used to identify model 
mkdir ${flag}_VCFs_${rundate}_blocks
out=${flag}_VCFs_${rundate}_blocks

                   
echo $j           
mkdir $out/group_${j}_${flag}_vcfs_${rundate}
gpout=$out/group_${j}_${flag}_vcfs_${rundate}                                                                         
ids= # population ids, e.g. yri ceu chb


len=100000

for id in $ids
do
echo $id
### need to only have 1 # before CHROM line. If two ## you lose individuals 
for i in {1..100}
do
echo "##fileformat=VCFv4.2" > $gpout/group_${j}_block_${i}.$flag.$id.output.vcf
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tIND1\tIND2\tIND3\tIND4\tIND5\tIND6\tIND7\tIND8\tIND9\tIND10" >> $gpout/group_${j}_block_${i}.$flag.$id.output.vcf

# grab sites from macs
# pos=sprintf("%0.f",($3*len*i)) : %0.f means round decimals, sprintf assigns output to a variable
# so multiply the position ($3) by length of window simulated and then *i
# so that each group goes from 1>> 10Mb instead of each window going from 1>>100kb
# end up with all 100 blocks in a group (chromosome) with non-overlapping positions
# split($5,geno,sep="") : use awk to split up the genotypes and then write them out separated by a | and tab between each individual (end up with 10 columns)

grep "SITE:" group_$j.${model}/*/group_${j}_block_${i}.${flag}.${id}.macs.macsFormat.OutputFile.$rundate.txt | \
awk -v len=$len -v chr=$j -v i=$i '{OFS="\t"; pos=sprintf("%0.f",($3*len)); split($5,geno,sep=""); print chr,pos,".","A","G",".","PASS",".","GT",geno[1]"|"geno[2],geno[3]"|"geno[4],geno[5]"|"geno[6],geno[7]"|"geno[8],geno[9]"|"geno[10],geno[11]"|"geno[12],geno[13]"|"geno[14],geno[15]"|"geno[16],geno[17]"|"geno[18],geno[19]"|"geno[20]}' >> \
$gpout/group_${j}_block_${i}.$flag.$id.output.vcf

done
done

sleep 10m