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

# Adapted by Philippe Cubry to deal with 163 generated individuals

for i in {1..12}
do
echo "##fileformat=VCFv4.2" > 163_NEW.chr$i.macs.output.vcf
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tIND1\tIND2\tIND3\tIND4\tIND5\tIND6\tIND7\tIND8\tIND9\tIND10\tIND11\tIND12\tIND13\tIND14\tIND15\tIND16\tIND17\tIND18\tIND19\tIND20\tIND21\tIND22\tIND23\tIND24\tIND25\tIND26\tIND27\tIND28\tIND29\tIND30\tIND31\tIND32\tIND33\tIND34\tIND35\tIND36\tIND37\tIND38\tIND39\tIND40\tIND41\tIND42\tIND43\tIND44\tIND45\tIND46\tIND47\tIND48\tIND49\tIND50\tIND51\tIND52\tIND53\tIND54\tIND55\tIND56\tIND57\tIND58\tIND59\tIND60\tIND61\tIND62\tIND63\tIND64\tIND65\tIND66\tIND67\tIND68\tIND69\tIND70\tIND71\tIND72\tIND73\tIND74\tIND75\tIND76\tIND77\tIND78\tIND79\tIND80\tIND81\tIND82" >>  163_NEW.chr$i.macs.output.vcf

len=20000000

# grab sites from macs
# pos=sprintf("%0.f",($3*len*i)) : %0.f means round decimals, sprintf assigns output to a variable
# so multiply the position ($3) by length of window simulated and then *i
# so that each group goes from 1>> 10Mb instead of each window going from 1>>100kb
# end up with all 100 blocks in a group (chromosome) with non-overlapping positions
# split($5,geno,sep="") : use awk to split up the genotypes and then write them out separated by a | and tab between each individual (end up with 10 columns)

grep "SITE:"  < 163_NEW.chr$i.macs.out| \
awk -v len=$len -v chr=Chr$i -v i=$i '{OFS="\t"; pos=sprintf("%0.f",($3*len)); split($5,geno,sep=""); print chr,pos,".","A","G",".","PASS",".","GT",geno[1]"|"geno[2],geno[3]"|"geno[4],geno[5]"|"geno[6],geno[7]"|"geno[8],geno[9]"|"geno[10],geno[11]"|"geno[12],geno[13]"|"geno[14],geno[15]"|"geno[16],geno[17]"|"geno[18],geno[19]"|"geno[20],geno[21]"|"geno[22],geno[23]"|"geno[24],geno[25]"|"geno[26],geno[27]"|"geno[28],geno[29]"|"geno[30],geno[31]"|"geno[32],geno[33]"|"geno[34],geno[35]"|"geno[36],geno[37]"|"geno[38],geno[39]"|"geno[40],geno[41]"|"geno[42],geno[43]"|"geno[44],geno[45]"|"geno[46],geno[47]"|"geno[48],geno[49]"|"geno[50],geno[51]"|"geno[52],geno[53]"|"geno[54],geno[55]"|"geno[56],geno[57]"|"geno[58],geno[59]"|"geno[60],geno[61]"|"geno[62],geno[63]"|"geno[64],geno[65]"|"geno[66],geno[67]"|"geno[68],geno[69]"|"geno[70],geno[71]"|"geno[72],geno[73]"|"geno[74],geno[75]"|"geno[76],geno[77]"|"geno[78],geno[79]"|"geno[80],geno[81]"|"geno[82],geno[83]"|"geno[84],geno[85]"|"geno[86],geno[87]"|"geno[88],geno[89]"|"geno[90],geno[91]"|"geno[92],geno[93]"|"geno[94],geno[95]"|"geno[96],geno[97]"|"geno[98],geno[99]"|"geno[100],geno[101]"|"geno[102],geno[103]"|"geno[104],geno[105]"|"geno[106],geno[107]"|"geno[108],geno[109]"|"geno[110],geno[111]"|"geno[112],geno[113]"|"geno[114],geno[115]"|"geno[116],geno[117]"|"geno[118],geno[119]"|"geno[120],geno[121]"|"geno[122],geno[123]"|"geno[124],geno[125]"|"geno[126],geno[127]"|"geno[128],geno[129]"|"geno[130],geno[131]"|"geno[132],geno[133]"|"geno[134],geno[135]"|"geno[136],geno[137]"|"geno[138],geno[139]"|"geno[140],geno[141]"|"geno[142],geno[143]"|"geno[144],geno[145]"|"geno[146],geno[147]"|"geno[148],geno[149]"|"geno[150],geno[151]"|"geno[152],geno[153]"|"geno[154],geno[155]"|"geno[156],geno[157]"|"geno[158],geno[159]"|"geno[160],geno[161]"|"geno[162],geno[163]"|"geno[16]}' >> \
163_NEW.chr$i.macs.output.vcf 

done
