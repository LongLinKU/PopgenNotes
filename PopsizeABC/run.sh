module load python2/2.7.18
# module load gsl/2.7.1
module load gsl/2.5 perl/5.26.3 bcftools
export LD_LIBRARY_PATH=/opt/software/gsl/2.7.1:/opt/software/gsl/2.7.1/lib
export PYTHONPATH=/projects/seqafrica/people/dmv161/muskox/popsizeABC/bin/py2.7.18.packages/
scriptFolder=/projects/seqafrica/people/dmv161/muskox/popsizeABC/bin
outputDir=/maps/projects/seqafrica/people/dmv161/elephant/popsizeABC/res210000/dp10/LO
vcf_file=/maps/projects/seqafrica/people/dmv161/elephant/popsizeABC/dataset/LoxAfr4_elephant_no1stdegree_variable_sites_nomultiallelics_noindels_10dp_3het.vcf.gz
ped_file=/maps/projects/seqafrica/people/dmv161/elephant/popsizeABC/pop_and_ped/dp10/elephant.dp10.ped

#ped file shall have the same order of individuals as in vcf.
check_ped=$(paste <(cut -f 2 $ped_file) <(bcftools query -l $vcf_file) | awk '$1!=$2{print $0}')
if [ -n "$check_ped" ]; then
    echo "The ped file has different order of samples in vcf, which is wrong."
    exit 1
fi

if [ ! -d "$outputDir" ]; then
    mkdir -p "$outputDir"
fi

#step 1
pop=LO
inds_file=/maps/projects/seqafrica/people/dmv161/elephant/popsizeABC/pop_and_ped/dp10/$pop.inds.list
num_lines=$(wc -l < "$inds_file")
# Calculate two times the number of lines
n=$((num_lines * 2))
tol=0.1
gen_time=31
runs=210000
segs=100
#python $scriptFolder/comp_stat1/stat_from_vcf_ex1.py $pop $inds_file $vcf_file $ped_file $outputDir

#step 2
echo "starts step2"
#python $scriptFolder/comp_stat1/simul_data_ex1.py $outputDir/ex2_$pop $runs $segs $n


#step3
echo "starts step3"
infile_params=$(ls $outputDir/ex2_${pop}*params)
infile_stat=$(ls $outputDir/ex2_${pop}*stat)
infile_obs=$(ls $outputDir/ex1*stat)
module load  gcc/11.2.0 R/4.2.2
Rscript $scriptFolder/estim/abc_ex1.R $infile_params $infile_stat $pop $n $infile_obs $tol $outputDir $gen_time

#step4
echo "starts step4"
echo Rscript $scriptFolder/estim/abc_cv_ex1.R $infile_params $infile_stat $pop $n $tol $outputDir
