# Simulation studies

We conducted simulation studies in order to assess the performance of different colocalization methods.

## Prepare genotype

We randomly sampled five 1-Mb loci from the genome and extracted their genotypes for 25,000 and 1,000 (non-overlapping) [UK Biobank](https://www.ukbiobank.ac.uk) European ancestry individuals for simulating trait 1 and trait 2. 
For each locus, the LD matrix was calculated using [PLINK](https://www.cog-genomics.org/plink/).

```
## select individuals
head -26000 ~/scratch/UKB_Geno/1.fam > UKB_A.fam

## sample genotypes
while read a b c d; do sed "s/LOCI/$a/g" LOCI.sh | sed "s/CHR/$b/g" | sed "s/ST/$c/g" | sed "s/ED/$d/g" > $a\.sh; done < random.gtf
for i in $(seq 1 5); do ./Locus$i\.sh; done
```

## Simulate traits and obtain summary statistics

We randomly sampled `KC` shared causal variants and `KS` specific causal variants for both traits in each locus. 
For example, when `KC = 0` and `KS = 1`, there was one causal variant for trait 1 and a different causal variant for trait 2; When `KC = 1` and `KS = 0`, there was one causal variant shared by both traits. 
With simulated traits, we performed GWAS using [GCTA](https://yanglab.westlake.edu.cn/software/gcta/#Overview) to obtain summary statistics. We repeated this process 50 times for each setting in each locus.

```
./simulate_traits.sh LOCI KC KS HQ HG
```

## Perform colocalization

With LD information and simulated summary statistics, 
we performed colocalization analysis with different methods using default settings. 

```
./run_SH.sh LOCI KC_KS_HQ_HG 1000 25000
./run_PWCoCo.sh LOCI KC_KS_HQ_HG
./run_coloc.sh LOCI KC_KS_HQ_HG
./run_ecaviar.sh LOCI KC_KS_HQ_HG
```

## Result visualization

SharePro achieved high power and well-controlled false discoveries across most simulation settings.

```
Rscript ../doc/plot_sim.R
```