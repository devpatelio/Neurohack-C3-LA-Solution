version 1.0
# Workflow Developer: Ramanandan Prabhakaran
# Step1: Annotate VCF using VEP tool
# Step2: Select only variants that are passed
# Step3: Prepare filter ready input file for VEP
# Step4: Filter by gene symbol list
# Step5: Filter by consequence like coding sequence variant or stop lost or start gain
# Step6: Filter by sift prediction like tolerated 
# Step7: Filter by sift score less than 0.5
# Step8: Filter by allele frequency less than 0.01
# VCF 2 BCF conversion
# BCF 2 PLINK genotype generation
workflow Annotate_Filter_VCF 
{
   input {
     File vcfFile
     File reference
     Int? alleleFrequency
     String? consequenceTerm
     String? siftPrediction
     String? siftScore
     File? geneList
   }
   call AnnotateVCF{
     input: vcfFile = vcfFile,
            referenceFasta = reference
   }
   call FilterPASS{
     input: vcfFile = AnnotateVCF.annotatedVCF
   }
   call VCF2TXT{
     input: vcfFile = Filter_PASS.Passed_variants
   }
   call FilterByGene{
     input: vepInput = VCF2TXT.Filter_ready_txt_file,
            gene = geneList
   }
   call FilterByConsequence{
     input: vepInput = VCF2TXT.Filter_ready_txt_file,
            consequence = consequenceTerm
   }
   call FilterByAllelFrequency{
     input: vepInput = VCF2TXT.Filter_ready_txt_file,
            alleleFrequency = alleleFrequency
   }
   call FilterBySiftPrediction{
     input: vepInput = VCF2TXT.Filter_ready_txt_file,
            siftPrediction = siftPrediction
   }
   call FilterBySiftScore{
     input: vepInput = VCF2TXT.Filter_ready_txt_file,
            siftScore = siftScore
   }
   call VCF2BCF {
     input: vcfFile = AnnotateVCF.annotatedVCF
   }
   call VCF2Plink {
     input: bcfFile = VCF2BCF.bcfFile,
            referenceFasta = reference
   }
  call PruneVariants {
    input: genotypeFile=VCF2Plink.genotypeFile
  }
}
   
task AnnotateVCF 
{
   input {
     File vcfFile
     File vcfIndex
     File referenceFasta
     String basename = basename("~{vcfFile}", ".vcf.gz")
   }
   command <<<
     vep --offline --dir ~{vepCacheDir} -i ~{vcfFile} --fasta ~{referenceFasta} \
         --species homo_sapiens --assembly ~{ncbiBuild} -o ~{basename}.vep.vcf.gz \
         --vcf --compress_output bgzip --ccds --uniprot --hgvs --symbol --numbers \ 
         --domains --gene_phenotype --canonical --protein --biotype --uniprot
   >>>
   runtime {
     docker: "ensemblorg/ensembl-vep:latest"
     memory: "16 GB"
     cpu : "2"
     retry : "2"
   }
   output {
     File annotatedVCF = "~{basename}.vep.vcf.gz"
   }
}

task FilterPASS 
{
   input {
     File vcfFile
     String basename = basename("~{vcfFile}", ".vcf.gz")
   }
   command <<<
     bcftools view -f PASS ~{vcfFile} > ~{basename}.PASS.vcf.gz
   >>>
   runtime {
     docker: "biocontainers/bcftools:latest"
     memory: "16 GB"
     cpu : "2"
     retry : "2"
   }
   output {
     File Passed_variants = "~{basename}.PASS.vcf.gz"
   }
}

task VCF2BCF
{
   input {
     File vcfFile
     String basename = basename("~{vcfFile}", ".vcf.gz")
   }
   command <<<
     vcftools --vcf ~{vcfFile} --recode-bcf-to-stream > ~{basename}.bcf
   >>>
   runtime {
     docker: "biocontainers/vcftools:latest"
     memory: "16 GB"
     cpu : "2"
     retry : "2"
   }
   output {
     File bcfFile = "~{basename}.bcf"
   }
}

task VCF2PLINK
{
   input {
     File bcfFile
     String basename = basename("~{bcfFile}", ".bcf")
   }
   command <<<
     plink --noweb --bcf ~{bcfFile} --keep-allele-order \
      --vcf-idspace-to _ \
      --const-fid \
      --allow-extra-chr 0 \
      --split-x b37 no-fail \
      --make-bed \
      --out ~{basename}.genotypes
   >>>
   runtime {
     docker: "asherkhb/plink:latest"
     memory: "16 GB"
     cpu : "2"
     retry : "2"
   }
   output {
     File genotypeFile = "~{basename}.genotypes"
   }
}

task PruneVariants
{
   input {
     File genotypeFile
     String basename = basename("~{genotypes}", ".genotype")
   }
   command <<<
     plink --noweb --bfile ~{genotypeFile} \
      --maf 0.10 --indep 50 5 1.5 \
      --out Pruned/~{basename}.genotypes
     plink --noweb --bfile ~{genotypeFile} \
      --extract Pruned/~{basename}.genotypes.prune.in \
      --make-bed \
      --out Pruned/~{basename}.genotypes
   >>>
   runtime {
     docker: "asherkhb/plink:latest"
     memory: "16 GB"
     cpu : "2"
     retry : "2"
   }
   output {
     File prunedGenotypeFile = "~{basename}.genotypes"
   }
}

task VCF2TXT
{
   input {
     File vcfFile
     String basename = basename("~{vcfFile}", ".vcf.gz")
}
   command <<<
     vep -i ~{vcfFile} -o ~{basename}.txt --cache --offline
   >>>
   runtime {
     docker: "ensemblorg/ensembl-vep:latest"
     memory: "16 GB"
     cpu : "2"
     retry : "2"
   }
   output
   {
     File Filter_ready_txt_file = "~{basename}.txt"
   }
}

task FilterByGene
{
   input {
     File filter_ready_vcf_file
     String basename = basename("~{filter_ready_vcf_file}", ".txt")
     File? geneList
   }
   command <<<
     filter_vep -i ~{filter_ready_vcf_file} -o ~{basename}.GeneList.txt --filter "SYMBOL in ~{geneList}"
   >>>
   runtime {
     docker: "ensemblorg/ensembl-vep:latest"
     memory: "16 GB"
     cpu : "2"
     retry : "2"
   }
   output
   {
     File Gene_specific_variants = "~{basename}.GeneList.txt"
   }
}

task FilterByConsequence
{
   input {
     File filter_ready_vcf_file
     String basename = basename("~{filter_ready_vcf_file}", ".txt")
     String? consequenceTerm
   }
   command <<<
     filter_vep -i ~{filter_ready_vcf_file} -o ~{basename}.~{consequenceTerm}.txt --filter "Consequence in ~{consequenceTerm}"
   >>>
   runtime {
     docker: "ensemblorg/ensembl-vep:latest"
     memory: "16 GB"
     cpu : "2"
     retry : "2"
   }
   output
   {
     File ConsequenceTerm_specific_variants = "~{basename}.~{consequenceTerm}.txt"
   }
}

task FilterBySiftPrediction
{
   input {
     File filter_ready_vcf_file
     String basename = basename("~{filter_ready_vcf_file}", ".txt")
     String? siftPrediction
   }
   command <<<
     filter_vep -i ~{filter_ready_vcf_file} -o ~{basename}.~{siftPrediction}.txt --filter "SIFT == ~{siftPrediction}"
   >>>
   runtime {
     docker: "ensemblorg/ensembl-vep:latest"
     memory: "16 GB"
     cpu : "2"
     retry : "2"
   }
   output
   {
     File SiftPrediction_specific_variants = "~{basename}.~{siftPrediction}.txt"
   }
}

task FilterBySiftScore
{
   input {
     File filter_ready_vcf_file
     String basename = basename("~{filter_ready_vcf_file}", ".txt")
     String? siftScore
   }
   command <<<
     filter_vep -i ~{filter_ready_vcf_file} -o ~{basename}.~{siftScore}.txt --filter "SIFT < ~{siftScore}"
   >>>
   runtime {
     docker: "ensemblorg/ensembl-vep:latest"
     memory: "16 GB"
     cpu : "2"
     retry : "2"
   }
   output
   {
     File SiftScore_specific_variants = "~{basename}.~{siftScore}.txt"
   }
}

task FilterByAlleleFrequency
{
   input {
     File filter_ready_vcf_file
     String basename = basename("~{filter_ready_vcf_file}", ".txt")
     Int? alleleFrequency
   }
   command <<<
     filter_vep -i ~{filter_ready_vcf_file} -o ~{basename}.~{alleleFrequency}.txt --filter "AF < ~{alleleFrequency}"
   >>>
   runtime {
     docker: "ensemblorg/ensembl-vep:latest"
     memory: "16 GB"
     cpu : "2"
     retry : "2"
   }
   output
   {
     File AlleleFrequency_specific_variants = "~{basename}.~{alleleFrequency}.txt"
   }
}

