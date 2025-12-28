library(EthSEQ)

# Directory of the VCF file to be analyzed
vcf_file <- "/storage2/jpitta/genomas_luydson/shareRodolfo/all_gatk_final.vcf.gz" 

# Definition of the desired EthSEQ model, you can list with getModelsList() and choose according to your interest
model_available <- "Gencode.Exome"
model_assembly <- "hg38"          
model_pop <- "All"                

# Output directory for the results
out_dir <- "EthSEQ_Results_3D"

# Running the analysis with the VCF as input, passing the compressed VCF file in target.vcf
ethseq.Analysis(
  target.vcf = vcf_file,
  model.available = model_available,
  model.assembly = model_assembly,
  model.pop = model_pop,
  out.dir = out_dir,
  verbose = TRUE,
  cores = 2,
  composite.model.call.rate = 1,
  space = "3D"
)
