import pandas as pd

vcf = pd.read_csv("test_parental.vcf",sep="\t",header=None)
vcf.columns = ['chr','pos','ref','alt','allele']

SNPs = []

for row in vcf.itertuples():
    chr = row[1]
    pos = row[2]
    ref = row[3]
    alt = row[4]
    
    if len(ref) != 1 or len(alt) != 1:
        continue
        
    alleles = row[5]
    
    ref_alt = [ref,alt]
    
    tmp = [alleles[start:start+3] for start in range(0, len(alleles), 3)]
    bay = "".join([allele for i,allele in enumerate(tmp) if i%2 == 0 ])
    sha = "".join([allele for i,allele in enumerate(tmp) if i%2 == 1 ])
    
    bay = sum(list(map(int, bay)))
    sha = sum(list(map(int, sha)))
    
    if bay > 2 and bay < 10:
        continue
        
    if sha > 2 and sha < 10:
        continue
    
    bay_allele = bay//10
    sha_allele = sha//10
    
    if bay_allele == sha_allele:
        continue
        
    SNPs.append([chr,pos,ref_alt[bay_allele],ref_alt[sha_allele]])
    
SNP_df = pd.DataFrame(SNPs,columns=['chr','pos','bay','sha'])

SNP_df.to_csv("bay_sha_SNPs.csv",index=False)
