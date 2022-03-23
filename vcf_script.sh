cd /mnt/scratch/xreij011/
for filename in /mnt/scratch/xreij011/parent/*.vcf ; do grep -v "${filename}" | awk '$6 > 20' > "${filename}" 
done
cd /mnt/scratch/xreij011/bin/
for filename in /mnt/scratch/xreij011/parent/*.vcf; do ./bcftools view "${filename}"  -Oz -o "${filename}".gz; done
for filename in /mnt/scratch/xreij011/parent/*.gz; do ./bcftools index  "${filename}"; done
./bcftools isec -n+10 /mnt/scratch/xreij011/parent/*.gz -o /mnt/scratch/xreij011/test_parental.vcf
