cd /mnt/scratch/xreij011/bin/
for filename in /mnt/scratch/xreij011/parent/*.vcf
 do ./bcftools view -i '%QUAL>=20' -v snps "${filename}"  -Oz -o "${filename}".gz
done

for filename in /mnt/scratch/xreij011/parent/*.gz
 do ./bcftools index  "${filename}"
done

./bcftools isec -n+10 /mnt/scratch/xreij011/parent/*.gz -o /mnt/scratch/xreij011/test_parental.vcf
