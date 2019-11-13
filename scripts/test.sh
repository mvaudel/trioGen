

java -Xmx16G -cp bin/triogen-0.1.0/triogen-0.1.0.jar no.uib.triogen.cmd.transmission.ExtractTransmission -g /mnt/archive/MOBAGENETICS/genotypes-base/imputed/all/vcf/22.vcf.gz -gf 1 -f /mnt/archive/helgeland/helgeland2020/sensitive/fulltriads-n24280-180919 -o /mnt/work/marc/moba/trioGen/tmp/test_h -test


java -Xmx16G -cp bin/triogen-0.1.0/triogen-0.1.0.jar no.uib.triogen.cmd.association.LinearModel -g /mnt/archive/MOBAGENETICS/genotypes-base/imputed/all/vcf/22.vcf.gz -gf 1 -f /mnt/archive/helgeland/helgeland2020/sensitive/fulltriads-n24280-180919 -p src/main/resources/transmission/phenos_linear_model.txt -pn pheno1,pheno2,pheno3,pheno4 -o /mnt/work/marc/moba/trioGen/tmp/test_lm.gz -test
