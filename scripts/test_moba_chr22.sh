

java -Xmx16G -cp bin/triogen-0.1.0/triogen-0.1.0.jar no.uib.triogen.cmd.transmission.ExtractTransmission -g /mnt/archive/MOBAGENETICS/genotypes-base/imputed/all/vcf/22.vcf.gz -gf 1 -f /mnt/archive/helgeland/helgeland2020/sensitive/fulltriads-n24280-180919 -o /mnt/work/marc/moba/trioGen/tmp/test_h -test


java -Xmx16G -cp bin/triogen-0.1.0/triogen-0.1.0.jar no.uib.triogen.cmd.association.LinearModel -g /mnt/archive/MOBAGENETICS/genotypes-base/imputed/all/vcf/22.vcf.gz -gf 1 -f /mnt/archive/helgeland/helgeland2020/sensitive/fulltriads-n24280-180919 -p tmp/lw.gz -pn length0,length1,length2,length3,length4,length5,length6,length7,length8,length9,length10,length11,weight0,weight1,weight2,weight3,weight4,weight5,weight6,weight7,weight8,weight9,weight10,weight11 -o /mnt/work/marc/moba/trioGen/tmp/test_lm.gz -test