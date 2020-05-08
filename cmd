find data/libraries/ -name *.fastq* -exec fastqc --noextract -q -t 12 {} \;
for i in * ; do [[ $i =~ RNA-m ]] && mv $i/*R2* ${i/RNA-m/RNA-adu-m}/${i/RNA-m/RNA-adu-m}_R2.fastq.gz ; done
for i in * ; do [[ $i =~ RNA-adu-f ]] && mv $i/*R2* ${i}/${i}_R2.fastq.gz ; done 
for i in * ; do [[ $i =~ RNA-m ]] && mv $i ${i/RNA-m/RNA-adu-m} ; done
