tagdust -arch arch.txt -t 80 ERR048988.fastq -o 15c_repeat1_lane1_extracted
tagdust -arch arch.txt -t 80 ERR048990.fastq -o 25c_repeat1_lane1_extracted
tagdust -arch arch.txt -t 80 ERR048992.fastq -o RNA_corresponding_to_1000_cells_lane1_extracted
tagdust -arch arch.txt -t 80 ERR048994.fastq -o RNA_corresponding_to_10_cells_lane1_extracted
tagdust -arch arch.txt -t 80 ERR048989.fastq -o 15c_repeat2_lane1_extracted
tagdust -arch arch.txt -t 80 ERR048991.fastq -o 25c_repeat1_lane2_extracted
tagdust -arch arch.txt -t 80 ERR048993.fastq -o RNA_corresponding_to_100_cells_lane1_extracted

for file in ./*
        do
                if [ -f $file ]; then 
                        if [[ $file =~ _extracted_BC_GACTT.fq$ ]]; then
                                if [ -f $$file.sorted.bam ]; then
                                        echo "$file already processed"
                                else
                                                echo "working on $file"
                                                sed  -e 's/[ \t]//g' $file > tmp.fq
                                                 bwa aln -t 80  $GENOME_DIR tmp.fq  | bwa samse $GENOME_DIR - tmp.fq | samtools view -Su - | samtools sort - $file.sorted

                                fi
                        fi
                fi
        done
samtools view -F 16 -q 20 15c_repeat1_lane1_extracted_BC_GACTT.fq.sorted.bam | awk '{x = split($1,a,";"); b = substr(a[2],4);print b , $3 ,$4 }' | sort -S 16G -k 2,2 -k 3,3n -k 1,1n | uniq -c > 15c_repeat1_lane1_extracted_processed.csv 

samtools view -F 16 -q 20 RNA_corresponding_to_1000_cells_lane1_extracted_BC_GACTT.fq.sorted.bam | awk '{x = split($1,a,";"); b = substr(a[2],4);print b , $3 ,$4 }'|sort -S 16G -k 2,2 -k 3,3n -k 1,1n | uniq -c > RNA_corresponding_to_1000_cells_lane1_extracted_processed.csv

 
samtools view -F 16 -q 20 15c_repeat2_lane1_extracted_BC_GACTT.fq.sorted.bam | awk '{x = split($1,a,";"); b = substr(a[2],4);print b , $3 ,$4 }'| sort -S 16G -k 2,2 -k 3,3n -k 1,1n | uniq -c > 15c_repeat2_lane1_extracted_processed.csv &

samtools view -F 16 -q 20 RNA_corresponding_to_100_cells_lane1_extracted_BC_GACTT.fq.sorted.bam | awk '{x = split($1,a,";"); b = substr(a[2],4);print b , $3 ,$4 }'| sort -S 16G -k 2,2 -k 3,3n -k 1,1n | uniq -c > RNA_corresponding_to_100_cells_lane1_extracted_processed.csv

 
samtools view -F 16 -q 20 25c_repeat1_lane1_extracted_BC_GACTT.fq.sorted.bam | awk '{x = split($1,a,";"); b = substr(a[2],4);print b , $3 ,$4 }'  | sort -S 16G -k 2,2 -k 3,3n -k 1,1n | uniq -c > 25c_repeat1_lane1_extracted_processed.csv
 
samtools view -F 16 -q 20 RNA_corresponding_to_10_cells_lane1_extracted_BC_GACTT.fq.sorted.bam | awk  '{x = split($1,a,";"); b = substr(a[2],4);print b , $3 ,$4 }' | sort -S 16G -k 2,2 -k 3,3n -k 1,1n | uniq -c > RNA_corresponding_to_10_cells_lane1_extracted_processed.csv

samtools view -F 16 -q 20 25c_repeat1_lane2_extracted_BC_GACTT.fq.sorted.bam | awk  '{x = split($1,a,";"); b = substr(a[2],4);print b , $3 ,$4 }'  | sort -S 16G -k 2,2 -k 3,3n -k 1,1n | uniq -c > 25c_repeat1_lane2_extracted_processed.csv 

cat 15c_repeat1_lane1_extracted_processed.csv | awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,$4}'  | groupBy  -g 3 -c 1 -o count > 15c_repeat1_lane1_UMI.csv 
cat 15c_repeat1_lane1_extracted_processed.csv | awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,$4}'  | groupBy  -g 3 -c 1 -o sum > 15c_repeat1_lane1_NOUMI.csv  

cat 15c_repeat2_lane1_extracted_processed.csv | awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,$4}'  | groupBy  -g 3 -c 1 -o count > 15c_repeat2_lane1_UMI.csv 
cat 15c_repeat2_lane1_extracted_processed.csv | awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,$4}'  | groupBy  -g 3 -c 1 -o sum > 15c_repeat2_lane1_NOUMI.csv 


cat 25c_repeat1_lane1_extracted_processed.csv | awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,$4}'  | groupBy  -g 3 -c 1 -o count > 25c_repeat1_lane1_UMI.csv 
cat 25c_repeat1_lane1_extracted_processed.csv | awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,$4}'  | groupBy  -g 3 -c 1 -o sum > 25c_repeat1_lane1_NOUMI.csv 

cat 25c_repeat1_lane2_extracted_processed.csv | awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,$4}'  | groupBy  -g 3 -c 1 -o count > 25c_repeat1_lane2_UMI.csv 
cat 25c_repeat1_lane2_extracted_processed.csv | awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,$4}'  | groupBy  -g 3 -c 1 -o sum > 25c_repeat1_lane2_NOUMI.csv  


cat RNA_corresponding_to_1000_cells_lane1_extracted_processed.csv | awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,$4}'  | groupBy  -g 3 -c 1 -o count > RNA_corresponding_to_1000_cells_UMI.csv  
cat RNA_corresponding_to_1000_cells_lane1_extracted_processed.csv | awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,$4}'  | groupBy  -g 3 -c 1 -o sum > RNA_corresponding_to_1000_cells_NOUMI.csv 

cat RNA_corresponding_to_100_cells_lane1_extracted_processed.csv | awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,$4}'  | groupBy  -g 3 -c 1 -o count > RNA_corresponding_to_100_cells_UMI.csv 
cat RNA_corresponding_to_100_cells_lane1_extracted_processed.csv | awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,$4}'  | groupBy  -g 3 -c 1 -o sum > RNA_corresponding_to_100_cells_NOUMI.csv 

cat RNA_corresponding_to_10_cells_lane1_extracted_processed.csv | awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,$4}'  | groupBy  -g 3 -c 1 -o count > RNA_corresponding_to_10_cells_UMI.csv 
cat RNA_corresponding_to_10_cells_lane1_extracted_processed.csv | awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,$4}'  | groupBy  -g 3 -c 1 -o sum > RNA_corresponding_to_10_cells_NOUMI.csv 

cat 15c_repeat1_lane1_NOUMI.csv 15c_repeat2_lane1_NOUMI.csv | sort -k 1,1  | groupBy  -g 1 -c 2 -o sum > 15c_combined_noumi.csv 
cat 15c_repeat1_lane1_UMI.csv  15c_repeat2_lane1_UMI.csv  | sort -k 1,1  | groupBy  -g 1 -c 2 -o sum > 15c_combined_UMI.csv 

cat 25c_repeat1_lane1_NOUMI.csv  25c_repeat1_lane2_NOUMI.csv | sort -k 1,1  |  groupBy  -g 1 -c 2 -o sum > 25c_combined_noumi.csv 
cat 25c_repeat1_lane1_UMI.csv   25c_repeat1_lane2_UMI.csv | sort -k 1,1  |  groupBy  -g 1 -c 2 -o sum > 25c_combined_UMI.csv

join -a1 -a2 -1 1 -2 1 -o 0 1.2 2.2 -e "0" 15c_repeat1_lane1_UMI.csv  15c_repeat2_lane1_UMI.csv | awk '{printf "%s\t%d\t%d\n",$1,$2,$3 }' > comp_15c_umi.csv

join -a1 -a2 -1 1 -2 1 -o 0 1.2 2.2 -e "0" 25c_repeat1_lane1_NOUMI.csv  25c_repeat1_lane2_NOUMI.csv | awk '{printf "%s\t%d\t%d\n",$1,$2,$3 }' > comp_25c_noumi.csv

join -a1 -a2 -1 1 -2 1 -o 0 1.2 2.2 -e "0" 25c_repeat1_lane1_UMI.csv  25c_repeat1_lane2_UMI.csv | awk '{printf "%s\t%d\t%d\n",$1,$2,$3 }' > comp_25c_umi.csv

join -a1 -a2 -1 1 -2 1 -o 0 1.2 2.2 -e "0" 15c_repeat1_lane1_NOUMI.csv  25c_repeat1_lane1_NOUMI.csv | awk '{printf "%s\t%d\t%d\n",$1,$2,$3 }' > comp_15c_repeat1_25c_repeat1_noumi.csv

join -a1 -a2 -1 1 -2 1 -o 0 1.2 2.2 -e "0" 15c_repeat1_lane1_UMI.csv  25c_repeat1_lane1_UMI.csv | awk '{printf "%s\t%d\t%d\n",$1,$2,$3 }' > comp_15c_repeat1_25c_repeat1_umi.csv

join -a1 -a2 -1 1 -2 1 -o 0 1.2 2.2 -e "0" 15c_repeat2_lane1_NOUMI.csv  25c_repeat1_lane2_NOUMI.csv | awk '{printf "%s\t%d\t%d\n",$1,$2,$3 }' > comp_15c_repeat2_25c_repeat1_lane2_noumi.csv

join -a1 -a2 -1 1 -2 1 -o 0 1.2 2.2 -e "0" 15c_repeat2_lane1_UMI.csv  25c_repeat1_lane2_UMI.csv | awk '{printf "%s\t%d\t%d\n",$1,$2,$3 }' > comp_15c_repeat2_25c_repeat1_lane2_umi.csv

join -a1 -a2 -1 1 -2 1 -o 0 1.2 2.2 -e "0" 15c_repeat1_lane1_NOUMI.csv  25c_repeat1_lane2_NOUMI.csv | awk '{printf "%s\t%d\t%d\n",$1,$2,$3 }' > comp_15c_repeat1_25c_repeat1_lane2_noumi.csv

join -a1 -a2 -1 1 -2 1 -o 0 1.2 2.2 -e "0" 15c_repeat1_lane1_UMI.csv  25c_repeat1_lane2_UMI.csv | awk '{printf "%s\t%d\t%d\n",$1,$2,$3 }' > comp_15c_repeat1_25c_repeat1_lane2_umi.csv

join -a1 -a2 -1 1 -2 1 -o 0 1.2 2.2 -e "0" 15c_repeat2_lane1_NOUMI.csv  25c_repeat1_lane1_NOUMI.csv | awk '{printf "%s\t%d\t%d\n",$1,$2,$3 }' > comp_15c_repeat2_25c_repeat1_lane1_noumi.csv

join -a1 -a2 -1 1 -2 1 -o 0 1.2 2.2 -e "0" 15c_repeat2_lane1_UMI.csv  25c_repeat1_lane1_UMI.csv | awk '{printf "%s\t%d\t%d\n",$1,$2,$3 }' > comp_15c_repeat2_25c_repeat1_lane1_umi.csv

join -a1 -a2 -1 1 -2 1 -o 0 1.2 2.2 -e "0"  15c_combined_noumi.csv 25c_combined_noumi.csv | awk '{printf "%s\t%d\t%d\n",$1,$2,$3 }' > 15_25_combined_noumi.csv
join -a1 -a2 -1 1 -2 1 -o 0 1.2 2.2 -e "0"  15c_combined_UMI.csv 25c_combined_UMI.csv | awk '{printf "%s\t%d\t%d\n",$1,$2,$3 }' > 15_25_combined_UMI.csv
