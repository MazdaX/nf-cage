
.PHONY:  message benchmark rplot

benchmark: Figure1.pdf Figure2.pdf
	@echo Done

Figure1.pdf: Figure1.tex
	pdflatex $<;

Figure2.pdf: barread_6nt_4r.tsv barread_4nt_4r.tsv 5barread3_4nt_4r.tsv 5barread3_6nt_4r.tsv
	R --slave --vanilla < plotting.R
	rm Rplots.pdf


barread_4nt_4r.tsv: barread_results.txt
	cat barread_results.txt | grep _4nt_ | awk 'BEGIN{print "Program\tsimerror\tbarcodes\tRecall\tSpecificity\tPrecision\tKappa"} {if(NR>0){x = split($$1,a,"_");printf "%s\t%f\t%f\t%f\t%f\t%f\t%f\n" ,  a[1],a[2],a[3],$$2,$$3,$$4,$$5 }}' > barread_4nt_4r.tsv
	
5barread3_4nt_4r.tsv: 5barread3_results.txt
	cat 5barread3_results.txt |  grep _4nt_ | awk 'BEGIN{print "Program\tsimerror\tbarcodes\tRecall\tSpecificity\tPrecision\tKappa"} {if(NR>0){x = split($$1,a,"_");printf "%s\t%f\t%f\t%f\t%f\t%f\t%f\n" ,  a[1],a[2],a[3],$$2,$$3,$$4,$$5 }}' > 5barread3_4nt_4r.tsv

barread_6nt_4r.tsv: barread_results.txt
	cat barread_results.txt | grep _6nt_ | awk 'BEGIN{print "Program\tsimerror\tbarcodes\tRecall\tSpecificity\tPrecision\tKappa"} {if(NR>0){x = split($$1,a,"_");printf "%s\t%f\t%f\t%f\t%f\t%f\t%f\n" ,  a[1],a[2],a[3],$$2,$$3,$$4,$$5 }}' > barread_6nt_4r.tsv
	
5barread3_6nt_4r.tsv: 5barread3_results.txt
	cat 5barread3_results.txt |  grep _6nt_ | awk 'BEGIN{print "Program\tsimerror\tbarcodes\tRecall\tSpecificity\tPrecision\tKappa"} {if(NR>0){x = split($$1,a,"_");printf "%s\t%f\t%f\t%f\t%f\t%f\t%f\n" ,  a[1],a[2],a[3],$$2,$$3,$$4,$$5 }}' > 5barread3_6nt_4r.tsv

barread_results.txt: all_arch.txt
	./barread.sh -b ../../dev/EDITTAG_4nt_ed_2.txt -r
	./barread.sh -b ../../dev/EDITTAG_6nt_ed_3.txt -r
	
5barread3_results.txt: all_arch.txt
	./5barread3.sh -b ../../dev/EDITTAG_4nt_ed_2.txt -r
	./5barread3.sh -b ../../dev/EDITTAG_6nt_ed_3.txt -r

all_arch.txt: 
	./barread.sh -b ../../dev/EDITTAG_4nt_ed_2.txt
	./barread.sh -b ../../dev/EDITTAG_6nt_ed_3.txt
	./5barread3.sh -b ../../dev/EDITTAG_4nt_ed_2.txt
	./5barread3.sh -b ../../dev/EDITTAG_6nt_ed_3.txt
	cat  *tagdust_arch.txt  | sort | uniq  > all_arch.txt
	
all: message

message:
	@echo To reproduce the figures run "make -f run.mk benchmark"
