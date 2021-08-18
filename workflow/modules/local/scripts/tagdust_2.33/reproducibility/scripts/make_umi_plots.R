library(ggplot2)

mat = read.table("15_25_combined_noumi.csv",header = F);
cor = cor(mat[,2],mat[,3]);
frame=ggplot(mat,aes(mat[,2],mat[,3]))
frame + geom_point(colour="black",size=1,shape=20,alpha=0.3) + labs(x="15 cycles ",y="25 cycles",title="without UMI") + annotate("text", label =paste("R=", round (cor,6), sep=""), x = min(mat[,2])+10, y = max(mat[,3]) / 10, size = 6, colour = "red") +scale_x_log10() + scale_y_log10()+theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"),title=element_text(size=18,face="bold"));
ggsave("15_25_combined_noumi.pdf")
mat = read.table("15_25_combined_UMI.csv",header = F);

cor = cor(mat[,2],mat[,3]);
frame=ggplot(mat,aes(mat[,2],mat[,3]))
frame + geom_point(colour="black",size=1,shape=20,alpha=0.3) + labs(x="15 cycles",y="25 cycles",title="with UMI") + annotate("text", label =paste("R=", round (cor,6), sep=""), x = min(mat[,2])+10, y = max(mat[,3]) / 10, size = 6, colour = "red") +scale_x_log10() + scale_y_log10()+theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"),title=element_text(size=18,face="bold"));
ggsave("15_25_combined_UMI.pdf")
