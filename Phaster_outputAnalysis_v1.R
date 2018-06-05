
################################################################################################################################

##################### PHASTER OUPUT ANALYSIS
################################################################################################################################

# before in basth
#1. mkdir Prophages
#for i in $(ls Vibrio_cholerae_*.fna); do echo $i ; wget --header="Content-type: multipart/form-data boundary=FILEUPLOAD" --post-file $i http://phaster.ca/phaster_api?contigs=1 -O Prophages/$(echo $i | cut -d'_' -f3).Job ;done

#2. Extract info about job ID
#for i in $(ls [0-9]*.Job); do cat $i| cut -d',' -f1 |cut -d':' -f2 | sed 's/\"//' | sed 's/^/wget \"htp:\/\/phaster.ca\/phaster_api?acc=/';  echo $i ; done | tr '\n' ' '  | sed 's/.Job/.txt\n/g' | sed 's/\" /\" -O Phaster_/' | sed 's/^ //'

#3. copy paster to get data

#4. modify results to import in R
#for i in $(ls Phaster_*.txt); do cat $i | sed 's/\\n/@/g' | tr '@' '\n' | grep -i -A 100 'gc%:' | grep -v 'gc%:' | grep -v '\-\-' | sed 's/  / /g' | sed 's/  / /g' | sed 's/  / /g' | grep -v '\"\}' > PhasterVF_$(echo $i | cut -d'_' -f2); done 


################################################################################################################################
library("pheatmap")
library(ggplot2)
library('stringr')
library(grid)
library(msa)
library(PopGenome)
options(digits=2)
# phylogeny
library(ape)
library("ggtree")
library('dendextend')

setwd('~/Documents/EPFL/Vibrio_Chromosomal_rearrangements/PopGenomics_V2/ALL_POPS/Prophages')


filesToProcess <- dir(pattern = "PhasterVF")  #files to pro# if event 3 merged
listOfFiles <- lapply(filesToProcess, function(x) tryCatch(read.delim(x, sep="", h=T,stringsAsFactors=F),
                                                           error= function (e) cbind.data.frame(V1="NA",V2="NA",V3="NA",
                                                                                                V4="NA",V5="NA",V6="NA",
                                                                                                V7=0,V8=0,V9="NA",
                                                                                                V10="NA",V11="NA",V12="NA")))
names(listOfFiles)<- gsub("\\.txt","", gsub("PhasterVF_","", filesToProcess) )
names(listOfFiles)

# Reformat DAta

Phago<-list()
for ( i in 1:length(listOfFiles)){
  pha<-listOfFiles[[i]]
  
  Straino<-rep( names(listOfFiles[i]),dim(pha)[1])
    Posas<-lapply(strsplit(pha$REGION_POSITION,split = ":") , function (x) cbind.data.frame(Contig=x[1], Position=x[2]))
  Posas1<-lapply(strsplit(pha$REGION_POSITION,split = ":") , function (x) x[1])
    Posas2<-lapply(Posas,function (x) unlist(strsplit(as.character(x$Position),split="-")))
    df <- data.frame(matrix(unlist(Posas2), nrow=length(Posas2), byrow=T))
  df<-cbind.data.frame(cbind.data.frame(Chromosome=unlist(Posas1) ),df)
  df<-cbind.data.frame(Strain=Straino,df,Length=as.numeric(levels(df[,3]))[df[,3]] - as.numeric(levels(df[,2]))[df[,2]], State=unlist(lapply(strsplit(pha[,3], split="\\(" ), function (x) x[1])), pha[,c(4,6:13,15)],
                       phagePerc=as.numeric(gsub("%","",pha[,16]))/100, gc=as.numeric(gsub("%","",pha[,17]))/100,
                                              PhageName=unlist(lapply(strsplit(unlist(lapply(strsplit(pha[,14],split=","), function (x) x[1]) ),split="\\(") , function (x) x[1]) ) )
  
  Phago[[i]]<-df
  names(Phago[i])<-names(listOfFiles[i])

}

colnam<-c("Strain","Chr","Start","End","Length","State","Parts","tRNA","Proteins","PhageProteins","HypotheticalProteins","Phage_HypProtsPerc","BacterialProts","AttachmentSite","PhageSpecies","MostCommonPhageNb","MostCommonPhagePerc","GCPerc","MostCommonPhage" )
Phago <- lapply(Phago, setNames, nm=colnam)
Phago


PhageTotalLength<-lapply( Phago,function (x) cbind.data.frame(unique(x[1]),ProphageLength=sum(x[5]), ProphageNb=dim(x)[1]) )


PhageTotalLength<- plyr::ldply(PhageTotalLength, data.frame)
PhageTotalLength$Strain<-gsub("-","_",PhageTotalLength$Strain)

rownames(PhageTotalLength)<-PhageTotalLength$Strain
PhageTotalLength


########################################################333




# Phages plus phylogeny

# PHYLOGENY PLUS HEATMAP 


Populations<-read.csv("~/Documents/EPFL/Vibrio_Chromosomal_rearrangements/PopGenomics_V2/ALL_POPS/Phylogenies_ALL/Strains_Pop_vf.csv", h=T)
Populations$strain<-gsub("-","_",Populations$strain)

tree <- read.tree('~/Documents/EPFL/Vibrio_Chromosomal_rearrangements/PopGenomics_V2/ALL_POPS/Phylogenies_ALL/aLRT_WAG/ALL6POP_OUT_Vibrio_CoreGenes_MafftAlign.phy_phyml_tree.txt')
tree$tip.label<-gsub("_faa","", gsub("Core_","", tree$tip.label)   )

tree  <- root(tree, outgroup = "Vf_ATCC33809",resolve.root = T)
tree <- drop.tip(tree, "Vf_ATCC33809")

popColor<-data.frame(strain=tree$tip.label)
popColor<-merge(Populations,popColor,by = "strain")
popColor<-popColor[match(tree$tip.label, popColor$strain),]


# add heatmap
pdf('~/Documents/EPFL/Vibrio_Chromosomal_rearrangements/PopGenomics_V2/ALL_POPS/Prophages/Phylogeny_Prophages_info.pdf',height = 16,width=12)
P1<-ggtree(tree) +geom_treescale(fontsize=6, linesize=1, offset=-10) + 
  geom_tiplab(size=3.5, color=popColor$color) 

p20<-facet_plot(P1, panel="Prophages total Nb", data=PhageTotalLength, geom=geom_segment, aes(x=0, xend=ProphageNb, y=y, yend=y), size=3, color='gray24',alpha=3/4)+ theme_tree2()

p21<-facet_plot(p20, panel="Prophages total length", data=PhageTotalLength, geom=geom_segment, aes(x=0, xend=ProphageLength, y=y, yend=y), size=3, color='gray24',alpha=3/4)+ theme_tree2()

plot(p21)

dev.off()


Phago


#########################################################################################################################################################
table(unlist(lapply(Phago, function (x) as.character(x[,19])) ) )
barplot(table(unlist(lapply(Phago, function (x) as.character(x[,19])) ) ), las=3)

# Get more insights about phages content.

# generalist phages? specific phages?

# Gene  content in phage?  carrying function?
# Clade specific phages? 
# can  clusterize phages genomes to see wich names can be put under the same label. download phages sequences. and clusterize.


#########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################

########################################################################################################################################################################


# Lex A involved in Prophage acquisition ? 
########################################################

lexA<-read.table("~/Documents/EPFL/Vibrio_Chromosomal_rearrangements/PopGenomics_V2/ALL_POPS/INSERTIONS_VF/Chr1_InsE/4_LexA_in_ALL_POP/LexA_copies_inPOP.txt",h=F,stringsAsFactors=F)

lexA[,1]<-gsub("-","_",lexA[,1])

plot(lexA[,2],PhageTotalLength$ProphageLength)
