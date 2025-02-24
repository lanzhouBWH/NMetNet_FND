####### Loading Required Libraries #############################################
library(irrCAC)
library(readxl)
library(openxlsx2)
library(exact2x2)
library(survival)
library(ggplot2)
library(ggfortify)
#library(meta)
library(tidyr)
library(pwr)
library(miceadds)
#library(stringr)
library(lme4)
library(lmerTest)
#library(emmeans)
library(rlist)
#library(abglasso)
library(igraph)
#library(stringr)
library(visNetwork)
library(abind)
################################################################################

Base="/Users/zhoulan/Partners HealthCare Dropbox/Zhou Lan/NeuroClinical_Review/Revisions/"### Please add you home direcory here

##### Defining Directories #####################################################
Input_Path=paste0(Base,"/NMetNet_FND/Input/")
Output_Path=paste0(Base,"NMetNet_FND/Output/")
################################################################################



##### General Information ######################################################
Regions=c("aDMN","pDMN","SMA")
MMM=c("NAA", "tCho" ,"mI", "Glx", "GSH","Cr")
Dict=expand.grid(MMM,Regions)
Dict$LL=as.vector(sapply(Regions, function(rr) paste(MMM,rr,sep=" at ")))
colnames(Dict)<-c("Metabolite","Region","Label")
################################################################################

##### Loading Results ##########################################################
load(paste0(Input_Path,"MCMC_FND-Others_vs_FND-Seizure.Rdata"))


Graph_Creation_Fud2<-function(Post,theta=.1,Dict){
  
  Niter=length(Post)
  Network_list=lapply(1:Niter, function(it) cov2cor(Post[[it]]))
  Network_cc=abind(Network_list, along = 3)
  Network=matrix(0,ncol(Post[[1]]),nrow(Post[[1]]))
  Network_weight=matrix(0,ncol(Post[[1]]),nrow(Post[[1]]))
  for(i in 1:(nrow(Post[[1]])-1) ){
    for(j in (i+1):nrow(Post[[1]]) ){
      Network[i,j]=mean(abs(Network_cc[i,j,])>0.001)*sign(mean(Network_cc[i,j,]))
      Network_weight[i,j]=(mean(Network_cc[i,j,]))
    }
  }
  diag(Network)<-0
  Network=Network+t(Network)
  
  abs=abs(Network)
  
  off_diag <- !diag(nrow(abs))
  off_diag_elements <- abs[off_diag]
  
  Network_weight[which(abs<=0.5)]=0
  Network_weight=Network_weight*Network
  Network<-Network_weight
  
  Network_graph=graph.adjacency(Network, weighted = TRUE,mode = "undirected")
  V(Network_graph)$ID<-1:nrow(Network)
  V(Network_graph)$Metabolite<-as.character(Dict$Metabolite)
  V(Network_graph)$Region<-as.character(Dict$Region)
  V(Network_graph)$Label<-Dict$Label
  
  deg <- which(degree(Network_graph, mode="all")>=1)
  sub=subgraph(Network_graph, deg)
  
  
  return(sub)
}



Post1=lapply(1:10000,function(i) Result1$Omega[,,i])
Post2=lapply(1:10000,function(i) Result2$Omega[,,i])
Graph1=Graph_Creation_Fud2(Post1,theta=.1,Dict)
Graph2=Graph_Creation_Fud2(Post2,theta=.1,Dict)

V(Graph1)$Type="FND-Seizure"
V(Graph2)$Type="FND-Others"

ff1=cluster_louvain(Graph1, weights = abs(E(Graph1)$weight))
plot(ff1,Graph1,vertex.label=V(Graph1)$Label)

ff2=cluster_louvain(Graph2, weights = abs(E(Graph2)$weight))
plot(ff2,Graph2,vertex.label=V(Graph2)$Label)



##### FND-Seizure #######

perm=unlist(sapply(c("aDMN","SMA","pDMN"),function(v) which(c(V(Graph1)$Region)==v)))
Graph1 <- permute(Graph1, perm)

nodes=data.frame(id=1:(length(V(Graph1)$Type)),
                 label=c(V(Graph1)$Label),
                 Metabolite=c(V(Graph1)$Metabolite),
                 Region=c(V(Graph1)$Region),
                 Type=c(V(Graph1)$Type),
                 Community=c(ff1$membership)
)
edges=as.data.frame(
  rbind(
    as_edgelist(Graph1)
  )
)

colnames(edges)<-c("from","to")
edges$weight=c(sign(E(Graph1)$weight))

nodes$shape <-c("triangle")
nodes$color <- factor(nodes$Region,labels = c("#ffb09c","#96F97A","lightskyblue"))
edges$color <- ifelse(edges$weight == 1, "red", "blue")
nodes$label<-nodes$Metabolite
edges$width_ori<-(c(abs(E(Graph1)$weight)))
edges$width<-(c(abs(E(Graph1)$weight)))*20+5

nodes$x=NA
nodes$y=1+rnorm(nrow(nodes), sd = 0.2)
nodes$y[nodes$Region=="SMA"]=nodes$y[nodes$Region=="SMA"]-0.6

for(rr in c("aDMN","SMA","pDMN")){
  nodes$x[nodes$Region==rr]=which(nodes$Region==rr)+ rnorm(sum(nodes$Region==rr), sd = 0.2)
}
nodes$label[nodes$label=="Ins"]="mI"

nodes$font.vadjust=-80
visNetwork(nodes, edges) %>%
  visNodes(font = list(size = 60, bold = TRUE))%>%
  visIgraphLayout(layoutMatrix =cbind(nodes$x,nodes$y) )
#%>%
 # visOptions(selectedBy = c("Community"))

Connection_Table=data.frame(Metabolite1=NULL,Region1=NULL,Metabolite2=NULL,Region2=NULL,Partial_Correlation=NULL)

for(i in 1:nrow(edges)){
  temp=data.frame(Metabolite1=nodes$Metabolite[edges$from[i]]
                  ,Region1=nodes$Region[edges$from[i]],
                  Metabolite2=nodes$Metabolite[edges$to[i]],
                  Region2=nodes$Region[edges$to[i]],Partial_Correlation=(edges$width_ori[i]))
  Connection_Table=rbind(Connection_Table,temp)
}
write_xlsx(Connection_Table,file=paste0(Output_Path,"ConTable_FND_Seizures.xlsx"))

##### FND-Others #######

perm=unlist(sapply(c("aDMN","SMA","pDMN"),function(v) which(c(V(Graph2)$Region)==v)))
Graph2 <- permute(Graph2, perm)

nodes=data.frame(id=1:(length(V(Graph2)$Type)),
                 label=c(V(Graph2)$Label),
                 Metabolite=c(V(Graph2)$Metabolite),
                 Region=c(V(Graph2)$Region),
                 Type=c(V(Graph2)$Type),
                 Community=c(ff2$membership)
)
edges=as.data.frame(
  rbind(
    as_edgelist(Graph2)
  )
)

colnames(edges)<-c("from","to")
edges$weight=c(sign(E(Graph2)$weight))

nodes$shape <- ifelse(nodes$Type == "Health", "square", "triangle")
nodes$color <- factor(nodes$Region,labels = c("red","green","lightskyblue") )
edges$color <- ifelse(edges$weight == 1, "red", "blue")
nodes$label<-nodes$Metabolite
edges$width_ori<-(c(abs(E(Graph2)$weight)))
edges$width<-(c(abs(E(Graph2)$weight)))*20+5

nodes$x=NA
nodes$y=1+rnorm(nrow(nodes), sd = 0.2)
nodes$y[nodes$Region=="SMA"]=nodes$y[nodes$Region=="SMA"]-0.6

for(rr in c("aDMN","SMA","pDMN")){
  nodes$x[nodes$Region==rr]=which(nodes$Region==rr)+ rnorm(sum(nodes$Region==rr), sd = 0.2)
}
nodes$label[nodes$label=="Ins"]="mI"
nodes$font.vadjust=-80
visNetwork(nodes, edges) %>%
  visNodes(font = list(size = 60, bold = TRUE))%>%
  visIgraphLayout(layoutMatrix =cbind(nodes$x,nodes$y) )
#%>%
#  visOptions(selectedBy = c("Community"))

Connection_Table=data.frame(Metabolite1=NULL,Region1=NULL,Metabolite2=NULL,Region2=NULL,Partial_Correlation=NULL)

for(i in 1:nrow(edges)){
  temp=data.frame(Metabolite1=nodes$Metabolite[edges$from[i]]
                  ,Region1=nodes$Region[edges$from[i]],
                  Metabolite2=nodes$Metabolite[edges$to[i]],
                  Region2=nodes$Region[edges$to[i]],Partial_Correlation=edges$width_ori[i] )
  Connection_Table=rbind(Connection_Table,temp)
}
write_xlsx(Connection_Table,file=paste0(Output_Path,"ConTable_FND_Others.xlsx"))








