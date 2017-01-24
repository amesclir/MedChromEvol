# Esta función toma los resultados de Chromevol y los resume en una tabla
# Se necesita el arbol original (Nexus tree) y los archivos *_mlAncestors.tree y *_expectations.txt 
# del mejor modelo estimado de ChromEvol
# La estructura de los archivos debe ser como estan en Drive. En la carpeta principal el tree
# en la carpeta OUT/1_Best el resto de archivos
# Tiene dos argumentos name: Nombre general del archivo (género) y path: la dirección 
# a la carpeta con los datos



#ChromR(name="Reseda",path="/home/fbalao/Datos/R/Rpackages/ChromTT/Reseda", best="1_Best_CONST_RATE_DEMI_EST")

ChromR<- function(name, path, bestmodel){
#Cargar paquetes necesarios. si no los tienes instalar con install.packages
  library(phytools)
  library(data.table)
  library(plyr)
#Fija el directorio
setwd(path)
#lee el arbol original para obtener las edades de los nodos
tree<-paste(name,".tree", sep="")
arbol<-read.nexus(tree)
arbol$tip.label<-gsub("\\'|\\]", "", arbol$tip.label)
ntip<-length(arbol$tip.label)
nnodes<-length(node.height(arbol))
ninodes<-nnodes-ntip
nodenames<-paste("N", 1:ninodes, sep="")
nodeage<- branching.times(arbol)


#Obtener la edad del stem asociado a cada nodo
stem1<-arbol$edge
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
agestem<-cbind(as.numeric.factor(revalue(factor(stem1[,1]), nodeage)), stem1[,2])
rootnode<-ntip+1
rootage<-nodeage[as.character(rootnode)]
agestem<-rbind(agestem, c(rootage,rootnode))
agestem<-agestem[order(agestem[,2]),]

#
names(nodeage)<-nodenames
#plot(arbol)
#nodelabels(nodenames)
time<-c(rep(0,ntip),nodeage)

#Obtener los números de cromosomas en cada nodo/tip *_mlAncestors.tree
path2<- paste(path,"/OUT/", bestmodel,"/", name, "_mlAncestors.tree", sep="")
chromtree<-readLines(path2)
chromtree2<-gsub("\\[|\\]", "", chromtree)
chromnumtree<-read.tree(text=chromtree2)
chromnum<-data.frame(name=c(chromnumtree$tip.label, chromnumtree$node.label))
chromnum$chrom<-lapply(strsplit(as.character(chromnum$name), "\\-"), "[", 2)
chromnum$name<-lapply(strsplit(as.character(chromnum$name), "\\-"), "[", 1)

chromtable<-data.frame(name=c(arbol$tip.label, nodenames), age=time, agestem=agestem[,1])

chromtable2<-merge(chromtable, chromnum, by="name")
chromtable2<-chromtable2[with(chromtable2, order(-age)), ]

#scatter.smooth(-chromtable2$age, chromtable2$chrom, pch=16, lwd=2, col=1)

# Leer la tabla de mutaciones *_expectations.txt con paquete data.table
path3<- paste(path,"/OUT/", bestmodel,"/", name, "_expectations.txt", sep="")
mutations <- fread(path3, 
                   skip="#ALL EVENTS EXPECTATIONS PER NODE", nrows=nnodes-1,  header=T) # No importa el error
colnames(mutations)[1]<- "name"
D1<-data.frame("N1",0,0,0,0,0)
colnames(D1)<- colnames(mutations)
mutations<-rbind(D1, mutations, fill=T)
mutations<-mutations[with(mutations, order(name)), ]


chromtable3<-merge(chromtable2, mutations, by="name")

#Redondeo de las mutaciones
round2 = function(x, n) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5
  z = trunc(z)
  z = z/10^n
  z*posneg
}


chromtable3[,5:dim(chromtable3)[2]]<-round2(chromtable3[,5:dim(chromtable3)[2]],0)
chromtable3

#Exportar la tabla a archivo tabulado
name2<-paste(name, "_summary.txt", sep="")
write.table(as.matrix(chromtable3), file=name2, sep="\t", fileEncoding = "utf8", row.names = F)
}


#plot(-chromtable2$age, chromtable2$chrom, type="p")
#cbind(-chromtable2$age, as.numeric(chromtable2$chrom))

#Chi2 test for the stem node
chischrom.stem<-function(file){
x<-read.table(file, header=T, strip.white = T)
preMed<-x[x$agestem>3.4,]
Med<-x[x$agestem<3.4,]

preMedt<-table(rowSums(preMed[,5:9], na.rm = T)==0)
Medt<-table(rowSums(Med[,5:9], na.rm = T)==0)    
chidata<-cbind(preMedt, Medt)
print(chisq.test(chidata))
}

#Chi2 test for the crown node
chischrom.crown<-function(file){
  x<-read.table(file, header=T, strip.white = T)
  preMed<-x[x$age>3.4,]
  Med<-x[x$age<3.4,]
  
  preMedt<-table(rowSums(preMed[,5:9], na.rm = T)==0)
  Medt<-table(rowSums(Med[,5:9], na.rm = T)==0)    
  chidata<-cbind(preMedt, Medt)
  print(chisq.test(chidata))
}

#Chi2 test for the half way
chischrom.halfway<-function(file){
  x<-read.table(file, header=T, strip.white = T)
  preMed<-x[((x$agestem - x$age) + x$age) >3.4,]
  Med<-x[((x$agestem - x$age) + x$age) <3.4,]
  
  preMedt<-table(rowSums(preMed[,5:9], na.rm = T)==0)
  Medt<-table(rowSums(Med[,5:9], na.rm = T)==0)    
  chidata<-cbind(preMedt, Medt)
  print(chisq.test(chidata))
}

#Chi2 test for the stem node and disploidy
chischrom.stem_disploidy<-function(file){
  x<-read.table(file, header=T, strip.white = T)
  preMed<-x[x$agestem>3.4,]
  Med<-x[x$agestem<3.4,]
  
  preMedt<-table(rowSums(preMed[,5:6], na.rm = T)==0)
  Medt<-table(rowSums(Med[,5:6], na.rm = T)==0)    
  chidata<-cbind(preMedt, Medt)
  print(chisq.test(chidata))
}

#Chi2 test for the crown node and disploidy
chischrom.crown_disploidy<-function(file){
  x<-read.table(file, header=T, strip.white = T)
  preMed<-x[x$age>3.4,]
  Med<-x[x$age<3.4,]
  
  preMedt<-table(rowSums(preMed[,5:6], na.rm = T)==0)
  Medt<-table(rowSums(Med[,5:6], na.rm = T)==0)    
  chidata<-cbind(preMedt, Medt)
  print(chisq.test(chidata))
}

#Chi2 test for the halfway and disploidy
chischrom.halfway_disploidy<-function(file){
  x<-read.table(file, header=T, strip.white = T)
  preMed<-x[((x$agestem - x$age) + x$age)>3.4,]
  Med<-x[((x$agestem - x$age) + x$age)<3.4,]
  
  preMedt<-table(rowSums(preMed[,5:6], na.rm = T)==0)
  Medt<-table(rowSums(Med[,5:6], na.rm = T)==0)    
  chidata<-cbind(preMedt, Medt)
  print(chisq.test(chidata))
}

#Chi2 test for the stem node and polyploidy
chischrom.stem_polyploidy<-function(file){
  x<-read.table(file, header=T, strip.white = T)
  preMed<-x[x$agestem>3.4,]
  Med<-x[x$agestem<3.4,]
  
  preMedt<-table(rowSums(preMed[,7:9], na.rm = T)==0)
  Medt<-table(rowSums(Med[,7:9], na.rm = T)==0)    
  chidata<-cbind(preMedt, Medt)
  print(chisq.test(chidata))
}

#Chi2 test for the crown node and polyploidy
chischrom.crown_polyploidy<-function(file){
  x<-read.table(file, header=T, strip.white = T)
  preMed<-x[x$age>3.4,]
  Med<-x[x$age<3.4,]
  
  preMedt<-table(rowSums(preMed[,7:9], na.rm = T)==0)
  Medt<-table(rowSums(Med[,7:9], na.rm = T)==0)    
  chidata<-cbind(preMedt, Medt)
  print(chisq.test(chidata))
}


#Chi2 test for the halfway and polyploidy
chischrom.halfway_polyploidy<-function(file){
  x<-read.table(file, header=T, strip.white = T)
  preMed<-x[((x$agestem - x$age) + x$age)>3.4,]
  Med<-x[((x$agestem - x$age) + x$age)<3.4,]
  
  preMedt<-table(rowSums(preMed[,7:9], na.rm = T)==0)
  Medt<-table(rowSums(Med[,7:9], na.rm = T)==0)    
  chidata<-cbind(preMedt, Medt)
  print(chisq.test(chidata))
}


#Loop para correr el test en todos los archivo y hacer una tabla STEM NODE
getwd()
setwd("/home/amesclir/Documents/Ubuntufiles/GITHUB/MedChromEvol/datasets/")
listfile<-dir()
results<-matrix(nrow=length(listfile), ncol=2)
for (i in 1:length(listfile)){
  res<-chischrom.stem(file=listfile[i])
  results.stem[i,1]<-res[[1]]
  results.stem[i,2]<-res[[3]]
}


results.stem<-as.data.frame(results.stem)
rownames(results.stem)<- listfile
colnames(results.stem)<-c("Chisq","p-value")
results.stem

#Loop para correr el test en todos los archivo y hacer una tabla CROWN NODE
results.crown<-matrix(nrow=length(listfile), ncol=2)
for (i in 1:length(listfile)){
  res<-chischrom.crown(file=listfile[i])
  results.crown[i,1]<-res[[1]]
  results.crown[i,2]<-res[[3]]
}


results.crown<-as.data.frame(results.crown)
rownames(results.crown)<- listfile
colnames(results.crown)<-c("Chisq","p-value")
results.crown

#Loop para correr el test en todos los archivo y hacer una tabla HALFWAY
results.halfway<-matrix(nrow=length(listfile), ncol=2)
for (i in 1:length(listfile)){
  res<-chischrom.halfway(file=listfile[i])
  results.halfway[i,1]<-res[[1]]
  results.halfway[i,2]<-res[[3]]
}


results.halfway<-as.data.frame(results.halfway)
rownames(results.halfway)<- listfile
colnames(results.halfway)<-c("Chisq","p-value")
results.halfway

#Loop para correr el test en todos los archivo y hacer una tabla DATOS con disploidía STEM NODE
resultsdisplo.stem<-matrix(nrow=length(listfile), ncol=2)
for (i in 1:length(listfile)){
  res<-chischrom.stem_disploidy(file=listfile[i])
  resultsdisplo.stem[i,1]<-res[[1]]
  resultsdisplo.stem[i,2]<-res[[3]]
}

resultsdisplo.stem<-as.data.frame(resultsdisplo.stem)
rownames(resultsdisplo.stem)<- listfile
colnames(resultsdisplo.stem)<-c("Chisq","p-value")
resultsdisplo.stem

#Loop para correr el test en todos los archivo y hacer una tabla DATOS con disploidía CROWN NODE
resultsdisplo.crown<-matrix(nrow=length(listfile), ncol=2)
for (i in 1:length(listfile)){
  res<-chischrom.crown_disploidy(file=listfile[i])
  resultsdisplo.crown[i,1]<-res[[1]]
  resultsdisplo.crown[i,2]<-res[[3]]
}

resultsdisplo.crown<-as.data.frame(resultsdisplo.crown)
rownames(resultsdisplo.crown)<- listfile
colnames(resultsdisplo.crown)<-c("Chisq","p-value")
resultsdisplo.crown

#Loop para correr el test en todos los archivo y hacer una tabla DATOS con disploidía HALF WAY
resultsdisplo.halfway<-matrix(nrow=length(listfile), ncol=2)
for (i in 1:length(listfile)){
  res<-chischrom.halfway_disploidy(file=listfile[i])
  resultsdisplo.halfway[i,1]<-res[[1]]
  resultsdisplo.halfway[i,2]<-res[[3]]
}

resultsdisplo.halfway<-as.data.frame(resultsdisplo.halfway)
rownames(resultsdisplo.halfway)<- listfile
colnames(resultsdisplo.halfway)<-c("Chisq","p-value")
resultsdisplo.halfway


#Loop para correr el test en todos los archivo y hacer una tabla DATOS con poliploidía STEM NODE
resultspoly.stem<-matrix(nrow=length(listfile), ncol=2)
for (i in 1:length(listfile)){
  res<-chischrom.stem_polyploidy(file=listfile[i])
  resultspoly.stem[i,1]<-res[[1]]
  resultspoly.stem[i,2]<-res[[3]]
}

resultspoly.stem<-as.data.frame(resultspoly.stem)
rownames(resultspoly.stem)<- listfile
colnames(resultspoly.stem)<-c("Chisq","p-value")
resultspoly.stem

#Loop para correr el test en todos los archivo y hacer una tabla DATOS con poliploidía CROWN NODE
resultspoly.crown<-matrix(nrow=length(listfile), ncol=2)
for (i in 1:length(listfile)){
  res<-chischrom.crown_polyploidy(file=listfile[i])
  resultspoly.crown[i,1]<-res[[1]]
  resultspoly.crown[i,2]<-res[[3]]
}

resultspoly.crown<-as.data.frame(resultspoly.crown)
rownames(resultspoly.crown)<- listfile
colnames(resultspoly.crown)<-c("Chisq","p-value")
resultspoly.crown

#Loop para correr el test en todos los archivo y hacer una tabla DATOS con poliploidía HALF WAY
resultspoly.halfway<-matrix(nrow=length(listfile), ncol=2)
for (i in 1:length(listfile)){
  res<-chischrom.halfway_polyploidy(file=listfile[i])
  resultspoly.halfway[i,1]<-res[[1]]
  resultspoly.halfway[i,2]<-res[[3]]
}

resultspoly.halway<-as.data.frame(resultspoly.halfway)
rownames(resultspoly.halfway)<- listfile
colnames(resultspoly.halfway)<-c("Chisq","p-value")
resultspoly.halfway


write.table(results.stem, file ="ChisqResults_stem_todo.txt", sep="\t", row.names=T)
write.table(resultspoly.stem, file ="ChisqResults_stem_todo_polyploidy.txt", sep="\t", row.names=T)
write.table(resultsdisplo.stem, file ="ChisqResults_stem_todo_disploidy.txt", sep="\t", row.names=T)

write.table(results.crown, file ="ChisqResults_crown_todo.txt", sep="\t", row.names=T)
write.table(resultspoly.crown, file ="ChisqResults_crown_todo_polyploidy.txt", sep="\t", row.names=T)
write.table(resultsdisplo.crown, file ="ChisqResults_crown_todo_disploidy.txt", sep="\t", row.names=T)

write.table(results.halfway, file ="ChisqResults_halfway_todo.txt", sep="\t", row.names=T)
write.table(resultspoly.halfway, file ="ChisqResults_halfway_todo_polyploidy.txt", sep="\t", row.names=T)
write.table(resultsdisplo.halfway, file ="ChisqResults_halfway_todo_disploidy.txt", sep="\t", row.names=T)


#Función para hacer los Chisq test en todos los nodos, solo tips y solo nodos internos STEM NODE
chischrom2.stem<-function(file){
  x<-read.table(file, header=T, strip.white = T)
  x$n<- substr(x$name, start=1, stop=2)
  nodes<-x[grep("N[1-9]", x$n),]
  tips<-x[-(grep("N[1-9]", x$n)),]
  preMed<-x[x$agestem>3.4,]
  Med<-x[x$agestem<3.4,]
  preMedt<-table(rowSums(preMed[,5:9], na.rm = T)==0)
  Medt<-table(rowSums(Med[,5:9], na.rm = T)==0)    
  chidata<-cbind(preMedt, Medt)
   restodo<-(chisq.test(chidata))
  
 #Chisqtest for tips STEM NODE
  preMedtips<-tips[tips$agestem>3.4,]
  Medtips<-tips[tips$agestem<3.4,]
  preMedtipst<-table(rowSums(preMedtips[,5:9], na.rm = T)==0)
  Medtipst<-table(rowSums(Medtips[,5:9], na.rm = T)==0)    
  chidatatips<-cbind(preMedtipst, Medtipst)
  restips<-(chisq.test(chidatatips))
  
  #Chisqtest for nodes STEM NODE 
  
  preMednodes<-nodes[nodes$agestem>3.4,]
  Mednodes<-nodes[nodes$agestem<3.4,]
  preMednodest<-table(rowSums(preMednodes[,5:9], na.rm = T)==0)
  Mednodest<-table(rowSums(Mednodes[,5:9], na.rm = T)==0)    
  chidatanodes<-cbind(preMednodest, Mednodest)
  resnodes<-(chisq.test(chidatanodes))
  
  chisq<-c(c(restodo[[1]],restodo[[3]]),c(restips[[1]],restips[[3]]), c(resnodes[[1]],resnodes[[3]]))
  chisq
    }

#Analisis para los stem con mutaciones 1/0 STEM NODE
setwd("/home/amesclir/Documents/Ubuntufiles/GITHUB/MedChromEvol/datasets/")
listfile<-dir()

resultschrom.stem<-matrix(nrow=length(listfile), ncol=6)
for (i in 1:length(listfile)){
    tryCatch(resultschrom.stem[i,]<-chischrom2.stem(file=listfile[i]), error=function(e) {
      print('Error')    })
 }

colnames(resultschrom.stem)<-c("Chisq_all","p-value_all","Chisq_tips","p-value_tips","Chisq_nodes","p-value_nodes")
row.names(resultschrom.stem)<-listfile
resultschrom.stem

setwd("/home/amesclir/Documents/Ubuntufiles/GITHUB/MedChromEvol/")
write.table(resultschrom.stem, file="/home/amesclir/Documents/Ubuntufiles/GITHUB/MedChromEvol/datasets/snalysis_stem_binarytransitions.txt", sep="\t")

#Función para hacer los Chisq test en todos los nodos, solo tips y solo nodos internos CROWN NODE
chischrom2.crown<-function(file){
  x<-read.table(file, header=T, strip.white = T)
  x$n<- substr(x$name, start=1, stop=2)
  nodes<-x[grep("N[1-9]", x$n),]
  tips<-x[-(grep("N[1-9]", x$n)),]
  preMed<-x[x$age>3.4,]
  Med<-x[x$ages<3.4,]
  preMedt<-table(rowSums(preMed[,5:9], na.rm = T)==0)
  Medt<-table(rowSums(Med[,5:9], na.rm = T)==0)    
  chidata<-cbind(preMedt, Medt)
  restodo<-(chisq.test(chidata))
  
  #Chisqtest for tips CROWN NODE
  preMedtips<-tips[tips$age>3.4,]
  Medtips<-tips[tips$age<3.4,]
  preMedtipst<-table(rowSums(preMedtips[,5:9], na.rm = T)==0)
  Medtipst<-table(rowSums(Medtips[,5:9], na.rm = T)==0)    
  chidatatips<-cbind(preMedtipst, Medtipst)
  restips<-(chisq.test(chidatatips))
  
  #Chisqtest for nodes CROWN NODE 
  
  preMednodes<-nodes[nodes$age>3.4,]
  Mednodes<-nodes[nodes$age<3.4,]
  preMednodest<-table(rowSums(preMednodes[,5:9], na.rm = T)==0)
  Mednodest<-table(rowSums(Mednodes[,5:9], na.rm = T)==0)    
  chidatanodes<-cbind(preMednodest, Mednodest)
  resnodes<-(chisq.test(chidatanodes))
  
  chisq<-c(c(restodo[[1]],restodo[[3]]),c(restips[[1]],restips[[3]]), c(resnodes[[1]],resnodes[[3]]))
  chisq
}

#Analisis para los stem con mutaciones 1/0 CROWN NODE
setwd("/home/amesclir/Documents/Ubuntufiles/GITHUB/MedChromEvol/datasets/")
listfile<-dir()

resultschrom.crown<-matrix(nrow=length(listfile), ncol=6)
for (i in 1:length(listfile)){
  tryCatch(resultschrom.crown[i,]<-chischrom2.crown(file=listfile[i]), error=function(e) {
    print('Error')    })
}

colnames(resultschrom.crown)<-c("Chisq_all","p-value_all","Chisq_tips","p-value_tips","Chisq_nodes","p-value_nodes")
row.names(resultschrom.crown)<-listfile
resultschrom.crown

setwd("/home/amesclir/Documents/Ubuntufiles/GITHUB/MedChromEvol/")
write.table(resultschrom.crown, file="/home/amesclir/Documents/Ubuntufiles/GITHUB/MedChromEvol/datasets/snalysis_crown_binarytransitions.txt", sep="\t")



#Función para hacer los Chisq test en todos los nodos, solo tips y solo nodos internos HALFWAY
chischrom2.halfway<-function(file){
  x<-read.table(file, header=T, strip.white = T)
  x$n<- substr(x$name, start=1, stop=2)
  nodes<-x[grep("N[1-9]", x$n),]
  tips<-x[-(grep("N[1-9]", x$n)),]
  preMed<-x[((x$agestem - x$age) + x$age) >3.4,]
  Med<-x[((x$agestem - x$age) + x$age)<3.4,]
  preMedt<-table(rowSums(preMed[,5:9], na.rm = T)==0)
  Medt<-table(rowSums(Med[,5:9], na.rm = T)==0)    
  chidata<-cbind(preMedt, Medt)
  restodo<-(chisq.test(chidata))
  
  #Chisqtest for tips HALFWAY
  preMedtips<-tips[((tips$agestem - tips$age) + tips$age) >3.4,]
  Medtips<-tips[((tips$agestem - tips$age) + tips$age) <3.4,]
  preMedtipst<-table(rowSums(preMedtips[,5:9], na.rm = T)==0)
  Medtipst<-table(rowSums(Medtips[,5:9], na.rm = T)==0)    
  chidatatips<-cbind(preMedtipst, Medtipst)
  restips<-(chisq.test(chidatatips))
  
  #Chisqtest for nodes HALFWAY NODE 
  
  preMednodes<-nodes[((nodes$agestem - nodes$age) + nodes$age) >3.4,]
  Mednodes<-nodes[((nodes$agestem - nodes$age) + nodes$age)<3.4,]
  preMednodest<-table(rowSums(preMednodes[,5:9], na.rm = T)==0)
  Mednodest<-table(rowSums(Mednodes[,5:9], na.rm = T)==0)    
  chidatanodes<-cbind(preMednodest, Mednodest)
  resnodes<-(chisq.test(chidatanodes))
  
  chisq<-c(c(restodo[[1]],restodo[[3]]),c(restips[[1]],restips[[3]]), c(resnodes[[1]],resnodes[[3]]))
  chisq
}

#Analisis para los stem con mutaciones 1/0 HALFWAY
setwd("/home/amesclir/Documents/Ubuntufiles/GITHUB/MedChromEvol/datasets/")
listfile<-dir()

resultschrom.halfway<-matrix(nrow=length(listfile), ncol=6)
for (i in 1:length(listfile)){
  tryCatch(resultschrom.halfway[i,]<-chischrom2.halfway(file=listfile[i]), error=function(e) {
    print('Error')    })
}

colnames(resultschrom.halfway)<-c("Chisq_all","p-value_all","Chisq_tips","p-value_tips","Chisq_nodes","p-value_nodes")
row.names(resultschrom.halfway)<-listfile
resultschrom.halfway

setwd("/home/amesclir/Documents/Ubuntufiles/GITHUB/MedChromEvol/")
write.table(resultschrom.halfway, file="/home/amesclir/Documents/Ubuntufiles/GITHUB/MedChromEvol/datasets/snalysis_halfway_binarytransitions.txt", sep="\t")
