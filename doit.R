library(data.table)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(ggplot2)
library(gtable)


load('todo.RData')
todo[,.N,sample][order(-N)][,sample]->sams
todo[,.N,gene][order(-N)][,gene]->genes

melt(dcast(todo,gene~sample,value.var='mutation'),id.var='gene')->todo
names(todo)<-c('gene','sample','mutation')
todo[is.na(todo)]<-'None'

todo[,mut:=mutation]
todo[grep("///",mut),mut:='multipe hit']


mut.Cate<-c("Substitution - Missense","Substitution - Nonsense",
"Insertion - In frame","Insertion - Frameshift",
"Deletion - In frame","Deletion - Frameshift",
"Nonstop extension","Whole gene deletion")

setNames(brewer.pal(n=10,'Paired')[0-5:6],mut.Cate)->mut.Cate
c(mut.Cate,'None'='gray','multipe hit'='red')->mut.Cate

todo[,mut:=ordered(mut,levels=names(mut.Cate))]


inter.x<-setNames(ppoints(length(sams),1/2),sams)
inter.y<-setNames(rev(ppoints(length(genes),1/2)),genes)

todo[,x:=inter.x[sample]]
todo[,y:=inter.y[gene]]
todo[,w:=1/length(sams)]
todo[,h:=1/length(genes)]


################################the ggplot graph
tu<-ggplot(todo,aes(x=x,y=y,height=h,width=w,fill=mut)) + 
    geom_tile(colour='white',show.legend=F)+
    scale_x_continuous(expand=expansion(),breaks=inter.x,labels=names(inter.x))+
    scale_y_continuous(expand=expansion(),breaks=inter.y,labels=names(inter.y),position = "right")+
    scale_fill_manual(values=mut.Cate) +
    theme_minimal() %+replace% theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position='top') +
    labs(x='',y='')


gai<-ggplotGrob(tu)

##########################the legends
leg<-ggplot(todo[mut!='None'],aes(x=x,y=y,height=h,width=w,fill=mut)) + 
     geom_tile()+
     scale_fill_manual(values=mut.Cate) +
     guides(fill= guide_legend(title='',nrow=1)) + theme(legend.position='top',legend.direction='horizontal')

leg<-ggplotGrob(leg)
leg<-leg$grobs[[grep('guide',leg$layout$name)]]


#######################################ok1
zz<-todo[mut=='multipe hit']

apply(zz[,.(x,y,w,h)],1,function(xxx)
{
segmentsGrob(x0=c(0,0),x1=c(1,1),y0=c(0,1),y1=c(1,0),gp=gpar(col='#2B1B17'),vp=viewport(xxx[1],xxx[2],xxx[3],xxx[4]))
})->rs

append(rs,list(rectGrob(zz$x,zz$y,zz$w,zz$h,gp=gpar(col='#2B1B17'))))->rs
addOn<-gTree(children=Reduce(gList,rs))
gtable_add_grob(gai,addOn,t=7,l=5)->ok1


##################################ok2
apply(zz[,.(x,y,w,h,mutation)],1,function(xxx)
{
strsplit(xxx[5]," /// ")[[1]]->jus
grid.rect(x=ppoints(length(jus),0.5),width=1/length(jus),vp=viewport(xxx[1],xxx[2],xxx[3],xxx[4]),gp=gpar(fill=mut.Cate[jus],col='#2B1B17'))
})->rs

addOn<-gTree(children=Reduce(gList,rs))
gtable_add_grob(gai,addOn,t=7,l=5)->ok2


###################################combination
grid.arrange(grobs=list(leg,gai,ok1,ok2),layout_matrix=rbind(c(NA,1,NA),c(2,3,4)),heights=c(1,8))
