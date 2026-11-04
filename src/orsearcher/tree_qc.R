args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=4L) stop('tree_qc.R tree.nwk id_map.tsv trimmed.faa output_dir')
suppressPackageStartupMessages(library(ape))
tree <- read.tree(args[1]); tree <- reorder(tree,'postorder')
ids <- read.delim(args[2],stringsAsFactors=FALSE,check.names=FALSE)
out <- args[4]; n <- length(tree$tip.label)
if(anyDuplicated(ids$tree_id) || anyDuplicated(tree$tip.label)) stop('Duplicate tree labels')
ids <- ids[match(tree$tip.label,ids$tree_id),]
if(anyNA(ids$tree_id)) stop('Unmapped tips')
refOR <- ids$state=='reference' & ids$family=='OR'
refGR <- ids$state=='reference' & ids$family=='GR'
if(!any(refOR) || !any(refGR)) stop('Both OR and GR references are required')
# Postorder accumulation is iterative, O(nodes), and independent of tree depth.
oc <- gc <- numeric(n+tree$Nnode)
oc[seq_len(n)] <- refOR; gc[seq_len(n)] <- refGR
for(i in seq_len(nrow(tree$edge))) {
    parent <- tree$edge[i,1]; child <- tree$edge[i,2]
    oc[parent] <- oc[parent]+oc[child]; gc[parent] <- gc[parent]+gc[child]
}
child <- tree$edge[,2]
errors <- cbind(gc[child]+sum(refOR)-oc[child],oc[child]+sum(refGR)-gc[child])
best <- which(errors==min(errors),arr.ind=TRUE)[1,]
children <- split(tree$edge[,2],tree$edge[,1])
stack <- tree$edge[best[1],2]; inside <- integer(0)
while(length(stack)) {
    node <- stack[length(stack)]; stack <- stack[-length(stack)]
    if(node<=n) inside <- c(inside,node)
    else stack <- c(stack,children[[as.character(node)]])
}
orSide <- rep(best[2]==2,n); orSide[inside] <- best[2]==1
fa <- readLines(args[3]); headers <- which(startsWith(fa,'>')); seqs <- list()
for(i in seq_along(headers)) {
    last <- if(i<length(headers)) headers[i+1]-1L else length(fa)
    label <- strsplit(sub('^>','',fa[headers[i]]),'[[:space:]]')[[1]][1]
    seqs[[label]] <- if(last>headers[i]) paste(fa[seq.int(headers[i]+1L,last)],collapse='') else ''
}
gap <- vapply(tree$tip.label,function(id) {
    s <- strsplit(seqs[[id]],'')[[1]]
    if(!length(s)) return(1)
    mean(s %in% c('-','X','?'))
},numeric(1))
qc <- data.frame(candidate_id=ids$original_id,tree_side=ifelse(orSide,'OR_side','GR_side'),gap_fraction=gap,
    reference_split_error_fraction=min(errors)/(sum(refOR)+sum(refGR)),stringsAsFactors=FALSE)
# A global OR/GR cut alone is not reliable identity evidence. Require both
# reference separation and three concordant nearest distinct reference tips,
# adequate aligned residues, and local/global agreement. Missing QC is unknown.
distance <- cophenetic.phylo(tree)
referenceTips <- which(ids$state=='reference')
qc$nearest_reference_family <- 'mixed'; qc$nearest_reference_ids <- ''
qc$identity_tree_vote <- 'unresolved'
for(i in which(ids$state=='candidate')) {
    near <- referenceTips[order(distance[i,referenceTips])][seq_len(min(3L,length(referenceTips)))]
    families <- ids$family[near]
    consensus <- if(length(near)==3L && length(unique(families))==1L) families[1] else 'mixed'
    qc$nearest_reference_family[i] <- consensus
    qc$nearest_reference_ids[i] <- paste(ids$original_id[near],collapse=',')
    aligned <- nchar(gsub('[-X?]','',seqs[[tree$tip.label[i]]]))
    if(consensus %in% c('OR','GR') && gap[i]<=.6 && aligned>=100L &&
       qc$reference_split_error_fraction[i]<=.05 &&
       qc$tree_side[i]==paste0(consensus,'_side')) qc$identity_tree_vote[i] <- consensus
}
write.table(qc[ids$state=='candidate',],file.path(out,'candidate_qc.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
reference <- data.frame(reference_id=ids$original_id,family=ids$family,tree_side=qc$tree_side,stringsAsFactors=FALSE)
write.table(reference[ids$state=='reference',],file.path(out,'reference_qc.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(c(paste('tips',n),paste('OR_references',sum(refOR)),paste('GR_references',sum(refGR)),
    paste('best_split_reference_errors',min(errors)),paste('candidates_on_GR_side',sum(!orSide & ids$state=='candidate'))),file.path(out,'summary.txt'))
colors <- ifelse(ids$state=='candidate',ifelse(orSide,'#117733','#CC3311'),ifelse(ids$family=='OR','#4477AA','#EE7733'))
# Graphics failure is a warning, not failure of biological-result export.
tryCatch({
    pdf(file.path(out,'OR_GR_tree.pdf'),width=12,height=12)
    plot(tree,type='unrooted',show.tip.label=FALSE,edge.width=.35,no.margin=TRUE)
    tiplabels(pch=19,col=colors,cex=.3)
    legend('topleft',c('reference OR','reference GR','candidate: OR side','candidate: GR side'),
        col=c('#4477AA','#EE7733','#117733','#CC3311'),pch=19,bty='n',cex=.8)
    title(main='OR/GR reference-guided unrooted tree QC',sub=paste('Reference split errors:',min(errors),'/',sum(refOR)+sum(refGR)))
    dev.off()
},error=function(e) {
    if(dev.cur()>1) dev.off()
    writeLines(conditionMessage(e),file.path(out,'plot_warning.txt'))
})
