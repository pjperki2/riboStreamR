shinyServer(function(input,output,session){
  options(shiny.maxRequestSize=1000*1024^2)
  obs_p1_add<-observe({
    print (input$go)
  })
  
  getGenes = reactive({
    withProgress(message = 'Updating Inputs', value = 0, {
      anno = input$ANNOfile
      txdb <- makeTxDbFromGFF(anno$datapath, format="gff3", organism="Arabidopsis")
      eByg <- exonsBy(txdb, by=c("gene"))
      names(eByg)
    })
  })
  getSamps = reactive({
    withProgress(message = 'Updating Inputs', value = 0, {
      rnafiles = input$RNAfile
      rpffiles = input$RPFfile
      c(rnafiles$name, rpffiles$name)
    })
  })
  
  output$inputvals1 <- renderUI({
    nums1 <- as.numeric(input$expers)
    nums2 <- as.numeric(input$treats)
    if (nums1 == 1){
      lapply(1:nums2, function(i) {
        selectInput(paste("Treat_", i, sep = ""),paste("Treatment ", i, sep = "" ), choices = getSamps(), multiple = T)
      })
    }
  })
  output$inputvals2 <- renderUI({
    nums1 <- as.numeric(input$expers)
    nums2 <- as.numeric(input$treats)
    n = nums1*nums2
    n2 = rep(c(1:nums2), times = 2)
    nams = rep(c("RNA ", "Ribo "), each = nums2)
    if(nums1 > 1){
      lapply(1:n, function(i) {
        selectInput(paste(nams[i],n2[i],sep = ""),paste(nams[i], " Treat ",n2[i] , sep = "" ), choices = getSamps(), multiple = T)
      })
    }
  })
  
  
  ##################################################
  #Set data set input choices
  ##################################################
  output$choose_dataset1 <- renderUI({
    selectInput("samples", "Data sets", getSamps(),getSamps(),TRUE)
  })
  output$choose_dataset2 <- renderUI({
    selectInput("samples2", "Data sets", getSamps(),getSamps(),TRUE)
  })
  output$choose_dataset3 <- renderUI({
    selectInput("samples3", "Data sets", getSamps(),getSamps(),TRUE)
  })
  output$choose_dataset4 <- renderUI({
    selectInput("samples4", "Data sets", getSamps(),getSamps(),TRUE)
  })
  output$choose_dataset5 <- renderUI({
    selectInput("samples5", "Data sets", getSamps(),getSamps(),TRUE)
  })
  output$choose_dataset6 <- renderUI({
    selectInput("samples6", "Data sets", getSamps(),getSamps(),TRUE)
  })
  output$choose_dataset7 <- renderUI({
    selectInput("samples7", "Data sets", getSamps(),getSamps(),TRUE)
  })
  output$choose_dataset8 <- renderUI({
    selectInput("samples8", "Data sets", getSamps(),getSamps(),TRUE)
  })
  
  ##################################################
  #
  ##################################################
  output$choose_dataset15 <- renderUI({
    selectInput("samples15", "Choose Control", c(getSamps(), "none"),"none",FALSE)
  })
  output$choose_dataset25 <- renderUI({
    selectInput("samples25", "Choose Control", c(getSamps(), "none"),"none",FALSE)
  })
  output$choose_dataset35 <- renderUI({
    selectInput("samples35", "Choose Control", c(getSamps(), "none"),"none",FALSE)
  })
  output$choose_dataset45 <- renderUI({
    selectInput("samples45", "Choose Control", c(getSamps(), "none"),"none",FALSE)
  })
  output$choose_dataset55 <- renderUI({
    selectInput("samples55", "Choose Control", c(getSamps(), "none"),"none",FALSE)
  })
  output$choose_dataset65 <- renderUI({
    selectInput("samples65", "Choose Control", c(getSamps(), "none"),"none",FALSE)
  })
  output$choose_dataset75 <- renderUI({
    selectInput("samples75", "Choose Control", c(getSamps(), "none"),"none",FALSE)
  })
  output$choose_dataset85 <- renderUI({
    selectInput("samples85", "Choose Control", c(getSamps(), "none"),"none",FALSE)
  })
  ##################################################
  #
  ##################################################
  ##################################################
  #Set gene input choices
  ##################################################
  output$choose_gene <- renderUI({
    selectInput("genes", "Genes", c("all",getGenes()), "all", T)
  })
  output$choose_gene2 <- renderUI({
    selectInput("genes2", "Genes", c("all",getGenes()), "all", T)
  })
  output$choose_gene3 <- renderUI({
    selectInput("genes3", "Genes", c("all",getGenes()), "all", T)
  })
  output$choose_gene4 <- renderUI({
    selectInput("genes4", "Genes", c("all",getGenes()), "all", T)
  })
  output$choose_gene5 <- renderUI({
    selectInput("genes5", "Genes", c("all",getGenes()), "all", T)
  })
  output$choose_gene6 <- renderUI({
    selectInput("genes6", "Genes", c("all",getGenes()), "AT3G09260", T)
  })
  output$choose_gene7 <- renderUI({
    selectInput("genes7", "Genes", c("all",getGenes()), "all", T)
  })
  output$choose_gene8 <- renderUI({
    selectInput("genes8", "Genes", c("all",getGenes()), "all", T)
  })
  
  ##################################################
  #Preprocess reactive
  ##################################################
  getTxdb = reactive({
    anno = input$ANNOfile
    txdb <- makeTxDbFromGFF(anno$datapath, format="gff3", organism="Arabidopsis")
    txdb
  })
  getFeat = reactive({
    txdb = getTxdb()
    feat <- genFeatures(txdb, featuretype="all", reduce_ranges=TRUE, upstream=1000, downstream=0,verbose=TRUE)
    feat
  })
  
  grang = reactive({
    rnafiles = input$RNAfile
    rpffiles = input$RPFfile
    ribofiles = c(rnafiles$datapath,rpffiles$datapath)
    ribonames = c(rnafiles$name, rpffiles$name)
    grl = list()
    withProgress(message = 'Preprocessing', value = 0, {
      txdb = getTxdb()
      feat = getFeat()
      for (i in 1:length(ribofiles)){
        ga = readGappedReads(ribofiles[i])
        gr = GRanges(ga)
        #m = sample(1:length(ga), 1000000)
        #ga = ga[m]
        gr = GRanges(ga)
        gr$sample = ribonames[i]
        gr$strands = as.character(gr@strand)
        gr$widths = width(gr)
        sbp = ScanBamParam(what = "mapq")
        mapq = unlist(scanBam(ribofiles[i], param = sbp))
        #mapq = mapq[m]
        gr$mapq = mapq
        gr$mapping = "Multi"
        gr$mapping[which(gr$mapq == 50)] = "Unique"
        reads <- qseq(ga)  
        gr$gc <- rowSums( alphabetFrequency( reads, as.prob = TRUE )[ , c( "C", "G" ) ] )
        gr$feature = "other"
        featNames = c("pseudogenic_transcript", "ncRNA","intergenic", "intron", 
                      "miRNA", "snRNA", "snoRNA","promoter", "mRNA", "exon", 
                      "fiveUTR", "threeUTR","cds", "tRNA", "rRNA")
        features = list(feat$pseudogenic_transcript, feat$ncRNA,feat$intergenic, feat$intron, 
                        feat$miRNA, feat$snRNA, feat$snoRNA,feat$promoter, feat$mRNA, feat$exon, 
                        feat$fiveUTR, feat$threeUTR,feat$cds, feat$tRNA, feat$rRNA)
        
        for (k in 1:length(features)){
          o = findOverlaps(gr, features[[k]], ignore.strand = T)
          p = queryHits(o)
          gr$feature[p] = featNames[k]
        }
        gr$gene = "other"
        eByg <- exonsBy(txdb, by=c("gene"))
        cByg <- cdsBy(txdb, by=c("gene"))#check to make sure this is right
        cb = cdsBy(txdb,"tx", use.names =T)
        tb = transcriptsBy(txdb,"gene")
        tbo = findOverlaps(tb,gr)
        cbo = findOverlaps(cb,gr)
        gr$gene[subjectHits(tbo)]  = names(tb[queryHits(tbo)])
        gr$transcript = "other"
        gr$transcript[subjectHits(cbo)]  = names(cb[queryHits(cbo)])
        tstart = as.vector(start(ranges(cb)))
        tstart2 = match(gr$transcript, names(tstart))
        tind = which(!is.na(tstart2))
        gr$frame = "none"
        gr$distA = "none"
        gr$distB = "none"
        gr$posA = NA
        gr$posB = NA
        gr$loc = 0
        gr$tag = "none"
        neg = which(gr[tind]$strands == "-")
        plus = which(gr[tind]$strands == "+")
        gr$distA[tind][plus] = (as.vector(start(ranges(gr[tind][plus]))) - min(start(ranges(cb[names(tstart[tstart2[tind][plus]])]))))
        gr$distA[tind][neg] = (max(end(ranges(cb[names(tstart[tstart2[tind][neg]])]))) - as.vector(end(ranges(gr[tind][neg]))))
        gr$distA = as.numeric(gr$distA) + 1
        gr$distB[tind][plus] = (max(end(ranges(cb[names(tstart[tstart2[tind][plus]])]))) - as.vector(start(ranges(gr[tind][plus]))))
        gr$distB[tind][neg] = (as.vector(end(ranges(gr[tind][neg]))) - min(start(ranges(cb[names(tstart[tstart2[tind][neg]])]))))
        gr$distB = as.numeric(gr$distB)
        mtt = mapToTranscripts(gr, cByg)
        gr$posA[mtt$xHits] = start(ranges(mtt))
        gr$posB[mtt$xHits] = sum(width(cByg[mtt$transcriptsHits])) - as.numeric(gr$posA[mtt$xHits])
        neg = which(gr$strands == "-")
        plus = which(gr$strands == "+")
        gr[plus]$loc = start(ranges(gr[plus]))
        gr[neg]$loc = end(ranges(gr[neg]))
        incProgress(1/length(ribofiles), detail = paste("File ", i, sep = ""))
        grl[i] = gr
      }
    })
    GRangesList(grl)
  })
  
  #toListen <- reactive({
  #  list(input$ANNOfile,input$RNAfile, input$RPFfile)
  #})
  gettab <- eventReactive(input$prepro,{
    #if(is.null(input$ANNOfile) || is.null(input$RNAfile) || is.null(input$RPFfile)){
    #  return()
    #}
    print("yes")
    uList = adjust()
    rnafiles = input$RNAfile
    rpffiles = input$RPFfile
    ribonames = c(rnafiles$name, rpffiles$name)
    rnaind = which(ribonames %in% rnafiles$name)
    rpfind = which(ribonames %in% rpffiles$name)
    types = c(rep("RNA-seq", length(rnaind)), rep("Ribo-seq", length(rpfind)))
    tabl = data.frame()
    for (i in 1:length(uList)){
      tabl[i,1] = ribonames[i]
      tabl[i,2] = types[i]
      tabl[i,3] = length(uList[[i]])
    }
    colnames(tabl) = c("Sample Names", "Type", "Number of alignments")
    tabl
  })
  
  output$starter = renderDataTable({
    print(gettab())
  })
  
  
  offs = reactive({
    rnafiles = input$RNAfile
    rpffiles = input$RPFfile
    ribonames = c(rnafiles$name, rpffiles$name)
    rpfind = which(ribonames %in% rpffiles$name)
    gr = grang()
    offset1 = 1:7
    offsets = list()
    withProgress(message = 'Getting Offsets', value = 0, {
      m = 1
      for(k in rpfind){
        exp1 = gr[[k]]
        for (i in 26:32){
          clen = exp1[which(exp1$widths == i)]
          clen = clen[which(clen$distA < 40)]
          clen = clen[which(clen$distA > 0)]
          clen$frame = clen$distA %% 3
          ftab = table(clen$frame)
          frame = as.numeric(names(which.max(ftab)))
          
          clen = exp1[which(exp1$widths == i)]
          clen = clen[which(clen$distA < 0)]
          ta = table(clen$distA)
          ta = ta[-1]
          off1 = names(which.max(ta))
          offframe = as.numeric(off1)%%3
          print(offframe)
          print(frame)
          if (offframe != frame){
            t2 = which(names(ta) == off1)
            print(as.numeric(names(ta[t2-1]))%%3)
            print(ta)
            print(t2)
            if (t2 != length(ta) & t2 != 1){
              if (as.numeric(names(ta[t2+1]))%%3 == frame){
                off1 = names(ta[t2+1])
              }
              if (as.numeric(names(ta[t2-1]))%%3 == frame){
                off1 = names(ta[t2-1])
              }
            }
            if (t2 == length(ta)){
              if (frame == 0){
                off1 = 0
              }
              if (as.numeric(names(ta[t2-1]))%%3 == frame){
                off1 = names(ta[t2-1])
              }
            }
            if (t2 == 1){
              if (as.numeric(names(ta[t2+1]))%%3 == frame){
                off1 = names(ta[t2+1])
              }
              if (as.numeric(names(ta[t2+1]))%%3 != frame){
                off1 = names(ta[t2])
              }
            }
          }
          offset1[i-25] = off1
        }
        offsets[[m]] = c(as.numeric(offset1))
        incProgress(1/m, detail = paste("File ", m, sep = ""))
        m = m+1
      }
      print(offsets)
    })
    print(offsets)
  })
  
  adjust = reactive({#input$pre,{
    rnafiles = input$RNAfile
    rpffiles = input$RPFfile
    ribonames = c(rnafiles$name, rpffiles$name)
    rpfind = which(ribonames %in% rpffiles$name)
    gr = grang()
    print(grang())
    offsets = offs()
    withProgress(message = 'Adjusting Positions', value = 0, {
      m = 1
      for (k in rpfind){
        gr[[k]]$distB = as.numeric(gr[[k]]$distB)
        os = offsets[[m]]
        ls = c(26:32)
        exp = gr[[k]]
        for ( i in 1:7){
          sub = which(exp$widths == ls[i])
          exp1 = exp[sub]
          gr[[k]][sub]$distA = exp1$distA - os[i]
          gr[[k]][sub]$posA = exp1$posA - os[i]
          gr[[k]][sub]$distB = as.numeric(exp1$distB) + os[i]
          gr[[k]][sub]$posB = exp1$posB + os[i]
        }
        m = m+1
      }
      for (j in 1:length(ribonames)){
        gr[[j]][which(gr[[j]]$distA < 25)]$posA = gr[[j]][which(gr[[j]]$distA < 25)]$distA
        endind = which(gr[[j]]$distB < 40)
        gr[[j]][endind]$frame = (((-1 * gr[[j]][endind]$distB)-1) %% 3)
        distind = which(!is.na(gr[[j]]$posA))
        gr[[j]][distind]$frame = gr[[j]][distind]$posA %% 3
      }
      incProgress(1/k, detail = paste("File ", k, sep = ""))
      gr
    }) 
  })
  
  
  
  getrld = eventReactive(input$submit,{
    withProgress(message = 'Generating Plot', value = 0, {
      uList = unlist(adjust())
      print(uList)
      sampInd = which(uList$sample %in% input$samples)
      strandInd = which(uList$strands %in% input$strand)
      mapInd = which(uList$mapping %in% input$mapping)
      featInd = which(uList$feature %in% input$features)
      print(input$genes)
      if(! "all" %in% input$genes){
        geneInd = which(uList$gene %in% input$genes)
      }
      if("all" %in% input$genes){
        geneInd = 1:length(uList)
      }
      frameInd = which(uList$frame %in% input$frame)
      print(length(sampInd))
      print(length(strandInd))
      print(length(frameInd))
      print(length(featInd))
      print(length(geneInd))
      print(length(mapInd))
      
      common = Reduce(intersect,list(sampInd,strandInd,mapInd,featInd, geneInd, frameInd))
      uList = uList[common]
      print(uList)
      if (length(input$sepFeats) >0){
        sepInd = which(uList$feature %in% input$sepFeats)
        v = 1:length(uList)
        v = v[-sepInd]
        if (length(v)>0){
          uList[v]$feature = "other"
        }
      }
      if(input$splitPlot != "none"){
        if(input$splitLine != "none"){
          uList = split(uList, f = mcols(uList)[which(names(mcols(uList)) == input$splitPlot)][,1])
          plotList = list()
          avec = vector()
          
          for (i in 1:length(uList)){
            uList2 = split(uList[[i]], f = mcols(uList[[i]])[which(names(mcols(uList[[i]])) == input$splitLine)][,1])
            w = vector()
            n = vector()
            ttl = mcols(uList[[i]])[which(names(mcols(uList[[i]])) == input$splitPlot)][1,1]
            for (k in 1:length(uList2)){
              w = c(w,uList2[[k]]$widths)
              n = c(n,mcols(uList2[[k]])[which(names(mcols(uList2[[k]])) == input$splitLine)][,1])
            }
            df = data.frame(w,n)
            gg = ggplot(df, aes(x=w)) + geom_density(aes(group=n,y = ..scaled.., color = n),bw=input$bw) +xlim(input$xmin,input$xmax) + ylim(0,1.1) + labs(title = ttl, x = "Length") +
              scale_color_brewer(palette=input$cols1)
            plotList[[i]] = gg
          }
          print(plot_grid(plotlist = plotList,ncol=floor(i/2)))
          
        }
        if(input$splitLine =="none"){
          uList = split(uList, f = mcols(uList)[which(names(mcols(uList)) == input$splitPlot)][,1])
          plotList = list()
          for (i in 1:length(uList)){
            ttl = mcols(uList[[i]])[which(names(mcols(uList[[i]])) == input$splitPlot)][1,1]
            w = uList[[i]]$widths
            #n = mcols(uList[[i]])[which(names(mcols(uList[[i]])) == input$splitPlot)][,1]
            df = data.frame(w)
            gg = ggplot(df, aes(x=w)) + geom_density(aes(y = ..scaled..),bw=input$bw) +xlim(input$xmin,input$xmax) + ylim(0,1.1) + labs(title = ttl, x = "Length") + scale_color_brewer(palette=input$cols1)
            plotList[[i]] = gg
          }
          print(plot_grid(plotlist = plotList,ncol=floor(i/2)))
        }
      }
      if (input$splitPlot == "none"){
        if(input$splitLine != "none"){
          uList2 = split(uList, f = mcols(uList)[which(names(mcols(uList)) == input$splitLine)][,1])
          w = vector()
          n = vector()
          for (k in 1:length(uList2)){
            w = c(w,uList2[[k]]$widths)
            n = c(n,mcols(uList2[[k]])[which(names(mcols(uList2[[k]])) == input$splitLine)][,1])
          }
          df = data.frame(w,n)
          print(ggplot(df, aes(x=w)) + geom_density(aes(group=n,y = ..scaled.., color = n),bw=input$bw) + xlim(input$xmin,input$xmax) + ylim(0,1.1) + labs(title = "Read Lengths", x = "Length") + scale_color_brewer(palette=input$cols1))
        }
        
        if(input$splitLine =="none"){
          print(length(uList))
          w = uList$widths
          df = data.frame(w)
          print(ggplot(df, aes(x=w)) + geom_density(aes(y = ..scaled..),bw=input$bw) +xlim(input$xmin,input$xmax) + ylim(0,1.1) + labs(title = "Read Lengths", x = "Length") + scale_color_brewer(palette=input$cols1))
        }
      }
    })
  })
  
  output$dens = renderPlot({
    getrld()
  })
  
  plotComp = eventReactive(input$submit2,{
    withProgress(message = 'Generating Plot', value = 0, {
      uList = unlist(adjust())
      sampInd = which(uList$sample %in% input$samples2)
      strandInd = which(uList$strands %in% input$strand2)
      mapInd = which(uList$mapping %in% input$mapping2)
      featInd = which(uList$feature %in% input$features2)
      if(input$genes2 != "all"){geneInd = which(uList$gene %in% input$genes2)}
      if(input$genes2 == "all"){geneInd = 1:length(uList)}
      frameInd = which(uList$frame %in% input$frame2)
      lengthInd = which(uList$widths %in% seq(input$length[1],input$length[2]))
      common = Reduce(intersect,list(sampInd,strandInd,mapInd,featInd, geneInd, frameInd, lengthInd))
      uList = uList[common]
      if (length(input$sepFeats2) >0){
        sepInd = which(uList$feature %in% input$sepFeats2)
        v = 1:length(uList)
        v = v[-sepInd]
        if (length(v)>0){
          uList[v]$feature = "other"
        }
      }
      if (input$splitGroup != "none"){
        nuList = split(uList, f = mcols(uList)[which(names(mcols(uList)) == input$splitBar)][,1])
        listy = list()
        m = 1
        for (i in 1:length(nuList)){
          uList2 = split(nuList[[i]], f = mcols(nuList[[i]])[which(names(mcols(nuList[[i]])) == input$splitGroup)][,1])
          for (k in 1:length(uList2)){
            f = mcols(uList2[[k]][,which(names(mcols(uList2[[k]])) == input$splitCol)])[,1]
            x = table(f)/length(f)
            listy[[m]] = x
            m = m+1
          }
        }
        m2 = do.call(rbind,listy)
        m = do.call(rbind, lapply(lapply(listy, unlist), "[",unique(unlist(c(sapply(listy,names))))))
        colnames(m) <- unique(unlist(c(sapply(listy,names))))
        m[is.na(m)] <- 0
        onenames = unique(mcols(uList[,which(names(mcols(uList)) == input$splitBar)])[,1])
        twonames = colnames(m)
        threenames = sort(unique(mcols(uList[,which(names(mcols(uList)) == input$splitGroup)])[,1]))
        dat = expand.grid(onenames, twonames, threenames)
        names(dat) = c(input$splitBar, input$splitCol, input$splitGroup)
        d1 = length(threenames)
        d2 = length(onenames) * length(twonames)
        dim(m) = c(d1,d2)
        m = t(m)
        dat$percentage = as.vector(m)
        print(ggplot() + geom_bar(data = dat, aes(y = percentage, x = dat[,which(names(dat) == input$splitBar)], fill = dat[,which(names(dat) == input$splitCol)]),stat = "identity",position = "stack") + theme_bw() + facet_grid(paste("~", input$splitGroup)) + scale_fill_discrete(name = input$splitCol) + xlab(input$splitBar) + scale_color_brewer(palette=input$cols2)) 
        
      }
      if(input$splitGroup == "none" ){
        if (input$splitCol != "none"){
          nuList = split(uList, f = mcols(uList)[which(names(mcols(uList)) == input$splitBar)][,1])
          listy = list()
          for (i in 1:length(nuList)){
            f = mcols(nuList[[i]][,which(names(mcols(nuList[[i]])) == input$splitCol)])[,1]
            x = table(f)/length(f)
            listy[[i]] = x
          }
          m = do.call(rbind, lapply(lapply(listy, unlist), "[",
                                    unique(unlist(c(sapply(listy,names))))))
          colnames(m) <- unique(unlist(c(sapply(listy,names))))
          m[is.na(m)] <- 0
          rownames(m) = sort(unique(mcols(uList[,which(names(mcols(uList)) == input$splitBar)])[,1]))
          dat = expand.grid(rownames(m), colnames(m))
          names(dat) = c(input$splitBar, input$splitCol)
          dat$percentage = as.vector(m)
          print(ggplot() + geom_bar(data = dat, aes(y = percentage, x = dat[,which(names(dat) == input$splitBar)], fill = dat[,which(names(dat) == input$splitCol)]), stat = "identity", position = "dodge") + scale_fill_discrete(name = input$splitCol) + xlab(input$splitBar) + scale_color_brewer(palette=input$cols2))
          
        }
      }
    })   
  })
  
  output$comp = renderPlot({
    plotComp()
  })
  
  plotGc = eventReactive(input$submit3,{
    withProgress(message = 'Generating Plot', value = 0, {
      uList = unlist(adjust())
      sampInd = which(uList$sample %in% input$samples3)
      strandInd = which(uList$strands %in% input$strand3)
      mapInd = which(uList$mapping %in% input$mapping3)
      if(input$genes3 != "all"){geneInd = which(uList$gene %in% input$genes3)}
      if(input$genes3 == "all"){geneInd = 1:length(uList)}
      frameInd = which(uList$frame %in% input$frame3)
      featInd = which(uList$feature %in% input$features3)
      lengthInd = which(uList$widths %in% seq(input$length3[1],input$length3[2]))
      common = Reduce(intersect,list(sampInd,strandInd,mapInd,featInd, geneInd,frameInd,lengthInd))
      uList = uList[common]
      if (length(input$sepFeats3) >0){
        sepInd = which(uList$feature %in% input$sepFeats3)
        v = 1:length(uList)
        v = v[-sepInd]
        if (length(v)>0){
          uList[v]$feature = "other"
        }
      }
      if(input$splitPlot3 != "none"){
        if(input$splitLine3 != "none"){
          uList = split(uList, f = mcols(uList)[which(names(mcols(uList)) == input$splitPlot3)][,1])
          plotList = list()
          for (i in 1:length(uList)){
            uList2 = split(uList[[i]], f = mcols(uList[[i]])[which(names(mcols(uList[[i]])) == input$splitLine3)][,1])
            w = vector()
            n = vector()
            ttl = mcols(uList[[i]])[which(names(mcols(uList[[i]])) == input$splitPlot3)][1,1]
            for (k in 1:length(uList2)){
              w = c(w,uList2[[k]]$gc)
              n = c(n,mcols(uList2[[k]])[which(names(mcols(uList2[[k]])) == input$splitLine3)][,1])
            }
            df = data.frame(w,n)
            gg = ggplot(df, aes(x=w)) + geom_density(aes(group=n,y = ..scaled.., color = n),bw=.01) +xlim(0,1) + ylim(0,1) + labs(title = ttl, x = "GC%") + scale_color_brewer(palette=input$cols3)
            plotList[[i]] = gg
          }
          print(plot_grid(plotlist = plotList,ncol=floor(i/2)))
        }
        if(input$splitLine3 =="none"){
          uList = split(uList, f = mcols(uList)[which(names(mcols(uList)) == input$splitPlot3)][,1])
          plotList = list()
          for (i in 1:length(uList)){
            ttl = mcols(uList[[i]])[which(names(mcols(uList[[i]])) == input$splitPlot3)][1,1]
            w = uList[[i]]$gc
            #n = mcols(uList[[i]])[which(names(mcols(uList[[i]])) == input$splitPlot)][,1]
            df = data.frame(w)
            gg = ggplot(df, aes(x=w)) + geom_density(aes(y = ..scaled..),bw=.01) +xlim(0,1) + ylim(0,1) + labs(title = ttl, x = "GC%") + scale_color_brewer(palette=input$cols3)
            plotList[[i]] = gg
          }
          print(plot_grid(plotlist = plotList,ncol=floor(i/2)))
        }
      }
      if (input$splitPlot3 == "none"){
        if(input$splitLine3 != "none"){
          uList2 = split(uList, f = mcols(uList)[which(names(mcols(uList)) == input$splitLine3)][,1])
          w = vector()
          n = vector()
          for (k in 1:length(uList2)){
            w = c(w,uList2[[k]]$gc)
            n = c(n,mcols(uList2[[k]])[which(names(mcols(uList2[[k]])) == input$splitLine3)][,1])
          }
          df = data.frame(w,n)
          print(ggplot(df, aes(x=w)) + geom_density(aes(group=n,y = ..scaled.., color = n),bw=.01) + xlim(0,1) + ylim(0,1) + labs(title = "GC", x = "GC%") + scale_color_brewer(palette=input$cols3))
        }
        
        if(input$splitLine3 =="none"){
          print(length(uList))
          w = uList$gc
          df = data.frame(w)
          print(ggplot(df, aes(x=w)) + geom_density(aes(y = ..scaled..),bw=.01) + xlim(0,1) + ylim(0,1) + labs(title = "GC", x = "GC%") + scale_color_brewer(palette=input$cols3))
        }
      }
    })
  })
  
  output$gc = renderPlot({
    plotGc()
  })
  
  plotTab = eventReactive(input$submit4,{
    withProgress(message = 'Generating Table', value = 0, {
      grr = adjust()
      uList = unlist(grr)
      sampInd = which(uList$sample %in% input$samples4)
      strandInd = which(uList$strands %in% input$strand4)
      mapInd = which(uList$mapping %in% input$mapping4)
      featInd = which(uList$feature %in% input$features4)
      lengthInd = which(uList$widths %in% seq(input$length4[1],input$length4[2]))
      if(! "all" %in% input$genes4){
        geneInd = which(uList$gene %in% input$genes4)
      }
      if("all" %in% input$genes4){
        geneInd = 1:length(uList)
      }
      frameInd = which(uList$frame %in% input$frame4)
      common = Reduce(intersect,list(sampInd,strandInd,mapInd,featInd, geneInd, frameInd, lengthInd))
      uList = uList[common]
      uList2 = split(uList, f = mcols(uList)[which(names(mcols(uList)) == "sample")][,1])
      tab = data.frame()
      print(uList2)
      for (i in 1:length(uList2)){
        tab[i,1] =  mcols(uList2[[i]])[which(names(mcols(uList2[[i]])) == "sample")][1,1]
        print(tab)
        tab[i,2] = length(uList2[[i]])
        print(tab)
        u = table(uList2[[i]]$mapping)
        tab[i,3] = u[which(names(u) == "Unique")] / sum(u)
        print(tab)
        uf = table(uList2[[i]]$feature) 
        if ("cds" %in% names(uf)){
          tab[i,4] = uf[which(names(uf) == "cds")] / sum(uf)
        }
        else{
          tab[i,4] = 0
        }
        if ("rRNA" %in% names(uf)){
          tab[i,5] = uf[which(names(uf) == "rRNA")] / sum(uf)
        }
        else{
          tab[i,5] = 0
        }
        if ("tRNA" %in% names(uf)){
          tab[i,6] = uf[which(names(uf) == "tRNA")] / sum(uf)
        }
        else{
          tab[i,6] = 0
        }
        if ("threeUTR" %in% names(uf)){
          tab[i,7] = uf[which(names(uf) == "threeUTR")] / sum(uf)
        }
        else{
          tab[i,7] = 0
        }
        if ("fiveUTR" %in% names(uf)){
          tab[i,8] = uf[which(names(uf) == "fiveUTR")] / sum(uf)
        }
        else{
          tab[i,8] = 0
        }
        if ("intergenic" %in% names(uf)){
          tab[i,9] = uf[which(names(uf) == "intergenic")] / sum(uf)
        }
        else{
          tab[i,9] = 0
        }
        if ("promoter" %in% names(uf)){
          tab[i,10] = uf[which(names(uf) == "promoter")] / sum(uf)
        }
        else{
          tab[i,10] = 0
        }
        #tab[i,5] = uf[which(names(uf) == "rRNA")] / sum(uf)
        #tab[i,6] = uf[which(names(uf) == "tRNA")] / sum(uf)
        #tab[i,7] = uf[which(names(uf) == "threeUTR")] / sum(uf)
        #tab[i,8] = uf[which(names(uf) == "fiveUTR")] / sum(uf)
        #tab[i,9] = uf[which(names(uf) == "intergenic")] / sum(uf)
        #tab[i,10] = uf[which(names(uf) == "intron")] / sum(uf)
        #tab[i,10] = uf[which(names(uf) == "promoter")] / sum(uf)
        tab[i,11] = length(unique(uList2[[i]])) / length(uList2[[i]])
        ps = vector()
        print(tab)
        t = table(grr[[i]]$widths)
        t = t / length(grr[[i]])
        l = names(which(t > .07))
        print(l)
        # l = 27:30
        for (k in l){
          g = grr[[i]][which(!is.na(grr[[i]]$posA))]
          g = g$frame[which(g$widths == k)]
          tabl = table(g)
          p = max(tabl[1:3]) / sum(tabl[1:3])
          ps = c(ps, p)
        }
        tab[i,12] = max(ps)
        cols = c("sample","Aligned Reads","Unique%", "cds%", "rRNA%", "tRNA%", "3UTR%", "5UTR%", "intergenic%", "promoter%", "complexity", "periodicity")
        colnames(tab) = cols
      }
      #print(tab[,input$cols])
      return(tab[input$cols])
    })
  })
  
  output$table = renderDataTable({
    #plotTab()
    print(plotTab())
  })
  
    
  #Generate metagene plots
  plotMeta = eventReactive(input$submit5,{
    withProgress(message = 'Generating Plot', value = 0, {
      uList = unlist(adjust())
      feat = getFeat()
      #uList = ul3
      mRNA = c("cds", "threeUTR", "fiveUTR")
      sampInd = which(uList$sample %in% input$samples5)
      strandInd = which(uList$strands %in% input$strand5)
      mapInd = which(uList$mapping %in% input$mapping5)
      featInd = which(uList$feature %in% mRNA)
      if(input$genes5 != "all"){geneInd = which(uList$gene %in% input$genes5)}
      if(input$genes5 == "all"){geneInd = 1:length(uList)}
      frameInd = which(uList$frame %in% input$frame5)
      lengthInd = which(uList$widths %in% seq(input$length5[1],input$length5[2]))
      common = Reduce(intersect,list(sampInd,strandInd,mapInd,featInd, geneInd, frameInd, lengthInd))
      uList = uList[common]
      uList$featStart = -1
      uList$featEnd = -1
      uList$dist = -1
      uList$pos = "none"
      fo = findOverlaps(uList, feat$fiveUTR, ignore.strand = T)
      fp = queryHits(fo)
      fs = subjectHits(fo)
      to = findOverlaps(uList, feat$threeUTR, ignore.strand = T)
      tp = queryHits(to)
      ts = subjectHits(to)
      uList$featStart[fp] = start(feat$fiveUTR[fs])
      uList$featEnd[fp] = end(feat$fiveUTR[fs])
      uList$featStart[tp] = start(feat$threeUTR[ts])
      uList$featEnd[tp] = end(feat$threeUTR[ts])
      c = which(uList$feature == "cds")
      t = which(uList$feature == "threeUTR")
      f = which(uList$feature == "fiveUTR")
      uList$dist[t] = start(uList[t]) - (uList$featStart[t] - 1)
      uList$dist[f] = start(uList[f]) - (uList$featEnd[f] + 1)
      pm = pmin(uList$distA[c], uList$distB[c])
      uList$dist[c] = pm
      uList$pos[c][which(pm == uList$distA[c])] = "a"
      uList$pos[c][which(pm == uList$distB[c])] = "b"
      #s = start(uList[c]) - uList$featStart[c]
      #e = start(uList[c]) - uList$featEnd[c]
      
      uList$pos[t] = "b"
      uList$pos[f] = "a"
      uList$dist[which(uList$pos == "b")] = uList$dist[which(uList$pos == "b")] * -1
      
      uList = uList[which(abs(uList$dist) < input$window5)]
      uList = split(uList, f = mcols(uList)[which(names(mcols(uList)) == "pos")][,1])
      plotList = list()
      ttl = c("Start", "End")
      for (i in 1:2){
        uList2 = split(uList[[i]], f = mcols(uList[[i]])[which(names(mcols(uList[[i]])) == input$splitLine5)][,1])
        w = vector()
        n = vector()
        for (k in 1:length(uList2)){
          w = c(w,uList2[[k]]$dist)
          n = c(n,mcols(uList2[[k]])[which(names(mcols(uList2[[k]])) == input$splitLine5)][,1])
        }
        df = data.frame(w,n)
        #gg = ggplot(df, aes(x=w, fill=n)) + geom_histogram( position="identity", binwidth = 1) + labs(title = ttl[i], x = "Position")
        gg = ggplot(df, aes(x=w)) + geom_density(aes(group=n,y = ..scaled.., color = n),bw=1) + labs(title = ttl[i], x = "Position") + scale_color_brewer(palette=input$cols4)
        plotList[[i]] = gg
      }
      print(plot_grid(plotlist = plotList,ncol=floor(i/2)))
    })  
  })
  
  output$meta = renderPlot({
    plotMeta()
  })
  
  plotGview = eventReactive(input$submit6,{
    withProgress(message = 'Generating Plot', value = 0, {
      
      uList = unlist(adjust())
      sampInd = which(uList$sample %in% input$samples6)
      strandInd = which(uList$strands %in% input$strand6)
      mapInd = which(uList$mapping %in% input$mapping6)
      featInd = which(uList$feature  == "cds")
      if(! "all" %in% input$genes6){
        geneInd = which(uList$gene %in% input$genes6)
      }
      if("all" %in% input$genes6){
        geneInd = 1:length(uList)
      }
      frameInd = which(uList$frame %in% input$frame6)
      lengthInd = which(uList$widths %in% seq(input$length6[1],input$length6[2]))
      common = Reduce(intersect,list(sampInd,strandInd,mapInd,featInd, geneInd, frameInd, lengthInd))
      cByg <- cdsBy(getTxdb(), by=c("gene"))#check to make sure this is right
      uList = uList[common]
      
      if (length(input$genes6) > 1){
        uList = split(uList, f = mcols(uList)[which(names(mcols(uList)) == "gene")][,1])
        plotList = list()
        for (i in 1:length(uList)){
          str = unique(uList[[i]]$strands)
          print(str)
          uList2 = split(uList[[i]], f = mcols(uList[[i]])[which(names(mcols(uList[[i]])) == input$splitLine6)][,1])
          w = vector()
          sample = vector()
          ttl = mcols(uList[[i]])[which(names(mcols(uList[[i]])) == "gene")][1,1]
          for (k in 1:length(uList2)){
            w = c(w,uList2[[k]]$distA)
            sample = c(sample,mcols(uList2[[k]])[which(names(mcols(uList2[[k]])) == input$splitLine6)][,1])
          }
          df = data.frame(w,sample)
          gg = ggplot(df, aes(x=w, fill=sample)) + geom_histogram( position="identity", binwidth = 1) + labs(title = ttl, x = "Position") + theme(plot.margin = unit(c(1,1,5,1), "lines"))
          s = start(cByg[ttl])
          e = end(cByg[ttl])
          e = e - min(s)
          s = s - min(s)
          s2 = (s - max(s)) * -1 + (max(e) - max(s))
          e2 = (e - max(s)) * -1 + (max(e) - max(s))
          s2=rev(s2[[1]])
          e2 = rev(e2[[1]])
          print(s)
          print(e)
          if (str == "+"){
            for(l in 1:length(s[[1]])){
              gg = gg + annotation_custom(grob = linesGrob(), xmin = s[[1]][l],xmax = e[[1]][l], ymin = 0, ymax = 0)
            }
          }
          if (str == "-"){
            for(l in 1:length(s2)){
              gg = gg + annotation_custom(grob = linesGrob(), xmin = s2[l],xmax = e2[l], ymin = 0, ymax = 0)
            }
          }
          plotList[[i]] = gg
        }
        print(plot_grid(plotlist = plotList,ncol=floor(i/2)))
      }
      if (length(input$genes6) == 1){
        str = unique(uList$strands)
        print(str)
        if(input$splitPlot6 != "none"){
          if (input$splitLine6 != "none"){
            uList = split(uList, f = mcols(uList)[which(names(mcols(uList)) == input$splitPlot6)][,1])
            plotList = list()
            for (i in 1:length(uList)){
              uList2 = split(uList[[i]], f = mcols(uList[[i]])[which(names(mcols(uList[[i]])) == input$splitLine6)][,1])
              w = vector()
              n = vector()
              ttl2 = mcols(uList[[i]])[which(names(mcols(uList[[i]])) == input$splitPlot6)][1,1]
              ttl = mcols(uList[[i]])[which(names(mcols(uList[[i]])) == "gene")][1,1]
              for (k in 1:length(uList2)){
                w = c(w,uList2[[k]]$distA)
                n = c(n,mcols(uList2[[k]])[which(names(mcols(uList2[[k]])) == input$splitLine6)][,1])
              }
              df = data.frame(w,n)
              df$n = as.character(df$n)
              gg = ggplot(df, aes(x=w, fill=n)) + geom_histogram( position="identity", binwidth = 1) + labs(title = ttl2, x = "Position")
              s = start(cByg[ttl])
              e = end(cByg[ttl])
              e = e - min(s)
              s = s - min(s)
              s2 = (s - max(s)) * -1 + (max(e) - max(s))
              e2 = (e - max(s)) * -1 + (max(e) - max(s))
              s2=rev(s2[[1]])
              e2 = rev(e2[[1]])
              print(s2)
              print(e2)
              if (str == "+"){
                for(l in 1:length(s)){
                  gg = gg + annotation_custom(grob = linesGrob(), xmin = s[[1]][l],xmax = e[[1]][l], ymin = 0, ymax = 0)
                }
              }
              if (str == "-"){
                for(l in 1:length(s2)){
                  gg = gg + annotation_custom(grob = linesGrob(), xmin = s2[l],xmax = e2[l], ymin = 0, ymax = 0)
                }
              }
              plotList[[i]] = gg
            }
            print(plot_grid(plotlist = plotList,ncol=floor(i/2)))
          }
          if (input$splitLine6 == "none"){
            uList = split(uList, f = mcols(uList)[which(names(mcols(uList)) == input$splitPlot6)][,1])
            plotList = list()
            for (i in 1:length(uList)){
              ttl2 = mcols(uList[[i]])[which(names(mcols(uList[[i]])) == input$splitPlot6)][1,1]
              ttl = mcols(uList[[i]])[which(names(mcols(uList[[i]])) == "gene")][1,1]
              w = uList[[i]]$distA
              n = mcols(uList[[i]])[which(names(mcols(uList[[i]])) == input$splitPlot6)][,1]
              df = data.frame(w,n)
              df$n = as.character(df$n)
              gg = ggplot(df, aes(x=w, fill=n)) + geom_histogram( position="identity", binwidth = 1) + labs(title = ttl2, x = "Position")
              s = start(cByg[ttl])
              e = end(cByg[ttl])
              e = e - min(s)
              s = s - min(s)
              s2 = (s - max(s)) * -1 + (max(e) - max(s))
              e2 = (e - max(s)) * -1 + (max(e) - max(s))
              s2=rev(s2[[1]])
              e2 = rev(e2[[1]])
              print(s2)
              print(e2)
              if (str == "+"){
                for(l in 1:length(s)){
                  gg = gg + annotation_custom(grob = linesGrob(), xmin = s[[1]][l],xmax = e[[1]][l], ymin = 0, ymax = 0)
                }
              }
              if (str == "-"){
                for(l in 1:length(s2)){
                  gg = gg + annotation_custom(grob = linesGrob(), xmin = s2[l],xmax = e2[l], ymin = 0, ymax = 0)
                }
              }
              plotList[[i]] = gg
            }  
            print(plot_grid(plotlist = plotList,ncol=floor(i/2)))
            
          }
        }
        if (input$splitPlot6 == "none"){
          if (input$splitLine6 != "none"){
            ttl = input$genes6
            w = uList$distA
            frame = mcols(uList)[which(names(mcols(uList)) == input$splitLine6)][,1]
            df = data.frame(w,frame)
            df$frame = as.character(df$frame)
            gg = ggplot(df, aes(x=w, fill=frame)) + geom_histogram( position="identity", binwidth = 1) + labs(title = ttl, x = "Position")+ theme(plot.margin = unit(c(1,1,5,1), "lines"))
            s = start(cByg[ttl])
            e = end(cByg[ttl])
            e = e - min(s)
            s = s - min(s)
            s2 = (s - max(s)) * -1 + (max(e) - max(s))
            e2 = (e - max(s)) * -1 + (max(e) - max(s))
            s2=rev(s2[[1]])
            e2 = rev(e2[[1]])
            print(s2)
            print(e2)
            if (str == "+"){
              for(l in 1:length(s)){
                gg = gg + annotation_custom(grob = linesGrob(), xmin = s[[1]][l],xmax = e[[1]][l], ymin = 0, ymax = 0)
              }
            }
            if (str == "-"){
              for(l in 1:length(s2)){
                gg = gg + annotation_custom(grob = linesGrob(), xmin = s2[l],xmax = e2[l], ymin = 0, ymax = 0)
              }
            }
            print(gg)
          }
          if (input$splitLine6 == "none"){
            ttl = input$genes6
            w = uList$distA
            df = data.frame(w)
            gg = ggplot(df, aes(x=w)) + geom_histogram( position="identity", binwidth = 1) + labs(title = ttl, x = "Position")
            s = start(cByg[ttl])
            e = end(cByg[ttl])
            e = e - min(s)
            s = s - min(s)
            s2 = (s - max(s)) * -1 + (max(e) - max(s))
            e2 = (e - max(s)) * -1 + (max(e) - max(s))
            s2=rev(s2[[1]])
            e2 = rev(e2[[1]])
            print(s2)
            print(e2)
            if (str == "+"){
              for(l in 1:length(s)){
                gg = gg + annotation_custom(grob = linesGrob(), xmin = s[[1]][l],xmax = e[[1]][l], ymin = 0, ymax = 0)
              }
            }
            if (str == "-"){
              for(l in 1:length(s2)){
                gg = gg + annotation_custom(grob = linesGrob(), xmin = s2[l],xmax = e2[l], ymin = 0, ymax = 0)
              }
            }
            print(gg)
          }
        }
      }
    })
  })
  
  output$gview = renderPlot({
    plotGview()
  })
  plotMetap = eventReactive(input$submit7,{    #uList = unlist(gr)
    withProgress(message = 'Generating Plot', value = 0, {
      
      uList = unlist(adjust())
      feat = getFeat()
      #uList = nr
      mRNA = c("cds", "threeUTR", "fiveUTR")
      sampInd = which(uList$sample %in% input$samples7)
      strandInd = which(uList$strands %in% input$strand7)
      mapInd = which(uList$mapping %in% input$mapping7)
      featInd = which(uList$feature %in% mRNA)
      if(input$genes7 != "all"){geneInd = which(uList$gene %in% input$genes7)}
      if(input$genes7 == "all"){geneInd = 1:length(uList)}
      #frameInd = which(uList$frame %in% input$frame5)
      lengthInd = which(uList$widths %in% seq(input$length7[1],input$length7[2]))
      common = Reduce(intersect,list(sampInd,strandInd,mapInd,featInd, geneInd, lengthInd))
      uList = uList[common]
      uList$featStart = -1
      uList$featEnd = -1
      uList$dist = -1
      uList$pos = "none"
      fo = findOverlaps(uList, feat$fiveUTR, ignore.strand = T)
      fp = queryHits(fo)
      fs = subjectHits(fo)
      to = findOverlaps(uList, feat$threeUTR, ignore.strand = T)
      tp = queryHits(to)
      ts = subjectHits(to)
      uList$featStart[fp] = start(feat$fiveUTR[fs])
      uList$featEnd[fp] = end(feat$fiveUTR[fs])
      uList$featStart[tp] = start(feat$threeUTR[ts])
      uList$featEnd[tp] = end(feat$threeUTR[ts])
      c = which(uList$feature == "cds")
      t = which(uList$feature == "threeUTR")
      f = which(uList$feature == "fiveUTR")
      uList$dist[t] = start(uList[t]) - (uList$featStart[t] - 1)
      uList$dist[f] = start(uList[f]) - (uList$featEnd[f] + 1)
      pm = pmin(uList$distA[c], uList$distB[c])
      uList$dist[c] = pm
      uList$pos[c][which(pm == uList$distA[c])] = "a"
      uList$pos[c][which(pm == uList$distB[c])] = "b"
      #s = start(uList[c]) - uList$featStart[c]
      #e = start(uList[c]) - uList$featEnd[c]
      
      uList$pos[t] = "b"
      uList$pos[f] = "a"
      uList$dist[which(uList$pos == "b")] = uList$dist[which(uList$pos == "b")] * -1
      
      #uList$dist[c[which(abs(e) < abs(s))]] = e[which(abs(e) < abs(s))]
      #uList$dist[c[which(abs(s) <= abs(e))]] = s[which(abs(s) <= abs(e))]
      #uList$pos[c[which(abs(s) <= abs(e))]] = "e"
      #uList$pos[c[which(abs(e) < abs(s))]] = "s"
      uList = uList[which(abs(uList$dist) < input$window7)]
      plotList = list()
      uList = split(uList, f = mcols(uList)[which(names(mcols(uList)) == input$splitPlot7)][,1])
      m = 1
      ttl = c("Start", "End")
      for (l in 1:length(uList)){
        uList2 = split(uList[[l]], f = mcols(uList[[l]])[which(names(mcols(uList[[l]])) == "pos")][,1])
        for (i in 1:2){
          w = uList2[[i]]$dist
          #if( i == 1){ n = as.character(w %%3)}
          #if( i == 2){ n = as.character(abs(w) %%3)}
          frame = as.character(uList2[[i]]$frame)
          df = data.frame(w,frame)
          tit = paste(ttl[i],mcols(uList[[l]])[which(names(mcols(uList[[l]])) == input$splitPlot7)][1,1], sep = " ")
          gg = ggplot(df, aes(x=w, fill=frame)) + geom_histogram( position="identity", binwidth = 1) + labs(title = tit, x = "Position") + scale_color_brewer(palette=input$cols5)
          #gg = ggplot(df, aes(x=w)) + geom_density(aes(group=n,y = ..scaled.., color = n),bw=0.5) + labs(title = ttl[i], x = "Position")
          plotList[[m]] = gg
          m = m + 1
        }
      }
      print(plot_grid(plotlist = plotList,ncol=2))
    })  
  })
  
  output$meta_perio = renderPlot({
    plotMetap()
  })
  
  plotPerio = eventReactive(input$submit8,{
    withProgress(message = 'Generating Plot', value = 0, {
      
      #print((length(input$samples8)))
      uList = unlist(adjust())
      sampInd = which(uList$sample %in% input$samples8)
      strandInd = which(uList$strands %in% input$strand8)
      mapInd = which(uList$mapping %in% input$mapping8)
      featInd = which(uList$feature %in% input$features8)
      lengthInd = which(uList$widths %in% seq(input$length8[1],input$length8[2]))
      if(input$genes8 != "all"){geneInd = which(uList$gene %in% input$genes8)}
      if(input$genes8 == "all"){geneInd = 1:length(uList)}
      frameInd = which(uList$frame %in% input$frame8)
      common = Reduce(intersect,list(sampInd,strandInd,mapInd,featInd, lengthInd, frameInd, geneInd))
      uList = uList[common]
      if(input$splitGroup8 == "none" || length(input$samples8) == 1){
        nuList = split(uList, f = mcols(uList)[which(names(mcols(uList)) == "widths")][,1])
        listy = list()
        for (i in 1:length(nuList)){
          f = mcols(nuList[[i]][,which(names(mcols(nuList[[i]])) == "frame")])[,1]
          x = table(f)/length(f)
          listy[[i]] = x
        }
        m = do.call(rbind, lapply(lapply(listy, unlist), "[",
                                  unique(unlist(c(sapply(listy,names))))))
        colnames(m) <- unique(unlist(c(sapply(listy,names))))
        m[is.na(m)] <- 0
        rownames(m) = sort(unique(mcols(uList[,which(names(mcols(uList)) == "widths")])[,1]))
        dat = expand.grid(rownames(m), colnames(m))
        names(dat) = c("widths", "frame")
        dat$percentage = as.vector(m)
        print(ggplot() + geom_bar(data = dat, aes(y = percentage, x = dat[,which(names(dat) == "widths")], fill = dat[,which(names(dat) == "frame")]), stat = "identity", position = "dodge") + scale_fill_discrete(name = "frame") + xlab("widths") + scale_color_brewer(palette=input$cols6))
        
      }
      if (input$splitGroup8 != "none" && length(input$samples8) > 1 ){
        nuList = split(uList, f = mcols(uList)[which(names(mcols(uList)) == "widths")][,1])
        listy = list()
        m = 1
        for (i in 1:length(nuList)){
          uList2 = split(nuList[[i]], f = mcols(nuList[[i]])[which(names(mcols(nuList[[i]])) == input$splitGroup8)][,1])
          for (k in 1:length(uList2)){
            f = mcols(uList2[[k]][,which(names(mcols(uList2[[k]])) == "frame")])[,1]
            x = table(f)/length(f)
            listy[[m]] = x
            m = m+1
          }
        }
        m2 = do.call(rbind,listy)
        m = do.call(rbind, lapply(lapply(listy, unlist), "[",unique(unlist(c(sapply(listy,names))))))
        colnames(m) <- unique(unlist(c(sapply(listy,names))))
        m[is.na(m)] <- 0
        onenames = unique(mcols(uList[,which(names(mcols(uList)) == "widths")])[,1])
        twonames = colnames(m)
        threenames = sort(unique(mcols(uList[,which(names(mcols(uList)) == input$splitGroup8)])[,1]))
        dat = expand.grid(onenames, twonames, threenames)
        names(dat) = c("widths", "frame", input$splitGroup8)
        d1 = length(threenames)
        d2 = length(onenames) * length(twonames)
        print(m)
        dim(m) = c(d1,d2)
        m = t(m)
        dat$percentage = as.vector(m)
        print(ggplot() + geom_bar(data = dat, aes(y = percentage, x = dat[,which(names(dat) == "widths")], fill = dat[,which(names(dat) == "frame")]),stat = "identity",position = "stack") + theme_bw() + facet_grid(paste("~", input$splitGroup8)) + scale_fill_discrete(name = "frame") + xlab("widths") + scale_color_brewer(palette=input$cols6)) 
        
      }
    })
  })
  
  output$perio <- renderPlot({
    plotPerio()
  })
      
  output$DownloadRld <- downloadHandler(
    filename = "RLD.png",
    content = function(file) {
      png(file)
      print(getrld())
      dev.off()
    })
  output$DownloadComp <- downloadHandler(
    filename = "comp.png",
    content = function(file) {
      png(file)
      print(plotComp())
      dev.off()
    })
  output$DownloadGc <- downloadHandler(
    filename = "gc.png",
    content = function(file) {
      png(file)
      print(plotGc())
      dev.off()
    })
  output$DownloadMeta <- downloadHandler(
    filename = "meta.png",
    content = function(file) {
      png(file)
      print(plotMeta())
      dev.off()
    })
  output$DownloadGview <- downloadHandler(
    filename = "gview.png",
    content = function(file) {
      png(file)
      print(plotGview())
      dev.off()
    })
  output$DownloadMetap <- downloadHandler(
    filename = "meta_perio.png",
    content = function(file) {
      png(file)
      print(plotMetap())
      dev.off()
    })
  output$DownloadPerio <- downloadHandler(
    filename = "perio.png",
    content = function(file) {
      png(file)
      print(plotPerio())
      dev.off()
    })
  
})


