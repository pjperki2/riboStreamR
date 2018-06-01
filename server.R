shinyServer(function(input,output,session){
  options(shiny.maxRequestSize=10000*1024^2)
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
    selectInput("genes6", "Genes", c("all",getGenes()), "AT3G20770", T)
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
      eByg <- exonsBy(txdb, by=c("gene"))
      cByg <- cdsBy(txdb, by=c("gene"))#check to make sure this is right
      for (i in 1:length(ribofiles)){
        ga = readGappedReads(ribofiles[i])
        #gr = GRanges(ga)
        o = findOverlaps(ga, cByg[c("AT3G20770", "AT5G25350")])
        #gr = GRanges(ga)
        m = sample(1:length(ga), 100000)
        #m = unique(c(m, queryHits(o)))
        ga = ga[m]
        gr = GRanges(ga)
        #if (ribonames[i] == "riboseq_root_2_1.bam"){
        #  seqlevels(gr) = c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "ChrC", "ChrM")
        #  seqlengths(gr) = c(30427671, 19698289, 23459830, 18585056, 26975502,   154478,   366924)
        #}
        #print(seqlengths(gr))
        seqlevels(gr) = c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "ChrC", "ChrM")
        gr$sample = ribonames[i]
        gr$experiment = "other"
        if (ribonames[i] %in% rnafiles$name){
          gr$experiment = "RNA"
        }
        if (ribonames[i] %in% rpffiles$name){
          gr$experiment = "RPF"
        }
        gr$strands = as.character(gr@strand)
        gr$widths = qwidth(ga)
        neg = which(gr$strands == "-")
        pos = which(gr$strands == "+")
        end(ranges(gr[pos])) = start(ranges(gr[pos])) + gr[pos]$widths
        start(ranges(gr[neg])) = end(ranges(gr[neg])) - gr[neg]$widths
        sbp = ScanBamParam(what = "mapq")
        mapq = unlist(scanBam(ribofiles[i], param = sbp))
        mapq = mapq[m]
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
        o = findOverlaps(gr, feat$cds_red, ignore.strand = T)
        cdss = feat$cds_red
        gr[queryHits(o)]$gene = cdss[subjectHits(o)]$feature_by
        cb = cdsBy(txdb,"tx", use.names =T)
        tb = transcriptsBy(txdb,"gene")
        #tbo = findOverlaps(tb,gr)
        #cbo = findOverlaps(cb,gr)
        #gr$gene[subjectHits(tbo)]  = names(tb[queryHits(tbo)])
        #gr$transcript = "other"
        #gr$transcript[subjectHits(cbo)]  = names(cb[queryHits(cbo)])
        #tstart = as.vector(start(ranges(cb)))
        #tstart2 = match(gr$transcript, names(tstart))
        #tind = which(!is.na(tstart2))
        tstart = as.vector(start(ranges(cByg)))
        tstart2 = match(gr$gene, names(tstart))
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
        #gr$distA[tind][plus] = (as.vector(start(ranges(gr[tind][plus]))) - min(start(ranges(cb[names(tstart[tstart2[tind][plus]])]))))
        #gr$distA[tind][neg] = (max(end(ranges(cb[names(tstart[tstart2[tind][neg]])]))) - as.vector(end(ranges(gr[tind][neg]))))
        #gr$distA = as.numeric(gr$distA) + 1
        #gr$distB[tind][plus] = (max(end(ranges(cb[names(tstart[tstart2[tind][plus]])]))) - as.vector(start(ranges(gr[tind][plus]))))
        #gr$distB[tind][neg] = (as.vector(end(ranges(gr[tind][neg]))) - min(start(ranges(cb[names(tstart[tstart2[tind][neg]])]))))
        #gr$distB = as.numeric(gr$distB)
        
        gr$distA[tind][plus] = (as.vector(start(ranges(gr[tind][plus]))) - min(start(ranges(cByg[names(tstart[tstart2[tind][plus]])]))))
        gr$distA[tind][neg] = (max(end(ranges(cByg[names(tstart[tstart2[tind][neg]])]))) - as.vector(end(ranges(gr[tind][neg]))))
        gr$distA = as.numeric(gr$distA) + 1
        gr$distB[tind][plus] = (max(end(ranges(cByg[names(tstart[tstart2[tind][plus]])]))) - as.vector(start(ranges(gr[tind][plus]))))
        gr$distB[tind][neg] = (as.vector(end(ranges(gr[tind][neg]))) - min(start(ranges(cByg[names(tstart[tstart2[tind][neg]])]))))
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
          #print("ftab")
          #print(ftab)
          frame = as.numeric(names(which.max(ftab)))
          clen = exp1[which(exp1$widths == i)]
          clen = clen[which(clen$distA < 0)]
          ta = table(clen$distA)
          ta = ta[which(as.numeric(names(ta)) > -20)]
          ta = ta[-1]
          #print("ta")
          #print(ta)
          off1 = names(which.max(ta))
          offframe = as.numeric(off1)%%3
          if(length(frame) > 0 & length(offframe) > 0){
            if (offframe != frame ){
              t2 = which(names(ta) == off1)
              #print(as.numeric(names(ta[t2-1]))%%3)
              print(ta)
              print(t2)
              print(frame)
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
                if (length(t2) > 1){
                  if (as.numeric(names(ta[t2-1]))%%3 == frame){
                    off1 = names(ta[t2-1])
                  }
                }
                if (length(t2) == 1){
                  off1 = 0
                }
              }
              if (t2 == 1){
                if (length(t2) > 1){
                  if (as.numeric(names(ta[t2+1]))%%3 == frame){
                  off1 = names(ta[t2+1])
                  }
                }
                if (length(t2) > 1){
                  if (as.numeric(names(ta[t2+1]))%%3 != frame){
                    off1 = names(ta[t2])
                  }
                }
                if (length(t2) > 1){
                  off1 = 0
                }
              }
            }
          }
          if( length(off1) > 0){
            offset1[i-25] = off1
          }
          if( length(off1) < 1){
            offset1[i-25] = 0
          }
        }
        offsets[[m]] = c(as.numeric(offset1))
        incProgress(1/m, detail = paste("File ", m, sep = ""))
        m = m+1
      }
      #print(offsets)
    })
    offsets
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
      expInd = which(uList$experiment %in% input$experiment)
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
      #print(length(sampInd))
      #print(length(strandInd))
      #print(length(frameInd))
      #print(length(featInd))
      #print(length(geneInd))
      #print(length(mapInd))
      
      common = Reduce(intersect,list(sampInd,strandInd,mapInd,featInd, geneInd, frameInd, expInd))
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
            #gg = ggplot(df, aes(x=w)) + stat_density(geom = "line",aes(group=n,y = ..scaled.., color = n),bw=input$bw) +xlim(input$xmin,input$xmax) + ylim(0,1.1) + labs(title = ttl, x = "Length") +scale_color_brewer(palette=input$cols1)
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
            #gg = ggplot(df, aes(x=w)) + stat_density(geom = "line", aes(y = ..scaled..),bw=input$bw) +xlim(input$xmin,input$xmax) + ylim(0,1.1) + labs(title = ttl, x = "Length") + scale_color_brewer(palette=input$cols1)
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
          #print(ggplot(df, aes(x=w)) + stat_density(geom = "line",aes(group=n,y = ..scaled.., color = n),bw=input$bw) + xlim(input$xmin,input$xmax) + ylim(0,1.1) + labs(title = "Read Lengths", x = "Length") + scale_color_brewer(palette=input$cols1))
          print(ggplot(df, aes(x=w)) + geom_density(aes(group=n,y = ..scaled.., color = n),bw=input$bw) + xlim(input$xmin,input$xmax) + ylim(0,1.1) + labs(title = "Read Lengths", x = "Length") + scale_color_brewer(palette=input$cols1))
        }
        
        if(input$splitLine =="none"){
          print(length(uList))
          w = uList$widths
          df = data.frame(w)
          #print(ggplot(df, aes(x=w)) + stat_density(geom = "line",aes(y = ..scaled..),bw=input$bw) +xlim(input$xmin,input$xmax) + ylim(0,1.1) + labs(title = "Read Lengths", x = "Length") + scale_color_brewer(palette=input$cols1))
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
      expInd = which(uList$experiment %in% input$experiment2)
      strandInd = which(uList$strands %in% input$strand2)
      mapInd = which(uList$mapping %in% input$mapping2)
      featInd = which(uList$feature %in% input$features2)
      if(input$genes2 != "all"){geneInd = which(uList$gene %in% input$genes2)}
      if(input$genes2 == "all"){geneInd = 1:length(uList)}
      frameInd = which(uList$frame %in% input$frame2)
      lengthInd = which(uList$widths %in% seq(input$length[1],input$length[2]))
      common = Reduce(intersect,list(sampInd,strandInd,mapInd,featInd, geneInd, frameInd, lengthInd, expInd))
      uList = uList[common]
      if (length(input$sepFeats2) >0){
        sepInd = which(uList$feature %in% input$sepFeats2)
        v = 1:length(uList)
        v = v[-sepInd]
        if (length(v)>0){
          uList[v]$feature = "other"
        }
      }
      
      #################
      if (input$splitPlotComp != "none"){
        if (input$splitGroup != "none"){
          plotList = list()
          tuList = split(uList, f = mcols(uList)[which(names(mcols(uList)) == input$splitPlotComp)][,1])
          #tuList = split(uList, f = mcols(uList)[which(names(mcols(uList)) == input$splitBar)][,1])
          for(co in 1:length(tuList)){  
            nuList = split(tuList[[co]], f = mcols(tuList[[co]])[which(names(mcols(tuList[[co]])) == input$splitBar)][,1])
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
            gg = ggplot() + geom_bar(data = dat, aes(y = percentage, x = dat[,which(names(dat) == input$splitBar)], fill = dat[,which(names(dat) == input$splitCol)]),stat = "identity",position = "stack") + theme_bw() + facet_grid(paste("~", input$splitGroup)) + scale_fill_discrete(name = input$splitCol) + xlab(input$splitBar) + scale_color_brewer(palette=input$cols2) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
            plotList[[co]] = gg
          }
          print(plot_grid(plotlist = plotList,ncol=1))
        }
        if(input$splitGroup == "none" ){
          plotList = list()
          tuList = split(uList, f = mcols(uList)[which(names(mcols(uList)) == input$splitPlotComp)][,1])
          for(co in 1:length(tuList)){ 
            if (input$splitCol != "none"){
              nuList = split(tuList[[co]], f = mcols(tuList[[co]])[which(names(mcols(tuList[[co]])) == input$splitBar)][,1])
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
              gg = ggplot() + geom_bar(data = dat, aes(y = percentage, x = dat[,which(names(dat) == input$splitBar)], fill = dat[,which(names(dat) == input$splitCol)]), stat = "identity", position = "dodge") + scale_fill_discrete(name = input$splitCol) + xlab(input$splitBar) + scale_color_brewer(palette=input$cols2) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
              plotList[[co]] = gg
            }
            print(plot_grid(plotlist = plotList,ncol=1))
          }
        }
      }
      ######################
      if (input$splitPlotComp == "none"){
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
          print(ggplot() + geom_bar(data = dat, aes(y = percentage, x = dat[,which(names(dat) == input$splitBar)], fill = dat[,which(names(dat) == input$splitCol)]),stat = "identity",position = "stack") + theme_bw() + facet_grid(paste("~", input$splitGroup)) + scale_fill_discrete(name = input$splitCol) + xlab(input$splitBar) + scale_color_brewer(palette=input$cols2) + theme(axis.text.x = element_text(angle = 90, hjust = 1))) 
          
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
            print(ggplot() + geom_bar(data = dat, aes(y = percentage, x = dat[,which(names(dat) == input$splitBar)], fill = dat[,which(names(dat) == input$splitCol)]), stat = "identity", position = "dodge") + scale_fill_discrete(name = input$splitCol) + xlab(input$splitBar) + scale_color_brewer(palette=input$cols2) + theme(axis.text.x = element_text(angle = 90, hjust = 1)))
            
          }
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
      expInd = which(uList$experiment %in% input$experiment3)
      strandInd = which(uList$strands %in% input$strand3)
      mapInd = which(uList$mapping %in% input$mapping3)
      if(input$genes3 != "all"){geneInd = which(uList$gene %in% input$genes3)}
      if(input$genes3 == "all"){geneInd = 1:length(uList)}
      frameInd = which(uList$frame %in% input$frame3)
      featInd = which(uList$feature %in% input$features3)
      lengthInd = which(uList$widths %in% seq(input$length3[1],input$length3[2]))
      common = Reduce(intersect,list(sampInd,strandInd,mapInd,featInd, geneInd,frameInd,lengthInd, expInd))
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
            #gg = ggplot(df, aes(x=w)) + stat_density(geom = "line",aes(group=n,y = ..scaled.., color = n),bw=.01) +xlim(0,1) + ylim(0,1) + labs(title = ttl, x = "GC%") + scale_color_brewer(palette=input$cols3)
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
            #gg = ggplot(df, aes(x=w)) + stat_density(geom = "line",aes(y = ..scaled..),bw=.01) +xlim(0,1) + ylim(0,1) + labs(title = ttl, x = "GC%") + scale_color_brewer(palette=input$cols3)
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
          #print(ggplot(df, aes(x=w)) + stat_density(geom = "line",aes(group=n,y = ..scaled.., color = n),bw=.01) + xlim(0,1) + ylim(0,1) + labs(title = "GC", x = "GC%") + scale_color_brewer(palette=input$cols3))
          print(ggplot(df, aes(x=w)) + geom_density(aes(group=n,y = ..scaled.., color = n),bw=.01) + xlim(0,1) + ylim(0,1) + labs(title = "GC", x = "GC%") + scale_color_brewer(palette=input$cols3))
        }
        
        if(input$splitLine3 =="none"){
          print(length(uList))
          w = uList$gc
          df = data.frame(w)
          #print(ggplot(df, aes(x=w)) + stat_density(geom = "line",aes(y = ..scaled..),bw=.01) + xlim(0,1) + ylim(0,1) + labs(title = "GC", x = "GC%") + scale_color_brewer(palette=input$cols3))
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
      expInd = which(uList$experiment %in% input$experiment4)
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
      common = Reduce(intersect,list(sampInd,strandInd,mapInd,featInd, geneInd, frameInd, lengthInd, expInd))
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
        if(length(which(names(u) == "Unique")) > 0){
          tab[i,3] = round((u[which(names(u) == "Unique")] / sum(u)) * 100, digits = 2)
        }
        if(length(which(names(u) == "Unique")) < 1){
          tab[i,3] = 0.21855341
        }
        print(tab)
        uf = table(uList2[[i]]$feature) 
        if ("cds" %in% names(uf)){
          tab[i,4] = round((uf[which(names(uf) == "cds")] / sum(uf)) * 100, digits = 2)
        }
        else{
          tab[i,4] = 0
        }
        if ("rRNA" %in% names(uf)){
          tab[i,5] = round((uf[which(names(uf) == "rRNA")] / sum(uf)) * 100, digits = 2)
        }
        else{
          tab[i,5] = 0
        }
        if ("tRNA" %in% names(uf)){
          tab[i,6] = round((uf[which(names(uf) == "tRNA")] / sum(uf)) * 100, digits = 2)
        }
        else{
          tab[i,6] = 0
        }
        if ("threeUTR" %in% names(uf)){
          tab[i,7] = round((uf[which(names(uf) == "threeUTR")] / sum(uf)) * 100, digits = 2)
        }
        else{
          tab[i,7] = 0
        }
        if ("fiveUTR" %in% names(uf)){
          tab[i,8] = round((uf[which(names(uf) == "fiveUTR")] / sum(uf)) * 100, digits = 2)
        }
        else{
          tab[i,8] = 0
        }
        if ("intergenic" %in% names(uf)){
          tab[i,9] = round((uf[which(names(uf) == "intergenic")] / sum(uf)) * 100, digits = 2)
        }
        else{
          tab[i,9] = 0
        }
        if ("promoter" %in% names(uf)){
          tab[i,10] = round((uf[which(names(uf) == "promoter")] / sum(uf)) * 100, digits = 2)
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
        tab[i,11] = round((length(unique(uList2[[i]])) / length(uList2[[i]])) * 100, digits = 2)
        ps = vector()
        print(tab)
        t = table(uList2[[i]]$widths)
        t = t / length(uList2[[i]])
        l = names(which(t > .07))
        print(l)
        # l = 27:30
        for (k in l){
          g = uList2[[i]][which(!is.na(uList2[[i]]$posA))]
          g = g$frame[which(g$widths == k)]
          tabl = table(g)
          p = max(tabl[1:3]) / sum(tabl[1:3])
          ps = c(ps, p)
        }
        tab[i,12] = round(max(ps) * 100, digits = 2) 
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
      expInd = which(uList$experiment %in% input$experiment5)
      strandInd = which(uList$strands %in% input$strand5)
      mapInd = which(uList$mapping %in% input$mapping5)
      featInd = which(uList$feature %in% mRNA)
      if(input$genes5 != "all"){geneInd = which(uList$gene %in% input$genes5)}
      if(input$genes5 == "all"){geneInd = 1:length(uList)}
      frameInd = which(uList$frame %in% input$frame5)
      lengthInd = which(uList$widths %in% seq(input$length5[1],input$length5[2]))
      common = Reduce(intersect,list(sampInd,strandInd,mapInd,featInd, geneInd, frameInd, lengthInd, expInd))
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
        gg = ggplot(df, aes(x=w)) + stat_density(geom = "line", aes(group=n,y = ..scaled.., color = n),bw=1) + labs(title = ttl[i], x = "Position") + scale_color_brewer(palette=input$cols4)
        #gg = ggplot(df, aes(x=w)) + geom_density(aes(group=n,y = ..scaled.., color = n),bw=1) + labs(title = ttl[i], x = "Position") + scale_color_brewer(palette=input$cols4)
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
      expInd = which(uList$experiment %in% input$experiment6)
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
      common = Reduce(intersect,list(sampInd,strandInd,mapInd,featInd, geneInd, frameInd, lengthInd, expInd))
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
          len = max(end(cByg[ttl])) - min(start(cByg[ttl]))
          #gg = ggplot(df, aes(x=w, fill=sample)) + stat_density( geom = "line",position="identity", binwidth = 1, bw = .5) + labs(title = ttl, x = "Position") + theme(plot.margin = unit(c(1,1,5,1), "lines")) + xlim(0, len)
          gg = ggplot(df, aes(x=w, fill=sample)) + geom_histogram( position="identity", binwidth = 30) + labs(title = ttl, x = "Position") + theme(plot.margin = unit(c(1,1,5,1), "lines")) + xlim(0, len)
          s = start(cByg[ttl])
          e = end(cByg[ttl])
          e = e - min(s)
          s = s - min(s)
          s2 = (s - max(s)) * -1 + (max(e) - max(s))
          e2 = (e - max(s)) * -1 + (max(e) - max(s))
          s2=rev(s2[[1]])
          e2 = rev(e2[[1]])
          #print(s)
          #print(e)
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
        print(plot_grid(plotlist = plotList,ncol=2))#floor(i/2)))
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
              #gg = ggplot(df, aes(x=w, fill=n)) + stat_density(geom = "line", position="identity", binwidth = 1, bw = .5) + labs(title = ttl2, x = "Position")
              gg = ggplot(df, aes(x=w, fill=n)) + geom_histogram( position="identity", binwidth = 30) + labs(title = ttl2, x = "Position")
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
            print(plot_grid(plotlist = plotList,ncol=2))#floor(i/2)))
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
              #gg = ggplot(df, aes(x=w, fill=n)) + stat_density(geom = "line", position="identity", binwidth = 1, bw = .5) + labs(title = ttl2, x = "Position")
              gg = ggplot(df, aes(x=w, fill=n)) + geom_histogram( position="identity", binwidth = 30) + labs(title = ttl2, x = "Position")
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
            #gg = ggplot(df, aes(x=w, fill=frame)) + stat_density( geom = "line",position="identity", binwidth = 1, bw = .5) + labs(title = ttl, x = "Position")+ theme(plot.margin = unit(c(1,1,5,1), "lines"))
            gg = ggplot(df, aes(x=w, fill=frame)) + geom_histogram( position="identity", binwidth = 30) + labs(title = ttl, x = "Position")+ theme(plot.margin = unit(c(1,1,5,1), "lines"))
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
            #df2 = data.frame(count(w))
            #gg = ggplot(df, aes(x=w)) + stat_density(geom = "line", position="identity", bw = 0.5, color = "red") + labs(title = ttl, x = "Position") 
            gg = ggplot(df, aes(x=w, fill=frame)) + geom_histogram( position="identity", binwidth = 30) + labs(title = ttl, x = "Position")+ theme(plot.margin = unit(c(1,1,5,1), "lines"))
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
      expInd = which(uList$experiment %in% input$experiment7)
      strandInd = which(uList$strands %in% input$strand7)
      mapInd = which(uList$mapping %in% input$mapping7)
      featInd = which(uList$feature %in% mRNA)
      if(input$genes7 != "all"){geneInd = which(uList$gene %in% input$genes7)}
      if(input$genes7 == "all"){geneInd = 1:length(uList)}
      #frameInd = which(uList$frame %in% input$frame5)
      lengthInd = which(uList$widths %in% seq(input$length7[1],input$length7[2]))
      common = Reduce(intersect,list(sampInd,strandInd,mapInd,featInd, geneInd, lengthInd, expInd))
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
          uList2[[i]]$frame = uList2[[i]]$dist %% 3
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
      expInd = which(uList$experiment %in% input$experiment8)
      strandInd = which(uList$strands %in% input$strand8)
      mapInd = which(uList$mapping %in% input$mapping8)
      featInd = which(uList$feature %in% input$features8)
      lengthInd = which(uList$widths %in% seq(input$length8[1],input$length8[2]))
      if(input$genes8 != "all"){geneInd = which(uList$gene %in% input$genes8)}
      if(input$genes8 == "all"){geneInd = 1:length(uList)}
      frameInd = which(uList$frame %in% input$frame8)
      common = Reduce(intersect,list(sampInd,strandInd,mapInd,featInd, lengthInd, frameInd, geneInd, expInd))
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
    nums1 <- length(input$expers)
    nums2 <- as.numeric(input$treats)
    if (nums1 == 1){
      lapply(1:nums2, function(i) {
        #selectInput(paste("Treat_", i, sep = ""),paste("Treatment ", i, sep = "" ), choices = getSamps(), multiple = T)
        ex = "RNA"
        if(input$expers == "Ribo-seq"){
          ex = "RPF"
        }
        fileInput(paste(ex, i, sep = ""),paste(ex," Treatment ", i, sep = "" ), multiple = T)
        
      })
    }
  })
  output$inputvals2 <- renderUI({
    print(names(input))
    nums1 <- length(input$expers)
    nums2 <- as.numeric(input$treats)
    n = nums1*nums2
    n2 = rep(c(1:nums2), times = 2)
    nams = rep(c("RNA", "RPF"), each = nums2)
    if(nums1 > 1){
      lapply(1:n, function(i) {
        #selectInput(paste(nams[i],n2[i],sep = ""),paste(nams[i], " Treat ",n2[i] , sep = "" ), choices = getSamps(), multiple = T)
        fileInput(paste(nams[i],n2[i],sep = ""),paste(nams[i]," Treatment ",n2[i] , sep = "" ),  multiple = T)
      })
    }
    
  })

  getEByg = reactive({
    txdb = getTxdb()
    eByg <- exonsBy(txdb, by=c("gene"))#check to make sure this is right
    feat = getFeat()
    nams = feat[[1]]$feature_by
    eByg = eByg[nams]
    eByg = reduce(eByg)
    
  })
  getCByg = reactive({
    txdb = getTxdb()
    cByg <- cdsBy(txdb, by=c("gene"))#check to make sure this is right
    feat = getFeat()
    nams = feat[[1]]$feature_by
    cByg = cByg[nams]
    cByg = reduce(cByg)
    
  })
  
  getRNAcounts = reactive({
    nums = 0
    print(nrow(input$RNA1))
    for(i in 1:input$treats){
      nums = nums + nrow(input[[paste("RNA",i,sep="")]])
    }
    print(nums)
    RNA_list = c(paste("RNA", 1:nums, sep = ""))
    print(RNA_list)
    eByg <- getEByg()
    RNA_ribofiles = vector()
    for (i in 1:length(RNA_list)){
      RNA_ribofiles = c(RNA_ribofiles, input[[RNA_list[i]]]$datapath)
    }
    #RNA_ribofiles <- c(input[[RNA_list]])
    # RNA_ribofiles <- RNA_ribofiles$datapath
    RNA_bfl = BamFileList(RNA_ribofiles)
    names(RNA_bfl) = RNA_list
    RNA_raw_count <- summarizeOverlaps(eByg, RNA_bfl,"Union", ignore.strand=FALSE,inter.feature = F,BPPARAM=SerialParam())
    RNA_count = assays(RNA_raw_count)$counts
    
  })
  
  getRPFcounts = reactive({
    nums = 0
    for(i in 1:input$treats){
      nums = nums + nrow(input[[paste("RPF",i,sep="")]])
    }
    print(nums)
    RPF_list = c(paste("RPF", 1:nums, sep = ""))
    cByg <- getCByg()
    RPF_ribofiles = vector()
    for (i in 1:length(RPF_list)){
      RPF_ribofiles = c(RPF_ribofiles, input[[RPF_list[i]]]$datapath)
    }
    #RPF_ribofiles <- c(get(RPF_list))
    #RPF_ribofiles <- RPF_ribofiles$datapath
    RPF_bfl = BamFileList(RPF_ribofiles)
    names(RPF_bfl) = RPF_list
    RPF_raw_count <- summarizeOverlaps(cByg, RPF_bfl,"Union", ignore.strand=FALSE,inter.feature = FALSE,BPPARAM=SerialParam())
    RPF_count = assays(RPF_raw_count)$counts
    
  })
  
  getRNArpkm = reactive({
    nums = vector()
    print(nrow(input$RNA1))
    for(i in 1:input$treats){
      nums = c(nums, nrow(input[[paste("RNA",i,sep="")]]))
    }
    groups = vector()
    for (i in 1:input$treats){
      groups = c(groups, rep(i, nums[i]))
    }
    print(groups)
    #indRNA = which(substr(names(input),1,3) == "RNA")
    #indRPF = which(substr(names(input),1,3) == "RPF")
    #groups = c(as.numeric(substr(names(input)[indRNA],4,4)))
    RNA_l = sum(width(getEByg()))
    print(RNA_l)
    RNA_dgel = DGEList(getRNAcounts(), group = groups)
    RNA_rpkm = rpkm(RNA_dgel, RNA_l)
    return(RNA_rpkm)
    
  })
  
  getRPFrpkm = reactive({
    nums = vector()
    print(nrow(input$RPF1))
    for(i in 1:input$treats){
      nums = c(nums, nrow(input[[paste("RPF",i,sep="")]]))
    }
    groups = vector()
    for (i in 1:input$treats){
      groups = c(groups, rep(i, nums[i]))
    }
    RPF_l = sum(width(getCByg()))
    RPF_dgel = DGEList(getRPFcounts(), group = groups)
    RPF_rpkm = rpkm(RPF_dgel, RPF_l)
    return(RPF_rpkm)
    
  })
  
  getRPKMs = reactive({
    if(length(input$expers) == 2){
      total_rpkm = cbind(getRNArpkm(), getRPFrpkm())
    }
    if(length(input$expers) == 1){
      if(input$expers == "RNA-seq"){
        total_rpkm = getRNArpkm()
      }
      if(input$expers == "Ribo-seq"){
        total_rpkm = getRPFrpkm()
      }
    }
    return(total_rpkm)
  })
  
  getCounts = reactive({
    print(input$expers)
    if(length(input$expers) == 2){
      total_count = cbind(getRNAcounts(), getRPFcounts())
    }
    if(length(input$expers) == 1){
      print("t2")
      if(input$expers == "RNA-seq"){
        print("t3")
        total_count = getRNAcounts()
      }
      if(input$expers == "Ribo-seq"){
        total_count = getRPFcounts()
      }
    }
    print(head(total_count))
    return(total_count)
  })
  
  doEdger = reactive({
    x1 = as.data.frame(getCounts())
    if(length(input$expers) == 1){
      group = vector()
      ex = "RNA"
      if(input$expers == "Ribo-seq"){
        ex = "RPF"
      }
      for (i in 1:input$treats){
        l = nrow(input[[paste(ex,i,sep = "")]])
        group = c(group, rep(i, l))
      }
      print(group)
      group <- factor(group)
      y <- DGEList(counts=x1,group=group)
      y <- calcNormFactors(y)
      design <- model.matrix(~group)
      y <- estimateDisp(y,design)
      fit <- glmQLFit(y,design)
      qlf <- glmQLFTest(fit,coef=2)
      topTags(qlf)
      fit <- glmFit(y,design)
      lrt <- glmLRT(fit,coef=2)
      #topTags(lrt)
      tt = as.data.frame(topTags(lrt, n = nrow(glmLRT)))
      return(tt)
    }
    if (length(input$expers) == 2){
      Exp = vector()
      numRNA = vector()
      for(i in 1:input$treats){
        numRNA = c(numRNA, nrow(input[[paste("RNA",i,sep="")]]))
      }
      numRPF = vector()
      for(i in 1:input$treats){
        numRPF = c(numRPF, nrow(input[[paste("RPF",i,sep="")]]))
      }
      witches = paste(rep(c("RNA", "RPF"),each = 2), rep(1:2, times = 2), sep = "")
      ind = which(names(input) %in% witches)
      Exp = c(rep("RNA", sum(numRNA)), rep('RPF', sum(numRPF)))
      print(ind)
      Treat = vector()
      for ( i in 1:length(ind)){
        Treat = c(Treat, rep(i, nrow(input[[names(input)[ind[i]]]])))
      }
      targets1 = data.frame(Exp, Treat)
      samps1 = c(paste("Sample_", 1:nrow(targets1), sep = ""))
      rownames(targets1) = samps1
      Group1 <- factor(paste(targets1$Exp,targets1$Treat,sep="."))
      targets1 = cbind(targets1,Group=Group1)
      #targets1$Treat <- relevel(targets1$Treat, ref="treat1")
      #targets1$Exp <- relevel(targets1$Exp, ref="RNA")
      group = factor(Treat)
      y1 <- DGEList(counts=x1,group=group)
      keep <- rowSums(cpm(y1)>1) >= 3
      y1 <- y1[keep, keep.lib.sizes=FALSE]
      y1 <- calcNormFactors(y1)
      design1 <- model.matrix(~Treat + Exp + Treat:Exp, data=targets1)
      y1 <- estimateDisp(y1,design1)
      fit1 <- glmFit(y1, design1)
      lrt1 <- glmLRT(fit1, coef=4)
      tt = as.data.frame(topTags(lrt1, n = nrow(glmLRT)))
      return(tt)
    } 
  })
  bigTab = reactive({
    tab = as.data.frame(getRPKMs())
    tab$RNA_treat1_ave = rowMeans(tab[,1:2])
    tab$RNA_treat2_ave = rowMeans(tab[,3:4])
    tab$RPF_treat1_ave = rowMeans(tab[,5:6])
    tab$RPF_treat2_ave = rowMeans(tab[,7:8])
    
    tab$RNA_treat1_sd = apply(tab[,1:2],1,sd,na.rm = T)
    tab$RNA_treat2_sd = apply(tab[,3:4],1,sd,na.rm = T)
    tab$RPF_treat1_sd = apply(tab[,5:6],1,sd,na.rm = T)
    tab$RPF_treat2_sd = apply(tab[,7:8],1,sd,na.rm = T)
    
    tab$treat1_TE_ave = tab$RPF_treat1_ave / tab$RNA_treat1_ave
    tab$treat2_TE_ave = tab$RPF_treat2_ave / tab$RNA_treat2_ave
    tab$foldChange_ave = (tab$treat2_TE_ave - tab$treat1_TE_ave) / tab$treat1_TE_ave 
    rownames(tab) = names(getEByg())
    tab
  })
  
  countTable = eventReactive(input$submit9,{
    #datatable( doEdger(), rownames = TRUE )
    er = doEdger()
    rpkms = as.data.frame(getRPKMs())
    #rpkms = bigTab()
    er = er[order(row.names(er)),]
    rpkms = rpkms[order(row.names(rpkms)),]
    rpkms = rpkms[which(row.names(rpkms)%in% row.names(er)),]
    print(head(rpkms))
    print(head(er))
    table = cbind(rpkms,er)
    table = table[order(table$FDR),]
    #print(table)
    return(table)
  })
  
  output$countTable2 = renderDataTable({
    #plotTab()
    ct = cbind(rownames(countTable()), countTable())
    colnames(ct)[1] = "Gene"
    print(ct)
  })
  
  
  output$gen = downloadHandler(
    filename = "markdown_test.html",
    content = function(file){
      rnafiles = input$RNAfile
      rpffiles = input$RPFfile
      RNAf1 = input$RNA1
      RNAf2 = input$RNA2
      RPFf1 = input$RPF1
      RPFf2 = input$RPF2
      ribonames = c(rnafiles$name, rpffiles$name)
      params= list(RNAfile = input$RNAfile, RPFfile = input$RPFfile, ribofiles = c(rnafiles$datapath,rpffiles$datapath), ribonames = c(rnafiles$name, rpffiles$name), rpfind = which(ribonames %in% rpffiles$name), rnaind = which(ribonames %in% rnafiles$name), plots = input$plots, samples = input$samples, strand = input$strand, mapping = input$mapping, features = input$features, sepFeats = input$sepFeats, genes = input$genes, frame = input$frame, splitLine = input$splitLine, splitPlot = input$splitPlot, bw = input$bw, xmin = input$xmin, xmax = input$xmax, 
                   samples2 = input$samples2, strand2 = input$strand2, mapping2 = input$mapping2, features2 = input$features2, sepFeats2 = input$sepFeats2, genes2 = input$genes2, frame2 = input$frame2, splitBar = input$splitBar, splitCol = input$splitCol, length = input$length, splitGroup = input$splitGroup,
                   samples7 = input$samples7, strand7 = input$strand7, mapping7 = input$mapping7, genes7 = input$genes7, length7 = input$length7, splitPlot7 = input$splitPlot7, window7 = input$window7,cols5 = input$cols5,
                   samples6 = input$samples6, strand6 = input$strand6, mapping6 = input$mapping6, genes6 = input$genes6, frame6 = input$frame6, length6 = input$length6, splitLine6 = input$splitLine6, splitPlot6 = input$splitPlot6,
                   samples4 = input$samples4, strand4 = input$strand4, mapping4 = input$mapping4, genes4 = input$genes4, frame4 = input$frame4, length4 = input$length4, features4 = input$features4, cols = input$cols,
                   samples3 = input$samples3, strand3 = input$strand3, mapping3 = input$mapping3, genes3 = input$genes3, frame3 = input$frame3, length3 = input$length3, features3 = input$features3, sepFeats3 = input$sepFeats3, splitLine3 = input$splitLine3, splitPlot3 = input$splitPlot3, cols3 = input$cols3,
                   samples5 = input$samples5, strand5 = input$strand5, mapping5 = input$mapping5, genes5 = input$genes5, frame5 = input$frame5, length5 = input$length5, window5 = input$window5, splitLine5 = input$splitLine5, cols4 = input$cols4,
                   samples8 = input$samples8, strand8 = input$strand8, mapping8 = input$mapping8, genes8 = input$genes8, frame8 = input$frame8, length8 = input$length8, splitGroup8 = input$splitGroup8, cols6 = input$cols6, features8 = input$features8,
                   select = input$select, expers = input$expers, treats = input$treats, RNA1 = input$RNA1, RNA2 = input$RNA2, RPF1 = input$RPF1, RPF2 = input$RPF2, RNAfiles2 = c(RNAf1$datapath, RNAf2$datapath), RPFfiles2 = c(RPFf1$datapath, RPFf2$datapath),
                   experiment = input$experiment,experiment2 = input$experiment2,experiment3 = input$experiment3,experiment4 = input$experiment4,experiment5 = input$experiment5,experiment6 = input$experiment6,experiment7 = input$experiment7,experiment8 = input$experiment8)
      tempReport = file.path(getwd(), "markdown_test_4.Rmd")
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
      
    }
  )
      
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


