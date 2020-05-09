library(shiny)
library(dplyr)
library(treemap)
library(dynamicTreeCut)
library(fastcluster)
library(ggplot2)
library(WGCNA)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  Fileread <- reactive({
    #read csv
    infile <- input$file
    if(is.null(infile)) return(NULL)
    read.csv(infile$datapath)
  })
  
  Rdataread <- reactive({
    #read rdata
    infile <- input$rdata
    load(infile$datapath)
  })
  
  queryfun <- reactive({
    #gene--mudule query
    geneData <- Fileread()
    gene_Name <- input$Target_Name
    tmp <- geneData[which(geneData$Gene==gene_Name),"Module"]
    res <- geneData[geneData$Module==tmp,]
    return(res)
  })
  
  TreeMap <- reactive({
    geneData <- Fileread()
    result <- table(cut(geneData$Module,c(0,10,20,30,40,50,60,70,80)))
    result_df <- data.frame(result)
    colnames(result_df) <- c("Module","Count")
    treemap(result_df,index = "Module",vSize = "Count",title = "")
  })
  
  outlierdetect <- reactive({
    Rdataread()
    gsg <- goodSamplesGenes(datExpr, verbose = 3)
    gsg$allOK
    sampleTree = hclust(dist(datExpr), method = "average")
    #sizeGrWindow(12,9)
    #par(cex = 0.6);
    #par(mar = c(0,4,2,0))
    plot(sampleTree, sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
  })
  
  sftthreshold <- reactive({
    Rdataread()
    powers = c(c(1:10), seq(from = 12, to=40, by=2))
    sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
    #sizeGrWindow(9, 5)
    par(mfrow = c(1,2));
    cex1 = 0.9;
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
         main = paste("Scale independence"));
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         labels=powers,cex=cex1,col="red");
    #abline(h=0.90,col="red")
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
         xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
         main = paste("Mean connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  })
  
  moduleindent <- reactive({
    Rdataread()
    softPower = 14
    #adjacency = adjacency(datExpr, power = softPower,type = 'unsigned')
    #TOM = TOMsimilarity(adjacency);
    #dissTOM = 1-TOM
    geneTree = hclust(as.dist(dissTOM), method = "average");
    #sizeGrWindow(12,9)
    #plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
    #     labels = FALSE, hang = 0.04)
    minModuleSize = 30;
    # Module identification using dynamic tree cut:
    dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                deepSplit = 2, pamRespectsDendro = FALSE,
                                minClusterSize = minModuleSize);
    table(dynamicMods)
    dynamicColors = labels2colors(dynamicMods)
    table(dynamicColors)
    # Plot the dendrogram and colors underneath
    #sizeGrWindow(8,6)
    plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05,
                        main = "Gene dendrogram and module colors")
  })
  
  eigene <- reactive({
    MEList = moduleEigengenes(datExpr, colors = dynamicColors)
    MEs = MEList$eigengenes
    # Calculate dissimilarity of module eigengenes
    MEDiss = 1-cor(MEs);
    # Cluster module eigengenes
    METree = hclust(as.dist(MEDiss), method = "average");
    # Plot the result
    #sizeGrWindow(7, 6)
    #plot(METree, main = "Clustering of module eigengenes",
    #     xlab = "", sub = "")
    MEDissThres = 0.05
    # Plot the cut line into the dendrogram
    #abline(h=MEDissThres, col = "red")
    # Call an automatic merging function
    merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
    # The merged module colors
    mergedColors = merge$colors;
    # Eigengenes of the new merged modules:
    mergedMEs = merge$newMEs
    #pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
    plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                        c("Dynamic Tree Cut", "Merged dynamic"),
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05)
  })
  
  Heatmap <- reactive({
    moduleColors = mergedColors
    # Construct numerical labels corresponding to the colors
    colorOrder = c("grey", standardColors(50));
    moduleLabels = match(moduleColors, colorOrder)-1;
    MEs = mergedMEs;
    
    table(mergedColors)
    plotEigengeneNetworks(MEs, 
                          "Eigengene adjacency heatmap", 
                          marHeatmap = c(3,4,2,2), 
                          plotDendrograms = FALSE, 
                          xLabelsAngle = 90) 
  })
  
  genefre <- reactive({
    module_data <- as.data.frame(table(mergedColors))
    ggplot(module_data)+
      geom_bar(aes(x=mergedColors,y=Freq),
               stat="identity",
               width=0.8)+
      theme_bw(base_size=15)+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.position = 'none',
            axis.text.x = element_text(angle = 50,hjust = 1,size=12))+
      geom_text(aes(x=mergedColors,y=Freq,label = paste(module_data$Freq)),
                size=2.2,vjust = -1,
                hjust = 0.5,
                position=position_dodge(0.8))+
      scale_fill_manual(values=rainbow(3))+
      # ylim(0,1000)+
      labs(y='Gene Number',x='')
  })
  
  output$contents <- DT::renderDataTable(DT::datatable({
    
    dat=Fileread()
    
    if(input$disp == "Head") {
      return(head(dat))
    }
    else {
      return(dat)}
  },rownames = FALSE))
    
  output$result <- DT::renderDataTable(DT::datatable({
    queryfun()
  },rownames = FALSE))
  
  output$download <- downloadHandler(
    filename = "query_result.csv",
    content =function(file){
      write.csv(queryfun(),file)
    })
  
  output$Treemap <- renderPlot({
    TreeMap()
  })
  
  output$outlier <- renderPlot({
    outlierdetect()
  })
  
  output$sft <- renderPlot({
    sftthreshold()
  })
  
  output$mod <- renderPlot({
    eigene()
  })
  
  output$heatmap <- renderPlot({
    Heatmap()
  })
  
  output$genefre <- renderPlot({
    genefre()
  })
})