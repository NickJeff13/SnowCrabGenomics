#Plotting treemix results
source("C:/Users/JEFFERYN/Downloads/treemix-1.13/treemix-1.13/src/plotting_funcs.R")

plot_tree2("Treemix_out_pops/snowcrab_lcWGS_100k_filter_pops_out.gz.8")

plot_resid("ZosteraTreemixOuput6",pop_order = "Treemixpoporder.txt")

pop_order = c("MASI","SAC","L3F","SUM","POK","PRJ","SAM","NAH",
              "SEPT","GRB","HEB","PORT", "NRIV","EBAY","POUL","BUCK","MELM","TAYH","PETI","RIM","JB33","JB38","TSW")



###########Modify Treemix plotting function############
plot_tree2 <- function(stem, o = NA, cex = 2, disp = 0.003, plus = 0.01, flip = vector(), 
                       arrow = 0.5, scale = T, ybar = 0.1, mbar = T, plotmig = T, plotnames = T, xmin = 0, lwd = 4, font = 2){
  d = paste(stem, ".vertices.gz", sep = "")
  e = paste(stem, ".edges.gz", sep = "")
  se = paste(stem, ".covse.gz", sep = "")
  d = read.table(gzfile(d), as.is = T, comment.char = "", quote = "")
  e = read.table(gzfile(e), as.is  = T, comment.char = "", quote = "")
  if (!is.na(o)){
    o = read.table(o, as.is = T, comment.char = "", quote = "")
  }
  e[,3] = e[,3]*e[,4]
  e[,3] = e[,3]*e[,4]
  
  se = read.table(gzfile(se), as.is = T, comment.char = "", quote = "")
  m1 = apply(se, 1, mean)
  m = mean(m1)
  #m = 0
  for(i in 1:length(flip)){
    d = flip_node(d, flip[i])
  }
  d$x = "NA"
  d$y = "NA"
  d$ymin = "NA"
  d$ymax = "NA"
  d$x = as.numeric(d$x)
  d$y = as.numeric(d$y)
  d$ymin = as.numeric(d$ymin)
  d$ymax = as.numeric(d$ymax)
  
  d = set_y_coords(d)
  d = set_x_coords(d, e)
  print(d)
  d = set_mig_coords(d, e)
  plot_tree_internal(d, e, o = o, cex = cex, xmin = xmin, disp = disp, plus = plus, arrow = arrow, 
                     ybar = ybar, mbar = mbar, mse = m, scale = scale, plotmig = plotmig, plotnames = plotnames, lwd = lwd, font = font)
  return(list( d= d, e = e))
}

plot_tree_internal<- function(d, e, o = NA, cex = 1, disp = 0.005, plus = 0.005, 
                              arrow = 0.05, ybar = 0.01, scale = T, mbar = T, 
                              mse = 0.01, plotmig = T, plotnames = T, xmin = 0, lwd = 2, font = 1){
  plot(d$x, d$y, axes = F, ylab = "", xlab = "Drift parameter", cex.lab=2, cex.axis=2, xlim = c(xmin, max(d$x)+plus), pch = "") #added cex.lab=2
  axis(1)
  mw = max(e[e[,5]=="MIG",4])
  mcols = rev(heat.colors(150))
  for(i in 1:nrow(e)){
    col = "black"
    if (e[i,5] == "MIG"){
      w = floor(e[i,4]*200)+50
      if (mw > 0.5){
        w = floor(e[i,4]*100)+50
      }
      col = mcols[w]
      if (is.na(col)){
        col = "blue"
      }
    }
    v1 = d[d[,1] == e[i,1],]
    v2 = d[d[,1] == e[i,2],]
    if (e[i,5] == "MIG"){
      if (plotmig){
        arrows( v1[1,]$x, v1[1,]$y, v2[1,]$x, v2[1,]$y, col = col, length = arrow, lwd=3)
      }
    }
    else{
      lines( c(v1[1,]$x, v2[1,]$x), c(v1[1,]$y, v2[1,]$y), col = col, lwd = 3)
    }
  }
  tmp = d[d[,5] == "TIP",]
  print(tmp$x)
  print(disp)
  if ( !is.na(o)){
    for(i in 1:nrow(tmp)){
      tcol = o[o[,1] == tmp[i,2],2]
      if(plotnames){
        #print(tmp[i,2])
        text(tmp[i,]$x+disp, tmp[i,]$y, labels = tmp[i,2], adj = 0, cex = cex, col  = tcol, font = font)
      }
    }
  }
  else{
    if (plotnames){
      text(tmp$x+disp, tmp$y, labels = tmp[,2], adj = 0, cex = cex, font = font)
    }
  }
  if (scale){
    print (paste("mse", mse))
    lines(c(0, mse*50), c(ybar, ybar),lwd=2) #changed mse*10 to mse*40 and added lwd=2
    text( 0, ybar - 0.04, lab = "50 s.e.", adj = 0, cex  = 1.5)
    lines( c(0, 0), c( ybar - 0.01, ybar+0.01))
    lines( c(mse*50, mse*50), c(ybar- 0.01, ybar+ 0.01),lwd=2)
  }
  if (mbar){
    mcols = rev( heat.colors(150) )
    mcols = mcols[50:length(mcols)]
    ymi = ybar+0.15
    yma = ybar+0.35
    l = 0.2
    w = l/100
    xma = max(d$x/20)
    rect( rep(0, 100), ymi+(0:99)*w, rep(xma, 100), ymi+(1:100)*w, col = mcols, border = mcols)
    text(xma+disp, ymi, lab = "0", adj = 0, cex = 1.2)
    if ( mw >0.5){ text(xma+disp, yma, lab = "1", adj = 0, cex = 1.2)}
    else{
      text(xma+disp, yma, lab = "0.5", adj = 0, cex =1.2)
    }
    text(0, yma+0.08, lab = "Migration", adj = 0 , cex = 2) #changed cex=2
    text(0, yma+0.03, lab = "weight", adj = 0 , cex = 2)
  }	
}
