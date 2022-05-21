###-*- coding: utf-8 -*-

#Copyright (C) Li Tao(litao0401@163.com)

#Description:
#********************************************************************************
#                                                                               *
#                                                                               *
#********************************************************************************

#Change Log:
#********************************************************************************
#(1).2020-01-20: Created and finished the first version of this program.        *
#                                                                               *
#                                                                               *
#********************************************************************************

version <- "0.0.1"
#--------------------------------------------------------------------------------
#得到网络空间布局
getpoints <- function(R,ox,oy,seq){
  point_x <- c()
  point_y <- c()
  radians <- (pi/180) * round(360/length(seq))
  for(i in 1:length(seq)){
    x = ox + R * sin(radians * i)
    y = ox + R * cos(radians * i)
    point_x <- append(point_x, x)
    point_y <- append(point_y, y)
  }
  points <- data.frame("point" = seq, "point_x" = point_x, "point_y" = point_y)
  return(points)
}
#--------------------------------------------------------------------------------
#按时间节点分割原始序列
cutdf <- function(df,start,end,len,step){
  times_df_list <- list()
  for(i in seq(start, end, step)){
    if(i >= len){
      times_df_list[[df$time[i]]] <- df[(i - len + 1):i, 2:length(df)]
    }else{
      times_df_list[[df$time[i]]] <- df[-(i+1):-(dim(df)[1]-len + i), 2:length(df)]
    }
  }
  return(times_df_list)
}
#--------------------------------------------------------------------------------
#计算数据的转移熵
caculate_TransferEntropy <- function(df){
  #install.packages("RTransferEntropy")
  library(RTransferEntropy)
  library(future)
  plan(multiprocess)
  x <- c()
  y <- c()
  x.y <- c()
  y.x <- c()
  x.y.p_value <- c()
  y.x.p_value <- c()
  for(i in 1:length(df)){
  #for(i in 1:1){
    for(j in 1:length(df)){
      if (i == j){
        next
      }
      x <- append(x, names(df)[i])
      y <- append(y, names(df)[j])
      transfer_entropy <- transfer_entropy(df[[i]], df[[j]], nboot = 100, quiet = T)
      x.y <- append(x.y, transfer_entropy[["coef"]][1,1])
      y.x <- append(y.x, transfer_entropy[["coef"]][2,1])
      x.y.p_value <- append(x.y.p_value, transfer_entropy[["coef"]][1,4])
      y.x.p_value <- append(y.x.p_value, transfer_entropy[["coef"]][2,4])
    }
  }
  TransferEntropy_df <- data.frame(x, y, x.y, y.x, x.y.p_value, y.x.p_value)
  return(TransferEntropy_df)
}
#--------------------------------------------------------------------------------
#得到单个基因对其他基因的转移熵
get_gene_TransferEntropy_df <- function(TransferEntropy_df, gene){
  #names(TransferEntropy_df) <- c("x", "y", "x.y", "y.x")
  gene_TransferEntropy_df <- subset(TransferEntropy_df, x == gene)
  return(gene_TransferEntropy_df)
}
#--------------------------------------------------------------------------------
#得到每个基因对除去自己以外基因的转移熵
get_genes_TransferEntropy_df <- function(TransferEntropy_df){
  genes_TransferEntropy_list <- list()
  genes <- unique(as.character(TransferEntropy_df[,1]))
  #genes_TransferEntropy_list <- lapply(TransferEntropy_df, get_gene_TransferEntropy_df, genes)
  #names(genes_TransferEntropy_list) <- genes
  for (i in 1:length(genes)){
    genes_TransferEntropy_list[[genes[i]]] <- get_gene_TransferEntropy_df(TransferEntropy_df, genes[i])
  }
  return(genes_TransferEntropy_list)
}
#--------------------------------------------------------------------------------
#得到单个基因对其他基因的网络图
get_gene_network_picture <- function(TransferEntropy_df, gene){
  R = 150
  ox = 250
  oy = 250
  library(ggplot2)
  windowsFonts(CA = windowsFont("Calibri"))
  points <- getpoints(R,ox,oy,unique(c(as.character(TransferEntropy_df[,1]),as.character(TransferEntropy_df[,2]))))
  gene_TransferEntropy_df <- get_gene_TransferEntropy_df(TransferEntropy_df, gene)
  data <- gene_TransferEntropy_df
  data$ratio <- data$x.y/data$y.x
  data$information <- data$x.y + data$y.x
  data$dirction_from <- ifelse(data$x.y > data$y.x, as.character(data$x), as.character(data$y))
  data$dirction_to <- ifelse(data$x.y > data$y.x, as.character(data$y), as.character(data$x))
  a <- merge(data, points,  by.x = "x", by.y = "point", all.x = T)
  names(a)[-(1:(length(a)-2))] <- c("x_point_x","x_point_y")
  b <- merge(a, points,  by.x = "y", by.y = "point", all.x = T )
  names(b)[-(1:(length(b)-2))] <- c("y_point_x","y_point_y")
  c <- merge(b, points, by.x = "dirction_from", by.y = "point", all.x = T)
  names(c)[-(1:(length(c)-2))] <- c("from_point_x","from_point_y")
  d <- merge(c, points, by.x = "dirction_to", by.y = "point", all.x = T)
  names(d)[-(1:(length(d)-2))] <- c("to_point_x","to_point_y")
  df <- subset(d, select = c('x','y','ratio','information','x_point_x', 'x_point_y', 'y_point_x', 'y_point_y', 'from_point_x', 'from_point_y', "to_point_x", "to_point_y" ))
  df$ratio_point_x <- (df$ratio * df$y_point_x + df$x_point_x)/(1 + df$ratio)
  df$ratio_point_y <- (df$ratio * df$y_point_y + df$x_point_y)/(1 + df$ratio)
  p <- ggplot() 
  p_dian <- p + geom_point(aes(x = point_x, y = point_y), points, size = 11, pch = 1) + 
                geom_text(aes(x = point_x, y = point_y + 15, label = point), subset(points, point_y > 250), family = "CA") + 
                geom_text(aes(x = point_x, y = point_y - 15, label = point), subset(points, point_y <= 250), family = "CA")
  p_red <- p_dian + geom_segment(aes(x = x_point_x, y = x_point_y, xend = y_point_x, yend = y_point_y, lwd = information), data = df, col = "red2", show.legend = F)
  p_blue <- p_red + geom_segment(aes(x = x_point_x, y = x_point_y, xend = ratio_point_x, yend = ratio_point_y, lwd = information), data = df, col = "#3399FF", show.legend = F)
  p_arrow <- p_blue + geom_segment(aes( x = from_point_x, y = from_point_y,  xend = ratio_point_x, yend = ratio_point_y), df, lwd = 0, lty = 0, arrow = arrow(length = unit(20, "points"),type = "closed"), alpha = 0.8)
  gene_network_picture <- p_arrow
  return(gene_network_picture)
}
#--------------------------------------------------------------------------------
#得到每个基因对除去自己以外基因的网络图
get_genes_network_picture <- function(TransferEntropy_df, R, ox, oy){
  genes_network_picture_list <- list()
  genes <- unique(as.charactor(TransferEntropy_df[,1]))
  for (i in 1:length(genes)){
    genes_network_picture_list[[genes[i]]] <- get_gene_network_picture(TransferEntropy_df, genes[i], R, ox, oy)
  }
  return(genes_network_picture_list)
}
#--------------------------------------------------------------------------------
#得到单个基因对其他基因的监控图谱数据
get_gene_monitor_data <- function(TransferEntropy_df, gene){
  gene_TransferEntropy_df <- get_gene_TransferEntropy_df(TransferEntropy_df, gene)
  gene_monitor_data <- gene_TransferEntropy_df
  return(gene_monitor_data)
}
#--------------------------------------------------------------------------------
#得到每个基因对其他基因的监控图谱数据
get_genes_monitor_data <- function(TransferEntropy_df){
  genes_monitor_data_list <- list()
  genes <- unique(as.charactor(TransferEntropy_df[,1]))
  for (i in 1:length(genes)){
    genes_monitor_data_list[[genes[i]]] <- get_gene_monitor_data(TransferEntropy_df, gene[i])
  }
  return(genes_monitor_data_list)
}
#--------------------------------------------------------------------------------
#实现单个基因对其他基因的所有功能
#gene_allfunction <- function(TransferEntropy_df, gene){
  #cell <- NULL
  #R = 150
  #ox = 250
  #oy = 250
  #gene_network_picture <- get_gene_network_picture(TransferEntropy_df, gene, R, ox, oy)
  #gene_monitor_data <- get_gene_monitor_data(TransferEntropy_df, gene)
  #cell <- list(gene_network_picture, gene_monitor_data)
  #return (cell)
#}
#--------------------------------------------------------------------------------
#实现单个基因对其他基因的影响随时间的演化规律(网络)
get_times_gene_network_development <- function(TransferEntropy_df_list, gene){
  library(ggplot2)
  windowsFonts(CA = windowsFont("Calibri"))
  theme1 <- theme_bw() + theme(legend.position = 'top', 
                               panel.border = element_blank(),
                               panel.grid.major = element_line(linetype = 'dashed'), 
                               panel.grid.minor = element_blank(), 
                               legend.text = element_text(size=9,color='#003087',family = "CA"),
                               plot.title = element_text(size=15,color="#003087",family = "CA"), 
                               legend.key = element_blank(), 
                               axis.text = element_text(size=10,color='#003087',family = "CA"),
                               strip.text = element_text(size=12,color="#EF0808",family = "CA"),
                               strip.background = element_blank())
  cell_development <- list()
  times_gene_network_development <- list()
  cell_development <- lapply(TransferEntropy_df_list, get_gene_network_picture, gene)
  for(i in 1:length(cell_development)){
    times_gene_network_development[[names(cell_development)[i]]] <-  cell_development[[i]] + 
                                                          theme1 +
                                                          ggtitle(names(cell_development)[i]) + 
                                                          scale_x_continuous(name = gene,      #y轴坐标名称
                                                                             breaks = NULL) +  #连续的标签和坐标轴 
                                                          scale_y_continuous(name = '',        #y轴坐标名称
                                                                             breaks = NULL)    #连续的标签和坐标轴 
  }
  return (times_gene_network_development)
}
#--------------------------------------------------------------------------------
#实现每个基因对除去自己以外基因的影响随时间的演化规律(网络)
get_times_genes_network_development <- function(TransferEntropy_df_list){
  times_genes_network_development <- list()
  genes <- unique(as.vector(TransferEntropy_df_list[[1]]$x))
  for (i in 1:length(genes)){
    times_genes_network_development[[genes[i]]] <- get_times_gene_network_development(TransferEntropy_df_list, genes[i])
  }
  return(times_genes_network_development)
}
#--------------------------------------------------------------------------------
#监控量变到质变图谱模块
monitor_quant_function <- function(times_gene_monitor_data, gene){
  library('ggplot2')
  library('dplyr')
  library('reshape2')
  library('purrr')
  library('ggthemr')
  #ggthemr("solarized")
  mav <- function(x, n = 5){stats::filter(x, rep(1/n, n), sides = 1)}
  windowsFonts(CA = windowsFont("Calibri"))
  theme1 <- theme_bw() + theme(legend.position = 'left', 
                               panel.border = element_blank(),
                               panel.grid.major = element_line(linetype = 'dashed'), 
                               panel.grid.minor = element_blank(), 
                               legend.text = element_text(size=9,color='#003087',family = "CA"),
                               plot.title = element_text(size=15,color="#003087",family = "CA"), 
                               legend.key = element_blank(), 
                               axis.text = element_text(size=10,color='#003087',family = "CA"),
                               axis.text.x = element_text(angle=90),
                               strip.text = element_text(size=12,color="#EF0808",family = "CA"),
                               strip.background = element_blank())
  times_gene_monitor_data_up <- subset(times_gene_monitor_data, select = c(time, y , x.y, x.y.p_value))
  times_gene_monitor_data_up$p_value <- ifelse(times_gene_monitor_data_up[["x.y.p_value"]] <= 0.05, "x->y(p < 0.05)", "x<->y(p > 0.05)")
  p.dat.up <- data.frame(time = times_gene_monitor_data_up[["time"]], y = times_gene_monitor_data_up[["y"]], x.y = times_gene_monitor_data_up[["x.y"]], p_value = times_gene_monitor_data_up[["p_value"]])
  p <- ggplot()
  p1 <- p + geom_bar(aes(x = time, y = mav(x.y), group = y, fill = p_value), p.dat.up, stat='identity', position = "stack", width = 1.0)
  p2 <- p1 + geom_area(aes(x = time, y = mav(x.y), group = y, col = y), p.dat.up, position = "stack", alpha = 0.0)
  #p1 <- p + geom_bar(aes(x = time, y = x.y, group = y, fill = p_value), p.dat.up, stat='identity', position = "stack", width = 1.0)
  #p2 <- p1 + geom_area(aes(x = time, y = x.y, group = y, col = y), p.dat.up, position = "stack", alpha = 0.0)
  p_up <- p2 
  
  times_gene_monitor_data_down <- subset(times_gene_monitor_data, select = c(time, y , y.x, y.x.p_value))
  times_gene_monitor_data_down$p_value <- ifelse(times_gene_monitor_data_down[["y.x.p_value"]] <= 0.05, "y->x(p < 0.05)", "x<->y(p > 0.05)")
  p.dat.down <- data.frame(time = times_gene_monitor_data_down[["time"]], y = times_gene_monitor_data_down[["y"]], y.x = -times_gene_monitor_data_down[["y.x"]], p_value = times_gene_monitor_data_down[["p_value"]])
  p1_1 <- geom_bar(aes(x = time, y = mav(y.x), group = y, fill = p_value), p.dat.down, stat = 'identity', position = "stack", width = 1.0)
  p2_1 <- geom_area(aes(x = time, y = mav(y.x), group = y, col = y), p.dat.down, stat = 'identity', position = "stack", alpha = 0.0)
  #p1_1 <- geom_bar(aes(x = time, y = y.x, group = y, fill = p_value), p.dat.down, stat = 'identity', position = "stack", width = 1.0)
  #p2_1 <- geom_area(aes(x = time, y = y.x, group = y, col = y), p.dat.down, stat = 'identity', position = "stack", alpha = 0.0)

  #cat(length(unique(p.dat.up$p_value)))
  
  if(length(unique(p.dat.up$p_value)) == 2){
    times_gene_monitor_quant_picture <- p_up + p1_1 + p2_1 +
      theme1 +
      ggtitle(gene) + 
      labs(fill = "p_value", col = "pollutants") +
      scale_x_discrete(name = '', #x轴坐标名称
                       breaks = levels(p.dat.down$time)[seq(1,length(levels(p.dat.down$time)),length(levels(p.dat.down$time))/12)]) +  #连续的标签和坐标轴 
      scale_y_continuous(name = 'TransferEntropy') + #y轴坐标名称
      #breaks = NULL) +  #连续的标签和坐标轴 
      scale_fill_manual(values = c("#FF0000", "#8DD3C7", "#0033FF"))
    return (times_gene_monitor_quant_picture)
  }
  if(length(unique(p.dat.up$p_value)) == 1){
    times_gene_monitor_quant_picture <- p_up + p1_1 + p2_1 +
      theme1 +
      ggtitle(gene) + 
      labs(fill = "p_value", col = "pollutants") +
      scale_x_discrete(name = '', #x轴坐标名称
                       breaks = levels(p.dat.down$time)[seq(1,length(levels(p.dat.down$time)),length(levels(p.dat.down$time))/12)]) +  #连续的标签和坐标轴 
      scale_y_continuous(name = 'TransferEntropy') + #y轴坐标名称
      #breaks = NULL) +  #连续的标签和坐标轴 
      scale_fill_manual(values = c("#8DD3C7", "#0033FF"))
    return (times_gene_monitor_quant_picture)
  }
}
#--------------------------------------------------------------------------------
#监控质变图谱模块
monitor_qualit_function <- function(times_gene_monitor_data, gene, R, ox, oy){
  library('ggplot2')
  windowsFonts(CA = windowsFont("Calibri"))
  theme1 <- theme_bw() + theme(legend.position = 'left', 
                               panel.border = element_blank(),
                               panel.grid.major = element_line(linetype = 'dashed'), 
                               panel.grid.minor = element_blank(), 
                               legend.text = element_text(size=9,color='#003087',family = "CA"),
                               plot.title = element_text(size=15,color="#003087",family = "CA"), 
                               legend.key = element_blank(), 
                               axis.text = element_text(size=10,color='#003087',family = "CA"),
                               axis.text.x = element_text(angle=90),
                               strip.text = element_text(size=12,color="#EF0808",family = "CA"),
                               strip.background = element_blank())
  times_gene_monitor_data$val <- times_gene_monitor_data$x.y + times_gene_monitor_data$y.x
  times_gene_monitor_data$col <- times_gene_monitor_data$x.y - times_gene_monitor_data$y.x
  data <- aggregate(cbind(val, col) ~ y, times_gene_monitor_data, sum)
  #points <- getpoints(R, ox, oy, unique(as.vector(times_gene_monitor_data$y)))
  points <- getpoints(R, ox, oy, unique(c(as.vector(times_gene_monitor_data$y), as.vector(times_gene_monitor_data$x))))
  a <- merge(data, points,  by.x = "y", by.y = "point", all.x = T)
  p <- ggplot() 
  p_dian <- p + geom_point(aes(x = point_x, y = point_y, size = val, fill = col), a,  pch = 21, show.legend = F) + 
    geom_text(aes(x = point_x, y = point_y + 15, label = y), subset(a,point_y > 250),  size = 4.5, family = "CA") + 
    geom_text(aes(x = point_x, y = point_y - 15, label = y), subset(a,point_y <= 250), size = 4.5, family = "CA")+ 
    geom_text(aes(x = 250, y = 250, label = gene), size = 10, alpha = 0.3, family = "CA") +
    theme1 +
    labs(x='', y='') + 
    scale_x_continuous(name = '', #x轴坐标名称
                       breaks = NULL) +  #连续的标签和坐标轴 
    scale_y_continuous(name = '', #y轴坐标名称
                       breaks = NULL) +  #连续的标签和坐标轴 )  
    scale_fill_gradient2(low = "#0033FF", mid = "#8DD3C7",
                         high = "#FF0000", midpoint = 0, space = "Lab",
                         na.value = "lightgrey")
  return(list(p_dian, a))
}
#--------------------------------------------------------------------------------
#监控量变到质变图谱模块
monitor_quant_function2 <- function(times_gene_monitor_data, gene){
  library('ggplot2')
  library('dplyr')
  library('reshape2')
  windowsFonts(CA = windowsFont("Calibri"))
  theme1 <- theme_bw() + theme(legend.position = 'left', 
                               panel.border = element_blank(),
                               panel.grid.major = element_line(linetype = 'dashed'), 
                               panel.grid.minor = element_blank(), 
                               legend.text = element_text(size=9,color='#003087',family = "CA"),
                               plot.title = element_text(size=15,color="#003087",family = "CA"), 
                               legend.key = element_blank(), 
                               axis.text = element_text(size=10,color='#003087',family = "CA"),
                               axis.text.x = element_text(angle=90),
                               strip.text = element_text(size=12,color="#EF0808",family = "CA"),
                               strip.background = element_blank())
  times_gene_monitor_data_up <- subset(times_gene_monitor_data, select = c(time, y , x.y))
  #p.dat.up <- times_gene_monitor_data_up %>%
  #dcast(y ~ time)       #数据格式调整，y为行，time为列
  #p.dat.up <- p.dat.up[order(p.dat.up$y, decreasing = F),]
  #p.dat.up$y <- factor(p.dat.up$y, levels(as.factor(as.vector(p.dat.up$y))))
  #link_dat_up <- p.dat.up %>% 
  #group_by(Site) %>% #按大组Site进行分组，后续作图以Site进行分页展示
  #mutate_if(is.numeric, cumsum) %>% #判断是否为数值型变量，才执行本操作
  #as.data.frame()
  #bar.width <- 0.5
  #link_dat_up <- link_dat_up[, c(1, 2, rep(3:(ncol(link_dat_up)-1),each=2),ncol(link_dat_up))] #对数值型变量部分，掐头去尾，对中间的变量各复制一遍
  #link_dat_up <- data.frame(y=t(matrix(t(link_dat_up[,-1]), nrow=2)))#除去非数值型变量后，转换为两列的数据框，分别作为连接线的两个端点的纵坐标
  #link_dat_up$x.1 <- 1:(length(TransferEntropy_df_list) - 1) + bar.width/2 #设定连接线的两个端点的横坐标，1:n，n = 小组数-1，类似于砍木头，4刀5段
  #link_dat_up$x.2 <- 1:(length(TransferEntropy_df_list) - 1) + (1-bar.width/2)#同上
  #p.dat.up.rev <- melt(p.dat.up)
  #names(p.dat.up.rev) <- c("y", "time", "x.y")
  #p.dat.up.rev$y <- factor(p.dat.up.rev$y, levels = rev(as.vector(p.dat.up$y)))
  p.dat.up.rev <- times_gene_monitor_data_up
  a  <- aggregate(x.y~time, p.dat.up.rev, mean) 
  names(a) <- c("time", "mean(x.y|time)")
  b  <- aggregate(x.y~time, p.dat.up.rev, sd)
  names(b) <- c("time", "sd(x.y|time)")
  c <- merge(p.dat.up.rev, a , by = "time", all.x = T)
  d <- merge(c, b, by = "time", all.x = T)
  d[["p_value"]] = ifelse(d[["x.y"]] >= d[["mean(x.y|time)"]] + 1.56*d[["sd(x.y|time)"]], "x->y(p < 0.05)", "x<->y(p > 0.05)")
  p.dat.up <- data.frame(time = d[["time"]], y = d[["y"]], x.y = d[["x.y"]], p_value = d[["p_value"]])
  p <- ggplot()
  p1 <- p + geom_histogram(aes(x = time, y = x.y, group = y, fill = p_value), p.dat.up, stat='identity', position = "stack")
  #p2 <- p1 + geom_segment(data=link_dat_up, aes(x=x.1, xend=x.2, y=y.1, yend=y.2), col='black', inherit.aes = F)
  p_up <- p1 #+
  #theme1 +
  #ggtitle(gene) + 
  #scale_fill_manual(values = c(3, 4))
  times_gene_monitor_data_down <- subset(times_gene_monitor_data, select = c(time, y , y.x))
  #times_gene_monitor_data_down$y.x <- -(times_gene_monitor_data_down$y.x)
  #p.dat.down <- times_gene_monitor_data_down %>%
  #dcast(y ~ time)       #数据格式调整，y为行，time为列
  #p.dat.down <- p.dat.down[order(p.dat.down$y, decreasing = F),]
  #p.dat.down$y <- factor(p.dat.down$y, levels(as.factor(as.vector(p.dat.down$y))))
  #link_dat_down <- p.dat.down %>% 
  #group_by(Site) %>% #按大组Site进行分组，后续作图以Site进行分页展示
  #mutate_if(is.numeric, cumsum) %>% #判断是否为数值型变量，才执行本操作
  #as.data.frame()
  #link_dat_down <- link_dat_down[, c(1, 2, rep(3:(ncol(link_dat_down)-1),each=2),ncol(link_dat_down))] #对数值型变量部分，掐头去尾，对中间的变量各复制一遍
  #link_dat_down <- data.frame(y=t(matrix(t(link_dat_down[,-1]), nrow=2)))#除去非数值型变量后，转换为两列的数据框，分别作为连接线的两个端点的纵坐标
  #link_dat_down$x.1 <- 1:(length(TransferEntropy_df_list) - 1) + bar.width/2 #设定连接线的两个端点的横坐标，1:n，n = 小组数-1，类似于砍木头，4刀5段
  #link_dat_down$x.2 <- 1:(length(TransferEntropy_df_list) - 1) + (1-bar.width/2)#同上
  #p.dat.down.rev <- melt(p.dat.down)
  #names(p.dat.down.rev) <- c("y", "time", "y.x")
  #p.dat.down.rev$y <- factor(p.dat.down.rev$y, levels = rev(as.vector(p.dat.down$y)))
  #p.dat.down.rev$y.x <- -p.dat.down.rev$y.x
  p.dat.down.rev <- times_gene_monitor_data_down
  a1  <- aggregate(y.x ~ time, p.dat.down.rev, mean) 
  names(a1) <- c("time", "mean(y.x|time)")
  b1  <- aggregate(y.x ~ time, p.dat.down.rev, sd)
  names(b1) <- c("time", "sd(y.x|time)")
  c1 <- merge(p.dat.down.rev, a1 , by = "time", all.x = T)
  d1 <- merge(c1, b1, by = "time", all.x = T)
  d1[["p_value"]] = ifelse(d1[["y.x"]] >= d1[["mean(y.x|time)"]] + 1.96 * d1[["sd(y.x|time)"]], "y->x(p < 0.05)", "x<->y(p > 0.05)")
  p.dat.down <- data.frame(time = d1[["time"]], y = d1[["y"]], y.x = -d1[["y.x"]], p_value = d1[["p_value"]])
  p1_1 <- geom_histogram(aes(x = time, y = y.x, group = y, fill = p_value), p.dat.down, stat = 'identity', position = "stack")
  #p2_1 <- geom_segment(data = link_dat_down, aes(x = x.1, xend = x.2, y = y.1, yend = y.2), col = 'black', inherit.aes = F)
  times_gene_monitor_quant_picture <- p_up + p1_1 + #p2_1 +
    theme1 +
    ggtitle(gene) + 
    scale_x_discrete(name = '', #x轴坐标名称
                     breaks = levels(p.dat.down$time)[seq(1,length(levels(p.dat.down$time)),length(levels(p.dat.down$time))/12)]) +  #连续的标签和坐标轴 
    scale_y_continuous(name = 'TransferEntropy') + #y轴坐标名称
    #breaks = NULL) +  #连续的标签和坐标轴 
    scale_fill_manual(values = c("#FF0000", "#8DD3C7", "#0033FF"))
  return (times_gene_monitor_quant_picture)
}
#--------------------------------------------------------------------------------
#实现单个基因对其他基因的影响随时间的演化规律(监控)
get_times_gene_monitor_picture <- function(TransferEntropy_df_list, gene){
  R = 150
  ox = 250
  oy = 250
  cell_development <- NULL
  times_gene_monitor_data <- NULL
  times_gene_monitor_quant_picture <- NULL
  times_gene_monitor_qualit_picture <- NULL
  cell_development <- lapply(TransferEntropy_df_list, get_gene_monitor_data, gene)
  #for(i in 1:length(cell_development)){
    #times_gene_monitor_data <- rbind(times_gene_monitor_data, cell_development[[i]])
  #}
  times_gene_monitor_data <- do.call(rbind, cell_development)
  times_gene_monitor_data$time <- rep(names(cell_development), each = length(unique(c(as.vector(TransferEntropy_df_list[[1]]$y), as.vector(TransferEntropy_df_list[[1]]$x))))-1)
  times_gene_monitor_quant_picture <- monitor_quant_function(times_gene_monitor_data, gene)
  times_gene_monitor_qualit_picture <- monitor_qualit_function(times_gene_monitor_data, gene, R ,ox, oy)
  times_gene_monitor_data <- list(times_gene_monitor_quant_picture, times_gene_monitor_qualit_picture)
  return(times_gene_monitor_data)
}
#--------------------------------------------------------------------------------
#实现每个基因对除去自己以外基因的影响随时间的演化规律(监控)
get_times_genes_monitor_picture <- function(TransferEntropy_df_list){
  times_genes_monitor_picture <- NULL
  genes <- unique(as.vector(TransferEntropy_df_list[[1]]$x))
  for (i in 1:length(genes)){
    times_genes_monitor_picture[[genes[i]]] <- get_times_gene_monitor_picture(TransferEntropy_df_list, genes[i])
  }
  return(times_genes_monitor_picture)
}
#--------------------------------------------------------------------------------
#自定义函数批量整合目录中符合条件的数据框
readCSV <- function(dir_dta, city){
  #install.packages("foreach")
  #install.packages("doParallel")
  #install.packages("parallel")
  library('doParallel')
  library("foreach")
  library('parallel')
  file_list <- list.files(path = dir_dta, full.names = T)
  varSave_func <- function(x){
    df <- read.csv(file = x, sep = ",", header = T, fileEncoding = "utf-8")
    if(dim(df)[1] == 360){
      table_x <- subset(df, select = c("date", "hour", "type", "city"))
    }
  }
  # system.time({
    # cl.cores <- detectCores()
    # cl <- makeCluster(cl.cores)
    # registerDoParallel(cl)
    # a <- foreach(x = file_list) %dopar% varSave_func(x)
    # b <- do.call(rbind, a)
    # stopCluster(cl)
  # })
  system.time({
    cl.cores <- detectCores()
    cl <- makeCluster(cl.cores)
    a <- parLapply(cl, file_list, varSave_func)
    b <- do.call(rbind, a)
    stopCluster(cl)
  })
  return(b)
}
#--------------------------------------------------------------------------------
doing_things <- function(infile, start, end, len, step){
  #install.packages("foreach")
  library('cowplot')
  library('ggplot2')
  library('reshape2')
  library('doParallel')
  library('foreach')
  ######################加工目标转移熵##########################
  start <- 1
  end <- 480
  len <- 40
  step <- 1
  infile <- "C:\\Users\\Administrator\\Desktop\\results_60806\\fpkm_t_arr_select_27.csv"
  dir_result <- paste(strsplit(infile, "\\.")[[1]][1], "_result\\", sep = "")
  if(file_test("-d", dir_result) == F){dir.create(dir_result)}
  df <- read.table(file = infile, sep = ',', head = T, dec = '.', fill = T, na.strings = 'NA', stringsAsFactors = F)
  
  #mydata <- as.data.frame(scale(df[,2:ncol(df)]))
  #mydata1 <- cbind(data.frame(time = df$time), mydata)
  #df <- mydata1
  
  df$time <- as.character(df$time)
  times_df_list <- cutdf(df, start, end, len, step)
  system.time({
    TransferEntropy_df_list <- lapply(times_df_list, caculate_TransferEntropy)
  })
  save(TransferEntropy_df_list, file = paste(dir_result, "(", start, ",", end, ",", len, ",", step, ")", "27TransferEntropy_df_list.txt", sep = ""))
  #write.csv(as.data.frame(TransferEntropy_df_list), file = paste(dir_result, "(", start, ",", end, ",", len, ",", step, ")", "TransferEntropy_df_list.csv", sep = ""), row.names = TRUE)
  #load( paste(dir_result, "(", start, ",", end, ",", len, ",", step, ")", "scaletpmTransferEntropy_df_list.txt", sep = ""))
  #####################转移熵片段选择###########################
  begin <- 120
  over <- 240
  TransferEntropy_df_list_select <- TransferEntropy_df_list[begin:over]
  #####################可视化目标转移熵###########################
  times_genes_network_development <- get_times_genes_network_development(TransferEntropy_df_list_select)
  times_genes_monitor_picture <- get_times_genes_monitor_picture(TransferEntropy_df_list_select)
  ##########################目标显示##############################
  genes <- unique(as.vector(TransferEntropy_df_list_select[[1]]$x))
  for(i in 1:length(genes)){
    png(file = paste(dir_result, "(", start, ",", end, ",", len, ",", step, ")", "_", genes[i], '_network_picture.png', sep = ""), width = 4000, height = 4000, res = 100)
    print(plot_grid(plotlist = times_genes_network_development[[genes[i]]]))
    dev.off()
  }
  for(i in 1:length(genes)){
    png(file = paste(dir_result, "(", begin, ",", over, ",", len, ",", step, ")", "_", genes[i], '_monitor_quant_picture.png', sep = ""), width = 4000, height = 2000, res = 300)
    print(times_genes_monitor_picture[[genes[i]]][[1]])
    dev.off()
  }
  for(i in 1:length(genes)){
    png(file = paste(dir_result, "(", begin, ",", over, ",", len, ",", step, ")", "_", genes[i], '_monitor_qualit_picture.png', sep = ""), width = 4000, height = 3000, res = 300)
    print(times_genes_monitor_picture[[genes[i]]][[2]])
    dev.off()
  }
  
  to_gene <- subset(times_genes_monitor_picture[[genes[1]]][[2]][[2]], val < 0, select = c("y", "val"))
  write.csv(to_gene, file = paste(dir_result, "(", start, ",", end, ",", len, ",", step, ")", "to_gene_df_list.csv", sep = ""), row.names = F, quote = F)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  start <- 168
  end <- 23304
  len <- 168
  step <- 24
  city <- "拉萨"
  dir_data <- "E:\\BaiduNetdiskDownload\\全国空气质量\\data\\"
  dir_result <- "E:\\BaiduNetdiskDownload\\全国空气质量\\test\\"
  ######################提取目标#########################
  result <- readCSV(dir_data, city)
  df1 <- subset(result,  select = c("date","hour", "type", "北京"))
  df1 <- dcast(df1, date + hour ~ type)
  df1$time <- paste(df1$date, paste(df1$hour, ":00", sep = ""),sep = " ")
  df <- subset(df1, select = c("date", "CO", "NO2", "O3", "PM2.5", "SO2"))
  df <- aggregate(cbind(CO, NO2,O3,PM2.5,SO2) ~ date, data = df, mean)
  #write.table(df, file = paste(dir_result, "(", city , ")20170101-20191231.csv", sep = ""), sep = ",", quote = F, row.names = F)
  ######################加工目标转移熵##########################
  df <- read.csv(file = paste(dir_result, "(", city , ")20170101-20191231.csv", sep = ""), sep = ',', head = T, dec = '.', fill = T, na.strings = 'NA',  stringsAsFactors = F)
  names(df) <- c("time", "gene1", "gene2", "gene3", "gene4", "gene5")
  df$time <- as.character(df$time)
  times_df_list <- cutdf(df, start, end, len, step)
  #system.time({
    #cl.cores <- detectCores()
    #cl <- makeCluster(cl.cores)
    #registerDoParallel(cl)
    #TransferEntropy_df_list = foreach(i = times_df_list) %dopar% caculate_TransferEntropy(i)
    #stopCluster(cl)
  #})
  #system.time({
    #cl.cores <- detectCores()
    #l <- makeCluster(cl.cores)
    #TransferEntropy_df_list <- parLapply(cl, times_df_list, caculate_TransferEntropy)
    #stopCluster(cl)
  #})
  system.time({
    TransferEntropy_df_list <- lapply(times_df_list, caculate_TransferEntropy)
  })
  save(TransferEntropy_df_list, file = paste(dir_result, "(", start, ",", end, ",", len, ",", step, ")", "(" ,city, ")20170101-20191231TransferEntropy_df_list.txt", sep = ""))
  #load(paste(dir_result, "(", start, ",", end, ",", len, ",", step, ")",  "(", city, ")20170101-20191231TransferEntropy_df_list.txt", sep = ""))
  #first <- as.data.frame(TransferEntropy_df_list)
  #write.csv(first, file = paste(dir_result, '1.csv', sep = ""), row.names = TRUE)
  #####################可视化目标转移熵###########################
  times_genes_network_development <- get_times_genes_network_development(TransferEntropy_df_list)
  times_genes_monitor_picture <- get_times_genes_monitor_picture(TransferEntropy_df_list)
  #save(times_genes_monitor_picture, file = paste(dir_result, "(", start, ",", end, ",", len, ",", step, ")", "(" ,city, ")20170101-20191231times_genes_monitor_picture", sep = ""))
  ######################目标显示##########################
  #load(paste(dir_result, "(", start, ",", end, ",", len, ",", step, ")", "(" ,city, ")20170101-20191231times_genes_monitor_picture", sep = ""))
  genes <- unique(as.vector(TransferEntropy_df_list[[1]]$y))
  for(i in 1:length(genes)){
    png(file = paste(dir_result, "(", start, ",", end, ",", len, ",", step, ")",  "(", city, ")",  "_", genes[i], '_network_picture1.png', sep = ""), width = 4000, height = 4000)
    print(plot_grid(plotlist = times_genes_network_development[[genes[i]]]))
    dev.off()
  }
  for(i in 1:length(genes)){
    #setEPS()
    #post
    #postscript(file = paste(genes[i],'_monitor_picture.ps'), width = 1000, height = 500)
    png(file = paste(dir_result, "(", start, ",", end, ",", len, ",", step, ")",  "(", city, ")", "_", genes[i], '_monitor_quant_picture1.png', sep = ""), width = 1500, height = 500)
    print(times_genes_monitor_picture[[genes[i]]][[1]])
    dev.off()
    #p <- times_genes_monitor_picture[[genes[i]]]
    #ggsave(file = paste(genes[i],'_monitor_picture.tiff'), width = 1000, height = 500, dpi = 600, limitsize = FALSE)
  }
  for(i in 1:length(genes)){
    png(file = paste(dir_result, "(", start, ",", end, ",", len, ",", step, ")",  "(", city, ")",  "_", genes[i], '_monitor_qualit_picture1.png', sep = ""), width = 800, height = 500)
    print(times_genes_monitor_picture[[genes[i]]][[2]])
    dev.off()
  }
}
#--------------------------------------------------------------------------------
main <- function(){
  library(optparse)
  if (TRUE){
    option_list <- list(
      make_option(c("-i", "--input"), type = "character", default = NULL,
                 help = "Please enter the input file! [default %default]", metavar = "character"),
      make_option(c("-s", "--start"), type = "integer", default = NULL,
                  help = "Please enter the start site! [default %default]", metavar = "number"),
      make_option(c("-e", "--end"), type = "integer", default = NULL,
                  help = "Please enter the end site! [default %default]", metavar = "number"),
      make_option(c("-l", "--length"), type = "integer", default = NULL,
                  help = "Please enter the window length! [default %default]", metavar = "number"),
      make_option(c("-s", "--step"), type = "integer", default = NULL,
                  help = "Please enter the step length! [default %default]", metavar = "number")
    );
    opts <- parse_args(OptionParser(option_list=option_list))
  }
  if (is.null(opts$input)){
    cat("\nNo input file to read! Please enter the input file!\n\n")
	cat("Please use -h or --help for useage information!\n\n\n")
    q(status = 1)
  }
  if (is.null(opts$start)){
    cat("\nNo input start site! Please enter the start site!\n\n")
	cat("Please use -h or --help for useage information!\n\n\n")
    q(status = 1)
  }
  if (is.null(opts$end)){
    cat("\nNo input end site! Please enter the end site!\n\n")
	cat("Please use -h or --help for useage information!\n\n\n")
    q(status = 1)
  }
  if (is.null(opts$length)){
    cat("\nNo input the length of window! Please enter the window length!\n\n")
	cat("Please use -h or --help for useage information!\n\n\n")
    q(status = 1)
  }
  if (is.null(opts$step)){
    cat("\nNo input step length! Please enter the step length!\n\n")
	cat("Please use -h or --help for useage information!\n\n\n")
    q(status = 1)
  }
  doing_things(opts$infile, opts$start, opts$end, opts$length, opts$step)
}
#--------------------------------------------------------------------------------
main()
#--------------------------------------------------------------------------------
