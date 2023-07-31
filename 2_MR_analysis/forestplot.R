###########################绘制MR分析结果的森林图################
########################Date:2023年6月7日#############################
#####################Author: Yidie########################

####代码参考：https://blog.csdn.net/dege857/article/details/127859291####

# 加载所需的R包
library(forestploter)
library(eoffice)
library(ggplot2)
library(grid)


dt <- read.table(paste0("/home/yidie/work_dir/analysis0607_nopleio/results/plot/forest_data/LAS.csv"), 
                 header = T, sep = ",", stringsAsFactors = FALSE)

####整理数据
# 添加空白列以显示 CI
# 用空格调整列宽
dt$` ` <- paste(rep(" ", 20), collapse = " ")

# 创建要显示的列，"%.2f"意为保留小数点后两位精度,去掉“NA”值
dt$`SNPs` <- ifelse(is.na(dt$pval), "", dt$No..of.SNPs)
dt$`SNPs` <- ifelse(is.na(dt$pval), "",
                      paste0("   ", dt$No..of.SNPs))
dt$`OR (95% CI)` <- ifelse(is.na(dt$pval), "",
                           sprintf("%.2f (%.2f to %.2f)",
                                   dt$or, dt$or_lci95, dt$or_uci95))
dt$`P-value` <- ifelse(is.na(dt$pval), "", dt$Pvalue)


# 定义简单的主题
tm <- forest_theme(base_size = 10,  #文本的大小
                   # Confidence interval point shape, line type/color/width
                   ci_pch = 15,   #可信区间点的形状
                   ci_col = "black",    #CI的颜色
                   ci_fill = "#4682B4",     #ci颜色填充
                   ci_alpha = 0.8,        #ci透明度
                   ci_lty = 1,            #CI的线型
                   ci_lwd = 1.5,          #CI的线宽
                   ci_Theight = 0.1, # Set an T end at the end of CI  ci的高度，默认是NULL
                   # Reference line width/type/color   参考线默认的参数，中间的竖的虚线
                   refline_lwd = 1,       #中间的竖的虚线
                   refline_lty = "dashed",
                   refline_col = "grey20",
                   # Vertical line width/type/color  垂直线宽/类型/颜色   可以添加一条额外的垂直线，如果没有就不显示
                   vertline_lwd = NULL,              #可以添加一条额外的垂直线，如果没有就不显示
                   vertline_lty = NULL,
                   vertline_col = NULL,
                   # Change summary color for filling and borders   更改填充和边框的摘要颜色
                   summary_fill = NULL,       #汇总部分大菱形的颜色
                   summary_col = NULL,
                   # Footnote font size/face/color  脚注字体大小/字体/颜色
                   footnote_cex = 0.6,
                   footnote_fontface = "italic",
                   footnote_col = "grey16",
                   padding = unit(c(4, 4), "mm")
                    )


p <- forest(dt[,c(1:2,15,14,16:17)],
            est = dt$or,
            lower = dt$or_lci95, 
            upper = dt$or_uci95,
            sizes = 0.4,
            ci_column = 4,
            xlim = c(0.2, 1.6),
            ticks_at = c(0.5, 0.8, 1.0, 1.4),
            ref_line = 1,
            theme = tm)


p_1 <- edit_plot(p, row = c(1, 2, 9, 16, 24, 25, 32, 39, 47, 48, 55, 62, 70,
                            71, 78, 85),
                 gp=gpar(fontface = "bold"))
p_2 <- edit_plot(p_1, row = c(1:91), which = "background",
                 gp = gpar(fill = "white"))
p_3 <- add_border(p_2, part = "header")


# 绘图
output_file <- paste0("/home/yidie/work_dir/analysis0607_nopleio/results/plot/forest_LAS.pdf")
ggsave(
  filename = output_file,
  p_3,
  width = 8,             
  height = 5.5,            
  units = "in",          
  dpi = 300              
)
output_file <- paste0("/home/yidie/work_dir/analysis0607_nopleio/results/plot/forest_LAS.tiff")
ggsave(
  filename = output_file,
  p_3,
  width = 8,             
  height = 5.5,            
  units = "in",          
  dpi = 300              
)

