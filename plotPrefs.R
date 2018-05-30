
theme_normal <- function(base_size = 12, base_family = "Helvetica")
{
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(text = element_text(size =12,family = 'Helvetica',colour="black",face='bold',hjust = 0.5, vjust = 0.5, angle = 0, lineheight =0,
                               margin = margin(), debug = FALSE),
          axis.ticks.length=unit(-0.15, "cm"), 
          axis.text.x = element_text(size =12,face = "bold", color = "black",margin=margin(0.3,0,0,0.3,"cm")),
          axis.title.x = element_text(size =12,face = "bold", color = "black",margin=margin(0.2,0,0,0.1,"cm")),
          axis.text.y = element_text(size =12,face = "bold", color = "black",margin=margin(0,0.3,0,0,"cm")),
          axis.title.y = element_text(size =12,face = "bold", color = "black",margin=margin(0,0.1,0,0,"cm"),angle=90),
          axis.title.y.right = element_blank(),
          axis.text.y.right = element_blank(),
          axis.title.x.top = element_blank(),
          axis.text.x.top = element_blank(),
          axis.line =  element_blank(),
          panel.border=element_rect(fill=NA,size=0.5),
          legend.title = element_blank(),legend.key=element_blank(),
          legend.text = element_text(size=12),
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          panel.grid = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
}

