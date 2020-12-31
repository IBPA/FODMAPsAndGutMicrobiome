#' Theme related utilities relating to plot colors, fonts, etc.

# Create the theme for ggplot figures
get_base_theme <- function(font_ratio){
  return(
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.box.background = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(fill=NA, size=1),
          legend.key.size = unit(0.5, "cm"),
          text = element_text(family="Helvetica", size=7*font_ratio),
          axis.text = element_text(size=8*font_ratio),
          axis.title = element_text(size=8*font_ratio, face='bold'),
          strip.text = element_text(size=8*font_ratio, vjust = 0)
    )
  )
}

my_base_theme <- get_base_theme(1.25)
colors_assigned <- c("A"="#619CFF", 
                     "B"="#F8766D", "C"="#E1BC29", "D" = "#00BA38",
                     "E"="#7768AE", "F"="#000000", "Transparent"="#ffffff00")