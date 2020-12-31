# Generic utility functions

lib.util.log <- function(message, obj = NULL){
  print(sprintf("** LOG ** %s", message))
  if (!is.null(obj)){
    print(obj)
  }
}

lib.util.ggsave <- function(filename, ...){
  lib.util.log(sprintf("Saved to into %s", filename))
  ggsave(filename, ...)
}

lib.util.load_fonts <- function(){
  # Next few steps are necessary for saving the pdf due to a font error
  library(showtext)
  showtext_auto()
  font_add(family="HiraKakuProN-W3", 
           file.path(Sys.getenv("HOME"), "/Library/Fonts/copyfonts.com_hirakakupro-w3-opentype.otf"))
  font_add(family="Helvetica",
           "/System/Library/Fonts/Helvetica.ttc")
}
