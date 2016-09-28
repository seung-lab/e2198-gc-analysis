

#Color palette from iwanthue
color_palette_iwh <- c(
  rgb(2,115,240, max=255),
  rgb(199,204,19, max=255),
  rgb(198,116,255, max=255),
  rgb(255,47,122, max=255),
  rgb(0,184,171, max=255),
  rgb(220,60,31, max=255),
  rgb(221,184,250, max=255),
  rgb(21,171,44, max=255),
  rgb(174,0,137, max=255),
  rgb(56,112,0, max=255),
  rgb(165,138,63, max=255),
  rgb(130,60,78, max=255),
  rgb(198,101,77, max=255),
  rgb(255,157,60, max=255),
  rgb(106,175,103, max=255),
  rgb(104,129,199, max=255),
  rgb(196,102,165, max=255)
)

color_palette_iwh2 <- c(
  #for 17028
  #rgb(255,255,255, max=255),
  #rgb(0,0,0),
  rgb(210,92,43, max=255),
  rgb(106,92,209, max=255),
  rgb(95,197,90, max=255),
  rgb(178,88,206, max=255),
  rgb(83,147,44, max=255),
  rgb(201,73,169, max=255),
  rgb(171,182,70, max=255),
  rgb(76,124,220, max=255),
  rgb(212,158,60, max=255),
  rgb(169,129,227, max=255),
  rgb(83,157,101, max=255),
  rgb(217,64,127, max=255),
  rgb(75,188,173, max=255),
  rgb(207,66,72, max=255),
  rgb(87,171,223, max=255),
  rgb(156,98,47, max=255),
  rgb(155,153,221, max=255),
  rgb(121,122,50, max=255),
  rgb(122,86,161, max=255),
  rgb(226,146,112, max=255),
  rgb(73,110,168, max=255),
  rgb(166,75,91, max=255),
  rgb(212,137,210, max=255),
  rgb(226,133,164, max=255),
  rgb(159,81,132, max=255),
  "gray80"
)

color_palette_iwh_somahist <- c(
  "#c57c3c",
  "#ab62c0",
  "#72a555",
  "#ca5670",
  "#638ccc",
  'gray40'
)

#not really used anymore, but could be useful later
colorpal <- function( counts ){
  
  cols <- c()
  cols <- c(cols, colorRampPalette(c("black","red")))
  cols <- c(cols, colorRampPalette(c("black","blue")))
  cols <- c(cols, colorRampPalette(c("black","green")))
  cols <- c(cols, colorRampPalette(c("black","yellow")))
  cols <- c(cols, colorRampPalette(c("black","purple")))
  cols <- c(cols, colorRampPalette(c("black","white")))
  
  res <- c()
  for( i in 1:length(counts) ){
    col_fn <- cols[[i]]
    res <- c(res, col_fn(counts[i]+2)[2:(counts[i]+1)])
  }
  res
}