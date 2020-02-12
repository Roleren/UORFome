if (file.exists("/export/valenfs/projects/Hakon/slackr.config")) {
  slackr::slackr_setup(config_file = "/export/valenfs/projects/Hakon/slackr.config")
  library(httr)
  ggslackr <- function (plot = last_plot(), channels = Sys.getenv("SLACK_CHANNEL"), 
            scale = 1, width = par("din")[1], height = par("din")[2], 
            units = c("in", "cm", "mm"), dpi = 300, limitsize = TRUE, 
            api_token = Sys.getenv("SLACK_API_TOKEN"), file = "ggplot", 
            device = "pdf", ...) 
  {
    loc <- Sys.getlocale("LC_CTYPE")
    Sys.setlocale("LC_CTYPE", "C")
    on.exit(Sys.setlocale("LC_CTYPE", loc))
    ftmp <- tempfile(file, fileext = paste0(".", device))
    ggsave(filename = ftmp, plot = plot, scale = scale, width = width, 
           height = height, units = units, dpi = dpi, limitsize = limitsize, 
           ...)
    modchan <- slackr_chtrans(channels)
    res <- POST(url = "https://slack.com/api/files.upload", 
                add_headers(`Content-Type` = "multipart/form-data"), 
                body = list(file = upload_file(ftmp), token = api_token, 
                            channels = modchan))
    invisible(res)
  }
  ggslackR <- function(plot = last_plot(), width = 300, height = 150) ggslackr(plot = plot, dpi = 300, width = width, height = height, units = "mm", device = "pdf")
} else {
  ggslackR <- function(plot = last_plot(), width = 300, height = 150) {
    message("You dont have slackR set up, so no plots will be outputted!")
    return(NULL)
  }
}
