#' Send plots to slack channel in one line
#'
#' You need a config file:
#' To make it, follow:
#' https://github.com/hrbrmstr/slackr
#' @param plot last_plot(), a ggplot sent to Rstudio frame viewer
#' @param width in mm, default (300)
#' @param heigth in mm, default (150)
#' @param config_file path to you config file
#' @import httr
#' @export
ggslackR <- function(plot = last_plot(), width = 300, height = 150,
                     config_file = "/export/valenfs/projects/Hakon/slackr.config") {
  if (file.exists(config_file)) {
    slackr::slackr_setup(config_file = config_file)
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
    return(ggslackr(plot = plot, dpi = 300, width = width, height = height, units = "mm", device = "pdf"))
  }
  message("You dont have slackR set up, so no plots will be outputted!")
  return(NULL)
}

