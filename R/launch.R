#' Launch the MultiplexIFGater
#'
#' Run the Shiny App MultiplexIFGater
#' @export
tSNEfier <- function () {
    shiny::runApp(system.file('shiny',package="tSNEfier"),launch.browser = T)
}
