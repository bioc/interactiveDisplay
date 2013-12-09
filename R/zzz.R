.onLoad <- function(libname, pkgname)
{
    suppressMessages({
        addResourcePath("js", system.file("www", "js", package="interactiveDisplay"))
        addResourcePath("css", system.file("www", "css", package="interactiveDisplay"))
    })
}