.onLoad <- function(libname, pkgname)
{
    suppressMessages({
        addResourcePath("js", system.file("www", "js", package="interactiveDisplayBase"))
        addResourcePath("css", system.file("www", "css", package="interactiveDisplayBase"))
    })
}