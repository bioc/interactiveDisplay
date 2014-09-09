.onLoad <- function(libname, pkgname)
{
    suppressMessages({
        addResourcePath("js-interactiveDisplayBase", system.file("www", "js", package="interactiveDisplayBase"))
        addResourcePath("css-interactiveDisplayBase", system.file("www", "css", package="interactiveDisplayBase"))
        addResourcePath("js-interactiveDisplay", system.file("www", "js", package="interactiveDisplay"))
        addResourcePath("css-interactiveDisplay", system.file("www", "css", package="interactiveDisplay"))
    })
}