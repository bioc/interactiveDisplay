.onLoad <- function(libname, pkgname)
{
    suppressMessages({
        addResourcePath("js-interactiveDisplayBase", system.file("www", "js", package="interactiveDisplayBase"))
        addResourcePath("css-interactiveDisplayBase", system.file("www", "css", package="interactiveDisplayBase"))
        addResourcePath("js-interactiveDisplay", system.file("www", "js", package="interactiveDisplay"))
        addResourcePath("css-interactiveDisplay", system.file("www", "css", package="interactiveDisplay"))
    })
}

.onAttach <- function(libname, pkgname) {
    msg <- sprintf(
        "Package '%s' is deprecated and will be removed from Bioconductor
         version %s", pkgname, "3.23")
    .Deprecated(msg=paste(strwrap(msg, exdent=2), collapse="\n"))
}
