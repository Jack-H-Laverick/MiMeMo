## Build pkgdown site

# Run to build the website
pkgdown::build_site()


# Run to place the hex sticker in the right places
usethis::use_logo("logo.svg")
pkgdown::build_favicons(pkg = ".", overwrite = T)
