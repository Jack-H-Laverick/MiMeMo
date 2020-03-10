
library(hexSticker)
library(desc)
library(ggplot2)

desc = desc::description$new()
fig_dir = file.path("man", "figures")
if (dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}

holder <- ggplot(x=1, y=1) + geom_point() + theme_void()

package = desc$get("Package")
hexSticker::sticker(
  subplot = holder, package = package, p_size = 30, s_x=1.22, s_y= 0.77, s_width = 0.900,
  s_height = 0.540, p_x = 0.85, h_fill = "#003399", h_color = "#FFCC33", dpi = 700,
  url = 'install_github("Jack-H-Laverick/MiMeMo.tools")', u_color = "white", u_size = 7,
  filename = file.path(fig_dir, "sticker.png"))

usethis::use_logo(file.path(fig_dir, "sticker.png"))
