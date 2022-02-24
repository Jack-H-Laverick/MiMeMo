# MiMeMo.tools 0.4.0 
<span style="color:grey;">24/02/2022</span>

* `characterise_flows` gains a precision argument to decrease the gap between sample points. A new check also tests whether the sampling points are in different polygons, when this fails the offending points are returned.
*`characterise_flows` bug fix swapping `st_line_sample` for `st_centroid`.
*`box_zoom` lets you control the ggplot plotting window with the bounding box of an sf object.
*`NM_volume_summary` gains a switch to ignore ice variables when ice_mod files weren't processed.
* Updated pkgdown site to V2 and styling to match group landing page.

# MiMeMo.tools 0.3.1 
<span style="color:grey;">26/12/2020</span>

* Initialised a test suite.

# MiMeMo.tools 0.2.1 
<span style="color:grey;">03/12/2020</span>

* Added functions to assist with unit conversions.

# MiMeMo.tools 0.2.0 
<span style="color:grey;">09/10/2020</span>

* Migrated NEMO-MEDUSA functions into their own dedicated package *nemomedusR*.
* Specified package dependencies so that the tidyverse and other useful packages are automatically loaded when calling MiMeMo.tools with `library()`

# MiMeMo.tools 0.1.0 
<span style="color:grey;">15/04/2020</span>

* Finished checking that existing scripts execute correctly using MiMiMo.tools. Lots of bug hunting!
* Added a scripts overview page to the website.
* Added a data sources page to the website.

# MiMeMo.tools 0.0.1 
<span style="color:grey;">16/03/2020</span>

* Finished documenting existing functions.
* Started a package website.
* Very importantly, made a hex sticker.

# MiMeMo.tools 0.0.0.1 
<span style="color:grey;">09/03/2020</span>

* Began work on MiMeMo.tools by transferring function files into a package.
