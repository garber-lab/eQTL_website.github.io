# eQTL_website.github.io

## Introduction

This repository is responsible for maintaining eQTL data, rendering it as an HTML file, and hosting it as a standalone static website via GitHub Pages.

The following are the primary files that are integral to the rendering and deployment of the website:

-   `index.Rmd` (R source code with the data to be rendered in HTML)
-   `index.html` (Rendered HTML of the website that will be deployed)
-   `_site.yml` (Site configuration file)

------------------------------------------------------------------------

## Description

To be included later....

### `_site.yml`

This file is integral for deploying the website because it can have several options that affect website output, including where it is written and what files are included and excluded from the website. View [this page](https://bookdown.org/yihui/rmarkdown/rmarkdown-site.html#site-configuration) for more info on the `_site.yml` file.

For this project specifically, I specified that the project root directory should be the output file for any site content that is generated when the HTML is rendered. This was done by specifying the field `output_dir: "."` in the file. Otherwise by default, all generated output would be stored in a `_site/` subdirectory, which wouldn't facilitate webpage deployment via GitHub Pages.


---

## Getting Started

The URL for the publicly hosted website is: https://garber-lab.github.io/eQTL_website.github.io/

To render the HTML website yourself, you can do either of the following: 

(1) Clone this repository and knit the `index.Rmd` file
(2) Clone this repository, navigate to the root directory in RStudio, and enter the `rmarkdown::render_site()` function


---

## Authors

The data and plots for this project comes from Shuo Shan (Crystal).

The project was rendered into HTML and deployed as a website by Thomas Jacob.

## Acknowledgements

Special thanks to our PI, Dr. Manuel Garber, for their assistance in the completion of this project. 

And special thanks to the UMass Chan Department of Genomics and Computational Biology for their resources.


