#
#	Rpackage.R
#Thu Sep  8 13:47:26 CEST 2022

# T: Lars, F: Stefan
if (F) {
	packagesDir = "C:\\Users\\lljvanderburg\\Desktop";
	sourceFiles = '..\\KIRhaplotypes\\src\\Functions\\Working_functions';
} else {
	packagesDir = '~/src/Rpackages';
	sourceFiles = '../KIRhaplotypes/src/Functions/Working_functions';
}


packageDefinition = list(
	name = 'haploemcnv',
	files = list.files(sourceFiles, full.names = TRUE, recursive = TRUE, pattern = '.R$'),
	instFiles = list(data = 'data-raw/HaploemCNV_data.rda'),
	#testing = list(
	#	doInstall = TRUE,
	#	tests = c('RtestsPackages/package/package.R')
	#),
	description = list(
		title = 'EM algorithms for reconstructing haplotypes in coplex genetic regions',
		# version to be documented in news section
		#version = '0.1-0',
		author = 'Lars van der Burg <lars.lj.vdburg@gmail.com>',
		description = 'EM algorithms for reconstructing haplotypes in coplex genetic regions inclouding copy number variations.',
		depends = c(
					'stringr', 'Matrix', 'glmnet', 'MASS', "rlang", "gridExtra", "philentropy", 'gtools', "penalized", "cluster", "survival", "dplyr"
		            ),
		enhances = c(
		            'ggplot2',  "ggnewscale", "RColorBrewer", "readxl", "mctest", "ggpubr", "cowplot", "scales", "ggtern", "bench", "ggthemes"
		            # not available on CRAN
		            #"gt", 
		),
		suggests = c(),
		news = "0.3-0	Documentation\n0.2-0	Installable version\n0.1-0	Initial release",
		license = 'LGPL-2'
		#vignettes = "vignettes/vignette-package.Rmd"
	),
	git = list(
		readme = '## Installation\n```{r}\nlibrary(devtools);\ninstall_github("sboehringer/haploemcnv")\n```\n',
		push = FALSE,
		pushOnNewVersion = TRUE,
		remote = 'https://github.com/sboehringer/haploemcnv.git'
	)
)

