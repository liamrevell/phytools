\name{anoletree}
\alias{anole.data}
\alias{anoletree}
\alias{flatworm.tree}
\alias{flatworm.data}
\alias{mammal.tree}
\alias{mammal.data}
\alias{salamanders}
\alias{sunfish.tree}
\alias{sunfish.data}
\alias{vertebrate.tree}
\alias{vertebrate.data}
\alias{wasp.trees}
\alias{wasp.data}
\title{Phylogenetic datasets}
\description{
	\code{anoletree} is a phylogeny of Greater Antillean anole species with a mapped discrete character - \emph{ecomorph class}. \code{anole.data} is a data frame of morphological characters. Data and tree are from Mahler et al. (2010).
	
	\code{flatworm.tree} and \code{flatworm.data} are a phylogeny and dataset of habitat preferences for flatworms from Benitez-Alvarez et al. (2000). \code{flatworm.tree} has been made ultrametric using penalized likelihood.
	
	\code{mammal.tree} and \code{mammal.data} are the phylogeny and dataset for mammal body size and home range size from Garland et al. (1992).
	
	\code{salamanders} is a phylogeny of \emph{Plethodon} salamanders from Highton and Larson (1979). According to Wikipedia, the genus \emph{Plethodon} contains 55 species in total.
	
	\code{sunfish.tree} and \code{sunfish.data} are the phylogeny and dataset for Centrarchidae and buccal morphology (respectively) from Revell and Collar (2009).
	
	\code{vertebrate.tree} is a time-calibrated phylogeny of vertebrates and \code{vertebrate.data} is a dataset of phenotypic traits. The phylogeny is from \href{http://www.timetree.org/}{http://www.timetree.org/} (Hedges et al. 2006).
	
	\code{wasp.trees} and \code{wasp.data} are the phylogeny and host-parasite associations from Lopez-Vaamonde et al. (2001).
}
\usage{
data(anole.data)
data(anoletree)
data(flatworm.tree)
data(flatworm.data)
data(mammal.data)
data(mammal.tree)
data(salamanders)
data(sunfish.data)
data(sunfish.tree)
data(vertebrate.tree)
data(vertebrate.data)
data(wasp.data)
data(wasp.trees)
}
\format{
	\code{anoletree} is an object of class \code{"simmap"}. \code{anole.data} is a data frame.
	
	\code{flatworm.tree} is an object of class \code{"phylo"}. \code{flatworm.data} is a data frame.
	
	\code{mammal.tree} is an object of class \code{"phylo"}. \code{mammal.data} is a data frame.
	
	\code{salamanders} is an object of class \code{"phylo"}.
	
	\code{sunfish.tree} is an object of class \code{"simmap"}. \code{sunfish.data} is a data frame.
	
	\code{vertebrate.tree} is an object of class \code{"phylo"}. \code{vertebrate.data} is a data frame.
	
	\code{wasp.trees} is an object of class \code{"multiPhylo"}. \code{wasp.data} is a data frame.
}
\source{
	Benitez-Alvarez, L., A. Maria Leal-Zanchet, A. Oceguera-Figueroa, R. Lopes Ferreira, D. de Medeiros Bento, J. Braccini, R. Sluys, and M. Riutort. (2020) Phylogeny and biogeography of the Cavernicola (Platyhelminthes: Tricladida): Relicts of an epigean group sheltering in caves? \emph{Molecular Phylogenetics and Evolution}, \bold{145}, 106709.

	Garland, T., Jr., P. H. Harvey, and A. R. Ives. (1992) Procedures for the analysis of comparative data using phylogenetically independent contrasts. \emph{Systematic Biology}, \bold{41}, 18-32.
	
	Hedges, S. B., J. Dudley, and S. Kumar. (2006) TimeTree: A public knowledgebase of divergence times among organisms. \emph{Bioinformatics}, \bold{22}, 2971-2972.

	Highton, R., and A. Larson. (1979) The genetic relationships of the salamanders of the genus \emph{Plethodon}. \emph{Systematic Zoology}, \bold{28}, 579-599.
	
	Lopez-Vaamonde, C., J. Y. Rasplus, G. D. Weiblen, and J. M. Cook. (2001) Molecular phylogenies of fig wasps: Partial cocladogenesis of pollinators and parasites. \emph{Molecular Phylogenetics and Evolution}, \bold{21}, 55-71.

	Mahler, D. L, L. J. Revell, R. E. Glor, and J. B. Losos. (2010) Ecological opportunity and the rate of morphological evolution in the diversification of Greater Antillean anoles. \emph{Evolution}, \bold{64}, 2731-2745.
	
	Revell, L. J., and D. C. Collar. (2009) Phylogenetic analysis of the evolutionary correlation using likelihood. \emph{Evolution}, \bold{63}, 1090-1100.
}
\keyword{datasets}
