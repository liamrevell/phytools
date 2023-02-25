\name{anoletree}
\alias{anole.data}
\alias{anoletree}
\alias{bonyfish.data}
\alias{bonyfish.tree}
\alias{butterfly.data}
\alias{butterfly.tree}
\alias{cordylid.data}
\alias{cordylid.tree}
\alias{flatworm.data}
\alias{flatworm.tree}
\alias{mammal.data}
\alias{mammal.tree}
\alias{primate.data}
\alias{primate.tree}
\alias{salamanders}
\alias{sunfish.data}
\alias{sunfish.tree}
\alias{tropidurid.data}
\alias{tropidurid.tree}
\alias{vertebrate.data}
\alias{vertebrate.tree}
\alias{wasp.data}
\alias{wasp.trees}
\title{Phylogenetic datasets}
\description{
	Various phylogenetic datasets for comparative analysis.
}
\details{
	\code{anoletree} is a phylogeny of Greater Antillean anole species with a mapped discrete character - \emph{ecomorph class}. \code{anole.data} is a data frame of morphological characters. Data and tree are from Mahler et al. (2010).
	
	\code{butterfly.tree} and \code{butterfly.data} are a phylogeny and dataset of habitat use in Mycalesina butterflies from Halali et al. (2020).
	
	\code{bonyfish.tree} and \code{bonyfish.data} are a phylogeny and dataset of spawning mode and parental care in bony fishes from Benun Sutton and Wilson (2019).
	
	\code{cordylid.tree} and \code{cordylid.data} are a phylogeny and dataset of morphological traits for three different principal components axes from Broeckhoven et al. (2016).

	\code{flatworm.tree} and \code{flatworm.data} are a phylogeny and dataset of habitat preferences for flatworms from Benitez-Alvarez et al. (2020). \code{flatworm.tree} has been made ultrametric using penalized likelihood.
	
	\code{mammal.tree} and \code{mammal.data} are the phylogeny and dataset for mammal body size and home range size from Garland et al. (1992).
	
	\code{primate.tree} and \code{primate.data} are a phylogeny and phenotypic trait dataset from Kirk and Kay (2004).

	\code{salamanders} is a phylogeny of \emph{Plethodon} salamanders from Highton and Larson (1979). According to Wikipedia, the genus \emph{Plethodon} contains 55 species in total.
	
	\code{sunfish.tree} and \code{sunfish.data} are the phylogeny and dataset for Centrarchidae and buccal morphology (respectively) from Revell and Collar (2009).

	\code{tropidurid.tree} and \code{tropidurid.data} are the phylogeny and bivariate quantitative trait dataset of tropidurid lizards from Revell et al. (2022).
	
	\code{vertebrate.tree} is a time-calibrated phylogeny of vertebrates and \code{vertebrate.data} is a dataset of phenotypic traits. The phylogeny is from \url{http://www.timetree.org/} (Hedges et al. 2006).
	
	\code{wasp.trees} and \code{wasp.data} are the phylogeny and host-parasite associations from Lopez-Vaamonde et al. (2001).
}
\usage{
data(anole.data)
data(anoletree)
data(bonyfish.data)
data(bonyfish.tree)
data(butterfly.data)
data(butterfly.tree)
data(cordylid.data)
data(cordylid.tree)
data(flatworm.data)
data(flatworm.tree)
data(mammal.data)
data(mammal.tree)
data(primate.data)
data(primate.tree)
data(salamanders)
data(sunfish.data)
data(sunfish.tree)
data(tropidurid.data)
data(tropidurid.tree)
data(vertebrate.data)
data(vertebrate.tree)
data(wasp.data)
data(wasp.trees)
}
\format{
	\code{anoletree} is an object of class \code{"simmap"}. \code{anole.data} is a data frame.
	
	\code{bonyfish.tree} is an object of class \code{"phylo"}. \code{bonyfish.data} is a data frame.
	
	\code{butterfly.tree} is an object of class \code{"phylo"}. \code{butterfly.data} is a data frame.
	
	\code{cordylid.tree} is an object of class \code{"phylo"}. \code{cordylid.data} is a data frame.
	
	\code{flatworm.tree} is an object of class \code{"phylo"}. \code{flatworm.data} is a data frame.
	
	\code{mammal.tree} is an object of class \code{"phylo"}. \code{mammal.data} is a data frame.
	
	\code{primate.tree} is an object of class \code{"phylo"}. \code{primate.data} is a data frame.

	\code{salamanders} is an object of class \code{"phylo"}.
	
	\code{sunfish.tree} is an object of class \code{"simmap"}. \code{sunfish.data} is a data frame.
	
	\code{tropidurid.tree} is an object of class \code{"simmap"}. \code{tropidurid.data} is a data frame.
	
	\code{vertebrate.tree} is an object of class \code{"phylo"}. \code{vertebrate.data} is a data frame.
	
	\code{wasp.trees} is an object of class \code{"multiPhylo"}. \code{wasp.data} is a data frame.
}
\source{
	Benitez-Alvarez, L., A. Maria Leal-Zanchet, A. Oceguera-Figueroa, R. Lopes Ferreira, D. de Medeiros Bento, J. Braccini, R. Sluys, and M. Riutort (2020) Phylogeny and biogeography of the Cavernicola (Platyhelminthes: Tricladida): Relicts of an epigean group sheltering in caves? \emph{Molecular Phylogenetics and Evolution}, \bold{145}, 106709.
	
	Benun Sutton, F., and A. B. Wilson (2019) Where are all the moms? External fertilization predicts the rise of male parental care in bony fishes. \emph{Evolution}, \bold{73}, 2451-2460.

	Broeckhoven, C., G. Diedericks, C. Hui, B. G. Makhubo, P. le Fras N. Mouton (2016) Enemy at the gates: Rapid defensive trait diversification in an adaptive radiation of lizards. \emph{Evolution}, \bold{70}, 2647-2656.
	
	Garland, T., Jr., P. H. Harvey, and A. R. Ives (1992) Procedures for the analysis of comparative data using phylogenetically independent contrasts. \emph{Systematic Biology}, \bold{41}, 18-32.
	
	Kirk, E. C., and R. F. Kay (2004) The evolution of high visual acuity in the Anthropoidea. In: Ross, C. F., Kay R. F. (Eds), \emph{Anthropoid Origins. Developments in Primatology: Progress and Prospects}, 539-602. Springer, Boston, MA.
	
	Halali, S., E. van Bergen, C. J. Breuker, P. M. Brakefield, and O. Brattstrom (2020) Seasonal environments drive convergent evolution of a faster pace-of-life in tropical butterflies. \emph{Ecology Letters}, \bold{24}, 102-112.
	
	Hedges, S. B., J. Dudley, and S. Kumar (2006) TimeTree: A public knowledgebase of divergence times among organisms. \emph{Bioinformatics}, \bold{22}, 2971-2972.

	Highton, R., and A. Larson (1979) The genetic relationships of the salamanders of the genus \emph{Plethodon}. \emph{Systematic Zoology}, \bold{28}, 579-599.
	
	Lopez-Vaamonde, C., J. Y. Rasplus, G. D. Weiblen, and J. M. Cook (2001) Molecular phylogenies of fig wasps: Partial cocladogenesis of pollinators and parasites. \emph{Molecular Phylogenetics and Evolution}, \bold{21}, 55-71.

	Mahler, D. L, L. J. Revell, R. E. Glor, and J. B. Losos (2010) Ecological opportunity and the rate of morphological evolution in the diversification of Greater Antillean anoles. \emph{Evolution}, \bold{64}, 2731-2745.
	
	Revell, L. J., and D. C. Collar (2009) Phylogenetic analysis of the evolutionary correlation using likelihood. \emph{Evolution}, \bold{63}, 1090-1100.

	Revell, L. J., K. S. Toyama, and D. L. Mahler (2022) A simple hierarchical model for heterogeneity in the evolutionary correlation a phylogenetic tree. \emph{PeerJ}, \bold{10}, e13910.
}
\keyword{datasets}
