\name{anoletree}
\alias{anole.data}
\alias{anoletree}
\alias{ant.geog}
\alias{ant.tree}
\alias{bat.tree}
\alias{bat_virus.data}
\alias{betaCoV.tree}
\alias{bonyfish.data}
\alias{bonyfish.tree}
\alias{butterfly.data}
\alias{butterfly.tree}
\alias{cordylid.data}
\alias{cordylid.tree}
\alias{darter.tree}
\alias{eel.data}
\alias{eel.tree}
\alias{elapidae.tree}
\alias{flatworm.data}
\alias{flatworm.tree}
\alias{liolaemid.data}
\alias{liolaemid.tree}
\alias{mammal.data}
\alias{mammal.geog}
\alias{mammal.tree}
\alias{primate.data}
\alias{primate.tree}
\alias{salamanders}
\alias{sunfish.data}
\alias{sunfish.tree}
\alias{tortoise.geog}
\alias{tortoise.tree}
\alias{tropidurid.data}
\alias{tropidurid.tree}
\alias{vertebrate.data}
\alias{vertebrate.tree}
\alias{wasp.data}
\alias{wasp.trees}
\alias{whale.tree}
\title{Phylogenetic datasets}
\description{
	Various phylogenetic datasets for comparative analysis.
}
\details{
	\code{anoletree} is a phylogeny of Greater Antillean anole species with a mapped discrete character - \emph{ecomorph class}. \code{anole.data} is a data frame of morphological characters. Data and tree are from Mahler et al. (2010).
	
	\code{ant.tree} is a phylogeny containing 84 species of fungus farming attine ants. \code{ant.geog} is a vector containing biogegraphic information for the species of the tree. Both the phylogeny and biogeographic information were originally published by Branstetter et al. (2017).
	
	\code{bat.tree} is a ML phylogeny of bats from the families Mormoopidae, Molossidae, and Phyllostomidae, based on GenBank accession numbers reported in Caraballo (2022) and estimated using \pkg{phangorn}. \code{betaCoV.tree} is a ML phylogenetic tree of betacoronaviruses based on GenBank accession numbers reported in Caraballo (2022) and estimated using \pkg{phangorn}. \code{bat_virus.data} is an association table of host and virus from Caraballo (2022).
	
	\code{butterfly.tree} and \code{butterfly.data} are a phylogeny and dataset of habitat use in Mycalesina butterflies from Halali et al. (2020). The phylogeny of that study was adapted from one originally published by Brattstrom et al. (2020).
	
	\code{bonyfish.tree} and \code{bonyfish.data} are a phylogeny and dataset of spawning mode and parental care in bony fishes from Benun Sutton and Wilson (2019). Their phylogenetic tree was adapted from a tree originally published by Betancur-R et al. (2017)
	
	\code{cordylid.tree} and \code{cordylid.data} are a phylogeny and dataset of morphological traits for three different principal components axes from Broeckhoven et al. (2016).
	
	\code{darter.tree} is a phylogeny of 201 species of darter fish (Etheostomatinae) from Near et al. (2011).

	\code{eel.tree} and \code{eel.data} are a phylogeny and dataset of feeding mode and maximum total length from 61 species of elopomorph eels from Collar et al. (2014).
	
	\code{elapidae.tree} is a reconstructed phylogeny containing 175 species of snakes of the family Elapidae from Lee et al. (2016).

	\code{flatworm.tree} and \code{flatworm.data} are a phylogeny and dataset of habitat preferences for flatworms from Benitez-Alvarez et al. (2020). \code{flatworm.tree} has been made ultrametric using penalized likelihood.
	
	\code{liolaemid.tree} and \code{liolaemid.data} are a phylogenetic tree and phenotypic trait dataset of lizards from the family Liolaemidae from Esquerre et al. (2019).
	
	\code{mammal.tree} and \code{mammal.data} are the phylogeny and dataset for mammal body size and home range size from Garland et al. (1992). \code{mammal.geog} is a matrix containing geographic coordinates for various of the species in \code{mammal.tree} pulled from citizen science observations registered on the \emph{iNaturalist} web platform.
	
	\code{primate.tree} and \code{primate.data} are a phylogeny and phenotypic trait dataset from Kirk and Kay (2004).

	\code{salamanders} is a phylogeny of \emph{Plethodon} salamanders from Highton and Larson (1979). According to Wikipedia, the genus \emph{Plethodon} contains 55 species in total.
	
	\code{sunfish.tree} and \code{sunfish.data} are the phylogeny and dataset for Centrarchidae and buccal morphology (respectively) from Revell and Collar (2009).
	
	\code{tortoise.tree} and \code{tortoise.geog} are a phylogeny and dataset of geographic coordinates for Galapagos tortoises. The geographic locality information was approximated from Figure 1 of Poulakakis et al. (2020), and the tree was estimated from nucleotide sequence data published with the same study.

	\code{tropidurid.tree} and \code{tropidurid.data} are the phylogeny and bivariate quantitative trait dataset of tropidurid lizards from Revell et al. (2022).
	
	\code{vertebrate.tree} is a time-calibrated phylogeny of vertebrates and \code{vertebrate.data} is a dataset of phenotypic traits. The phylogeny is from \url{http://www.timetree.org/} (Hedges et al. 2006).
	
	\code{wasp.trees} and \code{wasp.data} are the phylogeny and host-parasite associations from Lopez-Vaamonde et al. (2001).
	
	\code{whale.tree} is a phylogeny of cetaceans originally published in Steeman et al. (2009).
}
\usage{
data(anole.data)
data(anoletree)
data(ant.geog)
data(ant.tree)
data(bat.tree)
data(bat_virus.data)
data(betaCoV.tree)
data(bonyfish.data)
data(bonyfish.tree)
data(butterfly.data)
data(butterfly.tree)
data(cordylid.data)
data(cordylid.tree)
data(darter.tree)
data(eel.data)
data(eel.tree)
data(elapidae.tree)
data(flatworm.data)
data(flatworm.tree)
data(liolaemid.tree)
data(mammal.data)
data(mammal.geog)
data(mammal.tree)
data(primate.data)
data(primate.tree)
data(salamanders)
data(sunfish.data)
data(sunfish.tree)
data(tortoise.geog)
data(tortoise.tree)
data(tropidurid.data)
data(tropidurid.tree)
data(vertebrate.data)
data(vertebrate.tree)
data(wasp.data)
data(wasp.trees)
data(whale.tree)
}
\format{
	\code{anoletree} is an object of class \code{"simmap"}. \code{anole.data} is a data frame.
	
	\code{ant.tree} is an object of class \code{"phylo"}. \code{ant.geog} is a vector.
	
	\code{bat.tree} and \code{betaCoV.tree} are objects of class \code{"phylo"}. \code{bat_virus.data} is a data frame.
	
	\code{bonyfish.tree} is an object of class \code{"phylo"}. \code{bonyfish.data} is a data frame.
	
	\code{butterfly.tree} is an object of class \code{"phylo"}. \code{butterfly.data} is a data frame.
	
	\code{cordylid.tree} is an object of class \code{"phylo"}. \code{cordylid.data} is a data frame.
	
	\code{darter.tree} is an object of class \code{"phylo"}.
	
	\code{eel.tree} is an object of class \code{"phylo"}. \code{eel.data} is a data frame.
	
	\code{elapidae.tree} is an object of class \code{"phylo"}.
	
	\code{flatworm.tree} is an object of class \code{"phylo"}. \code{flatworm.data} is a data frame.
	
	\code{liolaemid.tree} is an object of class \code{"phylo"}. \code{liolaemid.data} is a data frame.
	
	\code{mammal.tree} is an object of class \code{"phylo"}. \code{mammal.data} is a data frame. \code{mammal.geog} is a numerical matrix.
	
	\code{primate.tree} is an object of class \code{"phylo"}. \code{primate.data} is a data frame.

	\code{salamanders} is an object of class \code{"phylo"}.
	
	\code{sunfish.tree} is an object of class \code{"simmap"}. \code{sunfish.data} is a data frame.
	
	\code{tortoise.tree} is an object of class \code{"phylo"}. \code{tortoise.geog} is a data frame.
	
	\code{tropidurid.tree} is an object of class \code{"simmap"}. \code{tropidurid.data} is a data frame.
	
	\code{vertebrate.tree} is an object of class \code{"phylo"}. \code{vertebrate.data} is a data frame.
	
	\code{wasp.trees} is an object of class \code{"multiPhylo"}. \code{wasp.data} is a data frame.
	
	\code{whale.tree} is an object of class \code{"phylo"}.
}
\source{
	Benitez-Alvarez, L., A. Maria Leal-Zanchet, A. Oceguera-Figueroa, R. Lopes Ferreira, D. de Medeiros Bento, J. Braccini, R. Sluys, and M. Riutort (2020) Phylogeny and biogeography of the Cavernicola (Platyhelminthes: Tricladida): Relicts of an epigean group sheltering in caves? \emph{Molecular Phylogenetics and Evolution}, \bold{145}, 106709.
	
	Benun Sutton, F., and A. B. Wilson (2019) Where are all the moms? External fertilization predicts the rise of male parental care in bony fishes. \emph{Evolution}, \bold{73}, 2451-2460.
	
	Betancur-R, R., E. O. Wiley, G. Arratia, A. Acero, N. Bailly, M. Miya, G. Lecointre, and G. Orti (2017) Phylogenetic classification of bony fishes. \emph{BMC Ecology and Evolution}, \bold{17}, 162.
	
	Branstetter, M. G., A. Jesovnik, J. Sosa-Calvo, M. W. Lloyd, B. C. Faircloth, S. G. Brady, and T. R. Schultz (2017) Dry habitats were crucibles of domestication in the evolution of agriculture in ants. \emph{Proceedings of the Royal Society B}, \bold{284}, 20170095.
	
	Brattstrom, O., K. Aduse-Poku, E. van Bergen, V. French, and P. M. Brakefield (2020) A release from developmental bias accelerates morphological diversification in butterfly eyespots. \emph{Proceedings of the National Academy of Sciences}, \bold{177}, 27474-27480.

	Broeckhoven, C., G. Diedericks, C. Hui, B. G. Makhubo, P. le Fras N. Mouton (2016) Enemy at the gates: Rapid defensive trait diversification in an adaptive radiation of lizards. \emph{Evolution}, \bold{70}, 2647-2656.
	
	Caraballo, D. A. (2022) Cross-species transmission of bat coronaviruses in the Americas: Contrasting patterns between alphacoronavirus and betacoronavirus. \emph{Microbiology Spectrum}, \bold{10}, e01411-22.
	
	Collar, D. C., P. C. Wainwright, M. E. Alfaro, L. J. Revell, and R. S. Mehta (2014) Biting disrupts integration to spur skull evolution in eels. \emph{Nature Communications}, \bold{5}, 5505.
	
	Esquerre, D., D. Ramirez-Alvarez, C. J. Pavon-Vazquez, J. Troncoso-Palacios, C. F. Garin, J. S. Keogh, and A. D. Leache (2019) Speciation across mountains: Phylogenomics, species delimitation and taxonomy of the \emph{Liolaemus leopardinus} clade (Squamata, Liolaemidae). \emph{Molecular Phylogenetics and Evolution}, \bold{139}, 106524.
	
	Garland, T., Jr., P. H. Harvey, and A. R. Ives (1992) Procedures for the analysis of comparative data using phylogenetically independent contrasts. \emph{Systematic Biology}, \bold{41}, 18-32.
	
	Kirk, E. C., and R. F. Kay (2004) The evolution of high visual acuity in the Anthropoidea. In: Ross, C. F., Kay R. F. (Eds), \emph{Anthropoid Origins. Developments in Primatology: Progress and Prospects}, 539-602. Springer, Boston, MA.
	
	Halali, S., E. van Bergen, C. J. Breuker, P. M. Brakefield, and O. Brattstrom (2020) Seasonal environments drive convergent evolution of a faster pace-of-life in tropical butterflies. \emph{Ecology Letters}, \bold{24}, 102-112.
	
	Hedges, S. B., J. Dudley, and S. Kumar (2006) TimeTree: A public knowledgebase of divergence times among organisms. \emph{Bioinformatics}, \bold{22}, 2971-2972.

	Highton, R., and A. Larson (1979) The genetic relationships of the salamanders of the genus \emph{Plethodon}. \emph{Systematic Zoology}, \bold{28}, 579-599.
	
	Lee, M. S. Y., K. L. Saunders, B. King, and A. Palci (2016) Diversification rates and phenotypic evolution in venomous snakes (Elapidae). \emph{Royal Society Open Science}, \bold{3}, 150277.
	
	Lopez-Vaamonde, C., J. Y. Rasplus, G. D. Weiblen, and J. M. Cook (2001) Molecular phylogenies of fig wasps: Partial cocladogenesis of pollinators and parasites. \emph{Molecular Phylogenetics and Evolution}, \bold{21}, 55-71.

	Mahler, D. L, L. J. Revell, R. E. Glor, and J. B. Losos (2010) Ecological opportunity and the rate of morphological evolution in the diversification of Greater Antillean anoles. \emph{Evolution}, \bold{64}, 2731-2745.
	
	Near, T. J., C. M. Bossu, G. S. Bradburd, R. L. Carlson, R. C. Harrington, P. R. Hollingsworth Jr., B. P. Keck, and D. A. Etnier. (2011) Phylogeny and temporal diversification of darters (Percidae: Etheostomatinae). \emph{Systematic Biology}, \bold{60}, 565-595.
	
	Poulakakis, N., J. M. Miller, E. L. Jensen, L. B. Beheregaray, M. A. Russello, S. Glaberman, J. Boore, and A. Caccone. (2020) Colonization history of Galapagos giant tortoises: Insights from mitogenomes support the progression rule. \emph{Journal of Zoological Systematics and Evolutionary Research}, \bold{58}, 1262-1275.
	
	Revell, L. J., and D. C. Collar (2009) Phylogenetic analysis of the evolutionary correlation using likelihood. \emph{Evolution}, \bold{63}, 1090-1100.

	Revell, L. J., K. S. Toyama, and D. L. Mahler (2022) A simple hierarchical model for heterogeneity in the evolutionary correlation a phylogenetic tree. \emph{PeerJ}, \bold{10}, e13910.
	
	Steeman, M. E., M. B. Hebsgaard, R. E. Fordyce, S. Y. W. Ho, D. L. Rabosky, R. Nielsen, C. Rahbek, H. Glenner, M. V. Sorensen, and E. Willerslev (2009) Radiation of extant cetaceans driven by restructuring of the oceans. \emph{Systematic Biology}, \bold{58}, 573-585.
}
\keyword{datasets}
