# Coral metabolism and contaminant loads track human land use
 <br />

 
Zachary A. Quinlan, Wesley J. Sparagon, Irina Koester, Sean Swift, Jacqueline Comstock, Rayna McClintock, Teegan Innis, Milou G.I. Arts,
Andreas Haas, Lihini Aluwihare, Craig A. Carlson, Hollie M. Putnam, Katie L. Barott, Pieter C. Dorrestein, Craig E. Nelson, Linda Wegley Kelly
*Corresponding author <br />

### Abstract:
Stress biomarkers released by coral into the surrounding water prior to and during bleaching can be used for the development of non-invasive monitoring tests of coral heat stress. With coral reef managers lacking the tools to accurately predict coral stress in response to elevated temperatures, seawater-based diagnostic tools may facilitate better assessment of reef-specific stress responses. Here, we used high-resolution tandem mass spectrometry to discover exometabolites that are exuded into seawater by thermally stressed corals. During a controlled mesocosm experiment, we observed a pulse of dissolved organic carbon being released from unbleached corals exposed to acute heat stress. The composition of this organic matter pool, which includes coral exometabolites, was profiled to evaluate putative biomarkers of heat stress. These stress biomarkers were further validated in situ during a marine heatwave and bleaching event in Hawaiʻi. A candidate biomarker that was abundant in the mesocosm experiment was also enriched adjacent to bleaching corals in situ and found in six public datasets comprising stressed coral exometabolites or tissue metabolites from both the Caribbean and Pacific. The characterization of this broad coral heat stress biomarker constitutes a critical first step for developing monitoring tools to accelerate coral reef management responses.
<br />


### Funding sources:
National Science Foundation Graduate Research Fellowship grant 2019252845 (ZAQ) <br />
National Science Foundation Postdoctoral Research Fellowship grant 2308400 (ZAQ) <br />
National Science Foundation grant 2516814 (ZAQ) <br />
National Science Foundation grant OCE-2023298 (C.E.N.) <br />
National Science Foundation OCE-2118618 (L.W.K) <br />
National Science Foundation OCE-2023509 (L.I.A & P.C.D.)  <br />
National Science Foudnation OCE-1637396 (Moorea Coral Reef Long-Term Ecological Research project) <br />
MGIA received support from Treub-Maatschappij, the Society for the Advancement of Research in the Tropics. <br />

We would like to acknowledge Christopher Suchocki and Sean Mahaffey for their assistance during fieldwork during depth-profile sampling in Kāne‘ohe Bay. <br />


*****
## Repository Overview:
In this repository you will find the code used to analyze this data and all raw data <br /> <br />
We have tried to clearly index each script and file below but if you have any questions please email me: zquinlan@gmail.com <br /> <br />

All raw data is in the subdirectory data/raw/

### Requirements:
R - Code was written and tested in 3.6.1 <br />
Libraries - tidyverse, broom, data.table, DescTools, readxl, vegan, ape, lme4, MuMIn, RColorBrewer <br />

### Installation and running the code:
1) Clone the repo to your desktop
2) install R and the required libraries using the install.libraries() function
3) Open the .R file in the code directory
4) The locations of all the files are hard coded in the second section of the code "Loading -- Datasets". Change the working directory at the beginning of the script to the location where you have cloned the respository manually
