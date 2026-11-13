This repository contains workflows for recalculating values provided in the BioShifts v2 database and BioShiftR R package.

Use workflow_script.R to sequentially run the following workflows on data from one example source article:
1. standardization of range shift rates as reported in source articles
2. taxonomic harmonization of species names as provided in source articles
3. the calculation of spatial summary values using the target polygon
4. the calculation of poleward climate velocity using the target polygon and gridded climate data
