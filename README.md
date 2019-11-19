Valpred program to predict transmembrane segments from amino acid sequence
==========================================================================

This program takes a FASTA amino acid sequence on input, and predicts the transmembrane segments.
The prediction algorithm winds the sequence into a helix, and calculates accessible surface areas,
assigning 'environments' consisting of how hydrophobic or hydrophilic the surface area is.
Areas that are compatible with hydrophobic solvent and incompatible with hydrophilic water
are assigned as being transmembrane segments (TMS).
There are 2 environment categories calculated for each residue in a sequence.
"Profiles 3D" is how compatible the residue is with a hydrophilic water environment.
"REPIMPS" is how compatible the resiude is with a hydrophobic solvent environment.

This perl program can be run on the command line in batch mode,
or on a web page as a cgi script on a web server.

The Valpred program was first implemented as a C++ program by Lawrence K. Lee.  
This perl program is a perl implementation of the C++ program.  
Valpred was presented in [Rath et al. 2013](https://www.ncbi.nlm.nih.gov/pubmed/?term=23530628).  
The REPIMPS algorithm was presented by [Dastmalchi et al. 2001](https://www.ncbi.nlm.nih.gov/pubmed/11468350).  
The Profile 3D algorithm was presented by Bowie et al. 1991 and Lüthy et al. 1992.  

### Citations

[A benchmark server using high resolution protein structure data, and benchmark results for membrane helix predictions.](https://www.ncbi.nlm.nih.gov/pubmed/?term=23530628)  
Rath EM, Tessier D, Campbell AA, Lee HC, Werner T, Salam NK, Lee LK, Church WB.  
BMC Bioinformatics. 2013 Mar 27;14:111. doi: 10.1186/1471-2105-14-111.  
PMID: [23530628](https://www.ncbi.nlm.nih.gov/pubmed/23530628) PMCID: PMC3620685 DOI: [10.1186/1471-2105-14-111](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-111)  

[Modeling of the structural features of integral-membrane proteins reverse-environment prediction of integral membrane protein structure (REPIMPS).](https://www.ncbi.nlm.nih.gov/pubmed/11468350)  
Dastmalchi S, Morris MB, Church WB.  
Protein Sci. 2001 Aug;10(8):1529-38.  
PMID: [11468350](https://www.ncbi.nlm.nih.gov/pubmed/11468350) PMCID: PMC2374085 DOI: [10.1110/ps.6301](https://onlinelibrary.wiley.com/doi/full/10.1110/ps.6301)  

[A method to identify protein sequences that fold into a known three-dimensional structure.](https://www.ncbi.nlm.nih.gov/pubmed/1853201)  
Bowie JU, Lüthy R, Eisenberg D.  
Science. 1991 Jul 12;253(5016):164-70.  
PMID: [1853201](https://www.ncbi.nlm.nih.gov/pubmed/1853201) DOI: [10.1126/science.1853201](https://doi.org/10.1126/science.1853201)  

[Assessment of protein models with three-dimensional profiles.](https://www.ncbi.nlm.nih.gov/pubmed/1538787)  
Lüthy R, Bowie JU, Eisenberg D.  
Nature. 1992 Mar 5;356(6364):83-5.  
PMID: [1538787](https://www.ncbi.nlm.nih.gov/pubmed/1538787) DOI: [10.1038/356083a0](https://www.nature.com/articles/356083a0)

