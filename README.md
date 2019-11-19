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
Valpred was presented in Rath et al. 2013.  
The REPIMPS algorithm was presented by Dastmalchi et al. 2001.  
The Profile 3D algorithm was presented by Bowie et al. 1991 and Lüthy et al. 1992.  

### Citations

A benchmark server using high resolution protein structure data, and benchmark results for membrane helix predictions.  
Rath EM, Tessier D, Campbell AA, Lee HC, Werner T, Salam NK, Lee LK, Church WB.  
BMC Bioinformatics. 2013 Mar 27;14:111. doi: 10.1186/1471-2105-14-111.  
PMID: [23530628](https://www.ncbi.nlm.nih.gov/pubmed/23530628) PMCID: PMC3620685 DOI: [10.1186/1471-2105-14-111](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-111)  

Modelling of the structural features of integral-membrane proteins using REPIMPS (Reverse-Environment Prediction of Integral Membrane Protein Structure)  
Siavoush Dastmalchi, Michael B. Morris and W. Bret Church (2001)  
Protein Science, 10, 1529-1538  

A Method to Identify Protein Sequences that Fold into a Known Three-Dimensional Structure  
James U. Bowie, Roland Lüthy, David Eisenberg (1991)  
Science, 253, 164-170  

Assessment of protein models with three-dimensional profiles  
Roland Lüthy, James U. Bowie, David Eisenberg (1992)  
Nature, 356, 83-85  

