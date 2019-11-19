#!/usr/bin/perl   -wT

# Valpred Windows Program Copyright © 2010 Lawrence K. Lee. All rights reserved.
# Valpred Perl Program Copyright © 2010 Emma M. Rath. All rights reserved.
# Soon this program will be released under an open source software license such as GNU General Public License or 
# Creative Commons license for Free Software Foundation's GNU General Public License at creativecommons.org

# This program takes a FASTA amino acid sequence on input, and predicts the transmembrane segments.
# The prediction algorithm winds the sequence into a helix, and calculates accessible surface areas,
# assigning 'environments' consisting of how hydrophobic or hydrophilic the surface area is.
# Areas that are compatible with hydrophobic solvent and incompatible with hydrophilic water
# are assigned as being transmembrane segments (TMS).
# There are 2 environment categories calculated for each residue in a sequence.
# "Profiles 3D" is how compatible the residue is with a hydrophilic water environment.
# "REMIMPS" is how compatible the resiude is with a hydrophobic solvent environment.

# There are 2 ways that this program can take the PROFILES3D and REPIMPS lines to predict where TMS are.
# Both methods use the moving average of the difference between the two lines.
# The method referred to as VALPRED or VALPRED1 takes 1 moving average, and has various limits on TMS size.
# The method referred to as VALPRED2 takes 2 moving averages of hydrophobicity,
# and then subtracts any areas of hydrophilicity inside that.
# VALPRED1 is the default, and performs better on the Rost/Kernytsky benchmark server.
# VALPRED2 is better at finding short helices belonging to potassium channels,
# but also tends to predict the short hydrophobic helices of soluble proteins as being TMS.

# This perl program can be run on the command line in batch mode,
# or on a web page as a cgi script on a web server.

# References :
#
#      Lawrence K. Lee, Noeris K. Salam, Hong Wing Lee, Emma M. Rath and W. Bret Church (in preparation)
#      'Transmembrane helix analysis using threading approaches'
#
#      James U. Bowie, Roland Lüthy, David Eisenberg (1991)
#      'A Method to Identify Protein Sequences that Fold into a Known Three-Dimensional Structure'
#      Science, 253, 164-170
#
#      Roland Lüthy, James U. Bowie, David Eisenberg (1992)
#      'Assessment of protein models with three-dimensional profiles'
#      Nature, 356, 83-85
# 
#      Siavoush Dastmalchi, Michael B. Morris and W. Bret Church (2001)
#      'Modelling of the structural features of integral-membrane proteins using REPIMPS
#      (Reverse-Environment Prediction of Integral Membrane Protein Structure)' 
#      Protein Science, 10, 1529-1538

# examples of calling this program in command line mode, 
# with 1 or more sequences, 1 line per sequence in the special valpred format, but each sequence has only 1 chain :
# perl perl_valpred_2d.pl -infile 1HML_human_ala_zinc_binding.pdb.fasta.1line
# perl perl_valpred_2d.pl -infile KCSA_upto50.fasta.1line -options valpred_optional_options_file_1.txt
# perl perl_valpred_2d.pl -infile KCSA_upto50.fasta.1line -options valpred_optional_options_file_3.txt
# perl perl_valpred_2d.pl -infile AAO76171.fasta.1line
# perl perl_valpred_2d.pl -infile Hol2438.fasta.1line
# perl perl_valpred_2d.pl -infile KCSA.fasta.1line -options valpred_optional_options_file_2.txt
# perl perl_valpred_2d.pl -infile bacteriorhodopsin.fasta.1line
# perl perl_valpred_2d.pl -infile cytochromeCC3.fasta.1line
# perl perl_valpred_2d.pl -infile HLAclassI.fasta.1line
# perl perl_valpred_2d.pl -infile IRE1.fasta.1line
# perl perl_valpred_2d.pl -infile 1BL8_kchan_slividans.fasta.1line
# perl perl_valpred_2d.pl -infile a_few.fasta.1line
# perl perl_valpred_2d.pl -infile a_few.fasta.1line -multiplegraphs
# perl perl_valpred_2d.pl -infile wang2000_holins.1linefasta -options 2outof3x3.txt -multiplegraphs
# perl perl_valpred_2d.pl -fastafile KCSA_STRLI_v1.fasta -options valpred_options_valpred1.txt
# perl perl_valpred_2d.pl -fastafile KCSA_STRLI_v2.fasta -options valpred_options_valpred2.txt
# perl perl_valpred_2d.pl -fastafile pdb_1BL8_kchan_chainA.pdb.fasta -tms pdb_1BL8_kchan_chainA.pdb.tms_structure.txt valpred_options_valpred1.txt
# perl perl_valpred_2d.pl -fastafile pdb_1BL8_kchan_chainA.pdb.fasta -tms pdb_1BL8_kchan_chainA.pdb.tms_structure.txt valpred_options_valpred2.txt
# perl perl_valpred_2d.pl -infile wang2000_holins_v1.1linefasta -multiplegraphs -options valpred_options_valpred1.txt
# perl perl_valpred_2d.pl -infile wang2000_holins_v2.1linefasta -multiplegraphs -options valpred_options_valpred2.txt
# perl perl_valpred_2d.pl -infile wang2000_holins_v1.1linefasta -multiplegraphs -version valpred1
# perl perl_valpred_2d.pl -infile wang2000_holins_v2.1linefasta -multiplegraphs -version valpred2

# examples of calling this program in command line mode, 
# with 1 sequence in standard fasta format, with 1 or more chains, each chain being a new fasta sequence,
# and can extend over 1 line:
# perl perl_valpred_2d.pl -fastachainfile KCSA.fasta -debug
# perl perl_valpred_2d.pl -fastachainfile 1BL8_kchan_slividans.fasta.txt
# perl perl_valpred_2d.pl -fastachainfile 2OAR_prokaryotic_mscl.fasta.txt

# examples of calling this program in command line mode, 
# with 1 or more sequences in standard fasta format which can extend over 1 line, each sequence with only 1 chains:
# perl perl_valpred_2d.pl -fastafile a_few.fasta
# perl perl_valpred_2d.pl -fastafile a_few.fasta -multiplegraphs
# perl perl_valpred_2d.pl -fastafile a_few2.fasta -multiplegraphs
# perl perl_valpred_2d.pl -fastafile some_kchans.fasta
# perl perl_valpred_2d.pl -fastafile some_kchans.fasta -multiplegraphs
# perl perl_valpred_2d.pl -fastafile pdb_1BL8_kchan_chainA.pdb.fasta -tms pdb_1BL8_kchan_chainA.pdb.tms_structure.txt
# perl perl_valpred_2d.pl -fastafile test1.fa -benchmark
# perl perl_valpred_2d.pl -fastafile rost_tmh_benchmark_data_tmheval.fa -benchmark

# examples of calling this program in command line mode, 
# with 1 or more sequences for which the valpred results have already been calculated,
# to regraph the results, perhaps with different graph options to the last time the graphs were produced :
# perl perl_valpred_2d.pl -graphpointsfile a_few.fasta.valpred_graphpoints.txt
# perl perl_valpred_2d.pl -graphpointsfile a_few2.fasta.valpred_graphpoints.txt -options valpred_optional_options_file_4.txt
# perl perl_valpred_2d.pl -graphpointsfile a_few2.fasta.valpred_graphpoints.txt -options valpred_optional_options_file_4.txt -multiplegraphs
# perl perl_valpred_2d.pl -graphpointsfile some_kchans.fasta.valpred_graphpoints.txt -options valpred_optional_options_file_4.txt -multiplegraphs

# examples of calling this program in command line mode,
# with the valpred calculation already done
# and doing the evaluate calculation to benchmark how good the valpred prediction is :
# perl perl_valpred_2d.pl -predictions 1BL8_Kchan_pdb_chainA_valpred2.fasta.valpred_segments.txt -tms 1BL8_Kchan.pdb.chainA.pdb_structure.txt
# perl perl_valpred_2d.pl -predictions many_valpred2.fasta.valpred_segments.txt -tms many.observed_tms.txt
# perl perl_valpred_2d.pl -predictions many_valpred1.fasta.valpred_segments.txt -tms many.observed_tms.txt
# perl perl_valpred_2d.pl -benchmarkpredictions many_tmhmm2.fasta.1line.benchmark.txt -tms many.observed_tms.txt
# perl perl_valpred_2d.pl -version valpred2 -inscoresfile wang2000_seq10.valpred2_scores.txt

# PLEASE NOTE THAT THE AMPHIPATHIC ALGORITHM IS NOT YET PROPERLY IMPLEMENTED IN THIS VERSION OF VALPRED.
# PLEASE DO NOT USE THIS VERSION OF VALPRED FOR PREDICTING AMPHIPATHIC HELICES.

# example of calling this program in command line mode, to predict amphipathic helices :
# perl perl_valpred_2d.pl -amphipathic -fastafile 2J58_A_TMH.fasta
# perl perl_valpred_2d.pl -amphipathic -fastafile 2J58_A.fasta

# in command line mode, 
#	this program will produce the following output files :
#  - output file xxx.valpred_summary.txt (eg. KCSA.fasta.1line.valpred_summary.txt)
#       which contains the valpred transmembrane summary for the first sequence in the input file
#  - output file xxx.valpred_graph.gif (eg. KCSA.fasta.1line.valpred_graph.gif)
#       which contains the valpred transmembrane graph for the first sequence in the input file
#  - output file xxx.valpred_calculations.txt (eg. KCSA.fasta.1line.valpred_calculations.txt)
#       which contains the valpred calculation results of accessibility and hydropathy environment for each residue
#  - output file xxx.valpred_segments.txt (eg. KCSA.fasta.1line.valpred_segments.txt)
#       which is the valpred transmembrane segments output file with one line per input sequence
#  - output file xxx.valpred_scores.txt (eg. KCSA.fasta.1line.valpred_segments.txt)
#       which is the valpred transmembrane scores output file with one line per input sequence,
#	containing the PROFILES3D and REPIMPS scores for each residue of the segment
# if the input file is -graphpointsfile, 
#	then only the xxx.valpred_graph.gif file is produced,
#	and xxx.valpred_summary.txt, xxx.valpred_calculations.txt, and xxx.valpred_segments.txt are not produced.

# in command line mode, if the -multiplegraphs option is present on the command line,
#	then this program will produce the following additional output files :
#  - output file xxx.valpred_graphpoints.txt (eg. KCSA.fasta.1line.valpred_graphpoints.txt)
#       which is the valpred transmembrane segments and graph points output file with one line per input sequence
#  - output file xxx.valpred_graphs.html (eg. KCSA.fasta.1line.valpred_valpred_graphs.html)
#       which is an html page containing multiple valpred graphs, 1 graph per input sequence
#  - output files xxx.valpred_graph_aaa.gif, xxx.vapred_graph_bbb.gif, etc.
#       which are the graph output files, one for each input sequence, 
#       where aaa, bbb, etc are the sequence-ids from the input file.
# if the input file is -graphpointsfile, 
#	then only the xxx.valpred_graphs.html and xxx.vapred_graph_aaa.gif files are produce,
#	 and xxx.valpred_graphpoints.txt is not produced.

# in command line mode, if the -benchmark option is present on the command line :
#  - output file xxx.valpred_benchmark.txt (eg. KCSA.fasta.1line.valpred_benchmark.txt)
#       which is the output file of TMS helices for the Rost TMH benchmark server at http://cubic.bioc.columbia.edu/services/tmh_benchmark/

# in command line mode, if the -tms option is present on the command line :
#  - output file xxx.valpred_evaluate.txt (eg. KCSA.fasta.1line.valpred_evaluate.txt)
#       which compares the position of TMS in the tms file to that predicted by valpred,
#	one line in the file per input protein.
#  - output file xxx.valpred_evaluate_summary.txt (eg. KCSA.fasta.1line.valpred_evaluate_summary.txt)
#	which gives a summary of the comparison of TMS in the tms file to that predicted by valpred
#	one summary for all the input proteins.

# in command line mode, if the -amphipathic option is present on the command line,
#	then this program will produce the following additional output files :
#  - output file xxx.valpred_amphipathic_graph.gif,
#       which is the graph output file is the (first sequence's) predicted amphipathic helices. 
#  - output file xxx.valpred_amphipathic_segments.txt
#  - output file xxx.valpred_amphipathic_residues.txt
# Note that the predicted amphipathic segments can be overlapping each other.

# in command line mode, if the -debug option is present on the command line :
#  - output file valpred.pdb, which is pdb file of helix built from sequence

# If this program is called with the environment variable $ENV{'REQUEST_METHOD'} set, 
# then it will assume that it is being run as a CGI script on a web server.
# If there are no input form fields, then the output will be an HTML form.
# If there are input form fields, then the output will be the HTML results screen.
# If the input form fields are a valpred request, then this program will run the valpred calculation
# and the HTML output will be the summary text results, the residues table, and a link to the valpred graph image result.
# If the input form fields are a valpred result, then the HTML output will be the valpred graph image result.

# If this program is called with $ENV{'REQUEST_METHOD'} not set, 
# then it will assume that it is being run from the command line.
# If an appropriate input file is not provided, 
# then the program will output instructions to the screen on how to run this program from the command line.

# If this program is being run as a CGI script on a web server, then only 1 valpred result will be calculated,
# and multiple fasta input will be considered to be multiple chains of the one protein.

# If this program is being run as a CGI script on a web server, then the input must be in standard fasta format, 
# with an optional first line starting with a '>' identifier, followed by one or more sequence lines.
# If this program is being run in command line mode, then the input must be in standard fasta format, in the -fastachainfile or in the -fastafile,
# or in the special valpred format that contains 1 line per protein, with ':::::' between fields, in the -infile,
# or in the special valpred graph_points formats, with ':::::' between fields, in the -graphpointsfile.

# Example of input file in the special valpred format :
# AAO76171:::::>gi|29338370|gb|AAO76171.1| Rhodopsin-like GPCR superfamily [Bacteroides thetaiotaomicron VPI-5482]:::::MWWKQGSKKRNDTIRYRYILPYESWMDDARVDVQRDECGCGEIQLMDVEPLGDIELERILVPYVVTPFFAYLQPKAEEVKSRDIQAECFLDFEVNKINIRPEYMNNPKELAKIRAMIDELKSDPSIKVNKLDIVGYASPEGSLANNKRLSEGRAMALRDYLASRYDFSRNQYYIIFGGENWDGLVKALDTIDFEYKDEALNIINDIPVEKGREAKLMQLRGGVPYRYMLKYIFPSLRVAICKVNYEIKNFNLDEAKEIIKTRPQNLSLNEMFMVANSYPKGSQEFIDVFETAVRMYPKDEIASINAAAAALSRNDLVSAERYLNMVNVNKQLPEYSNAMGVLMLLKGEYEHAEEYLKAAAKSGLQAAGQNLEELAKKKTNAAEIEKIENRDK
# AAA26845:::::>gi|153550|gb|AAA26845.1| penicillin-binding protein (PBPB2):::::AVIASISKEMPGISISTSWDRKILETSLSSIVGSVSSEKAGLPAEEAETYLKKGYSLNDRVGTSYLEKQYEETLQGKRSVKEIHLDKYGNMESVENIEDGTKGNNIKLTIDLSFQDSVDALLKSYFNSELGNGGAKYSEGVYAVALNPKTGAVLSMSGIKHDLKTGELTPDSLGTVTNVFVPGSVVKAATISSGWENGVLSGNQTLTDQSIVFQGSAPINSWYTAFSRPMPITAVQALEYSSNAYMVQTALGLMGQTYQPNMFVGTSNLESAMGKLRSTFGEYGLGSATGIDLPDESTGFIPKEYSFANFITNAFGQFDNYTPMQLAQYVATIANDGVRVAPRIVEGIYGNNDKGGLGGLIQQLQPTEMNKVNISDSDMSVLHQGFYQVAHGTSGLTTGRAFSNGAAVSISGKTGTAESYVAGGQEANNTNAVAYAPSDNPQIAVAVVFPHNTNLTNGVGPSIARDIINLYNQHHPMN
# AAA26872:::::>gi|153613|gb|AAA26872.1| DpnII DNA methylase:::::MKIKEIKKVTLQPFTKWTGGKRQLLPVIRELIPKTYNRYFEPFVGGGALFFDLAPKDAVINDFNAELINCYQQIKDNPQELIEILKVHQEYNSKEYYLDLRSADRDERIDMMSEVQRAARILYMLRVNFNGLYRVNSKNQFNVPYGRYKNPKIVDEELISAISVYINNNQLEIKVGDFEKAIVDVRTGDFVYFDPPYIPLSETSAFTSYTHEGFSFADQVRLRDAFKRLSDTGAYVMLSNSSSALVEELYKDFNIHYVEATRTNGAKSSSRGKISEIIVTNYEK

# If this program is being called as a cgi-script, 
# then multiple fasta input will be treated as multiple chains, eg.
# 
# >2OAR:E|PDBID|CHAIN|SEQUENCE
# MGHHHHHHHHHHSSGHIDDDDKHMLKGFKEFLARGNIVDLAVAVVIGTAFTALVTKFTDSIITPLINRIGVNAQSDVGIL
# RIGIGGGQTIDLNVLLSAAINFFLIAFAVYFLVVLPYNTLRKKGEVEQPGDTQVVLLTEIRDLLAQTNGDSPGRHGGRGT
# PSPTDGPRASTESQ
# >2OAR:D|PDBID|CHAIN|SEQUENCE
# MGHHHHHHHHHHSSGHIDDDDKHMLKGFKEFLARGNIVDLAVAVVIGTAFTALVTKFTDSIITPLINRIGVNAQSDVGIL
# RIGIGGGQTIDLNVLLSAAINFFLIAFAVYFLVVLPYNTLRKKGEVEQPGDTQVVLLTEIRDLLAQTNGDSPGRHGGRGT
# PSPTDGPRASTESQ
# >2OAR:C|PDBID|CHAIN|SEQUENCE
# MGHHHHHHHHHHSSGHIDDDDKHMLKGFKEFLARGNIVDLAVAVVIGTAFTALVTKFTDSIITPLINRIGVNAQSDVGIL
# RIGIGGGQTIDLNVLLSAAINFFLIAFAVYFLVVLPYNTLRKKGEVEQPGDTQVVLLTEIRDLLAQTNGDSPGRHGGRGT
# PSPTDGPRASTESQ
# >2OAR:B|PDBID|CHAIN|SEQUENCE
# MGHHHHHHHHHHSSGHIDDDDKHMLKGFKEFLARGNIVDLAVAVVIGTAFTALVTKFTDSIITPLINRIGVNAQSDVGIL
# RIGIGGGQTIDLNVLLSAAINFFLIAFAVYFLVVLPYNTLRKKGEVEQPGDTQVVLLTEIRDLLAQTNGDSPGRHGGRGT
# PSPTDGPRASTESQ
# >2OAR:A|PDBID|CHAIN|SEQUENCE
# MGHHHHHHHHHHSSGHIDDDDKHMLKGFKEFLARGNIVDLAVAVVIGTAFTALVTKFTDSIITPLINRIGVNAQSDVGIL
# RIGIGGGQTIDLNVLLSAAINFFLIAFAVYFLVVLPYNTLRKKGEVEQPGDTQVVLLTEIRDLLAQTNGDSPGRHGGRGT
# PSPTDGPRASTESQ
# 
# will produce :
# 
# Predicted Transmembrane Domains:
#            Chain                   Residues
#              1                     30  -  58
#              1                     86  -  112
#              2                     208  -  232
#              2                     260  -  286
#              3                     440  -  460
#              4                     556  -  580
#              4                     608  -  634
#              5                     782  -  808

# Example of the valpred_scores.txt file, produced on output, and can be read on input in the command line -inscoresfile :
# seq3:::::>APSE-1 P11, gamma family, class 1, from bacteriophages of gram-negative host:::::MSEFCKPLLDILRHQGTCAALAFIMALLRARYHRKDFYRSLLDALMCAMLGGVAHELLQFLGLKADYSWLASVAIGYLGVDRIGNWLKKKTGKL:::::PROFILES3D;;;;;-0.9,0.16,0.6,-1.35,-0.17,0.07,-0.25,-0.3,-0.46,0.29,-0.59,-0.46,0.56,-0.06,0.62,0.63,-0.2,-0.17,0.44,0.44,-0.46,0.44,-1.35,-0.06,-0.27,0.44,-0.3,-0.46,0.56,0.44,-0.11,-0.55,-0.06,-0.51,0.07,0.29,-1.35,-0.55,0.56,0.16,-0.3,-0.3,0.44,0.44,-0.46,-0.27,-0.17,0.44,-0.9,-1.37,0.63,0.63,-0.62,0.44,-0.06,0.2,-0.46,-0.46,0.29,-0.85,-0.46,0.63,-0.46,-0.5,0.44,0.44,-0.55,0.16,-1.09,-0.46,0.44,0.16,0.3,0.44,-0.59,0.63,-0.55,-0.46,0.63,-1.25,-0.28,0.56,-0.59,0.63,-0.02,0.86,-0.46,0.07,0.07,-0.5,-0.2,0.63,0.07,-1.37,:::::REPIMPS;;;;;1.26,0.47,-2.15,1.28,0.95,-1.37,-1.56,1.3,1.3,-0.28,1.11,1.3,-1.8,-0.34,-1.38,-0.46,0.39,0.95,0.76,0.76,1.3,0.76,1.28,1.11,1.26,0.76,1.3,1.3,-1.8,0.76,-1.8,0.27,-0.34,-1.8,-1.37,-0.28,1.28,0.27,-1.8,0.47,1.3,1.3,-0.28,0.76,1.3,1.26,0.95,0.76,1.26,1.3,-0.46,-0.46,0.74,0.76,-0.34,-2.15,1.3,1.3,-1.38,1.28,1.3,-0.46,1.3,-1.37,0.76,-0.28,0.27,0.47,1.11,1.3,0.76,0.47,0.74,0.76,1.11,-0.46,0.27,1.3,-0.46,0.74,-0.28,-1.8,1.11,-0.46,-1.76,1.11,1.3,-1.37,-1.37,-1.37,0.39,-0.46,-1.37,1.3,:::::
# seq4:::::>ES18 13, gamma family, class 1 holin, from bacteriophages of gram-negative host:::::MKKMPEKHDLLTAMMAAKEQGIGAILAFAMAYLRGRYNGGAFKKTLIDATMCAIIAWFIRDLLVFAGLSSNLAYIASVFIGYIGTDSIGSLIKRFAAKKAGVDDANQQ:::::PROFILES3D;;;;;-0.9,0.07,0.07,-0.27,-0.25,0.2,0.07,0.17,0.44,-1.37,-0.46,-0.2,0.44,-0.27,-0.27,0.44,0.44,0.07,0.2,0.29,0.63,-0.59,0.63,0.44,-2.36,-0.46,0.44,-1.81,0.44,-0.27,0.44,-0.55,-0.46,0.56,0.63,-0.51,-1.7,0.32,0.63,0.63,0.44,-1.35,0.07,0.07,-0.2,-0.3,-0.59,0.44,0.44,-0.2,-0.27,-0.17,0.44,-0.59,-0.59,0.44,-1.09,-0.85,-0.06,-0.51,-0.28,-0.3,-0.46,-1.25,-1.35,0.44,0.63,-0.46,0.47,0.16,-0.02,-0.46,0.44,-0.55,-0.59,0.44,0.16,-1.25,-1.35,-2.36,0.63,-0.55,-0.06,0.63,-0.2,-0.28,0.16,-0.59,0.63,-0.38,-0.46,-0.59,0.07,-0.51,-0.85,0.44,0.44,0.07,-0.5,0.44,0.63,-1.25,0.44,0.44,0.44,0.32,0.29,0.29,:::::REPIMPS;;;;;1.26,-1.37,-1.37,1.26,-1.56,-2.15,-1.37,-0.34,-0.28,1.3,1.3,0.39,0.76,1.26,1.26,0.76,0.76,-1.37,-2.15,-1.38,-0.46,1.11,-0.46,0.76,1.11,1.3,0.76,1.28,0.76,1.26,0.76,0.27,1.3,-1.8,-0.46,-1.8,0.27,-1.76,-0.46,-0.46,0.76,1.28,-1.37,-1.37,0.39,1.3,1.11,-0.28,0.76,0.39,1.26,0.95,0.76,1.11,1.11,0.76,1.11,1.28,1.11,-1.8,-0.28,1.3,1.3,0.74,1.28,0.76,-0.46,1.3,0.47,0.47,-1.76,1.3,0.76,0.27,1.11,0.76,0.47,0.74,1.28,1.11,-0.46,0.27,1.11,-0.46,0.39,-0.28,0.47,1.11,-0.46,0.47,1.3,1.11,-1.37,-1.8,1.28,0.76,0.76,-1.37,-1.37,0.76,-0.46,0.74,-0.28,-0.28,0.76,-1.76,-1.38,-1.38,:::::
# seq5:::::>HK022 S, gamma family, class 1 holin, from bacteriophages of gram-negative host:::::MKMPEKNDLLAAILAAKEQGIGAILAFAMAYLRGRYNGGAFTKTVIDATMCAIIAWFIRDLLGFAGLSSNLAYITSVFIGYIGTDSIGSLIKRFAAKKAGVEDGGNQ:::::PROFILES3D;;;;;-0.9,0.07,-0.9,-0.25,0.2,0.66,-0.02,0.44,-0.46,-0.46,0.44,0.44,-0.06,-0.46,0.44,0.44,0.66,0.2,0.29,0.63,-0.59,0.63,0.44,-2.36,-0.46,0.44,-1.81,0.44,-0.27,0.44,-0.55,-0.46,0.56,0.63,-0.51,-1.7,0.32,0.63,0.63,0.44,-1.81,-0.2,0.07,-0.2,-0.62,-0.59,0.44,0.44,-0.2,-0.27,-0.17,0.44,-0.59,-0.59,0.44,-1.09,-1.35,-0.06,-0.51,-0.28,-0.3,-1.37,0.63,-1.35,0.44,0.63,-0.46,0.47,0.16,-0.02,-0.46,0.44,-0.55,-0.59,-0.2,0.16,-1.25,-1.35,-0.59,0.63,-0.55,-0.06,0.63,-0.2,-0.28,0.16,-0.59,0.63,-0.38,-0.46,-2.36,0.07,-0.51,-0.85,0.44,0.44,0.07,-0.5,0.44,0.63,-1.25,0.6,0.44,0.63,0.63,0.32,0.29,:::::REPIMPS;;;;;1.26,-1.37,1.26,-1.56,-2.15,-1.37,-1.76,-0.28,1.3,1.3,0.76,0.76,1.11,1.3,0.76,0.76,-1.37,-2.15,-1.38,-0.46,1.11,-0.46,0.76,1.11,1.3,0.76,1.28,0.76,1.26,0.76,0.27,1.3,-1.8,-0.46,-1.8,0.27,-1.76,-0.46,-0.46,0.76,1.28,0.39,-1.37,0.39,0.74,1.11,-0.28,0.76,0.39,1.26,0.95,0.76,1.11,1.11,0.76,1.11,1.28,1.11,-1.8,-0.28,1.3,1.3,-0.46,1.28,0.76,-0.46,1.3,0.47,0.47,-1.76,1.3,0.76,0.27,1.11,0.39,0.47,0.74,1.28,1.11,-0.46,0.27,1.11,-0.46,0.39,-0.28,0.47,1.11,-0.46,0.47,1.3,1.11,-1.37,-1.8,1.28,0.76,0.76,-1.37,-1.37,0.76,-0.46,0.74,-2.15,-0.28,-0.46,-0.46,-1.76,-1.38,:::::

# Example of valpred calculations output file :
# Index,Residue,Chain,No.,Area Buried,Fraction Polar,REPIMPS FP,Category,REPIMPS Cat,Main SASA,Side SASA,SASA,Solvent Score,REPIMPS Score,Overall Score,is transmembrane,
# MET,1,0,11.1328979071137,0.826815642458101,-0.108521954274218,profEnv,repimpsEnv,18.4241942495527,161.036668092886,-0.9,1.26,-0.9,0,
# PRO,1,1,37.9770085093625,0.689542483660131,-0.0323420129374274,profEnv,repimpsEnv,2.34572251468038,98.5742014906374,0.05,-1.56,0.05,0,
# PRO,1,2,39.7063804653386,0.765472312703583,0.05625245639721,profEnv,repimpsEnv,6.64621379159441,96.8448295346614,0.05,-1.56,0.05,0,
# MET,1,3,31.595761389371,0.80225988700565,-0.0142247442597214,profEnv,repimpsEnv,0.78190750489346,140.573804610629,-0.9,1.26,-0.9,0,
# LEU,1,4,51.4199930921647,0.717325227963526,0.0391145508895036,profEnv,repimpsEnv,0.39095375244673,108.373975907835,-0.46,1.3,-0.46,0,

# Example of valpred segments output file :
# 1BL8:A:::::>1BL8:A|PDBID|CHAIN|SEQUENCE:::::ALHWRAAGAATVLLVIVLLAGSYLAVLAERGAPGAQLITYPRALWWSVETATTVGYGDLYPVTLWGRCVAVVVMVAGITSFGLVTAALATWFVGREQ:::::0-20,49-88,:::::
# 2OAR:E:::::>2OAR:E|PDBID|CHAIN|SEQUENCE:::::MGHHHHHHHHHHSSGHIDDDDKHMLKGFKEFLARGNIVDLAVAVVIGTAFTALVTKFTDSIITPLINRIGVNAQSDVGILRIGIGGGQTIDLNVLLSAAINFFLIAFAVYFLVVLPYNTLRKKGEVEQPGDTQVVLLTEIRDLLAQTNGDSPGRHGGRGTPSPTDGPRASTESQ:::::30-58,86-112,:::::
# 1HML:::::>1HML CALCIUM-BINDING PROTEIN:::::KQFTKCELSQLLKDIDGYGGIALPELICTMFHTSGYDTQAIVENNESTEYGLFQISNKLWCKSSQVPQSRNICDISCDKFLDDDITDDIMCAKKILDIKGIDYWLAHKALCTEKLEQWLCEKL:::::-:::::

# Example of valpred graphpoints file :
# (is output file when -multiplegraphs option is present, is input file for -graphpointsfile option) :
# SEQ_1BL8A:::::>1BL8:A|PDBID|CHAIN|1BL8:A:::::ALHWRAAGAATVLLVIVLLAGSYLAVLAERGAPGAQLITYPRALWWSVETATTVGYGDLYPVTLWGRCVAVVVMVAGITSFGLVTAALATWFVGREQ:::::0-20,49-88,:::::PROFILES3D,,,,,-0.518,-0.358333333333333,-0.244285714285714,-0.135,-0.0711111111111112,-0.02,-0.084,-0.072,-0.112,-0.14,-0.059,-0.109,-0.278,-0.478,-0.568,-0.568,-0.485,-0.344,-0.353,-0.262,-0.248,-0.304,-0.225,-0.044,0.062,-0.033,-0.033,-0.00500000000000003,0.025,0.134,0.134,0.225,0.134,0.031,-0.049,-0.053,-0.141,-0.129,-0.06,-0.169,-0.322,-0.46,-0.361,-0.427,-0.45,-0.356,-0.287,-0.363,-0.427,-0.506,-0.334,-0.28,-0.179,-0.01,-0.013,-0.107,-0.176,-0.281,-0.281,-0.202,-0.374,-0.256,-0.37,-0.431,-0.355,-0.256,-0.356,-0.356,-0.461,-0.505,-0.521,-0.54,-0.426,-0.645,-0.695,-0.723,-0.779,-0.591,-0.603,-0.638,-0.533,-0.533,-0.552,-0.362,-0.298,-0.334,-0.262,-0.46,-0.448,-0.26,-0.291,-0.275,-0.29,-0.271111111111111,-0.36,-0.382857142857143,-0.265,:::::REPIMPS,,,,,0.206,0.298333333333333,0.364285714285714,0.26125,0.316666666666667,0.361,0.324,0.268,0.432,0.451,0.705,0.74,0.738,0.914,0.968,0.968,0.883,0.856,0.753,0.753,0.755,0.718,0.774,0.72,0.375,0.119,0.119,0.148,-0.0350000000000001,-0.211,-0.211,-0.423,-0.423,-0.388,-0.134,0.073,-0.037,-0.293,-0.061,0.115,0.15,0.399,0.316,0.279,0.025,0.037,0.269,0.488,0.451,0.395,0.238,0.154,0.061,-0.041,0.304,0.292,0.06,0.095,0.095,0.151,0.308,0.235,0.101,0.224,0.168,0.217,0.447,0.447,0.482,0.478,0.441,0.563,0.697,0.713,0.678,0.649,0.703,0.583,0.639,0.587,0.552,0.552,0.674,0.693,0.73,0.722,0.705,0.879,0.823,0.703,0.484,0.193,-0.021,-0.167777777777778,-0.28375,-0.38,-0.628333333333333,:::::TMS,,,,,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,:::::
# SEQ_2OARE:::::>2OAR:E|PDBID|CHAIN|2OAR:E:::::MGHHHHHHHHHHSSGHIDDDDKHMLKGFKEFLARGNIVDLAVAVVIGTAFTALVTKFTDSIITPLINRIGVNAQSDVGILRIGIGGGQTIDLNVLLSAAINFFLIAFAVYFLVVLPYNTLRKKGEVEQPGDTQVVLLTEIRDLLAQTNGDSPGRHGGRGTPSPTDGPRASTESQ:::::30-58,86-112,:::::PROFILES3D,,,,,-0.044,-0.00833333333333333,0.0171428571428572,0.03625,0.0511111111111111,0.04,0.124,0.055,0.077,0.099,0.145,0.122,0.046,0.073,0.085,0.135,0.113,0.185,0.186,0.08,-0.029,0.043,0.165,-0.014,-0.036,-0.06,-0.167,-0.263,-0.236,-0.197,-0.088,-0.122,-0.244,-0.234,-0.197,-0.354,-0.175,-0.207,-0.207,-0.281,-0.469,-0.737,-0.615,-0.51,-0.51,-0.508,-0.572,-0.466,-0.556,-0.493,-0.388,-0.145,-0.343,-0.284,-0.284,-0.133,-0.119,-0.169,-0.143,-0.106,-0.132,-0.198,-0.065,-0.155,-0.435,-0.388,-0.507,-0.469,-0.405,-0.351,-0.289,-0.186,-0.246,-0.132,0.045,-0.064,0.00999999999999999,-0.258,-0.239,-0.327,-0.28,-0.261,-0.136,-0.17,-0.131,-0.144,-0.049,0.141,0.076,0.01,-0.099,-0.299,-0.346,-0.331,-0.267,-0.267,-0.313,-0.448,-0.627,-0.548,-0.561,-0.38,-0.531,-0.531,-0.7,-0.696,-0.779,-0.644,-0.588,-0.512,-0.483,-0.552,-0.472,-0.518,-0.354,-0.345,-0.311,-0.258,-0.067,-0.034,0.016,-0.084,0.031,0.062,-0.002,0.107,0.202,0.175,0.138,-0.05,-0.195,-0.116,-0.222,-0.264,-0.219,-0.288,-0.276,-0.212,-0.287,-0.208,-0.039,0.036,0.062,0.107,0.15,0.2,0.106,0.037,0.146,0.248,0.221,0.255,0.338,0.362,0.362,0.298,0.311,0.352,0.264,0.195,0.207,0.207,0.119,0.119,0.1,0.136,0.141,0.185,0.226,0.268,0.265555555555556,0.22,0.287142857142857,0.241666666666667,:::::REPIMPS,,,,,-0.044,-0.0933333333333333,-0.128571428571429,-0.155,-0.175555555555556,-0.192,-0.352,-0.34,-0.259,-0.178,-0.19,-0.19,-0.045,-0.039,-0.033,-0.027,-0.021,-0.124,-0.205,-0.126,0.05,-0.053,-0.21,-0.0540000000000001,-0.163,-0.35,-0.194,0.073,0.183,-0.123,-0.299,-0.338,-0.181,-0.235,-0.126,0.219,0.167,0.111,0.111,0.365,0.485,0.772,0.615,0.58,0.684,0.682,0.645,0.647,0.701,0.701,0.666,0.418,0.592,0.592,0.488,0.407,0.479,0.514,0.423,0.193,0.284,0.532,0.228,0.00899999999999999,0.148,0.0550000000000001,0.0180000000000001,-0.269,-0.232,-0.214,-0.297,-0.436,-0.186,-0.052,-0.052,0.124,-0.13,0.157,0.035,0.284,0.191,0.173,0.053,-0.039,-0.111,-0.13,0.022,0.041,-0.089,-0.126,0.05,0.226,0.319,0.533,0.57,0.57,0.422,0.42,0.724,0.78,0.761,0.707,0.788,0.788,0.786,0.702,1.006,1.008,0.954,0.898,0.917,0.685,0.584,0.332,0.297,0.4,0.092,-0.175,-0.386,-0.506,-0.851,-0.621,-0.863,-0.825,-1.02,-1.196,-1.044,-0.868,-0.869,-0.749,-0.46,-0.404,-0.059,0.118,0.059,0.216,0.0640000000000001,-0.00299999999999998,0.265,0.321,0.323,0.0550000000000001,-0.036,-0.251,-0.0819999999999999,-0.221,0.00600000000000001,-0.122,-0.298,-0.608,-0.718,-0.626,-0.711,-0.715,-0.715,-0.648,-0.851,-0.648,-0.758,-0.539,-0.533,-0.533,-0.643,-0.643,-0.521,-0.513,-0.318,-0.58,-0.377,-0.554,-0.584444444444444,-0.6,-0.462857142857143,-0.24,:::::TMS,,,,,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,:::::
# SEQ_3:::::>1HML CALCIUM-BINDING PROTEIN:::::KQFTKCELSQLLKDIDGYGGIALPELICTMFHTSGYDTQAIVENNESTEYGLFQISNKLWCKSSQVPQSRNICDISCDKFLDDDITDDIMCAKKILDIKGIDYWLAHKALCTEKLEQWLCEKL:::::-:::::PROFILES3D,,,,,-0.316,-0.291666666666667,-0.311428571428571,-0.31,-0.317777777777778,-0.257,-0.31,-0.385,-0.197,-0.133,-0.199,-0.138,-0.032,-0.057,0.044,0.078,0.065,0.155,0.102,0.033,0.112,0.022,-0.1,-0.062,-0.145,-0.298,-0.42,-0.47,-0.444,-0.403,-0.36,-0.369,-0.266,-0.269,-0.22,-0.086,0.036,-0.02,0.02,0.036,-0.029,0.046,0.018,0.018,0.009,-0.09,0.032,0.048,-0.107,-0.171,-0.175,-0.179,-0.253,-0.226,-0.292,-0.346,-0.426,-0.314,-0.163,-0.084,-0.049,-0.19,-0.157,-0.135,-0.073,-0.015,1.38777878078145e-18,-0.125,-0.158,-0.161,-0.249,-0.108,-0.1,-0.085,-0.094,-0.178,-0.222,-0.119,-0.058,-0.13,-0.13,-0.166,-0.105,-0.105,-0.118,-0.01,0.019,0.019,0.041,0.076,0.129,0.103,0.103,2.22044604925031e-17,0.072,0.162,0.12,0.12,-0.000999999999999984,-0.137,-0.177,-0.087,-0.222,-0.213,-0.235,-0.344,-0.302,-0.366,-0.291,-0.155,-0.246,-0.27,-0.15,-0.209,-0.283,-0.254,-0.177,-0.15,-0.307,-0.348888888888889,-0.22125,-0.281428571428571,-0.376666666666667,:::::REPIMPS,,,,,-0.49,-0.25,-0.521428571428571,-0.29375,-0.208888888888889,-0.326,-0.059,0.209,-0.056,-0.123,0.125,0.00200000000000002,0.171,0.068,-0.025,0.067,0.048,-0.00599999999999999,0.261,0.133,-0.193,-0.035,0.122,0.19,0.275,0.447,0.464,0.354,0.263,0.466,0.635,0.532,0.393,0.337,0.16,0.11,0.093,0.201,-0.053,-0.276,-0.406,-0.648,-0.573,-0.573,-0.65,-0.699,-0.856,-0.8,-0.457,-0.419,-0.132,0.13,-0.093,-0.269,0.076,0.16,0.301,0.034,-0.047,0.138,-0.111,-0.084,-0.005,-0.00600000000000001,-0.089,-0.38,-0.651,-0.403,-0.355,-0.43,-0.181,-0.208,-0.016,0.094,-0.09,0.218,0.524,0.385,0.262,0.262,0.262,0.254,0.131,0.131,0.379,0.377,0.342,0.446,0.337,0.228,0.228,0.319,0.319,0.458,0.21,0.038,0.054,-0.05,0.114,0.362,0.381,0.327,0.321,0.073,0.286,0.462,0.446,0.513,0.271,0.023,0.023,-0.268,-0.372,-0.124,-0.07,-0.105,-0.415,-0.591,-0.246,-0.121111111111111,-0.29875,-0.0342857142857143,0.19,:::::TMS,,,,,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,:::::
# SEQ_AAO76171:::::>gi|29338370|gb|AAO76171.1| Rhodopsin-like GPCR superfamily [Bacteroides thetaiotaomicron VPI-5482]:::::MWWKQGSKKRNDTIRYRYILPYESWMDDARVDVQRDECGCGEIQLMDVEPLGDIELERILVPYVVTPFFAYLQPKAEEVKSRDIQAECFLDFEVNKINIRPEYMNNPKELAKIRAMIDELKSDPSIKVNKLDIVGYASPEGSLANNKRLSEGRAMALRDYLASRYDFSRNQYYIIFGGENWDGLVKALDTIDFEYKDEALNIINDIPVEKGREAKLMQLRGGVPYRYMLKYIFPSLRVAICKVNYEIKNFNLDEAKEIIKTRPQNLSLNEMFMVANSYPKGSQEFIDVFETAVRMYPKDEIASINAAAAALSRNDLVSAERYLNMVNVNKQLPEYSNAMGVLMLLKGEYEHAEEYLKAAAKSGLQAAGQNLEELAKKKTNAAEIEKIENRDK:::::-:::::PROFILES3D,,,,,-0.512,-0.321666666666667,-0.33,-0.28,-0.241111111111111,-0.161,-0.073,0.082,0.171,0.105,0.099,-0.019,-0.032,-0.094,-0.16,-0.246,-0.269,-0.386,-0.346,-0.271,-0.436,-0.423,-0.343,-0.244,-0.141,-0.055,-0.092,0.025,-0.12,-0.074,0.091,0.177,0.168,0.107,0.126,0.053,0.178,0.169,0.235,0.202,0.00899999999999999,-0.062,-0.038,-0.146,-0.149,-0.157,-0.266,-0.223,-0.12,-0.208,-0.051,-0.07,-0.094,-0.02,-0.139,-0.16,-0.239,-0.327,-0.426,-0.492,-0.637,-0.611,-0.672,-0.756,-0.832,-0.742,-0.705,-0.71,-0.687,-0.587,-0.512,-0.448,-0.387,-0.232,-0.159,-0.196,-0.092,-0.113,-0.037,-0.071,0.008,0.00800000000000001,0.00800000000000001,-0.029,-0.102,-0.155,-0.127,-0.211,-0.298,-0.301,-0.332,-0.426,-0.505,-0.49,-0.414,-0.419,-0.488,-0.333,-0.345,-0.31,-0.31,-0.262,-0.244,-0.176,-0.097,-0.092,-0.023,-0.036,0.013,0.096,0.142,0.117,0.152,0.13,0.17,0.17,0.133,0.142,0.12,0.039,0.011,-0.021,-0.00799999999999998,-0.114,-0.176,-0.064,-0.208,-0.18,-0.388,-0.488,-0.441,-0.437,-0.4,-0.322,-0.345,-0.454,-0.254,-0.282,-0.092,0.077,0.012,0.065,0.028,0.068,0.047,0.052,0.049,0.096,0.198,0.198,0.173,0.219,0.075,-0.032,0.058,0.041,-0.065,-0.084,-0.124,-0.348,-0.376,-0.391,-0.389,-0.376,-0.471,-0.384,-0.37,-0.469,-0.54,-0.596,-0.6,-0.764,-0.566,-0.465,-0.394,-0.484,-0.561,-0.534,-0.416,-0.317,-0.228,-0.086,-0.105,-0.214,-0.19,-0.152,-0.102,-0.03,-0.228,-0.071,-0.156,-0.156,-0.228,-0.162,-0.162,-0.188,-0.131,-0.234,-0.158,-0.18,-0.081,-0.094,-0.107,-0.252,-0.276,-0.223,-0.158,-0.15,-0.071,-0.025,-0.062,-0.102,-0.088,0.066,-0.091,-0.149,-0.149,-0.035,-0.117,-0.186,-0.248,-0.253,-0.281,-0.337,-0.23,-0.172,-0.29,-0.359,-0.432,-0.432,-0.33,-0.325,-0.321,-0.356,-0.282,-0.348,-0.31,-0.297,-0.287,-0.32,-0.422,-0.356,-0.364,-0.295,-0.397,-0.473,-0.458,-0.511,-0.357,-0.279,-0.18,-0.193,-0.114,-0.18,-0.181,-0.039,-0.057,-0.062,-0.132,-0.123,-0.169,-0.222,-0.226,-0.213,-0.156,-0.143,-0.15,-0.234,-0.22,-0.374,-0.328,-0.34,-0.34,-0.464,-0.487,-0.5,-0.41,-0.228,-0.172,0.013,-0.212,-0.213,-0.185,-0.14,-0.25,-0.237,-0.313,-0.316,-0.315,-0.319,-0.165,-0.161,-0.23,-0.098,0.081,0.018,-0.028,-0.028,-0.042,-0.334,-0.275,-0.176,-0.107,-0.07,-0.07,0.017,-0.061,-0.089,-0.049,0.185,0.197,0.016,-0.0900000000000001,-0.172,-0.172,-0.196,-0.003,-0.074,-0.176,-0.176,-0.247,-0.08,-0.02,-0.044,-0.09,-0.044,-0.071,-0.062,-0.041,-0.019,-0.047,-0.061,-0.027,0.079,0.039,0.036,-0.055,-0.055,-0.12,-0.277,-0.268,-0.277,-0.246,-0.27,-0.283,-0.326,-0.27,-0.18,-0.133,0.024,-0.1,-0.153,-0.209,-0.185,-0.086,-0.062,-0.049,-0.077,0.029,-0.128,0.071,0.161,0.198,0.217,0.202,0.156,0.103,0.107,0.064,0.171,0.186,0.149,0.112,0.056,0.007,0.007,0.097,0.121,0.121,0.092,0.068,0.068,0.002,0.055,0.107,0.089,0.089,0.052,0.0355555555555556,0.11375,0.101428571428571,0.106666666666667,:::::REPIMPS,,,,,0.146,0.045,0.105714285714286,-0.07875,-0.222222222222222,-0.38,-0.682,-0.821,-0.893,-0.645,-0.687,-0.614,-0.841,-0.677,-0.429,-0.119,-0.099,-0.0439999999999999,-0.298,-0.362,-0.071,0.028,0.18,0.125,0.09,-0.22,0.00999999999999999,-0.045,0.244,0.059,-0.232,-0.386,-0.573,-0.45,-0.572,-0.297,-0.417,-0.604,-0.567,-0.567,-0.257,-0.103,0.084,0.0630000000000001,-0.106,-0.357,-0.181,-0.0119999999999999,-0.151,0.098,-0.247,-0.243,-0.43,-0.684,-0.358,-0.072,-0.128,-0.238,-0.183,-0.22,0.069,-0.0219999999999999,0.037,0.345,0.362,0.308,0.261,0.547,0.382,0.152,-0.0589999999999999,-0.0220000000000001,-0.081,-0.424,-0.478,-0.691,-0.671,-0.981,-0.871,-0.604,-0.605,-0.605,-0.605,-0.295,-0.241,0.026,-0.049,0.259,0.072,0.035,-0.00299999999999994,-0.216,0.11,-0.161,-0.178,-0.488,-0.616,-0.959,-0.717,-0.665,-0.665,-0.704,-0.971,-0.932,-1.258,-0.948,-0.716,-0.638,-0.554,-0.86,-0.608,-0.306,-0.039,0.07,0.07,0.07,-0.143,0.041,-0.098,-0.074,-0.103,-0.118,-0.366,-0.264,-0.225,-0.492,-0.225,-0.3,-0.161,0.069,-0.0239999999999999,-0.108,0.105,0.078,0.098,0.02,-0.156,-0.081,-0.062,-0.06,-0.19,-0.393,-0.606,-0.833,-0.547,-0.285,-0.454,-0.547,-0.857,-0.857,-0.555,-0.303,-0.0359999999999999,-0.036,-0.194,-0.214,0.131,0.253,0.48,0.224,0.125,0.021,0.019,0.246,0.094,-0.109,-0.377,-0.426,-0.446,-0.155,-0.071,0.0850000000000001,-0.089,-0.182,-0.217,-0.217,0.032,-0.023,-0.096,-0.0769999999999999,-0.114,-0.379,-0.257,-0.081,0.106,0.321,0.321,0.321,0.495,0.15,0.103,0.103,-0.000999999999999979,-0.346,-0.242,-0.151,-0.438,-0.299,-0.316,-0.277,-0.332,-0.084,-0.212,0.077,-0.214,-0.481,-0.351,-0.642,-0.968,-0.716,-0.825,-0.806,-0.524,-0.736,-0.391,-0.434,-0.434,-0.3,-0.011,-0.243,-0.079,-0.389,-0.488,-0.224,-0.224,-0.181,-0.108,0.049,0.103,0.162,0.182,0.492,0.285,0.233,0.179,0.427,0.495,0.247,0.193,0.114,0.094,-0.251,0.04,-0.171,-0.423,-0.406,-0.677,-0.41,-0.512,-0.551,-0.502,-0.424,-0.75,-0.502,-0.215,-0.48,-0.265,-0.575,-0.703,-0.626,-0.878,-0.611,-0.349,-0.33,-0.617,-0.695,-0.608,-0.3,-0.018,0.194,0.446,0.14,0.14,0.037,0.057,0.135,-0.037,-0.118,-0.382,-0.671,-0.619,-0.332,-0.407,-0.36,-0.0759999999999999,-0.154,-0.069,-0.0399999999999999,0.172,0.207,0.205,0.121,-0.00700000000000001,-0.218,-0.374,-0.374,-0.302,-0.302,-0.329,-0.038,-0.34,-0.291,-0.059,0.154,0.258,0.549,0.568,0.539,0.312,0.025,0.173,0.227,0.225,0.196,0.196,-0.095,-0.405,-0.425,-0.115,-0.115,0.039,-0.0170000000000001,-0.267,-0.24,-0.492,-0.414,-0.372,-0.269,-0.555,-0.594,-0.693,-0.72,-0.72,-0.718,-0.416,-0.325,-0.113,-0.113,0.169,0.514,0.617,0.433,0.563,0.272,0.173,0.004,-0.104,-0.158,-0.499,-0.844,-0.947,-0.68,-0.771,-0.48,-0.431,-0.14,-0.243,-0.272,-0.103,0.242,0.077,0.023,0.236,0.114,-0.1,-0.352,-0.085,-0.347,-0.516,-0.516,-0.302,-0.515,-0.728,-0.819,-0.642,-0.642,-0.696,-0.405,-0.405,-0.424,-0.715,-0.715,-0.467,-0.545,-0.76,-0.764,-0.868,-1.081,-0.962222222222222,-1.22125,-1.08857142857143,-1.04166666666667,:::::TMS,,,,,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,:::::

# Example of input tms file :
# SEQ_14:::::>2JQY:A:::::HELIX;;;;;:::::TMS;;;;;:::::
# SEQ_15:::::>1PHO:A:::::HELIX;;;;;106-109,198-200,:::::TMS;;;;;106-109,198-200,:::::INSIDE;;;;;112-116,:::::
# SEQ_16:::::>2MPR:A:::::HELIX;;;;;113-115,237-239,:::::TMS;;;;;113-115,:::::OUTSIDE;;;;;118-121,237-239,:::::
# SEQ_17:::::>1MAL:A:::::HELIX;;;;;113-115,160-162,398-400,:::::TMS;;;;;113-115,160-162,398-400,:::::

# If this program is being called as a cgi-script with multiple fasta input,
# then the output sequence identifier displayed will be that of the first fasta input.

# If this program is being run in command line mode, 
# then one of the following command line options and corresponding input file must be present :
#  -infile : input file must be in the special valpred format, with one line per sequence
#  -fastafile : input file must be in fasta format, with one or more sequences, 
#		each sequence taking up two or more lines, the first line being the identifier starting with a >
#  -fastachainfile : input file must be in fasta format, with one or more chains for the one sequence, 
#		each chain taking up two or more lines, the first line being the identifier starting with a >
#  -inscoresfile : input file must contains the PROFILES3D and REPIMPS scores, 1 line per sequence.
#		This input file is the valpred_scores file produced by this program.
#  -graphpointsfile : input file must be in valpred graph-points output format, with one line per sequence
#		with the valpred results already calculated and turned into streams of numbers to be plotted on a graph
#  -predictions : input file, 1 line per protein, contains previously done valpred predictions,
#		is the valpred_segments.txt output file from running this program with -fastafile, -fastachainfile or -infile
#  -benchmarkpredictions : input file, 1 line per protein, contains previously done predictions,
#		is the valpred_benchmark.txt output file from running this program with the -benchmark option, with -fastafile, -fastachainfile or -infile

# Running this program in command line mode with the -fastafile or -infile parameter, 
# to produce the valpred segments output file,
# is a handy way to run this program on a whole genome in order to find out 
# which open reading frame inputs have possible transmembrane segments.

# Running this program in command line mode with the -version parameter set to valpred1 or valpred2
# will make this program set parameters and output graph lines to the defaults for the chosen method : valpred1 or valpred2,
# before reading the -options files (if specified) for fine-tuning/changing of those defaults.

# The following data files are opened and read by this program.
# If those files are not able to be read by this program,
# then this program will use a hard-coded version of the data that is hard-coded into this program.
#  - ./Amino_Acids/ala.pdb, arg.pdb, asn.pdb, ... trp.pdb, tyr.pdb, val.pdb
#  - ./Sphere_Points/maxvol.3.252.txt
#  - ./Valpred_Files/ChartOptions.opt

# Some default options are hard-coded into this program.
# If the ./Valpred_Files/ChartOptions.opt file exists, then that file's contents override the default options.
# In command line mode, if the command line has an -options file, then that file's contents override the options.
# In web server CGI script mode, if input fields for options are received from the HTML form, then they override the options.

# Example of the program options input file, this file can be specified on the command line :
# SASARES:::::2:::::
# GRAPHRESDIUEWIDTH:::::10:::::
# GRAPHHEIGHT:::::600:::::
# VALPRED1_OR_VALPRED2:::::1:::::
# MINAREADIFF_FOR_VALPRED1:::::10:::::
# MINAREADIFF_FOR_VALPRED2:::::10:::::
# MVEAVE:::::10:::::
# MOV_AVG_1:::::10:::::
# MOV_AVG_2:::::5:::::
# MOV_AVG_3:::::5:::::
# MOV_AVG_4:::::6:::::
# MINAVESASA:::::0:::::
# TMLENGTHMIN_FOR_VALPRED1:::::12:::::
# TMLENGTHMAX_FOR_VALPRED1:::::40:::::
# TMLENGTHMIN_FOR_VALPRED2:::::6:::::
# TMLENGTHMAX_FOR_VALPRED2:::::38:::::
# AVEREPRANGEMIN:::::0.3:::::
# AVEREPRANGEMAX:::::0.8:::::
# AVEPROFRANGEMIN:::::-0.6:::::
# AVEPROFRANGEMAX:::::0:::::
# LENDIVAREAMIN:::::0.:::::
# LENDIVAREAMAX:::::1.8:::::
# MINSCOREDIFF_FOR_VALPRED1:::::0.2:::::
# MINSCOREDIFF_MOVAVG1_FOR_VALPRED2:::::0.2:::::
# MINSCOREDIFF_MOVAVG2_FOR_VALPRED2:::::0.2:::::
# MIN_AVG_AREA_MOVAVG1:::::0.2:::::
# MIN_AVG_AREA_MOVAVG2:::::0.2:::::
# MAX_NONTMS_LENGTH:::::3:::::
# MIN_NONTMS_LENGTH:::::1:::::
# MIN_NONTMS_SCORE_DIFF:::::0:::::
# MIN_NONTMS_AREA:::::1.2:::::
# TWILIGHT_AREA_PER_RESIDUE_LOWER_LIMIT:::::0.8:::::
# TWILIGHT_AREA_PER_RESIDUE_UPPER_LIMIT:::::1.3:::::
# IGNORE_TWILIGHT_AREA:::::0:::::
# SEARCH_AREA_FOR_NEIGHBOUR_TMS:::::60:::::
# HELIX_LENGTH_MIN:::::6:::::
# MIN_SCORE_DIFF_FOR_TMS_ENDS:::::0.6:::::
# PROFILES3D-LINE:::::1:::::1:::::1:::::0:::::0:::::255:::::blue:::::Profiles 3D:::::0:::::0:::::
# REPIMPS-LINE:::::1:::::1:::::1:::::255:::::0:::::0:::::red:::::REPIMPS:::::0:::::0:::::
# OVERALL:::::0:::::1:::::10:::::0:::::0:::::0:::::black:::::Overall:::::0:::::1:::::
# PROFILES3D-MOV-AVG-1:::::1:::::1:::::10:::::0:::::0:::::255:::::blue:::::Profiles 3D (mv.av.1):::::0:::::1:::::
# REPIMPS-MOV-AVG-1:::::1:::::1:::::10:::::255:::::0:::::0:::::red:::::REPIMPS (mv.av.1):::::0:::::1:::::
# TMS-MOV-AVG-1:::::1:::::1:::::0:::::0:::::0:::::0:::::red:::::TMS (mv.av.1):::::1.6:::::0:::::
# PROFILES3D-MOV-AVG-2:::::1:::::1:::::5:::::0:::::0:::::255:::::lblue:::::Profiles 3D (mv.av.2):::::0:::::0:::::
# REPIMPS-MOV-AVG-2:::::1:::::1:::::5:::::255:::::0:::::0:::::lred:::::REPIMPS (mv.av.2):::::0:::::0:::::
# TMS-MOV-AVG-2:::::1:::::1:::::0:::::0:::::0:::::0:::::lred:::::TMS (mv.av.2):::::1.58:::::0:::::
# PROFILES3D-MOV-AVG-3:::::0:::::1:::::5:::::0:::::0:::::255:::::lblue:::::Profiles 3D (mv.av.3):::::0:::::0:::::
# REPIMPS-MOV-AVG-3:::::0:::::1:::::5:::::255:::::0:::::0:::::lred:::::REPIMPS (mv.av.3):::::0:::::0:::::
# HYDROPHILIC-MOV-AVG-3:::::1:::::1:::::0:::::0:::::0:::::0:::::lblue:::::Hydrophilic area (mv.av.3):::::1.56:::::0:::::
# PROFILES3D-MOV-AVG-4:::::1:::::1:::::10:::::0:::::0:::::255:::::cyan:::::Profiles 3D (mv.av.4):::::0:::::1:::::
# REPIMPS-MOV-AVG-4:::::1:::::1:::::10:::::255:::::0:::::0:::::pink:::::REPIMPS (mv.av.4):::::0:::::1:::::
# HYDROPHILIC-MOV-AVG-4:::::1:::::1:::::0:::::0:::::0:::::0:::::cyan:::::Hydrophilic area (mv.av.4):::::1.54:::::0:::::
# TMS:::::1:::::1:::::0:::::0:::::0:::::0:::::black:::::Transmembrane Segments:::::1.7:::::0:::::
# TRUE-TMS:::::0:::::1:::::10:::::0:::::0:::::0:::::green:::::true TMS:::::1.9:::::0:::::
# TRUE-INSIDE:::::0:::::1:::::10:::::0:::::0:::::0:::::lgreen:::::true inside:::::1.86:::::0:::::
# TRUE-OUTSIDE:::::0:::::1:::::10:::::0:::::0:::::0:::::dgreen:::::true outside:::::1.84:::::0:::::

# The default SASA resolution (SASARes) is 2.
# Here are the permitted SASA resolution values and the files used for each resolution :
# 0 ./maxvol.3.92.txt
# 1 ./maxvol.3.162.txt
# 2 ./maxvol.3.252.txt
# 3 ./maxvol.3.392.txt
# 4 ./maxvol.3.482.txt
# 5 ./maxvol.3.1002.txt

# Example of the *.valpred_benchmark.txt output file of TMS helices 
# for the Rost TMH benchmark server at http://cubic.bioc.columbia.edu/services/tmh_benchmark/
# >0
# 18,35
# 43,59
# >1
# >2
# 66,90
# >3
# 22,39
# 43,62
# 69,88

# Example of the *.valpred_evaluate.txt output file
# SEQ_1:::::>1BL8:A:::::GGGGGGGGGGGGGGGGGGGGGGALHWRAAGAATVLLVIVLLAGSYLAVLAERGAPGAQLITYPRALWWSVETATTVGYGDLYPVTLWGRCVAVVVMVAGITSFGLVTAALATWFVGREQ:::::27-51,62-74,86-112,:::::23-43,72-111,:::::toplogy_correct=0,num_tms_observed=3,num_tms_observed_correctly_predicted=2,num_tms_predicted=2,num_tms_predicted_correctly_predicted=2,num_residues=119,num_observed_tms_residues=65,num_observed_non_tms_residues=54,num_predicted_tms_residues=61,num_predicted_non_tms_residues=58,num_correctly_predicted_residues=85,num_correctly_predicted_tms_residues=46,num_correctly_predicted_non_tms_residues=39,:::::
# SEQ_2:::::>1AP9:::::GGGGGGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLR:::::10-30,37-61,83-93,105-126,134-152,166-169,172-189,202-220,:::::8-29,38-67,77-155,166-190,197-224,:::::toplogy_correct=0,num_tms_observed=8,num_tms_observed_correctly_predicted=5,num_tms_predicted=5,num_tms_predicted_correctly_predicted=5,num_residues=225,num_observed_tms_residues=139,num_observed_non_tms_residues=86,num_predicted_tms_residues=184,num_predicted_non_tms_residues=41,num_correctly_predicted_residues=176,num_correctly_predicted_tms_residues=137,num_correctly_predicted_non_tms_residues=39,:::::
# SEQ_3:::::>1AR1:A:::::GGGGGGGGGGGGGGGGGFFTRWFMSTNHKDIGILYLFTAGIVGLISVCFTVYMRMELQHPGVQYMCLEGARLIADASAECTPNGHLWNVMITYHGVLMMFFVVIPALFGGFGNYFMPLHIGAPDMAFPRLNNLSYWMYVCGVALGVASLLAPGGNDQMGSGVGWVLYPPLSTTEAGYSMDLAIFAVHVSGASSILGAINIITTFLNMRAPGMTLFKVPLFAWSVFITAWLILLSLPVLAGAITMLLMDRNFGTQFFDPAGGGDPVLYQHILWFFGHPEVYIIILPGFGIISHVISTFAKKPIFGYLPMVLAMAAIGILGFVVWAHHMYTAGMSLTQQAYFMLATMTIAVPTGIKVFSWIATMWGGSIEFKTPMLWAFGFLFLFTVGGVTGVVLSQAPLDRVYHDTYYVVAHFHYVMSLGAVFGIFAGVYYWIGKMSGRQYPEWAGQLHFWMMFIGSNLIFFPQHFLGRQGMPRRYIDYPVEFAYWNNISSIGAYISFASFLFFIGIVFYTLFAGKRVNVPNYWNEHADTLEWTLPSPPPEHTFET:::::28-58,84-101,103-120,128-150,178-205,219-250,257-260,264-298,305-218,323-326,328-330,334-363,371-394,396-402,406-416,420-436,442-460,462-466,480-513,530-532,:::::25-51,82-117,129-146,174-245,263-297,299-328,331-363,365-394,397-430,442-461,474-511,:::::toplogy_correct=0,num_tms_observed=20,num_tms_observed_correctly_predicted=11,num_tms_predicted=11,num_tms_predicted_correctly_predicted=11,num_residues=545,num_observed_tms_residues=346,num_observed_non_tms_residues=199,num_predicted_tms_residues=373,num_predicted_non_tms_residues=172,num_correctly_predicted_residues=430,num_correctly_predicted_tms_residues=302,num_correctly_predicted_non_tms_residues=128,:::::



use warnings;
use strict;
use diagnostics;
use Getopt::Long;
use Math::Complex;
use Math::Trig;
use GD;
use GD::Text;
use GD::Graph::lines;
use CGI;
use Data::Dumper;
$Data::Dumper::Purity = 1;
# print "var x = " . Dumper( $x ) . "\n"; # debug
# my $wait_for_key_press = <STDIN>; # debug

# global constants

				# in Residue.h, enum SecStruct{
use constant SHEET => 0;	# SHEET = 0,
use constant HELIX => 1;	# HELIX = 1,	
use constant OTHER => 2;	# OTHER = 2,
				# in Residue.h, enum ResCode{
use constant ALA => 7;		# ALA = 7,
use constant ARG => 19;		# ARG = 19,
use constant ASN => 14;		# ASN = 14,
use constant ASP => 16;		# ASP = 16,
use constant CYS => 10;		# CYS = 10,
use constant GLN => 13;		# GLN = 13,
use constant GLU => 15;		# GLU = 15,
use constant GLY => 8;		# GLY = 8,
use constant HIS => 17;		# HIS = 17,
use constant ILE => 4;		# ILE = 4,
use constant LEU => 3;		# LEU = 3,
use constant LYS => 18;		# LYS = 18,
use constant MET => 6;		# MET = 6,
use constant PHE => 1;		# PHE = 1,
use constant PRO => 9;		# PRO = 9,
use constant SER => 12;		# SER = 12,
use constant THR => 11;		# THR = 11,
use constant TRP => 0;		# TRP = 0,
use constant TYR => 2;		# TYR = 2,
use constant VAL => 5;		# VAL = 5,
use constant UNK => 21;		# UNK = 21
use constant PI => 3.1415926535897932384626433832795; 	# define PI 3.1415926535897932384626433832795f;
							# in Atom.h, enum AtomType{
use constant C_Ali => 0;				# C_Ali = 0,		// Aliphatic C
use constant C_Car => 1;				# C_Car = 1,		// Carboxyl or Carbonyl
use constant C_Aro => 2;				# C_Aro = 2,		// Aromatic C
use constant N_All => 3;				# N_All = 3,		// All Nitrogen
use constant O_All => 4;				# O_All = 4,		// All Oxygen
use constant S_Oxy => 5;				# S_Oxy = 5,		// Sulfur and thiol S of S-S bridge
use constant S_Red => 6;				# S_Red = 6		// SH or S-CH3 etc
my %charToRes;
set_charToRes();
my %assignResType;
set_assignResType();
my %Residue_SecStruct;
my %Residue_ResCode;
my %Residue_EnvCode;
set_Residue_constants();
my @SCGlyXGly_array;
my $SCGlyXGly = \@SCGlyXGly_array;
set_Value_constants();
my $atom_default_values = Atom_default_values();
my %amino_acid_pdbs_hash;
my $amino_acid_pdbs = \%amino_acid_pdbs_hash;
my @tmDomains_array;
my $tmDomains = \@tmDomains_array; # vector<Range> tmDomains;
my @tmDomain_start;
my @amphipathic_helices_array;
my $amphipathic_helices = \@amphipathic_helices_array;
my @amphipathic_helices_chains_array;
my $amphipathic_helices_chains = \@amphipathic_helices_chains_array;
my $path_to_files_Sphere_Points = 'Sphere_Points';
my $path_to_files_Amino_Acids = 'Amino_Acids';
my $path_to_files_Valpred_Files = 'Valpred_Files';
my $program_options;
my $TMPred2D;

# global variables

my $analyser; #C++ Analyser
my $structure; #C++ Structure
my $helices; #C++ Range
my $multiplegraphs_data;
my $multiplegraphs_filename;
my $input_fasta_id = '';
my $input_sequence_id = '';
my $input_aa_sequence = '';
my $input_infile;
my $input_inscoresfile;
my $input_flag_multiplegraphs;
my $input_flag_benchmark;
my $input_flag_version;
my $input_flag_amphipathic;
my $input_flag_help;
my $input_flag_h;
my $input_debug_flag;
my $input_fasta_chain_file;
my $input_fasta_file;
my $input_predictions_file;
my $input_benchmarkpredictions_file;
my $input_graph_points_file;
my $input_options_file;
my $input_file_name_for_output_files = 'valpred_perl_program';
my $input_tms_file;
my $input_tms_line = '';
my $input_inside_line = '';
my $input_outside_line = '';
my $input_TMS_line = '';
my $input_TMS_HELIX_line = '';
my $input_HELIX_line = '';
my $input_HELIX_INSIDE_line = '';
my $input_HELIX_OUTSIDE_line = '';
my $input_HELIX_MEMBRANE_line = '';
my @input_tms_array_fasta_id;
my @input_tms_array_tms_line;
my @input_tms_array_inside_line;
my @input_tms_array_outside_line;
my @input_tms_array_TMS_line;
my @input_tms_array_TMS_HELIX_line;
my @input_tms_array_HELIX_line;
my @input_tms_array_HELIX_INSIDE_line;
my @input_tms_array_HELIX_OUTSIDE_line;
my @input_tms_array_HELIX_MEMBRANE_line;
my $segment_output_file;
my $scores_output_file;
my $amphipathic_output_file;
my $benchmark_output_file;
my $html_graphs_output_file;
my $evaluate_output_file;
my $stats;
my %stats_hash;
$stats = \%stats_hash;



#===============================================================================
#==========               MAIN PROGRAM LOGIC           =========================
#===============================================================================

my $program_mode = '';
my $q = new CGI;

init();

program_get_input();

if ($program_mode eq 'cgi-display-form') {

	program_display_html_form();

} elsif ($program_mode eq 'cgi-output-graph') {

	program_output_html_graph();

} elsif ($program_mode eq 'cgi-output') {

	program_output_html_results();

} else { # $program_mode eq 'command-line'

	program_output_command_line();
}

#===============================================================================

sub adjust_options {

	# If the VALPRED1 option has be chosen, 
	# then don't try to use or print scores and lines from the VALPRED2 option.

	if ($TMPred2D->{'valpred1_or_valpred2'} == 1) {
		for ( my $ln = 0; $ln <= $program_options->{'option_max_index'}; $ln++ ) {
			if (	($ln == $program_options->{'lnOpt_index'}->{'PROFILES3D-LINE'}) ||
				($ln == $program_options->{'lnOpt_index'}->{'REPIMPS-LINE'}) ||
				($ln == $program_options->{'lnOpt_index'}->{'TMS-MOV-AVG-1'}) ||
				($ln == $program_options->{'lnOpt_index'}->{'PROFILES3D-MOV-AVG-2'}) ||
				($ln == $program_options->{'lnOpt_index'}->{'REPIMPS-MOV-AVG-2'}) ||
				($ln == $program_options->{'lnOpt_index'}->{'TMS-MOV-AVG-2'}) ||
				($ln == $program_options->{'lnOpt_index'}->{'PROFILES3D-MOV-AVG-3'}) ||
				($ln == $program_options->{'lnOpt_index'}->{'REPIMPS-MOV-AVG-3'}) ||
				($ln == $program_options->{'lnOpt_index'}->{'HYDROPHILIC-MOV-AVG-3'}) ||
				($ln == $program_options->{'lnOpt_index'}->{'PROFILES3D-MOV-AVG-4'}) ||
				($ln == $program_options->{'lnOpt_index'}->{'REPIMPS-MOV-AVG-4'}) ||
				($ln == $program_options->{'lnOpt_index'}->{'HYDROPHILIC-MOV-AVG-4'}) ) {

				$program_options->{'lnOpt'}->[$ln]->{'showGraph'} = 0;
			}
		}
	}

	# If the VALPRED2 option has been chosen, 
	# then don't try to use or print scores and lines from the VALPRED2 option.

	if ($TMPred2D->{'valpred1_or_valpred2'} == 2) {
		for ( my $ln = 0; $ln <= $program_options->{'option_max_index'}; $ln++ ) {
			if ($ln == $program_options->{'lnOpt_index'}->{'BOTH-LINE'}) {

				$program_options->{'lnOpt'}->[$ln]->{'showGraph'} = 0;
			}
		}
	}

	# If this program is being run as a cgi-script rather than from the command line,
	# then there will not be a file containing the true-tms/true-inside/true-outside residues,
	# so don't try to display those lines.

	if ($program_mode eq 'cgi-output') {
		for ( my $ln = 0; $ln <= $program_options->{'option_max_index'}; $ln++ ) {
			if (	($ln == $program_options->{'lnOpt_index'}->{'TRUE-TMS'}) ||
				($ln == $program_options->{'lnOpt_index'}->{'TRUE-INSIDE'}) ||
				($ln == $program_options->{'lnOpt_index'}->{'TRUE-OUTSIDE'}) ) {

				$program_options->{'lnOpt'}->[$ln]->{'showGraph'} = 0;
			}
		}
	}

	# None of the TMS lines should be smoothed. They are straight lines.
	# They can also be discontinuous lines, which would cause the smoothing algorithm to crash.

	for ( my $ln = 0; $ln <= $program_options->{'option_max_index'}; $ln++ ) {
		if (	($ln == $program_options->{'lnOpt_index'}->{'TMS'}) ||
			($ln == $program_options->{'lnOpt_index'}->{'TRUE-TMS'}) ||
			($ln == $program_options->{'lnOpt_index'}->{'TRUE-INSIDE'}) ||
			($ln == $program_options->{'lnOpt_index'}->{'TRUE-OUTSIDE'}) ||
			($ln == $program_options->{'lnOpt_index'}->{'TMS-MOV-AVG-1'}) ||
			($ln == $program_options->{'lnOpt_index'}->{'TMS-MOV-AVG-2'}) ||
			($ln == $program_options->{'lnOpt_index'}->{'HYDROPHILIC-MOV-AVG-3'}) ||
			($ln == $program_options->{'lnOpt_index'}->{'HYDROPHILIC-MOV-AVG-4'}) ) {

			$program_options->{'lnOpt'}->[$ln]->{'smooth'} = 0;
		}
	}

	# If the VALPRED1 option has be chosen and the input is scores, not residues,
	# then the side ASA (side chain accessible surface area) is not calculated
	# and is not available for the VALPRED1 algorithm.
	if ($TMPred2D->{'valpred1_or_valpred2'} == 1) {
		if (defined ($input_inscoresfile)) {
			$TMPred2D->{'ASA_are_available'} = 0;
			print "\nResides are not available to calculate side-ASA (side chain accessible surface area) for VALPRED1 algorithm.\n";
			print "Thus, side-ASA restrictions will be ignored in the VALPRED1 algorithm.\n\n";
		}
	}
}



# // Allocates environment categories and scores to residues
sub Analyser_doProfile { #C++ void Analyser::doProfile(void)

	my @show_aacode = ("TRP","PHE","TYR","LEU","ILE","VAL","MET","ALA","GLY","PRO","CYS","THR","SER","GLN","ASN","GLU","ASP","HIS","LYS","ARG","  ","UNK");
	my @show_envcode = ("B1a","B1b","B1","B2a","B2b","B2","B3a","B3b","B3","P1a","P1b","P1","P2a","P2b","P2","Ea","Eb","E","unknown");

	$structure = $analyser->{'st'};
	my $st_dot_allResidues_dot_size = $#{$structure->{'allResidues'}} + 1;
	for (my $i = 0; $i < $st_dot_allResidues_dot_size; $i++ ) {					# for(unsigned int i = 0;i < st.allResidues.size(); i++){
		# //soluble analysis
		if ( $structure->{'allResidues'}->[$i]->{'areaBuried'} < $analyser->{'ABLow'} ) {	# if(st.allResidues[i].areaBuried < ABLow){
			if ( $structure->{'allResidues'}->[$i]->{'SS'} == HELIX ) {			# if(st.allResidues[i].SS == HELIX){
				$structure->{'allResidues'}->[$i]->{'profEnv'} = $Residue_EnvCode{'Ea'};# st.allResidues[i].profEnv = Ea;
			} elsif ( $structure->{'allResidues'}->[$i]->{'SS'} == SHEET) {			# else if(st.allResidues[i].SS == SHEET){
				$structure->{'allResidues'}->[$i]->{'profEnv'} = $Residue_EnvCode{'Eb'};# st.allResidues[i].profEnv = Eb;
			} else {
				$structure->{'allResidues'}->[$i]->{'profEnv'} = $Residue_EnvCode{'E'};	# st.allResidues[i].profEnv = E;
			}
		} elsif ( $structure->{'allResidues'}->[$i]->{'areaBuried'} < $analyser->{'ABHigh'} ) {	# else if(st.allResidues[i].areaBuried < ABHigh){

			if ( $structure->{'allResidues'}->[$i]->{'fractionPolar'} < $analyser->{'FPPartial'} ) { # if(st.allResidues[i].fractionPolar <FPPartial){
				if ( $structure->{'allResidues'}->[$i]->{'SS'} == HELIX ) {		# if(st.allResidues[i].SS == HELIX){
					$structure->{'allResidues'}->[$i]->{'profEnv'} = $Residue_EnvCode{'P1a'}; # st.allResidues[i].profEnv = P1a;
				} elsif ( $structure->{'allResidues'}->[$i]->{'SS'} == SHEET ) {	# else if(st.allResidues[i].SS == SHEET){
					$structure->{'allResidues'}->[$i]->{'profEnv'} = $Residue_EnvCode{'P1b'}; # st.allResidues[i].profEnv = P1b;
				} else {
					$structure->{'allResidues'}->[$i]->{'profEnv'} = $Residue_EnvCode{'P1'}; # st.allResidues[i].profEnv = P1;
				}
			}

			if ( $structure->{'allResidues'}->[$i]->{'fractionPolar'} >= $analyser->{'FPPartial'} ) { # if(st.allResidues[i].fractionPolar >=FPPartial){
				if ( $structure->{'allResidues'}->[$i]->{'SS'} == HELIX ) {		# if(st.allResidues[i].SS == HELIX){
					$structure->{'allResidues'}->[$i]->{'profEnv'} = $Residue_EnvCode{'P2a'}; # st.allResidues[i].profEnv = P2a;
				} elsif ( $structure->{'allResidues'}->[$i]->{'SS'} == SHEET ) {	# else if(st.allResidues[i].SS == SHEET){
					$structure->{'allResidues'}->[$i]->{'profEnv'} = $Residue_EnvCode{'P2b'}; # st.allResidues[i].profEnv = P2b;
				} else {
					$structure->{'allResidues'}->[$i]->{'profEnv'} = $Residue_EnvCode{'P2'}; # st.allResidues[i].profEnv = P2;
				}
			}

		} else {
			if ( $structure->{'allResidues'}->[$i]->{'fractionPolar'} < $analyser->{'FPBLow'} ) { # if(st.allResidues[i].fractionPolar <FPBLow){
				if ( $structure->{'allResidues'}->[$i]->{'SS'} == HELIX ) {		# if(st.allResidues[i].SS == HELIX){
					$structure->{'allResidues'}->[$i]->{'profEnv'} = $Residue_EnvCode{'B1a'}; # st.allResidues[i].profEnv = B1a;
				} elsif ( $structure->{'allResidues'}->[$i]->{'SS'} == SHEET ) {	# else if(st.allResidues[i].SS == SHEET){
					$structure->{'allResidues'}->[$i]->{'profEnv'} = $Residue_EnvCode{'B1b'}; # st.allResidues[i].profEnv = B1b;
				} else { 
					$structure->{'allResidues'}->[$i]->{'profEnv'} = $Residue_EnvCode{'B1'}; # st.allResidues[i].profEnv = B1;
				}
			} elsif ( $structure->{'allResidues'}->[$i]->{'fractionPolar'} < $analyser->{'FPBHigh'} ) { # else if(st.allResidues[i].fractionPolar <FPBHigh){
				if ( $structure->{'allResidues'}->[$i]->{'SS'} == HELIX ) {		# if(st.allResidues[i].SS == HELIX){
					$structure->{'allResidues'}->[$i]->{'profEnv'} = $Residue_EnvCode{'B2a'}; # st.allResidues[i].profEnv = B2a;
				} elsif ( $structure->{'allResidues'}->[$i]->{'SS'} == SHEET ) {	# else if(st.allResidues[i].SS == SHEET){
					$structure->{'allResidues'}->[$i]->{'profEnv'} = $Residue_EnvCode{'B2b'}; # st.allResidues[i].profEnv = B2b;
				} else {
					$structure->{'allResidues'}->[$i]->{'profEnv'} = $Residue_EnvCode{'B2'}; # st.allResidues[i].profEnv = B2;
				}
			} else {
				if ( $structure->{'allResidues'}->[$i]->{'SS'} == HELIX ) {		# if(st.allResidues[i].SS == HELIX){
					$structure->{'allResidues'}->[$i]->{'profEnv'} = $Residue_EnvCode{'B3a'}; # st.allResidues[i].profEnv = B3a;
				} elsif ( $structure->{'allResidues'}->[$i]->{'SS'} == SHEET ) {	# else if(st.allResidues[i].SS == SHEET){
					$structure->{'allResidues'}->[$i]->{'profEnv'} = $Residue_EnvCode{'B3b'}; # st.allResidues[i].profEnv = B3b;
				} else {
					$structure->{'allResidues'}->[$i]->{'profEnv'} = $Residue_EnvCode{'B3'}; # st.allResidues[i].profEnv = B3;
				}
			}
		}

		# repimps analysis
		my $temp1 = $structure->{'allResidues'}->[$i]->{'aaCode'};
		if ( $SCGlyXGly->[$temp1] < $analyser->{'ABLow'} ) {					# if(SCGlyXGly[st.allResidues[i].aaCode] < ABLow){
			if ( $structure->{'allResidues'}->[$i]->{'SS'} == HELIX ) {			# if(st.allResidues[i].SS == HELIX){
				$structure->{'allResidues'}->[$i]->{'repEnv'} = $Residue_EnvCode{'Ea'};	# st.allResidues[i].repEnv = Ea;
			} elsif ( $structure->{'allResidues'}->[$i]->{'SS'} == SHEET ) {		# else if(st.allResidues[i].SS == SHEET){
				$structure->{'allResidues'}->[$i]->{'repEnv'} = $Residue_EnvCode{'Eb'};	# st.allResidues[i].repEnv = Eb;
			} else {
				$structure->{'allResidues'}->[$i]->{'repEnv'} = $Residue_EnvCode{'E'};	# st.allResidues[i].repEnv = E;
			}

		} elsif ( $SCGlyXGly->[$temp1] < $analyser->{'ABHigh'} ) {				# else if(SCGlyXGly[st.allResidues[i].aaCode] < ABHigh){

			if ( $structure->{'allResidues'}->[$i]->{'repimpsFP'} < 0.67 ) {		# if(st.allResidues[i].repimpsFP <0.67){
				if ($structure->{'allResidues'}->[$i]->{'SS'} == HELIX ) {		# if(st.allResidues[i].SS == HELIX){
					$structure->{'allResidues'}->[$i]->{'repEnv'} = $Residue_EnvCode{'P1a'}; # st.allResidues[i].repEnv = P1a;
				} elsif ( $structure->{'allResidues'}->[$i]->{'SS'} == SHEET ) {	# else if(st.allResidues[i].SS == SHEET){
					$structure->{'allResidues'}->[$i]->{'repEnv'} = $Residue_EnvCode{'P1b'}; # st.allResidues[i].repEnv = P1b;
				} else {
					$structure->{'allResidues'}->[$i]->{'repEnv'} = $Residue_EnvCode{'P1'}; # st.allResidues[i].repEnv = P1;
				}
			}

			if ( $structure->{'allResidues'}->[$i]->{'repimpsFP'} >= $analyser->{'FPPartial'} ) { # if(st.allResidues[i].repimpsFP >=FPPartial){
				if ( $structure->{'allResidues'}->[$i]->{'SS'} == HELIX ) {		# if(st.allResidues[i].SS == HELIX){
					$structure->{'allResidues'}->[$i]->{'repEnv'} = $Residue_EnvCode{'P2a'}; # st.allResidues[i].repEnv = P2a;
				} elsif ( $structure->{'allResidues'}->[$i]->{'SS'} == SHEET ) {	# else if(st.allResidues[i].SS == SHEET){
					$structure->{'allResidues'}->[$i]->{'repEnv'} = $Residue_EnvCode{'P2b'}; # st.allResidues[i].repEnv = P2b;
				} else {
					$structure->{'allResidues'}->[$i]->{'repEnv'} = $Residue_EnvCode{'P2'}; # st.allResidues[i].repEnv = P2;
				}
			}

		} else {
			if ( $structure->{'allResidues'}->[$i]->{'repimpsFP'} < $analyser->{'FPBLow'} ) { # if(st.allResidues[i].repimpsFP <FPBLow){
				if ( $structure->{'allResidues'}->[$i]->{'SS'} == HELIX ) {		# if(st.allResidues[i].SS == HELIX){
					$structure->{'allResidues'}->[$i]->{'repEnv'} = $Residue_EnvCode{'B1a'}; # st.allResidues[i].repEnv = B1a;
				} elsif ( $structure->{'allResidues'}->[$i]->{'SS'} == SHEET ) {	# else if(st.allResidues[i].SS == SHEET){
					$structure->{'allResidues'}->[$i]->{'repEnv'} = $Residue_EnvCode{'B1b'}; # st.allResidues[i].repEnv = B1b;
				} else {
					$structure->{'allResidues'}->[$i]->{'repEnv'} = $Residue_EnvCode{'B1'}; # st.allResidues[i].repEnv = B1;
				}
			} elsif ( $structure->{'allResidues'}->[$i]->{'repimpsFP'} < $analyser->{'FPBHigh'} ) { # else if(st.allResidues[i].repimpsFP <FPBHigh){
				if ( $structure->{'allResidues'}->[$i]->{'SS'} == HELIX ) {		# if(st.allResidues[i].SS == HELIX){
					$structure->{'allResidues'}->[$i]->{'repEnv'} = $Residue_EnvCode{'B2a'}; # st.allResidues[i].repEnv = B2a;
				} elsif ( $structure->{'allResidues'}->[$i]->{'SS'} == SHEET ) {	# else if(st.allResidues[i].SS == SHEET){
					$structure->{'allResidues'}->[$i]->{'repEnv'} = $Residue_EnvCode{'B2b'}; # st.allResidues[i].repEnv = B2b;
				} else {
					$structure->{'allResidues'}->[$i]->{'repEnv'} = $Residue_EnvCode{'B2'}; # st.allResidues[i].repEnv = B2;
				}
			} else {
				if ( $structure->{'allResidues'}->[$i]->{'SS'} == HELIX ) {		# if(st.allResidues[i].SS == HELIX){
					$structure->{'allResidues'}->[$i]->{'repEnv'} = $Residue_EnvCode{'B3a'}; # st.allResidues[i].repEnv = B3a;
				} elsif ( $structure->{'allResidues'}->[$i]->{'SS'} == SHEET ) {	# else if(st.allResidues[i].SS == SHEET){
					$structure->{'allResidues'}->[$i]->{'repEnv'} = $Residue_EnvCode{'B3b'}; # st.allResidues[i].repEnv = B3b;
				} else {
					$structure->{'allResidues'}->[$i]->{'repEnv'} = $Residue_EnvCode{'B3'}; # st.allResidues[i].repEnv = B3;
				}
			}
		}
		if ( $structure->{'allResidues'}->[$i]->{'aaCode'} == UNK ) {				# if (st.allResidues[i].aaCode == UNK){
			$structure->{'allResidues'}->[$i]->{'scores'}->[0] = 0.0;			# st.allResidues[i].scores[0] = 0.0f;
			$structure->{'allResidues'}->[$i]->{'scores'}->[1] = 0.0;			# st.allResidues[i].scores[1] = 0.0f;
			$structure->{'allResidues'}->[$i]->{'repEnv'} = $Residue_EnvCode{'unknown'};		# st.allResidues[i].repEnv = unknown;
		} else {
			# st.allResidues[i].scores[0] = scores[st.allResidues[i].profEnv][st.allResidues[i].aaCode];
			my $temp1 = $structure->{'allResidues'}->[$i]->{'profEnv'};
			my $temp2 = $structure->{'allResidues'}->[$i]->{'aaCode'};
			$structure->{'allResidues'}->[$i]->{'scores'}->[0] = $analyser->{'scores'}->[$temp1][$temp2];
			# st.allResidues[i].scores[1] = scores[st.allResidues[i].repEnv][st.allResidues[i].aaCode];
			my $temp3 = $structure->{'allResidues'}->[$i]->{'repEnv'};
			my $temp4 = $structure->{'allResidues'}->[$i]->{'aaCode'};
			$structure->{'allResidues'}->[$i]->{'scores'}->[1] = $analyser->{'scores'}->[$temp3][$temp4];
		}
		# //change when writing overall score handler!!
		# //st.aveResScore += st.allResidues[i].scores[0];

	}
	# //st.aveResScore /= st.allResidues.size();
}



sub Analyser_findCloseAtoms { #C++ void Analyser::findCloseAtoms()

	# ///recalculate xyzextremities to include hetatoms
	$structure = $analyser->{'st'};
	my $st_dot_hetAtoms_dot_size = $#{$structure->{'hetAtoms'}} + 1;
	if ($analyser->{'ignoreHet'} == 0) { # if(!ignoreHet){
		for (my $h = 0; $h < $st_dot_hetAtoms_dot_size; $h++) { # for(unsigned int h = 0; h < st.hetAtoms.size(); h++){
			if ($structure->{'hetAtoms'}->[$h]->{'coords'}->{'x'} < $structure->{'xyzExtremities'}->[0]->{'f1'}) { # if(st.hetAtoms[h].coords.x < st.xyzExtremities[0].f1){
				$structure->{'xyzExtremities'}->[0]->{'f1'} = $structure->{'hetAtoms'}->[$h]->{'coords'}->{'x'}; # st.xyzExtremities[0].f1 = st.hetAtoms[h].coords.x;
			}
			if ($structure->{'hetAtoms'}->[$h]->{'coords'}->{'x'} > $structure->{'xyzExtremities'}->[0]->{'f2'}) { # if(st.hetAtoms[h].coords.x > st.xyzExtremities[0].f2){
				$structure->{'xyzExtremities'}->[0]->{'f2'} = $structure->{'hetAtoms'}->[$h]->{'coords'}->{'x'}; # st.xyzExtremities[0].f2 = st.hetAtoms[h].coords.x;
			}
			if ($structure->{'hetAtoms'}->[$h]->{'coords'}->{'y'} < $structure->{'xyzExtremities'}->[1]->{'f1'}) { # if(st.hetAtoms[h].coords.y < st.xyzExtremities[1].f1){
				$structure->{'xyzExtremities'}->[1]->{'f1'} = $structure->{'hetAtoms'}->[$h]->{'coords'}->{'y'}; # st.xyzExtremities[1].f1 = st.hetAtoms[h].coords.y;
			}
			if ($structure->{'hetAtoms'}->[$h]->{'coords'}->{'y'} > $structure->{'xyzExtremities'}->[1]->{'f2'}) { # if(st.hetAtoms[h].coords.y > st.xyzExtremities[1].f2){
				$structure->{'xyzExtremities'}->[1]->{'f2'} = $structure->{'hetAtoms'}->[$h]->{'coords'}->{'y'}; # st.xyzExtremities[1].f2 = st.hetAtoms[h].coords.y;
			}
			if ($structure->{'hetAtoms'}->[$h]->{'coords'}->{'z'} < $structure->{'xyzExtremities'}->[2]->{'f1'}) { # if(st.hetAtoms[h].coords.z < st.xyzExtremities[2].f1){
				$structure->{'xyzExtremities'}->[2]->{'f1'} = $structure->{'hetAtoms'}->[$h]->{'coords'}->{'z'}; # st.xyzExtremities[2].f1 = st.hetAtoms[h].coords.z;
			}
			if ($structure->{'hetAtoms'}->[$h]->{'coords'}->{'z'} > $structure->{'xyzExtremities'}->[2]->{'f2'}) { # if(st.hetAtoms[h].coords.z > st.xyzExtremities[2].f2){
				$structure->{'xyzExtremities'}->[2]->{'f2'} = $structure->{'hetAtoms'}->[$h]->{'coords'}->{'z'}; # st.xyzExtremities[2].f2 = st.hetAtoms[h].coords.z;
			}
		}
	}

	# //cubing the atoms for faster processing//
	# noOfXCubes = (int)((st.xyzExtremities[0].f2 - st.xyzExtremities[0].f1) /cubeSize) + 1;
	$analyser->{'noOfXCubes'} = int (($structure->{'xyzExtremities'}->[0]->{'f2'} - $structure->{'xyzExtremities'}->[0]->{'f1'}) / $analyser->{'cubeSize'}) + 1;
	# noOfYCubes = (int)((st.xyzExtremities[1].f2 - st.xyzExtremities[1].f1) / cubeSize) + 1;
	$analyser->{'noOfYCubes'} = int (($structure->{'xyzExtremities'}->[1]->{'f2'} - $structure->{'xyzExtremities'}->[1]->{'f1'}) / $analyser->{'cubeSize'}) + 1;
	# noOfZCubes = (int)((st.xyzExtremities[2].f2 - st.xyzExtremities[2].f1) / cubeSize) + 1;
	$analyser->{'noOfZCubes'} = int (($structure->{'xyzExtremities'}->[2]->{'f2'} - $structure->{'xyzExtremities'}->[2]->{'f1'}) / $analyser->{'cubeSize'}) + 1;

	my $cubes; # vector<Atom*> *cubes = new vector<Atom*>[noOfXCubes * noOfYCubes * noOfZCubes];
	my $XYZ = $analyser->{'noOfXCubes'} * $analyser->{'noOfYCubes'} * $analyser->{'noOfZCubes'};
	for (my $init_index = 0; $init_index < $XYZ; $init_index++) {
		my @atom_vector;
		$cubes->[$init_index] = \@atom_vector;
	}

	if ($analyser->{'ignoreHet'} == 0) { # if(!ignoreHet){
		$st_dot_hetAtoms_dot_size = $#{$structure->{'hetAtoms'}} + 1;
		for (my $e = 0; $e < $st_dot_hetAtoms_dot_size; $e++) { # for (unsigned int e = 0 ; e < st.hetAtoms.size(); e++){
			# cubes[xyzToIndex(	(int)((st.hetAtoms[e].coords.x  - st.xyzExtremities[0].f1)/ cubeSize),
			#					(int)((st.hetAtoms[e].coords.y  - st.xyzExtremities[1].f1)/ cubeSize),
			#					(int)((st.hetAtoms[e].coords.z  - st.xyzExtremities[2].f1)/ cubeSize) )].push_back(&(st.hetAtoms[e]));
			my $v1 = int (($structure->{'hetAtoms'}->[$e]->{'coords'}->{'x'} - $structure->{'xyzExtremities'}->[0]->{'f1'}) / $analyser->{'cubeSize'});
			my $v2 = int (($structure->{'hetAtoms'}->[$e]->{'coords'}->{'y'} - $structure->{'xyzExtremities'}->[1]->{'f1'}) / $analyser->{'cubeSize'});
			my $v3 = int (($structure->{'hetAtoms'}->[$e]->{'coords'}->{'z'} - $structure->{'xyzExtremities'}->[2]->{'f1'}) / $analyser->{'cubeSize'});
			my $v4 = Analyser_xyzToIndex( $v1, $v2, $v3 );
			my $push_back_index = $#{$cubes->[$v4]} + 1;
			$cubes->[$v4]->[$push_back_index] = $structure->{'hetAtoms'}->[$e];
		}
	}

	my $st_dot_allAtoms_dot_size = $#{$structure->{'allAtoms'}} + 1;
	for (my $e = 0; $e < $st_dot_allAtoms_dot_size; $e++) { # for (unsigned int e = 0 ; e < st.allAtoms.size(); e++){
		# cubes[xyzToIndex(	(int)((st.allAtoms[e].coords.x  - st.xyzExtremities[0].f1)/ cubeSize),
		#					(int)((st.allAtoms[e].coords.y  - st.xyzExtremities[1].f1)/ cubeSize),
		#					(int)((st.allAtoms[e].coords.z  - st.xyzExtremities[2].f1)/ cubeSize) )].push_back(&(st.allAtoms[e]));
		my $v1 = int (($structure->{'allAtoms'}->[$e]->{'coords'}->{'x'} - $structure->{'xyzExtremities'}->[0]->{'f1'}) / $analyser->{'cubeSize'});
		my $v2 = int (($structure->{'allAtoms'}->[$e]->{'coords'}->{'y'} - $structure->{'xyzExtremities'}->[1]->{'f1'}) / $analyser->{'cubeSize'});
		my $v3 = int (($structure->{'allAtoms'}->[$e]->{'coords'}->{'z'} - $structure->{'xyzExtremities'}->[2]->{'f1'}) / $analyser->{'cubeSize'});
		my $v4 = Analyser_xyzToIndex( $v1, $v2, $v3 );
		my $push_back_index = $#{$cubes->[$v4]} + 1;
		$cubes->[$v4]->[$push_back_index] = $structure->{'allAtoms'}->[$e];
	}

	for (my $x = 0; $x < $analyser->{'noOfXCubes'}; $x++) { # for (unsigned int x = 0; x < noOfXCubes; x++){
		for (my $y = 0; $y < $analyser->{'noOfYCubes'}; $y++) { # for (unsigned int y = 0; y < noOfYCubes; y++){
			for (my $z = 0; $z < $analyser->{'noOfZCubes'}; $z++) { # for (unsigned int z = 0; z < noOfZCubes; z++){
				my $closeCubes = Analyser_getSurroundingCubes( $x, $y, $z, $cubes ); # vector<Atom*> closeCubes = getSurroundingCubes(x,y,z,cubes);
				my $xyzToIndex = Analyser_xyzToIndex( $x, $y, $z );
				my $print_close_size = $#{$closeCubes} + 1;
				my $currentCube = $cubes->[ $xyzToIndex ]; # vector<Atom*> currentCube = cubes[xyzToIndex(x,y,z)];
				my $currentCube_dot_size = $#{$currentCube} + 1;
				for (my $i = 0; $i < $currentCube_dot_size; $i++) { # for (unsigned int i = 0; i < currentCube.size(); i++){
					my $closeCubes_dot_size = $#{$closeCubes} + 1;
					for (my $j = 0; $j < $closeCubes_dot_size; $j++) { # for (unsigned int j = 0; j < closeCubes.size(); j++){
						# float xyzDistance = (currentCube[i]->coords - closeCubes[j]->coords).lengthSqr();
						my $v1 = Vec3_subtract_Vec3( $currentCube->[$i]->{'coords'}, $closeCubes->[$j]->{'coords'} );
						my $xyzDistance = Vec3_lengthSqr( $v1 );
						# float solventDistance =(currentCube[i]->solventRadii + 
						#						closeCubes[j]->solventRadii) *
						#						(currentCube[i]->solventRadii + 
						#						closeCubes[j]->solventRadii);
						my $solventDistance = ($currentCube->[$i]->{'solventRadii'} + $closeCubes->[$j]->{'solventRadii'})
							* ($currentCube->[$i]->{'solventRadii'} + $closeCubes->[$j]->{'solventRadii'});

						if ($xyzDistance < $solventDistance) { # if (xyzDistance <   solventDistance){
							# currentCube[i]->closeAtoms.push_back(closeCubes[j]);
							my $push_back_index = $#{$currentCube->[$i]->{'closeAtoms'}} + 1;
							$currentCube->[$i]->{'closeAtoms'}->[$push_back_index]->{'atomID'} = $closeCubes->[$j]->{'atomID'};
							# currentCube[i]->distances.push_back(xyzDistance);
							$push_back_index = $#{$currentCube->[$i]->{'distances'}} + 1;
							$currentCube->[$i]->{'distances'}->[$push_back_index] = $xyzDistance;
						}
					}
					$currentCube_dot_size = $#{$currentCube} + 1;
					for (my $k = $i + 1; $k < $currentCube_dot_size; $k++) { # for (unsigned int k = i + 1; k < currentCube.size(); k++){
						# float xyzDistance = (currentCube[i]->coords - currentCube[k]->coords).lengthSqr();
						my $v1 = Vec3_subtract_Vec3( $currentCube->[$i]->{'coords'}, $currentCube->[$k]->{'coords'} );
						my $xyzDistance = Vec3_lengthSqr( $v1 );
						# float solventDistance =(currentCube[i]->solventRadii + 
						#						currentCube[k]->solventRadii) *
						#						(currentCube[i]->solventRadii + 
						#						currentCube[k]->solventRadii);						
						my $solventDistance = ($currentCube->[$i]->{'solventRadii'} + $currentCube->[$k]->{'solventRadii'})
							* ($currentCube->[$i]->{'solventRadii'} + $currentCube->[$k]->{'solventRadii'});
						if ($xyzDistance < $solventDistance) { # if (xyzDistance <   solventDistance){
							# currentCube[i]->closeAtoms.push_back(currentCube[k]);
							my $push_back_index = $#{$currentCube->[$i]->{'closeAtoms'}} + 1;
							$currentCube->[$i]->{'closeAtoms'}->[$push_back_index]->{'atomID'} = $currentCube->[$k]->{'atomID'};
							# currentCube[i]->distances.push_back(xyzDistance);
							$push_back_index = $#{$currentCube->[$i]->{'distances'}} + 1;
							$currentCube->[$i]->{'distances'}->[$push_back_index] = $xyzDistance;
							# currentCube[k]->closeAtoms.push_back(currentCube[i]);
							$push_back_index = $#{$currentCube->[$k]->{'closeAtoms'}} + 1;
							$currentCube->[$k]->{'closeAtoms'}->[$push_back_index]->{'atomID'} = $currentCube->[$i]->{'atomID'};
							# currentCube[k]->distances.push_back(xyzDistance);
							$push_back_index = $#{$currentCube->[$k]->{'distances'}} + 1;
							$currentCube->[$k]->{'distances'}->[$push_back_index] = $xyzDistance;
						}
					}
					my $v = $currentCube->[$i]->{'closeAtoms'};
					my $d = $currentCube->[$i]->{'distances'};
					my $currentCube_i_distances_dot_size = $#{$currentCube->[$i]->{'distances'}} + 1;
					# quickSort(currentCube[i]->closeAtoms,
					#		currentCube[i]->distances,
					#		0,
					#		(int)currentCube[i]->distances.size()-1);
					Analyser_quickSort( $v, $d, 0, $currentCube_i_distances_dot_size - 1 );
					$currentCube->[$i]->{'closeAtoms'} = $v;
					$currentCube->[$i]->{'distances'} = $d;	
				}
			}
		}
	}

	$analyser->{'st'} = $structure;
}



sub Analyser_getSurroundingCubes { #C++ vector<Atom*> Analyser::getSurroundingCubes(int x, int y, int z, vector<Atom*> *a){

	my $x = shift;
	my $y = shift;
	my $z = shift;
	my $a = shift;
	my $cubeIndices;				# vector<Atom*> cubeIndices;

	# //x series
	# //if (xyzToIndex(x, y, z) > 0){
	# //	for (unsigned int i = 0; i < a[xyzToIndex(x, y, z)].size(); i++){
	# //		cubeIndices.push_back(a[xyzToIndex(x, y, z)][i]);
	# //	}
	# //}
	my $xyzToIndex = Analyser_xyzToIndex( $x, $y, $z+1 );
	if ($xyzToIndex > 0) { # if (xyzToIndex(x, y, z+1) > 0){
		my $a_xyzToIndex_dot_size = $#{$a->[$xyzToIndex]} + 1;
		for (my $i = 0; $i < $a_xyzToIndex_dot_size; $i++ ) { # for (unsigned int i = 0; i < a[xyzToIndex(x, y, z+1)].size(); i++){
			my $push_back_index = $#{$cubeIndices} + 1;
			$cubeIndices->[$push_back_index] = $a->[$xyzToIndex]->[$i]; # cubeIndices.push_back(a[xyzToIndex(x, y, z+1)][i]);}
		}
	}
	$xyzToIndex = Analyser_xyzToIndex( $x, $y, $z-1 );
	if ($xyzToIndex > 0) { # if (xyzToIndex(x , y , z-1) > 0){
		my $a_xyzToIndex_dot_size = $#{$a->[$xyzToIndex]} + 1;
		for (my $i = 0; $i < $a_xyzToIndex_dot_size; $i++ ) { # for (unsigned int i = 0; i < a[xyzToIndex(x, y, z-1)].size(); i++){
			my $push_back_index = $#{$cubeIndices} + 1;
			$cubeIndices->[$push_back_index] = $a->[$xyzToIndex]->[$i]; # cubeIndices.push_back(a[xyzToIndex(x, y, z-1)][i]);}
		}
	}
	$xyzToIndex = Analyser_xyzToIndex( $x, $y+1, $z );
	if ($xyzToIndex > 0) { # if (xyzToIndex(x, y+1, z) > 0){
		my $a_xyzToIndex_dot_size = $#{$a->[$xyzToIndex]} + 1;
		for (my $i = 0; $i < $a_xyzToIndex_dot_size; $i++ ) { # for (unsigned int i = 0; i < a[xyzToIndex(x, y+1, z)].size(); i++){
			my $push_back_index = $#{$cubeIndices} + 1;
			$cubeIndices->[$push_back_index] = $a->[$xyzToIndex]->[$i]; # cubeIndices.push_back(a[xyzToIndex(x, y+1, z)][i]);}
		}
	}
	$xyzToIndex = Analyser_xyzToIndex( $x, $y+1, $z+1 );
	if ($xyzToIndex > 0) { # if (xyzToIndex(x, y+1, z+1) > 0){
		my $a_xyzToIndex_dot_size = $#{$a->[$xyzToIndex]} + 1;
		for (my $i = 0; $i < $a_xyzToIndex_dot_size; $i++ ) { # for (unsigned int i = 0; i < a[xyzToIndex(x, y+1, z+1)].size(); i++){
			my $push_back_index = $#{$cubeIndices} + 1;
			$cubeIndices->[$push_back_index] = $a->[$xyzToIndex]->[$i]; # cubeIndices.push_back(a[xyzToIndex(x, y+1, z+1)][i]);}
		}
	}
	$xyzToIndex = Analyser_xyzToIndex( $x, $y+1, $z-1 );
	if ($xyzToIndex > 0) { # if (xyzToIndex(x, y+1, z-1) > 0){
		my $a_xyzToIndex_dot_size = $#{$a->[$xyzToIndex]} + 1;
		for (my $i = 0; $i < $a_xyzToIndex_dot_size; $i++ ) { # for (unsigned int i = 0; i < a[xyzToIndex(x, y+1, z-1)].size(); i++){
			my $push_back_index = $#{$cubeIndices} + 1;
			$cubeIndices->[$push_back_index] = $a->[$xyzToIndex]->[$i]; # cubeIndices.push_back(a[xyzToIndex(x, y+1, z-1)][i]);}
		}
	}
	$xyzToIndex = Analyser_xyzToIndex( $x, $y-1, $z );
	if ($xyzToIndex > 0) { # if (xyzToIndex(x, y-1, z) > 0){
		my $a_xyzToIndex_dot_size = $#{$a->[$xyzToIndex]} + 1;
		for (my $i = 0; $i < $a_xyzToIndex_dot_size; $i++ ) { # for (unsigned int i = 0; i < a[xyzToIndex(x, y-1, z)].size(); i++){
			my $push_back_index = $#{$cubeIndices} + 1;
			$cubeIndices->[$push_back_index] = $a->[$xyzToIndex]->[$i]; # cubeIndices.push_back(a[xyzToIndex(x, y-1, z)][i]);}
		}
	}
	$xyzToIndex = Analyser_xyzToIndex( $x, $y-1, $z+1 );
	if ($xyzToIndex > 0) { # if (xyzToIndex(x, y-1, z+1) > 0){
		my $a_xyzToIndex_dot_size = $#{$a->[$xyzToIndex]} + 1;
		for (my $i = 0; $i < $a_xyzToIndex_dot_size; $i++ ) { # for (unsigned int i = 0; i < a[xyzToIndex(x, y-1, z+1)].size(); i++){
			my $push_back_index = $#{$cubeIndices} + 1;
			$cubeIndices->[$push_back_index] = $a->[$xyzToIndex]->[$i]; # cubeIndices.push_back(a[xyzToIndex(x, y-1, z+1)][i]);}
		}
	}
	$xyzToIndex = Analyser_xyzToIndex( $x, $y-1, $z-1 );
	if ($xyzToIndex > 0) { # if (xyzToIndex(x, y-1, z-1) > 0){
		my $a_xyzToIndex_dot_size = $#{$a->[$xyzToIndex]} + 1;
		for (my $i = 0; $i < $a_xyzToIndex_dot_size; $i++ ) { # for (unsigned int i = 0; i < a[xyzToIndex(x, y-1, z-1)].size(); i++){
			my $push_back_index = $#{$cubeIndices} + 1;
			$cubeIndices->[$push_back_index] = $a->[$xyzToIndex]->[$i]; # cubeIndices.push_back(a[xyzToIndex(x, y-1, z-1)][i]);}
		}
	}
	# //x+1 series
	$xyzToIndex = Analyser_xyzToIndex( $x+1, $y, $z );
	if ($xyzToIndex > 0) { # if (xyzToIndex(x+1, y, z) > 0){
		my $a_xyzToIndex_dot_size = $#{$a->[$xyzToIndex]} + 1;
		for (my $i = 0; $i < $a_xyzToIndex_dot_size; $i++ ) { # for (unsigned int i = 0; i < a[xyzToIndex(x+1, y, z)].size(); i++){
			my $push_back_index = $#{$cubeIndices} + 1;
			$cubeIndices->[$push_back_index] = $a->[$xyzToIndex]->[$i]; # cubeIndices.push_back(a[xyzToIndex(x+1, y, z)][i]);}
		}
	}
	$xyzToIndex = Analyser_xyzToIndex( $x+1, $y, $z+1 );
	if ($xyzToIndex > 0) { # if (xyzToIndex(x+1, y, z+1) > 0){
		my $a_xyzToIndex_dot_size = $#{$a->[$xyzToIndex]} + 1;
		for (my $i = 0; $i < $a_xyzToIndex_dot_size; $i++ ) { # for (unsigned int i = 0; i < a[xyzToIndex(x+1, y, z+1)].size(); i++){
			my $push_back_index = $#{$cubeIndices} + 1;
			$cubeIndices->[$push_back_index] = $a->[$xyzToIndex]->[$i]; # cubeIndices.push_back(a[xyzToIndex(x+1, y, z+1)][i]);}
		}
	}
	$xyzToIndex = Analyser_xyzToIndex( $x+1, $y, $z-1 );
	if ($xyzToIndex > 0) { # if (xyzToIndex(x+1 , y , z-1) > 0){
		my $a_xyzToIndex_dot_size = $#{$a->[$xyzToIndex]} + 1;
		for (my $i = 0; $i < $a_xyzToIndex_dot_size; $i++ ) { # for (unsigned int i = 0; i < a[xyzToIndex(x+1 , y , z-1)].size(); i++){
			my $push_back_index = $#{$cubeIndices} + 1;
			$cubeIndices->[$push_back_index] = $a->[$xyzToIndex]->[$i]; # cubeIndices.push_back(a[xyzToIndex(x+1 , y , z-1)][i]);}
		}
	}
	$xyzToIndex = Analyser_xyzToIndex( $x+1, $y+1, $z );
	if ($xyzToIndex > 0) { # if (xyzToIndex(x+1, y+1, z) > 0){
		my $a_xyzToIndex_dot_size = $#{$a->[$xyzToIndex]} + 1;
		for (my $i = 0; $i < $a_xyzToIndex_dot_size; $i++ ) { # for (unsigned int i = 0; i < a[xyzToIndex(x+1, y+1, z)].size(); i++){
			my $push_back_index = $#{$cubeIndices} + 1;
			$cubeIndices->[$push_back_index] = $a->[$xyzToIndex]->[$i]; # cubeIndices.push_back(a[xyzToIndex(x+1, y+1, z)][i]);}
		}
	}
	$xyzToIndex = Analyser_xyzToIndex( $x+1, $y+1, $z+1 );
	if ($xyzToIndex > 0) { # if (xyzToIndex(x+1, y+1, z+1) > 0){
		my $a_xyzToIndex_dot_size = $#{$a->[$xyzToIndex]} + 1;
		for (my $i = 0; $i < $a_xyzToIndex_dot_size; $i++ ) { # for (unsigned int i = 0; i < a[xyzToIndex(x+1, y+1, z+1)].size(); i++){
			my $push_back_index = $#{$cubeIndices} + 1;
			$cubeIndices->[$push_back_index] = $a->[$xyzToIndex]->[$i]; # cubeIndices.push_back(a[xyzToIndex(x+1, y+1, z+1)][i]);}
		}
	}
	$xyzToIndex = Analyser_xyzToIndex( $x+1, $y+1, $z-1 );
	if ($xyzToIndex > 0) { # if (xyzToIndex(x+1, y+1, z-1) > 0){
		my $a_xyzToIndex_dot_size = $#{$a->[$xyzToIndex]} + 1;
		for (my $i = 0; $i < $a_xyzToIndex_dot_size; $i++ ) { # for (unsigned int i = 0; i < a[xyzToIndex(x+1, y+1, z-1)].size(); i++){
			my $push_back_index = $#{$cubeIndices} + 1;
			$cubeIndices->[$push_back_index] = $a->[$xyzToIndex]->[$i]; # cubeIndices.push_back(a[xyzToIndex(x+1, y+1, z-1)][i]);}
		}
	}
	$xyzToIndex = Analyser_xyzToIndex( $x+1, $y-1, $z );
	if ($xyzToIndex > 0) { # if (xyzToIndex(x+1, y-1, z) > 0){
		my $a_xyzToIndex_dot_size = $#{$a->[$xyzToIndex]} + 1;
		for (my $i = 0; $i < $a_xyzToIndex_dot_size; $i++ ) { # for (unsigned int i = 0; i < a[xyzToIndex(x+1, y-1, z)].size(); i++){
			my $push_back_index = $#{$cubeIndices} + 1;
			$cubeIndices->[$push_back_index] = $a->[$xyzToIndex]->[$i]; # cubeIndices.push_back(a[xyzToIndex(x+1, y-1, z)][i]);}
		}
	}
	$xyzToIndex = Analyser_xyzToIndex( $x+1, $y-1, $z+1 );
	if ($xyzToIndex > 0) { # if (xyzToIndex(x+1, y-1, z+1) > 0){
		my $a_xyzToIndex_dot_size = $#{$a->[$xyzToIndex]} + 1;
		for (my $i = 0; $i < $a_xyzToIndex_dot_size; $i++ ) { # for (unsigned int i = 0; i < a[xyzToIndex(x+1, y-1, z+1)].size(); i++){
			my $push_back_index = $#{$cubeIndices} + 1;
			$cubeIndices->[$push_back_index] = $a->[$xyzToIndex]->[$i]; # cubeIndices.push_back(a[xyzToIndex(x+1, y-1, z+1)][i]);}
		}
	}
	$xyzToIndex = Analyser_xyzToIndex( $x+1, $y-1, $z-1 );
	if ($xyzToIndex > 0) { # if (xyzToIndex(x+1, y-1, z-1) > 0){
		my $a_xyzToIndex_dot_size = $#{$a->[$xyzToIndex]} + 1;
		for (my $i = 0; $i < $a_xyzToIndex_dot_size; $i++ ) { # for (unsigned int i = 0; i < a[xyzToIndex(x+1, y-1, z-1)].size(); i++){
			my $push_back_index = $#{$cubeIndices} + 1;
			$cubeIndices->[$push_back_index] = $a->[$xyzToIndex]->[$i]; # cubeIndices.push_back(a[xyzToIndex(x+1, y-1, z-1)][i]);}
		}
	}
	# //x-1 series
	$xyzToIndex = Analyser_xyzToIndex( $x-1, $y, $z );
	if ($xyzToIndex > 0) { # if (xyzToIndex(x-1, y, z) > 0){
		my $a_xyzToIndex_dot_size = $#{$a->[$xyzToIndex]} + 1;
		for (my $i = 0; $i < $a_xyzToIndex_dot_size; $i++ ) { # for (unsigned int i = 0; i < a[xyzToIndex(x-1, y, z)].size(); i++){
			my $push_back_index = $#{$cubeIndices} + 1;
			$cubeIndices->[$push_back_index] = $a->[$xyzToIndex]->[$i]; # cubeIndices.push_back(a[xyzToIndex(x-1, y, z)][i]);}
		}
	}
	$xyzToIndex = Analyser_xyzToIndex( $x-1, $y, $z+1 );
	if ($xyzToIndex > 0) { # if (xyzToIndex(x-1, y, z+1) > 0){
		my $a_xyzToIndex_dot_size = $#{$a->[$xyzToIndex]} + 1;
		for (my $i = 0; $i < $a_xyzToIndex_dot_size; $i++ ) { # for (unsigned int i = 0; i < a[xyzToIndex(x-1, y, z+1)].size(); i++){
			my $push_back_index = $#{$cubeIndices} + 1;
			$cubeIndices->[$push_back_index] = $a->[$xyzToIndex]->[$i]; # cubeIndices.push_back(a[xyzToIndex(x-1, y, z+1)][i]);}
		}
	}
	$xyzToIndex = Analyser_xyzToIndex( $x-1, $y, $z-1 );
	if ($xyzToIndex > 0) { # if (xyzToIndex(x-1 , y , z-1) > 0){
		my $a_xyzToIndex_dot_size = $#{$a->[$xyzToIndex]} + 1;
		for (my $i = 0; $i < $a_xyzToIndex_dot_size; $i++ ) { # for (unsigned int i = 0; i < a[xyzToIndex(x-1 , y , z-1)].size(); i++){
			my $push_back_index = $#{$cubeIndices} + 1;
			$cubeIndices->[$push_back_index] = $a->[$xyzToIndex]->[$i]; # cubeIndices.push_back(a[xyzToIndex(x-1 , y , z-1)][i]);}
		}
	}
	$xyzToIndex = Analyser_xyzToIndex( $x-1, $y+1, $z );
	if ($xyzToIndex > 0) { # if (xyzToIndex(x-1, y+1, z) > 0){
		my $a_xyzToIndex_dot_size = $#{$a->[$xyzToIndex]} + 1;
		for (my $i = 0; $i < $a_xyzToIndex_dot_size; $i++ ) { # for (unsigned int i = 0; i < a[xyzToIndex(x-1, y+1, z)].size(); i++){
			my $push_back_index = $#{$cubeIndices} + 1;
			$cubeIndices->[$push_back_index] = $a->[$xyzToIndex]->[$i]; # cubeIndices.push_back(a[xyzToIndex(x-1, y+1, z)][i]);}
		}
	}
	$xyzToIndex = Analyser_xyzToIndex( $x-1, $y+1, $z+1 );
	if ($xyzToIndex > 0) { # if (xyzToIndex(x-1, y+1, z+1) > 0){
		my $a_xyzToIndex_dot_size = $#{$a->[$xyzToIndex]} + 1;
		for (my $i = 0; $i < $a_xyzToIndex_dot_size; $i++ ) { # for (unsigned int i = 0; i < a[xyzToIndex(x-1, y+1, z+1)].size(); i++){
			my $push_back_index = $#{$cubeIndices} + 1;
			$cubeIndices->[$push_back_index] = $a->[$xyzToIndex]->[$i]; # cubeIndices.push_back(a[xyzToIndex(x-1, y+1, z+1)][i]);}
		}
	}
	$xyzToIndex = Analyser_xyzToIndex( $x-1, $y+1, $z-1 );
	if ($xyzToIndex > 0) { # if (xyzToIndex(x-1, y+1, z-1) > 0){
		my $a_xyzToIndex_dot_size = $#{$a->[$xyzToIndex]} + 1;
		for (my $i = 0; $i < $a_xyzToIndex_dot_size; $i++ ) { # for (unsigned int i = 0; i < a[xyzToIndex(x-1, y+1, z-1)].size(); i++){
			my $push_back_index = $#{$cubeIndices} + 1;
			$cubeIndices->[$push_back_index] = $a->[$xyzToIndex]->[$i]; # cubeIndices.push_back(a[xyzToIndex(x-1, y+1, z-1)][i]);}
		}
	}
	$xyzToIndex = Analyser_xyzToIndex( $x-1, $y-1, $z );
	if ($xyzToIndex > 0) { # if (xyzToIndex(x-1, y-1, z) > 0){
		my $a_xyzToIndex_dot_size = $#{$a->[$xyzToIndex]} + 1;
		for (my $i = 0; $i < $a_xyzToIndex_dot_size; $i++ ) { # for (unsigned int i = 0; i < a[xyzToIndex(x-1, y-1, z)].size(); i++){
			my $push_back_index = $#{$cubeIndices} + 1;
			$cubeIndices->[$push_back_index] = $a->[$xyzToIndex]->[$i]; # cubeIndices.push_back(a[xyzToIndex(x-1, y-1, z)][i]);}
		}
	}
	$xyzToIndex = Analyser_xyzToIndex( $x-1, $y-1, $z+1 );
	if ($xyzToIndex > 0) { # if (xyzToIndex(x-1, y-1, z+1) > 0){
		my $a_xyzToIndex_dot_size = $#{$a->[$xyzToIndex]} + 1;
		for (my $i = 0; $i < $a_xyzToIndex_dot_size; $i++ ) { # for (unsigned int i = 0; i < a[xyzToIndex(x-1, y-1, z+1)].size(); i++){
			my $push_back_index = $#{$cubeIndices} + 1;
			$cubeIndices->[$push_back_index] = $a->[$xyzToIndex]->[$i]; # cubeIndices.push_back(a[xyzToIndex(x-1, y-1, z+1)][i]);}
		}
	}
	$xyzToIndex = Analyser_xyzToIndex( $x-1, $y-1, $z-1 );
	if ($xyzToIndex > 0) { # if (xyzToIndex(x-1, y-1, z-1) > 0){
		my $a_xyzToIndex_dot_size = $#{$a->[$xyzToIndex]} + 1;
		for (my $i = 0; $i < $a_xyzToIndex_dot_size; $i++ ) { # for (unsigned int i = 0; i < a[xyzToIndex(x-1, y-1, z-1)].size(); i++){
			my $push_back_index = $#{$cubeIndices} + 1;
			$cubeIndices->[$push_back_index] = $a->[$xyzToIndex]->[$i]; # cubeIndices.push_back(a[xyzToIndex(x-1, y-1, z-1)][i]);}
		}
	}
	return $cubeIndices; # return cubeIndices;
}



sub Analyser_new {

	# // Construction/Destruction

	$analyser->{'cubeSize'} = 7;				# float Analyser::cubeSize = 7.0f;  
	$analyser->{'ABHigh'} = 114;				# float Analyser::ABHigh = 114.0f;
	$analyser->{'ABLow'} = 40;				# float Analyser::ABLow = 40.0f;
	$analyser->{'FPBHigh'} = 0.58;				# float Analyser::FPBHigh = 0.58f;
	$analyser->{'FPBLow'} = 0.45;				# float Analyser::FPBLow = 0.45f;
	$analyser->{'FPPartial'} = 0.67;			# float Analyser::FPPartial = 0.67f;

	# //profiles scoring table
	# float Analyser::scores[18][20] = {
	#	{1.11f,1.28f,0.27f,1.30f,1.11f,0.74f,1.26f,-0.77f,-2.22f,-1.56f,-0.43f,-1.72f,-2.43f,-1.38f,-1.76f,-2.15f,-2.48f,-0.34f,-1.37f,-1.80},	//B1alpha
	#	{0.92f,0.96f,0.17f,1.07f,1.50f,1.18f,0.51f,-1.05f,-2.35f,-0.77f,-0.45f,-1.27f,-2.56f,-2.03f,-2.18f,-1.59f,-1.80f,-2.26f,-3.04f,-1.52},	//B1beta
	#	{0.96f,1.40f,0.52f,1.06f,0.93f,1.00f,0.91f,-0.54f,-2.78f,0.59f,-0.59f,-1.41f,-2.99f,-0.84f,-2.61f,-2.01f,-2.63f,-0.61f,-2.78f,-2.35},	//B1
	#
	#	{1.01f,0.87f,0.86f,0.71f,0.55f,0.41f,1.02f,-0.65f,-2.04f,-0.97f,0.15f,-0.67f,-1.33f,0.16f,-0.48f,-0.58f,-0.80f,0.82f,-0.94f,-0.11},	//B2alpha
	#	{0.83f,1.32f,1.30f,0.36f,1.08f,0.81f,0.49f,-1.52f,-2.22f,-0.86f,-0.72f,-1.14f,-0.82f,-0.79f,-0.26f,-0.20f,-2.08f,-0.05f,-0.83f,-0.41},	//B2beta
	#	{1.62f,1.04f,1.14f,0.77f,0.81f,0.66f,1.00f,-0.81f,-1.71f,-0.07f,-0.62f,-1.03f,-1.23f,-0.87f,-0.56f,-1.13f,-1.7f,0.54f,-2.12f,-0.44},	//B2
	#
	#	{0.86f,-0.22f,0.50f,0.16f,-0.02f,-0.29f,0.87f,-0.44f,-1.09f,-1.11f,-1.38f,-0.69f,-1.01f,0.16f,-0.07f,0.09f,-0.43f,0.61f,0.56f,1.10},	//B3alpha
	#	{0.07f,0.37f,1.09f,0.14f,0.26f,0.16f,-0.68f,-1.08f,-2.29f,-0.01f,-0.79f,-0.10f,-0.71f,0.52f,-0.33f,-0.42f,-0.76f,0.80f,0.35f,0.84},	//B3beta
	#	{1.12f,0.71f,1.25f,0.29f,-0.54f,-0.40f,0.23f,-0.87f,-0.61f,-0.11f,-0.98f,-0.48f,-0.61f,0.10f,0.09f,-0.46f,-0.83f,1.04f,0.08f,0.71},	//B3
	#
	#	{-1.29f,-0.85f,-0.88f,-0.30f,-0.06f,0.30f,-0.42f,0.76f,-0.46f,-0.41f,0.95f,0.39f,0.47f,-0.32f,-0.58f,-0.43f,-0.28f,-0.91f,-0.50f,-0.51},	//P1alpha
	#	{0.34f,-0.61f,-0.09f,-0.81f,0.09f,0.44f,-0.40f,0.59f,-0.22f,-0.65f,1.28f,0.95f,0.49f,-2.38f,-0.92f,-0.68f,-0.61f,-0.53f,-2.01f,-0.89},		//P1beta
	#	{-1.25f,-1.29f,-1.40f,-0.33f,-0.28f,-0.09f,-0.90f,0.49f,-0.39f,0.64f,1.29f,0.55f,0.59f,-0.57f,-0.26f,-0.59f,0.34f,-1.21f,-0.72f,-0.88},	//P1
	#
	#	{-1.09f,-1.35f,-0.55f,-0.46f,-0.59f,-0.62f,-0.27f,-0.02f,-0.58f,-0.25f,-0.70f,-0.13f,-0.38f,0.62f,-0.02f,0.2f,0.29f,0.17f,0.66f,0.56},	//P2alpha
	#	{-0.71f,-0.56f,-0.30f,-1.33f,-0.35f,0.08f,-0.76f,-0.52f,-0.87f,-1.01f,-0.87f,0.79f,0.49f,0.10f,0.00f,0.41f,-0.03f,-0.49f,0.55f,0.19},	//P2beta
	#	{-0.42f,-0.84f,-0.43f,-0.68f,-0.94f,-0.74f,-0.83f,-0.25f,-0.42f,0.44f,-0.81f,0.08f,0.17f,0.25f,0.51f,0.28f,0.51f,0.20f,0.47f,0.24},	//P2
	#
	#	{-1.26f,-1.81f,-1.70f,-1.37f,-2.36f,-1.25f,-0.90f,0.44f,0.63f,0.05f,-0.17f,-0.20f,0.16f,0.29f,0.32f,0.60f,0.44f,-0.06f,0.07f,-0.20},		//Ealpha
	#	{0.81f,-0.83f,-0.03f,-1.60f,-1.39f,-1.66f,-0.62f,0.14f,1.75f,-0.88f,-0.04f,-0.17f,0.65f,-0.12f,0.01f,-0.37f,-0.30f,-0.76f,-1.54f,-1.12},	//Ebeta
	#	{-2.06f,-1.63f,-1.04f,-1.14f,-1.63f,0.80f,-1.30f,0.14f,1.10f,0.25f,-0.35f,0.08f,0.34f,-0.03f,0.41f,0.04f,0.23f,-0.41f,-0.10f,-0.41}
	# };
	$analyser->{'scores'} = [
		[1.11,1.28,0.27,1.30,1.11,0.74,1.26,-0.77,-2.22,-1.56,-0.43,-1.72,-2.43,-1.38,-1.76,-2.15,-2.48,-0.34,-1.37,-1.80],	# //B1alpha
		[0.92,0.96,0.17,1.07,1.50,1.18,0.51,-1.05,-2.35,-0.77,-0.45,-1.27,-2.56,-2.03,-2.18,-1.59,-1.80,-2.26,-3.04,-1.52],	# //B1beta
		[0.96,1.40,0.52,1.06,0.93,1.00,0.91,-0.54,-2.78,0.59,-0.59,-1.41,-2.99,-0.84,-2.61,-2.01,-2.63,-0.61,-2.78,-2.35],	# //B1

		[1.01,0.87,0.86,0.71,0.55,0.41,1.02,-0.65,-2.04,-0.97,0.15,-0.67,-1.33,0.16,-0.48,-0.58,-0.80,0.82,-0.94,-0.11],	# //B2alpha
		[0.83,1.32,1.30,0.36,1.08,0.81,0.49,-1.52,-2.22,-0.86,-0.72,-1.14,-0.82,-0.79,-0.26,-0.20,-2.08,-0.05,-0.83,-0.41],	# //B2beta
		[1.62,1.04,1.14,0.77,0.81,0.66,1.00,-0.81,-1.71,-0.07,-0.62,-1.03,-1.23,-0.87,-0.56,-1.13,-1.7,0.54,-2.12,-0.44],	# //B2

		[0.86,-0.22,0.50,0.16,-0.02,-0.29,0.87,-0.44,-1.09,-1.11,-1.38,-0.69,-1.01,0.16,-0.07,0.09,-0.43,0.61,0.56,1.10],	# //B3alpha
		[0.07,0.37,1.09,0.14,0.26,0.16,-0.68,-1.08,-2.29,-0.01,-0.79,-0.10,-0.71,0.52,-0.33,-0.42,-0.76,0.80,0.35,0.84],	# //B3beta
		[1.12,0.71,1.25,0.29,-0.54,-0.40,0.23,-0.87,-0.61,-0.11,-0.98,-0.48,-0.61,0.10,0.09,-0.46,-0.83,1.04,0.08,0.71],	# //B3

		[-1.29,-0.85,-0.88,-0.30,-0.06,0.30,-0.42,0.76,-0.46,-0.41,0.95,0.39,0.47,-0.32,-0.58,-0.43,-0.28,-0.91,-0.50,-0.51],	# //P1alpha
		[0.34,-0.61,-0.09,-0.81,0.09,0.44,-0.40,0.59,-0.22,-0.65,1.28,0.95,0.49,-2.38,-0.92,-0.68,-0.61,-0.53,-2.01,-0.89],	# //P1beta
		[-1.25,-1.29,-1.40,-0.33,-0.28,-0.09,-0.90,0.49,-0.39,0.64,1.29,0.55,0.59,-0.57,-0.26,-0.59,0.34,-1.21,-0.72,-0.88],	# //P1

		[-1.09,-1.35,-0.55,-0.46,-0.59,-0.62,-0.27,-0.02,-0.58,-0.25,-0.70,-0.13,-0.38,0.62,-0.02,0.2,0.29,0.17,0.66,0.56],	# //P2alpha
		[-0.71,-0.56,-0.30,-1.33,-0.35,0.08,-0.76,-0.52,-0.87,-1.01,-0.87,0.79,0.49,0.10,0.00,0.41,-0.03,-0.49,0.55,0.19],	# //P2beta
		[-0.42,-0.84,-0.43,-0.68,-0.94,-0.74,-0.83,-0.25,-0.42,0.44,-0.81,0.08,0.17,0.25,0.51,0.28,0.51,0.20,0.47,0.24],	# //P2

		[-1.26,-1.81,-1.70,-1.37,-2.36,-1.25,-0.90,0.44,0.63,0.05,-0.17,-0.20,0.16,0.29,0.32,0.60,0.44,-0.06,0.07,-0.20],	# //Ealpha
		[0.81,-0.83,-0.03,-1.60,-1.39,-1.66,-0.62,0.14,1.75,-0.88,-0.04,-0.17,0.65,-0.12,0.01,-0.37,-0.30,-0.76,-1.54,-1.12],	# //Ebeta
		[-2.06,-1.63,-1.04,-1.14,-1.63,0.80,-1.30,0.14,1.10,0.25,-0.35,0.08,0.34,-0.03,0.41,0.04,0.23,-0.41,-0.10,-0.41]
	];

	# private:
	$analyser->{'SASARes'} = 0;				# unsigned int SASARes; # //records resolution of SASA calculation 

	# ////for cubing atoms///
	$analyser->{'noOfXCubes'} = 0;			# unsigned int noOfXCubes;
	$analyser->{'$noOfYCubes'} = 0;			# unsigned int noOfYCubes;
	$analyser->{'$noOfZCubes'} = 0;			# unsigned int noOfZCubes;

	$analyser->{'st'} = $structure;			# public: Structure st;

	return;
}



sub Analyser_partition { #C++ int Analyser::partition(vector<Atom*> & v,vector<float> &d, int start, int end)

	my $v = shift;
	my $d = shift;
	my $start = shift;
	my $end = shift;

	use constant RAND_MAX => 32767;
	my $r = rand( RAND_MAX );			# int r = rand();
	$r = $r % ($end - $start);			# r = r%(end - start);
	$r = $r + $start;				# r = r + start;
	my $partition = $d->[$r];			# float partition = d[r];
	my $buffer;					# float buffer;
	my $atomBuffer;					# Atom* atomBuffer;

	my $s = $start - 1;				# register int s = start-1; 
	my $e = $end + 1;				# register int e = end+1;
	while (1) {					# while(true){ 

		$e--;					# do { e--; }
		while ($d->[$e] > $partition) {		# while(d[e] > partition);
			$e--;
		}

		$s++;					# do { s++; } 
		while ($d->[$s] < $partition) {		# while(d[s] < partition);
			$s++; 
		}

		if ($s < $e) {				# if(s < e){
			$buffer = $d->[$s];		# buffer = d[s];
			$d->[$s] = $d->[$e];		# d[s] = d[e];
			$d->[$e] = $buffer;		# d[e] = buffer;
			$atomBuffer = $v->[$s];		# atomBuffer = v[s];
			$v->[$s] = $v->[$e];		# v[s] = v[e];
			$v->[$e] = $atomBuffer;		# v[e] = atomBuffer;
		} else { 
			return $e;			# return e; 
		}
	}
	return -1;					#return -1;
}



# bool Analyser::predTM2D(unsigned int mveAve, unsigned int TMLengthMin, unsigned int TMLengthMax,
#						float minScoreDiff, float minAreaDiff, float aveProfMin, float aveProfMax,
#						float aveRepMin, float aveRepMax, float lenDivAreaMin, float lenDivAreaMax,
#						float aveSASAMin)
sub Analyser_predTM2D_valpred1 { #C++ void Analyser::Analyser::predTM2D(...)

	$structure = $analyser->{'st'};

	my $mveAve = $TMPred2D->{'mveAve'};
	my $TMLengthMin = $TMPred2D->{'TMLengthMin_for_valpred1'};
	my $TMLengthMax = $TMPred2D->{'TMLengthMax_for_valpred1'};
	my $minScoreDiff = $TMPred2D->{'minScoreDiff_for_valpred1'};
	my $minAreaDiff = $TMPred2D->{'minAreaDiff_for_valpred1'};
	my $aveProfMin = $TMPred2D->{'aveProfRangeMin'};
	my $aveProfMax = $TMPred2D->{'aveProfRangeMax'};
	my $aveRepMin = $TMPred2D->{'aveRepRangeMin'};
	my $aveRepMax = $TMPred2D->{'aveRepRangeMax'};
	my $lenDivAreaMin = $TMPred2D->{'lenDivAreaMin'};
	my $lenDivAreaMax = $TMPred2D->{'lenDivAreaMax'};
	my $aveSASAMin = $TMPred2D->{'minAveSasa'};
	my $ASA_are_available = $TMPred2D->{'ASA_are_available'};

	@tmDomain_start = (); # initialise

	my $res = $structure->{'allResidues'}; 							# vector<Residue> &res = st.allResidues;
	# ///fill residue moving average scores
	my $res_dot_size = $#{$res} + 1;
	for (my $i = 0; $i < $res_dot_size; $i++) {						# for (unsigned int i = 0; i < res.size(); i++){
		$res->[$i]->{'mov_avg_1_score'}->[0] = 0;						# res[i].aveScores[0] = 0.0f;
		$res->[$i]->{'mov_avg_1_score'}->[1] = 0;						# res[i].aveScores[1] = 0.0f;
		my $j = 0;									# unsigned int j = 0;
		while ($j < $mveAve) {								# while ( j < mveAve){
			if (!( ($j+$i) < $res_dot_size )) {					# if(!((j + i) < res.size())){
				$j++;								# j++;
												# continue;
			} else {
				$res->[$i]->{'mov_avg_1_score'}->[0] += $res->[$i+$j]->{'scores'}->[0]; # res[i].aveScores[0] += res[i+j].scores[0];
				$res->[$i]->{'mov_avg_1_score'}->[1] += $res->[$i+$j]->{'scores'}->[1]; # res[i].aveScores[1] += res[i+j++].scores[1];
				$j++;
			}
		}
		$res->[$i]->{'mov_avg_1_score'}->[0] = $res->[$i]->{'mov_avg_1_score'}->[0] / $j; 		# res[i].aveScores[0] = res[i].aveScores[0]/j;
		$res->[$i]->{'mov_avg_1_score'}->[1] = $res->[$i]->{'mov_avg_1_score'}->[1] / $j;		# res[i].aveScores[1] = res[i].aveScores[1]/j;
		$res->[$i]->{'isTransMem'} = 0;							# res[i].isTransMem = false;	# //reset tm predictions
	}

	# ////////find TM domains///////////

	# /////start with the largest TM length and decrease
	for ( my $length = $TMLengthMax; $length >= $TMLengthMin; $length-- ) {			# for(unsigned int length = TMLengthMax; length >= TMLengthMin; length--){
		my $st_dot_chains_dot_size = $#{$structure->{'chains'}} + 1;
		for ( my $c = 0; $c < $st_dot_chains_dot_size; $c++ ) {				# for (unsigned int c = 0; c < st.chains.size(); c++){
			my $residues = Structure_getResChain_int( $c );				# vector<Residue*> residues = st.getResChain(c);
			# ////for each length test all criteria on window at every residue
			my $residues_dot_size = $#{$residues} + 1;
			for ( my $i = 0; $i < $residues_dot_size; $i++ ) {			# for(unsigned int i = 0; i < residues.size(); i++){
				my $area;			# vector<float> area; # //vector to record the area between the curves
				my $alwaysGreater = 1;		# bool alwaysGreater = true; # //test that the difference between every residue in the domain is greater than the minimum distance allowed///
				my $isRec = 0;			# bool isRec = false; # //Make sure the same residue is not tested twice
				
				# //pre run test - make sure no residues have already been dubbed transmembrane
				my $j;				# unsigned int j;
				# for(j = 0; j < length && i + length + 1 < residues.size() && !isRec; j++){
				for ( $j = 0; ($j < $length) && (($i + $length + 1) < $residues_dot_size) && (!$isRec); $j++ ) {
					if ( $residues->[$i + $j]->{'isTransMem'} ) { 		# if (residues[i+j]->isTransMem){
						$isRec = 1;					# isRec = true;
					}
				}
				# ////test a window of the TMLength
				# for(j = 0; j < length && i + length + 1 < residues.size() && !isRec; j++){
				for ( $j = 0; ($j < $length) && (($i + $length + 1) < $residues_dot_size) && (!$isRec); $j++ ) { 
					# ////calculate the area between the curves
					# float repScoreArea = (residues[i+j]->aveScores[1] +
					#					 residues[i+j+1]->aveScores[1])/2;
					my $repScoreArea = ( $residues->[$i+$j]->{'mov_avg_1_score'}->[1] + $residues->[$i+$j+1]->{'mov_avg_1_score'}->[1] ) / 2;
					# float profScoreArea = (residues[i+j]->aveScores[0] +
					#					 residues[i+j+1]->aveScores[0])/2;
					my $profScoreArea = ( $residues->[$i+$j]->{'mov_avg_1_score'}->[0] + $residues->[$i+$j+1]->{'mov_avg_1_score'}->[0] ) / 2;
					my $push_back_index = $#{$area} + 1;
					$area->[$push_back_index] = $repScoreArea - $profScoreArea; # area.push_back(repScoreArea - profScoreArea);
					# if(residues[i+j]->aveScores[1] - residues[i+j]->aveScores[0] <	minScoreDiff){
					if ( ($residues->[$i+$j]->{'mov_avg_1_score'}->[1] - $residues->[$i+$j]->{'mov_avg_1_score'}->[0]) < $minScoreDiff ) {
						$alwaysGreater = 0;				# alwaysGreater = false;
					}
				}
				# ///sum the total area between the 2 curves
				my $totalArea = 0;						# float totalArea = 0.0f;
				my $area_dot_size = $#{$area} + 1;
				for ( my $s = 0; $s < $area_dot_size; $s++ ) {			# for (unsigned int s = 0; s < area.size(); s++){
					$totalArea += $area->[$s];				# totalArea+=area[s];
				}
				# //////calculate average profiles, repimps scores & SASA///////
				my $aveRepScore = 0;						# float aveRepScore = 0;
				my $aveProfScore = 0;						# float aveProfScore = 0;
				my $aveSASA = 0;						# float aveSASA = 0;
				$residues_dot_size = $#{$residues} + 1;
				# for(int q = i; q <= i+length && q < residues.size(); q++){
				for ( my $q = $i; ($q <= ($i + $length)) && ($q < $residues_dot_size); $q++ ) {
					$aveRepScore += $residues->[$q]->{'mov_avg_1_score'}->[1];	# aveRepScore += residues[q]->aveScores[1];
					$aveProfScore += $residues->[$q]->{'mov_avg_1_score'}->[0];	# aveProfScore += residues[q]->aveScores[0];
					$aveSASA += $residues->[$q]->{'side_ASA'}; 		# aveSASA += residues[q]->side_ASA;
				}

				if ($length != 0) {
					$aveSASA = $aveSASA / $length;					# aveSASA = aveSASA / length;
					$aveRepScore = $aveRepScore / $length;				# aveRepScore = aveRepScore / length;
					$aveProfScore = $aveProfScore / $length;			# aveProfScore = aveProfScore / length;
				}

				# //////calculate TMlength / sum of area difference///////
				my $lDivA = 0;
				if ($totalArea != 0) {
					$lDivA = $length / $totalArea;				# float lDivA = length/totalArea;
				}

				# if(	alwaysGreater &&
				#	totalArea > minAreaDiff &&
				#	aveProfScore > aveProfMin && aveProfScore < aveProfMax	&&
				#	aveRepScore > aveRepMin && aveRepScore < aveRepMax		&&
				#	lDivA > lenDivAreaMin && lDivA < lenDivAreaMax		&&
				#	aveSASA > aveSASAMin){
				if ( $alwaysGreater && ($totalArea > $minAreaDiff) &&
					($aveProfScore > $aveProfMin) && ($aveProfScore < $aveProfMax) &&
					($aveRepScore > $aveRepMin) && ($aveRepScore < $aveRepMax) &&
					($lDivA > $lenDivAreaMin) && ($lDivA < $lenDivAreaMax) ) {

					my $is_trans_mem = 0;
					if ($ASA_are_available == 1) { # which means that aveSASA was calculated from building a helix out of the residues
						if ($aveSASA > $aveSASAMin) {
							$is_trans_mem = 1;
						}
					} else { # which means that the PROFILES3D and REPIMPS scores were read in from the input valpred_scores.txt file
							# which means that aveSASA is not available, so it will be ignored.
							# thus, it is possible that the VALPRED1 algorithm will give different results
							# depending on whether the input is the residues or the scores previously computed from the residues.
						$is_trans_mem = 1;
					}

					if ($is_trans_mem == 1) {
						for ( my $tm = $i; $tm <= ($length + $i); $tm++ ) {	# for (unsigned int tm = i; tm <= length+i; tm++){
							$residues->[$tm]->{'isTransMem'} = 1;		# residues[tm]->isTransMem = true;
						}
						push( @tmDomain_start, $i ); # for valpred1, if 2 predicted TMS are stuck together, separate them, don't predict 1 long TMS
					}
				}
			}
		}
	}
	return;		# return false;
}



sub Analyser_predTM2D_valpred2 { #C++ void Analyser::Analyser::predTM2D(...)

	$structure = $analyser->{'st'};

	my $mov_avg_1 = $TMPred2D->{'mov_avg_1'};
	my $mov_avg_2 = $TMPred2D->{'mov_avg_2'};
	my $mov_avg_3 = $TMPred2D->{'mov_avg_3'};
	my $mov_avg_4 = $TMPred2D->{'mov_avg_4'};
	my $TMLengthMin = $TMPred2D->{'TMLengthMin_for_valpred2'};
	my $TMLengthMax = $TMPred2D->{'TMLengthMax_for_valpred2'};
	my $minScoreDiff_movavg1 = $TMPred2D->{'minScoreDiff_movavg1_for_valpred2'};
	my $minScoreDiff_movavg2 = $TMPred2D->{'minScoreDiff_movavg2_for_valpred2'};
	my $minAreaDiff = $TMPred2D->{'minAreaDiff_for_valpred2'};
	my $min_avg_area_movavg1 = $TMPred2D->{'min_avg_area_movavg1'};
	my $min_avg_area_movavg2 = $TMPred2D->{'min_avg_area_movavg2'};
	my $max_nonTMS_length = $TMPred2D->{'max_nonTMS_length'};
	my $min_nonTMS_length = $TMPred2D->{'min_nonTMS_length'};
	my $min_nonTMS_score_diff = $TMPred2D->{'min_nonTMS_score_diff'};
	my $min_nonTMS_area = $TMPred2D->{'min_nonTMS_area'};
	my $twilight_area_per_residue_lower_limit = $TMPred2D->{'twilight_area_per_residue_lower_limit'};
	my $twilight_area_per_residue_upper_limit = $TMPred2D->{'twilight_area_per_residue_upper_limit'};
	my $ignore_twilight_area = $TMPred2D->{'ignore_twilight_area'};
	my $search_area_for_neighbour_tms = $TMPred2D->{'search_area_for_neighbour_tms'};
	my $helix_length_min = $TMPred2D->{'helix_length_min'};
	my $min_score_diff_for_TMS_ends = $TMPred2D->{'min_score_diff_for_TMS_ends'};

	my $res = $structure->{'allResidues'}; 							# vector<Residue> &res = st.allResidues;
	my $res_dot_size = $#{$res} + 1;

	# Algorithm for finding transmembrane domains : "AREA ALGORITHM"
	# Scores are given for each residue's propensity for the REPIMPS membrane environment and the PROFILES3D watery anti-membrane environment.
	# Scores are smoothed with a moving average.
	# Areas between the positive REPIMPS and negative PROFLES3D scores are calculated.
	# A large enough area for long enough, will mark those residues as transmembrane.

	# Calculate the moving-average-1 and moving-average-2 scores
	# If defaults are being used, then moving-average-1 = 10 and moving-average-2 = 5.
	# Also calculate a small moving-average (default = 5) for places inside a putative transmembrane segment,
	# where there is an area where PROFILES3D is greater than REPIMPS. 
	# This will be a break in the putative transmembrane helix.

	my $st_dot_chains_dot_size = $#{$structure->{'chains'}} + 1;
	for ( my $c = 0; $c < $st_dot_chains_dot_size; $c++ ) {	
		my $residues = Structure_getResChain_int( $c );
		my $residues_dot_size = $#{$residues} + 1;
		for (my $i = 0; $i < $residues_dot_size; $i++) {

			# Calculate the moving-average-1 scores, moving-average-1 = 10

			$residues->[$i]->{'mov_avg_1_score'}->[0] = 0; # PROFILES3D
			$residues->[$i]->{'mov_avg_1_score'}->[1] = 0; # REPIMPS
			my $mov_avg_1_half = int($mov_avg_1/2);
			my $mov_avg_1_count = 0;
			for (my $j = ($i-$mov_avg_1_half); $j <= ($i+$mov_avg_1_half); $j++) {
				if (($j >= 0) && ($j < $residues_dot_size)) {
					$residues->[$i]->{'mov_avg_1_score'}->[0] += $residues->[$j]->{'scores'}->[0];
					$residues->[$i]->{'mov_avg_1_score'}->[1] += $residues->[$j]->{'scores'}->[1];
					$mov_avg_1_count++;
				}
			}
			if (($mov_avg_1 % 2) == 0) {
				my $j = $i - $mov_avg_1_half - 1;
				if ($j >= 0) {
					$residues->[$i]->{'mov_avg_1_score'}->[0] += $residues->[$j]->{'scores'}->[0];
					$residues->[$i]->{'mov_avg_1_score'}->[1] += $residues->[$j]->{'scores'}->[1];
					$mov_avg_1_count = $mov_avg_1_count + 0.5;
				}
				$j = $i + $mov_avg_1_half + 1;
				if ($j < $residues_dot_size) {
					$residues->[$i]->{'mov_avg_1_score'}->[0] += $residues->[$j]->{'scores'}->[0];
					$residues->[$i]->{'mov_avg_1_score'}->[1] += $residues->[$j]->{'scores'}->[1];
					$mov_avg_1_count = $mov_avg_1_count + 0.5;
				}
			}
			$residues->[$i]->{'mov_avg_1_score'}->[0] = $residues->[$i]->{'mov_avg_1_score'}->[0] / $mov_avg_1_count;
			$residues->[$i]->{'mov_avg_1_score'}->[1] = $residues->[$i]->{'mov_avg_1_score'}->[1] / $mov_avg_1_count;

			# Calculate the moving-average-2 scores, moving-average-2 = 5

			$residues->[$i]->{'mov_avg_2_score'}->[0] = 0; # PROFILES3D
			$residues->[$i]->{'mov_avg_2_score'}->[1] = 0; # REPIMPS
			my $mov_avg_2_half = int($mov_avg_2/2);
			my $mov_avg_2_count = 0;
			for (my $j = ($i-$mov_avg_2_half); $j <= ($i+$mov_avg_2_half); $j++) {
				if (($j >= 0) && ($j < $residues_dot_size)) {
					$residues->[$i]->{'mov_avg_2_score'}->[0] += $residues->[$j]->{'scores'}->[0];
					$residues->[$i]->{'mov_avg_2_score'}->[1] += $residues->[$j]->{'scores'}->[1];
					$mov_avg_2_count++;
				}
			}
			if (($mov_avg_2 % 2) == 0) {
				my $j = $i - $mov_avg_2_half - 1;
				if ($j >= 0) {
					$residues->[$i]->{'mov_avg_2_score'}->[0] += $residues->[$j]->{'scores'}->[0];
					$residues->[$i]->{'mov_avg_2_score'}->[1] += $residues->[$j]->{'scores'}->[1];
					$mov_avg_2_count = $mov_avg_2_count + 0.5;
				}
				$j = $i + $mov_avg_2_half + 1;
				if ($j < $residues_dot_size) {
					$residues->[$i]->{'mov_avg_2_score'}->[0] += $residues->[$j]->{'scores'}->[0];
					$residues->[$i]->{'mov_avg_2_score'}->[1] += $residues->[$j]->{'scores'}->[1];
					$mov_avg_2_count = $mov_avg_2_count + 0.5;
				}
			}
			$residues->[$i]->{'mov_avg_2_score'}->[0] = $residues->[$i]->{'mov_avg_2_score'}->[0] / $mov_avg_2_count;
			$residues->[$i]->{'mov_avg_2_score'}->[1] = $residues->[$i]->{'mov_avg_2_score'}->[1] / $mov_avg_2_count;

			# Calculate the moving-average-3 scores, moving-average-3 = 5

			$residues->[$i]->{'mov_avg_3_score'}->[0] = 0; # PROFILES3D
			$residues->[$i]->{'mov_avg_3_score'}->[1] = 0; # REPIMPS
			my $mov_avg_3_half = int($mov_avg_3/2);
			my $mov_avg_3_count = 0;
			for (my $j = ($i-$mov_avg_3_half); $j <= ($i+$mov_avg_3_half); $j++) {
				if (($j >= 0) && ($j < $residues_dot_size)) {
					$residues->[$i]->{'mov_avg_3_score'}->[0] += $residues->[$j]->{'scores'}->[0];
					$residues->[$i]->{'mov_avg_3_score'}->[1] += $residues->[$j]->{'scores'}->[1];
					$mov_avg_3_count++;
				}
			}
			if (($mov_avg_3 % 2) == 0) {
				my $j = $i - $mov_avg_3_half - 1;
				if ($j >= 0) {
					$residues->[$i]->{'mov_avg_3_score'}->[0] += $residues->[$j]->{'scores'}->[0];
					$residues->[$i]->{'mov_avg_3_score'}->[1] += $residues->[$j]->{'scores'}->[1];
					$mov_avg_3_count = $mov_avg_3_count + 0.5;
				}
				$j = $i + $mov_avg_3_half + 1;
				if ($j < $residues_dot_size) {
					$residues->[$i]->{'mov_avg_3_score'}->[0] += $residues->[$j]->{'scores'}->[0];
					$residues->[$i]->{'mov_avg_3_score'}->[1] += $residues->[$j]->{'scores'}->[1];
					$mov_avg_3_count = $mov_avg_3_count + 0.5;
				}
			}
			$residues->[$i]->{'mov_avg_3_score'}->[0] = $residues->[$i]->{'mov_avg_3_score'}->[0] / $mov_avg_3_count;
			$residues->[$i]->{'mov_avg_3_score'}->[1] = $residues->[$i]->{'mov_avg_3_score'}->[1] / $mov_avg_3_count;

			# Calculate the moving-average-4 scores, moving-average-4 = 10

			$residues->[$i]->{'mov_avg_4_score'}->[0] = 0; # PROFILES3D
			$residues->[$i]->{'mov_avg_4_score'}->[1] = 0; # REPIMPS
			my $mov_avg_4_half = int($mov_avg_4/2);
			my $mov_avg_4_count = 0;
			for (my $j = ($i-$mov_avg_4_half); $j <= ($i+$mov_avg_4_half); $j++) {
				if (($j >= 0) && ($j < $residues_dot_size)) {
					$residues->[$i]->{'mov_avg_4_score'}->[0] += $residues->[$j]->{'scores'}->[0];
					$residues->[$i]->{'mov_avg_4_score'}->[1] += $residues->[$j]->{'scores'}->[1];
					$mov_avg_4_count++;
				}
			}
			if (($mov_avg_4 % 2) == 0) {
				my $j = $i - $mov_avg_4_half - 1;
				if ($j >= 0) {
					$residues->[$i]->{'mov_avg_4_score'}->[0] += $residues->[$j]->{'scores'}->[0];
					$residues->[$i]->{'mov_avg_4_score'}->[1] += $residues->[$j]->{'scores'}->[1];
					$mov_avg_4_count = $mov_avg_4_count + 0.5;
				}
				$j = $i + $mov_avg_4_half + 1;
				if ($j < $residues_dot_size) {
					$residues->[$i]->{'mov_avg_4_score'}->[0] += $residues->[$j]->{'scores'}->[0];
					$residues->[$i]->{'mov_avg_4_score'}->[1] += $residues->[$j]->{'scores'}->[1];
					$mov_avg_4_count = $mov_avg_4_count + 0.5;
				}
			}
			$residues->[$i]->{'mov_avg_4_score'}->[0] = $residues->[$i]->{'mov_avg_4_score'}->[0] / $mov_avg_4_count;
			$residues->[$i]->{'mov_avg_4_score'}->[1] = $residues->[$i]->{'mov_avg_4_score'}->[1] / $mov_avg_4_count;
		}

		my $residues_indices = Structure_getResChain_indices( $c );
		my $index = $residues_indices->{'min'};
		for (my $i = 0; $i < $residues_dot_size; $i++) {
			$res->[$index]->{'mov_avg_1_score'}->[0] = $residues->[$i]->{'mov_avg_1_score'}->[0];
			$res->[$index]->{'mov_avg_1_score'}->[1] = $residues->[$i]->{'mov_avg_1_score'}->[1];
			$res->[$index]->{'mov_avg_2_score'}->[0] = $residues->[$i]->{'mov_avg_2_score'}->[0];
			$res->[$index]->{'mov_avg_2_score'}->[1] = $residues->[$i]->{'mov_avg_2_score'}->[1];
			$res->[$index]->{'mov_avg_3_score'}->[0] = $residues->[$i]->{'mov_avg_3_score'}->[0];
			$res->[$index]->{'mov_avg_3_score'}->[1] = $residues->[$i]->{'mov_avg_3_score'}->[1];
			$res->[$index]->{'mov_avg_4_score'}->[0] = $residues->[$i]->{'mov_avg_4_score'}->[0];
			$res->[$index]->{'mov_avg_4_score'}->[1] = $residues->[$i]->{'mov_avg_4_score'}->[1];
			$res->[$index]->{'isTransMem'} = 0;
			$res->[$index]->{'isTransMem_area_mov_avg_1'} = 0;
			$res->[$index]->{'isTransMem_area_mov_avg_2'} = 0;
			$res->[$index]->{'is_hydrophilic_area_mov_avg_3'} = 0;
			$res->[$index]->{'is_hydrophilic_area_mov_avg_4'} = 0;
			$index++;
		}
	}

	# find TM domains using moving-average-1 scores, moving-average-1 = 10

	for ( my $c = 0; $c < $st_dot_chains_dot_size; $c++ ) {
		my $residues = Structure_getResChain_int( $c );
		my $residues_dot_size = $#{$residues} + 1;
		for ( my $i = 0; $i < $residues_dot_size; $i++ ) {
			my $totalArea = 0;
			my $tms_length = 0;
			my $alwaysGreater = 1;
			for ( my $j = 0; (($i+$j) < $residues_dot_size) && ($alwaysGreater == 1); $j++ ) { 
				# ////calculate the area between the curves
				my $repScoreArea = $residues->[$i+$j]->{'mov_avg_1_score'}->[1];
				my $profScoreArea = $residues->[$i+$j]->{'mov_avg_1_score'}->[0];
				my $area_for_this_residue = $repScoreArea - $profScoreArea;
				if ( $area_for_this_residue >= $minScoreDiff_movavg1 ) {
					$tms_length++;
					$totalArea += $area_for_this_residue;
				} else {
					$alwaysGreater = 0;
				}
			}

			if ($tms_length > 0) {
				my $avg_area = $totalArea / $tms_length;
				if (($totalArea > $minAreaDiff) && ($avg_area > $min_avg_area_movavg1)) {
					for ( my $tm = $i; $tm <= ($tms_length + $i); $tm++ ) {
						$residues->[$tm]->{'isTransMem_area_mov_avg_1'} = 1;
					}
				}
			}
		}
		my $residues_indices = Structure_getResChain_indices( $c );
		my $index = $residues_indices->{'min'};
		for (my $i = 0; $i < $residues_dot_size; $i++) {
			$res->[$index]->{'isTransMem_area_mov_avg_1'} = $residues->[$i]->{'isTransMem_area_mov_avg_1'};
			$index++;
		}
	}

	# find TM domains using moving-average-2 scores, moving-average-2 = 5

	for ( my $c = 0; $c < $st_dot_chains_dot_size; $c++ ) {
		my $residues = Structure_getResChain_int( $c );
		my $residues_dot_size = $#{$residues} + 1;
		for ( my $i = 0; $i < $residues_dot_size; $i++ ) {
			my $totalArea = 0;
			my $tms_length = 0;
			my $alwaysGreater = 1;
			for ( my $j = 0; (($i+$j) < $residues_dot_size) && ($alwaysGreater == 1); $j++ ) { 
				# ////calculate the area between the curves
				my $repScoreArea = $residues->[$i+$j]->{'mov_avg_2_score'}->[1];
				my $profScoreArea = $residues->[$i+$j]->{'mov_avg_2_score'}->[0];
				my $area_for_this_residue = $repScoreArea - $profScoreArea;
				if ( $area_for_this_residue >= $minScoreDiff_movavg2 ) {
					$tms_length++;
					$totalArea += $area_for_this_residue;
				} else {
					$alwaysGreater = 0;
				}
			}

			if ($tms_length > 0) {
				my $avg_area = $totalArea / $tms_length;
				if (($totalArea > $minAreaDiff) && ($avg_area > $min_avg_area_movavg2)) {
					for ( my $tm = $i; $tm <= ($tms_length + $i); $tm++ ) {
						$residues->[$tm]->{'isTransMem_area_mov_avg_2'} = 1;
					}
				}
			}
		}
		my $residues_indices = Structure_getResChain_indices( $c );
		my $index = $residues_indices->{'min'};
		for (my $i = 0; $i < $residues_dot_size; $i++) {
			$res->[$index]->{'isTransMem_area_mov_avg_2'} = $residues->[$i]->{'isTransMem_area_mov_avg_2'};
			$index++;
		}
	}

	# Decide the TMS by looking at the results of the two moving averages algorithms,
	# before removing hydrophilic areas.

	for ( my $c = 0; $c < $st_dot_chains_dot_size; $c++ ) {
		my $residues_indices = Structure_getResChain_indices( $c );
		my $index = $residues_indices->{'min'};
		for (my $index = $residues_indices->{'min'}; $index <= $residues_indices->{'max'}; $index++) {
			if (($res->[$index]->{'isTransMem_area_mov_avg_1'} == 1) || ($res->[$index]->{'isTransMem_area_mov_avg_2'} == 1)) {
				$res->[$index]->{'isTransMem'} = 1;
			} else {
				$res->[$index]->{'isTransMem'} = 0;
			}
		}
	}

	# find NON-TM hydrophilic domains inside TM domains, 
	# wherever moving-average-3 PROFILES3D line is above the moving-average-3 REPIMPS line
	# (moving average = 5)

	for ( my $c = 0; $c < $st_dot_chains_dot_size; $c++ ) {
		my $residues = Structure_getResChain_int( $c );
		my $residues_dot_size = $#{$residues} + 1;
		for ( my $i = 0; $i < $residues_dot_size; $i++ ) {
			my $in_TMS = 0;
			for ( my $j = ($i - $max_nonTMS_length + 1); $j < ($i + $max_nonTMS_length - 1); $j++ ) {
				if (($j >= 0) && ($j < $residues_dot_size)) {
					if ($residues->[$j]->{'isTransMem'} == 1) {
						$in_TMS = 1;
					}
				}
			}
			if ($in_TMS == 1) {
				for ( my $j = ($i - $max_nonTMS_length + 1); $j <= $i; $j++ ) {
					my $totalArea = 0;
					my $non_tms_length = 0;
					my $alwaysGreater = 1;
					for ( my $k = 0; ($k < $max_nonTMS_length) && ($alwaysGreater == 1); $k++ ) { 
						if ((($i+$k) >= 0) && (($i+$k) < $residues_dot_size)) {
							# ////calculate the area between the curves
							my $repScoreArea = $residues->[$i+$k]->{'mov_avg_3_score'}->[1];
							my $profScoreArea = $residues->[$i+$k]->{'mov_avg_3_score'}->[0];
							my $area_for_this_residue = $profScoreArea - $repScoreArea;
							if ( $area_for_this_residue >= $min_nonTMS_score_diff ) {
								$non_tms_length++;
								$totalArea += $area_for_this_residue;
							} else {
								$alwaysGreater = 0;
							}
						}
					}
					# //////calculate TMlength / sum of area difference///////
					my $lDivA = 0;
					if ($totalArea != 0) {
						$lDivA = $non_tms_length / $totalArea;
					}

					if ($totalArea > $min_nonTMS_area) {
						for ( my $tm = $i; $tm <= ($i + $non_tms_length - 1); $tm++ ) {
							$residues->[$tm]->{'is_hydrophilic_area_mov_avg_3'} = 1;
						}
					}
				}
			}
		}
		my $residues_indices = Structure_getResChain_indices( $c );
		my $index = $residues_indices->{'min'};
		for (my $i = 0; $i < $residues_dot_size; $i++) {
			$res->[$index]->{'is_hydrophilic_area_mov_avg_3'} = $residues->[$i]->{'is_hydrophilic_area_mov_avg_3'};
			$index++;
		}
	}

	# Removing the identified hydrophilic areas inside hydrophobic helices.

	for ( my $i = 0; $i < $res_dot_size; $i++ ) {
		if ($res->[$i]->{'isTransMem'} == 1) {
			if ($res->[$i]->{'is_hydrophilic_area_mov_avg_3'} == 1) {
				$res->[$i]->{'isTransMem'} = 0;
			}
		}
	}

	# Shorten the ends of TMS.
	# A large window for MOV_AVG_1 means that there could be weakly hydrophobic areas 
	# in a TMS in the middle of stronger hydrophobic areas.
	# This may be where the TMS was broken into 2 pieces.
	# So trim the weakly hydrophobic areas at the ends of TMS 
	# (if these weakly hydrophobic areas are now the ends of 2 TMS, not in the middle of 1 TMS).
	# Otherwise, they will make the TMS too long for no good reason,
	# and then it will be vulnerable to being broken in 2 in the middle.

	for ( my $c = 0; $c < $st_dot_chains_dot_size; $c++ ) {
		my $residues = Structure_getResChain_int( $c );
		my $residues_dot_size = $#{$residues} + 1;
		my @start_tms;
		my @end_tms;
		my $in_tms = 0;
		for ( my $i = 0; $i < $residues_dot_size; $i++ ) {
			if ($residues->[$i]->{'isTransMem'} == 1) {
				if ($in_tms == 0) {
					push( @start_tms, $i );
					$in_tms = 1;
				}
			} else { # $residues->[$i]->{'isTransMem'} == 1
				if ($in_tms == 1) {
					push( @end_tms, $i );
					$in_tms = 0;
				}
			}
		}
		if ($in_tms == 1) {
			my $i = $residues_dot_size - 1;
			push( @end_tms, $i );
		}
		for ( my $i = 0; $i < @start_tms; $i++ ) {
			my $continue_removing = 1;
			for ( my $j = $start_tms[$i]; $j <= $end_tms[$i]; $j++ ) {
				if ($continue_removing == 1) {
					my $repScoreArea = $residues->[$j]->{'mov_avg_1_score'}->[1];
					my $profScoreArea = $residues->[$j]->{'mov_avg_1_score'}->[0];
					my $area_for_this_residue = $repScoreArea - $profScoreArea;
					if ( $area_for_this_residue < $min_score_diff_for_TMS_ends ) {
						$residues->[$j]->{'isTransMem'} = 0;
					} else {
						$continue_removing = 0;
					}
				}
			}
		}
		for ( my $i = 0; $i < @start_tms; $i++ ) {
			my $continue_removing = 1;
			for ( my $j = $end_tms[$i]; $j >= $start_tms[$i]; $j-- ) {
				if ($continue_removing == 1) {
					my $repScoreArea = $residues->[$j]->{'mov_avg_1_score'}->[1];
					my $profScoreArea = $residues->[$j]->{'mov_avg_1_score'}->[0];
					my $area_for_this_residue = $repScoreArea - $profScoreArea;
					if ( $area_for_this_residue < $min_score_diff_for_TMS_ends ) {
						$residues->[$j]->{'isTransMem'} = 0;
					} else {
						$continue_removing = 0;
					}
				}
			}
		}

		# record that the trimmed ends are not TMS

		my $residues_indices = Structure_getResChain_indices( $c );
		my $index = $residues_indices->{'min'};
		for (my $i = 0; $i < $residues_dot_size; $i++) {
			$res->[$index]->{'isTransMem'} = $residues->[$i]->{'isTransMem'};
			$index++;
		}
	}

	# break up TMS that are too long.
	# break them at the "weakest" point where 
	# the area between the moving-average-4 REPIMPS and PROFILES3D lines is the least
	# (moving average = 6)

	for ( my $c = 0; $c < $st_dot_chains_dot_size; $c++ ) {
		my $residues = Structure_getResChain_int( $c );
		my $residues_dot_size = $#{$residues} + 1;

		my $all_TMS_respect_max_length = 0;
		while ($all_TMS_respect_max_length == 0) {
			my $found_TMS_too_long = 0;

			for ( my $i = 0; $i < $residues_dot_size; $i++ ) {
				my $previous_residue_is_tms = 0;
				if ($i > 0) {
					if ($residues->[$i-1]->{'isTransMem'} == 1) {
						$previous_residue_is_tms = 1;
					} else {
						$previous_residue_is_tms = 0;
					}
				}
				if (($residues->[$i]->{'isTransMem'} == 1) && ($previous_residue_is_tms == 0)) {
					my $this_tms_length = 0;
					my $reached_end_of_this_tms = 0;
					my $j = $i;
					while ($reached_end_of_this_tms == 0) {
						if ($residues->[$j]->{'isTransMem'} == 1) {
							$this_tms_length++;
							$j++;
							if ($j >= $residues_dot_size) {
								$reached_end_of_this_tms = 1;
							}
						} else {
							$reached_end_of_this_tms = 1;
						}
					}
					my $this_tms_start = $i;
					my $this_tms_end = $i + $this_tms_length - 1;

					if ($this_tms_length > $TMLengthMax) {
						$found_TMS_too_long = 1;

						my $broke_up_this_tms = 0;
						my @this_area_array;
						for (my $k = $this_tms_start; $k <= $this_tms_end; $k++) {
							if (($k >= 0) && ($k < $residues_dot_size)) {
								my $repScoreArea = $residues->[$k]->{'mov_avg_4_score'}->[1];
								my $profScoreArea = $residues->[$k]->{'mov_avg_4_score'}->[0];
								$this_area_array[$k] = $repScoreArea - $profScoreArea;
							}
						}

						# try to break up this long TMS at a minima point that is nearest the middle of the TMS,
						# but may not be the most minimum point

						if ($broke_up_this_tms == 0) {
							my $finished_looking_for_a_minima = 0;
							my $k = $this_tms_start + int(($this_tms_end - $this_tms_start) / 2);
							my $k_left = $k;
							my $k_right = $k;
							my $go_which_way = 'left';
							while ($finished_looking_for_a_minima == 0) {
								my $this_residue_is_a_minima = 0;
								if ((($k - 1) >= 0) && (($k - 1) >= $this_tms_start) && (($k + 1) < $this_tms_end) && (($k + 1) < $residues_dot_size)) {
									if (($this_area_array[$k - 1] < $this_area_array[$k]) && ($this_area_array[$k] < $this_area_array[$k + 1])) {
										$this_residue_is_a_minima = 1;
									}
								}
								if ($this_residue_is_a_minima == 1) {
									$residues->[$k]->{'isTransMem'} = 0;
									$residues->[$k]->{'is_hydrophilic_area_mov_avg_4'} = 1;
									$broke_up_this_tms = 1;
									$finished_looking_for_a_minima = 1;
								}
								if ($go_which_way eq 'left') {
									$k = $k_left - 1;
									$k_left = $k;
									$go_which_way = 'right';
								} else { # $go_which_way eq 'right'
									$k = $k_right + 1;
									$k_right = $k;
									$go_which_way = 'left';
								}
								if (($k < 0) || ($k < $this_tms_start) || ($k >= $this_tms_end) || ($k >= $residues_dot_size)) {
									$finished_looking_for_a_minima = 1;
								}
							}
						}

						# if didn't find a local minima,
						# then just break up the TMS at the minimum point

						if ($broke_up_this_tms == 0) {
							my $finished_looking_for_the_minimum_point = 0;
							my $k = $this_tms_start + int(($this_tms_end - $this_tms_start) / 2);
							my $k_left = $k;
							my $k_right = $k;
							my $go_which_way = 'left';
							my $min_area = $this_area_array[$k];
							my $min_area_index = $k;
							while ($finished_looking_for_the_minimum_point == 0) {
								my $this_area = $this_area_array[$k];
								if ($this_area < $min_area) {
									$min_area = $this_area;
									$min_area_index = $k;
								}
								if ($go_which_way eq 'left') {
									$k = $k_left - 1;
									$k_left = $k;
									$go_which_way = 'right';
								} else { # $go_which_way eq 'right'
									$k = $k_right + 1;
									$k_right = $k;
									$go_which_way = 'left';
								}
								if (($k < 0) || ($k < $this_tms_start) || ($k >= $this_tms_end) || ($k >= $residues_dot_size)) {
									$finished_looking_for_the_minimum_point = 1;
								}
							}
							$residues->[$min_area_index]->{'isTransMem'} = 0;
							$residues->[$min_area_index]->{'is_hydrophilic_area_mov_avg_4'} = 1;
							if (($min_area_index - 1) >= 0) {
								$residues->[$min_area_index - 1]->{'isTransMem'} = 0;
								$residues->[$min_area_index - 1]->{'is_hydrophilic_area_mov_avg_4'} = 1;
								$broke_up_this_tms = 1;
							}
							if (($min_area_index + 1) < $residues_dot_size) {
								$residues->[$min_area_index + 1]->{'isTransMem'} = 0;
								$residues->[$min_area_index + 1]->{'is_hydrophilic_area_mov_avg_4'} = 1;
								$broke_up_this_tms = 1;
							}
						}
					}
				}
			}

			if ($found_TMS_too_long == 0) {
				$all_TMS_respect_max_length = 1;
			}
		}
		my $residues_indices = Structure_getResChain_indices( $c );
		my $index = $residues_indices->{'min'};
		for (my $i = 0; $i < $residues_dot_size; $i++) {
			$res->[$index]->{'is_hydrophilic_area_mov_avg_4'} = $residues->[$i]->{'is_hydrophilic_area_mov_avg_4'};
			$res->[$index]->{'isTransMem'} = $residues->[$i]->{'isTransMem'};
			$index++;
		}
	}

	# Removing the hydrophilic area in hydrophobic helices may have created some orphan bits of hydrophobic helix.
	# Remove any hydrophobic helices that are less than 5 residues long.
	# (eg. _XXXX_, _XXX_, _XX_, _X_)

	for ( my $c = 0; $c < $st_dot_chains_dot_size; $c++ ) {

		my $residues = Structure_getResChain_int( $c );
		my $residues_dot_size = $#{$residues} + 1;

		for ( my $i = 0; $i < $residues_dot_size; $i++ ) {

			my $previous_residue_is_tms = 0;
			if ($i > 0) {
				if ($residues->[$i-1]->{'isTransMem'} == 1) {
					$previous_residue_is_tms = 1;
				} else {
					$previous_residue_is_tms = 0;
				}
			}

			if (($residues->[$i]->{'isTransMem'} == 1) && ($previous_residue_is_tms == 0)) {

				my $this_tms_length = 0;
				my $reached_end_of_this_tms = 0;
				my $j = $i;

				while ($reached_end_of_this_tms == 0) {
					if ($residues->[$j]->{'isTransMem'} == 1) {
						$this_tms_length++;
						$j++;
						if ($j >= $residues_dot_size) {
							$reached_end_of_this_tms = 1;
						}
					} else {
						$reached_end_of_this_tms = 1;
					}
				}

				if ($this_tms_length < $helix_length_min) {
					for (my $k = $i; $k <= ($i + $this_tms_length - 1); $k++) {
						$residues->[$k]->{'isTransMem'} = 0;
					}
				}
			}
		}
		my $residues_indices = Structure_getResChain_indices( $c );
		my $index = $residues_indices->{'min'};
		for (my $i = 0; $i < $residues_dot_size; $i++) {
			$res->[$index]->{'isTransMem'} = $residues->[$i]->{'isTransMem'};
			$index++;
		}
	}

	# For each TMS, what is the total area (between REPIMPS and PROFILES3D lines), and the area per residues?
	# Try to distinguish between TMS helices (larger area) and hydrophobic soluble protein helices (less area).

	my $total_area_mov_avg_1 = 0;
	my $total_area_mov_avg_2 = 0;
	my $total_residues = 0;
	my $start_tms;
	for ( my $c = 0; $c < $st_dot_chains_dot_size; $c++ ) {

		my $residues = Structure_getResChain_int( $c );
		my $residues_dot_size = $#{$residues} + 1;

		my @isTransMem_new = ();
		my @isTransMem_above_min_tms = ();
		my @sum_of_total_areas_per_residue = ();

		for ( my $i = 0; $i < $residues_dot_size; $i++ ) {
			$isTransMem_new[$i] = $residues->[$i]->{'isTransMem'};
		}

		for ( my $i = 0; $i < $residues_dot_size; $i++ ) {

			my $previous_residue_is_tms = 0;
			if ($i > 0) {
				if ($residues->[$i - 1]->{'isTransMem'} == 1) {
					$previous_residue_is_tms = 1;
				} else {
					$previous_residue_is_tms = 0;
				}
			}

			if (($residues->[$i]->{'isTransMem'} == 1) && ($previous_residue_is_tms == 0)) {
				my $mov_avg_1_area = $residues->[$i]->{'mov_avg_1_score'}->[1] - $residues->[$i]->{'mov_avg_1_score'}->[0];
				my $mov_avg_2_area = $residues->[$i]->{'mov_avg_2_score'}->[1] - $residues->[$i]->{'mov_avg_2_score'}->[0];
				$total_area_mov_avg_1 = $mov_avg_1_area;
				$total_area_mov_avg_2 = $mov_avg_2_area;
				$total_residues = 1;
				$start_tms = $i
			} elsif ($residues->[$i]->{'isTransMem'} == 1) {
				my $mov_avg_1_area = $residues->[$i]->{'mov_avg_1_score'}->[1] - $residues->[$i]->{'mov_avg_1_score'}->[0];
				my $mov_avg_2_area = $residues->[$i]->{'mov_avg_2_score'}->[1] - $residues->[$i]->{'mov_avg_2_score'}->[0];
				$total_area_mov_avg_1 += $mov_avg_1_area;
				$total_area_mov_avg_2 += $mov_avg_2_area;
				$total_residues += 1;
			} elsif (($residues->[$i]->{'isTransMem'} == 0) && ($previous_residue_is_tms == 1)) {
				my $total_area_mov_avg_1_per_residue = $total_area_mov_avg_1 / $total_residues;
				my $total_area_mov_avg_2_per_residue = $total_area_mov_avg_2 / $total_residues;
				my $end_tms = $i - 1;
				for ( my $j = $start_tms; $j <= $end_tms; $j++ ) {
					$sum_of_total_areas_per_residue[$j] = $total_area_mov_avg_1_per_residue + $total_area_mov_avg_2_per_residue;
				}
				$total_area_mov_avg_1 = 0;
				$total_area_mov_avg_2 = 0;
				$total_residues = 0;
			}
		}
		if ($residues_dot_size > 0) {
			if ($residues->[$residues_dot_size - 1]->{'isTransMem'} == 1) {
				my $total_area_mov_avg_1_per_residue = $total_area_mov_avg_1 / $total_residues;
				my $total_area_mov_avg_2_per_residue = $total_area_mov_avg_2 / $total_residues;
				my $end_tms = $residues_dot_size - 1;
				for ( my $j = $start_tms; $j <= $end_tms; $j++ ) {
					$sum_of_total_areas_per_residue[$j] = $total_area_mov_avg_1_per_residue + $total_area_mov_avg_2_per_residue;
				}
			}
		}

		# for each TMS whose areas were calculated above,
		# see if the area is big enough for a TMS.
		# if it is not big enough, then assume it is a hydrophobic helix in a soluble protein, not a TMS,
		# and unmark it as being TMS.

		if ($ignore_twilight_area == 0) {
			my $i = 0;
			while ( $i < $residues_dot_size ) {
				if ($isTransMem_new[$i] == 1) {
					if ($sum_of_total_areas_per_residue[$i] < $twilight_area_per_residue_lower_limit) {
						$isTransMem_new[$i] = 0; # if area per residue is less than a certain amount, then unmark this tms, it is not a tms
						$i++;
					} else {
						if (($sum_of_total_areas_per_residue[$i] >= $twilight_area_per_residue_lower_limit) && ($sum_of_total_areas_per_residue[$i] < $twilight_area_per_residue_upper_limit)) {
							# if area per residue is in the twilight zone, then see if it is surrounded by other TMS.
							# if so, keep it being marked as a TMS. 
							# if not, then assume it is a hydrophobic part of a water soluble protein, and unmark it as TMS.
							my $is_tms_protein = 0;
							my $start_of_this_tms = $i;
							my $start_search_on_left = $start_of_this_tms - $search_area_for_neighbour_tms + 1;
							my $end_search_on_left = $start_of_this_tms;
							for ( my $j = $start_search_on_left; $j < $end_search_on_left; $j++ ) {
								if (($j >= 0) && ($j < $residues_dot_size)) {
									if ($isTransMem_new[$j] == 1) {
										if ($sum_of_total_areas_per_residue[$j] >= $twilight_area_per_residue_upper_limit) {
											$is_tms_protein = 1;
										}
									}
								}
							}
							my $still_inside_this_tms = 1;
							while ($still_inside_this_tms == 1) {
								if ($isTransMem_new[$i] == 1) {
									$i++;
								} else {
									$still_inside_this_tms = 0;
								}
								if ($i >= $residues_dot_size) {
									$still_inside_this_tms = 0;
								}
							}
							my $end_of_this_tms = $i - 1;
							my $start_search_on_right = $end_of_this_tms;
							my $end_search_on_right = $end_of_this_tms + $search_area_for_neighbour_tms - 1;
							for ( my $j = $start_search_on_right; $j < $end_search_on_right; $j++ ) {
								if (($j >= 0) && ($j < $residues_dot_size)) {
									if ($isTransMem_new[$j] == 1) {
										if ($sum_of_total_areas_per_residue[$j] >= $twilight_area_per_residue_upper_limit) {
											$is_tms_protein = 1;
										}
									}
								}
							}
							if ($is_tms_protein == 0) {
								for ( my $j = $start_of_this_tms; $j <= $end_of_this_tms; $j++ ) {
									# this tms is not a tms, so unmark it as tms
									$isTransMem_new[$j] = 0;
								}
							}
							if (($start_search_on_left < 0) && ($end_search_on_right > $residues_dot_size)) {
								$is_tms_protein = 1; # if there is not enough area around this tms for there to be other tms, then assume that this really is a tms
							}
						} else { # area per residue is greater than $twilight_area_per_residue_upper_limit, so do keep it marked as TMS
							$i++;
						}
					}
				} else {
					$i++;
				}
			}
		}

		for ( my $i = 0; $i < $residues_dot_size; $i++ ) {
			$residues->[$i]->{'isTransMem'} = $isTransMem_new[$i];
		}

		my $residues_indices = Structure_getResChain_indices( $c );
		my $index = $residues_indices->{'min'};
		for (my $i = 0; $i < $residues_dot_size; $i++) {
			$res->[$index]->{'isTransMem'} = $residues->[$i]->{'isTransMem'};
			$index++;
		}
	}

	# For each TMS, if it is shorter than the minimum number of residues allowed for a TMS, then mark it as not TMS.

	for ( my $c = 0; $c < $st_dot_chains_dot_size; $c++ ) {

		my $residues = Structure_getResChain_int( $c );
		my $residues_dot_size = $#{$residues} + 1;

		my $start_of_this_tms = 0;
		my $length_of_this_tms = 0;
		for ( my $i = 0; $i < $residues_dot_size; $i++ ) {

			my $previous_residue_is_tms = 0;
			if ($i > 0) {
				if ($residues->[$i - 1]->{'isTransMem'} == 1) {
					$previous_residue_is_tms = 1;
				} else {
					$previous_residue_is_tms = 0;
				}
			}

			if (($residues->[$i]->{'isTransMem'} == 1) && ($previous_residue_is_tms == 0)) {
				$start_of_this_tms = $i;
				$length_of_this_tms = 1;
			} elsif ($residues->[$i]->{'isTransMem'} == 1) {
				$length_of_this_tms++;
			} elsif (($residues->[$i]->{'isTransMem'} == 0) && ($previous_residue_is_tms == 1)) {
				if ($length_of_this_tms < $TMLengthMin) {
					for ( my $j = $start_of_this_tms; $j < $i; $j++ ) {
						$residues->[$j]->{'isTransMem'} = 0;
					}
				}
			}
		}

		my $residues_indices = Structure_getResChain_indices( $c );
		my $index = $residues_indices->{'min'};
		for (my $i = 0; $i < $residues_dot_size; $i++) {
			$res->[$index]->{'isTransMem'} = $residues->[$i]->{'isTransMem'};
			$index++;
		}
	}

	return;		# return false;
}



sub Analyser_quickSort { #C++ void Analyser::quickSort(vector<Atom*> &v, vector<float> &d, int start, int end)

	my $v = shift;
	my $d = shift;
	my $start = shift;
	my $end = shift;

	if ($start < $end) {							# if( start < end )
		my $boundary = Analyser_partition( $v, $d, $start, $end );	# int boundary = partition(v,d,start,end);
		Analyser_quickSort( $v, $d, $start, $boundary );		# quickSort( v, d, start, boundary );
		Analyser_quickSort( $v, $d, $boundary + 1, $end );		# quickSort( v, d, boundary + 1, end  );
	}
}



sub Analyser_SASA { #C++ void Analyser::SASA(int resolution)

# Description of the PROFILES3D algorithm
# =======================================
#
# PROFILES3D is an algorithm that measures the compatibility of an amino acid sequence with a 3D protein structure, 
# producing an overall compatibility score (Bowie et al. 1991). 
# There is an extended version, called VERIFY3D, that plots a sliding window of PROFILES3D scores 
# for the residues in a sequence being compared to its proposed 3D structure, 
# in order to assess that 3D protein model (Lüthy et al. 1992). 
# It is assumed that the protein structure is for a water-soluble protein found in an aqueous environment. 
# The steps involved in calculating PROFILES3D scores are explained in Bowie et al. (1991) 
# and are briefly explained below.
#
# STEP 1 : For each residue in the 3D structure: 
# (i) calculate the total area of the side chain that is buried by other protein atoms; 
# (ii) calculate the fraction of the side chain area that is covered by polar atoms or water; and 
# (iii) read what the secondary structure (alpha helix, beta sheet or coil) is for that residue. 
# These calculations involve finding the solvent accessible surface area (SASA) 
# of the residue's side chain atoms in this particular 3D structure. 
# They also involve calculating the surface area of the side chain that is buried, 
# which is calculated by comparing the SASA in this protein 
# to what it is in the Gly-X-Gly tripeptides of Eisenberg et al. (1989).
#
# STEP 2 : Assign an environment class to each residue, 
# given the surface area and hydropathy of exposed and buried side chains, 
# as calculated in step 1, and given its secondary structure. 
# This converts a 3D protein structure into a 1D sequence of environment classes in place of amino acid residues.
#	ENVIRONMENT CLASS CODE	ENVIRONMENT CLASS DESCRIPTION
#		B1		buried; covered by non-polar atoms
#		B2		buried; covered by polar and non-polar atoms
#		B3		buried; covered by polar atoms or water
#		P1		partially exposed; covered by non-polar atoms
#		P2		partially exposed; covered by polar atoms or water
#		E		exposed
#
# STEP 3 : Make the PROFILES3D 3D structure profile for this 3D protein structure 
# (Gribskov et al. 1987; Gribskov et al. 1990). 
# For each residue in the protein, use its environment class to assign a score to each of the 20 amino acids. 
# The score is read from the precompiled PROFILES3D 3D-1D scoring table. 
# The score indicates how likely it is that the amino acid could exist in this position in the 3D structure. 
# For example, it is unlikely that a charged residue will be found buried in a nonpolar environment. 
# This produces a table that is similar to the PAM mutation data matrix of Dayhoff et al. (1978), 
# with a column for each of the 20 amino acids, 
# but where the rows represent positions in a 3D structure rather than being a list of the 20 amino acids.
#
# STEP 4 : For a given amino acid sequence for which you wish to see how compatible it is with the 3D structure 
# for which the PROFILES3D table has just been built, calculate the compatibility. 
# For each residue in the PROFILES3D table, observe what amino acid from the sequence is found in this residue position, 
# read off the PROFILES3D score for this amino acid, and add the score to the running total compatibility score.  
# The PROFILES3D method allows for insertions and deletions, 
# and so the dynamic programming algorithm (Needleman and Wunsch 1970; Smith and Waterman 1981) 
# needs to be used to find the optimal alignment of the sequence with the PROFILES3D table.
#
# THE PROFILES 3D-1D SCORING TABLE : 
# The 3D-1D scoring table contains the generally observed probabilities of finding each of the 20 amino acids 
# in each of the PROFILES3D environments classes and protein secondary structure types. 
# These were derived from 16 protein structures that were known as of 1991, 
# and aligned homologous sequences available at that time (Bowie et al. 1991; Lüthy et al. 1991).

# Description of the REPIMPS algorithm
# ====================================
#
# REPIMPS = reverse-environment prediction of integral membrane protein structure
#
# REPIMPS is an algorithm that measures the compatibility of an amino acid sequence 
# with a 3D protein structure of an integral membrane protein, 
# producing an overall compatibility score (Dastmalchi et al. 2001). 
# It works in the same way as the PROFILES3D method, 
# except that the exposed surface area of a residue is considered to be exposed to a hydrophobic solvent 
# instead of being exposed to a hydrophilic aqueous environment. 
# Dastmalchi et al. (2001) describe an algorithm for correcting the buried and aqueous phase side chains 
# calculated for PROFILES3D so that they can be used in the REPIMPS model of hydrophobic solvent 
# to calculate compatibility of a sequence to a membrane protein 3D structure.
# The REPIMPS method also describes how to predict the position of transmembrane segments in an amino acid sequence 
# (Dastmalchi et al. 2001). The sequence is folded into a single, long alpha helix, 
# with side chains automatically placed in a favorable rotamer. 
# The compatibility of each residue to a lipid alpha-helical environment is scored. 
# A sliding window of residue compatilibility scores with a smoothing function 
# is used to smooth the consecutive residue scores to a curve. 
# The particular combination of functions and parameters described in the Dastmalchi et al. (2001) paper 
# is a power spectrum of the fourier analysis at 100°, which corresponds to the periodicity of an alpha-helix, 
# and using a sliding window of 17 residues, 
# which is an average residue length of the transmembrane portion of transmembrane helices. 
# Transmembrane helices are predicted to occur where there are peaks in the curves 
# of smoothed residue compatibility scores, 
# as these are the areas with a high compatibility to the hydrophobic lipid environoment of the REPIMPS model.

# Solvent Accessible Surface Area (SASA) algorithm
# ================================================
#
# A method for calculating the van der Waals surface of a protein molecule 
# was described by Lee and Richards (1971). In their calculations, 
# for each atom a sphere is centred at the atom's position, 
# and the atom is considered to have a radius equal to the sum of that of the atom and that of the solvent molecule. 
# The protein molecule's surface area is computed as the locus of the center of a solvent molecule 
# as it rolls over the protein making the maximum contact possible in the given 3D structure. 
# The Lee and Richards method involves modelling protein molecule slices. 
# Each slice traces the line of the molecule's surface at that cross-section slice. 
# The Shrake and Rupley method (Shrake and Rupley 1973) describes having “test points” 
# on the surface of each atom of each residue, 
# and testing whether the test point falls inside the radius of another atom. 
# Shrake and Rupley (1973) use a set of 92 test points per atom, 
# that are nearly uniformly distributed around the surface area of the atom. 
# Hydrogen atoms are ignored (because there are so many of them, and because they are so small 
# that they are considered to make only a small dint in the surface area of the atom they are bonded to).

# Cubing algorithm
# ================
#
# The solvent accessible surface area (SASA) calculations for a 3D protein 
# require calculating how much of each atom is buried by other atoms in the protein. 
# If there are n atoms, then comparing each atom's coordinates to every other atom's coordinates 
# will take on the order of n2 amount of time 
# (or more precisely, the time is proportional n(n-1)/2 (Katz and Levinthal 1972)). 
# Using a cubing algorithm will reduce the number of comparisons necessary and thus reduce the processing time to do it. 
# The 3D protein is divided into 3D cubes. To identify buried surface areas, 
# each atom is compared only to atoms in its own cube and to atoms in the 26 neighbouring cubes to this atom's cube. 
# This "cube-testing" procedure was first developed at the Massachusetts Institute of Technology (Levinthal 1966) 
# and is used extensively in computer programs dealing with distances and energies of atoms in molecules.
# To further speed up the calculations, these neighbouring atoms and neighbouring cube atoms are sorted by distance. 
# The Quicksort algorithm is a sorting method that could be used. 
# On average, Quicksort takes n log n amount of time, depending on how unsorted or reverse sorted the set of data is.

# References
# ==========
#
# Bowie,J.U., Lüthy,R. and Eisenberg,D. (1991) A method to identify protein sequences that fold into a known three-dimensional structure. Science, 253, 164-170.
# Dastmalchi,S., Morris,M.B. and Church,W.B. (2001) Modelling of the structural features of integral-membrane proteins using REPIMPS. Protein Science 10, 1529-1538.
# Dayhoff,M.O., Schwartz,R.M. and Orcutt,B.C. (1978) A model of evolutionary change in proteins. In Dayhoff,M.O. (ed) Atlas of Protein Sequence and Structure vol 5 suppl 3 (pp.345-352) Washington, DC: National Biomedical Research Foundation.
# Eisenberg,D., Wesson,M. and Yamashita,M. (1989) Interpretation of protein folding and binding with atomic solvation parameters. Chemica Scripta 29A, 217-221.
# Gribskov,M., McLachlan,A.D. and Eisenberg,D. (1987) Profile analysis: Detection of distantly related proteins. Biochemistry 84, 4355-4358.
# Gribskov,M., Lüthy,R. and Eisenberg,D. (1990) Profile Analysis. Methods in Enzymology 183, 146-159.
# Lee,B. and Richards,F.M. (1971) The Interpretation of Protein Structures: Estimation of Static Accessibility. Journal of Molecular Biology 55, 379-490. 
# Levinthal,C. (1966) Molecular model-building by computer. Scientific American. 214, 42-52.
# Luthy,R., McLachlan,A.D. and Eisenberg,D. (1991) Secondary structure-based profiles: Use of structure-conserving scoring tables in searching protein sequence databases for structural similarities. Proteins: Structure, Function, and Genetics 10, 229-239. 
# Luthy,R., Bowie,J.U. and Eisenberg,D. (1992) Assessment of protein models with three-dimensional profiles. Nature, 356, 83-85. 
# Needleman,S.B. and Wunsch,C.D. (1970) A general method applicable to the search for similarities in the amino acid sequence of two proteins. Journal of Molecular Biology 48, 443. 
# Shrake,A. and Rupley,J.A. (1973) Environment and exposure to solvent of protein atoms lysozyme and insulin. Journal of Molecular Biology 79, 361-371.
# Smith,T.F. and Waterman,M.S. (1981) Identification of common molecular subsequences. Advances in Applied Mathematics 2, 482. 

	# for each point on the surface of a residue:
	#    if the point is not buried, then it is a surface point and it is considered to be a polar point
	#    if the point is buried, then it is a surface point, 
	#				and it is considered to be a polar point only if it really is polar (truly-polar)
	# for each residue :
	#    'solventSurface' solvent surface area = the surface area of this residue when not in a structure/helix
	#				= add up the surface area of each atom in this residue when in this particular amino acid
	#    surfPoints = surface point, whether or not it is buried in this particular helix or 3D protein structure
	#    ASA = accessible surface area = ASA for this residue in this particular helix or 3D protein structure
	#		= ( (num points - num buried points) / num points ) * solvent surface area
	#		= ( num not-buried points / num points ) * 'solventSurface' solvent surface area
	#    exp_SC_points = exposed-to-polar-environment side chain points 
	#				= buried polar points + all not buried points (because assumed to be in polar environment0
	#    all_SC_points = all side chain points = CA.surface_points + surface_points of each side chain atom
	#				= all not buried points + all buried points
	#    side_ASA = side chain ASA = CA.ASA + ASA of each side chain atom, 
	#				in this particular helix or 3D protein structure
	#    areaBuried = side chain GlyXGly - side chain ASA
	#    fractionPolar = fraction polar of side chain = exposed side chain points / all side chain points 
	#				= exposed-to-polar-environment side chain points / all side chain surface (whether buried or not, whether polar all not)
	#				= (buried polar side chain points + all not buried points) / all side chain points
	#    repimpsFP = repimps fraction polar = fraction polar - (side chain ASA / side chain GlyXGly)
	#				= fraction polar - fraction of non-buried ASA compared to a free standing residue having no other residue atoms in the 3D structure or helix covering this residue's side chains
	#				= ( (buried polar side chain points + all not buried points) / all side chain points)
	#						- fraction of non-buried ASA compared to a free standing residue
	#				repimpsFP : should = exposed-to-polar-environment side chain points / all side chain surface (whether buried or not, whether polar all not)
	#						= buried polar side chain points / all side chain points
	#				repimpsFP numerator does not include any of the non-buried side chain points, 
	#				because they are now assumed to be exposed to non-polar solvent and are not polar, 
	#				and so not included in the polar fraction for repimps.
	#				(side chain ASA / side chain GlyXGly) is the fraction that is exposed to solvent.
	#				In fractionPolar, it is assumed to be polar - exposed to polar solvent.
	#				But in repimpsFP, it is assumed to be non-polar - exposed to non-polar solvent.
	#				So remove this now-non-polar fraction from the fractionPolar total fraction,
	#				to give the remaining reimpsFP polar fraction
	#				as being only the fraction of the side chain that is buried AND polar.

	# Here is the key for how PROFILES3D and REPIMPS assume the side chain residues look like :
	#
	#	=	:	a buried residue
	#	.	:	a non-buried residue
	#	+	:	a polar residue
	#	 	:	a non-polar residue
	#
	#
	# Here is how PROFILES3D assumes the side chain residues look like :
	#
	# 	buried vs. non-buried :	==============.....................................======================
	#	polar  vs. non-polar  :       ++      +++++++++++++++++++++++++++++++++++++   ++         ++      
	#
	#
	# Here is how REPIMPS assumes the side chain residues look like :
	#
	# 	buried vs. non-buried :	==============.....................................======================
	#	polar  vs. non-polar  :       ++                                              ++         ++      
	#
	#
	# As shown above, polar buried vs. non-polar buried residues are the same.
	# But for the solvent-accessible residues, 
	# PROFILES3D assumes they are polar and REPIMPS assumes they are non-polar.

	# Calculations for the algorithm to find amphipathic helices :
	# ============================================================
	#
	# all not-buried polar points vs. all not-buried points 
	#				(which gives all not-buried hydrophobic points)
	# if a point is not buried, add it to amphipathic-calc not-buried points count
	# if a point is not buried, see if it is truly polar. 
	#				if so, add it to the amphipathic-calc not-buried polar points.
	#
	# 'solventSurface' solvent surface area = the surface area of this residue when not in a structure/helix
	#				= add up the surface area of each atom in this residue when in this particular amino acid
	# ASA = accessible surface area = ASA for this residue in this particular helix or 3D protein structure
	#		= ( (num points - num buried points) / num points ) * solvent surface area
	#		= ( num not-buried points / num points ) * solvent surface area
	# side_ASA = same as above = side chain ASA = CA.ASA + ASA of each side chain atom, 
	#			in this particular helix or 3D protein structure
	# amphipathic_calc_side_chain_not_buried_polar_ASA = side_ASA * 
	#			( amphipathic_calc_not_buried_is_polar_points count / amphipathic_calc_not_buried_points count) 

	my $resolution = shift;

	$structure = $analyser->{'st'};

	# ////////////local variables//////////////////
	my $genVec; #C++ vector<Vec3f>					# vector<Vec3f> genVec;		//vector template
	my $isBuried = 0;						# bool isBuried = false;

									# Analyser::logFile<<"Creating structure "<<fileName<<"..."<<endl;

	# /////////seting resolution/////////
	my $path = '';
	if ($resolution == 0) {						# if (resolution == 0){
		$path = $path_to_files_Sphere_Points . '/maxvol.3.92.txt';
		$genVec = SpherePoints_generate( $genVec, $path );	# SpherePoints::generate(genVec, currentDir + "/Sphere Points/maxvol.3.92.txt");
	} elsif ($resolution == 1) {					# else if(resolution == 1){	
		$path = $path_to_files_Sphere_Points . '/maxvol.3.162.txt';
		$genVec = SpherePoints_generate( $genVec, $path );	# SpherePoints::generate(genVec, currentDir + "/Sphere Points/maxvol.3.162.txt");
	} elsif ($resolution == 2) {					# else if(resolution == 2){	//default SASA resolution
		$path = $path_to_files_Sphere_Points . '/maxvol.3.252.txt';
		$genVec = SpherePoints_generate( $genVec, $path );	# SpherePoints::generate(genVec, currentDir + "/Sphere Points/maxvol.3.252.txt");
	} elsif ($resolution == 3) {					# else if(resolution == 3){
		$path = $path_to_files_Sphere_Points . '/maxvol.3.392.txt';
		$genVec = SpherePoints_generate( $genVec, $path );	# SpherePoints::generate(genVec, currentDir + "/Sphere Points/maxvol.3.392.txt");
	} elsif ($resolution == 4) {					# else if(resolution == 4){
		$path = $path_to_files_Sphere_Points . '/maxvol.3.482.txt';
		$genVec = SpherePoints_generate( $genVec, $path );	# SpherePoints::generate(genVec, currentDir + "/Sphere Points/maxvol.3.482.txt");
	} elsif ($resolution == 5) {					# else if(resolution == 5){
		$path = $path_to_files_Sphere_Points . '/maxvol.3.1002.txt';
		$genVec = SpherePoints_generate( $genVec, $path );	# SpherePoints::generate(genVec, currentDir + "/Sphere Points/maxvol.3.1002.txt");
	}

	# //SpherePoints::generate(genVec, "C:/Documents and Settings/johann/Desktop/Lawrence/PhD 2003/Program Dev/Valpred dev/spherical code/maxvol.3.252.txt");
	# /////////House-keeping seting resolution etc/////////
	# //if (resolution == 0){
	# //	SASARes = resolution;
	# //	genVec = v.lowRes();	//creates list of generic vectors
	# //}
	# //else if(resolution == 1){	
	# //	SASARes = resolution;
	# //	genVec = v.medRes();
	# //}
	# //else if(resolution == 2){	//default SASA resolution
	# //	SASARes = resolution;
	# //	genVec = v.highRes();
	# //}
	# //else if(resolution == 3){
	# //	SASARes = resolution;
	# //	genVec = v.veryHighRes();
	# //}

	# reset ASA values;
	my $st_dot_allResidues_dot_size = $#{$structure->{'allResidues'}} + 1;
	for (my $r = 0; $r < $st_dot_allResidues_dot_size; $r++ ) {		# for (unsigned int r = 0; r < st.allResidues.size(); r++){
		$structure->{'allResidues'}->[$r]->{'main_ASA'} = 0;		# st.allResidues[r].main_ASA = 0.0;
		$structure->{'allResidues'}->[$r]->{'side_ASA'} = 0;		# st.allResidues[r].side_ASA = 0.0;
		$structure->{'allResidues'}->[$r]->{'all_SC_points'} = 0;	# st.allResidues[r].all_SC_points = 0;
		$structure->{'allResidues'}->[$r]->{'exp_SC_points'} = 0;	# st.allResidues[r].exp_SC_points = 0;
		$structure->{'allResidues'}->[$r]->{'amphipathic_calc_side_chain_not_buried_polar_ASA'} = 0;
		$structure->{'allResidues'}->[$r]->{'amphipathic_calc_not_buried_is_polar_points'} = 0;
		$structure->{'allResidues'}->[$r]->{'amphipathic_calc_not_buried_points'} = 0;
	}
	$structure->{'ASA'} = 0;						# st.ASA = 0.0;

	my $genVec_dot_size = $#{$genVec} + 1;
	my $no_of_points = $genVec_dot_size;					# const unsigned int no_of_points = (int)genVec.size();
										# ProgDlg prog;
										# prog.Create(IDD_PROGRESS);
										# prog.ShowWindow(SW_SHOWNORMAL);
										# prog.setText("Calculating solvent accessible surface area...");

	# ///////calculating SASA///////////

	my $st_dot_allAtoms_dot_size = $#{$structure->{'allAtoms'}} + 1;
	for (my $i = 0; $i < $st_dot_allAtoms_dot_size; $i++ ) {		# for (unsigned int i = 0; i < st.allAtoms.size(); i++){
		my $buried = 0;							# int buried = 0;
		$structure->{'allAtoms'}->[$i]->{'polarPoints'} = 0;		# st.allAtoms[i].polarPoints = 0;
		$structure->{'allAtoms'}->[$i]->{'surfPoints'} = 0;		# st.allAtoms[i].surfPoints = 0;
		$structure->{'allAtoms'}->[$i]->{'amphipathic_calc_not_buried_points'} = 0;
		$structure->{'allAtoms'}->[$i]->{'amphipathic_calc_not_buried_is_polar_points'} = 0;

		my $genVec_dot_size = $#{$genVec} + 1;	
		for (my $v = 0; $v < $genVec_dot_size; $v++ ) {			# for (unsigned int v = 0; v < genVec.size(); v++){
			# // to specify point l(t) = C + bt vector in parametric form
			# // C is the centre of the atom, b is the generic vector
			# // t is the time or radial length
			# //////////////////////////////////////////////////////////
			my $temp1 = Vec3_multiply_number( $genVec->[$v], $structure->{'allAtoms'}->[$i]->{'solventRadii'} );
			my $p = Vec3_add_Vec3( $structure->{'allAtoms'}->[$i]->{'coords'}, $temp1 ); # Vec3f p = st.allAtoms[i].coords + 
													# 	(genVec[v] * st.allAtoms[i].solventRadii); 
			my $st_dot_allAtoms_dot_closeAtoms_dot_size = $#{$structure->{'allAtoms'}->[$i]->{'closeAtoms'}} + 1;
			for (my $j = 0; ($j < $st_dot_allAtoms_dot_closeAtoms_dot_size) && (!$isBuried); $j++ ) { # for(unsigned int j = 0; j < st.allAtoms[i].closeAtoms.size()
														#	&& !isBuried; j++){
				my $term1_atomID = $structure->{'allAtoms'}->[$i]->{'closeAtoms'}->[$j]->{'atomID'};
				my $term1 = $structure->{'allAtoms'}->[$term1_atomID]->{'coords'};
				my $temp2 = Vec3_subtract_Vec3( $p, $term1 );
				my $temp2sqr = Vec3_lengthSqr( $temp2 );
				my $term2_atomID = $structure->{'allAtoms'}->[$i]->{'closeAtoms'}->[$j]->{'atomID'};
				my $term2 = $structure->{'allAtoms'}->[$term2_atomID]->{'solventRadii'};
				my $temp3 = $term2 * $term2;
				if ( $temp2sqr < $temp3 ) {					# if ((p - st.allAtoms[i].closeAtoms[j]->coords).lengthSqr()
												#	< st.allAtoms[i].closeAtoms[j]->solventRadii *
												#	st.allAtoms[i].closeAtoms[j]->solventRadii){
					$buried++;						# buried++;
					$isBuried = 1;						# isBuried = true;
					my $term4_atomID = $structure->{'allAtoms'}->[$i]->{'closeAtoms'}->[$j]->{'atomID'};
					if ( $structure->{'allAtoms'}->[$i]->{'resSeq'} != $structure->{'allAtoms'}->[$term4_atomID]->{'resSeq'} ) { # if( st.allAtoms[i].resSeq !=
												#	st.allAtoms[i].closeAtoms[j]->resSeq ){
						$structure->{'allAtoms'}->[$i]->{'surfPoints'} += 1;	#st.allAtoms[i].surfPoints++;
						my $term5_atomID = $structure->{'allAtoms'}->[$i]->{'closeAtoms'}->[$j]->{'atomID'};
						if ( $structure->{'allAtoms'}->[$term5_atomID]->{'isPolar'} ) { # if (st.allAtoms[i].closeAtoms[j]->isPolar){
							$structure->{'allAtoms'}->[$i]->{'polarPoints'} += 1; # st.allAtoms[i].polarPoints++;
						}
					}
				}
			}
			if (!$isBuried) { 							# if (!isBuried){
				$structure->{'allAtoms'}->[$i]->{'polarPoints'} += 1;		# st.allAtoms[i].polarPoints++;
				$structure->{'allAtoms'}->[$i]->{'surfPoints'} += 1;		# st.allAtoms[i].surfPoints++;
				$structure->{'allAtoms'}->[$i]->{'amphipathic_calc_not_buried_points'} += 1;
				if ( $structure->{'allAtoms'}->[$i]->{'isPolar'} ) {
					# this point is not-buried and is truly polar (not just assumed to be polar by PROFILES3D due to not being buried)
					$structure->{'allAtoms'}->[$i]->{'amphipathic_calc_not_buried_is_polar_points'} += 1;
				}
			}
			$isBuried = 0;								# isBuried = false;
		}
		# st.allAtoms[i].ASA = ( ((float)no_of_points - (float)buried) / (float)no_of_points ) 
		#								* st.allAtoms[i].solventSurface;
		$structure->{'allAtoms'}->[$i]->{'ASA'} = ( ($no_of_points - $buried) / $no_of_points ) * $structure->{'allAtoms'}->[$i]->{'solventSurface'};
	}
	my $atomIndex;										# unsigned int atomIndex;		//copy of index for convenience

	$st_dot_allResidues_dot_size = $#{$structure->{'allResidues'}} + 1;
	for (my $t = 0; $t < $st_dot_allResidues_dot_size; $t++ ) {				# for (unsigned int t = 0; t < st.allResidues.size(); t++){
		my $haha = $structure->{'allResidues'}->[$t];					# Residue haha = st.allResidues[t];
		# //////main chain ASA////
		if ( Residue_checkMainChain( $structure->{'allResidues'}->[$t] ) ) {		# if(st.allResidues[t].checkMainChain()){
			my $temp1 = $structure->{'allResidues'}->[$t]->{'C'};
			$structure->{'allResidues'}->[$t]->{'main_ASA'} += $structure->{'allAtoms'}->[$temp1]->{'ASA'}; # st.allResidues[t].main_ASA += st.allAtoms[st.allResidues[t].C].ASA;
			my $temp2 = $structure->{'allResidues'}->[$t]->{'O'};
			$structure->{'allResidues'}->[$t]->{'main_ASA'} += $structure->{'allAtoms'}->[$temp2]->{'ASA'}; #st.allResidues[t].main_ASA += st.allAtoms[st.allResidues[t].O].ASA;
			my $temp3 = $structure->{'allResidues'}->[$t]->{'N'};
			$structure->{'allResidues'}->[$t]->{'main_ASA'} += $structure->{'allAtoms'}->[$temp3]->{'ASA'}; # st.allResidues[t].main_ASA += st.allAtoms[st.allResidues[t].N].ASA;
			# //alpha carbon is taken as part of a side chain
			my $temp4 = $structure->{'allResidues'}->[$t]->{'CA'};
			$structure->{'allResidues'}->[$t]->{'side_ASA'} += $structure->{'allAtoms'}->[$temp4]->{'ASA'}; # st.allResidues[t].side_ASA += st.allAtoms[st.allResidues[t].CA].ASA;
			$structure->{'allResidues'}->[$t]->{'all_SC_points'} += $structure->{'allAtoms'}->[$temp4]->{'surfPoints'}; # st.allResidues[t].all_SC_points += st.allAtoms[st.allResidues[t].CA].surfPoints;
			$structure->{'allResidues'}->[$t]->{'exp_SC_points'} += $structure->{'allAtoms'}->[$temp4]->{'polarPoints'}; # st.allResidues[t].exp_SC_points += st.allAtoms[st.allResidues[t].CA].polarPoints;
			$structure->{'allResidues'}->[$t]->{'amphipathic_calc_not_buried_points'} += $structure->{'allAtoms'}->[$temp4]->{'amphipathic_calc_not_buried_points'};
			$structure->{'allResidues'}->[$t]->{'amphipathic_calc_not_buried_is_polar_points'} += $structure->{'allAtoms'}->[$temp4]->{'amphipathic_calc_not_buried_is_polar_points'};
		} else {
			# //print error in log file here, missing main chain atoms
			# logFile << "Warning!! - residue "<<st.allResidues[t].name<<st.allResidues[t].resSeq<<" missing main chain atoms. Scores may be incorrect"<<endl;
			print "ERROR: Warning!! - residue " . $structure->{'allResidues'}->[$t]->{'name'} . " missing main chain atoms. Scores may be incorrect\n";
			$structure->{'allResidues'}->[$t]->{'all_SC_points'} = 1;		# st.allResidues[t].all_SC_points = 1;
			$structure->{'allResidues'}->[$t]->{'exp_SC_points'} = 1;		# st.allResidues[t].exp_SC_points = 1;
		}
		$haha = $structure->{'allResidues'}->[$t];					# haha = st.allResidues[t];
		# ///////handling side chain atoms//////////
		my $st_dot_allResidues_dot_sideChains_dot_size = $#{$structure->{'allResidues'}->[$t]->{'sideChains'}} + 1;
		for (my $s = 0; $s < $st_dot_allResidues_dot_sideChains_dot_size; $s++ ) {	# for (unsigned int s = 0; s < st.allResidues[t].sideChains.size(); s++){
			my $atomIndex = $structure->{'allResidues'}->[$t]->{'sideChains'}->[$s]; # atomIndex = st.allResidues[t].sideChains[s];
			$structure->{'allResidues'}->[$t]->{'side_ASA'} += $structure->{'allAtoms'}->[$atomIndex]->{'ASA'}; # st.allResidues[t].side_ASA += st.allAtoms[atomIndex].ASA;
			$structure->{'allResidues'}->[$t]->{'exp_SC_points'} += $structure->{'allAtoms'}->[$atomIndex]->{'polarPoints'}; # st.allResidues[t].exp_SC_points += st.allAtoms[atomIndex].polarPoints;
			$structure->{'allResidues'}->[$t]->{'all_SC_points'} += $structure->{'allAtoms'}->[$atomIndex]->{'surfPoints'}; # st.allResidues[t].all_SC_points += st.allAtoms[atomIndex].surfPoints;
			$structure->{'allResidues'}->[$t]->{'amphipathic_calc_not_buried_points'} += $structure->{'allAtoms'}->[$atomIndex]->{'amphipathic_calc_not_buried_points'};
			$structure->{'allResidues'}->[$t]->{'amphipathic_calc_not_buried_is_polar_points'} += $structure->{'allAtoms'}->[$atomIndex]->{'amphipathic_calc_not_buried_is_polar_points'};
		}
		if ($structure->{'allResidues'}->[$t]->{'aaCode'} == UNK) {			# if (st.allResidues[t].aaCode == UNK){
			$structure->{'allResidues'}->[$t]->{'areaBuried'} = 0;			# st.allResidues[t].areaBuried = 0.0f;
			$structure->{'allResidues'}->[$t]->{'fractionPolar'} = 0;		# st.allResidues[t].fractionPolar = 0.0f;
			$structure->{'allResidues'}->[$t]->{'repimpsFP'} = 0;			# st.allResidues[t].repimpsFP = 0.0f;
			# //print log file error here, unknown residue type
			# logFile << "Warning!! Unknown residue type, residue: "<<st.allResidues[t].resSeq<<" chain: "<<st.allResidues[t].chainID<<endl;
			print "ERROR: Warning!! Unknown residue type, residue: " . $structure->{'allResidues'}->[$t]->{'resSeq'} . " chain: " . $structure->{'allResidues'}->[$t]->{'chainID'} . "\n";
		}
		else{
			# st.allResidues[t].areaBuried = SCGlyXGly[st.allResidues[t].aaCode] - 
			#	st.allResidues[t].side_ASA ;
			my $temp4 = $structure->{'allResidues'}->[$t]->{'aaCode'};
			$structure->{'allResidues'}->[$t]->{'areaBuried'} = $SCGlyXGly->[$temp4] - $structure->{'allResidues'}->[$t]->{'side_ASA'};
			# st.allResidues[t].fractionPolar = (float)st.allResidues[t].exp_SC_points / 
			#	(float)st.allResidues[t].all_SC_points  ; 
			$structure->{'allResidues'}->[$t]->{'fractionPolar'} = $structure->{'allResidues'}->[$t]->{'exp_SC_points'} / $structure->{'allResidues'}->[$t]->{'all_SC_points'};
			# st.allResidues[t].repimpsFP =	st.allResidues[t].fractionPolar - 
			#	st.allResidues[t].side_ASA / 
			#	SCGlyXGly[st.allResidues[t].aaCode];

			$temp4 = $structure->{'allResidues'}->[$t]->{'aaCode'};
			$structure->{'allResidues'}->[$t]->{'repimpsFP'} = $structure->{'allResidues'}->[$t]->{'fractionPolar'} - $structure->{'allResidues'}->[$t]->{'side_ASA'} / $SCGlyXGly->[$temp4];
		}

		# amphipathic_calc_side_chain_not_buried_polar_ASA = side_ASA * 
		#			( amphipathic_calc_not_buried_is_polar_points count / amphipathic_calc_not_buried_points count) 
		if ($structure->{'allResidues'}->[$t]->{'amphipathic_calc_not_buried_points'} == 0) {
			$structure->{'allResidues'}->[$t]->{'amphipathic_calc_side_chain_not_buried_polar_ASA'} = 0;
		} else {
			$structure->{'allResidues'}->[$t]->{'amphipathic_calc_side_chain_not_buried_polar_ASA'} = 
				$structure->{'allResidues'}->[$t]->{'side_ASA'} *
					($structure->{'allResidues'}->[$t]->{'amphipathic_calc_not_buried_is_polar_points'} /
					$structure->{'allResidues'}->[$t]->{'amphipathic_calc_not_buried_points'});
		}

		$haha = $structure->{'allResidues'}->[$t];					# haha = st.allResidues[t];
		$structure->{'ASA'} += $structure->{'allResidues'}->[$t]->{'main_ASA'};		# st.ASA += st.allResidues[t].main_ASA;
		$structure->{'ASA'} += $structure->{'allResidues'}->[$t]->{'side_ASA'};		# st.ASA += st.allResidues[t].side_ASA;
	}
												# prog.DestroyWindow();
	$analyser->{'st'} = $structure;

	return;
}



sub Analyser_setOverallScores { #C++ void Analyser::setOverallScores(void)

	$structure = $analyser->{'st'};

	$structure->{'scoreTotal'} = 0;								# st.scoreTotal = 0.0f;
	my $st_dot_allResidues_dot_size = $#{$structure->{'allResidues'}} + 1;
	for ( my $i = 0; $i < $st_dot_allResidues_dot_size; $i++ ) {				# for(unsigned int i = 0; i < st.allResidues.size(); i++){
		if ( $structure->{'allResidues'}->[$i]->{'isTransMem'} ) {	# if(st.allResidues[i].isTransMem){
			# st.allResidues[i].scores[2] = st.allResidues[i].scores[1];
			$structure->{'allResidues'}->[$i]->{'scores'}->[2] = $structure->{'allResidues'}->[$i]->{'scores'}->[1];
		} else {
			# st.allResidues[i].scores[2] = st.allResidues[i].scores[0];
			$structure->{'allResidues'}->[$i]->{'scores'}->[2] = $structure->{'allResidues'}->[$i]->{'scores'}->[0];
		}
		$structure->{'scoreTotal'} += $structure->{'allResidues'}->[$i]->{'scores'}->[2]; # st.scoreTotal+=st.allResidues[i].scores[2];
	}
	$structure->{'aveResScore'} = $structure->{'scoreTotal'} / $st_dot_allResidues_dot_size; # st.aveResScore = st.scoreTotal/st.allResidues.size();

	$analyser->{'st'} = $structure;
	return;
}



sub Analyser_xyzToIndex { #C++ int Analyser::xyzToIndex(int x, int y, int z){

	my $x = shift;
	my $y = shift;
	my $z = shift;

	if (($x < 0) || ($y < 0) || ($z < 0) || # if (x < 0 || y < 0 || z < 0 ||
		($x >= $analyser->{'noOfXCubes'}) || ($y >= $analyser->{'noOfYCubes'}) || ($z >= $analyser->{'noOfZCubes'})) { # x >= noOfXCubes || y >= noOfYCubes ||z >= noOfZCubes){
		return -1;			# return -1;
	}
	# return z*noOfXCubes*noOfYCubes + noOfXCubes*y + x;
	my $rtn_var = $z * $analyser->{'noOfXCubes'} * $analyser->{'noOfYCubes'} + $analyser->{'noOfXCubes'} * $y + $x;
	return $rtn_var;
}



sub Atom_default_values {

	my $atom;

	# //initialise Wan Der Waals radii with default values
	$atom->{ 'CARadii' } = 2;			# float Atom::CARadii = 2.0f;
	$atom->{ 'CRadii' } = 1.5;			# float Atom::CRadii = 1.5f;
	$atom->{ 'ORadii' } = 1.4;			# float Atom::ORadii = 1.4f;
	$atom->{ 'NRadii' } = 1.5;			# float Atom::NRadii = 1.5f;
	$atom->{ 'SolventRadii' } = 1.4;		# float Atom::SolventRadii = 1.4f;

	$atom->{ 'S_ALL' } = 1.85;			# float Atom::S_ALL = 1.85f;
	$atom->{ 'C_ARO' } = 1.85;			# float Atom::C_ARO = 1.85f;

	# //Initialise Solvent Sphere Surface with Default Values 4*pi*r^2
	my $pi4 = 4 * PI;
	# float Atom::CASurface = 4.0f * pi * (CARadii + SolventRadii) * (CARadii + SolventRadii);
	$atom->{ 'CASurface' } = $pi4 * ($atom->{ 'CARadii' } + $atom->{ 'SolventRadii' }) * ($atom->{ 'CARadii' } + $atom->{ 'SolventRadii' });
	# float Atom::CSurface = 4.0f * pi * (CRadii + SolventRadii) * (CRadii + SolventRadii);
	$atom->{ 'CSurface' } = $pi4 * ($atom->{ 'CRadii' } + $atom->{ 'SolventRadii' }) * ($atom->{ 'CRadii' } + $atom->{ 'SolventRadii' });
	# float Atom::OSurface = 4.0f * pi * (ORadii + SolventRadii) * (ORadii + SolventRadii);
	$atom->{ 'OSurface' } = $pi4 * ($atom->{ 'ORadii' } + $atom->{ 'SolventRadii' }) * ($atom->{ 'ORadii' } + $atom->{ 'SolventRadii' });
	# float Atom::NSurface = 4.0f * pi * (NRadii + SolventRadii) * (NRadii + SolventRadii);
	$atom->{ 'NSurface' } = $pi4 * ($atom->{ 'NRadii' } + $atom->{ 'SolventRadii' }) * ($atom->{ 'NRadii' } + $atom->{ 'SolventRadii' });
	# float Atom::C_ARO_Surface = 4.0f * pi * (C_ARO + SolventRadii) * (C_ARO + SolventRadii);
	$atom->{ 'C_ARO_Surface' } = $pi4 * ($atom->{ 'C_ARO' } + $atom->{ 'SolventRadii' }) * ($atom->{ 'C_ARO' } + $atom->{ 'SolventRadii' });
	# float Atom::S_ALL_Surface = 4.0f * pi * (S_ALL + SolventRadii) * (S_ALL + SolventRadii);
	$atom->{ 'S_ALL_Surface' } = $pi4 * ($atom->{ 'S_ALL' } + $atom->{ 'SolventRadii' }) * ($atom->{ 'S_ALL' } + $atom->{ 'SolventRadii' });

	return $atom;
}



sub Atom_new { #C++ Atom::Atom()

	my $atom;

	$atom->{ 'polarPoints' } = 0;			# polarPoints = 0;
	$atom->{ 'surfPoints' } = 0;			# surfPoints = 0;
	$atom->{ 'amphipathic_calc_not_buried_points' } = 0;
	$atom->{ 'amphipathic_calc_not_buried_is_polar_points' } = 0;
	$atom->{ 'coords' }->{ 'x' } = 0;		# coords.x = 0.0f;
	$atom->{ 'coords' }->{ 'y' } = 0;
	$atom->{ 'coords' }->{ 'z' } = 0;
	$atom->{ 'isHet' } = 0;				# isHet = false;
	$atom->{ 'name' } = '';
	$atom->{ 'aType' } = '';
	$atom->{ 'aaCode' } = '';

	# //initialise Wan Der Waals radii with default values
	$atom->{ 'CARadii' } = $atom_default_values->{ 'CARadii' };		# float Atom::CARadii = 2.0f;
	$atom->{ 'CRadii' } = $atom_default_values->{ 'CRadii' };		# float Atom::CRadii = 1.5f;
	$atom->{ 'ORadii' } = $atom_default_values->{ 'ORadii' };		# float Atom::ORadii = 1.4f;
	$atom->{ 'NRadii' } = $atom_default_values->{ 'NRadii' };		# float Atom::NRadii = 1.5f;
	$atom->{ 'SolventRadii' } = $atom_default_values->{ 'SolventRadii' };	# float Atom::SolventRadii = 1.4f;

	$atom->{ 'S_ALL' } = $atom_default_values->{ 'S_ALL' };			# float Atom::S_ALL = 1.85f;
	$atom->{ 'C_ARO' } = $atom_default_values->{ 'C_ARO' };			# float Atom::C_ARO = 1.85f;

	# //Initialise Solvent Sphere Surface with Default Values 4*pi*r^2
	# float Atom::CASurface = 4.0f * pi * (CARadii + SolventRadii) * (CARadii + SolventRadii);
	$atom->{ 'CASurface' } = $atom_default_values->{ 'CASurface' };
	# float Atom::CSurface = 4.0f * pi * (CRadii + SolventRadii) * (CRadii + SolventRadii);
	$atom->{ 'CSurface' } = $atom_default_values->{ 'CSurface' };
	# float Atom::OSurface = 4.0f * pi * (ORadii + SolventRadii) * (ORadii + SolventRadii);
	$atom->{ 'OSurface' } = $atom_default_values->{ 'OSurface' };
	# float Atom::NSurface = 4.0f * pi * (NRadii + SolventRadii) * (NRadii + SolventRadii);
	$atom->{ 'NSurface' } = $atom_default_values->{ 'NSurface' };
	# float Atom::C_ARO_Surface = 4.0f * pi * (C_ARO + SolventRadii) * (C_ARO + SolventRadii);
	$atom->{ 'C_ARO_Surface' } = $atom_default_values->{ 'C_ARO_Surface' };
	# float Atom::S_ALL_Surface = 4.0f * pi * (S_ALL + SolventRadii) * (S_ALL + SolventRadii);
	$atom->{ 'S_ALL_Surface' } = $atom_default_values->{ 'S_ALL_Surface' };

	return $atom; #C++ Atom
}



# //sets the atom type and the corresponding radii and surface areas
sub Atom_setAtomType { #C++ void Atom::setAtomType(Residue &r, unsigned int i) 

	my $atom = shift; #C++ Atom
	my $r = shift; #C++ Residue
	my $i = shift; #C++ unsigned int

	if ($atom->{'name'} eq '') {				# if (name.empty()){
		return;						# return;
	}

	# ////main chain atoms///////
	if ($atom->{'name'} eq 'N') {				# if ( name == "N" ){
		$atom->{'aType'} = N_All;			# aType = N_All;
		$r->{'N'} = $i;					# r.N = i;

	} elsif ($atom->{'name'} eq 'CA') {			# else if ( name == "CA" ){
		$atom->{'aType'} = C_Ali;			# aType = C_Ali;
		$r->{'CA'} = $i;				# r.CA = i;

	} elsif ($atom->{'name'} eq 'C') {			# else if ( name == "C" ){
		$atom->{'aType'} = C_Car;			# aType = C_Car;
		$r->{'C'} = $i;					# r.C = i;

	} elsif ($atom->{'name'} eq 'O') {			# else if ( name == "O" ){
		$atom->{'aType'} = O_All;			# aType = O_All;
		$r->{'O'} = $i;					# r.O = i;

	# ///////side chain atoms//////////
	} else {						# else {
		
		my $push_back_index = $#{$r->{'sideChains'}} + 1;
		$r->{'sideChains'}->[$push_back_index] = $i; 	# r.sideChains.push_back(i);
		if ( substr($atom->{'name'},0,1) eq 'O') {	# if (name[0] == 'O'){
			$atom->{'aType'} = O_All;		# aType = O_All;
		} elsif ( substr($atom->{'name'},0,1) eq 'N') { # else if (name[0] == 'N'){
			$atom->{'aType'} = N_All;		# aType = N_All;
		} else {					# else{
			if ($atom->{'aType'} eq '') {		# if (aType == NULL){
				$atom->{'aType'} = C_Car; 	# aType = C_Car;
			}
			Atom_setSCType( $atom, $r );		# setSCType(r);
		}

	}
	Atom_setRadii( $atom );					# setRadii();
}



sub Atom_setRadii { #C++ void Atom::setRadii()

	my $atom = shift; #C++ Atom
							# switch (aType){
	if ($atom->{'aType'} eq C_Ali) {		# case C_Ali:
		$atom->{'solventRadii'} = $atom->{'CARadii'} + $atom->{'SolventRadii'}; # solventRadii = CARadii + SolventRadii;
		$atom->{'solventSurface'} = $atom->{'CASurface'}; # solventSurface = CASurface; # //4.pi.r^2
		$atom->{'isPolar'} = 0;			# isPolar = false;
							# break;
	} elsif ($atom->{'aType'} eq C_Car) {		# case C_Car:
		$atom->{'solventRadii'} = $atom->{'CRadii'} + $atom->{'SolventRadii'}; # solventRadii = CRadii + SolventRadii;
		$atom->{'solventSurface'} = $atom->{'CSurface'}; # solventSurface = CSurface; # //4.pi.r^2
		$atom->{'isPolar'} = 0;			# isPolar = false;
							# break; # //edited
	} elsif ($atom->{'aType'} eq C_Aro) {		# case C_Aro:
		$atom->{'solventRadii'} = $atom->{'C_ARO'} + $atom->{'SolventRadii'}; # solventRadii = C_ARO + SolventRadii;
		$atom->{'solventSurface'} = $atom->{'C_ARO_Surface'}; # solventSurface = C_ARO_Surface; # //4.pi.r^2
		$atom->{'isPolar'} = 0;			# isPolar = false;
							# break;
	} elsif ($atom->{'aType'} eq N_All) {		# case N_All:
		$atom->{'solventRadii'} = $atom->{'NRadii'} + $atom->{'SolventRadii'}; # solventRadii = NRadii + SolventRadii;
		$atom->{'solventSurface'} = $atom->{'NSurface'}; # solventSurface = NSurface; # //4.pi.r^2
		$atom->{'isPolar'} = 1;			# isPolar = true;
							# break;
	} elsif ($atom->{'aType'} eq O_All) {		# case O_All:
		$atom->{'solventRadii'} = $atom->{'ORadii'} + $atom->{'SolventRadii'}; # solventRadii = ORadii + SolventRadii;
		$atom->{'solventSurface'} = $atom->{'OSurface'}; # solventSurface = OSurface; # //4.pi.r^2
		$atom->{'isPolar'} = 1;			# isPolar = true;
							# break;
	} elsif ($atom->{'aType'} eq S_Oxy) {		# case S_Oxy:
		$atom->{'solventRadii'} = $atom->{'S_ALL'} + $atom->{'SolventRadii'}; # solventRadii = S_ALL + SolventRadii;
		$atom->{'solventSurface'} = $atom->{'S_ALL_Surface'}; # solventSurface = S_ALL_Surface; # //4.pi.r^2
		$atom->{'isPolar'} = 0;			# isPolar = false;
							# break;
	} elsif ($atom->{'aType'} eq S_Red) {		# case S_Red:
		$atom->{'solventRadii'} = $atom->{'S_ALL'} + $atom->{'SolventRadii'}; # solventRadii = S_ALL + SolventRadii;
		$atom->{'solventSurface'} = $atom->{'S_ALL_Surface'}; # solventSurface = S_ALL_Surface; # //4.pi.r^2
		$atom->{'isPolar'} = 0;			# isPolar = false;
							# break; # //edited
	} else {					# default:
		$atom->{'solventRadii'} = 3.4;		# solventRadii = 3.4f;
		$atom->{'solventSurface'} = 4 * PI * $atom->{'solventRadii'} * $atom->{'solventRadii'}; # solventSurface = 4.0f*pi*solventRadii*solventRadii;
							# break;
	}
}



sub Atom_setSCType { #C++ void Atom::setSCType(Residue &r)

	my $atom = shift; #C++ Atom
	my $r = shift; #C++ Residue
							# switch(r.aaCode){
							# /* case ALA:
							# aType = C_Ali;
							# break;*/
	if ($r->{'aaCode'} eq ARG) {			# case ARG:
		if ($atom->{'name'} eq 'CZ') {		# if(name == "CZ"){
			$atom->{'aType'} = C_Aro;	# aType = C_Aro;
		} else {				# else{
			$atom->{'aType'} = C_Ali;	# aType = C_Ali;
		}					# break;
	} elsif ($r->{'aaCode'} eq ASN) {		# case ASN:
		if ($atom->{'name'} eq 'CG') {		# if (name == "CG"){
			$atom->{'aType'} = C_Car;	# aType = C_Car;
		} else {				# else{
			$atom->{'aType'} = C_Ali;	# aType = C_Ali;
		}					# break;
	} elsif ($r->{'aaCode'} eq ASP) {		# case ASP:
		if ($atom->{'name'} eq 'CG') {		# if (name == "CG"){
			$atom->{'aType'} = C_Car;	# aType = C_Car;
		} else {				# else{
			$atom->{'aType'} = C_Ali;	# aType = C_Ali;
		}					# break;
	} elsif ($r->{'aaCode'} eq CYS) {		# case CYS:
		if ($atom->{'name'} eq 'SG') {		# if (name == "SG"){
			$atom->{'aType'} = S_Red;	# aType = S_Red;
		} else {				# else{
			$atom->{'aType'} = C_Ali;	# aType = C_Ali;
		}					# break;
	} elsif ($r->{'aaCode'} eq GLN) {		# case GLN:
		if ($atom->{'name'} eq 'CD') {		# if (name == "CD"){
			$atom->{'aType'} = C_Car;	# aType = C_Car;
		} else {				# else{
			$atom->{'aType'} = C_Ali;	# aType = C_Ali;
		}					# break;
	} elsif ($r->{'aaCode'} eq GLU) {		# case GLU:
		if ($atom->{'name'} eq 'CD') {		# if (name == "CD"){
			$atom->{'aType'} = C_Car;	# aType = C_Car;
		} else {				# else{
			$atom->{'aType'} = C_Ali;	# aType = C_Ali;
		}					# break;
	} elsif ($r->{'aaCode'} eq GLY) {		# case GLY:
		my $break = 1;				# break;
	} elsif ($r->{'aaCode'} eq HIS) {		# case HIS:
		if ($atom->{'name'} eq 'CB') {		# if (name == "CB"){
			$atom->{'aType'} = C_Ali;	# aType = C_Ali;
		} else {				# else{
			$atom->{'aType'} = C_Aro;	# aType = C_Aro;
		}					# break;
	} elsif ($r->{'aaCode'} eq ILE) {		# case ILE:
		$atom->{'aType'} = C_Ali;		# aType = C_Ali;
							# break;
	} elsif ($r->{'aaCode'} eq LEU) {		# case LEU:
		$atom->{'aType'} = C_Ali;		# aType = C_Ali;
							# break;
	} elsif ($r->{'aaCode'} eq LYS) {		# case LYS:
		$atom->{'aType'} = C_Ali;		# aType = C_Ali;
							# break;
	} elsif ($r->{'aaCode'} eq MET) {		# case MET:
		if ($atom->{'name'} eq 'SD') {		# if (name == "SD"){
			$atom->{'aType'} = S_Oxy;	# aType = S_Oxy;
		} else {				# else{
			$atom->{'aType'} = C_Ali;	# aType = C_Ali;
		}					# break;
	} elsif ($r->{'aaCode'} eq PHE) {		# case PHE:
		if ($atom->{'name'} eq 'CB') {		# if (name == "CB"){
			$atom->{'aType'} = C_Ali;	# aType = C_Ali;
		} else {				# else{
			$atom->{'aType'} = C_Aro;	# aType = C_Aro;
		}					# break;
	} elsif ($r->{'aaCode'} eq PRO) {		# case PRO:
		$atom->{'aType'} = C_Ali;		# aType = C_Ali;
							# break
	} elsif ($r->{'aaCode'} eq SER) {		# case SER:
		$atom->{'aType'} = C_Ali;		# aType = C_Ali;
							# break
	} elsif ($r->{'aaCode'} eq THR) {		# case THR:
		$atom->{'aType'} = C_Ali;		# aType = C_Ali;
							# break
	} elsif ($r->{'aaCode'} eq TRP) {		# case TRP:
		if ($atom->{'name'} eq 'CB') {		# if (name == "CB"){
			$atom->{'aType'} = C_Ali;	# aType = C_Ali;
		} else {				# else{
			$atom->{'aType'} = C_Aro;	# aType = C_Aro;
		}					# break;
	} elsif ($r->{'aaCode'} eq TYR) {		# case TYR:
		if ($atom->{'name'} eq 'CB') {		# if (name == "CB"){
			$atom->{'aType'} = C_Ali;	# aType = C_Ali;
		} else {				# else{
			$atom->{'aType'} = C_Aro;	# aType = C_Aro;
		}					# break;
	} elsif ($r->{'aaCode'} eq VAL) {		# case VAL:
		$atom->{'aType'} = C_Ali;		# aType = C_Ali;
							# break
	} else {					# default:
		$atom->{'aType'} = C_Car;		# aType = C_Car;
							# break
	}
}



sub CChartView_OnDraw { #C++ void CChartView::OnDraw(CDC* pDC)

	# The $chart_coords string that this subroutine returns, 
	# contains the scores and TMS data information, 
	# for the lines that were specified in the input options (default is PROFILES3D, REPIMPS, not BOTH, and TMS),
	# ready for a graph to be created from it.
	# The data is all in one long line, so that it can be passed as a hidden field in an HTML form.
	# It looks like this :
	# OPTIONS,0,1,1,10,blue,Profiles 3D,:::::OPTIONS,1,1,1,10,red,REPIMPS,:::::OPTIONS,2,0,1,10,black,Overall,:::::OPTIONS,3,0,3,10,black,Transmembrane Segments,:::::NUM_RESIDUES,160,:::::NUM_CHAINS,1,:::::CHAIN,0,160,0,:::::LINE,0,:::::DATA,-0.432,-0.333333333333333,-0.195714285714286,-0.22875,-0.254444444444444,-0.185,-0.146,-0.197,-0.264,-0.167,-0.167,-0.229,-0.338,-0.229,-0.234,-0.261,-0.147,-0.085,0.021,-0.032,0.008,-0.055,-0.06,-0.079,0.016,0.062,0.043,0.071,0.00699999999999998,-0.072,-0.112,-0.14,-0.059,-0.109,-0.278,-0.478,-0.568,-0.568,-0.485,-0.344,-0.353,-0.262,-0.248,-0.304,-0.225,-0.044,0.062,-0.033,-0.033,-0.00500000000000003,0.025,0.134,0.134,0.225,0.134,0.031,-0.049,-0.053,-0.141,-0.129,-0.06,-0.169,-0.322,-0.46,-0.307,-0.373,-0.396,-0.302,-0.233,-0.309,-0.373,-0.452,-0.28,-0.226,-0.179,-0.01,-0.104,-0.198,-0.267,-0.372,-0.372,-0.293,-0.465,-0.347,-0.461,-0.551,-0.384,-0.285,-0.385,-0.385,-0.49,-0.534,-0.487,-0.506,-0.392,-0.582,-0.632,-0.66,-0.716,-0.528,-0.54,-0.575,-0.533,-0.533,-0.552,-0.362,-0.298,-0.334,-0.262,-0.46,-0.385,-0.197,-0.357,-0.381,-0.363,-0.297,-0.392,-0.423,-0.251,-0.099,-0.172,-0.36,-0.124,-0.127,-0.227,-0.227,-0.11,-0.015,-0.034,0.00900000000000003,0.164,0.333,0.222,0.192,0.286,0.246,0.219,0.124,0.124,0.018,0.015,-0.00900000000000001,-0.005,-0.127,-0.154,-0.185,-0.27,-0.262,-0.357,-0.338,-0.401,-0.392,-0.297,-0.13,-0.179,-0.148,-0.113333333333333,-0.07375,-0.0114285714285714,0.0316666666666667,:::::LINE,1,:::::DATA,0.14,0.195,0.101428571428571,0.25125,0.367777777777778,0.407,0.101,0.387,0.617,0.354,0.354,0.437,0.613,0.437,0.127,0.017,0.151,0.0680000000000001,0.07,0.337,0.173,0.154,-0.156,-0.034,0.222,0.21,0.332,0.361,0.324,0.268,0.432,0.451,0.705,0.74,0.738,0.914,0.968,0.968,0.883,0.856,0.753,0.753,0.755,0.718,0.774,0.72,0.375,0.119,0.119,0.148,-0.0350000000000001,-0.211,-0.211,-0.423,-0.423,-0.388,-0.134,0.073,-0.037,-0.293,-0.061,0.115,0.15,0.399,0.316,0.279,0.025,0.037,0.269,0.488,0.451,0.395,0.238,0.154,0.061,-0.041,0.304,0.292,0.06,0.095,0.095,0.151,0.308,0.235,0.101,0.259,0.203,0.252,0.482,0.482,0.517,0.513,0.476,0.598,0.732,0.713,0.678,0.649,0.703,0.583,0.639,0.587,0.552,0.552,0.674,0.693,0.73,0.722,0.705,0.879,0.823,0.703,0.484,0.193,-0.021,-0.366,-0.622,-0.841,-0.998,-1.16,-1.106,-0.986,-0.986,-0.805,-0.62,-0.62,-0.577,-0.321,-0.199,-0.38,-0.723,-0.721,-0.514,-0.441,-0.668,-0.414,-0.238,-0.494,-0.494,-0.149,0.032,-0.259,-0.466,-0.377,-0.225,-0.444,-0.353,-0.388,-0.644,-0.648,-0.484,-0.297,-0.145,-0.449,-0.601,-0.601,-0.812222222222222,-0.645,-0.48,-0.77,:::::LINE,3,:::::DATA,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,:::::TITLE,,,,,>KCSA_STRLI Voltage-gated potassium channel:::::
	# OPTIONS,0,1,1,10,blue,Profiles 3D,:::::OPTIONS,1,1,1,10,red,REPIMPS,:::::OPTIONS,2,0,1,10,black,Overall,:::::OPTIONS,3,0,3,10,black,Transmembrane Segments,:::::NUM_RESIDUES,160,:::::NUM_CHAINS,1,:::::CHAIN,0,160,0,:::::LINE,0,:::::DATA,-0.4,-0.3,:::::LINE,1,:::::DATA,0.14,0.195,0.1,:::::LINE,3,:::::DATA,0,1,1,:::::TITLE,,,,,>KCSA_STRLI Voltage-gated potassium channel:::::

	# The $multiplegraphs_data global variable that this subroutine fills in
	# contains the scores and TMS data informations for the lines PROFILES3D, REPIMPS and TMS,
	# in the same format as the $chart_coords variable that is passed from the valpred-calculations cgi-output to the graph-drawing cgi-output,
	# except that the ::::: are changed to ;;;;;

	my $chart_coords = '';
	my $output_data = '';
	$multiplegraphs_filename = '';

	$chart_coords .= 'GRAPH_RESIDUE_WIDTH,' . $program_options->{'graph_residue_width'} . ':::::';
	$chart_coords .= 'GRAPH_HEIGHT,' . $program_options->{'graph_height'} . ':::::';
	$chart_coords .= 'VALPRED1_OR_VALPRED2,' . $TMPred2D->{'valpred1_or_valpred2'} . ':::::';

	# my @option_name = ('line_id','showGraph','weight','mov_avg','colour0','colour1','colour2','colour_name','legend_name','y_value','smooth');
	for ( my $ln = 0; $ln <= $program_options->{'option_max_index'}; $ln++ ) {
		if ($program_options->{'lnOpt'}->[$ln]->{'showGraph'} == 1) {
			$chart_coords .= 'OPTIONS,' . $ln . ',';
			for ( my $n = 0; $n <= $#{$program_options->{'option_name'}}; $n++ ) {
				my $n_value = $program_options->{'option_name'}->[$n];
				$chart_coords .= $program_options->{'lnOpt'}->[$ln]->{$n_value} . ',';
			}
			$chart_coords .= ':::::';
		}
	}

	if (defined ($input_graph_points_file)) {

		my $chart_coords = $multiplegraphs_data;
		$chart_coords =~ s/\;\;\;\;\;/\:\:\:\:\:/g;

	} else {

		$structure = $analyser->{'st'};
		my $chains = $structure->{'chains'};

		# calculate total number of residues to be graphed
		my $noOfRes = 0;
		my $chains_dot_size = $#{$chains} + 1;
		for ( my $r = 0; $r < $chains_dot_size; $r++ ) {
			my $getResChain = Structure_getResChain_Range( $chains->[$r] );
			my $getResChain_dot_size = $#{$getResChain} + 1;
			$noOfRes += $getResChain_dot_size;
		}
		$chart_coords .= 'NUM_RESIDUES,' . $noOfRes . ',:::::';

		my @graph_data;
		my @graph_data_x_axis;
		for ( my  $d = 0; $d < $noOfRes; $d++ ) {
			$graph_data_x_axis[$d] = $d;
		}
		$graph_data[0] = \@graph_data_x_axis;
		my $graph_data_line_index = 0;
		my @graph_colors;
		my @graph_legend;
		my $graph_colors_index = -1;
		my $graph_legend_index = -1;
		my @chain_length;
		my @chain_start_position;

		# Build up the contents of $chart_coords, which will be passed to the actual graph-drawing subroutine.
		# Here is an example of the contents of $chart_coords :
		# OPTIONS,0,1,1,10,blue,Profiles 3D,:::::
		# OPTIONS,1,1,1,10,red,REPIMPS,:::::
		# OPTIONS,2,0,1,10,black,Both,:::::
		# OPTIONS,3,1,1,10,black,Transmembrane Segments,:::::
		# NUM_RESIDUES,160,:::::
		# NUM_CHAINS,2,:::::
		# CHAIN,0,80,0,:::::
		# LINE,0,:::::DATA,-0.432,-0.333333333333333,-0.195714285714286,-0.22875,-0.254444444444444,-0.185,-0.146,-0.197,-0.264,-0.167,-0.167,-0.229,-0.338,-0.229,-0.234,-0.261,-0.147,-0.085,0.021,-0.032,0.008,-0.055,-0.06,-0.079,0.016,0.062,0.043,0.071,0.00699999999999998,-0.072,-0.112,-0.14,-0.059,-0.109,-0.278,-0.478,-0.568,-0.568,-0.485,-0.344,-0.353,-0.262,-0.248,-0.304,-0.225,-0.044,0.062,-0.033,-0.033,-0.00500000000000003,0.025,0.134,0.134,0.225,0.134,0.031,-0.049,-0.053,-0.141,-0.129,-0.06,-0.169,-0.322,-0.46,-0.307,-0.373,-0.396,-0.302,-0.233,-0.309,-0.373,-0.452,-0.28,-0.226,-0.179,-0.01,0.0366666666666667,-0.00750000000000001,-0.0714285714285714,-0.05,:::::
		# LINE,1,:::::DATA,0.14,0.195,0.101428571428571,0.25125,0.367777777777778,0.407,0.101,0.387,0.617,0.354,0.354,0.437,0.613,0.437,0.127,0.017,0.151,0.0680000000000001,0.07,0.337,0.173,0.154,-0.156,-0.034,0.222,0.21,0.332,0.361,0.324,0.268,0.432,0.451,0.705,0.74,0.738,0.914,0.968,0.968,0.883,0.856,0.753,0.753,0.755,0.718,0.774,0.72,0.375,0.119,0.119,0.148,-0.0350000000000001,-0.211,-0.211,-0.423,-0.423,-0.388,-0.134,0.073,-0.037,-0.293,-0.061,0.115,0.15,0.399,0.316,0.279,0.025,0.037,0.269,0.488,0.451,0.395,0.238,0.154,0.061,-0.041,0.193333333333333,0.16875,0.0842857142857143,0.0333333333333333,:::::
		# LINE,3,:::::DATA,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,:::::
		# CHAIN,1,80,80,:::::
		# LINE,0,:::::DATA,-0.204,-0.225,-0.371428571428571,-0.1875,-0.264444444444444,-0.271,-0.166,-0.012,0.00400000000000001,0.00400000000000003,0.076,-0.021,0.095,-0.00099999999999999,0.197,0.0670000000000001,0.084,0.104,-0.139,-0.109,-0.303,-0.0929999999999999,-0.0759999999999999,-0.076,-0.172,-0.0419999999999999,-0.036,-0.062,-0.024,-0.263,-0.158,-0.128,-0.371,-0.444,-0.515,-0.541,-0.643,-0.739,-0.504,-0.496,-0.616,-0.646,-0.499,-0.561,-0.445,-0.445,-0.429,-0.327,-0.423,-0.298,-0.228,-0.294,-0.346,-0.17,-0.317,-0.25,-0.123,-0.225,-0.225,-0.262,-0.324,-0.397,-0.345,-0.529,-0.407,-0.503,-0.591,-0.562,-0.664,-0.721,-0.633,-0.523,-0.412,-0.242,-0.317,-0.27,-0.263333333333333,-0.2225,-0.128571428571429,0,:::::
		# LINE,1,:::::DATA,0.744,0.796666666666667,0.82,0.66875,0.333333333333333,0.406,0.4,0.397,0.438,0.438,0.483,0.468,0.472,0.56,0.756,0.743,0.698,0.708,0.748,0.609,0.615,0.624,0.579,0.579,0.667,0.68,0.674,0.67,0.626,0.805,0.799,0.66,0.37,0.12,-0.0129999999999999,-0.32,-0.604,-0.894,-1.029,-1.23,-1.19,-1.051,-1.051,-0.911,-0.768,-0.768,-0.811,-0.527,-0.439,-0.579,-0.92,-0.971,-0.684,-0.568,-0.862,-0.606,-0.273,-0.557,-0.557,-0.25,-0.11,-0.36,-0.647,-0.562,-0.293,-0.583,-0.532,-0.498,-0.782,-0.797,-0.63,-0.395,-0.126,-0.527,-0.796,-0.796,-1.00222222222222,-0.87625,-0.665714285714286,-0.928333333333333,:::::
		# LINE,3,:::::DATA,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,:::::
		# TITLE,,,,,>chain1:::::
	
		# ///Draw Charts///
		$chains_dot_size = $#{$chains} + 1;
		$chart_coords .= 'NUM_CHAINS,' . $chains_dot_size . ',:::::';
		for ( my $c = 0; $c < $chains_dot_size; $c++ ) {

			my $res = Structure_getResChain_Range( $chains->[$c] );
			my $res_dot_size = $#{$res} + 1;

			$chain_length[$c] = $#{$res} + 1;
			if ($c == 0) {
				$chain_start_position[$c] = 0;
			} else {
				$chain_start_position[$c] = $chain_start_position[$c - 1] + $chain_length[$c - 1];
			}

			$chart_coords .= 'CHAIN,' . $c . ',' . $chain_length[$c] . ',' . $chain_start_position[$c] . ',:::::';

			for ( my $ln = 0; $ln <= $program_options->{'option_max_index'}; $ln++ ) { # for (unsigned int ln = 0; ln < 3; ln++){
				if ($program_options->{'lnOpt'}->[$ln]->{'showGraph'} == 1) {

					$chart_coords .= 'LINE,' . $ln . ',:::::DATA,';

					if ($ln == $program_options->{'lnOpt_index'}->{'PROFILES3D-LINE'}) {
						for ( my $i = 0; $i < $res_dot_size; $i++ ) {
							$chart_coords .= $res->[$i]->{'scores'}->[0] . ',';
						}
					} elsif ($ln == $program_options->{'lnOpt_index'}->{'REPIMPS-LINE'}) {
						for ( my $i = 0; $i < $res_dot_size; $i++ ) {
							$chart_coords .= $res->[$i]->{'scores'}->[1] . ',';
						}
					} elsif ($ln == $program_options->{'lnOpt_index'}->{'BOTH-LINE'}) {
						for ( my $i = 0; $i < $res_dot_size; $i++ ) {
							my $both_line = $res->[$i]->{'scores'}->[1] - $res->[$i]->{'scores'}->[0];
							$chart_coords .= $both_line . ',';
						}
					} elsif ($ln == $program_options->{'lnOpt_index'}->{'PROFILES3D-MOV-AVG-1'}) {
						for ( my $i = 0; $i < $res_dot_size; $i++ ) {
							$chart_coords .= $res->[$i]->{'mov_avg_1_score'}->[0] . ',';
						}
					} elsif ($ln == $program_options->{'lnOpt_index'}->{'REPIMPS-MOV-AVG-1'}) {
						for ( my $i = 0; $i < $res_dot_size; $i++ ) {
							$chart_coords .= $res->[$i]->{'mov_avg_1_score'}->[1] . ',';
						}
					} elsif ($ln == $program_options->{'lnOpt_index'}->{'TMS-MOV-AVG-1'}) {
						for (my $i = 0; $i < $res_dot_size; $i++) {
							if ($res->[$i]->{'isTransMem_area_mov_avg_1'} == 1) {
								$chart_coords .= '1,';
							} else {
								$chart_coords .= '0,';
							}
						}
					} elsif ($ln == $program_options->{'lnOpt_index'}->{'PROFILES3D-MOV-AVG-2'}) {
						for ( my $i = 0; $i < $res_dot_size; $i++ ) {
							$chart_coords .= $res->[$i]->{'mov_avg_2_score'}->[0] . ',';
						}
					} elsif ($ln == $program_options->{'lnOpt_index'}->{'REPIMPS-MOV-AVG-2'}) {
						for ( my $i = 0; $i < $res_dot_size; $i++ ) {
							$chart_coords .= $res->[$i]->{'mov_avg_2_score'}->[1] . ',';
						}
					} elsif ($ln == $program_options->{'lnOpt_index'}->{'TMS-MOV-AVG-2'}) {
						for (my $i = 0; $i < $res_dot_size; $i++) {
							if ($res->[$i]->{'isTransMem_area_mov_avg_2'} == 1) {
								$chart_coords .= '1,';
							} else {
								$chart_coords .= '0,';
							}
						}
					} elsif ($ln == $program_options->{'lnOpt_index'}->{'PROFILES3D-MOV-AVG-3'}) {
						for ( my $i = 0; $i < $res_dot_size; $i++ ) {
							$chart_coords .= $res->[$i]->{'mov_avg_3_score'}->[0] . ',';
						}
					} elsif ($ln == $program_options->{'lnOpt_index'}->{'REPIMPS-MOV-AVG-3'}) {
						for ( my $i = 0; $i < $res_dot_size; $i++ ) {
							$chart_coords .= $res->[$i]->{'mov_avg_3_score'}->[1] . ',';
						}
					} elsif ($ln == $program_options->{'lnOpt_index'}->{'HYDROPHILIC-MOV-AVG-3'}) {
						for (my $i = 0; $i < $res_dot_size; $i++) {
							if ($res->[$i]->{'is_hydrophilic_area_mov_avg_3'} == 1) {
								$chart_coords .= '1,';
							} else {
								$chart_coords .= '0,';
							}
						}
					} elsif ($ln == $program_options->{'lnOpt_index'}->{'PROFILES3D-MOV-AVG-4'}) {
						for ( my $i = 0; $i < $res_dot_size; $i++ ) {
							$chart_coords .= $res->[$i]->{'mov_avg_4_score'}->[0] . ',';
						}
					} elsif ($ln == $program_options->{'lnOpt_index'}->{'REPIMPS-MOV-AVG-4'}) {
						for ( my $i = 0; $i < $res_dot_size; $i++ ) {
							$chart_coords .= $res->[$i]->{'mov_avg_4_score'}->[1] . ',';
						}
					} elsif ($ln == $program_options->{'lnOpt_index'}->{'HYDROPHILIC-MOV-AVG-4'}) {
						for (my $i = 0; $i < $res_dot_size; $i++) {
							if ($res->[$i]->{'is_hydrophilic_area_mov_avg_4'} == 1) {
								$chart_coords .= '1,';
							} else {
								$chart_coords .= '0,';
							}
						}
					} elsif ($ln == $program_options->{'lnOpt_index'}->{'TMS'}) {
						for (my $i = 0; $i < $res_dot_size; $i++) {
							if ($res->[$i]->{'isTransMem'} == 1) {
								$chart_coords .= '1,';
							} else {
								$chart_coords .= '0,';
							}
						}
					} elsif ($ln == $program_options->{'lnOpt_index'}->{'TRUE-TMS'}) {
						if ($input_tms_line ne '') {
							my @true_array;
							for (my $i = 0; $i < $res_dot_size; $i++) {
								$true_array[$i] = 0;
							}
							my @bits = split(/,/, $input_tms_line);
							foreach my $true_segment (@bits) {
								my @bits2 = split(/-/, $true_segment);
								if (exists($bits2[1])) {
									my $start = $bits2[0] - 1;
									my $end = $bits2[1] - 1;
									for (my $t = $start; $t <= $end; $t++) {
										$true_array[$t] = 1;
									}
								}
							}
							for (my $i = 0; $i < $res_dot_size; $i++) {
								$chart_coords .= $true_array[$i] . ',';
							}
						}
					} elsif ($ln == $program_options->{'lnOpt_index'}->{'TRUE-INSIDE'}) {
						if ($input_inside_line ne '') {
							my @true_array;
							for (my $i = 0; $i < $res_dot_size; $i++) {
								$true_array[$i] = 0;
							}
							my @bits = split(/,/, $input_inside_line);
							foreach my $true_segment (@bits) {
								my @bits2 = split(/-/, $true_segment);
								if (exists($bits2[1])) {
									my $start = $bits2[0] - 1;
									my $end = $bits2[1] - 1;
									for (my $t = $start; $t <= $end; $t++) {
										$true_array[$t] = 1;
									}
								}
							}
							for (my $i = 0; $i < $res_dot_size; $i++) {
								$chart_coords .= $true_array[$i] . ',';
							}
						}
					} elsif ($ln == $program_options->{'lnOpt_index'}->{'TRUE-OUTSIDE'}) {
						if ($input_outside_line ne '') {
							my @true_array;
							for (my $i = 0; $i < $res_dot_size; $i++) {
								$true_array[$i] = 0;
							}
							my @bits = split(/,/, $input_outside_line);
							foreach my $true_segment (@bits) {
								my @bits2 = split(/-/, $true_segment);
								if (exists($bits2[1])) {
									my $start = $bits2[0] - 1;
									my $end = $bits2[1] - 1;
									for (my $t = $start; $t <= $end; $t++) {
										$true_array[$t] = 1;
									}
								}
							}
							for (my $i = 0; $i < $res_dot_size; $i++) {
								$chart_coords .= $true_array[$i] . ',';
							}
						}
					}

					$chart_coords .= ':::::';
				}
			}
		}

		my $graph_title = $input_fasta_id;
		$graph_title =~ s/\r//g;
		$graph_title =~ s/\n//g;
		$chart_coords .= 'TITLE,,,,,' . $graph_title . ':::::'; # has ,,,,, instead of , as the delimiter in case there is a comma in the actual title

		$multiplegraphs_data = $chart_coords;
		$multiplegraphs_data =~ s/\:\:\:\:\:/\;\;\;\;\;/g;
	}

	return $chart_coords;
}



sub CProtSumView_OnDraw { #C++ void CProtSumView::OnDraw(CDC* pDC)

	$structure = $analyser->{'st'};								# Structure *s = &(pDoc->a.st);

	my $output_line = '';

	$output_line .= "Protein summary :\n";							# lineOut("Protein summary for " + CString(pDoc->a.fileName.c_str()) + ":",pDC);
	$output_line .= "\n";
	if ($input_fasta_id ne '') {
		$output_line .= "Sequence identifier : $input_fasta_id\n";
		$output_line .= "\n";
	}
	my $s_arrow_chains_dot_size = $#{$structure->{'chains'}} + 1;
	$output_line .= "Number of chains: " . $s_arrow_chains_dot_size . "\n";			# lineOut("Number of chains: " + iToStr(s->chains.size()),pDC);
	my $s_arrow_allResidues_dot_size = $#{$structure->{'allResidues'}} + 1;
	$output_line .= "Number of Residues: " . $s_arrow_allResidues_dot_size . "\n";		# lineOut("Number of Residues: " + iToStr(s->allResidues.size()),pDC);
	$output_line .= "Total Solvent Accessible Surface Area: " . $structure->{'ASA'} . "\n";	# lineOut("Total Solvent Accessible Surface Area: " + fToStr(s->ASA),pDC);
	$output_line .= "Average Residue Score: " . $structure->{'aveResScore'} . "\n";		# lineOut("Average Residue Score: " + fToStr(s->aveResScore),pDC);
	$output_line .= "\n";
	$output_line .= "Predicted Transmembrane Domains:\n";					# lineOut("Predicted Transmembrane Domains:",pDC);
	$output_line .= "           Chain                   Residues\n";					# lineOut("           Chain           Residues",pDC);
	Structure_getTMDomains();							# vector<Range> tmDomains = s->getTMDomains();
	my $tmDomains_dot_size = $#{$tmDomains} + 1;
	for ( my $i = 0; $i < $tmDomains_dot_size; $i++ ) {				# for (unsigned int i = 0; i < tmDomains.size(); i++){
		# lineOut("             " + CString(tmDomains[i].ChainID) +"                     "+ 
		# 	iToStr(tmDomains[i].min) + "  -  " + iToStr(tmDomains[i].max),pDC);
		my $tms_start = $tmDomains->[$i]->{'min'} + 1;
		my $tms_end = $tmDomains->[$i]->{'max'} + 1;
		$output_line .= "             " . $tmDomains->[$i]->{'chainID'} .	"                     " . 
			$tms_start . "  -  " . $tms_end . "\n";
	}
	$output_line .= "\n\n";

	if (defined($input_flag_amphipathic)) {

		$output_line .= "Predicted Amphipathic Helices:\n";
		$output_line .= "           Chain                   Residues\n";
		get_amphipathic_helices();

		my $amphipathic_helices_dot_size = $#{$amphipathic_helices} + 1;
		for ( my $i = 0; $i < $amphipathic_helices_dot_size; $i++ ) {
			my $amphipathic_helix_start = $amphipathic_helices->[$i]->{'start_residue'};
			my $amphipathic_helix_end = $amphipathic_helices->[$i]->{'end_residue'};
			my $chain_id = ' ';
			$output_line .= "             " . $chain_id .	"                     " . 
				$amphipathic_helix_start . "  -  " . $amphipathic_helix_end . "\n";
		}
		$output_line .= "\n\n";
	}

	if ($program_mode eq 'cgi-output') {

		print $output_line;

	} else { # $program_mode eq 'command-line'

		$output_line .= "Input parameters that were used :\n";
		$output_line .= "\tSASA RESOLUTION (Index for Solvent Accessible Surface Area (SASA) resolution) = " . $program_options->{'SASARes'} . "\n";
		$output_line .= "\tTMS PREDICTION ALGORITHM (Algorithm for predicting transmembrane segments (TMS) from PROFILES3D and REPIMPS scores (VALPRED1 or VALPRED2)) = VALPRED" . $TMPred2D->{'valpred1_or_valpred2'} . "\n";
		if ($TMPred2D->{'valpred1_or_valpred2'} == 1) {
			$output_line .= "\tMINIMUM AREA DIFFERENCE = " . $TMPred2D->{'minAreaDiff_for_valpred1'} . "\n";
			$output_line .= "\tMOVING AVERAGE = " . $TMPred2D->{'mveAve'} . "\n";
			$output_line .= "\tMIN AV SASA = " . $TMPred2D->{'minAveSasa'} . "\n";
			$output_line .= "\tMIN TMS LENGTH = " . $TMPred2D->{'TMLengthMin_for_valpred1'} . "\n";
			$output_line .= "\tMAX TMS LENGTH = " . $TMPred2D->{'TMLengthMax_for_valpred1'} . "\n";
			$output_line .= "\tMIN AV PROFILES3D RANGE = " . $TMPred2D->{'aveRepRangeMin'} . "\n";
			$output_line .= "\tMAX AV PROFILES3D RANGE = " . $TMPred2D->{'aveRepRangeMax'} . "\n";
			$output_line .= "\tMIN AV REPIMPS RANGE = " . $TMPred2D->{'aveProfRangeMin'} ."\n";
			$output_line .= "\tMAX AV REPIMPS RANGE = " . $TMPred2D->{'aveProfRangeMax'} . "\n";
			$output_line .= "\tMIN LENGTH DIVIDED BY AREA = " . $TMPred2D->{'lenDivAreaMin'} . "\n";
			$output_line .= "\tMAX LENGTH DIVIDED BY AREA = " . $TMPred2D->{'lenDivAreaMax'} . "\n";
			$output_line .= "\tMIN SCORE DIFFERENCE = " . $TMPred2D->{'minScoreDiff_for_valpred1'} . "\n";
			$output_line .= "\tASA were available (surface areas were calculated from residues to give PROFILES3D and REPIMPS scores, the scores were not read in on input) = ";
			if ($TMPred2D->{'ASA_are_available'} == 1) {
				$output_line .= "Yes" . "\n";
			} else { # $TMPred2D->{'ASA_are_available'} == 0
				$output_line .= "No" . "\n";
			}
		} else { # $TMPred2D->{'valpred1_or_valpred2'} == 2
			$output_line .= "\tMINIMUM AREA DIFFERENCE = " . $TMPred2D->{'minAreaDiff_for_valpred2'} . "\n";
			$output_line .= "\tMV.AV.1 = " . $TMPred2D->{'mov_avg_1'} . "\n";
			$output_line .= "\tMV.AV.2 = " . $TMPred2D->{'mov_avg_2'} . "\n";
			$output_line .= "\tMV.AV.3 = " . $TMPred2D->{'mov_avg_3'} . "\n";
			$output_line .= "\tMV.AV.4 = " . $TMPred2D->{'mov_avg_4'} . "\n";
			$output_line .= "\tMIN TMS LENGTH = " . $TMPred2D->{'TMLengthMin_for_valpred2'} . "\n";
			$output_line .= "\tMAX TMS LENGTH = " . $TMPred2D->{'TMLengthMax_for_valpred2'} . "\n";
			$output_line .= "\tMAX LENGTH NON-TMS = " . $TMPred2D->{'max_nonTMS_length'} . "\n";
			$output_line .= "\tMIN LENGTH NON-TMS = " . $TMPred2D->{'min_nonTMS_length'} . "\n";
			$output_line .= "\tMIN NON-TMS SCORE DIFFERENCE = " . $TMPred2D->{'min_nonTMS_score_diff'} . "\n";
			$output_line .= "\tMIN NON-TMS AREA = " . $TMPred2D->{'min_nonTMS_area'} . "\n";
			$output_line .= "\tTWILIGHT TMS ZONE AREA PER RESIDUE - LOWER LIMIT = " . $TMPred2D->{'twilight_area_per_residue_lower_limit'} . "\n";
			$output_line .= "\tTWILIGHT TMS ZONE AREA PER RESIDUE - UPPER LIMIT = " . $TMPred2D->{'twilight_area_per_residue_upper_limit'} . "\n";
			$output_line .= "\tIGNORE TWILIGHT TMS ZONE = ";
			if ($TMPred2D->{'ignore_twilight_area'} == 1) {
				$output_line .= "Yes" . "\n";
			} else { # $TMPred2D->{'ignore_twilight_area'} == 0
				$output_line .= "No" . "\n";
			}
			$output_line .= "\tNO. RESIDUES SEARCH AREA FOR TMS NEIGHBOURS = " . $TMPred2D->{'search_area_for_neighbour_tms'} . "\n";
			$output_line .= "\tMIN HELIX LENGTH = " . $TMPred2D->{'helix_length_min'} . "\n";
			$output_line .= "\tMIN SCORE DIFF FOR TMS ENDS (MOV.AVG.1) = " . $TMPred2D->{'min_score_diff_for_TMS_ends'} . "\n";
			$output_line .= "\tMIN SCORE DIFF. MOV.AVG.1 = " . $TMPred2D->{'minScoreDiff_movavg1_for_valpred2'} . "\n";
			$output_line .= "\tMIN SCORE DIFF. MOV.AVG.2 = " . $TMPred2D->{'minScoreDiff_movavg2_for_valpred2'} . "\n";
			$output_line .= "\tMIN AV AREA MOV.AVG.1 = " . $TMPred2D->{'min_avg_area_movavg1'} . "\n";
			$output_line .= "\tMIN AV AREA MOV.AVG.2 = " . $TMPred2D->{'min_avg_area_movavg2'} . "\n";
		}
		$output_line .= "\n";

		# open output file
		my $summary_output_file = "$input_file_name_for_output_files.valpred_summary.txt";
		open( SUMOUTFILE, ">$summary_output_file") or
			die "Cannot open $summary_output_file for writing : $!\n";

		print SUMOUTFILE $output_line;

		close SUMOUTFILE;
	}
}



sub create_output_graphs_html_file {

	my $what_to_do_flag = shift;

	my $output_line = '';

	if ($what_to_do_flag == 1) {

		my $title = 'VALPRED GRAPH RESULTS';
		$output_line .= "<html><head><title>$title</title></head><body>\n";
		$output_line .= html_code_for_graphics_header($title);
		$output_line .= "<b>VALPRED GRAPH RESULTS for $input_file_name_for_output_files</b><br><br>\n";
		print HTMLGRAPHSOUTFILE $output_line;

	} elsif ($what_to_do_flag == 2) {

		my $display_fasta_id = $input_fasta_id;
		$display_fasta_id =~ s/>/\>/g;
		$display_fasta_id =~ s/</\</g;
		$output_line .= "$input_sequence_id&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$display_fasta_id<br>\n";
		$output_line .= "<a href='$multiplegraphs_filename' target='_blank'>Click here to see graph in a new window</a><br>\n";
		$output_line .= "<img src='$multiplegraphs_filename' border='0' height='250'><br><br>\n";
		print HTMLGRAPHSOUTFILE $output_line;

		close HTMLGRAPHSOUTFILE;
		open( HTMLGRAPHSOUTFILE, ">>$html_graphs_output_file") or
			die "Cannot open $html_graphs_output_file for writing : $!\n";

	} elsif ($what_to_do_flag == 3) {

		$output_line .= "</td></tr></table>\n";
		$output_line .= "</body></html>\n";
		print HTMLGRAPHSOUTFILE $output_line;
	}
}



sub create_structure_from_input_scores { 

	my $profiles3d_scores_string = shift;
	my $repimps_scores_string = shift;
	my $aaseq_string = shift;

	my %code_1to3;
	$code_1to3{'A'} = 'ALA';
	$code_1to3{'R'} = 'ARG';
	$code_1to3{'N'} = 'ASN';
	$code_1to3{'D'} = 'ASP';
	$code_1to3{'C'} = 'CYS';
	$code_1to3{'Q'} = 'GLN';
	$code_1to3{'E'} = 'GLU';
	$code_1to3{'G'} = 'GLY';
	$code_1to3{'H'} = 'HIS';
	$code_1to3{'I'} = 'ILE';
	$code_1to3{'L'} = 'LEU';
	$code_1to3{'K'} = 'LYS';
	$code_1to3{'M'} = 'MET';
	$code_1to3{'F'} = 'PHE';
	$code_1to3{'P'} = 'PRO';
	$code_1to3{'S'} = 'SER';
	$code_1to3{'T'} = 'THR';
	$code_1to3{'W'} = 'TRP';
	$code_1to3{'Y'} = 'TYR';
	$code_1to3{'V'} = 'VAL';
	$code_1to3{'X'} = 'UNK';

	my @profiles3d_array = split( '\,', $profiles3d_scores_string );
	my @repimps_array = split( '\,', $repimps_scores_string );

	my $residues;

	for ( my $i = 0; $i < @profiles3d_array; $i++ ) {
		$residues->[$i]->{'scores'}->[0] = $profiles3d_array[$i];
		$residues->[$i]->{'scores'}->[1] = $repimps_array[$i];
		my $code1 = substr( $aaseq_string, $i, 1 );
		my $code3 = $code_1to3{$code1};
		$residues->[$i]->{'name'} = $code3;
		$residues->[$i]->{'resSeq'} = $i;
		$residues->[$i]->{'chainID'} = 0;
		$residues->[$i]->{'side_ASA'} = -1;
		$residues->[$i]->{'amphipathic_calc_side_chain_not_buried_polar_ASA'} = -1;
	}

	$structure->{'allResidues'} = $residues;
	$structure->{'chains'}->[0]->{'min'} = 0;
	$structure->{'chains'}->[0]->{'max'} = @profiles3d_array - 1;
	$structure->{'ASA'} = -1;
	$structure->{'aveResScore'} = -1;
	$analyser->{'st'} = $structure;
}



sub CTableView_fillTable { #C++ void CTableView::fillTable(void)

	# if $program_mode eq 'command-line'
	my $field_start = "";
	my $field_end = ",";
	my $line_start = "";
	my $line_end = "\n";
	my $table_start = "";
	my $table_end = "";
	if ($program_mode eq 'cgi-output') {
		$field_start = "<td>";
		$field_end = "</td>";
		$line_start = "<tr>";
		$line_end = "</tr>\n";
		$table_start = "<table cellspacing=0 cellpadding=2 border=0>";
		$table_end = "</table>\n";
	}
	my $calc_output = $table_start;

	# void CTableView::OnInitialUpdate()
	my $calc_hdr = $line_start;
	$calc_hdr .= $field_start . "Index" . $field_end;				# ListCtrl.InsertColumn(0,"Index",LVCFMT_LEFT,30,0 );
	$calc_hdr .= $field_start . "Residue" . $field_end;				# ListCtrl.InsertColumn(1,"Residue" ,LVCFMT_LEFT,40,1 );
	$calc_hdr .= $field_start . "Chain" . $field_end;				# ListCtrl.InsertColumn(2,"Chain" ,LVCFMT_LEFT,40,1 );
	$calc_hdr .= $field_start . "No." . $field_end;					# ListCtrl.InsertColumn(3,"No." ,LVCFMT_LEFT,60,2 );
	$calc_hdr .= $field_start . "Area Buried" . $field_end;				# ListCtrl.InsertColumn(4,"Area Buried" ,LVCFMT_LEFT,80,3 );
	$calc_hdr .= $field_start . "Fraction Polar" . $field_end;			# ListCtrl.InsertColumn(5,"Fraction Polar" ,LVCFMT_LEFT,80,4 );
	$calc_hdr .= $field_start . "REPIMPS FP" . $field_end;				# ListCtrl.InsertColumn(6,"REPIMPS FP" ,LVCFMT_LEFT,70,5 );
	$calc_hdr .= $field_start . "Category" . $field_end;				# ListCtrl.InsertColumn(7,"Category" ,LVCFMT_LEFT,60 ,6);
	$calc_hdr .= $field_start . "REPIMPS Cat" . $field_end;				# ListCtrl.InsertColumn(8,"REPIMPS Cat" ,LVCFMT_LEFT,70 ,7);
	$calc_hdr .= $field_start . "Main SASA" . $field_end;				# ListCtrl.InsertColumn(9,"Main SASA" ,LVCFMT_LEFT,70 ,8);
	$calc_hdr .= $field_start . "Side SASA" . $field_end;				# ListCtrl.InsertColumn(10,"Side SASA" ,LVCFMT_LEFT,70 ,9);
	$calc_hdr .= $field_start . "Side not-buried polar ASA" . $field_end;
	$calc_hdr .= $field_start . "SASA" . $field_end;				# ListCtrl.InsertColumn(11,"SASA" ,LVCFMT_LEFT,60 ,10);
	$calc_hdr .= $field_start . "Solvent Score" . $field_end;			# ListCtrl.InsertColumn(12,"Solvent Score" ,LVCFMT_LEFT,60 ,11);
	$calc_hdr .= $field_start . "REPIMPS Score" . $field_end;			# ListCtrl.InsertColumn(13,"REPIMPS Score" ,LVCFMT_LEFT,60 ,12);
	$calc_hdr .= $field_start . "Overall Score" . $field_end;			# ListCtrl.InsertColumn(14,"Overall Score" ,LVCFMT_LEFT,60 ,13);
	$calc_hdr .= $field_start . "is transmembrane" . $field_end;			# ListCtrl.InsertColumn(15,"is transmembrane," ,LVCFMT_LEFT,60 ,13);
	$calc_hdr .= $line_end;

	$calc_output .= $calc_hdr;

	my $r = $analyser->{'st'}->{'allResidues'};					# vector<Residue> &r = GetDocument()->a.st.allResidues;

	my $i;										# unsigned int i;
	my $r_dot_size = $#{$r} + 1;
	for ( $i = 0; $i < $r_dot_size; $i++ ) {					# for (i = 0; i < r.size();  i++){
		my $calc_line = $line_start;
		$calc_line .= $field_start . $i . $field_end;				# buff.Format("%d",i);
		$calc_line .= $field_start . $r->[$i]->{'name'} . $field_end;		# ListCtrl.SetItemText(i, 1, r[i].name.c_str());
		$calc_line .= $field_start . $r->[$i]->{'chainID'} . $field_end;	# ListCtrl.SetItemText(i, 2, CString(r[i].chainID));
		$calc_line .= $field_start . $r->[$i]->{'resSeq'} . $field_end;		# buffer.Format("%d",r[i].resSeq);
		$calc_line .= $field_start . $r->[$i]->{'areaBuried'} . $field_end;	# buffer.Format("%f",r[i].areaBuried);
		$calc_line .= $field_start . $r->[$i]->{'fractionPolar'} . $field_end;	# buffer.Format("%f",r[i].fractionPolar);
		$calc_line .= $field_start . $r->[$i]->{'repimpsFP'} . $field_end; 	# buffer.Format("%f",r[i].repimpsFP);
		$calc_line .= $field_start . "profEnv" . $field_end; 			# buffer = "profEnv";
		$calc_line .= $field_start . "repimpsEnv" . $field_end; 		# buffer = "repimpsEnv";
		$calc_line .= $field_start . $r->[$i]->{'main_ASA'} . $field_end;	# buffer.Format("%f",r[i].main_ASA);
		$calc_line .= $field_start . $r->[$i]->{'side_ASA'} . $field_end;	# buffer.Format("%f",r[i].side_ASA);
		#$calc_line .= $field_start . $r->[$i]->{'amphipathic_calc_side_chain_not_buried_polar_ASA'} . $field_end;
		my $sasa = $r->[$i]->{'main_ASA'} + $r->[$i]->{'side_ASA'}; 		# buffer.Format("%f",r[i].main_ASA + r[i].side_ASA);
		$calc_line .= $field_start . $sasa . $field_end;
		$calc_line .= $field_start . $r->[$i]->{'scores'}->[0] . $field_end;	# buffer.Format("%f",r[i].scores[0]);
		$calc_line .= $field_start . $r->[$i]->{'scores'}->[1] . $field_end;	# buffer.Format("%f",r[i].scores[1]);
		$calc_line .= $field_start . $r->[$i]->{'scores'}->[2] . $field_end;	# buffer.Format("%f",r[i].scores[2]);
		$calc_line .= $field_start . $r->[$i]->{'isTransMem'} . $field_end;	# buffer.Format("%d",(int)(r[i].isTransMem));
		$calc_line .= $line_end;
		$calc_output .= $calc_line;
	}
	$calc_output .= $table_end;
	$calc_output .= "\n\n";

	if ($program_mode eq 'cgi-output') {

		print $calc_output;

	} else { # $program_mode eq 'command-line'

		# open output file
		my $calculations_output_file = "$input_file_name_for_output_files.valpred_calculations.txt";
		open( CALCOUTFILE, ">$calculations_output_file") or
			die "Cannot open $calculations_output_file for writing : $!\n";

		print CALCOUTFILE $calc_output;

		close CALCOUTFILE;
	}
}



sub CValpredDoc_OnToolsTmpred2d {

	my $ref_aa_seq_chains = shift; #C++ string
	my $seq_num = shift;

	if (!defined ($input_graph_points_file)) {
		# input file is $input_infile, $input_fasta_chain_file or $input_fasta_file
		# if input file is $input_graph_points_file, then valpred points 
		#	(scores and moving averages for plotting on a graph) have already been received on input
		# if input file is $input_inscoresfile, then the valpred scores have already been received on input (and have been built into the structure as if it came from a helix)
		#	(don't need to build a helix from residues, but do need to calculate moving averages, find TMS, and create graphs from scratch)

		if (!defined ($input_inscoresfile)) {

			my $genOpts = GeneralOpts_new();						# GeneralOpts genOpts;

			# for this FASTA sequence, add chains to sequence

			my $structSeq = $ref_aa_seq_chains; #C++ vector<string> 			# structSeq.push_back(sequence);

			Analyser_new();									# Analyser an;

			$helices = Range_new(); #C++ Range
			CValpredDoc_seqToHelix( $structSeq ); #C++ Structure				# Structure s = seqToHelix(structSeq);

													# if (s.allAtoms.size() == 0){
													# AfxMessageBox("Error: Could not process FASTA file");
													# return;

			$analyser->{'st'} = $structure;							# an.st  = s;
													# an.fileName = LPCTSTR(newOpenBox.GetFileName());
													# //an.logFilePath = currentDir + "LogFile.log";
			$analyser->{'ignoreHet'} = $genOpts->{'IgnoreHetAtoms'};			# an.ignoreHet = genOpts.IgnoreHetAtoms;
			Analyser_findCloseAtoms();							# an.findCloseAtoms();

			Analyser_SASA( $genOpts->{'SASARes'} );						# an.SASA(genOpts.SASARes);

			Analyser_doProfile();								# an.doProfile();
		}

		# my $TMPred2D = TMPred2DParam_TMPred2DParam();
		# an.predTM2D(TMPred2D.mveAve,TMPred2D.TMLengthMin,TMPred2D.TMLengthMax,TMPred2D.minScoreDiff,TMPred2D.minAreaDiff,
		#		TMPred2D.aveProfRangeMin,TMPred2D.aveProfRangeMax,TMPred2D.aveRepRangeMin,TMPred2D.aveRepRangeMax,
		#		TMPred2D.lenDivAreaMin,TMPred2D.lenDivAreaMax,TMPred2D.minAveSasa);

		if ($TMPred2D->{'valpred1_or_valpred2'} == 2) {

			Analyser_predTM2D_valpred2();

		} else { # $TMPred2D->{'valpred1_or_valpred2'} == 1

			Analyser_predTM2D_valpred1();
		}

		if ($TMPred2D->{'predict_amphipathic_helices'} == 1) {

			predict_amphipathic_helices();
		}

		Analyser_setOverallScores();							# an.setOverallScores();

		if (($program_mode eq 'command-line') && (defined($input_flag_benchmark))) {
			output_benchmark_file();
		}

		if (($program_mode eq 'command-line') && (defined($input_tms_file))) {
			output_evaluate_file();
		}
	}

	if ( ($seq_num == 1) ||
		(($program_mode eq 'command-line') && (defined($input_flag_multiplegraphs))) ) {

		# for command-line mode if $seq_num = 1,
		#  - output file xxx.valpred_summary.txt (eg. KCSA.fasta.1line.valpred_summary.txt)
		#       which contains the valpred transmembrane summary for the first sequence in the input file
		# for cgi-script mode (which only does one sequence whose $seq_num = 1), output it to the screen
		if (($seq_num == 1) && (!defined ($input_graph_points_file))) {
			if ($program_mode eq 'cgi-output') {
				print "<pre>\n";
			}

			CProtSumView_OnDraw();

			if ($program_mode eq 'cgi-output') {
				print "</pre>\n";
			}
		}

		# for command-line mode, if $seq_num = 1, 
		#  - output file xxx.valpred_graph.gif (eg. KCSA.fasta.1line.valpred_graph.gif)
		#       which contains the valpred transmembrane graph for the first sequence in the input file
		# if -multiplegraphs option, output file xxx.valpred_graph_aaa.gif
		#	and add this graph to output file xxx.valpred_graphs.html
		# for cgi-script mode, output a link to the graph instead of outputting the graph itself
		my $chart_coords = CChartView_OnDraw();

		if ($program_mode eq 'command-line') {
			if (($seq_num == 1) && (defined($input_flag_multiplegraphs))) { # output xxx.valpred_graph.gif and xxx.vapred_graph_aaa.gif
				draw_chart_from_coords( $chart_coords, 3 );
			} elsif ($seq_num == 1) { # output xxx.valpred_graph.gif only
				draw_chart_from_coords( $chart_coords, 1 );
			} elsif (defined($input_flag_multiplegraphs)) { # output xxx.valpred_graph_aaa.gif only
				draw_chart_from_coords( $chart_coords, 2 );
			}
			if (defined($input_flag_multiplegraphs)) {
				create_output_graphs_html_file( 2 );
			}
		} else { # $program_mode eq 'cgi-script'
			output_chart_coords_link( $chart_coords );
		}

		# for command-line mode, if $seq_num = 1, 
		#  - output file xxx.valpred_calculations.txt (eg. KCSA.fasta.1line.valpred_calculations.txt)
		#       which contains the valpred calculation results of accessibility and hydropathy environment for each residue
		# for cgi-script mode (which only does one sequence whose $seq_num = 1), output it to the screen
		if (($seq_num == 1) && (!defined ($input_graph_points_file)) && (!defined ($input_inscoresfile))) {
			CTableView_fillTable();
		}
	}

	if (($program_mode eq 'command-line') && (!defined ($input_graph_points_file))) {
		#  - output file xxx.valpred_segments.txt (eg. KCSA.fasta.1line.valpred_segments.txt)
		#       which is the valpred transmembrane segments output file with one line per input sequence
		write_output_segments_and_graphpoints_files();
	}

	if (($program_mode eq 'command-line') && (defined ($input_flag_amphipathic))) {
		#  - output file xxx.valpred_amphipathic.txt (eg. KCSA.fasta.1line.valpred_amphipathic.txt)
		#       which is the valpred amphipathic helix segments output file with one line per input sequence
		write_output_amphipathic_helices();
	}
}



sub CValpredDoc_seqToHelix { # Structure CValpredDoc::seqToHelix(vector<string> chains){

	my $chains = shift; #C++ vector<string>

	my $p = PolymerBuilder_new(); #C++ PolymerBuilder				# PolymerBuilder p;

	PolymerBuilder_buildStructure( $p, $chains, -62, -41, 180 );			# return p.buildStructure(chains,-62, -41, 180);

	return;
}



sub draw_chart_from_coords { 

	my $chart_coords = shift;
	my $control_flag_for_output_files = shift;

	my $noOfRes = 0;
	my $chains_dot_size = 0;
	my @chain_length;
	my @chain_start_position;
	my $graph_title = '';

	my @chart_coords_lines = split( /:::::/, $chart_coords );
	my $c = -1;
	my $ln = -1;
	my @chart_coords_scores;
	foreach my $chart_coords_line (@chart_coords_lines) {
		my @bits = split( /,/, $chart_coords_line );
		if ($bits[0] eq 'NUM_RESIDUES') {
			$noOfRes = trim($bits[1]);
		} elsif ($bits[0] eq 'NUM_CHAINS') {
			$chains_dot_size = trim($bits[1]);
		} elsif ($bits[0] eq 'CHAIN') {
			$c = trim($bits[1]);
			$chain_length[$c] = trim($bits[2]);
			$chain_start_position[$c] = trim($bits[3]);
		} elsif ($bits[0] eq 'LINE') {
			$ln = trim($bits[1]);
		} elsif ($bits[0] eq 'DATA') {
			for ( my $i2 = 1; $i2 < @bits; $i2++ ) {
				my $i = $i2 - 1;
				$chart_coords_scores[$c][$ln][$i] = trim($bits[$i2]);
			}
		} elsif ($bits[0] eq 'TITLE') {
			my @bits2 = split( /,,,,,/, $chart_coords_line ); # has ,,,,, instead of , as the delimiter in case there is a comma in the actual title
			$graph_title = trim($bits2[1]);
		} elsif ($bits[0] eq 'GRAPH_RESIDUE_WIDTH') {
			$program_options->{'graph_residue_width'} = trim($bits[1]);
		} elsif ($bits[0] eq 'GRAPH_HEIGHT') {
			$program_options->{'graph_height'} = trim($bits[1]);
		} elsif ($bits[0] eq 'VALPRED1_or_VALPRED2') {
			$TMPred2D->{'valpred1_or_valpred2'} = trim($bits[1]);
		}
	}

	my @graph_data;
	my @graph_data_line_id;
	my @graph_data_line_smooth;
	my @graph_data_x_axis;
	for (my $d = 0; $d < $noOfRes; $d++) {
		$graph_data_x_axis[$d] = $d + 1;
	}

	my $graph_data_line_index = 0;
	$graph_data[0] = \@graph_data_x_axis;
	$graph_data_line_id[0] = 'X-AXIS';
	$graph_data_line_smooth[0] = 0;
	my @graph_colors;
	my @graph_legend;

	# ///Draw Charts///	
	for ( $c = 0; $c < $chains_dot_size; $c++ ) {

		my @continuous_lines = ('PROFILES3D-LINE','REPIMPS-LINE','BOTH-LINE','PROFILES3D-MOV-AVG-1','REPIMPS-MOV-AVG-1','PROFILES3D-MOV-AVG-2','REPIMPS-MOV-AVG-2',
			'PROFILES3D-MOV-AVG-3','REPIMPS-MOV-AVG-3','PROFILES3D-MOV-AVG-4','REPIMPS-MOV-AVG-4');

		for ( my $o = 0; $o < @continuous_lines; $o++ ) {

			my $ln_name = $continuous_lines[$o];
			my $ln = $program_options->{'lnOpt_index'}->{$ln_name};

			if ( $program_options->{'lnOpt'}->[$ln]->{'showGraph'} == 1 ) {

				$graph_data_line_index++;
				my $graph_data_residue_index = $chain_start_position[$c];
				$graph_colors[$graph_data_line_index-1] = $program_options->{'lnOpt'}->[$ln]->{'colour_name'};
				if ($c == 0) {
					$graph_legend[$graph_data_line_index-1] = $program_options->{'lnOpt'}->[$ln]->{'legend_name'};
				}
				$graph_data_line_id[$graph_data_line_index] = $program_options->{'lnOpt'}->[$ln]->{'line_id'};
				$graph_data_line_smooth[$graph_data_line_index] = $program_options->{'lnOpt'}->[$ln]->{'smooth'};

				for ( my $i = 0; $i < $chain_length[$c]; $i++ ) {

					$graph_data[$graph_data_line_index][$graph_data_residue_index] = $chart_coords_scores[$c][$ln][$i];
					$graph_data_residue_index++;
				}
			}
		}


		my @discontinuous_lines = ('TMS-MOV-AVG-1','TMS-MOV-AVG-2','HYDROPHILIC-MOV-AVG-3','HYDROPHILIC-MOV-AVG-4',
			'TMS','TRUE-TMS','TRUE-INSIDE','TRUE-OUTSIDE');

		for ( my $o = 0; $o < @discontinuous_lines; $o++ ) {

			my $ln_name = $discontinuous_lines[$o];
			my $ln = $program_options->{'lnOpt_index'}->{$ln_name};

			my $have_started_a_new_graph_data_line = 0;
			my $have_done_the_legend = 0;
			my $res_dot_size = $chain_length[$c];
			for (my $i = 0; $i < $res_dot_size; $i++) {
				if ($program_options->{'lnOpt'}->[$ln]->{'showGraph'} == 1) {
					if ($chart_coords_scores[$c][$ln][$i] == 1) {
						# TMS have to be printed over more than one graph line
						# because the graph module does not handle discontinuous lines.
						# Each segment of a discontinuous line is printed as a new graph line.
						if ($have_started_a_new_graph_data_line == 0) {
							$graph_data_line_index++;
							$graph_colors[$graph_data_line_index-1] = $program_options->{'lnOpt'}->[$ln]->{'colour_name'};
							$have_started_a_new_graph_data_line = 1;
							if (($c == 0) && ($have_done_the_legend == 0)) {
								$graph_legend[$graph_data_line_index-1] = $program_options->{'lnOpt'}->[$ln]->{'legend_name'};
								$have_done_the_legend = 1;
							}
							$graph_data_line_id[$graph_data_line_index] = $program_options->{'lnOpt'}->[$ln]->{'line_id'};
							$graph_data_line_smooth[$graph_data_line_index] = $program_options->{'lnOpt'}->[$ln]->{'smooth'};
							for (my $j = 0; $j < $noOfRes; $j++) {
								my $graph_data_residue_index = $chain_start_position[$c] + $j;
								$graph_data[$graph_data_line_index][$graph_data_residue_index] = undef;
							}
						}
						my $graph_data_residue_index = $chain_start_position[$c] + $i;
						$graph_data[$graph_data_line_index][$graph_data_residue_index] = $program_options->{'lnOpt'}->[$ln]->{'y_value'};
					} else {
						$have_started_a_new_graph_data_line = 0;
					}
				}
			}
		}
	}

	# 1 residue is 10 pixels wide
	my $graph_width = $noOfRes * $program_options->{'graph_residue_width'};
	if ($graph_width < 600) {	# make sure the graph plot is wider than the legend,
		$graph_width = 600; 	# otherwise GD/Graph/axestype.pm will crash.
	}
	my $graph_height = $program_options->{'graph_height'};
	# $graph_width = $noOfRes * 20; 	# 10 	or 5
	# $graph_height = 800; 		# 600 	or 200
	my $graph_object = new GD::Graph::lines( $graph_width, $graph_height );

	my $graph_y_max_value = 2; # 2
	my $graph_y_min_value = -2; # -2
	my $num_duplications = 1;

	# make the graph lines look a more like a histogram than pointy,
	# by duplicating all values 4 times

	my @new_graph_data;
	$graph_data_line_index = -1;
	foreach my $graph_line_array_ref (@graph_data) {
		$graph_data_line_index++;
		my @graph_line_array = @$graph_line_array_ref;
		my @new_graph_line_array;
		my $line_id = $graph_data_line_id[$graph_data_line_index];
		if ($graph_data_line_smooth[$graph_data_line_index] == 1 ) {
			my @temp_graph_line_array = @graph_line_array;
			$temp_graph_line_array[@temp_graph_line_array] = $graph_line_array[(@graph_line_array-1)];
			for (my $i = 0; $i < @graph_line_array; $i++) {
				if ($i < ($noOfRes-1)) {
					if (defined($temp_graph_line_array[$i])) {
						my $diff = ($temp_graph_line_array[$i+1] - $temp_graph_line_array[$i]) / 4;
						$new_graph_line_array[@new_graph_line_array] = $temp_graph_line_array[$i];
						$new_graph_line_array[@new_graph_line_array] = $temp_graph_line_array[$i] + $diff;
						$new_graph_line_array[@new_graph_line_array] = $temp_graph_line_array[$i] + ($diff * 2);
						$new_graph_line_array[@new_graph_line_array] = $temp_graph_line_array[$i] + ($diff * 3);
					} else { # the value is undef, due to making room for chains before this chain
						$new_graph_line_array[@new_graph_line_array] = $temp_graph_line_array[$i];
						$new_graph_line_array[@new_graph_line_array] = $temp_graph_line_array[$i];
						$new_graph_line_array[@new_graph_line_array] = $temp_graph_line_array[$i];
						$new_graph_line_array[@new_graph_line_array] = $temp_graph_line_array[$i];
					}
				} else { # the value is undef, that was padded out to stop the graph from being less wide than the legend, which would cause a crash
					$new_graph_line_array[@new_graph_line_array] = $temp_graph_line_array[$i];
					$new_graph_line_array[@new_graph_line_array] = $temp_graph_line_array[$i];
					$new_graph_line_array[@new_graph_line_array] = $temp_graph_line_array[$i];
					$new_graph_line_array[@new_graph_line_array] = $temp_graph_line_array[$i];
				}
			}
		} else {
			foreach my $graph_line_point (@graph_line_array) {
				$new_graph_line_array[@new_graph_line_array] = $graph_line_point;
				$new_graph_line_array[@new_graph_line_array] = $graph_line_point;
				$new_graph_line_array[@new_graph_line_array] = $graph_line_point;
				$new_graph_line_array[@new_graph_line_array] = $graph_line_point;
			}
		}
		$new_graph_data[@new_graph_data] = \@new_graph_line_array;
	}
	@graph_data = @new_graph_data;
	$num_duplications = 4;

	my $requested_x_label_skip = 2;
	my $graph_x_label_skip = $requested_x_label_skip * $num_duplications;

	$graph_object->set( 
		# x_label => '',
		# y_label => '',
		title => $graph_title,
		y_max_value => $graph_y_max_value,
		y_min_value => $graph_y_min_value,
		y_tick_number => 4,
		box_axis => 0,
		line_width => 3,
		zero_axis_only => 1,
		x_label_position => 0,
		y_label_position => 0,
		x_label_skip => $graph_x_label_skip,
		# x_tick_offset => 0,
		# y_label_skip => 0,
		transparent => 0,
		legend_placement => 'B',
	);

	$graph_object->set_legend(@graph_legend);

	$graph_object->set(dclrs=>\@graph_colors);

	$graph_object->set_legend_font(['verdana', 'arial', gdGiantFont], 12);

	$graph_object->plot(\@graph_data);

	my $ext = $graph_object->export_format;

	if ($program_mode eq 'cgi-output-graph') {

		# Apparently we need to do this on Windows or the image will be garbled, and it doesn't hurt on Unix/Linux/etc.  
		binmode STDOUT;

		print $q->header("Content-type: image/gif");
		print $graph_object->gd->$ext();

	} else { # $program_mode eq 'command-line'

		if (($control_flag_for_output_files == 1) || ($control_flag_for_output_files == 3)) {

			my $graph_output_file_1 = "$input_file_name_for_output_files.valpred_graph.$ext";

			open( GRFOUTFILE1, ">$graph_output_file_1") or
				die "Cannot open $graph_output_file_1 for writing : $!\n";
			binmode GRFOUTFILE1;
			print GRFOUTFILE1 $graph_object->gd->$ext();
			close GRFOUTFILE1;
		}

		if (($control_flag_for_output_files == 2) || ($control_flag_for_output_files == 3)) {

			my $filename_bit = $input_sequence_id;
			$filename_bit =~ s/\W//g;
			if (substr($filename_bit,0,4) ne 'SEQ_') {
				$filename_bit = "SEQ_$filename_bit";
			}
			my $graph_output_file_2 = "$input_file_name_for_output_files.valpred_graph_$filename_bit.$ext";
			$multiplegraphs_filename = $graph_output_file_2;

			open( GRFOUTFILE2, ">$graph_output_file_2") or
				die "Cannot open $graph_output_file_2 for writing : $!\n";
			binmode GRFOUTFILE2;
			print GRFOUTFILE2 $graph_object->gd->$ext();
			close GRFOUTFILE2;
		}
	}
}



sub GeneralOpts_new {

	my $genOpts;

	$genOpts->{'IgnoreHetAtoms'} = 0;		# , IgnoreHetAtoms(FALSE)
	$genOpts->{'SASARes'} = $program_options->{'SASARes'}; # , SASARes(2)

	return $genOpts; #C++ GeneralOpts
}



# // find the angle, in radian between 2 vectors.
sub geometry_angle { #C++ inline Type angle(Vec3<Type> v1, Vec3<Type> v2)

	my $v1 = shift;
	my $v2 = shift;

	my $bit1 = Vec3_dot_Vec3( $v1, $v2 );
	my $bit2 = Vec3_length( $v1 );

	my $bit3 = Vec3_length( $v2 );

	my $cosAngle = $bit1 / $bit2 / $bit3;		# Type cosAngle = v1.dot(v2)/v1.length()/v2.length();

	if ($cosAngle >= 1) {				# if (cosAngle >= 1) { // sometimes fail due to rounding errors
		return 0;				# return 0;
	}

	my $acos = acos( $cosAngle );

	return $acos;					# return acos(cosAngle);
}



# // find the torsion angle, in radian between v1 and v2 with respect to axis u
# // according to the right hand rule
sub geometry_torsion { #C++ Type torsion(Vec3<Type> v1, Vec3<Type> v2, Vec3<Type> u)

	my $v1 = shift;
	my $v2 = shift;
	my $u = shift;

	$u = Vec3_normalize( $u );			# u.normalize();
	$v1 = Vec3_cross( $v1, $u );			# v1 = v1.cross(u);
	$v2 = Vec3_cross( $v2, $u );			# v2 = v2.cross(u);

	my $a = geometry_angle( $v1, $v2 );			# double a = angle(v1, v2);

	# // ensures that right hand rule is obeyed
	my $cross = Vec3_cross( $v1, $v2 );		# Vec3<Type> cross = v1.cross(v2);
	$cross = Vec3_normalize( $cross );		# cross.normalize();

	$cross = Vec3_subtract_Vec3( $cross, $u );	# cross -= u;

	my $cross_length = Vec3_length( $cross );
	if ($cross_length > 1) {			# if (cross.length() > 1) 
		$a = $a * -1;				# a = -a;
	}

	return $a;					# return a;

}



# // transformation Matrix representing the transformation relationship between
# // model 3 points to the actual 3 points
sub geometry_xformMat { #C++ Mat44<Type> xformMat(const Vec3<Type>& from_p1, const Vec3<Type>& from_p2, const Vec3<Type>& from_p3,
			#			const Vec3<Type>& p1, const Vec3<Type>& p2, const Vec3<Type>& p3)

	my $from_p1 = shift; #C++ const Vec3<Type>&
	my $from_p2 = shift; #C++ const Vec3<Type>&
	my $from_p3 = shift; #C++ const Vec3<Type>&
	my $p1 = shift; #C++ const Vec3<Type>&
	my $p2 = shift; #C++ const Vec3<Type>&
	my $p3 = shift; #C++ const Vec3<Type>&

	my $p1_p2 = Vec3_subtract_Vec3( $p2, $p1 );			# Vec3<Type> p1_p2 = p2 - p1;
	my $from_p1_p2 = Vec3_subtract_Vec3( $from_p2, $from_p1 );	# Vec3<Type> from_p1_p2 = from_p2 - from_p1;

	my $trans = Mat44_new();					# Mat44<Type> trans;
	my $param1 = Vec3_multiply_number( $from_p1, -1 );
	$trans = mat44impl_translate( $trans, $param1 );		# trans.translate(-from_p1);

	my $axis = Vec3_cross( $from_p1_p2, $p1_p2 );			# Vec3<Type> axis = from_p1_p2.cross(p1_p2);

	my $m = Mat44_new();						# Mat44<Type> m;
	my $param2 = geometry_torsion( $from_p1_p2, $p1_p2, $axis );
	$m = Mat44_rotate( $m, $param2, $axis );			# m.rotate(torsion(from_p1_p2, p1_p2, axis), axis);
	$trans = Mat44_multiply_Mat44( $m, $trans );			# trans = m*trans;

	my $from_p1_p3 = Mat44_multVec3d( $trans, $from_p3 );		# Vec3<Type> from_p1_p3 = trans.multVec3d(from_p3);

	my $p1_p3 = Vec3_subtract_Vec3( $p3, $p1 );			# Vec3<Type> p1_p3 = p3 - p1;
	$p1_p2 = Vec3_normalize( $p1_p2 );				# p1_p2.normalize();
	my $param3 = geometry_torsion( $from_p1_p3, $p1_p3, $p1_p2 );
	$m = Mat44_rotate( $m, $param3, $p1_p2 );			# m.rotate(torsion(from_p1_p3, p1_p3, p1_p2), p1_p2);
	$trans = Mat44_multiply_Mat44( $m, $trans );			# trans = m*trans;

	$m = mat44impl_translate( $m, $p1 );				# m.translate(p1);
	$trans = Mat44_multiply_Mat44( $m, $trans );			# trans = m*trans;

	return $trans;							# return trans;
}



sub get_2nd_struc_info {

	$input_tms_line = '';
	$input_inside_line = '';
	$input_outside_line = '';
	$input_TMS_line = '';
	$input_TMS_HELIX_line = '';
	$input_HELIX_line = '';
	$input_HELIX_INSIDE_line = '';
	$input_HELIX_OUTSIDE_line = '';
	$input_HELIX_MEMBRANE_line = '';

	if (defined($input_tms_file)) {
		my $found_this_fasta_id = 0;
		my $i = 0;
		while (($i < @input_tms_array_fasta_id) && ($found_this_fasta_id == 0)) {
			if ($input_fasta_id eq $input_tms_array_fasta_id[$i]) {
				$found_this_fasta_id = 1;
				if (exists($input_tms_array_tms_line[$i])) {
					if ($input_tms_array_tms_line[$i] ne '-') {
						$input_tms_line = $input_tms_array_tms_line[$i];
						my $ln = $program_options->{'lnOpt_index'}->{'TRUE-TMS'};
						$program_options->{'lnOpt'}->[$ln]->{'showGraph'} = 1;
					}
				}
				if (exists($input_tms_array_inside_line[$i])) {
					if ($input_tms_array_inside_line[$i] ne '-') {
						$input_inside_line = $input_tms_array_inside_line[$i];
						my $ln = $program_options->{'lnOpt_index'}->{'TRUE-INSIDE'};
						$program_options->{'lnOpt'}->[$ln]->{'showGraph'} = 1;
					}
				}
				if (exists($input_tms_array_outside_line[$i])) {
					if ($input_tms_array_outside_line[$i] ne '-') {
						$input_outside_line = $input_tms_array_outside_line[$i];
						my $ln = $program_options->{'lnOpt_index'}->{'TRUE-OUTSIDE'};
						$program_options->{'lnOpt'}->[$ln]->{'showGraph'} = 1;
					}
				}
				if (exists($input_tms_array_TMS_line[$i])) {
					if ($input_tms_array_TMS_line[$i] ne '-') {
						$input_TMS_line = $input_tms_array_TMS_line[$i];
					}
				}
				if (exists($input_tms_array_TMS_HELIX_line[$i])) {
					if ($input_tms_array_TMS_HELIX_line[$i] ne '-') {
						$input_TMS_HELIX_line = $input_tms_array_TMS_HELIX_line[$i];
					}
				}
				if (exists($input_tms_array_HELIX_line[$i])) {
					if ($input_tms_array_HELIX_line[$i] ne '-') {
						$input_HELIX_line = $input_tms_array_HELIX_line[$i];
					}
				}
				if (exists($input_tms_array_HELIX_INSIDE_line[$i])) {
					if ($input_tms_array_HELIX_INSIDE_line[$i] ne '-') {
						$input_HELIX_INSIDE_line = $input_tms_array_HELIX_INSIDE_line[$i];
					}
				}
				if (exists($input_tms_array_HELIX_OUTSIDE_line[$i])) {
					if ($input_tms_array_HELIX_OUTSIDE_line[$i] ne '-') {
						$input_HELIX_OUTSIDE_line = $input_tms_array_HELIX_OUTSIDE_line[$i];
					}
				}
				if (exists($input_tms_array_HELIX_MEMBRANE_line[$i])) {
					if ($input_tms_array_HELIX_MEMBRANE_line[$i] ne '-') {
						$input_HELIX_MEMBRANE_line = $input_tms_array_HELIX_MEMBRANE_line[$i];
					}
				}
				# if no TMS information was provided (probably manually created for displaying on graph),
				# then see if the TMS_HELIX was provided
				# (from running a program on the OPM Orientations of Proteins in Membrane database PDB files to find observed helices and their inside/membrane/outside position),
				if (($input_tms_line eq '-') or ($input_tms_line eq '')) {
					if ($input_TMS_HELIX_line ne '-') {
						$input_tms_line = $input_TMS_HELIX_line;
					}
				}
			} else {
				$i++;
			}
		}
	}
}



sub get_amphipathic_helices {

	$amphipathic_helices = (); # initialize this data structure for each new sequence

	# my $st_dot_chains_dot_size = $#{$structure->{'chains'}} + 1;
	# for ( my $c = 0; $c < $st_dot_chains_dot_size; $c++ ) {	 # this is only called in command-line mode, so should only have 1 chain
	my $c = 0;
	$amphipathic_helices = $amphipathic_helices_chains->[$c];

	return;
}



sub get_cgi_script_program_options {

	set_program_options_defaults();
	# read_program_options_defaults_file(); # defaults file from windows version
	get_html_program_options();
	adjust_options();
}



sub get_command_line_program_options {

	set_program_options_defaults();
	# read_program_options_defaults_file(); # defaults file from windows version
	read_program_options_input_file();
	adjust_options();
}



sub get_html_aaseq {

	my $params = $q->Vars;
	my @input_lines_1 = trim($params->{'search_text'});
	my $input_lines_2 = '';
	foreach my $input_line_1 (@input_lines_1) {
		$input_lines_2 .= $input_line_1;
	}
	my @input_lines_3 = split( /\n/, $input_lines_2 );
	my @aa_seq_chains;
	my $aaseq = '';
	my $i = 0;
	foreach my $input_line_3 (@input_lines_3) {
		if (substr($input_line_3,0,1) eq '>') {
			$aaseq =~ s/\W//g;
			$aaseq =~ s/\_//g;
			if ($aaseq ne '') {
				$aa_seq_chains[$i] = $aaseq;
				$i++;
				$aaseq = '';
			}
			if (($input_fasta_id eq '') or ($input_fasta_id eq '-')) {
				$input_fasta_id = trim($input_line_3);
			}
		} else {
			chomp($input_line_3);
			$aaseq .= trim($input_line_3);
		}
	}
	$aaseq =~ s/\W//g;
	$aaseq =~ s/\_//g;
	$aaseq =~ s/\s//g;
	if ($aaseq ne '') {
		$aa_seq_chains[$i] = $aaseq;
		$i++;
	}
	$input_aa_sequence .= $aaseq;

	my $ref_aa_seq_chains = \@aa_seq_chains;
	return $ref_aa_seq_chains;
}



sub get_html_program_options {

	my $params = $q->Vars;

	if (defined($params->{'option_SASARes'})) {
		$program_options->{'SASARes'} = trim($params->{'option_SASARes'});
	}
	if (defined($params->{'option_graph_residue_width'})) {
		$program_options->{'graph_residue_width'} = trim($params->{'option_graph_residue_width'});
	}
	if (defined($params->{'option_graph_height'})) {
		$program_options->{'graph_height'} = trim($params->{'option_graph_height'});
	}
	if (defined($params->{'option_valpred1_or_valpred2'})) {
		$TMPred2D->{'valpred1_or_valpred2'} = trim($params->{'option_valpred1_or_valpred2'});
	}
	if (defined($params->{'option_predict_amphipathic_helices'})) {
		$TMPred2D->{'predict_amphipathic_helices'} = trim($params->{'option_predict_amphipathic_helices'});
	}
	if (defined($params->{'option_mveAve'})) {
		$TMPred2D->{'mveAve'} = trim($params->{'option_mveAve'});
	}
	if (defined($params->{'option_mov_avg_1'})) {
		$TMPred2D->{'mov_avg_1'} = trim($params->{'option_mov_avg_1'});
	}
	if (defined($params->{'option_mov_avg_2'})) {
		$TMPred2D->{'mov_avg_2'} = trim($params->{'option_mov_avg_2'});
	}
	if (defined($params->{'option_mov_avg_3'})) {
		$TMPred2D->{'mov_avg_3'} = trim($params->{'option_mov_avg_3'});
	}
	if (defined($params->{'option_mov_avg_4'})) {
		$TMPred2D->{'mov_avg_4'} = trim($params->{'option_mov_avg_4'});
	}
	if (defined($params->{'option_minAveSasa'})) {
		$TMPred2D->{'minAveSasa'} = trim($params->{'option_minAveSasa'});
	}
	if (defined($params->{'option_TMLengthMin_for_valpred1'})) {
		$TMPred2D->{'TMLengthMin_for_valpred1'} = trim($params->{'option_TMLengthMin_for_valpred1'});
	}
	if (defined($params->{'option_TMLengthMax_for_valpred1'})) {
		$TMPred2D->{'TMLengthMax_for_valpred1'} = trim($params->{'option_TMLengthMax_for_valpred1'});
	}
	if (defined($params->{'option_TMLengthMin_for_valpred2'})) {
		$TMPred2D->{'TMLengthMin_for_valpred2'} = trim($params->{'option_TMLengthMin_for_valpred2'});
	}
	if (defined($params->{'option_TMLengthMax_for_valpred2'})) {
		$TMPred2D->{'TMLengthMax_for_valpred2'} = trim($params->{'option_TMLengthMax_for_valpred2'});
	}
	if (defined($params->{'option_aveRepRangeMin'})) {
		$TMPred2D->{'aveRepRangeMin'} = trim($params->{'option_aveRepRangeMin'});
	}
	if (defined($params->{'option_aveRepRangeMax'})) {
		$TMPred2D->{'aveRepRangeMax'} = trim($params->{'option_aveRepRangeMax'});
	}
	if (defined($params->{'option_aveProfRangeMin'})) {
		$TMPred2D->{'aveProfRangeMin'} = trim($params->{'option_aveProfRangeMin'});
	}
	if (defined($params->{'option_aveProfRangeMax'})) {
		$TMPred2D->{'aveProfRangeMax'} = trim($params->{'option_aveProfRangeMax'});
	}
	if (defined($params->{'option_lenDivAreaMin'})) {
		$TMPred2D->{'lenDivAreaMin'} = trim($params->{'option_lenDivAreaMin'});
	}
	if (defined($params->{'option_lenDivAreaMax'})) {
		$TMPred2D->{'lenDivAreaMax'} = trim($params->{'option_lenDivAreaMax'});
	}
	if (defined($params->{'option_minScoreDiff_for_valpred1'})) {
		$TMPred2D->{'minScoreDiff_for_valpred1'} = trim($params->{'option_minScoreDiff_for_valpred1'});
	}
	if (defined($params->{'option_minScoreDiff_movavg1_for_valpred2'})) {
		$TMPred2D->{'minScoreDiff_movavg1_for_valpred2'} = trim($params->{'option_minScoreDiff_movavg1_for_valpred2'});
	}
	if (defined($params->{'option_minScoreDiff_movavg2_for_valpred2'})) {
		$TMPred2D->{'minScoreDiff_movavg2_for_valpred2'} = trim($params->{'option_minScoreDiff_movavg2_for_valpred2'});
	}
	if (defined($params->{'option_min_avg_area_movavg1'})) {
		$TMPred2D->{'min_avg_area_movavg1'} = trim($params->{'option_min_avg_area_movavg1'});
	}
	if (defined($params->{'option_min_avg_area_movavg2'})) {
		$TMPred2D->{'min_avg_area_movavg2'} = trim($params->{'option_min_avg_area_movavg2'});
	}
	if (defined($params->{'option_max_nonTMS_length'})) {
		$TMPred2D->{'max_nonTMS_length'} = trim($params->{'option_max_nonTMS_length'});
	}
	if (defined($params->{'option_min_nonTMS_length'})) {
		$TMPred2D->{'min_nonTMS_length'} = trim($params->{'option_min_nonTMS_length'});
	}
	if (defined($params->{'option_min_nonTMS_score_diff'})) {
		$TMPred2D->{'min_nonTMS_score_diff'} = trim($params->{'option_min_nonTMS_score_diff'});
	}
	if (defined($params->{'option_min_nonTMS_area'})) {
		$TMPred2D->{'min_nonTMS_area'} = trim($params->{'option_min_nonTMS_area'});
	}
	if (defined($params->{'option_twilight_area_per_residue_lower_limit'})) {
		$TMPred2D->{'twilight_area_per_residue_lower_limit'} = trim($params->{'option_twilight_area_per_residue_lower_limit'});
	}
	if (defined($params->{'option_twilight_area_per_residue_upper_limit'})) {
		$TMPred2D->{'twilight_area_per_residue_upper_limit'} = trim($params->{'option_twilight_area_per_residue_upper_limit'});
	}
	if (defined($params->{'option_ignore_twilight_area'})) {
		$TMPred2D->{'ignore_twilight_area'} = trim($params->{'option_ignore_twilight_area'});
	}
	if (defined($params->{'option_search_area_for_neighbour_tms'})) {
		$TMPred2D->{'search_area_for_neighbour_tms'} = trim($params->{'option_search_area_for_neighbour_tms'});
	}
	if (defined($params->{'option_helix_length_min'})) {
		$TMPred2D->{'helix_length_min'} = trim($params->{'option_helix_length_min'});
	}
	if (defined($params->{'option_min_score_diff_for_TMS_ends'})) {
		$TMPred2D->{'min_score_diff_for_TMS_ends'} = trim($params->{'option_min_score_diff_for_TMS_ends'});
	}

	my @option_lines = ('PROFILES3D-LINE','REPIMPS-LINE','BOTH-LINE','PROFILES3D-MOV-AVG-1','REPIMPS-MOV-AVG-1','TMS-MOV-AVG-1',
		'PROFILES3D-MOV-AVG-2','REPIMPS-MOV-AVG-2','TMS-MOV-AVG-2','PROFILES3D-MOV-AVG-3','REPIMPS-MOV-AVG-3','HYDROPHILIC-MOV-AVG-3',
		'PROFILES3D-MOV-AVG-4','REPIMPS-MOV-AVG-4','HYDROPHILIC-MOV-AVG-4',
		'TMS');
		# 'TMS','TRUE-TMS','TRUE-INSIDE','TRUE-OUTSIDE');
	my @option_attributes = ('showGraph','colour_name','legend_name','smooth');

	for ( my $o = 0; $o < @option_lines; $o++ ) {

		my $option_name = $option_lines[$o];
		my $ln = $program_options->{'lnOpt_index'}->{$option_name};
		my $html_display_option_name = $program_options->{'lnOpt'}->[$ln]->{'legend_name'} . ' :&nbsp;';

		for ( my $n = 0; $n < @option_attributes; $n++ ) {

			my $option_attribute = $option_attributes[$n];
			my $param_name = 'option_' . $ln . '_' . $option_attribute;

			if (defined($params->{$param_name})) {
				my $param_value = trim($params->{$param_name});
				my $param_ok = 0;
				if ($n eq 'showGraph') {
					if (($param_value == 0) || ($param_value == 1)) { # is it boolean?
						$param_ok = 1;
					}
				} elsif ($n eq 'weight') {
					if ($param_value =~ /^(\d+\.?\d*|\.\d+)$/) { # is it a +ve numeric?
						$param_ok = 1;
					}
				} elsif ($n eq 'mov_avg') {
					if ($param_value =~ /^(\d+)$/) { # is it a +ve integer?
						$param_ok = 1;
					}
				} else {
					$param_ok = 1;
				}
				if ($param_ok == 1) {
					$program_options->{'lnOpt'}->[$ln]->{$option_attribute} = $param_value;
				}
			}
		}
	}

	if ( (defined($params->{'option_show_all_valpred2_lines'})) && (defined($params->{'option_valpred1_or_valpred2'})) ) {
		if ( ($params->{'option_show_all_valpred2_lines'} == 1) && ($params->{'option_valpred1_or_valpred2'} == 2) ) {
			my @option_lines = ('PROFILES3D-LINE','REPIMPS-LINE','PROFILES3D-MOV-AVG-1','REPIMPS-MOV-AVG-1','TMS-MOV-AVG-1',
				'PROFILES3D-MOV-AVG-2','REPIMPS-MOV-AVG-2','TMS-MOV-AVG-2','PROFILES3D-MOV-AVG-3','REPIMPS-MOV-AVG-3','HYDROPHILIC-MOV-AVG-3',
				'PROFILES3D-MOV-AVG-4','REPIMPS-MOV-AVG-4','HYDROPHILIC-MOV-AVG-4',
				'TMS');
			for ( my $o = 0; $o < @option_lines; $o++ ) {
				my $option_name = $option_lines[$o];
				my $ln = $program_options->{'lnOpt_index'}->{$option_name};
				$program_options->{'lnOpt'}->[$ln]->{'showGraph'} = 1;
			}
		}
	}
}



sub hard_coded_amino_acid_pdb {

	my $amino_acid = shift;
	$amino_acid = uc $amino_acid;
	my @pdb_lines;

	if ($amino_acid eq 'ALA') {

		@pdb_lines = (	'ATOM      1   N  ALA     1       1.458   0.000   0.000  1.00  0.00',
				'ATOM      2  CA  ALA     1       0.000   0.000   0.000  1.00  0.00',
				'ATOM      3   C  ALA     1      -0.551   0.000  -1.422  1.00  0.00',
				'ATOM      3   O  ALA     1      -1.535  -0.689  -1.687  1.00  0.00',
				'ATOM      5  CB  ALA     1      -0.540   1.211   0.775  1.00  0.00',
				'TER');

	} elsif ($amino_acid eq 'ARG') {

		@pdb_lines = (	'ATOM      1   N  ARG     1       1.458   0.000   0.000  1.00  0.00',
				'ATOM      2  CA  ARG     1       0.000   0.000   0.000  1.00  0.00',
				'ATOM      3   C  ARG     1      -0.551   0.000  -1.422  1.00  0.00',
				'ATOM      3   O  ARG     1      -1.535  -0.689  -1.687  1.00  0.00',
				'ATOM      5  CB  ARG     1      -0.536   1.207   0.772  1.00  0.00',
				'ATOM      6  CG  ARG     1      -0.246   1.169   2.263  1.00  0.00',
				'ATOM      7  CD  ARG     1      -1.062   2.210   3.012  1.00  0.00',
				'ATOM      8  NE  ARG     1      -0.692   3.571   2.633  1.00  0.00',
				'ATOM      9  CZ  ARG     1       0.313   4.249   3.176  1.00  0.00',
				'ATOM     10  NH1 ARG     1       1.054   3.693   4.124  1.00  0.00',
				'ATOM     11  NH2 ARG     1       0.575   5.484   2.768  1.00  0.00',
				'TER ');

	} elsif ($amino_acid eq 'ASN') {

		@pdb_lines = (	'ATOM      1   N  ASN     1       1.458   0.000   0.000  1.00  0.00',
				'ATOM      2  CA  ASN     1       0.000   0.000   0.000  1.00  0.00',
				'ATOM      3   C  ASN     1      -0.551   0.000  -1.422  1.00  0.00',
				'ATOM      3   O  ASN     1      -1.535  -0.689  -1.687  1.00  0.00',
				'ATOM      5  CB  ASN     1      -0.536   1.207   0.772  1.00  0.00',
				'ATOM      6  CG  ASN     1      -0.292   1.100   2.264  1.00  0.00',
				'ATOM      7  OD1 ASN     1      -0.082   0.009   2.791  1.00  0.00',
				'ATOM      8  ND2 ASN     1      -0.317   2.236   2.949  1.00  0.00',
				'TER'); 

	} elsif ($amino_acid eq 'ASP') {

		@pdb_lines = (	'ATOM      1   N  ASP     1       1.458   0.000   0.000  1.00  0.00',
				'ATOM      2  CA  ASP     1       0.000   0.000   0.000  1.00  0.00',
				'ATOM      3   C  ASP     1      -0.551   0.000  -1.422  1.00  0.00',
				'ATOM      3   O  ASP     1      -1.535  -0.689  -1.687  1.00  0.00',
				'ATOM      5  CB  ASP     1      -0.536   1.207   0.772  1.00  0.00',
				'ATOM      6  CG  ASP     1      -0.292   1.100   2.264  1.00  0.00',
				'ATOM      7  OD1 ASP     1       0.119   0.013   2.723  1.00  0.00',
				'ATOM      8  OD2 ASP     1      -0.510   2.103   2.976  1.00  0.00',
				'TER');

	} elsif ($amino_acid eq 'CYS') {

		@pdb_lines = (	'ATOM      1   N  CYS     1       1.458   0.000   0.000  1.00  0.00',
				'ATOM      2  CA  CYS     1       0.000   0.000   0.000  1.00  0.00',
				'ATOM      3   C  CYS     1      -0.551   0.000  -1.422  1.00  0.00',
				'ATOM      3   O  CYS     1      -1.535  -0.689  -1.687  1.00  0.00',
				'ATOM      5  CB  CYS     1      -0.536   1.207   0.772  1.00  0.00',
				'ATOM      6  SG  CYS     1      -0.143   1.198   2.551  1.00  0.00',
				'TER');

	} elsif ($amino_acid eq 'GLN') {

		@pdb_lines = (	'ATOM      1   N  GLN     1       1.458   0.000   0.000  1.00  0.00',
				'ATOM      2  CA  GLN     1       0.000   0.000   0.000  1.00  0.00',
				'ATOM      3   C  GLN     1      -0.551   0.000  -1.422  1.00  0.00',
				'ATOM      3   O  GLN     1      -1.535  -0.689  -1.687  1.00  0.00',
				'ATOM      5  CB  GLN     1      -0.536   1.207   0.772  1.00  0.00',
				'ATOM      6  CG  GLN     1      -0.246   1.168   2.264  1.00  0.00',
				'ATOM      7  CD  GLN     1      -0.791   2.380   2.994  1.00  0.00',
				'ATOM      8  OE1 GLN     1      -0.973   3.444   2.403  1.00  0.00',
				'ATOM      9  NE2 GLN     1      -1.055   2.221   4.286  1.00  0.00',
				'TER');

	} elsif ($amino_acid eq 'GLU') {

		@pdb_lines = (	'ATOM      1   N  GLU     1       1.458   0.000   0.000  1.00  0.00',
				'ATOM      2  CA  GLU     1       0.000   0.000   0.000  1.00  0.00',
				'ATOM      3   C  GLU     1      -0.551   0.000  -1.422  1.00  0.00',
				'ATOM      3   O  GLU     1      -1.535  -0.689  -1.687  1.00  0.00',
				'ATOM      5  CB  GLU     1      -0.536   1.207   0.772  1.00  0.00',
				'ATOM      6  CG  GLU     1      -0.246   1.168   2.264  1.00  0.00',
				'ATOM      7  CD  GLU     1      -0.791   2.380   2.994  1.00  0.00',
				'ATOM      8  OE1 GLU     1      -1.225   3.337   2.319  1.00  0.00',
				'ATOM      9  OE2 GLU     1      -0.785   2.373   4.244  1.00  0.00',
				'TER');

	} elsif ($amino_acid eq 'GLY') {

		@pdb_lines = (	'ATOM      1   N  GLY     1       1.458   0.000   0.000  1.00  0.00',
				'ATOM      2  CA  GLY     1       0.000   0.000   0.000  1.00  0.00',
				'ATOM      3   C  GLY     1      -0.551   0.000  -1.422  1.00  0.00',
				'ATOM      3   O  GLY     1      -1.535  -0.689  -1.687  1.00  0.00',
				'TER');

	} elsif ($amino_acid eq 'HIS') {

		@pdb_lines = (	'ATOM      1   N  HIS     1       1.458   0.000   0.000  1.00  0.00',
				'ATOM      2  CA  HIS     1       0.000   0.000   0.000  1.00  0.00',
				'ATOM      3   C  HIS     1      -0.551   0.000  -1.422  1.00  0.00',
				'ATOM      3   O  HIS     1      -1.535  -0.689  -1.687  1.00  0.00',
				'ATOM      5  CB  HIS     1      -0.536   1.207   0.772  1.00  0.00',
				'ATOM      6  CG  HIS     1      -0.205   1.187   2.232  1.00  0.00',
				'ATOM      7  ND1 HIS     1      -0.760   0.280   3.108  1.00  0.00',
				'ATOM      8  CD2 HIS     1       0.624   1.962   2.970  1.00  0.00',
				'ATOM      9  CE1 HIS     1      -0.286   0.500   4.322  1.00  0.00',
				'ATOM     10  NE2 HIS     1       0.555   1.512   4.266  1.00  0.00',
				'TER');

	} elsif ($amino_acid eq 'ILE') {

		@pdb_lines = (	'ATOM      1   N  ILE     1       1.458   0.000   0.000  1.00  0.00',
				'ATOM      2  CA  ILE     1       0.000   0.000   0.000  1.00  0.00',
				'ATOM      3   C  ILE     1      -0.551   0.000  -1.422  1.00  0.00',
				'ATOM      3   O  ILE     1      -1.535  -0.689  -1.687  1.00  0.00',
				'ATOM      5  CB  ILE     1      -0.564   1.215   0.759  1.00  0.00',
				'ATOM      6  CG1 ILE     1      -0.195   1.136   2.242  1.00  0.00',
				'ATOM      7  CG2 ILE     1      -2.080   1.261   0.641  1.00  0.00',
				'ATOM      8  CD1 ILE     1      -0.500   2.399   3.016  1.00  0.00',
				'TER');

	} elsif ($amino_acid eq 'LEU') {

		@pdb_lines = (	'ATOM      1   N  LEU     1       1.458   0.000   0.000  1.00  0.00',
				'ATOM      2  CA  LEU     1       0.000   0.000   0.000  1.00  0.00',
				'ATOM      3   C  LEU     1      -0.551   0.000  -1.422  1.00  0.00',
				'ATOM      3   O  LEU     1      -1.535  -0.689  -1.687  1.00  0.00',
				'ATOM      5  CB  LEU     1      -0.536   1.207   0.772  1.00  0.00',
				'ATOM      6  CG  LEU     1      -0.230   1.244   2.271  1.00  0.00',
				'ATOM      7  CD1 LEU     1      -0.690   2.559   2.881  1.00  0.00',
				'ATOM      8  CD2 LEU     1      -0.943   0.111   2.992  1.00  0.00',
				'TER');

	} elsif ($amino_acid eq 'LYS') {

		@pdb_lines = (	'ATOM      1   N  LYS     1       1.458   0.000   0.000  1.00  0.00',
				'ATOM      2  CA  LYS     1       0.000   0.000   0.000  1.00  0.00',
				'ATOM      3   C  LYS     1      -0.551   0.000  -1.422  1.00  0.00',
				'ATOM      3   O  LYS     1      -1.535  -0.689  -1.687  1.00  0.00',
				'ATOM      5  CB  LYS     1      -0.536   1.207   0.772  1.00  0.00',
				'ATOM      6  CG  LYS     1      -0.246   1.169   2.263  1.00  0.00',
				'ATOM      7  CD  LYS     1      -0.805   2.395   2.967  1.00  0.00',
				'ATOM      8  CE  LYS     1      -0.514   2.357   4.458  1.00  0.00',
				'ATOM      9  NZ  LYS     1      -1.057   3.553   5.160  1.00  0.00',
				'TER');

	} elsif ($amino_acid eq 'MET') {

		@pdb_lines = (	'ATOM      1   N  MET     1       1.458   0.000   0.000  1.00  0.00',
				'ATOM      2  CA  MET     1       0.000   0.000   0.000  1.00  0.00',
				'ATOM      3   C  MET     1      -0.551   0.000  -1.422  1.00  0.00',
				'ATOM      3   O  MET     1      -1.535  -0.689  -1.687  1.00  0.00',
				'ATOM      5  CB  MET     1      -0.536   1.207   0.772  1.00  0.00',
				'ATOM      6  CG  MET     1      -0.246   1.168   2.264  1.00  0.00',
				'ATOM      7  SD  MET     1      -0.893   2.607   3.136  1.00  0.00',
				'ATOM      8  CE  MET     1      -2.643   2.222   3.171  1.00  0.00',
				'TER');

	} elsif ($amino_acid eq 'PHE') {

		@pdb_lines = (	'ATOM      1   N  PHE     1       1.458   0.000   0.000  1.00  0.00',
				'ATOM      2  CA  PHE     1       0.000   0.000   0.000  1.00  0.00',
				'ATOM      3   C  PHE     1      -0.551   0.000  -1.422  1.00  0.00',
				'ATOM      3   O  PHE     1      -1.535  -0.689  -1.687  1.00  0.00',
				'ATOM      5  CB  PHE     1      -0.536   1.207   0.772  1.00  0.00',
				'ATOM      6  CG  PHE     1      -2.033   1.242   0.879  1.00  0.00',
				'ATOM      7  CD1 PHE     1      -2.693   0.470   1.820  1.00  0.00',
				'ATOM      8  CD2 PHE     1      -2.784   2.046   0.041  1.00  0.00',
				'ATOM      9  CE1 PHE     1      -4.071   0.501   1.919  1.00  0.00',
				'ATOM     10  CE2 PHE     1      -4.162   2.079   0.136  1.00  0.00',
				'ATOM     11  CZ  PHE     1      -4.806   1.307   1.077  1.00  0.00',
				'TER');

	} elsif ($amino_acid eq 'PRO') {

		@pdb_lines = (	'ATOM      1   N  PRO     1       1.458   0.000   0.000  1.00  0.00',
				'ATOM      2  CA  PRO     1       0.000   0.000   0.000  1.00  0.00',
				'ATOM      3   C  PRO     1      -0.551   0.000  -1.422  1.00  0.00',
				'ATOM      3   O  PRO     1      -1.535  -0.689  -1.687  1.00  0.00',
				'ATOM      5  CB  PRO     1      -0.344   1.256   0.803  1.00  0.00',
				'ATOM      6  CG  PRO     1       0.835   1.435   1.733  1.00  0.00',
				'ATOM      7  CD  PRO     1       2.018   1.062   0.857  1.00  0.00',
				'TER');

	} elsif ($amino_acid eq 'SER') {

		@pdb_lines = (	'ATOM      1   N  SER     1       1.458   0.000   0.000  1.00  0.00',
				'ATOM      2  CA  SER     1       0.000   0.000   0.000  1.00  0.00',
				'ATOM      3   C  SER     1      -0.551   0.000  -1.422  1.00  0.00',
				'ATOM      3   O  SER     1      -1.535  -0.689  -1.687  1.00  0.00',
				'ATOM      5  CB  SER     1      -0.536   1.207   0.772  1.00  0.00',
				'ATOM      6  OG  SER     1      -0.192   1.129   2.144  1.00  0.00',
				'TER');

	} elsif ($amino_acid eq 'THR') {

		@pdb_lines = (	'ATOM      1   N  THR     1       1.458   0.000   0.000  1.00  0.00',
				'ATOM      2  CA  THR     1       0.000   0.000   0.000  1.00  0.00',
				'ATOM      3   C  THR     1      -0.551   0.000  -1.422  1.00  0.00',
				'ATOM      3   O  THR     1      -1.535  -0.689  -1.687  1.00  0.00',
				'ATOM      5  CB  THR     1      -0.564   1.215   0.759  1.00  0.00',
				'ATOM      6  OG1 THR     1      -0.209   1.123   2.145  1.00  0.00',
				'ATOM      7  CG2 THR     1      -2.080   1.261   0.641  1.00  0.00',
				'TER');

	} elsif ($amino_acid eq 'TRP') {

		@pdb_lines = (	'ATOM      1   N  TRP     1       1.458   0.000   0.000  1.00  0.00',
				'ATOM      2  CA  TRP     1       0.000   0.000   0.000  1.00  0.00',
				'ATOM      3   C  TRP     1      -0.551   0.000  -1.422  1.00  0.00',
				'ATOM      3   O  TRP     1      -1.535  -0.689  -1.687  1.00  0.00',
				'ATOM      5  CB  TRP     1      -0.536   1.207   0.772  1.00  0.00',
				'ATOM      6  CG  TRP     1      -2.030   1.237   0.876  1.00  0.00',
				'ATOM      7  CD1 TRP     1      -2.794   0.683   1.862  1.00  0.00',
				'ATOM      8  CD2 TRP     1      -2.940   1.855  -0.042  1.00  0.00',
				'ATOM      9  NE1 TRP     1      -4.126   0.916   1.618  1.00  0.00',
				'ATOM     10  CE2 TRP     1      -4.241   1.635   0.453  1.00  0.00',
				'ATOM     11  CE3 TRP     1      -2.781   2.571  -1.231  1.00  0.00',
				'ATOM     12  CZ2 TRP     1      -5.379   2.105  -0.202  1.00  0.00',
				'ATOM     13  CZ3 TRP     1      -3.909   3.037  -1.879  1.00  0.00',
				'ATOM     14  CH2 TRP     1      -5.196   2.804  -1.365  1.00  0.00',
				'TER');

	} elsif ($amino_acid eq 'TYR') {

		@pdb_lines = (	'ATOM      1   N  TYR     1       1.458   0.000   0.000  1.00  0.00',
				'ATOM      2  CA  TYR     1       0.000   0.000   0.000  1.00  0.00',
				'ATOM      3   C  TYR     1      -0.551   0.000  -1.422  1.00  0.00',
				'ATOM      3   O  TYR     1      -1.535  -0.689  -1.687  1.00  0.00',
				'ATOM      5  CB  TYR     1      -0.536   1.207   0.772  1.00  0.00',
				'ATOM      6  CG  TYR     1      -2.044   1.244   0.882  1.00  0.00',
				'ATOM      7  CD1 TYR     1      -2.708   0.472   1.826  1.00  0.00',
				'ATOM      8  CD2 TYR     1      -2.799   2.052   0.041  1.00  0.00',
				'ATOM      9  CE1 TYR     1      -4.084   0.499   1.932  1.00  0.00',
				'ATOM     10  CE2 TYR     1      -4.176   2.093   0.133  1.00  0.00',
				'ATOM     11  CZ  TYR     1      -4.818   1.316   1.080  1.00  0.00',
				'ATOM     12  OH  TYR     1      -6.189   1.349   1.180  1.00  0.00',
				'TER');

	} elsif ($amino_acid eq 'VAL') {

		@pdb_lines = (	'ATOM      1   N  VAL     1       1.458   0.000   0.000  1.00  0.00',
				'ATOM      2  CA  VAL     1       0.000   0.000   0.000  1.00  0.00',
				'ATOM      3   C  VAL     1      -0.551   0.000  -1.422  1.00  0.00',
				'ATOM      3   O  VAL     1      -1.535  -0.689  -1.687  1.00  0.00',
				'ATOM      5  CB  VAL     1      -0.564   1.215   0.759  1.00  0.00',
				'ATOM      6  CG1 VAL     1      -2.080   1.261   0.641  1.00  0.00',
				'ATOM      7  CG2 VAL     1      -0.199   1.138   2.234  1.00  0.00',
				'TER');
	}

	my $pdb_lines_ref = \@pdb_lines;
	return $pdb_lines_ref;
}



sub hard_coded_sphere_points {

	# hard coded from resolution = 2, file = ./Sphere_Points/maxvol.3.252.txt

	my $points_string = "";
	$points_string .= "252\n";
	$points_string .= "0.364416027907\n0.839162952267\n0.403740632272\n0.659046753775\n0.741815371905\n0.123965036782\n0.585553996265\n0.731343352236\n0.349661863232\n0.679618729744\n";
	$points_string .= "0.555611950811\n0.478971546438\n0.525731112119\n0.850650808352\n0.000000000000\n0.000000000000\n0.525731112119\n-0.850650808352\n-0.525731112119\n0.850650808352\n";
	$points_string .= "0.000000000000\n0.000000000000\n0.525731112119\n0.850650808352\n0.850650808352\n0.000000000000\n0.525731112119\n-0.850650808352\n0.000000000000\n0.525731112119\n";
	$points_string .= "0.000000000000\n-0.525731112119\n0.850650808352\n-0.850650808352\n0.000000000000\n-0.525731112119\n-0.525731112119\n-0.850650808352\n0.000000000000\n0.850650808352\n";
	$points_string .= "0.000000000000\n-0.525731112119\n0.525731112119\n-0.850650808352\n0.000000000000\n0.000000000000\n-0.525731112119\n-0.850650808352\n0.403740632174\n0.364416027652\n";
	$points_string .= "-0.839162952425\n0.000000001184\n0.334502077164\n-0.942395012918\n0.216102916071\n0.235892133096\n-0.947446268244\n0.219380291097\n0.000000000079\n-0.975639425135\n";
	$points_string .= "-0.364416027907\n0.839162952267\n-0.403740632272\n-0.659046753775\n0.741815371905\n-0.123965036782\n-0.585553996265\n0.731343352236\n-0.349661863232\n-0.679618729744\n";
	$points_string .= "0.555611950811\n-0.478971546438\n-0.403740632174\n0.364416027652\n0.839162952425\n-0.000000001184\n0.334502077164\n0.942395012918\n-0.216102916071\n0.235892133096\n";
	$points_string .= "0.947446268244\n-0.219380291097\n0.000000000079\n0.975639425135\n0.768156660082\n0.589637518944\n0.249525433225\n0.659046754960\n0.741815371173\n-0.123965034866\n";
	$points_string .= "0.801656912336\n0.597784405036\n0.000000000039\n0.898999020841\n0.420027474451\n-0.124006778982\n-0.114890594426\n0.993378151216\n-0.000000000158\n-0.458467112030\n";
	$points_string .= "0.865780407955\n-0.200579641014\n-0.235892133057\n0.947446268268\n-0.216102916008\n-0.124006779060\n0.898999020889\n-0.420027474324\n-0.249525433481\n0.768156659924\n";
	$points_string .= "0.589637519041\n-0.200579641746\n0.458467113946\n0.865780406770\n-0.349661863208\n0.585553996328\n0.731343352197\n-0.555611950684\n0.478971546516\n0.679618729793\n";
	$points_string .= "0.589637519041\n-0.249525433481\n0.768156659924\n0.865780406770\n-0.200579641746\n0.458467113946\n0.731343352197\n-0.349661863208\n0.585553996328\n0.679618729793\n";
	$points_string .= "-0.555611950684\n0.478971546516\n-0.993378151216\n-0.000000000158\n0.114890594426\n-0.865780407955\n-0.200579641014\n0.458467112030\n-0.947446268268\n-0.216102916008\n";
	$points_string .= "0.235892133057\n-0.898999020889\n-0.420027474324\n0.124006779060\n0.768156659924\n0.589637519041\n-0.249525433481\n0.458467113946\n0.865780406770\n-0.200579641746\n";
	$points_string .= "0.585553996328\n0.731343352197\n-0.349661863208\n0.478971546516\n0.679618729793\n-0.555611950684\n0.249525433225\n0.768156660082\n0.589637518944\n-0.123965034866\n";
	$points_string .= "0.659046754960\n0.741815371173\n0.000000000039\n0.801656912336\n0.597784405036\n-0.124006778982\n0.898999020841\n0.420027474451\n-0.403740632272\n-0.364416027907\n";
	$points_string .= "0.839162952267\n-0.123965036782\n-0.659046753775\n0.741815371905\n-0.349661863232\n-0.585553996265\n0.731343352236\n-0.478971546438\n-0.679618729744\n0.555611950811\n";
	$points_string .= "0.839162952425\n-0.403740632174\n0.364416027652\n0.942395012918\n-0.000000001184\n0.334502077164\n0.947446268244\n-0.216102916071\n0.235892133096\n0.975639425135\n";
	$points_string .= "-0.219380291097\n0.000000000079\n0.403740632272\n0.364416027907\n0.839162952267\n0.123965036782\n0.659046753775\n0.741815371905\n0.349661863232\n0.585553996265\n";
	$points_string .= "0.731343352236\n0.478971546438\n0.679618729744\n0.555611950811\n-0.000000000158\n-0.114890594426\n0.993378151216\n-0.200579641014\n-0.458467112030\n0.865780407955\n";
	$points_string .= "-0.216102916008\n-0.235892133057\n0.947446268268\n-0.420027474324\n-0.124006779060\n0.898999020889\n0.000000000158\n0.114890594426\n0.993378151216\n0.200579641014\n";
	$points_string .= "0.458467112030\n0.865780407955\n0.216102916008\n0.235892133057\n0.947446268268\n0.420027474324\n0.124006779060\n0.898999020889\n-0.364416027652\n0.839162952425\n";
	$points_string .= "0.403740632174\n-0.334502077164\n0.942395012918\n0.000000001184\n-0.235892133096\n0.947446268244\n0.216102916071\n-0.000000000079\n0.975639425135\n0.219380291097\n";
	$points_string .= "-0.589637518944\n0.249525433225\n-0.768156660082\n-0.741815371173\n-0.123965034866\n-0.659046754960\n-0.597784405036\n0.000000000039\n-0.801656912336\n-0.420027474451\n";
	$points_string .= "-0.124006778982\n-0.898999020841\n-0.839162952267\n-0.403740632272\n0.364416027907\n-0.741815371905\n-0.123965036782\n0.659046753775\n-0.731343352236\n-0.349661863232\n";
	$points_string .= "0.585553996265\n-0.555611950811\n-0.478971546438\n0.679618729744\n-0.768156659924\n0.589637519041\n0.249525433481\n-0.458467113946\n0.865780406770\n0.200579641746\n";
	$points_string .= "-0.585553996328\n0.731343352197\n0.349661863208\n-0.478971546516\n0.679618729793\n0.555611950684\n-0.839162952267\n0.403740632272\n-0.364416027907\n-0.741815371905\n";
	$points_string .= "0.123965036782\n-0.659046753775\n-0.731343352236\n0.349661863232\n-0.585553996265\n-0.555611950811\n0.478971546438\n-0.679618729744\n-0.364416027652\n-0.839162952425\n";
	$points_string .= "-0.403740632174\n-0.334502077164\n-0.942395012918\n-0.000000001184\n-0.235892133096\n-0.947446268244\n-0.216102916071\n-0.000000000079\n-0.975639425135\n-0.219380291097\n";
	$points_string .= "0.589637518944\n0.249525433225\n0.768156660082\n0.741815371173\n-0.123965034866\n0.659046754960\n0.597784405036\n0.000000000039\n0.801656912336\n0.420027474451\n";
	$points_string .= "-0.124006778982\n0.898999020841\n0.839162952267\n-0.403740632272\n-0.364416027907\n0.741815371905\n-0.123965036782\n-0.659046753775\n0.731343352236\n-0.349661863232\n";
	$points_string .= "-0.585553996265\n0.555611950811\n-0.478971546438\n-0.679618729744\n0.364416027652\n0.839162952425\n-0.403740632174\n0.334502077164\n0.942395012918\n-0.000000001184\n";
	$points_string .= "0.235892133096\n0.947446268244\n-0.216102916071\n0.000000000079\n0.975639425135\n-0.219380291097\n0.839162952267\n0.403740632272\n0.364416027907\n0.741815371905\n";
	$points_string .= "0.123965036782\n0.659046753775\n0.731343352236\n0.349661863232\n0.585553996265\n0.555611950811\n0.478971546438\n0.679618729744\n0.993378151216\n-0.000000000158\n";
	$points_string .= "-0.114890594426\n0.865780407955\n-0.200579641014\n-0.458467112030\n0.947446268268\n-0.216102916008\n-0.235892133057\n0.898999020889\n-0.420027474324\n-0.124006779060\n";
	$points_string .= "-0.000000000158\n0.114890594426\n-0.993378151216\n-0.200579641014\n0.458467112030\n-0.865780407955\n-0.216102916008\n0.235892133057\n-0.947446268268\n-0.420027474324\n";
	$points_string .= "0.124006779060\n-0.898999020889\n0.364416027652\n-0.839162952425\n0.403740632174\n0.334502077164\n-0.942395012918\n0.000000001184\n0.235892133096\n-0.947446268244\n";
	$points_string .= "0.216102916071\n0.000000000079\n-0.975639425135\n0.219380291097\n-0.839162952425\n0.403740632174\n0.364416027652\n-0.942395012918\n0.000000001184\n0.334502077164\n";
	$points_string .= "-0.947446268244\n0.216102916071\n0.235892133096\n-0.975639425135\n0.219380291097\n0.000000000079\n-0.249525433225\n-0.768156660082\n0.589637518944\n0.123965034866\n";
	$points_string .= "-0.659046754960\n0.741815371173\n-0.000000000039\n-0.801656912336\n0.597784405036\n0.124006778982\n-0.898999020841\n0.420027474451\n-0.589637519041\n0.249525433481\n";
	$points_string .= "0.768156659924\n-0.865780406770\n0.200579641746\n0.458467113946\n-0.731343352197\n0.349661863208\n0.585553996328\n-0.679618729793\n0.555611950684\n0.478971546516\n";
	$points_string .= "-0.768156659924\n-0.589637519041\n-0.249525433481\n-0.458467113946\n-0.865780406770\n-0.200579641746\n-0.585553996328\n-0.731343352197\n-0.349661863208\n-0.478971546516\n";
	$points_string .= "-0.679618729793\n-0.555611950684\n0.249525433481\n-0.768156659924\n0.589637519041\n0.200579641746\n-0.458467113946\n0.865780406770\n0.349661863208\n-0.585553996328\n";
	$points_string .= "0.731343352197\n0.555611950684\n-0.478971546516\n0.679618729793\n-0.589637518944\n-0.249525433225\n0.768156660082\n-0.741815371173\n0.123965034866\n0.659046754960\n";
	$points_string .= "-0.597784405036\n-0.000000000039\n0.801656912336\n-0.420027474451\n0.124006778982\n0.898999020841\n0.114890594426\n-0.993378151216\n-0.000000000158\n0.458467112030\n";
	$points_string .= "-0.865780407955\n-0.200579641014\n0.235892133057\n-0.947446268268\n-0.216102916008\n0.124006779060\n-0.898999020889\n-0.420027474324\n0.993378151216\n0.000000000158\n";
	$points_string .= "0.114890594426\n0.865780407955\n0.200579641014\n0.458467112030\n0.947446268268\n0.216102916008\n0.235892133057\n0.898999020889\n0.420027474324\n0.124006779060\n";
	$points_string .= "0.403740632174\n-0.364416027652\n0.839162952425\n0.000000001184\n-0.334502077164\n0.942395012918\n0.216102916071\n-0.235892133096\n0.947446268244\n0.219380291097\n";
	$points_string .= "-0.000000000079\n0.975639425135\n0.589637518944\n-0.249525433225\n-0.768156660082\n0.741815371173\n0.123965034866\n-0.659046754960\n0.597784405036\n-0.000000000039\n";
	$points_string .= "-0.801656912336\n0.420027474451\n0.124006778982\n-0.898999020841\n-0.768156660082\n-0.589637518944\n0.249525433225\n-0.659046754960\n-0.741815371173\n-0.123965034866\n";
	$points_string .= "-0.801656912336\n-0.597784405036\n0.000000000039\n-0.898999020841\n-0.420027474451\n-0.124006778982\n0.768156659924\n-0.589637519041\n0.249525433481\n0.458467113946\n";
	$points_string .= "-0.865780406770\n0.200579641746\n0.585553996328\n-0.731343352197\n0.349661863208\n0.478971546516\n-0.679618729793\n0.555611950684\n-0.249525433481\n-0.768156659924\n";
	$points_string .= "-0.589637519041\n-0.200579641746\n-0.458467113946\n-0.865780406770\n-0.349661863208\n-0.585553996328\n-0.731343352197\n-0.555611950684\n-0.478971546516\n-0.679618729793\n";
	$points_string .= "0.249525433481\n0.768156659924\n-0.589637519041\n0.200579641746\n0.458467113946\n-0.865780406770\n0.349661863208\n0.585553996328\n-0.731343352197\n0.555611950684\n";
	$points_string .= "0.478971546516\n-0.679618729793\n-0.589637519041\n-0.249525433481\n-0.768156659924\n-0.865780406770\n-0.200579641746\n-0.458467113946\n-0.731343352197\n-0.349661863208\n";
	$points_string .= "-0.585553996328\n-0.679618729793\n-0.555611950684\n-0.478971546516\n-0.768156660082\n0.589637518944\n-0.249525433225\n-0.659046754960\n0.741815371173\n0.123965034866\n";
	$points_string .= "-0.801656912336\n0.597784405036\n-0.000000000039\n-0.898999020841\n0.420027474451\n0.124006778982\n-0.249525433225\n0.768156660082\n-0.589637518944\n0.123965034866\n";
	$points_string .= "0.659046754960\n-0.741815371173\n-0.000000000039\n0.801656912336\n-0.597784405036\n0.124006778982\n0.898999020841\n-0.420027474451\n0.403740632272\n-0.364416027907\n";
	$points_string .= "-0.839162952267\n0.123965036782\n-0.659046753775\n-0.741815371905\n0.349661863232\n-0.585553996265\n-0.731343352236\n0.478971546438\n-0.679618729744\n-0.555611950811\n";
	$points_string .= "-0.839162952425\n-0.403740632174\n-0.364416027652\n-0.942395012918\n-0.000000001184\n-0.334502077164\n-0.947446268244\n-0.216102916071\n-0.235892133096\n-0.975639425135\n";
	$points_string .= "-0.219380291097\n-0.000000000079\n-0.403740632272\n0.364416027907\n-0.839162952267\n-0.123965036782\n0.659046753775\n-0.741815371905\n-0.349661863232\n0.585553996265\n";
	$points_string .= "-0.731343352236\n-0.478971546438\n0.679618729744\n-0.555611950811\n0.839162952425\n0.403740632174\n-0.364416027652\n0.942395012918\n0.000000001184\n-0.334502077164\n";
	$points_string .= "0.947446268244\n0.216102916071\n-0.235892133096\n0.975639425135\n0.219380291097\n-0.000000000079\n0.249525433225\n-0.768156660082\n-0.589637518944\n-0.123965034866\n";
	$points_string .= "-0.659046754960\n-0.741815371173\n0.000000000039\n-0.801656912336\n-0.597784405036\n-0.124006778982\n-0.898999020841\n-0.420027474451\n0.589637519041\n0.249525433481\n";
	$points_string .= "-0.768156659924\n0.865780406770\n0.200579641746\n-0.458467113946\n0.731343352197\n0.349661863208\n-0.585553996328\n0.679618729793\n0.555611950684\n-0.478971546516\n";
	$points_string .= "-0.993378151216\n0.000000000158\n-0.114890594426\n-0.865780407955\n0.200579641014\n-0.458467112030\n-0.947446268268\n0.216102916008\n-0.235892133057\n-0.898999020889\n";
	$points_string .= "0.420027474324\n-0.124006779060\n0.000000000158\n-0.114890594426\n-0.993378151216\n0.200579641014\n-0.458467112030\n-0.865780407955\n0.216102916008\n-0.235892133057\n";
	$points_string .= "-0.947446268268\n0.420027474324\n-0.124006779060\n-0.898999020889\n0.364416027907\n-0.839162952267\n-0.403740632272\n0.659046753775\n-0.741815371905\n-0.123965036782\n";
	$points_string .= "0.585553996265\n-0.731343352236\n-0.349661863232\n0.679618729744\n-0.555611950811\n-0.478971546438\n-0.364416027907\n-0.839162952267\n0.403740632272\n-0.659046753775\n";
	$points_string .= "-0.741815371905\n0.123965036782\n-0.585553996265\n-0.731343352236\n0.349661863232\n-0.679618729744\n-0.555611950811\n0.478971546438\n-0.114890594426\n-0.993378151216\n";
	$points_string .= "0.000000000158\n-0.458467112030\n-0.865780407955\n0.200579641014\n-0.235892133057\n-0.947446268268\n0.216102916008\n-0.124006779060\n-0.898999020889\n0.420027474324\n";
	$points_string .= "0.114890594426\n0.993378151216\n0.000000000158\n0.458467112030\n0.865780407955\n0.200579641014\n0.235892133057\n0.947446268268\n0.216102916008\n0.124006779060\n";
	$points_string .= "0.898999020889\n0.420027474324\n-0.403740632174\n-0.364416027652\n-0.839162952425\n-0.000000001184\n-0.334502077164\n-0.942395012918\n-0.216102916071\n-0.235892133096\n";
	$points_string .= "-0.947446268244\n-0.219380291097\n-0.000000000079\n-0.975639425135\n0.768156660082\n-0.589637518944\n-0.249525433225\n0.659046754960\n-0.741815371173\n0.123965034866\n";
	$points_string .= "0.801656912336\n-0.597784405036\n-0.000000000039\n0.898999020841\n-0.420027474451\n0.124006778982\n";

	my @points_lines = split(/\n/, $points_string);
	my $points_lines_ref = \@points_lines;
	return $points_lines_ref;
}



sub html_code_for_graphics_header {

	my $title = shift;
	my $html_code = '';

	$html_code .= "<table width='100%' cellspacing='0' cellpadding='0' border='0'><tr><td align='left'>\n" .
		"<table height='160' cellspacing='0' cellpadding='0' border='0' align='left' bgcolor='#a4a5e9'>\n" .
		#"	style='background-image: url(\"bluegold.gif\");\n" .
		#"		background-attachment: fixed;\n" .
		#"		background-position: top right;\n" .
		#"		background-repeat: no-repeat;'>\n" .
		"	<tr>\n" .
		"		<td height='160' valign='center'>\n" .
		"			<table height='160' cellspacing='0' cellpadding='0' border='0'>\n" .
		"				<tr>\n" .
		"					<td width='20'></td>\n" .
		"					<td width='740' valign='center'>\n" .
		"						<font size='5' color='black'>$title</font>\n" .
		"					</td>\n" .
		#"					<td width='560' valign='center'>\n" .
		#"						<font size='5' color='black'>$title</font>\n" .
		#"					</td>\n" .
		#"					<td width='180' valign='center'>\n" .
		#"						<img src='bluequoll.gif' width='180' height='100' border='1'>\n" .
		#"					</td>\n" .
		"					<td width='240'>\n" .
		"						<img src='bluegold.gif' width='240' height='160' border='0'>\n" .
		"					</td>\n" .
		"				</tr>\n" .
		"			</table>\n" .
		"		</td>\n" .
		"	</tr>\n" .
		"</table>\n" .
		"</td></tr><tr><td align='left'>\n" .
		"<table height='20' cellspacing='0' cellpadding='0' border='0'><tr><td></td></tr></table>\n<br>";

	return $html_code;
}



sub html_header {

	print $q->header();
	print	#$q->start_html(-title=>$text_title),
		"<html><head><title>VALPRED 2D</title>\n",
		'<style type="text/css">' . "\n",
		'<!--' . "\n",
		'body {' . "\n",
		'	background-image: url("bluequollbg.gif");' . "\n",
		'	background-position: left bottom;' . "\n",
		'	background-repeat: no-repeat;' . "\n",
		'}' . "\n",
		'.basestyle {' . "\n",
		'	color: #FFFFFF;' . "\n",
		'	font-size: 12px;' . "\n",
		'	font-family: Arial, Helvetica, Verdana, sans-serif;' . "\n",
		'}' . "\n",
		'h1 {' . "\n",
		'	font-size: 16px;' . "\n",
		'	font-weight: bold;' . "\n",
		'	font-family: Arial, Helvetica,Verdana, sans-serif;' . "\n",
		'}' . "\n",
		'body, td, p, input {' . "\n",
		'	font-size: 12px;' . "\n",
		'	font-family: Arial, Helvetica, Verdana, sans-serif;' . "\n",
		'	color: #333333;' . "\n",
		'}' . "\n",
		'.navlinks {' . "\n",
		'	font-size: 14px;' . "\n",
		'	font-family: Arial, Helvetica, Verdana, sans-serif;' . "\n",
		'	font-weight: bold;' . "\n",
		'	color: #FFFFFF;' . "\n",
		'}' . "\n",
		'-->' . "\n",
		'</style>' . "\n",
		"</head><body>\n";
		my $title = 'VALPRED 2D';
		my $output_line = html_code_for_graphics_header($title);
		print $output_line;
}



sub init {

	# initialise the displayed statistics totalled for all sequences in the prediction method

	$stats->{'qok'} = 0;
	$stats->{'qhtm_obs'} = 0;
	$stats->{'qhtm_prd'} = 0;
	$stats->{'q2'} = 0;
	$stats->{'q2t_obs'} = 0;
	$stats->{'q2t_prd'} = 0;
	$stats->{'q2n_obs'} = 0;
	$stats->{'q2n_prd'} = 0;

	# consider a tms to be correctly predicted if 3 residues overlap between the predicted and observed/true segments
	$stats->{'minimum_overlap'} = 3;

	# initialise the overall counts that will be used to calculate the displayed statistics

	$stats->{'stats_num_seqs'} = 0;
	$stats->{'stats_num_obs_tms_in_this_method'} = 0;
	$stats->{'stats_num_prd_tms_in_this_method'} = 0;
	$stats->{'stats_num_correct_obs_tms_in_this_method'} = 0;
	$stats->{'stats_num_correct_prd_tms_in_this_method'} = 0;
	$stats->{'stats_num_correct_obs_and_prd_tms_in_this_method'} = 0;
	$stats->{'stats_num_obs_tms_residues_in_this_method'} = 0;
	$stats->{'stats_num_prd_tms_residues_in_this_method'} = 0;
	$stats->{'stats_num_obs_nontms_residues_in_this_method'} = 0;
	$stats->{'stats_num_prd_nontms_residues_in_this_method'} = 0;
	$stats->{'stats_num_correct_obs_tms_residues_in_this_method'} = 0;
	$stats->{'stats_num_correct_prd_tms_residues_in_this_method'} = 0;
	$stats->{'stats_num_correct_obs_nontms_residues_in_this_method'} = 0;
	$stats->{'stats_num_correct_prd_nontms_residues_in_this_method'} = 0;
	$stats->{'stats_ratio_correct_residues_in_this_method'} = 0;
}



sub Mat44_multiply_Mat44 { #C++ Mat44<Type> Mat44<Type>::operator * (const Mat44& A) const  // MULTIPLICATION (*)
	my $M = shift;
	my $A = shift;

	# Mat44<Type> NewM( M[0]*A.M[0] + M[4]*A.M[1] + M[8]*A.M[2] + M[12]*A.M[3],      // ROW 1
			# M[0]*A.M[4] + M[4]*A.M[5] + M[8]*A.M[6] + M[12]*A.M[7],
			# M[0]*A.M[8] + M[4]*A.M[9] + M[8]*A.M[10] + M[12]*A.M[11],
			# M[0]*A.M[12] + M[4]*A.M[13] + M[8]*A.M[14] + M[12]*A.M[15],
	my $NewM = Mat44_new();
	$NewM->[0] = $M->[0]*$A->[0] + $M->[4]*$A->[1] + $M->[8]*$A->[2] + $M->[12]*$A->[3];
	$NewM->[4] = $M->[0]*$A->[4] + $M->[4]*$A->[5] + $M->[8]*$A->[6] + $M->[12]*$A->[7];
	$NewM->[8] = $M->[0]*$A->[8] + $M->[4]*$A->[9] + $M->[8]*$A->[10] + $M->[12]*$A->[11];
	$NewM->[12] = $M->[0]*$A->[12] + $M->[4]*$A->[13] + $M->[8]*$A->[14] + $M->[12]*$A->[15];

			# M[1]*A.M[0] + M[5]*A.M[1] + M[9]*A.M[2] + M[13]*A.M[3],      // ROW 2
			# M[1]*A.M[4] + M[5]*A.M[5] + M[9]*A.M[6] + M[13]*A.M[7],
			# M[1]*A.M[8] + M[5]*A.M[9] + M[9]*A.M[10] + M[13]*A.M[11],
			# M[1]*A.M[12] + M[5]*A.M[13] + M[9]*A.M[14] + M[13]*A.M[15],
	$NewM->[1] = $M->[1]*$A->[0] + $M->[5]*$A->[1] + $M->[9]*$A->[2] + $M->[13]*$A->[3];
	$NewM->[5] = $M->[1]*$A->[4] + $M->[5]*$A->[5] + $M->[9]*$A->[6] + $M->[13]*$A->[7];
	$NewM->[9] = $M->[1]*$A->[8] + $M->[5]*$A->[9] + $M->[9]*$A->[10] + $M->[13]*$A->[11];
	$NewM->[13] = $M->[1]*$A->[12] + $M->[5]*$A->[13] + $M->[9]*$A->[14] + $M->[13]*$A->[15];

			# M[2]*A.M[0] + M[6]*A.M[1] + M[10]*A.M[2] + M[14]*A.M[3],     // ROW 3
			# M[2]*A.M[4] + M[6]*A.M[5] + M[10]*A.M[6] + M[14]*A.M[7],
			# M[2]*A.M[8] + M[6]*A.M[9] + M[10]*A.M[10] + M[14]*A.M[11],
			#  M[2]*A.M[12] + M[6]*A.M[13] + M[10]*A.M[14] + M[14]*A.M[15],
	$NewM->[2] = $M->[2]*$A->[0] + $M->[6]*$A->[1] + $M->[10]*$A->[2] + $M->[14]*$A->[3];
	$NewM->[6] = $M->[2]*$A->[4] + $M->[6]*$A->[5] + $M->[10]*$A->[6] + $M->[14]*$A->[7];
	$NewM->[10] = $M->[2]*$A->[8] + $M->[6]*$A->[9] + $M->[10]*$A->[10] + $M->[14]*$A->[11];
	$NewM->[14] = $M->[2]*$A->[12] + $M->[6]*$A->[13] + $M->[10]*$A->[14] + $M->[14]*$A->[15];

			# M[3]*A.M[0] + M[7]*A.M[1] + M[11]*A.M[2] + M[15]*A.M[3],     // ROW 4
			# M[3]*A.M[4] + M[7]*A.M[5] + M[11]*A.M[6] + M[15]*A.M[7],
			# M[3]*A.M[8] + M[7]*A.M[9] + M[11]*A.M[10] + M[15]*A.M[11],
			# M[3]*A.M[12] + M[7]*A.M[13] + M[11]*A.M[14] + M[15]*A.M[15] );
	$NewM->[3] = $M->[3]*$A->[0] + $M->[7]*$A->[1] + $M->[11]*$A->[2] + $M->[15]*$A->[3];
	$NewM->[7] = $M->[3]*$A->[4] + $M->[7]*$A->[5] + $M->[11]*$A->[6] + $M->[15]*$A->[7];
	$NewM->[11] = $M->[3]*$A->[8] + $M->[7]*$A->[9] + $M->[11]*$A->[10] + $M->[15]*$A->[11];
	$NewM->[15] = $M->[3]*$A->[12] + $M->[7]*$A->[13] + $M->[11]*$A->[14] + $M->[15]*$A->[15];

	return $NewM; # return(NewM);
}



# // MAT-POINT MULTIPLICATION _WITHOUT_ PERSP DIV 
# // (for transforming points in space)
# // Assumes matrix is affine, i.e. bottom row is 0,0,0,1
sub Mat44_multVec3d { #C++ Vec3<Type> Mat44<Type>::multVec3d(const Vec3<Type>& P) const

	my $M = shift;
	my $P = shift;
	# Vec3<Type> NewP( (M[0]*P.x + M[4]*P.y + M[8]*P.z  + M[12]),
	#                  (M[1]*P.x + M[5]*P.y + M[9]*P.z  + M[13]),
	#                  (M[2]*P.x + M[6]*P.y + M[10]*P.z + M[14]) );
	my $NewP = Vec3_new_set3( 0, 0, 0 );
	$NewP->{'x'} = ($M->[0] * $P->{'x'}) + ($M->[4] * $P->{'y'}) + ($M->[8] * $P->{'z'}) + $M->[12];
	$NewP->{'y'} = ($M->[1] * $P->{'x'}) + ($M->[5] * $P->{'y'}) + ($M->[9] * $P->{'z'}) + $M->[13];
	$NewP->{'z'} = ($M->[2] * $P->{'x'}) + ($M->[6] * $P->{'y'}) + ($M->[10] * $P->{'z'}) + $M->[14];

	return $NewP; #C++ Vec3f or Vec3d					# return (NewP);
}



sub Mat44_new { # typedef Mat44<float> Mat44f;

	my $mat44;

	## void Mat44<Type>::identity()
	## M[0]=M[5]=M[10]=M[15]=1;
	## M[1]=M[2]=M[3]=M[4]=M[6]=M[7]=M[8]=M[9]=M[11]=M[12]=M[13]=M[14]=0;

	$mat44->[0] = 1;
	$mat44->[1] = 0;
	$mat44->[2] = 0;
	$mat44->[3] = 0;
	$mat44->[4] = 0;
	$mat44->[5] = 1;
	$mat44->[6] = 0;
	$mat44->[7] = 0;
	$mat44->[8] = 0;
	$mat44->[9] = 0;
	$mat44->[10] = 1;
	$mat44->[11] = 0;
	$mat44->[12] = 0;
	$mat44->[13] = 0;
	$mat44->[14] = 0;
	$mat44->[15] = 1;

	return $mat44; #C++ Mat44f or Mat44d
}



sub Mat44_rotate { #C++ void Mat44<Type>::rotate(Type RadAng, const Vec3<Type>& Axis)

	my $M = shift;
	my $RadAng = shift;
	my $Axis = shift;

	my $ca = cos( $RadAng ); # Type ca=(Type)cos(RadAng),
	my $sa = sin( $RadAng ); # sa=(Type)sin(RadAng);

	if (($Axis->{'x'} == 1) && ($Axis->{'y'} == 0) && ($Axis->{'z'} == 0)) { # if (Axis.x==1 && Axis.y==0 && Axis.z==0)  // ABOUT X-AXIS
		$M->[0] = 1;				# M[0]=1; M[4]=0;  M[8]=0;   M[12]=0;
		$M->[4] = 0;
		$M->[8] = 0;
		$M->[12] = 0;
		$M->[1] = 0;				# M[1]=0; M[5]=ca; M[9]=-sa; M[13]=0;
		$M->[5] = $ca;
		$M->[9] = -1 * $sa;
		$M->[13] = 0;
		$M->[2] = 0;				# M[2]=0; M[6]=sa; M[10]=ca; M[14]=0;
		$M->[6] = $sa;
		$M->[10] = $ca;
		$M->[14] = 0;
		$M->[3] = 0;				# M[3]=0; M[7]=0;  M[11]=0;  M[15]=1;
		$M->[7] = 0;
		$M->[11] = 0;
		$M->[15] = 1;

	} elsif (($Axis->{'x'} == 0) && ($Axis->{'y'} == 1) && ($Axis->{'z'} == 0)) { # else if (Axis.x==0 && Axis.y==1 && Axis.z==0)  // ABOUT Y-AXIS
		$M->[0] = $ca;				# M[0]=ca;  M[4]=0; M[8]=sa;  M[12]=0;
		$M->[4] = 0;
		$M->[8] = $sa;
		$M->[12] = 0;
		$M->[1] = 0;				# M[1]=0;   M[5]=1; M[9]=0;   M[13]=0;
		$M->[5] = 1;
		$M->[9] = 0;
		$M->[13] = 0;
		$M->[2] = -1 * $sa;			# M[2]=-sa; M[6]=0; M[10]=ca; M[14]=0;
		$M->[6] = 0;
		$M->[10] = $ca;
		$M->[14] = 0;
		$M->[3] = 0;				# M[3]=0;   M[7]=0; M[11]=0;  M[15]=1;
		$M->[7] = 0;
		$M->[11] = 0;
		$M->[15] = 1;

	} elsif (($Axis->{'x'} == 0) && ($Axis->{'y'} == 0) && ($Axis->{'z'} == 1)) { # else if (Axis.x==0 && Axis.y==0 && Axis.z==1)  // ABOUT Z-AXIS
		$M->[0] = $ca;				# M[0]=ca; M[4]=-sa; M[8]=0;  M[12]=0;
		$M->[4] = -1 * $sa;
		$M->[8] = 0;
		$M->[12] = 0;
		$M->[1] = $sa;				# M[1]=sa; M[5]=ca;  M[9]=0;  M[13]=0;
		$M->[5] = $ca;
		$M->[9] = 0;
		$M->[13] = 0;
		$M->[2] = 0;				# M[2]=0;  M[6]=0;   M[10]=1; M[14]=0;
		$M->[6] = 0;
		$M->[10] = 1;
		$M->[14] = 0;
		$M->[3] = 0;				# M[3]=0;  M[7]=0;   M[11]=0; M[15]=1;
		$M->[7] = 0;
		$M->[11] = 0;
		$M->[15] = 1;

	} else { # // ARBITRARY AXIS

		my $l = Vec3_lengthSqr( $Axis );	# Type l = Axis.lengthSqr();
		my $x = $Axis->{'x'};			# Type x, y, z;
		my $y = $Axis->{'y'};			# x=Axis.x, y=Axis.y, z=Axis.z;
		my $z = $Axis->{'z'};

		if ((($l > 1.0001) || ($l < 0.9999)) && ($l != 0)) { # if (l > Type(1.0001) || l < Type(0.9999) && l!=0)
			# // needs normalization
			$l = 1 / sqrt($l);		# l=Type(1.0)/sqrt(l);
			$x = $x * $l;			# x*=l; y*=l; z*=l;
			$y = $y * $l;
			$z = $z * $l;
		}
		my $x2 = $x * $x;			# Type x2=x*x, y2=y*y, z2=z*z;
		my $y2 = $y * $y;
		my $z2 = $z * $z;

		$M->[0] = $x2 + $ca * (1 - $x2);				# M[0]=x2+ca*(1-x2);
		$M->[4] = ($x * $y) + $ca * (-1 * $x * $y) + $sa * (-1 * $z);	# M[4]=(x*y)+ca*(-x*y)+sa*(-z); 
		$M->[8] = ($x * $z) + $ca * (-1 * $x * $z) + $sa * $y;		# M[8]=(x*z)+ca*(-x*z)+sa*y;
		$M->[1] = ($x * $y) + $ca * (-1 * $x * $y) + $sa * $z;		# M[1]=(x*y)+ca*(-x*y)+sa*z; 
		$M->[5] = $y2 + $ca * (1 - $y2);				# M[5]=y2+ca*(1-y2); 
		$M->[9] = ($y * $z) + $ca * (-1 * $y * $z) + $sa * (-1 * $x);	# M[9]=(y*z)+ca*(-y*z)+sa*(-x);
		$M->[2] = ($x * $z) + $ca * (-1 * $x * $z) + $sa * (-1 * $y);	# M[2]=(x*z)+ca*(-x*z)+sa*(-y); 
		$M->[6] = ($y * $z) + $ca * (-1 * $y * $z) + $sa * $x;		# M[6]=(y*z)+ca*(-y*z)+sa*x; 
		$M->[10] = $z2 + $ca * (1 - $z2);				# M[10]=z2+ca*(1-z2);
		$M->[12] = 0;							# M[12]=M[13]=M[14]=M[3]=M[7]=M[11]=0;
		$M->[13] = 0;
		$M->[14] = 0;
		$M->[3] = 0;
		$M->[7] = 0;
		$M->[11] = 0;
		$M->[15] = 1;							# M[15]=1;
	}

	return $M; #C++ Mat44f or Mat44d 
}



# //============================================================================
# // mat44.hpp : 4x4 OpenGL-style matrix template.
# //============================================================================
#
# // from GLUV
# //Author: 
# //Kenny Hoff <hoff at cs.unc.edu> , Bill Baxter <baxter at billbaxter.com> ,
# //Wes Hunt <hunt at cs.unc.edu> , Mark Harris <harrism at cs.unc.edu> 
#
# //----------------------------------------------------------------------------
# // M[16] = [ 0  4  8 12 ]   |  16 floats were used instead of the normal [4][4]
# //         [ 1  5  9 13 ]   |  to be compliant with OpenGL. OpenGL uses
# //         [ 2  6 10 14 ]   |  premultiplication with column vectors. These
# //         [ 3  7 11 15 ]   |  matrices can be fed directly into the OpenGL
# //                          |  matrix stack with glLoadMatrix or glMultMatrix.
# //
# // [ x' y' z' w' ] = [ 0  4  8 12 ]   [ x ]
# //                   [ 1  5  9 13 ] * [ y ] 
# //                   [ 2  6 10 14 ]   [ z ]
# //                   [ 3  7 11 15 ]   [ w ]
# // 
# // Loading a [4][4] format matrix directly into the matrix stack (assuming
# // premult/col vecs) results in a transpose matrix. M[0]=M[0][0], but
# // M[1]!=M[1][0] since M[0][1] would cast to M[1].
# //
# // However, if we assumed postmult/row vectors we could use [4][4] format,
# // but all transformations in this module would be transposed.
# //----------------------------------------------------------------------------
sub mat44impl_translate { # void Mat44<Type>::translate(const Vec3<Type>& T)

	my $M = shift;
	my $T = shift;

	# M[0]=1; M[4]=0;  M[8]=0;  M[12]=T.x;
	# M[1]=0; M[5]=1;  M[9]=0;  M[13]=T.y;
	# M[2]=0; M[6]=0;  M[10]=1; M[14]=T.z;
	# M[3]=0; M[7]=0;  M[11]=0; M[15]=1;
	$M->[0] = 1;
	$M->[1] = 0;
	$M->[2] = 0;
	$M->[3] = 0;
	$M->[4] = 0;
	$M->[5] = 1;
	$M->[6] = 0;
	$M->[7] = 0;
	$M->[8] = 0;
	$M->[9] = 0;
	$M->[10] = 1;
	$M->[11] = 0;
	$M->[12] = $T->{'x'};
	$M->[13] = $T->{'y'};
	$M->[14] = $T->{'z'};
	$M->[15] = 1;

	return $M;
}



sub output_benchmark_file {

	my $output_line = $input_fasta_id;
	print BENCHOUTFILE "$output_line\n";

	Structure_getTMDomains();
	my $tmDomains_dot_size = $#{$tmDomains} + 1;
	for ( my $i = 0; $i < $tmDomains_dot_size; $i++ ) {
		my $tms_start = $tmDomains->[$i]->{'min'} + 1;
		my $tms_end = $tmDomains->[$i]->{'max'} + 1;
		$output_line = $tms_start . "," . $tms_end;
		print BENCHOUTFILE "$output_line\n";
	}

	close BENCHOUTFILE;
	open( BENCHOUTFILE, ">>$benchmark_output_file") or
		die "Cannot reopen $benchmark_output_file for writing : $!\n";
}



sub output_evaluate_file {

	my $predicted_tms_array_start_ref;
	my $predicted_tms_array_end_ref;
	if ((defined($input_predictions_file)) || (defined($input_benchmarkpredictions_file))) {
		$predicted_tms_array_start_ref = shift;
		$predicted_tms_array_end_ref = shift;
	}

	my $output_line = $input_sequence_id . ':::::' . $input_fasta_id . ':::::' . $input_aa_sequence . ':::::';

	# get the number of residues length of this sequence
	# program loops and calculated totals and percentages depend on this

	my $stats_num_residues = length($input_aa_sequence);
	if ($input_aa_sequence eq '-') {
		if (($input_tms_line ne '') && ($input_tms_line ne '-')) {
			my @bits = split(/\,/, $input_tms_line);
			foreach my $this_tms (@bits) {
				my @bits2 = split(/\-/, $this_tms);
				my $start = $bits2[0];
				my $end = $bits2[1];
				if ($end > $stats_num_residues) {
					$stats_num_residues = $end;
				}
			}
		}
		my @temp_predicted_tms_array_end = @$predicted_tms_array_end_ref;		
		foreach my $this_end (@temp_predicted_tms_array_end) {
			if ($this_end > $stats_num_residues) {
				$stats_num_residues = $this_end;
			}
		}
	}

	my $output_tms_line = $input_tms_line;
	if ($output_tms_line eq '') {
		$output_tms_line = '-';
	}
	$output_line .= $output_tms_line . ':::::';

	# get the tms prediction that was predicted by valpred from an amino acid sequence

	my @predicted_tms_array;
	my $predicted_num_tms = 0;
	my @predicted_tms_array_start;
	my @predicted_tms_array_end;
	my $predicted_tms_array_index = -1;
	my $predicted_segments = '';
	for ( my $i = 1; $i <= $stats_num_residues; $i++ ) {
		$predicted_tms_array[$i] = 0;
	}
	if ((defined($input_predictions_file)) || (defined($input_benchmarkpredictions_file))) {
		@predicted_tms_array_start = @$predicted_tms_array_start_ref;
		@predicted_tms_array_end = @$predicted_tms_array_end_ref;
	} else {
		Structure_getTMDomains();
		my $tmDomains_dot_size = $#{$tmDomains} + 1;
		for ( my $i = 0; $i < $tmDomains_dot_size; $i++ ) {
			my $tms_start = $tmDomains->[$i]->{'min'} + 1;
			my $tms_end = $tmDomains->[$i]->{'max'} + 1;
			$predicted_tms_array_index++;
			$predicted_tms_array_start[$predicted_tms_array_index] = $tms_start;
			$predicted_tms_array_end[$predicted_tms_array_index] = $tms_end;
		}
	}
	for ( my $a = 0; $a < @predicted_tms_array_start; $a++ ) {
		my $tms_start = $predicted_tms_array_start[$a];
		my $tms_end = $predicted_tms_array_end[$a];
		for ( my $i = $tms_start; $i <= $tms_end; $i++ ) {
			$predicted_tms_array[$i] = 1;
		}
		$predicted_num_tms++;
		$predicted_segments .= $tms_start . '-' . $tms_end . ',';
	}
	if ($predicted_segments eq '') {
		$predicted_segments = '-';
	}
	$output_line .= $predicted_segments . ':::::';

	# The input tms line may have continguous or overlapping segments, 
	# start and end may be before position 1 or after the end position of the sequence,
	# so normalise the true tms segments before comparing to the predicted for accuracy.

	if (($input_tms_line ne '-') and ($input_tms_line ne '')) {
		my @bits = split(/\,/, $input_tms_line);
		foreach my $this_tms (@bits) {
			# could be 3-23, -1-18, -6--1
			$this_tms =~ /(-?\d+)-(-?\d+)/;
			my $start = $1;
			my $end = $2;
			# ignore residues before residue number 1
			if ($start <= 0) {
				$start = 1;
			}
			# ignore helices before residue number 1
			if ($end > 0) {
				if ($start > $stats_num_residues) {
					$start = $stats_num_residues;
				}
				if ($end > $stats_num_residues) {
					$end = $stats_num_residues;
				}
			}
		}
	}

	# now compare the observed tms list to the predicted tms list,
	# to calculate how many observed tms were predicted to exist.
	# an observed tms is predicted to exist if there is an overlap of 3 residues between them.
	# each tms is counted only once, so if a predicted tms straddles 2 observed tms, 
	# then only one is counted as being predicted because the prediction has predicted only 1 long tms, not 2 tms.
	# If the observed tms is parallel to the membrane, and it half-in half-out of the membrane,
	# ie, if eg. IIMIMMIMMIIMI or IIMIMMIMIIIIIIIII or OOMOMMOMMOOMO or OOMOMMOMOOOOOOOOO,
	# ie, if there are never more than 2 M next to each other,
	# then this tms is to be counted as being correctly predict regardless of if it was predicted or not.
	# If part of the observed tms is in the membrane and part of it is non-membrane, 
	# ie. iIIIIMIIMMIMm or oOOOOMMMMMMMMMMm,
	# the consider the parts that are non-membrane to be correctly predicted residues even if predicted as membrane.

	# 1jb0_A:::::>1jb0_A:::::
	# HELIX;;;;;24-31,44-56,64-98,99-105,122-128,142-151,154-184,192-230,238-244,245-252,261-270,293-316,324-333,341-350,351-378,386-420,428-435,436-469,484-501,532-560,590-623,637-644,645-663,667-692,693-712,723-755,:::::
	# HELIX;;;;;24-31,44-56,64-98,99-105,122-128,142-151,154-184,192-230,238-244,245-252,261-270,293-316,324-333,341-350,351-378,386-420,428-435,436-469,484-501,532-560,590-623,637-644,645-663,667-692,693-712,723-755,:::::
	# HELIX_INSIDE;;;;;24-31,44-48,50-51,54-54,64-71,180-184,192-192,324-324,326-327,329-331,333-333,341-350,351-352,414-420,428-435,436-438,558-560,590-591,695-695,698-699,702-703,705-712,723-724,:::::
	# HELIX_MEMBRANE;;;;;49-49,52-53,55-56,72-93,159-179,193-215,267-270,299-316,325-325,328-328,332-332,353-374,392-413,439-462,464-464,484-484,486-487,536-557,592-612,646-646,674-692,693-694,696-697,700-701,704-704,725-745,:::::
	# HELIX_OUTSIDE;;;;;94-98,99-105,122-128,142-151,154-158,216-230,238-244,245-252,261-266,293-298,375-378,386-391,463-463,465-469,485-485,488-501,532-535,613-623,637-644,645-645,647-663,667-673,746-755,:::::
	#	HELIX			= list of entire helices
	#	HELIX_INSIDE		= list of the parts of helices that fall in the inside of the cell membrane
	#	HELIX_MEMBRANE		= list of the parts of helices that fall in the membrane
	#	HELIX_OUTSIDE		= list of the parts of helices that fall in the outside of the cell membrane	# >1jb0_A
	# Photosystem I
	# SEQUENCE: GGGGGGGGGGGGRVVVDNDPVPTSFEKWAKPGHFDRTLARGPQTTTWIWNLHALAHDFDTHTSDLEDISRKIFSAHFGHLAVVFIWLSGMYFHGAKFSNYEAWLADPTGIKPSAQVVWPIVGQGILNGDVGGGFHGIQITSGLFQLWRASGITNEFQLYCTAIGGLVMAGLMLFAGWFHYHKRAPKLEWFQNVESMLNHHLAGLLGLGSLAWAGHQIHVSLPINKLLDAGVAAKDIPLPHEFILNPSLMAELYPKVDWGFFSGGGPFFTFNWAAYSDFLTFNGGLNPVTGGLWLSDTAHHHLAIAVLFIIAGHMYRTNWGIGHSLKEILEAHKGPFTGAGHKGLYEVLTTSWHAQLAINLAMMGSLSIIVAQHMYAMPPYPYLATDYPTQLSLFTHHMWIGGFLVVGGAAHGAIFMVRDYDPAMNQNNVLDRVLRHRDAIISHLNWVCIFLGFHSFGLYVHNDTMRAFGRPQDMFSDTGIQLQPVFAQWVQNLHTLAPGGTAPNAAATASVAFGGDVVAVGGKVAMMPIVLGTADFMVHHIHAFTIHVTVLILLKGVLFARSSRLIPDKANLGFRFPCDGPGRGGTCQVSGWDHVFLGLFWMYNCISVVIFHFSWKMQSDVWGTVAPDGTVSHITGGNFAQSAITINGWLRDFLWAQASQVIGSYGSALSAYGLLFLGAHFIWAFSLMFLFSGRGYWQELIESIVWAHNKLKVAPAIQPRALSIIQGRAVGVAHYLLGGIATTWAFFLARIISVG
	# OBSERVED: _______________________iIIIIIIi____________iIIIIMIIMMIMm_______iIIIIIIIMMMMMMMMMMMMMMMMMMMMMMOOOOooOOOOOo________________oOOOOOo_____________oOOOOOOOOo__oOOOOMMMMMMMMMMMMMMMMMMMMMIIIIi_______iMMMMMMMMMMMMMMMMMMMMMMMOOOOOOOOOOOOOOo_______oOOOOOooOOOOOOo________oOOOOOMMMm______________________oOOOOOMMMMMMMMMMMMMMMMMm_______iMIIMIIIMi_______iIIIIIIIIiiIMMMMMMMMMMMMMMMMMMMMMMOOOo_______oOOOOOMMMMMMMMMMMMMMMMMMMMMMIIIIIIi_______iIIIIIIiiIIMMMMMMMMMMMMMMMMMMMMMMMMOMOOOOo______________mOMMOOOOOOOOOOOOOo______________________________oOOOMMMMMMMMMMMMMMMMMMMMMMIIi_____________________________iIMMMMMMMMMMMMMMMMMMMMMOOOOOOOOOOo_____________oOOOOOOooMOOOOOOOOOOOOOOOOo___oOOOOOOMMMMMMMMMMMMMMMMMMmmMIMMIIMMIIMIIIIIIIi__________iIMMMMMMMMMMMMMMMMMMMMMOOOOOOOOOo
	# VALPRED2: ____________________________________________mMMMMMMMMMMMMMMMMMMMm_____mMMMMMMMMMMMMMMMMMMMMMMMm_________________mMMMMMMMMMm_________mMMMMMMMMMMMMMMMMMm___mMMMMMMMMMMMMMMMMMMMMMMMMMm________________mMMMMMMMMMMMMMMMMMMMMm________________mMMMMMMMMMMMMMMMMMMMMMMMMMm___mMMMMMMMMMMMMMMMMm___mMMMMMMMMMMMMMMMMMMMMMMMMMMMMMm___________________________mMMMMMMMMMMMMMMMMMMMMMMMMMMMm__________mMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMm_____________mMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMm__________________________________________mMMMMMMMMm________mMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMm_________________________mMMMMMMMMMMMMMMMMMMMMMMMMMMMMm______________________mMMMMMMMMMMMMMMMMMMMMMMMm___mMMMMMMMMMMMMMMMMMMMMMMMMm______________________________mMMMMMMMMMMMMMMMMMMMMMMMMMMMMm_
	# VALPRED1: __________________________________________________________________mMMMMMMMMMMMMMMMMMMMMMMMMm_____________________________________________________________mMMMMMMMMMMMMMMMMMMMMMm___________________________________________________________________________mMMMMMMMMMMMMMMMMMMMMMMMm________mMMMMMMMMMMMMMMMMMMMMMMMm_______________________________mMMMMMMMMMMMMMMMMMMMMMMMMMMMMMm__________mMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMm___________________mMMMMMMMMMMMMMMMMMMMm_________________________________________________________mMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMm_____________________________mMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMm____________________________mMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMm____________________________________mMMMMMMMMMMMMMMMMMMMMMMMMMm__
	# TMHMM2  : _______________________________________________________________________mMMMMMMMMMMMMMMMMMMMMMm_____________________________________________________________mMMMMMMMMMMMMMMMMMMMMMm_______________________________________________________________________________mMMMMMMMMMMMMMMMMm______________mMMMMMMMMMMMMMMMMMMMMMm____________________________________________mMMMMMMMMMMMMMMMMMMm___________________mMMMMMMMMMMMMMMMMMMMMMm____________________mMMMMMMMMMMMMMMMMMMMMMm___________________________________________________________________________mMMMMMMMMMMMMMMMMMMMMMm____________________________________mMMMMMMMMMMMMMMMMMMm_____________________________________________________mMMMMMMMMMMMMMMMMMMMMMm______________________________________mMMMMMMMMMMMMMMMMMMMMMm___

	my $this_reference_segments = $input_tms_line;
	my $this_method_segments = $predicted_segments;

	my $stats_qok_for_this_seq = 0;
	my $stats_num_obs_tms_in_this_seq = 0;
	my $stats_num_prd_tms_in_this_seq = 0;
	my $stats_num_correct_obs_tms_in_this_seq = 0;
	my $stats_num_correct_prd_tms_in_this_seq = 0;

	# list the observed helices for this sequence, observed by the benchmark reference
	my @observed_helices;
	my @observed_helices_start;
	my @observed_helices_end;
	my @observed_helices_already_counted;
	my @observed_residues;
	for ( my $i = 1; $i <= $stats_num_residues; $i++ ) {
		$observed_residues[$i] = 0;
	}
	if ($this_reference_segments ne '-') {
		@observed_helices = split(/\,/, $this_reference_segments);
		$stats_num_obs_tms_in_this_seq = @observed_helices;
		$stats->{'stats_num_obs_tms_in_this_method'} += $stats_num_obs_tms_in_this_seq;
		for ( my $i = 0; $i < @observed_helices; $i++ ) {
			my $observed_helix = $observed_helices[$i];
			$observed_helix =~ /(-?\d+)-(-?\d+)/;
			my $start = $1;
			my $end = $2;
			$observed_helices_start[$i] = $start;
			$observed_helices_end[$i] = $end;
			$observed_helices_already_counted[$i] = 0;
			for ( my $j = $start; $j <= $end; $j++ ) {
				$observed_residues[$j] = 1;
			}
		}
	}

	# list the predicted helices for this sequence, predicted by the method to be benchmarked
	my @predicted_helices;
	my @predicted_helices_start;
	my @predicted_helices_end;
	my @predicted_helices_already_counted;
	my @predicted_residues;
	for ( my $i = 1; $i <= $stats_num_residues; $i++ ) {
		$predicted_residues[$i] = 0;
	}
	if ($this_method_segments ne '-') {
		@predicted_helices = split(/\,/, $this_method_segments);
		$stats_num_prd_tms_in_this_seq = @predicted_helices;
		$stats->{'stats_num_prd_tms_in_this_method'} += $stats_num_prd_tms_in_this_seq;
		for ( my $i = 0; $i < @predicted_helices; $i++ ) {
			my $predicted_helix = $predicted_helices[$i];
			$predicted_helix =~ /(-?\d+)-(-?\d+)/;
			my $start = $1;
			my $end = $2;
			$predicted_helices_start[$i] = $start;
			$predicted_helices_end[$i] = $end;
			$predicted_helices_already_counted[$i] = 0;
			for ( my $j = $start; $j <= $end; $j++ ) {
				$predicted_residues[$j] = 1;
			}
		}
	}

	# for each observed membrane segment in this sequence, find out if it was predicted
	for ( my $i = 0; $i < @observed_helices; $i++ ) {
		my $observed_start = $observed_helices_start[$i];
		my $observed_end = $observed_helices_end[$i];
		my $found_overlapping_tms = 0;
		my $found_a_predicted_tms_for_this_observed_tms = 0;
		my $minimum_overlap_for_this_segment = $stats->{'minimum_overlap'};
		my $length_of_this_segment = $observed_end - $observed_start + 1;
		if ($length_of_this_segment < $minimum_overlap_for_this_segment) {
			$minimum_overlap_for_this_segment = $length_of_this_segment;
		}
		for ( my $j = 0; $j < @predicted_helices; $j++ ) {
			my $predicted_start = $predicted_helices_start[$j];
			my $predicted_end = $predicted_helices_end[$j];
			if (($found_overlapping_tms == 0) && ($predicted_helices_already_counted[$j] == 0)) {
				if ( (($observed_start <= $predicted_start) && ($predicted_start <= $observed_end)) ||
					(($observed_start <= $predicted_end) && ($predicted_end <= $observed_end)) ||
					(($predicted_start <= $observed_start) && ($observed_end <= $predicted_end)) ) {
					my $num_overlap = 0;
					for ( my $k = $observed_start; $k <= $observed_end; $k++ ) {
						if (($k >= 1) && ($k <= $stats_num_residues)) {
							if ($predicted_residues[$k] == 1) {
								$num_overlap++;
							}
						}
					}
					if ($num_overlap >= $minimum_overlap_for_this_segment) {
						$stats_num_correct_obs_tms_in_this_seq++;
						$stats->{'stats_num_correct_obs_tms_in_this_method'} = $stats->{'stats_num_correct_obs_tms_in_this_method'} + 1;
						$found_a_predicted_tms_for_this_observed_tms = 1;
						$predicted_helices_already_counted[$j] = 1;
						$found_overlapping_tms = 1;
					}
				}
			}
		}
	}

	# for each predicted membrane segment in this sequence, find out if it was observed
	for ( my $i = 0; $i < @predicted_helices; $i++ ) {
		my $predicted_start = $predicted_helices_start[$i];
		my $predicted_end = $predicted_helices_end[$i];
		my $found_overlapping_tms = 0;
		my $found_an_observed_tms_for_this_predicted_tms = 0;
		my $minimum_overlap_for_this_segment = $stats->{'minimum_overlap'};
		my $length_of_this_segment = $predicted_end - $predicted_start + 1;
		if ($length_of_this_segment < $minimum_overlap_for_this_segment) {
			$minimum_overlap_for_this_segment = $length_of_this_segment;
		}
		for ( my $j = 0; $j < @observed_helices; $j++ ) {
			my $observed_start = $observed_helices_start[$j];
			my $observed_end = $observed_helices_end[$j];
			if (($found_overlapping_tms == 0) && ($observed_helices_already_counted[$j] == 0)) {
				if ( (($predicted_start <= $observed_start) && ($observed_start <= $predicted_end)) ||
					(($predicted_start <= $observed_end) && ($observed_end <= $predicted_end)) ||
					(($observed_start <= $predicted_start) && ($predicted_end <= $observed_end)) ) {
					my $num_overlap = 0;
					for ( my $k = $predicted_start; $k <= $predicted_end; $k++ ) {
						if (($k >= 1) && ($k <= $stats_num_residues)) {
							if ($observed_residues[$k] == 1) {
								$num_overlap++;
							}
						}
					}
					if ($num_overlap >= $minimum_overlap_for_this_segment) {
						$stats_num_correct_prd_tms_in_this_seq++;
						$stats->{'stats_num_correct_prd_tms_in_this_method'} = $stats->{'stats_num_correct_prd_tms_in_this_method'} + 1;
						$found_an_observed_tms_for_this_predicted_tms = 1;
						$observed_helices_already_counted[$j] = 1;
						$found_overlapping_tms = 1;
					}
				}
			}
		}
	}

	# find out if the segments prediction for this sequence matches the observed segments
	if (($stats_num_obs_tms_in_this_seq == $stats_num_prd_tms_in_this_seq) && ($stats_num_obs_tms_in_this_seq == $stats_num_correct_prd_tms_in_this_seq)) {
		$stats_qok_for_this_seq = 1;
		$stats->{'stats_num_correct_obs_and_prd_tms_in_this_method'} += 1;
	}

	# count how many residues are correctly predicted as being tms or not tms

	my $stats_num_obs_tms_residues_in_this_seq = 0;
	my $stats_num_prd_tms_residues_in_this_seq = 0;
	my $stats_num_correct_obs_tms_residues_in_this_seq = 0;
	my $stats_num_correct_prd_tms_residues_in_this_seq = 0;
	my $stats_num_obs_nontms_residues_in_this_seq = 0;
	my $stats_num_prd_nontms_residues_in_this_seq = 0;
	my $stats_num_correct_obs_nontms_residues_in_this_seq = 0;
	my $stats_num_correct_prd_nontms_residues_in_this_seq = 0;
	my $stats_num_correct_residues_in_this_seq = 0;

	for ( my $i = 1; $i <= $stats_num_residues; $i++ ) {
		if ($observed_residues[$i] == 1) {
			$stats_num_obs_tms_residues_in_this_seq++;
			if ($predicted_residues[$i] == 1) {
				$stats_num_correct_obs_tms_residues_in_this_seq++;
			}
		} else { # $observed_residues[$i] == 0
			$stats_num_obs_nontms_residues_in_this_seq++;
			if ($predicted_residues[$i] == 0) {
				$stats_num_correct_obs_nontms_residues_in_this_seq++;
			}
		}
		if ($predicted_residues[$i] == 1) {
			$stats_num_prd_tms_residues_in_this_seq++;
			if ($observed_residues[$i] == 1) {
				$stats_num_correct_prd_tms_residues_in_this_seq++;
			}
		} else { # $predicted_residues[$i] == 0
			$stats_num_prd_nontms_residues_in_this_seq++;
			if ($observed_residues[$i] == 0) {
				$stats_num_correct_prd_nontms_residues_in_this_seq++;
			}
		}
		if (($observed_residues[$i] == 1) && ($predicted_residues[$i] == 1)) {
			$stats_num_correct_residues_in_this_seq++;
		} elsif (($observed_residues[$i] == 0) && ($predicted_residues[$i] == 0)) {
			$stats_num_correct_residues_in_this_seq++;
		}
	}

	# add the counts of correct residues to the overall statistics

	if ($stats_num_residues > 0) {
		$stats->{'stats_ratio_correct_residues_in_this_method'} += ($stats_num_correct_residues_in_this_seq / $stats_num_residues);
	}
	if ($stats_num_obs_tms_residues_in_this_seq > 0) {
		$stats->{'stats_num_correct_obs_tms_residues_in_this_method'} += $stats_num_correct_obs_tms_residues_in_this_seq;
		$stats->{'stats_num_obs_tms_residues_in_this_method'} += $stats_num_obs_tms_residues_in_this_seq;
	}
	if ($stats_num_prd_tms_residues_in_this_seq > 0) {
		$stats->{'stats_num_correct_prd_tms_residues_in_this_method'} += $stats_num_correct_prd_tms_residues_in_this_seq;
		$stats->{'stats_num_prd_tms_residues_in_this_method'} += $stats_num_prd_tms_residues_in_this_seq;
	}
	if ($stats_num_obs_nontms_residues_in_this_seq > 0) {
		$stats->{'stats_num_correct_obs_nontms_residues_in_this_method'} += $stats_num_correct_obs_nontms_residues_in_this_seq;
		$stats->{'stats_num_obs_nontms_residues_in_this_method'} += $stats_num_obs_nontms_residues_in_this_seq;
	}
	if ($stats_num_prd_nontms_residues_in_this_seq > 0) {
		$stats->{'stats_num_correct_prd_nontms_residues_in_this_method'} += $stats_num_correct_prd_nontms_residues_in_this_seq;
		$stats->{'stats_num_prd_nontms_residues_in_this_method'} += $stats_num_prd_nontms_residues_in_this_seq;
	}

	$stats->{'stats_num_seqs'} = $stats->{'stats_num_seqs'} + 1;

	# The measures for accuracy are those used by the Rost/Kernytsky benchmark server at http://cubic.bioc.columbia.edu/services/tmh_benchmark/ 
	# Andrew Kernytsky and Burkhard Rost, "Static benchmarking of membrane helix predictions", Nucleic Acids Research (2003), 31:13:, pp. 3642–3644
	# The measures are actually described in
	# Chien Peter Chen, Andrew Kernytsky, Burkhard Rost, "Transmembrane helix predictions revisited", Protein Science (2002), 11: pp. 2774–2791

	# measures outputted in the *.valpred_evaluate.txt file :
	# toplogy_correct ( = 1 or 0, for yes or no )
	# num_tms_observed, num_tms_observed_correctly_predicted, num_tms_predicted, num_tms_predicted_correctly_predicted
	# num_residues, num_observed_tms_residues, num_observed_non_tms_residues, num_correctly_predicted_residues, 
	# num_predicted_tms_residues, num_predicted_non_tms_residues, num_correctly_predicted_tms_residues, num_correctly_predicted_non_tms_residues

	# measures outputted in the *.valpred_evaluate_summary.txt file :
	# accuracy_topology = number of proteins having all observed tms are predicted and vice versa
	#				ie, having correct number of tms in roughly the correct places
	# accuracy_segments_observed = % of all true/observed tms helices that were predicted.
	#				minimum overlap of 3 residues.
	#				if 1 long predicted helix overlaps 2 observed helices, only count it once.
	# accuracy_segments_predicted = % of all predicted tms helices are really are true/observed helices.
	#				minimum overlap of 3 residues.
	#				if 1 long observed helix overlaps 2 different predicted helices, only count it once.
	# evalsum_accuracy_residues_predicted = Per residue accuracy
	#				(sum for each protein( no. of observed tms residues predicted correctly / no. residues in this protein) / no. proteins * 100%)
	# evalsum_accuracy_tms_residues_observed = Residues correctly predicted in membrane segments
	#				(sum for each protein( no. of residues correctly predicted in TM helices / no. residues observed in TM helices) * 100%)
	#				is a measure of overprediction
	# evalsum_accuracy_tms_residues_predicted = Accuracy of predicted residues
	#				(sum for each protein( no. of predicted tms residues predicted correctly / no. residues observed in TM helices) * 100%)
	#				is a measure of overprediction
	# evalsum_accuracy_non_tms_residues_observed = Residues correctly predicted to not be in membrane segments
	#				(sum for each protein( no. of residues correctly predicted to not be in TM helices / no. residues observed in to not be in TM helices) * 100%)
	#				is a measure of overprediction
	# evalsum_accuracy_non_tms_residues_predicted = Accuracy of predicted residues that were predicted to not be in TM helices
	#				(sum for each protein( no. of predicted non-tms residues predicted correctly / no. residues observed to not be in TM helices) * 100%)
	#				is a measure of overprediction

	# $output_evaluate_string = '';
	# $output_evaluate_string .= "toplogy_correct=$topology_is_correct,";
	# $output_evaluate_string .= "num_tms_observed=$observed_num_tms,";
	# $output_evaluate_string .= "num_tms_observed_correctly_predicted=$observed_num_tms_predicted_correctly,";
	# $output_evaluate_string .= "num_tms_predicted=$predicted_num_tms,";
	# $output_evaluate_string .= "num_tms_predicted_correctly_predicted=$predicted_num_tms_predicted_correctly,";
	# $output_evaluate_string .= "num_residues=$num_residues,";
	# $output_evaluate_string .= "num_observed_tms_residues=$num_observed_tms_residues,";
	# $output_evaluate_string .= "num_observed_non_tms_residues=$num_observed_non_tms_residues,";
	# $output_evaluate_string .= "num_predicted_tms_residues=$num_predicted_tms_residues,";
	# $output_evaluate_string .= "num_predicted_non_tms_residues=$num_predicted_non_tms_residues,";
	# $output_evaluate_string .= "num_correctly_predicted_residues=$num_correctly_predicted_residues,";
	# $output_evaluate_string .= "num_correctly_predicted_tms_residues=$num_correctly_predicted_tms_residues,";
	# $output_evaluate_string .= "num_correctly_predicted_non_tms_residues=$num_correctly_predicted_non_tms_residues,";
	# $output_line .= $output_evaluate_string . ':::::';

	my $output_evaluate_string = '';
	$output_evaluate_string .= "stats_qok_for_this_seq=$stats_qok_for_this_seq,";
	$output_evaluate_string .= "stats_num_obs_tms_in_this_seq=$stats_num_obs_tms_in_this_seq,";
	$output_evaluate_string .= "stats_num_prd_tms_in_this_seq=$stats_num_prd_tms_in_this_seq,";
	$output_evaluate_string .= "stats_num_correct_obs_tms_in_this_seq=$stats_num_correct_obs_tms_in_this_seq,";
	$output_evaluate_string .= "stats_num_correct_prd_tms_in_this_seq=$stats_num_correct_prd_tms_in_this_seq,";
	$output_evaluate_string .= "stats_num_obs_tms_residues_in_this_seq=$stats_num_obs_tms_residues_in_this_seq,";
	$output_evaluate_string .= "stats_num_prd_tms_residues_in_this_seq=$stats_num_prd_tms_residues_in_this_seq,";
	$output_evaluate_string .= "stats_num_correct_obs_tms_residues_in_this_seq=$stats_num_correct_obs_tms_residues_in_this_seq,";
	$output_evaluate_string .= "stats_num_correct_prd_tms_residues_in_this_seq=$stats_num_correct_prd_tms_residues_in_this_seq,";
	$output_evaluate_string .= "stats_num_obs_nontms_residues_in_this_seq=$stats_num_obs_nontms_residues_in_this_seq,";
	$output_evaluate_string .= "stats_num_prd_nontms_residues_in_this_seq=$stats_num_prd_nontms_residues_in_this_seq,";
	$output_evaluate_string .= "stats_num_correct_obs_nontms_residues_in_this_seq=$stats_num_correct_obs_nontms_residues_in_this_seq,";
	$output_evaluate_string .= "stats_num_correct_prd_nontms_residues_in_this_seq=$stats_num_correct_prd_nontms_residues_in_this_seq,";
	$output_evaluate_string .= "stats_num_correct_residues_in_this_seq=$stats_num_correct_residues_in_this_seq,";

	$output_line .= $output_evaluate_string . ':::::';

	print EVALOUTFILE "$output_line\n";

	close EVALOUTFILE;
	open( EVALOUTFILE, ">>$evaluate_output_file") or
		die "Cannot open $evaluate_output_file for writing : $!\n";
}



sub output_evaluate_summary_file {

	# do the final calculations that will be displayed for this method's statistics, and format the figures for display

	if ($stats->{'stats_num_seqs'} > 0 ) {
		$stats->{'qok'} = $stats->{'stats_num_correct_obs_and_prd_tms_in_this_method'} * 100 / $stats->{'stats_num_seqs'};
		$stats->{'qok'} = sprintf("%0.f", $stats->{'qok'});
	} else {
		if ($stats->{'stats_num_correct_obs_and_prd_tms_in_this_method'} == 0) {
			$stats->{'qok'} = 100;
		} else {
			$stats->{'qok'} = 0;
		}
	}

	if ($stats->{'stats_num_obs_tms_in_this_method'} > 0) {
		$stats->{'qhtm_obs'} = $stats->{'stats_num_correct_obs_tms_in_this_method'} * 100 / $stats->{'stats_num_obs_tms_in_this_method'};
		$stats->{'qhtm_obs'} = sprintf("%0.f", $stats->{'qhtm_obs'});
	} else {
		if ($stats->{'stats_num_correct_obs_tms_in_this_method'} == 0) {
			$stats->{'qhtm_obs'} = 100;
		} else {
			$stats->{'qhtm_obs'} = 0;
		}
	}

	if ($stats->{'stats_num_prd_tms_in_this_method'} > 0) {
		$stats->{'qhtm_prd'} = $stats->{'stats_num_correct_prd_tms_in_this_method'} * 100 / $stats->{'stats_num_prd_tms_in_this_method'};
		$stats->{'qhtm_prd'} = sprintf("%0.f", $stats->{'qhtm_prd'});
	} else {
		if ($stats->{'stats_num_correct_prd_tms_in_this_method'} == 0) {
			$stats->{'qhtm_prd'} = 100;
		} else {
			$stats->{'qhtm_prd'} = 0;
		}
	}

	if ($stats->{'stats_num_seqs'} > 0 ) {
		$stats->{'q2'} = $stats->{'stats_ratio_correct_residues_in_this_method'} * 100 / $stats->{'stats_num_seqs'};
		$stats->{'q2'} = sprintf("%0.f", $stats->{'q2'});
	} else {
		if ($stats->{'stats_ratio_correct_residues_in_this_method'} == 0) {
			$stats->{'q2'} = 100;
		} else {
			$stats->{'q2'} = 0;
		}
	}

	if ($stats->{'stats_num_obs_tms_residues_in_this_method'} > 0) {
		$stats->{'q2t_obs'} = ($stats->{'stats_num_correct_prd_tms_residues_in_this_method'} * 100) / $stats->{'stats_num_obs_tms_residues_in_this_method'};
	} else {
		if ($stats->{'stats_num_correct_prd_tms_residues_in_this_method'} == 0) {
			$stats->{'q2t_obs'} = 100;
		} else {
			$stats->{'q2t_obs'} = 0;
		}
	}
	$stats->{'q2t_obs'} = sprintf("%0.f", $stats->{'q2t_obs'});

	if ($stats->{'stats_num_prd_tms_residues_in_this_method'} > 0) {
		$stats->{'q2t_prd'} = ($stats->{'stats_num_correct_prd_tms_residues_in_this_method'} * 100) / $stats->{'stats_num_prd_tms_residues_in_this_method'};
	} else {
		if ($stats->{'stats_num_correct_prd_tms_residues_in_this_method'} == 0) {
			$stats->{'q2t_prd'} = 100;
		} else {
			$stats->{'q2t_prd'} = 0;
		}
	}
	$stats->{'q2t_prd'} = sprintf("%0.f", $stats->{'q2t_prd'});

	if ($stats->{'stats_num_obs_nontms_residues_in_this_method'} > 0) {
		$stats->{'q2n_obs'} = ($stats->{'stats_num_correct_prd_nontms_residues_in_this_method'} * 100) / $stats->{'stats_num_obs_nontms_residues_in_this_method'};
	} else {
		if ($stats->{'stats_num_correct_prd_nontms_residues_in_this_method'} == 0) {
			$stats->{'q2n_obs'} = 100;
		} else {
			$stats->{'q2n_obs'} = 0;
		}
	}
	$stats->{'q2n_obs'} = sprintf("%0.f", $stats->{'q2n_obs'});

	if ($stats->{'stats_num_prd_nontms_residues_in_this_method'} > 0) {
		$stats->{'q2n_prd'} = ($stats->{'stats_num_correct_prd_nontms_residues_in_this_method'} * 100) / $stats->{'stats_num_prd_nontms_residues_in_this_method'};
	} else {
		if ($stats->{'stats_num_correct_prd_nontms_residues_in_this_method'} == 0) {
			$stats->{'q2n_prd'} = 100;
		} else {
			$stats->{'q2n_prd'} = 0;
		}
	}
	$stats->{'q2n_prd'} = sprintf("%0.f", $stats->{'q2n_prd'});

	my $evaluate_summary_output_file = "$input_file_name_for_output_files.valpred_evaluate_summary.txt";
	open( EVALSUMOUTFILE, ">$evaluate_summary_output_file") or
		die "Cannot open $evaluate_summary_output_file for writing : $!\n";

	my $outline = "\n\n";;
	$outline .= "\tInput file : $input_file_name_for_output_files\n\n";
	$outline .= "\tNumber of amino acid residue sequences : " . $stats->{'stats_num_seqs'} . "\n\n";
	$outline .= "\tTopology accuracy : " . $stats->{'stats_num_correct_obs_and_prd_tms_in_this_method'} . "\n";
	$outline .= "\t\t\tnumber of proteins having all observed membrane helices are predicted and vice versa,\n";
	$outline .= "\t\t\tie, having correct number of membrane helices in roughly the correct places\n\n";
	$outline .= "\tTopology accuracy : " . $stats->{'qok'} . "%\n\n";
	$outline .= "\tAccuracy of observed membrane helices : " . $stats->{'qhtm_obs'} . "%\n";
	$outline .= "\t\t\t% of all true/observed membrane helices that were predicted.\n";
	$outline .= "\t\t\tminimum overlap of 3 residues.\n";
	$outline .= "\t\t\tif 1 long predicted membrane helix overlaps 2 observed membrane helices, it is only counted once.\n\n";
	$outline .= "\tAccuracy of predicted membrane helices : " . $stats->{'qhtm_prd'}. "%\n";
	$outline .= "\t\t\t% of all predicted membrane helices helices that really are true/observed membrane helices.\n";
	$outline .= "\t\t\tminimum overlap of 3 residues.\n";
	$outline .= "\t\t\tif 1 long observed membrane helix overlaps 2 different predicted membrane helices, it is only counted once.\n\n";
	$outline .= "\tPer residue accuracy : " . $stats->{'q2'} . "%\n";
	$outline .= "\t\t\tsum for each protein( no. of observed helical membrane residues predicted correctly / no. residues in this protein) / no. proteins * 100%\n\n";
	$outline .= "\tResidues correctly predicted in membrane segments : " . $stats->{'q2t_obs'} . "%\n";
	$outline .= "\t\t\tsum for each protein( no. of residues correctly predicted in membrane helices / no. residues observed in membrane helices) * 100%\n";
	$outline .= "\t\t\tis a measure of overprediction.\n";
	$outline .= "\t\t\tambiguous residues (non-membrane residues in a helix that also contains membrane residues) are counted as being helical membrane residues that were predicted correctly, regardless of the prediction.\n\n";
	$outline .= "\tAccuracy of predicted residues : " . $stats->{'q2t_prd'}. "%\n";
	$outline .= "\t\t\tsum for each protein( no. of predicted helical membrane residues predicted correctly / no. residues observed in membrane helices) * 100%\n";
	$outline .= "\t\t\tis a measure of overprediction.\n";
	$outline .= "\t\t\tambiguous residues (non-membrane residues in a helix that also contains membrane residues) are counted as being predictions that were predicted correctly, regardless of the prediction.\n\n";
	$outline .= "\tResidues correctly predicted to not be in membrane helices : " . $stats->{'q2n_obs'} . "%\n";
	$outline .= "\t\t\tsum for each protein( no. of residues correctly predicted to not be in membrane helices / no. residues observed to not be in membrane helices) * 100%\n";
	$outline .= "\t\t\tis a measure of overprediction.\n";
	$outline .= "\t\t\tambiguous residues (non-membrane residues in a helix that also contains membrane residues) are counted as being correctly predicted.\n\n";
	$outline .= "\tAccuracy of predicted residues that were predicted to not be in membrane helices: " . $stats->{'q2n_prd'} . "%\n";
	$outline .= "\t\t\tsum for each protein( no. of predicted non-helical-membrane residues predicted correctly / no. residues observed to not be in membrane helices) * 100%\n";
	$outline .= "\t\t\tis a measure of overprediction.\n";
	$outline .= "\t\t\tambiguous residues (non-membrane residues in a helix that also contains membrane residues) are counted as being correctly predicted.\n\n";

	print EVALSUMOUTFILE $outline;

	$outline = "\n\n\n\n";
	$outline .= "High Resolution Integral Membrane Proteins\n";
	$outline .= "\t\tPer-segment-accuracy\tPer-residue-accuracy................\n";
	$outline .= "Method\t\tQok\tQhtm\tQhtm\tQ2\tQ2T\tQ2T\tQ2N\tQ2N\n";
	$outline .= "\t\t\%\t\%obs\t\%prd\t\%\t\%obs\t\%prd\t\%obs\t\%prd\n\n";
	$outline .= "\t\t" . $stats->{'qok'}. "\t" . $stats->{'qhtm_obs'} . "\t" . $stats->{'qhtm_prd'} . "\t";
	$outline .= $stats->{'q2'} . "\t" . $stats->{'q2t_obs'} . "\t" . $stats->{'q2t_prd'} . "\t" . $stats->{'q2n_obs'} . "\t" . $stats->{'q2n_prd'} . "\n";
	$outline .= "\n\n\n";

	print EVALSUMOUTFILE $outline;

	$outline = "\n\n";
	$outline .= "Qok \% : Percentage of proteins for which all membrane helices are predicted correctly\n";
	$outline .= "Qhtm \%obs : Percentage of all observed membrane helices that are predicted correctly\n";
	$outline .= "Qhtm \%prd : Percentage of all predicted membrane helices that are predicted correctly\n";
	$outline .= "Q2 \% : Percentage of correctly predicted residues in two-states: membrane helix / nonmembrane helix\n";
	$outline .= "Q2T \%obs : Percentage of all observed membrane helix residues that are predicted correctly\n";
	$outline .= "Q2T \%prd : Percentage of all predicted membrane helix residues that are predicted correctly\n";
	$outline .= "Q2N \%obs : Percentage of all observed non-membrane helix residues that are predicted correctly\n";
	$outline .= "Q2N \%prd : Percentage of all predicted non-membrane helix residues that are predicted correctly\n";
	$outline .= "\n\n\n";

	print EVALSUMOUTFILE $outline;

	close EVALSUMOUTFILE;
}



sub output_chart_coords_link {

	my $chart_coords = shift;

	my $request_uri = $ENV{'REQUEST_URI'};
	my @bits = split(/\//, $request_uri);
	my $pgm_name = $bits[$#bits];

	my $output_line;
	$output_line .= "<form action='$pgm_name' method='post'>";
	$output_line .= "<table cellspacing='0' cellpadding='0' border='0'><tr><td valign='top'>";
	$output_line .= "<input type='hidden' id='chart_coords' name='chart_coords' value='$chart_coords'>";
	$output_line .= "<input type='submit' name='submit' value='Click here to see chart'>";
	$output_line .= "</td><td width='10'></td><td valign='top'>";
	$output_line .= "<== PLEASE NOTE : If you are viewing this from a Windows machine,<br>";
	$output_line .= "your browser may not allow you to view this dynamically generated chart conveniently in the browser.<br>";
	$output_line .= "Instead, it may ask you to download the graph as a file and save it to disk.<br>";
	$output_line .= "When you do this, you may need to rename the ending file extension of the downloaded file to gif (eg. rename the file to valpred.gif),<br>";
	$output_line .= "so that double-clicking this saved file on your disk will display this graph that has been dynamically generated for you.";
	$output_line .= "</td></tr></table>\n\n";
	$output_line .= "</form>\n\n";

	print $output_line;
}



sub PolymerBuilder_addResidue { # void PolymerBuilder::addResidue(Residue &r, float phi, float psi, float omega, Structure &s)

	my $r = shift; #C++ Residue
	my $phi = shift; #C++ float
	my $psi = shift; #C++ float
	my $omega = shift; #C++ float

	PolymerBuilder_assignResType( $r );				# assignResType(r);
	my $C = Vec3_new(); #C++ Vec3f					# Vec3f C,CA,N,O;
	my $CA = Vec3_new(); #C++ Vec3f
	my $N = Vec3_new(); #C++ Vec3f
	my $O = Vec3_new(); #C++ Vec3f
	my $s_dot_allResidues_dot_size = $#{$structure->{'allResidues'}} + 1;
	if ($s_dot_allResidues_dot_size == 0) {				# if(s.allResidues.size() == 0){
									# //CA.set(0.0f, 0.0f, 0.0f); 
									# //C.set(-0.551f, 0.0f, -1.422f); 
									# ////O.set(-1.535f, -0.689f, -1.687f);
									# //O.set(-1.337f,  -0.836f,  -1.799f );
									# //
									# //	//-1.535  -0.689  -1.687
									# //N.set(1.458f, 0.0f, 0.0f);	
		$CA = Vec3_set3( $CA, -0.560, 1.090, -3.601 );		# CA.set(-0.560f,   1.090f,  -3.601f); 
		$C = Vec3_set3( $C, -0.147, -0.145, -4.392 );		# C.set(-0.147f,  -0.145f,  -4.392f); 
									# //O.set(-1.535f, -0.689f, -1.687f);
		$O = Vec3_set3( $O, -0.920, -0.650, -5.206 );		# O.set(-0.920f,  -0.650f,  -5.206f );
									# //-1.535  -0.689  -1.687
		$N = Vec3_set3( $N, -0.115, 0.971, -2.217 );		# N.set(-0.115f,   0.971f,  -2.217f);

	} else {

		my $temp_r;
		my $temp_s;
		# Vec3f preO = s.allAtoms[s.allResidues[s.allResidues.size() - 1].O].coords; // O from previous residue
		$temp_r = $structure->{'allResidues'}->[$s_dot_allResidues_dot_size - 1]->{'O'};
		$temp_s = $structure->{'allAtoms'}->[$temp_r]->{'coords'};
		my $preO = Vec3_new_set1( $temp_s );
		# Vec3f preCA = s.allAtoms[s.allResidues[s.allResidues.size() - 1].CA].coords; // CA from previous residue
		$temp_r = $structure->{'allResidues'}->[$s_dot_allResidues_dot_size - 1]->{'CA'};
		$temp_s = $structure->{'allAtoms'}->[$temp_r]->{'coords'};
		my $preCA = Vec3_new_set1( $temp_s );
		# Vec3f preC = s.allAtoms[s.allResidues[s.allResidues.size() - 1].C].coords; // C from previous residue
		$temp_r = $structure->{'allResidues'}->[$s_dot_allResidues_dot_size - 1]->{'C'};
		$temp_s = $structure->{'allAtoms'}->[$temp_r]->{'coords'};
		my $preC = Vec3_new_set1( $temp_s );

		my $C_CA = Vec3_subtract_Vec3( $preCA, $preC );		# Vec3f C_CA = preCA - preC;
		my $C_O = Vec3_subtract_Vec3( $preO, $preC );		# Vec3f C_O = preO - preC;

		my $norAmide = Vec3_cross( $C_CA, $C_O );		# Vec3f norAmide = C_CA.cross(C_O); // normal to the amide plane

		my $m = Mat44_new();					# Mat44f m;

		# // rotate the C-CA from previous residue to form N
		# // the angle is 116.1 degree - psi
		my $deg_to_rad = PolymerBuilder_degToRad( -116.1 );
		$m = Mat44_rotate( $m, $deg_to_rad, $norAmide );	# m.rotate(degToRad(-116.1f), norAmide);
		$C_CA = Vec3_normalize( $C_CA );
		$N = Mat44_multVec3d( $m, $C_CA );			# N = m.multVec3d(C_CA.normalize());
		$N = Vec3_multiply_number( $N, 1.329 );			# N *= 1.329f;
		$N = Vec3_add_Vec3( $N, $preC );			# N += preC;

		# // rotate the amidebond to form CA
		# // angle is 121.88 - phi
		$deg_to_rad = PolymerBuilder_degToRad( 121.88 );
		$m = Mat44_rotate( $m, $deg_to_rad, $norAmide );	# m.rotate(degToRad(121.88f), norAmide);
		my $param1 = Vec3_subtract_Vec3( $preC, $N );
		$CA = Mat44_multVec3d( $m, $param1 );			# CA = m.multVec3d(preC-N);
		$CA = Vec3_normalize( $CA );				# CA.normalize();
		$CA = Vec3_multiply_number( $CA, 1.458 );		# CA *= 1.458f;

		# // rotate by omega - 180
		$deg_to_rad = PolymerBuilder_degToRad( $omega - 180 );
		$param1 = Vec3_subtract_Vec3( $N, $preC );
		$m = Mat44_rotate( $m, $deg_to_rad, $param1 );		# m.rotate(degToRad(omega-180.0f), N-preC);
		$CA = Mat44_multVec3d( $m, $CA );			# CA = m.multVec3d(CA);
		$CA = Vec3_add_Vec3( $CA, $N );				# CA += N;

		# // position Vector of C when phi=0
		$deg_to_rad = PolymerBuilder_degToRad( 51.98 );
		$param1 = Vec3_subtract_Vec3( $preC, $N );
		my $param2 = Vec3_subtract_Vec3( $CA, $N );
		$param1 = Vec3_cross( $param1, $param2 );
		$m = Mat44_rotate( $m, $deg_to_rad, $param1 );		# m.rotate(degToRad(51.98f), (preC-N).cross(CA-N));
		$param1 = Vec3_subtract_Vec3( $preC, $N );
		$C = Mat44_multVec3d( $m, $param1 );			# C = m.multVec3d(preC-N);
		$C = Vec3_normalize( $C );				# C.normalize();
		$C = Vec3_multiply_number( $C, 1.525 );			# C *= 1.525f;

		# // rotate by phi
		$deg_to_rad = PolymerBuilder_degToRad( -1 * $phi );
		$param1 = Vec3_subtract_Vec3( $N, $CA );
		$m = Mat44_rotate( $m, $deg_to_rad, $param1 );		# m.rotate(degToRad(-phi), N-CA);
		$C = Mat44_multVec3d( $m, $C );				# C = m.multVec3d(C);
		$C = Vec3_add_Vec3( $C, $CA );				# C += CA;

		# // rotate the N-CA to form O, psi = 180
		$deg_to_rad = PolymerBuilder_degToRad( 50.5 );
		$param1 = Vec3_subtract_Vec3( $N, $CA );
		$param2 = Vec3_subtract_Vec3( $C, $CA );
		$param1 = Vec3_cross( $param1, $param2 );
		$m = Mat44_rotate( $m, $deg_to_rad, $param1 );		# m.rotate(degToRad(50.5f), (N-CA).cross(C-CA));
		$param1 = Vec3_subtract_Vec3( $N, $CA );
		$O = Mat44_multVec3d( $m, $param1 );			# O = m.multVec3d(N-CA);

		# // rotate by psi
		$deg_to_rad = PolymerBuilder_degToRad( $psi - 180 );
		$param1 = Vec3_subtract_Vec3( $C, $CA );
		$m = Mat44_rotate( $m, $deg_to_rad, $param1 );		# m.rotate(degToRad(psi-180.0f), C-CA);
		$O = Mat44_multVec3d( $m, $O );				# O = m.multVec3d(O);
		$O = Vec3_normalize( $O );				# O.normalize();
		$O = Vec3_multiply_number( $O, 1.231 );			# O *= 1.231f;
		$O = Vec3_add_Vec3( $O, $C );				# O += C;
	}

	my $i;
	my @newAtoms_array;
	my $newAtoms = \@newAtoms_array;
	for ( $i = 0; $i < 4; $i++ ) {
		$newAtoms->[$i] = Atom_new();				# Atom newAtoms[4];
	}
	$newAtoms->[0]->{'coords'} = $N;				# newAtoms[0].coords = N;
	$newAtoms->[0]->{'name'} = "N";					# newAtoms[0].name = "N";
	$newAtoms->[1]->{'coords'} = $CA;				# newAtoms[1].coords = CA;
	$newAtoms->[1]->{'name'} = "CA";				# newAtoms[1].name = "CA";
	$newAtoms->[2]->{'coords'} = $C;				# newAtoms[2].coords = C;
	$newAtoms->[2]->{'name'} = "C";					# newAtoms[2].name = "C";
	$newAtoms->[3]->{'coords'} = $O;				# newAtoms[3].coords = O;
	$newAtoms->[3]->{'name'} = "O";					# newAtoms[3].name = "O";

	my $push_back_index;
	for ( $i = 0; $i < 4; $i++ ) {					# for (unsigned int i = 0; i < 4; i++){
		my $s_dot_allAtoms_dot_size = $#{$structure->{'allAtoms'}} + 1;
		$newAtoms->[$i]->{'atomID'} = $s_dot_allAtoms_dot_size;	# newAtoms[i].atomID = s.allAtoms.size();
		$newAtoms->[$i]->{'index'} = $s_dot_allAtoms_dot_size;	# newAtoms[i].index = s.allAtoms.size();
		$newAtoms->[$i]->{'chainID'} = 	$r->{'chainID'};	# newAtoms[i].chainID = r.chainID;
		my $s_dot_allResidues_dot_size = $#{$structure->{'allResidues'}} + 1;
		$newAtoms->[$i]->{'resSeq'} = $s_dot_allResidues_dot_size; # newAtoms[i].resSeq = (unsigned int)s.allResidues.size();
		Atom_setAtomType( $newAtoms->[$i], $r, $s_dot_allAtoms_dot_size ); # newAtoms[i].setAtomType(r,(unsigned int)s.allAtoms.size());
		$push_back_index = $#{$structure->{'allAtoms'}} + 1;
		$structure->{'allAtoms'}->[$push_back_index] = $newAtoms->[$i]; # s.allAtoms.push_back(newAtoms[i]);
	}

	PolymerBuilder_createSideChain( $r );				# createSideChain(r,s);
	$s_dot_allResidues_dot_size = $#{$structure->{'allResidues'}} + 1;
	$r->{'resSeq'} = $s_dot_allResidues_dot_size;			# r.resSeq = s.allResidues.size();
	$push_back_index = $#{$structure->{'allResidues'}} + 1;
	$structure->{'allResidues'}->[$push_back_index] = $r;		# s.allResidues.push_back(r);
}



sub PolymerBuilder_assignResType { # void PolymerBuilder::assignResType(Residue &r)

	# // Checks whether the amino acid lies within the HELIX or 
	# // SHEET ranges. Otherwise allocates their secondary structure
	# // to OTHER
	# ///////////////////////////////////////////////////////////////
	# //for (unsigned int i = 0;	i< s.HelixInfo.size();i++){
	#	if( r.resSeq >= helices.min &&
	#		r.resSeq <= helices.max &&
	#		helices.ChainID == r.chainID){
	#		r.SS = HELIX;
	#	}
	# //}
	# //for (i = 0;	i< s.SheetInfo.size(); i++){
	# //	if( r.resSeq >= s.SheetInfo[i].min &&
	# //		r.resSeq <= s.SheetInfo[i].max &&
	# //		s.SheetInfo[i].ChainID == r.chainID){
	# //		r.SS = SHEET;
	# //	}
	# //}

	# Checks whether the amino acid lies within the HELIX or SHEET ranges. 
	# Otherwise allocates their secondary structure to OTHER.

	my $r = shift; #C++ Residue

	$r->{'SS'} = OTHER;

	if (($r->{'resSeq'} >= $helices->{'min'}) && ($r->{'resSeq'} <= $helices->{'max'}) && ($helices->{'chainID'} eq $r->{'chainID'})) {
									# if( r.resSeq >= helices.min &&
									# r.resSeq <= helices.max &&
									# helices.ChainID == r.chainID){
		$r->{'SS'} = HELIX;					# r.SS = HELIX;
	}

	# if (($r->{'SS'} != HELIX) and ($r->{'SS'} != SHEET)) {	# if (r.SS != HELIX && r.SS != SHEET){
	# 	$r->{'SS'} = OTHER;					# r.SS = OTHER;
	# }

	if ($r->{'name'} eq "ALA") {					# if (r.name == "ALA"){
		$r->{'aaCode'} = ALA;					# r.aaCode = ALA;
	} elsif ($r->{'name'} eq "ARG") {				# else if (r.name == "ARG"){
		$r->{'aaCode'} = ARG;					# r.aaCode = ARG;
	} elsif ($r->{'name'} eq "ASN") {				# else if (r.name == "ASN"){
		$r->{'aaCode'} = ASN;					# r.aaCode = ASN;
	} elsif ($r->{'name'} eq "ASP") {				# else if (r.name == "ASP"){
		$r->{'aaCode'} = ASP;					# r.aaCode = ASP;
	} elsif ($r->{'name'} eq "CYS") {				# else if (r.name == "CYS"){
		$r->{'aaCode'} = CYS;					# r.aaCode = CYS;
	} elsif ($r->{'name'} eq "GLN") {				# else if (r.name == "GLN"){
		$r->{'aaCode'} = GLN;					# r.aaCode = GLN;
	} elsif ($r->{'name'} eq "GLU") {				# else if (r.name == "GLU"){
		$r->{'aaCode'} = GLU;					# r.aaCode = GLU;
	} elsif ($r->{'name'} eq "GLY") {				# else if (r.name == "GLY"){
		$r->{'aaCode'} = GLY;					# r.aaCode = GLY;
	} elsif ($r->{'name'} eq "HIS") {				# else if (r.name == "HIS"){
		$r->{'aaCode'} = HIS;					# r.aaCode = HIS;
	} elsif ($r->{'name'} eq "ILE") {				# else if (r.name == "ILE"){
		$r->{'aaCode'} = ILE;					# r.aaCode = ILE;
	} elsif ($r->{'name'} eq "LEU") {				# else if (r.name == "LEU"){
		$r->{'aaCode'} = LEU;					# r.aaCode = LEU;
	} elsif ($r->{'name'} eq "LYS") {				# else if (r.name == "LYS"){
		$r->{'aaCode'} = LYS;					# r.aaCode = LYS;
	} elsif ($r->{'name'} eq "MET") {				# else if (r.name == "MET"){
		$r->{'aaCode'} = MET;					# r.aaCode = MET;
	} elsif ($r->{'name'} eq "PHE") {				# else if (r.name == "PHE"){
		$r->{'aaCode'} = PHE;					# r.aaCode = PHE;
	} elsif ($r->{'name'} eq "PRO") {				# else if (r.name == "PRO"){
		$r->{'aaCode'} = PRO;					# r.aaCode = PRO;
	} elsif ($r->{'name'} eq "SER") {				# else if (r.name == "SER"){
		$r->{'aaCode'} = SER;					# r.aaCode = SER;
	} elsif ($r->{'name'} eq "THR") {				# else if (r.name == "THR"){
		$r->{'aaCode'} = THR;					# r.aaCode = THR;
	} elsif ($r->{'name'} eq "TRP") {				# else if (r.name == "TRP"){
		$r->{'aaCode'} = TRP;					# r.aaCode = TRP;
	} elsif ($r->{'name'} eq "TYR") {				# else if (r.name == "TYR"){
		$r->{'aaCode'} = TYR;					# r.aaCode = TYR;
	} elsif ($r->{'name'} eq "VAL") {				# else if (r.name == "VAL"){
		$r->{'aaCode'} = VAL;					# r.aaCode = VAL;
	} else {							# else{
		$r->{'aaCode'} = UNK;					# r.aaCode = UNK;	//Unknown residue
	}
}



sub PolymerBuilder_buildStructure { # Structure PolymerBuilder::buildStructure(vector<string> chains, float phi, float psi, float omega)

	my $self = shift; #C++ PolymerBuilder
	my $chains = shift; #C++ vector<string>
	my $phi = shift; #C++ float
	my $psi = shift; #C++ float
	my $omega = shift; #C++ float
	my $push_back_index;

	Structure_new(); #C++ Structure					# Structure s;
	# my $chain = 'A';							# char chain = 'A';
	my $chain = 1;
	my $i; #C++ float							# unsigned int i;
	my $chains_dot_size = $#{$chains} + 1;
	my $debug_output = ''; # debug
	for ($i = 0; $i < $chains_dot_size; $i++) {				# for (i = 0; i < chains.size(); i++){
		# create new chain						# create new chain
		my $r = Range_new(); #C++ Range					# Range r;
		# $r->{'chainID'} = ord $chain;					# r.ChainID = chain++;
		$r->{'chainID'} = $chain;
		# $chain = chr ((ord $chain) + 1);
		$chain += 1;
		my $s_dot_allResidues_dot_size = $#{$structure->{'allResidues'}} + 1;
		$r->{'min'} = $s_dot_allResidues_dot_size;			# r.min = s.allResidues.size();
		my $chains_i_dot_length = length($chains->[$i]);
		$r->{'max'} = $r->{'min'} + $chains_i_dot_length - 1;		# r.max = r.min + (unsigned int)chains[i].length() - 1;
		# specify as HELIX						# specify as HELIX
		$helices->{'chainID'} = $r->{'chainID'};			# helices.ChainID = r.ChainID;
		$helices->{'min'} = $r->{'min'}; #C++ Range			# helices.min = r.min;
		$helices->{'max'} = $r->{'max'}; #C++ Range			# helices.max = r.max;
		# add chain to structure
		$push_back_index = $#{$structure->{'chains'}} + 1;
		$structure->{'chains'}->[$push_back_index] = $r; #C++ Range	# s.chains.push_back(r);
		$push_back_index = $#{$structure->{'HelixInfo'}} + 1;
		$structure->{'HelixInfo'}->[$push_back_index] = $helices; #C++ Range # s.HelixInfo.push_back(helices);
		my $j; #C++ unsigned int
		$chains_i_dot_length = length($chains->[$i]);
		for ($j = 0; $j < $chains_i_dot_length; $j++) {			# for (unsigned int j = 0; j < chains[i].length(); j++){
			my $res = Residue_new(); #C++ Residue			# Residue res;
			my $char_for_charToRes = substr( $chains->[$i], $j, 1 );
			$res->{'name'} = PolymerBuilder_charToRes($char_for_charToRes); # res.name = charToRes(chains[i][j]);
			$res->{'chainID'} = $r->{'chainID'};			# res.chainID = r.ChainID;
			$res->{'resSeq'} = $j;					# res.resSeq = j;
			PolymerBuilder_addResidue( $res, $phi, $psi, $omega ); 	# addResidue(res, phi, psi, omega, s);
		}
	}

	# find xyz extremities
	my $s_dot_allAtoms_dot_size = $#{$structure->{'allAtoms'}} + 1;
	for ($i = 0; $i < $s_dot_allAtoms_dot_size; $i++) {									# for(i = 0; i < s.allAtoms.size(); i++){
		if ($structure->{'allAtoms'}->[$i]->{'coords'}->{'x'} < $structure->{'xyzExtremities'}->[0]->{'f1'}) {		# if(s.allAtoms[i].coords.x < s.xyzExtremities[0].f1){
			$structure->{'xyzExtremities'}->[0]->{'f1'} = $structure->{'allAtoms'}->[$i]->{'coords'}->{'x'};	# s.xyzExtremities[0].f1 = s.allAtoms[i].coords.x;
		}
		if ($structure->{'allAtoms'}->[$i]->{'coords'}->{'x'} > $structure->{'xyzExtremities'}->[0]->{'f2'}) {		# if(s.allAtoms[i].coords.x > s.xyzExtremities[0].f2){
			$structure->{'xyzExtremities'}->[0]->{'f2'} = $structure->{'allAtoms'}->[$i]->{'coords'}->{'x'};	# s.xyzExtremities[0].f2 = s.allAtoms[i].coords.x;
		}
		if ($structure->{'allAtoms'}->[$i]->{'coords'}->{'y'} < $structure->{'xyzExtremities'}->[1]->{'f1'}) {		# if(s.allAtoms[i].coords.y < s.xyzExtremities[1].f1){
			$structure->{'xyzExtremities'}->[1]->{'f1'} = $structure->{'allAtoms'}->[$i]->{'coords'}->{'y'};	# s.xyzExtremities[1].f1 = s.allAtoms[i].coords.y;
		}
		if ($structure->{'allAtoms'}->[$i]->{'coords'}->{'y'} > $structure->{'xyzExtremities'}->[1]->{'f2'}) {		# if(s.allAtoms[i].coords.y > s.xyzExtremities[1].f2){
			$structure->{'xyzExtremities'}->[1]->{'f2'} = $structure->{'allAtoms'}->[$i]->{'coords'}->{'y'};	# s.xyzExtremities[1].f2 = s.allAtoms[i].coords.y;
		}
		if ($structure->{'allAtoms'}->[$i]->{'coords'}->{'z'} < $structure->{'xyzExtremities'}->[2]->{'f1'}) {		# if(s.allAtoms[i].coords.z < s.xyzExtremities[2].f1){
			$structure->{'xyzExtremities'}->[2]->{'f1'} = $structure->{'allAtoms'}->[$i]->{'coords'}->{'z'};	# s.xyzExtremities[2].f1 = s.allAtoms[i].coords.z;
		}
		if ($structure->{'allAtoms'}->[$i]->{'coords'}->{'z'} > $structure->{'xyzExtremities'}->[2]->{'f2'}) {		# if(s.allAtoms[i].coords.z > s.xyzExtremities[2].f2){
			$structure->{'xyzExtremities'}->[2]->{'f2'} = $structure->{'allAtoms'}->[$i]->{'coords'}->{'z'};	# s.xyzExtremities[2].f2 = s.allAtoms[i].coords.z;
		}
	}

	if (defined($input_debug_flag)) {
		my $output_pdb_file = 'valpred.pdb';
		open( PDBFILE, ">$output_pdb_file") or
			die "Cannot open $output_pdb_file for writing : $!\n";
		my $output_lines = '';
		$output_lines .= "HEADER    VALPRED PROTEIN                         01-JAN-00   1234              \n";
		$output_lines .= "TITLE     STRUCTURE PRODUCED FROM VALPRED                                       \n";

		my $st_dot_allAtoms_dot_size = $#{$structure->{'allAtoms'}} + 1;
		for (my $e = 0; $e < $st_dot_allAtoms_dot_size; $e++) { # for (unsigned int e = 0 ; e < st.allAtoms.size(); e++){

			my $atomID = $structure->{'allAtoms'}->[$e]->{'atomID'};
			my $atom_name = $structure->{'allAtoms'}->[$e]->{'name'};
			my $resSeq = $structure->{'allAtoms'}->[$e]->{'resSeq'};
			my $residue_number = $resSeq + 1;
			my $res_name = $structure->{'allResidues'}->[$resSeq]->{'name'};
			my $coords_x = $structure->{'allAtoms'}->[$e]->{'coords'}->{'x'};
			my $coords_y = $structure->{'allAtoms'}->[$e]->{'coords'}->{'y'};
			my $coords_z = $structure->{'allAtoms'}->[$e]->{'coords'}->{'z'};
			my $atom_type = substr( $atom_name, 0, 1 );

			my $print_ATOM = sprintf( "%-6s", 'ATOM' );
			my $print_atomID = sprintf( "% 5d", $atomID );
			my $print_atom_name = sprintf( "%-4s", $atom_name );
			my $print_res_name = sprintf( "%-3s", $res_name );
			my $print_resSeq = sprintf( "% 4d", $residue_number );
			my $print_coords_x = sprintf( "% 8.3f", $coords_x );
			my $print_coords_y = sprintf( "% 8.3f", $coords_y );
			my $print_coords_z = sprintf( "% 8.3f", $coords_z );
			my $print_bit = '  1.00 33.00           ';
			my $print_atom_type = sprintf( "%-1s", $atom_type );

			my $this_line = $print_ATOM . $print_atomID . '  ' . $print_atom_name . $print_res_name . ' A' . $print_resSeq . '    ' .
				$print_coords_x . $print_coords_y . $print_coords_z . $print_bit . $print_atom_type;

			$output_lines .= "$this_line\n";
		}

		print PDBFILE $output_lines;
		close PDBFILE;
	}
	# ///create PDB file for debugging
	# //ofstream outFile;
	# //outFile.open("./helixtest.pdb",ios::out);
	# //if (outFile.fail()){
	# //	AfxMessageBox("debug error: Could not open ./helixtest.pdb");
	# //	return s;
	# //}
	# //for (unsigned int i = 0; i < s.allAtoms.size(); i++){
	# //	outFile<<formatLine(&(s.allAtoms[i]),&s);
	# //}
	# ///create PDB for debugging//
	# string PolymerBuilder::formatLine(Atom *a, Structure *s){
	# //	ATOM      1   N  ASN     1       1.458   0.000   0.000  1.00  0.00
	# //	ATOM   2220  CG2 VAL B 137      52.871   5.236  16.888  1.00 12.02           C  

	#	string output,strBuf;
	#	char *cStarBuf = new char[15];
	#	output += "ATOM";
	#	sprintf(cStarBuf,"%d",a->atomID);
	#	strBuf =  cStarBuf;
	#	while(strBuf.length() < 7){
	#		strBuf = " " + strBuf;
	#	}
	#	output+= strBuf + "  ";
	#	strBuf = a->name;
	#	while(strBuf.length() < 4){
	#		strBuf += " ";
	#	}
	#	output+= strBuf;
	#	output += s->allResidues[a->resSeq].name + "  ";
	#	sprintf(cStarBuf,"%d",a->resSeq);
	#	strBuf =  cStarBuf;
	#	while(strBuf.length() < 4){
	#		strBuf = " " + strBuf;
	#	}
	#	output+= strBuf + "    ";
	#
	#	unsigned int dec = 0;
	#	sprintf(cStarBuf,"%f",a->coords.x);
	#	strBuf =  cStarBuf;
	#	while(strBuf[dec]!= '.'){
	#		dec++;
	#	}
	#	while(dec < strBuf.length() - 4){
	#		strBuf = strBuf.substr(0,strBuf.length() - 1);
	#	}
	#	while(strBuf.length() < 8){
	#		strBuf = " " + strBuf;
	#	}
	#	output+= strBuf;
	#
	#	sprintf(cStarBuf,"%f",a->coords.y);
	#	strBuf =  cStarBuf;
	#	dec = 0;
	#	while(strBuf[dec]!= '.'){
	#		dec++;
	#	}
	#	while(dec < strBuf.length() - 4){
	#		strBuf = strBuf.substr(0,strBuf.length() - 1);
	#	}
	#	while(strBuf.length() < 8){
	#		strBuf = " " + strBuf;
	#	}
	#	output+= strBuf;
	#
	#	sprintf(cStarBuf,"%f",a->coords.z);
	#	strBuf =  cStarBuf;
	#	dec = 0;
	#	while(strBuf[dec]!= '.'){
	#		dec++;
	#	}
	#	while(dec < strBuf.length() - 4){
	#		strBuf = strBuf.substr(0,strBuf.length() - 1);
	#	}
	#	while(strBuf.length() < 8){
	#		strBuf = " " + strBuf;
	#	}
	#	output+= strBuf;
	#
	#	output += "\n";
	#	delete [] cStarBuf;
	#	return output;
	#}

	return;			# return s;
}



sub PolymerBuilder_charToRes { # string PolymerBuilder::charToRes(char c){

	my $c = shift; #C++ char
	my $rtn;
	if (defined $charToRes{$c}) {
		$rtn = $charToRes{$c};
	} else {
		$rtn = 'UNK';
	}

	return $rtn; #C++ string
}



sub PolymerBuilder_createSideChain { #C++ void PolymerBuilder::createSideChain(Residue &r, Structure &s){

	my $r = shift; #C++ Residue

	my $m = Mat44_new();							# Mat44f m;

	my $modelN = Vec3_new(); #C++ Vec3f					# Vec3f modelN;
	my $modelCA = Vec3_new(); #C++ Vec3f					# Vec3f modelCA;
	my $modelC = Vec3_new(); #C++ Vec3f					# Vec3f modelC;
	my $modelO = Vec3_new(); #C++ Vec3f					# Vec3f modelO;

	my $a = Atom_new(); #C++ Atom

	my $r_name = $r->{'name'};						# filepath = currentDir + "/Amino Acids/" + r.name + ".pdb";
	my $num_amino_acid_pdb_lines = $#{$amino_acid_pdbs->{$r_name}};
	for ( my $i = 0; $i <= $num_amino_acid_pdb_lines; $i++ ) {		# while(!(f.eof())){

		my $name = $amino_acid_pdbs->{$r_name}->[$i]->{'name'};		# string name = buffer.substr(12,4);
		my $x = $amino_acid_pdbs->{$r_name}->[$i]->{'x'};		# float x = (float)(atof(buffer.substr(30, 8).c_str()));
		my $y = $amino_acid_pdbs->{$r_name}->[$i]->{'y'};		# float y = (float)(atof(buffer.substr(38, 8).c_str()));
		my $z = $amino_acid_pdbs->{$r_name}->[$i]->{'z'};		# float z = (float)(atof(buffer.substr(46, 8).c_str()));

		if ($name eq 'N') {						# if (name=="N") 
			$modelN = Vec3_set3( $modelN, $x, $y, $z );		# modelN.set(x, y, z); continue;
		} elsif ($name eq 'CA') {					# else if (name=="CA") 
			$modelCA = Vec3_set3( $modelCA, $x, $y, $z );		# modelCA.set(x, y, z); continue;
		} elsif ($name eq 'C') {					# else if (name=="C")
			$modelC = Vec3_set3( $modelC, $x, $y, $z );		# modelC.set(x, y, z); continue;
		} elsif ($name eq 'O') {					# else if (name=="O")
			$modelO = Vec3_set3( $modelO, $x, $y, $z );		# modelO.set(x, y, z);
			my $r_N = $r->{'N'};
			my $r_CA = $r->{'CA'};
			my $r_C = $r->{'C'};
			# m = xformMat(modelN, modelCA, modelC, s.allAtoms[r.N].coords, s.allAtoms[r.CA].coords,s.allAtoms[r.C].coords);
			$m = geometry_xformMat( $modelN, $modelCA, $modelC, $structure->{'allAtoms'}->[$r_N]->{'coords'}, $structure->{'allAtoms'}->[$r_CA]->{'coords'}, $structure->{'allAtoms'}->[$r_C]->{'coords'} );
										# continue;
		} else {							# else
			my $coord = Vec3_new();					# Vec3f coord;
			$coord = Vec3_set3( $coord, $x, $y, $z );		# coord.set(x, y, z);
			$coord = Mat44_multVec3d( $m, $coord );			# coord = m.multVec3d(coord);
			$a = Atom_new(); #C++ Atom				# Atom a;
			$a->{'name'} = $name;					# a.name = name;
			$a->{'coords'} = $coord;				# a.coords = coord;
			my $s_dot_allAtoms_dot_size = $#{$structure->{'allAtoms'}} + 1;
			$a->{'atomID'} = $s_dot_allAtoms_dot_size;		# a.atomID = (unsigned int)s.allAtoms.size();
			$a->{'index'} = $s_dot_allAtoms_dot_size;		# a.index = (unsigned int)s.allAtoms.size();
			$a->{'chainID'} = ' ';					# a.chainID = ' ';
			$a->{'resSeq'} = $r->{'resSeq'};			# a.resSeq = r.resSeq;
			Atom_setAtomType( $a, $r, $s_dot_allAtoms_dot_size );	# a.setAtomType(r,(unsigned int)(s.allAtoms.size()));
			my $push_back_index = $#{$structure->{'allAtoms'}} + 1;
			$structure->{'allAtoms'}->[$push_back_index] = $a;	# s.allAtoms.push_back(a);
		}
	}

	return;
}



sub PolymerBuilder_degToRad { # float PolymerBuilder::degToRad(float deg){

	my $deg = shift;
	my $deg_to_rad = $deg / 180 * PI;		# return deg/180*PI;
	return $deg_to_rad;
}



sub PolymerBuilder_new {

	my $p;

	return $p; #C++ PolymerBuilder
}



sub predict_amphipathic_helices {

	# Calculations for the algorithm to find amphipathic helices :
	# ============================================================
	#
	# all not-buried polar points vs. all not-buried points 
	#				(which gives all not-buried hydrophobic points)
	# if a point is not buried, add it to amphipathic-calc not-buried points count
	# if a point is not buried, see if it is truly polar. 
	#				if so, add it to the amphipathic-calc not-buried polar points.
	#
	# 'solventSurface' solvent surface area = the surface area of this residue when not in a structure/helix
	#				= add up the surface area of each atom in this residue when in this particular amino acid
	# ASA = accessible surface area = ASA for this residue in this particular helix or 3D protein structure
	#		= ( (num points - num buried points) / num points ) * solvent surface area
	#		= ( num not-buried points / num points ) * solvent surface area
	# side_ASA = same as above = side chain ASA = CA.ASA + ASA of each side chain atom, 
	#			in this particular helix or 3D protein structure
	# amphipathic_calc_side_chain_not_buried_polar_ASA = side_ASA * 
	#			( amphipathic_calc_not_buried_is_polar_points count / amphipathic_calc_not_buried_points count) 

	# In an amphipathic helix, we expect one side to be (almost) completely non-polar, 
	# to have (almost) no polar not-buried side chain surface area.
	# And we expect the other side to have some (or a lot) of polar not-buried side chain surface area.
	# Below, 'SCpolarSASA' refers to the polar not-buried side chain surface area = 'amphipathic_calc_side_chain_not_buried_polar_ASA'

	##### emma, put these values into parameters/options

	# my @residues_per_turn_array = (3, 3.6, 4.1); # 3 = 3-10 helix, 3.6 = alpha helix, 4.1 = pi-helix

	# my @residues_per_turn_array = (2.9, 3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4, 4.1, 4.2); 
	# my @residue_starting_point_array = (2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9);

	# my @residues_per_turn_array = (3.7, 3.71, 3.72, 3.73, 3.74, 3.75, 3.76, 3.77, 3.78, 3.79, 3.8); 
	# my @residue_starting_point_array = ( 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 
	#	2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 
	#	3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 
	#	4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 
	#	5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9,
	#	6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 
	#	7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 
	#	8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, );

	my @residues_per_turn_array = (2.9, 3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4, 4.1, 4.2); 
	my @residue_starting_point_array = (2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9);

	my @algorithm_residue_point_weights = (0.3, 0.4, 0.3); # :::::-::::: 
	# eg. (0.1, 0.8, 0.1) means SCpolarSASA for residue n is (0.1 x SCpolarSASA_of_n-1) + (0.8 x SCpolarSASA_of_n) + (0.1 x SCpolarSASA_of_n+1)

	my $algorithm_sliding_window = 5;

	my $max_SCpolarSASA_for_hydrophobic_side = 20; # also called the amphipathic_side
	# my $min_SCpolarSASA_for_not_hydrophobic_side = 1; # also called the opposite_side

	# my $min_num_residues_for_amphipathic_helix = 15;
	# my $avg_num_residues_per_turn = 3.6;
	# my $min_num_turns_for_amphipathic_helix = $min_num_residues_for_amphipathic_helix / $avg_num_residues_per_turn;

	my $min_num_turns_for_amphipathic_helix = 5;

	my $max_hydrophobic_side_over_opposite_side = 0.8;

	my $st_dot_chains_dot_size = $#{$structure->{'chains'}} + 1;
	for ( my $c = 0; $c < $st_dot_chains_dot_size; $c++ ) {	
		my $residues = Structure_getResChain_int( $c );
		my $residues_dot_size = $#{$residues} + 1;
		# for (my $i = 0; $i < $residues_dot_size; $i++) {

		my @amphipathic_helices_array;
		my $amphipathic_helices = \@amphipathic_helices_array;

		my @residue_is_amphipathic;
		for (my $r = 1; $r <= $residues_dot_size; $r++) {
			$residue_is_amphipathic[$r] = 0;
		}

print "\n";#####
for (my $ii = 0; $ii < $residues_dot_size; $ii++) {#####
my $rr = $ii + 1;#####
my $rr_SASA = $structure->{'allResidues'}->[$ii]->{'amphipathic_calc_side_chain_not_buried_polar_ASA'};#####
$rr = sprintf( "% 5d", $rr);#####
$rr_SASA = sprintf( "% 5d", $rr_SASA);#####
print "rr=$rr, amphipathic_calc_side_chain_not_buried_polar_ASA=$rr_SASA\n";#####
}#####

		my @start_end_array_for_this_chain;

		foreach my $this_residues_per_turn (@residues_per_turn_array) {

			foreach my $this_residue_starting_point (@residue_starting_point_array) {
print "\nthis_residues_per_turn=$this_residues_per_turn, this_residue_starting_point=$this_residue_starting_point\n\n";#####

				# Starting from this starting point in the helix (which is a residue),
				# see if this helix is amphipathic, by seeing if
				# going up this side of the helix is hydrophobic and 
				# going up the opposite side of the helix is not hydrophobic.

				my @residue_points_of_the_amphipathic_face;	# 2.0, 5.6, 9.2, 12.8, 16.4, 20.0, -1
				my @residue_points_of_the_opposite_face;	# -1,  3.8, 7.4, 11.0, 14.6, 18.2, 21.8
				# corresponding opposite_face element is below the corresponding amphipathic_face element, 
				# so has 1 element more for last opposite_face element above the last amphipathic_face element.
				# a value of -1 means that this residue point doesn't exist.
				my @residue_points_weighted_SCpolarSASA_of_the_amphipathic_face;
				my @residue_points_weighted_SCpolarSASA_of_the_opposite_face;
				my @residue_points_sliding_window_weighted_SCpolarSASA_of_the_amphipathic_face;
				my @residue_points_sliding_window_weighted_SCpolarSASA_of_the_opposite_face;
				my @residue_points_amphipathic_face_is_hydrophobic;
				my @residue_points_opposite_face_is_hydrophobic;
				my @residue_points_is_amphipathic_helix;

				# First, list all the residues points to be inspected.
				# A residue point may be a residue or may be part way between 2 residues.
				# list all the residue points to be inspected on the investigate-if-hydrophobic side 
				# (called the amphipathic_face),
				# and list all the residue points to be inspected on the investigate-if-not-hydrophobic side
				# (called the opposite_face).

				my $got_to_end_of_sequence = 0;
				my $this_amphipathic_face_residue_upto = $this_residue_starting_point;
				my $this_opposite_face_residue_upto = $this_residue_starting_point - ($this_residues_per_turn / 2);
				while ($got_to_end_of_sequence == 0) {

					if ($this_amphipathic_face_residue_upto >= 1) {
						push( @residue_points_of_the_amphipathic_face, $this_amphipathic_face_residue_upto );
					} else {
						push( @residue_points_of_the_amphipathic_face, -1 );
					}
					if ($this_opposite_face_residue_upto >= 1) {
						push( @residue_points_of_the_opposite_face, $this_opposite_face_residue_upto );
					} else {
						push( @residue_points_of_the_opposite_face, -1 );
					}

					$this_amphipathic_face_residue_upto += $this_residues_per_turn;
					$this_opposite_face_residue_upto += $this_residues_per_turn;
					
					if ($this_amphipathic_face_residue_upto > $residues_dot_size) {

						$got_to_end_of_sequence = 1;

						push( @residue_points_of_the_amphipathic_face, -1 );
						if ($this_opposite_face_residue_upto <= $residues_dot_size) {
							push( @residue_points_of_the_opposite_face, $this_opposite_face_residue_upto );
						} else {
							push( @residue_points_of_the_opposite_face, -1 );
						}
					}
				}

				# Next, for each residue point to be inspected,
				# inspect it to see if it is hydrophobic or not hydrophobic.
				# The hydrophobicity is measured by the area of the side chain not-buried surface area 
				# that is polar (called SCpolarSASA here).
				# But don't just use this one residue (and we might be partway between 2 residues anyway).
				# Use the polar side chain SASA (SCpolarSASA) of this residue 
				# and a tiny bit of that of its neighbours.
				# If we are partway between 2 residues, use part of the SCpolarSASA of each one,
				# according to how close each residue is to this residue point.
				# (Eg. if this is residue point is 2.3, then use more of residue 2's polar side chain SASA (SCpolarSASA)
				# than of residue 3's polar side chain SASA (SCpolarSASA).
				# Record this residue point's polar side chain SASA (SCpolarSASA) for each residue point
				# going up the helix on the inspect-if-hydrophic-side (called amphipathic_face)
				# and going up the helix on the inpsect-if-not-hydrophobic-side (called opposite_face).

				for ( my $i = 0; $i < @residue_points_of_the_amphipathic_face; $i++ ) {
					my $this_amphipathic_face_residue_point_upto = $residue_points_of_the_amphipathic_face[$i]; # eg. = 2.1

					if (($this_amphipathic_face_residue_point_upto >= 1) && ($this_amphipathic_face_residue_point_upto <= $residues_dot_size)) {
						my $pt_1 = $this_amphipathic_face_residue_point_upto - 1;
						my $pt_2 = $this_amphipathic_face_residue_point_upto;
						my $pt_3 = $this_amphipathic_face_residue_point_upto + 1;
						if ($pt_1 < 1) {
							$pt_1 = 1;
						}
						if ($pt_3 > $residues_dot_size) {
							$pt_3 = $residues_dot_size;
						}
						my @neighbouring_residue_points = ( $pt_1, $pt_2, $pt_3 );
						my @neighbouring_residue_points_unweighted_SCpolarSASA;
						for ( my $j = 0; $j < @neighbouring_residue_points; $j++ ) {
							my $this_residue_point = $neighbouring_residue_points[$j]; # eg. 2.1
							my $this_residue_point_weight = $algorithm_residue_point_weights[$j];

							my $lower_residue = int($this_residue_point); # eg. = 2
							my $upper_residue = $lower_residue + 1; # eg. = 3
							if ($lower_residue == $this_residue_point) {
								$upper_residue = $this_residue_point;
							}
							if ($lower_residue < 1) {
								$lower_residue = 1;
							}
							if ($upper_residue > $residues_dot_size) {
								$upper_residue = $residues_dot_size;
							}
							# when accessing $structure->{'allResidues'}-> array, use index of $residue_num - 1, 
							# because that array goes from 0..($residues_dot_size - 1)
							# whereas this algorithm goes from 1..$residues_dot_size
							# my $lower_residue_SCpolarSASA = $structure->{'allResidues'}->[$lower_residue - 1]->{'amphipathic_calc_side_chain_not_buried_polar_ASA'};
							# my $upper_residue_SCpolarSASA = $structure->{'allResidues'}->[$upper_residue - 1]->{'amphipathic_calc_side_chain_not_buried_polar_ASA'};
							my $lower_residue_SCpolarSASA = $residues->[$lower_residue - 1]->{'amphipathic_calc_side_chain_not_buried_polar_ASA'};
							my $upper_residue_SCpolarSASA = $residues->[$upper_residue - 1]->{'amphipathic_calc_side_chain_not_buried_polar_ASA'};
							my $lower_weight = $upper_residue - $this_residue_point; # eg. = 3 - 2.1 = 0.9
							my $upper_weight = $this_residue_point - $lower_residue; # eg. = 2.1 - 2 = 0.1
							if ($lower_weight == 0) {
								$upper_weight = 1;
							}
							my $unweighted_SCpolarSASA = ($lower_residue_SCpolarSASA * $lower_weight) + ($upper_residue_SCpolarSASA * $upper_weight);
							# my $unweighted_SCpolarSASA = $lower_residue_SCpolarSASA;
							# if ($upper_residue_SCpolarSASA < $lower_residue_SCpolarSASA) {
							#	$unweighted_SCpolarSASA = $upper_residue_SCpolarSASA;
							# }
							$neighbouring_residue_points_unweighted_SCpolarSASA[$j] = $unweighted_SCpolarSASA;
						}
						my $this_amphipathic_face_residue_point_upto_SCpolarSASA = 
							($neighbouring_residue_points_unweighted_SCpolarSASA[0] * $algorithm_residue_point_weights[0]) +
							($neighbouring_residue_points_unweighted_SCpolarSASA[1] * $algorithm_residue_point_weights[1]) +
							($neighbouring_residue_points_unweighted_SCpolarSASA[2] * $algorithm_residue_point_weights[2]);
						$residue_points_weighted_SCpolarSASA_of_the_amphipathic_face[$i] = $this_amphipathic_face_residue_point_upto_SCpolarSASA;
					} else {
						$residue_points_weighted_SCpolarSASA_of_the_amphipathic_face[$i] = 0;
					}
#if ($residue_points_of_the_amphipathic_face[$i] > 20) {#####
#my $print_i = sprintf( "%-5s", $i);#####
#my $print_point = sprintf( "% 7.1f", $residue_points_of_the_amphipathic_face[$i]);#####
#my $print_SASA = sprintf( "% 5d", $residue_points_weighted_SCpolarSASA_of_the_amphipathic_face[$i]);#####
#print "     i=$print_i, amphipathic_face_residue_point=" . $print_point . ", weighted_SCpolarSASA=" . $print_SASA . "\n";#####
#}#####
				}
#print "     ==================================//==================================\n";#####

				for ( my $i = 0; $i < @residue_points_of_the_opposite_face; $i++ ) {
					my $this_opposite_face_residue_point_upto = $residue_points_of_the_opposite_face[$i]; # eg. = 2.1

					if (($this_opposite_face_residue_point_upto >= 1) && ($this_opposite_face_residue_point_upto <= $residues_dot_size)) {
						my $pt_1 = $this_opposite_face_residue_point_upto - 1;
						my $pt_2 = $this_opposite_face_residue_point_upto;
						my $pt_3 = $this_opposite_face_residue_point_upto + 1;
						if ($pt_1 < 1) {
							$pt_1 = 1;
						}
						if ($pt_3 > $residues_dot_size) {
							$pt_3 = $residues_dot_size;
						}
						my @neighbouring_residue_points = ( $pt_1, $pt_2, $pt_3 );
						my @neighbouring_residue_points_unweighted_SCpolarSASA;
						for ( my $j = 0; $j < @neighbouring_residue_points; $j++ ) {
							my $this_residue_point = $neighbouring_residue_points[$j]; # eg. 28
							my $this_residue_point_weight = $algorithm_residue_point_weights[$j];

							my $lower_residue = int($this_residue_point); # eg. = 28
							my $upper_residue = $lower_residue + 1; # eg. = 29
							if ($lower_residue == $this_residue_point) {
								$upper_residue = $this_residue_point; # eg. = 28
							}
							if ($lower_residue < 1) {
								$lower_residue = 1;
							}
							if ($upper_residue > $residues_dot_size) {
								$upper_residue = $residues_dot_size;
							}
							# when accessing $structure->{'allResidues'}-> array, use index of $residue_num - 1, 
							# because that array goes from 0..($residues_dot_size - 1)
							# whereas this algorithm goes from 1..$residues_dot_size
							my $lower_residue_SCpolarSASA = $structure->{'allResidues'}->[$lower_residue - 1]->{'amphipathic_calc_side_chain_not_buried_polar_ASA'};
							my $upper_residue_SCpolarSASA = $structure->{'allResidues'}->[$upper_residue - 1]->{'amphipathic_calc_side_chain_not_buried_polar_ASA'};
							my $lower_weight = $upper_residue - $this_residue_point; # eg. = 28 - 28 = 0
							my $upper_weight = $this_residue_point - $lower_residue; # eg. = 28 - 28 = 0
							if ($lower_weight == 0) {
								$upper_weight = 1; # eg. = 1
							}
							my $unweighted_SCpolarSASA = ($lower_residue_SCpolarSASA * $lower_weight) + ($upper_residue_SCpolarSASA * $upper_weight);
							# my $unweighted_SCpolarSASA = $lower_residue_SCpolarSASA;
							# if ($upper_residue_SCpolarSASA < $lower_residue_SCpolarSASA) {
							#	$unweighted_SCpolarSASA = $upper_residue_SCpolarSASA;
							# }
							$neighbouring_residue_points_unweighted_SCpolarSASA[$j] = $unweighted_SCpolarSASA
						}
						my $this_opposite_face_residue_point_upto_SCpolarSASA = 
							($neighbouring_residue_points_unweighted_SCpolarSASA[0] * $algorithm_residue_point_weights[0]) +
							($neighbouring_residue_points_unweighted_SCpolarSASA[1] * $algorithm_residue_point_weights[1]) +
							($neighbouring_residue_points_unweighted_SCpolarSASA[2] * $algorithm_residue_point_weights[2]);
						$residue_points_weighted_SCpolarSASA_of_the_opposite_face[$i] = $this_opposite_face_residue_point_upto_SCpolarSASA;
					} else {
						$residue_points_weighted_SCpolarSASA_of_the_opposite_face[$i] = 0;
					}
#if ($residue_points_of_the_opposite_face[$i] > 20) {#####
#my $print_i = sprintf( "%-5s", $i);#####
#my $print_point = sprintf( "% 7.1f", $residue_points_of_the_opposite_face[$i]);#####
#my $print_SASA = sprintf( "% 5d", $residue_points_weighted_SCpolarSASA_of_the_opposite_face[$i]);#####
#print "     i=$print_i, opposite_face_residue_point=" . $print_point . ", weighted_SCpolarSASA=" . $print_SASA . "\n";#####
#}#####
				}

				# Now that we have the polar side chain SASA (SCpolarSASA : polar side chain not-buried surface area)
				# for each residue point going up both sides of the helix,
				# don't just use this SCpolarSASA to decide if one side is hydrophobic and the other not-hydrophobic.
				# Maybe there is a bit of polar surface area on the inspect-if-hydrophobic side
				# and it is tolerated in an amphiphatic helix's hydrophobic side
				# because it is surrounded by all hydrophobic surface area below and above.
				# So, use a sliding window for the inspect-if-hydrophobic side (called amphipathic_face)
				# and for the inspect-if-not-hydrophobic side (called opposite_face).
				# Store the found amphipathic helix segments (stored in @start_end_array_for_this_chain)

				for ( my $i = 0; $i < @residue_points_weighted_SCpolarSASA_of_the_amphipathic_face; $i++ ) {
					my $this_amphipathic_face_residue_point_upto_weighted_SCpolarSASA = $residue_points_weighted_SCpolarSASA_of_the_amphipathic_face[$i];
				
					my $sliding_window_sum = 0;
				
					my $start_window_index = $i - int($algorithm_sliding_window / 2);
					my $end_window_index = $i + int($algorithm_sliding_window / 2);
					my $use_half_of_start_and_end_index = 0;
					if ($algorithm_sliding_window == int($algorithm_sliding_window / 2)) {
						$use_half_of_start_and_end_index = 1;
					}
				
					for ( my $j = $start_window_index; $j <= $end_window_index; $j++ ) {
						my $this_SCpolarSASA = 0;
						if (($j >= 0) && ($j < @residue_points_weighted_SCpolarSASA_of_the_amphipathic_face)) {
							$this_SCpolarSASA = $residue_points_weighted_SCpolarSASA_of_the_amphipathic_face[$j];
						}
						if ($use_half_of_start_and_end_index == 1) {
							if (($j == $start_window_index) || ($j == $end_window_index)) {
								$this_SCpolarSASA = $this_SCpolarSASA / 2;
							}
						}
						
						$sliding_window_sum += $this_SCpolarSASA;
					}
				
					$residue_points_sliding_window_weighted_SCpolarSASA_of_the_amphipathic_face[$i] = $sliding_window_sum / $algorithm_sliding_window;
				}
				
				for ( my $i = 0; $i < @residue_points_weighted_SCpolarSASA_of_the_opposite_face; $i++ ) {
					my $this_opposite_face_residue_point_upto_weighted_SCpolarSASA = $residue_points_weighted_SCpolarSASA_of_the_opposite_face[$i];
				
					my $sliding_window_sum = 0;
				
					my $start_window_index = $i - int($algorithm_sliding_window / 2);
					my $end_window_index = $i + int($algorithm_sliding_window / 2);
					my $use_half_of_start_and_end_index = 0;
					if ($algorithm_sliding_window == int($algorithm_sliding_window / 2)) {
						$use_half_of_start_and_end_index = 1;
					}
				
					for ( my $j = $start_window_index; $j <= $end_window_index; $j++ ) {
						my $this_SCpolarSASA = 0;
						if (($j >= 0) && ($j < @residue_points_weighted_SCpolarSASA_of_the_opposite_face)) {
							$this_SCpolarSASA = $residue_points_weighted_SCpolarSASA_of_the_opposite_face[$j];
						}
						if ($use_half_of_start_and_end_index == 1) {
							if (($j == $start_window_index) || ($j == $end_window_index)) {
								$this_SCpolarSASA = $this_SCpolarSASA / 2;
							}
						}
						
						$sliding_window_sum += $this_SCpolarSASA;
					}
				
					$residue_points_sliding_window_weighted_SCpolarSASA_of_the_opposite_face[$i] = $sliding_window_sum / $algorithm_sliding_window;
				}

				# Go up the helix on the inspect-if-hydrophobic side and inpsect-if-not-hydrophobic side
				# to see if the sliding window value for respects the max/min polar values allowed 
				# for amphipathic helices.
				# Below, SCpolarSASA refers to the sliding window polar not-buried surface area value at a residue point.
				# And remember, for these corresponding amphipathic_face/opposite_face pairs,
				# the opposite_face residue point is half a turn below the amphipathic_face residue point in the helix.

				my @residue_point_is_amphipathic_helix;
				for ( my $r = 0; $r < @residue_points_of_the_amphipathic_face; $r++ ) {
					$residue_point_is_amphipathic_helix[$r] = 0;
				}
				
				for ( my $i = 0; $i < @residue_points_sliding_window_weighted_SCpolarSASA_of_the_amphipathic_face; $i++ ) {
				#for ( my $i = 0; $i < @residue_points_weighted_SCpolarSASA_of_the_amphipathic_face; $i++ ) {
				
					my $this_amphipathic_SCpolarSASA = $residue_points_sliding_window_weighted_SCpolarSASA_of_the_amphipathic_face[$i];
				
					if ($this_amphipathic_SCpolarSASA <= $max_SCpolarSASA_for_hydrophobic_side) {
						my $opposite_side_is_not_hydrophobic = 0;
						my $start_j = int($i - ($min_num_turns_for_amphipathic_helix / 2)) + 1;
						$start_j += 1; # opposite_face is half a turn below the corresponding amphipathic_face, so start half a turn above the corresponding amphipathic_face instead
						my $end_j = $i + ($min_num_turns_for_amphipathic_helix / 2);
						if ($end_j == int($end_j)) {
							$end_j = int($end_j) - 1;
						} else {
							$end_j = int($end_j);
						}
						my $avg_opposite_SCpolarSASA = 0;
						my $avg_opposite_SCpolarSASA_count = 0;
						for ( my $j = $start_j; $j <= $end_j; $j++ ) {
							if (($j >= 0) && ($j < @residue_points_sliding_window_weighted_SCpolarSASA_of_the_amphipathic_face)) {
							#if (($j >= 0) && ($j < @residue_points_weighted_SCpolarSASA_of_the_amphipathic_face)) {
								my $this_opposite_SCpolarSASA = $residue_points_sliding_window_weighted_SCpolarSASA_of_the_opposite_face[$j];
								#my $this_opposite_SCpolarSASA = $residue_points_weighted_SCpolarSASA_of_the_opposite_face[$j];
								if ($this_opposite_SCpolarSASA > $max_SCpolarSASA_for_hydrophobic_side) {
									$opposite_side_is_not_hydrophobic = 1;
#print "     i=$i, j=$j, r=" . $residue_points_of_the_amphipathic_face[$i] . ", this_opposite_SCpolarSASA=$this_opposite_SCpolarSASA, max_SCpolarSASA_for_hydrophobic_side=$max_SCpolarSASA_for_hydrophobic_side\n";#####
								}
								$avg_opposite_SCpolarSASA += $this_opposite_SCpolarSASA;
								$avg_opposite_SCpolarSASA_count += 1;
							}
						}
						if ($avg_opposite_SCpolarSASA_count > 0) {
							$avg_opposite_SCpolarSASA = $avg_opposite_SCpolarSASA / $avg_opposite_SCpolarSASA_count;
						} else {
							$avg_opposite_SCpolarSASA = 0;
						}
						my $this_hydrophobic_side_over_opposite_side = 0;
						if ($avg_opposite_SCpolarSASA > 0) {
							$this_hydrophobic_side_over_opposite_side = $this_amphipathic_SCpolarSASA / $avg_opposite_SCpolarSASA;
						}
						if ($opposite_side_is_not_hydrophobic == 1) {
							if ($this_hydrophobic_side_over_opposite_side <= $max_hydrophobic_side_over_opposite_side) {
my $print_SASA = sprintf( "% 7.1f", $this_amphipathic_SCpolarSASA );#####
my $print_div = sprintf( "% 7.2f", $this_hydrophobic_side_over_opposite_side );#####
print "i=$i, r=" . $residue_points_of_the_amphipathic_face[$i] . ", this_amphipathic_SCpolarSASA=$print_SASA, div=$print_div <===> IS AMPHIPATHIC\n";#####
								$residue_point_is_amphipathic_helix[$i] = 1;
							} else {#####
my $print_SASA = sprintf( "% 7.1f", $this_amphipathic_SCpolarSASA );#####
my $print_div = sprintf( "% 7.1f", $this_hydrophobic_side_over_opposite_side );#####
print "i=$i, r=" . $residue_points_of_the_amphipathic_face[$i] . ", this_amphipathic_SCpolarSASA=$print_SASA, div=$print_div\n";#####
							}
						} else {#####
my $print_SASA = sprintf( "% 7.1f", $this_amphipathic_SCpolarSASA );#####
print "i=$i, r=" . $residue_points_of_the_amphipathic_face[$i] . ", this_amphipathic_SCpolarSASA=$print_SASA\n";#####
						}
					} else {#####
my $print_SASA = sprintf( "% 7.1f", $this_amphipathic_SCpolarSASA );#####
print "i=$i, r=" . $residue_points_of_the_amphipathic_face[$i] . ", this_amphipathic_SCpolarSASA=$print_SASA\n";#####
					}
				}
				
				my $in_an_amphipathic_helix = 0;
				my $start_of_this_amphipathic_helix_i = 0;
				my $start_of_this_amphipathic_helix_residue_num = 0;
				for (my $r = 0; $r < @residue_point_is_amphipathic_helix; $r++) {
					if ($residue_point_is_amphipathic_helix[$r] == 1) {
						my $this_residue_point = $residue_points_of_the_amphipathic_face[$r];
						if ($in_an_amphipathic_helix == 0) {
							$in_an_amphipathic_helix = 1;
							$start_of_this_amphipathic_helix_i = $r;
							$start_of_this_amphipathic_helix_residue_num = int($this_residue_point);
						}
					} else { # $residue_point_is_amphipathic_helix[$r] == 0
						if ($in_an_amphipathic_helix == 1) {
							$in_an_amphipathic_helix = 0;
							my $num_turns_in_this_amphipathic_helix = $r - $start_of_this_amphipathic_helix_i;
							if ($num_turns_in_this_amphipathic_helix >= $min_num_turns_for_amphipathic_helix) {
								my $this_residue_point = $residue_points_of_the_amphipathic_face[$r];
								my $end_of_this_amphipathic_helix_residue_num = int($this_residue_point);
								if ($end_of_this_amphipathic_helix_residue_num != int($this_residue_point)) {
									$end_of_this_amphipathic_helix_residue_num = int($this_residue_point) + 1;
								}
print "     !!!!! GOT A HELIX !!!!! start=$start_of_this_amphipathic_helix_residue_num, end=$end_of_this_amphipathic_helix_residue_num\n";#####
								my %this_amphipathic_helix_hash;
								my $this_amphipathic_helix = \%this_amphipathic_helix_hash;
								$this_amphipathic_helix->{'start_residue'} = $start_of_this_amphipathic_helix_residue_num;
								$this_amphipathic_helix->{'end_residue'} = $end_of_this_amphipathic_helix_residue_num;
								my $push_back_index = $#{$amphipathic_helices} + 1;
								$amphipathic_helices->[$push_back_index] = $this_amphipathic_helix;
							}
						}
					}
				}
				if ($in_an_amphipathic_helix == 1) {
					$in_an_amphipathic_helix = 0;
					my $r = @residue_point_is_amphipathic_helix - 1;
					my $num_turns_in_this_amphipathic_helix = $r - $start_of_this_amphipathic_helix_i + 1;
					if ($num_turns_in_this_amphipathic_helix >= $min_num_turns_for_amphipathic_helix) {
						my $this_residue_point = $residue_points_of_the_amphipathic_face[$r];
						my $end_of_this_amphipathic_helix_residue_num = int($this_residue_point);
						if ($end_of_this_amphipathic_helix_residue_num != int($this_residue_point)) {
							$end_of_this_amphipathic_helix_residue_num = int($this_residue_point) + 1;
						}
print "     !!!!! GOT A HELIX !!!!! start=$start_of_this_amphipathic_helix_residue_num, end=$end_of_this_amphipathic_helix_residue_num\n";#####
						my %this_amphipathic_helix_hash;
						my $this_amphipathic_helix = \%this_amphipathic_helix_hash;
						$this_amphipathic_helix->{'start_residue'} = $start_of_this_amphipathic_helix_residue_num;
						$this_amphipathic_helix->{'end_residue'} = $end_of_this_amphipathic_helix_residue_num;
						my $push_back_index = $#{$amphipathic_helices} + 1;
						$amphipathic_helices->[$push_back_index] = $this_amphipathic_helix;
					}
				}
			}
		}

		# The helices will be overlapping, 
		# because this algorithm will find the same helix more than once.

		my @final_amphipathic_helices_array;
		my $final_amphipathic_helices = \@final_amphipathic_helices_array;
		my @final_residue_is_amphipathic;
		for (my $r = 1; $r <= $residues_dot_size; $r++) {
			$final_residue_is_amphipathic[$r] = 0;
		}

		my $amphipathic_helices_dot_size = $#{$amphipathic_helices} + 1;
		for ( my $i = 0; $i < $amphipathic_helices_dot_size; $i++ ) {
			my $amphipathic_helix_start = $amphipathic_helices->[$i]->{'start_residue'};
			my $amphipathic_helix_end = $amphipathic_helices->[$i]->{'end_residue'};
			for ( my $j = $amphipathic_helix_start; $j <= $amphipathic_helix_end; $j++ ) {
				$final_residue_is_amphipathic[$j] = 1;
			}
		}

print "final_residue_is_amphipathic = " . Dumper(@final_residue_is_amphipathic) . "\n";#####
		my $in_helix = 0;
		my $start_helix = 0;
		for (my $r = 1; $r <= $residues_dot_size; $r++) {
			if ($final_residue_is_amphipathic[$r] == 1) {
				if ($in_helix == 0) {
					$in_helix = 1;
					$start_helix = $r;
				}
			} else { # $final_residue_is_amphipathic[$r] == 0
				if ($in_helix == 1) {
					$in_helix = 0;
					my $end_helix = $r - 1;
					my %start_end_hash;
					my $start_end = \%start_end_hash;
					$start_end->{'start_residue'} = $start_helix;
					$start_end->{'end_residue'} = $end_helix;
					my $push_back_index = $#{$final_amphipathic_helices} + 1;
					$final_amphipathic_helices->[$push_back_index] = $start_end;
print "START=$start_helix, END=$end_helix\n";#####
				}
			}
		}
		if ($in_helix == 1) {
			$in_helix = 0;
			my $end_helix = $residues_dot_size;
			my %start_end_hash;
			my $start_end = \%start_end_hash;
			$start_end->{'start_residue'} = $start_helix;
			$start_end->{'end_residue'} = $end_helix;
			my $push_back_index = $#{$final_amphipathic_helices} + 1;
			$final_amphipathic_helices->[$push_back_index] = $start_end;
print "START=$start_helix, END=$end_helix\n";#####
		}

		##### emma, this is just for debugging
		##### emma, have to output this in output files
		print "\n\nFINAL AMPHIPATHIC HELICES :\n\n";#####
		my $final_amphipathic_helices_dot_size = $#{$final_amphipathic_helices} + 1;
		for ( my $i = 0; $i < $final_amphipathic_helices_dot_size; $i++ ) {
			my $amphipathic_helix_start = $final_amphipathic_helices->[$i]->{'start_residue'};
			my $amphipathic_helix_end = $final_amphipathic_helices->[$i]->{'end_residue'};
			my $output_line = $amphipathic_helix_start . "  -  " . $amphipathic_helix_end;
			print "$output_line\n";
		}

		$amphipathic_helices_chains->[$c] = $final_amphipathic_helices;
	}

	return;
}



sub program_display_html_form {

	html_header();

	my $request_uri = $ENV{'REQUEST_URI'};
	my @bits = split(/\//, $request_uri);
	my $pgm_name = $bits[$#bits];
	my $got_errors = 0;
	my $err_msg = '';

print "	<form action='$pgm_name' method='post'>",
"	<table cellspacing=0 cellpadding=2 border=0>";

	if ($got_errors == 1) {
print "		<tr>\n",
"			<td valign='top' align='left' colspan='3'>",
"				<p><font color='red'>$err_msg<br><b>Didn't do a search because of input error(s).</b></font></p>",
"			</td>",
"		</tr>\n";
	}

print "		<tr>\n",
"			<td valign='top' align='left'>\n",
"				<table cellspacing=0 cellpadding=0 border=0>\n",
"					<tr>\n",
"						<td valign='top' align='left'>\n",
							"Enter an amino acid sequence in FASTA format.<br>\n",
							"Then press the GO button to predict transmembrane segments (TMS)<br>\n",
							"using the VALPRED algorithm, which uses the PROFILES3D and REPIMPS<br>\n",
							"algorithms of solvent-accessible surface area.<br><br>\n",
"						</td>\n",
"					</tr>\n",
"					<tr>\n",
"						<td valign='top' align='left'>\n",
"							<textarea name='search_text' rows='10' cols='50'></textarea>\n",
"						</td>\n",
"					</tr>\n",
"					<tr>\n",
"						<td valign='top' align='left'>\n",
							"<input type='submit' name='submit' value='&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;GO&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'>\n",
"						</td>\n",
"					</tr>\n",
"					<tr>\n",
"						<td valign='top' align='left'>\n",
							"<br>eg. Enter the FASTA sequence as :<br><br>\n",
							">1h2s_B Sensory rhodopsin II<br>",
							"AVFIFVGALTVLFGAIAYGEVTAAAATGDAAAVQEAAVSA<br>",
							"ILGLIILLGINLGLVAATL<br>",
							"<br>or as :<br><br>\n",
							"AVFIFVGALTVLFGAIAYGEVTAAAATGDAAAVQEAAVSA<br>",
							"ILGLIILLGINLGLVAATL<br>",
	
							"<br>Here is an example output graph from VALPRED :<br>\n",
							"<img src='graph_valpred1.gif' border='0' alt='output from VALPRED'><br>\n",

							"<br>Here is an example output graph from VALPRED2 :<br>\n",
							"<img src='graph_valpred2.gif' border='0' alt='output from VALPRED2'><br>\n",

							"<br>Here is another example of input :<br><br>\n",
							">KCSA_STRLI Voltage-gated potassium channel<br>",
							"MPPMLSGLLARLVKLLLGRHGSALHWRAAGAATVLLVIVL<br>",
							"LAGSYLAVLAERGAPGAQLITYPRALWWSVETATTVGYGD<br>",
							"LYPVTLWGRLVAVVVMVAGITSFGLVTAALATWFVGREQE<br>",
							"RRGHFVRHSEKAAEEAYTRTTRALHERFDRLERMLDDNRR<br>",
							"<br><br>\n",

							"<b>Notes about running this web page :</b><br>\n",
							"<br>\n",
							"If your amino acid sequence is too long for this server (around over 600 residues long), ",
							"this web page will time out and not display the result (it will display a blank page, or take a long time to process and then eventually time out.).<br><br>",
							"In that case, you can either download this Perl program ",
							"from <a href='http://www.canoz.com/valpred/perl_valpred_2d_pl.txt'>http://www.canoz.com/valpred/perl_valpred_2d_pl.txt</a><br>",
							"and run it on your own cgi-script perl-enabled web server<br>",
							"or run it in command line mode in your own Perl installation<br>",
							"(Linux Ubuntu 8.04 was used for testing).<br><br>",
							"Your own web server may allow the program to run for longer before timing out (eg. you may be able to process an 800 residue long sequence and get a result without timing out, which is why no length validation has been added to this program so that you can use the program as is on a more powerful computer.).<br><br>",
							"In command line mode, depending on your computer's memory size, you will probably be able to process longer sequences (our Linux testing machine was able to process 1800 residue long sequences).<br><br>",
							"You can also download the Windows version of this program ",
							"from <a href='http://www.canoz.com/valpred/valpred_execution_files_jun2009.zip'>http://www.canoz.com/valpred/valpred_execution_files_jun2009.zip</a> ",
							"or from <a href='http://www.canoz.com/valpred/valpred_execution_files_jun2009.rar'>http://www.canoz.com/valpred/valpred_execution_files_jun2009.rar</a> ",
							"and run it on your Windows PC (Windows XP Home Edition was used for testing).<br><br>",
							"When you need to process sequences that are too long for the machine you are using, then a valid strategy is to break up your sequence into overlapping pieces, run each separately, and visually piece them together.<br><br>",
							"The Windows version has implemented VALPRED, whereas the Perl version has implemented VALPRED and VALPRED2 implemented.<br><br>",
							"The Perl command line version allows you to specify and display the true TMS in addition to the predicted TMS.<br><br>",
							"The Perl command line version allows you to run multiple sequences in batch mode.<br><br>",
							"Incidently, the Windows version has the 'VALPRED 3D' implementation of the 3D verification of a protein structure with its REPIMPS and PROFILES3D environments, as described by Dastmalchi et al. 2001. ",
							"(The 'VALPRED 3D' implementation may need some final attention before general use.) This Perl program does not have the 'VALPRED 3D' implementation.<br><br>",
							"<br>\n",

							"<br><br><b>References :</b><br>\n",
							"<br>\n",
							"Lawrence K. Lee, Noeris K. Salam, Hong Wing Lee, Emma M. Rath and W. Bret Church (in preparation)<br>\n",
							"'Transmembrane helix analysis using threading approaches'<br>\n",
							"<br>\n",
							"James U. Bowie, Roland L&uuml;thy, David Eisenberg (1991)<br>\n",
							"'A Method to Identify Protein Sequences that Fold into a Known Three-Dimensional Structure'<br>\n",
							"Science, 253, 164-170<br>\n",
							"<br>\n",
							"Roland L&uuml;thy, James U. Bowie, David Eisenberg (1992)<br>\n",
							"'Assessment of protein models with three-dimensional profiles'<br>\n",
							"Nature, 356, 83-85<br>\n",
							"<br>\n",
							"Siavoush Dastmalchi, Michael B. Morris and W. Bret Church (2001)<br>\n",
							"'Modelling of the structural features of integral-membrane proteins using REPIMPS<br>\n",
							"(Reverse-Environment Prediction of Integral Membrane Protein Structure)'<br>\n",
							"Protein Science, 10, 1529-1538<br>\n",

							"<br><br><b>Benchmarking :</b><br>\n",
							"<br>\n",
							"VALPRED and VALPRED2 program have been benchmarked at the<br>\n",
							"<a href='http://www.canoz.com/benchmark/benchmark.pl'>Benchmark of Membrane Helix Predictions from Sequence</a> server.<br>\n",

							"<br><br><b>Source Code :</b> <a href='http://www.canoz.com/valpred/perl_valpred_2d_pl.txt'>perl_valpred_2d_pl.txt</a><br>\n",
							"<br>\n",

							"<br><br><b>COPYRIGHT :</b><br>\n",
							"Valpred Windows Program Copyright &copy; 2010 Lawrence K. Lee. All rights reserved.<br>\n",
							"Valpred Perl Program Copyright &copy; 2010 Emma M. Rath. All rights reserved.<br>\n",
							"Soon this program will be released under an open source software license such as GNU General Public License or<br>\n",
							"Creative Commons license for Free Software Foundation's GNU General Public License at creativecommons.org<br>\n",
							"<br>\n",

"						</td>\n",
"					</tr>\n",
"				</table>\n",
"			</td>\n",
"			<td width='20'></td>\n",
"			<td valign='top' align='left'>\n",
				"Please note that if your sequence is too long for this server (around over 600 residues long), this server will not give an answer and will either display a blank page or time out. ",
				"To know what can be done if this happens, please see the explanations below left.<br>\n",
				"<br>\n",
				"OPTIONS :<br><br>\n",
"				<table cellspacing=0 cellpadding=0 border=0>\n";

print "					<tr>\n",
"						<td valign='top' align='left'>SASA RESOLUTION :&nbsp;<br>(Solvent Assessible Surface Area)</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<select name='option_SASARes'><option value='0'>0</option><option value='1'>1</option><option value='2' selected='selected'>2</option><option value='3'>3</option><option value='4'>4</option><option value='5'>5</option></select>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>GRAPH RESIDUE WIDTH :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_graph_residue_width' value='" . $program_options->{'graph_residue_width'} . "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>GRAPH HEIGHT :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_graph_height' value='" . $program_options->{'graph_height'} . "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>TMS PREDICTION ALGORITHM :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<select name='option_valpred1_or_valpred2'><option value='1' selected='selected'>VALPRED</option><option value='2'>VALPRED2</option></select>\n",
"						</td>\n",
"					</tr>\n";
#print "					<tr>\n",
#"						<td valign='top' align='left'>PREDICT AMPHIPATHIC HELICES? :&nbsp;</td>\n",
#"						<td valign='top' align='left' colspan='2'>\n",
#							"<select name='option_predict_amphipathic_helices'><option value='0' selected='selected'>No</option><option value='1'>Yes</option></select>\n",
#"						</td>\n",
#"					</tr>\n";

print "					<tr><td height='20' colspan='3'></td></tr>\n";



print "					<tr>\n",
"						<td colspan='3'>OPTIONS FOR VALPRED ALGORITHM :</td>\n",
"					</tr>\n";
print "					<tr><td height='20' colspan='3'></td></tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>MINIMUM AREA DIFFERENCE :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_minAreaDiff_for_valpred1' value='" . $TMPred2D->{'minAreaDiff_for_valpred1'} . "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>MOVING AVERAGE :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_mveAve' value='" . $TMPred2D->{'mveAve'} . "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>MIN AV SASA :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_minAveSasa' value='" . $TMPred2D->{'minAveSasa'}. "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>MIN TMS LENGTH :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_TMLengthMin_for_valpred1' value='" . $TMPred2D->{'TMLengthMin_for_valpred1'}. "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>MAX TMS LENGTH :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_TMLengthMax_for_valpred1' value='" . $TMPred2D->{'TMLengthMax_for_valpred1'} . "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>MIN AV PROFILES3D RANGE :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_aveProfRangeMin' value='" . $TMPred2D->{'aveProfRangeMin'}. "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>MAX AV PROFILES3D RANGE :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_aveProfRangeMax' value='" . $TMPred2D->{'aveProfRangeMax'} . "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>MIN AV REPIMPS RANGE :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_aveRepRangeMin' value='" . $TMPred2D->{'aveRepRangeMin'} . "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>MAX AV REPIMPS RANGE :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_aveRepRangeMax' value='" . $TMPred2D->{'aveRepRangeMax'} . "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>MIN LENGTH DIVIDED BY AREA :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_lenDivAreaMin' value='" . $TMPred2D->{'lenDivAreaMin'} . "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>MAX LENGTH DIVIDED BY AREA :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_lenDivAreaMax' value='" . $TMPred2D->{'lenDivAreaMax'}. "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>MIN SCORE DIFFERENCE :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_minScoreDiff_for_valpred1' value='" . $TMPred2D->{'minScoreDiff_for_valpred1'}. "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr><td height='20' colspan='3'></td></tr>\n";



print "					<tr>\n",
"						<td colspan='3'>OPTIONS FOR VALPRED 2 ALGORITHM :</td>\n",
"					</tr>\n";
print "					<tr><td height='20' colspan='3'></td></tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>IF VALPRED2 USED,<br>SHOW ALL LINES USED BY VALPRED2 :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<select name='option_show_all_valpred2_lines'><option value='1' selected='selected'>Yes</option><option value='0'>No</option></select>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>MINIMUM AREA DIFFERENCE :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_minAreaDiff_for_valpred2' value='" . $TMPred2D->{'minAreaDiff_for_valpred2'} . "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>MV.AV.1 :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_mov_avg_1' value='" . $TMPred2D->{'mov_avg_1'} . "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>MV.AV.2 :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_mov_avg_2' value='" . $TMPred2D->{'mov_avg_2'} . "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>MV.AV.3 :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_mov_avg_3' value='" . $TMPred2D->{'mov_avg_3'}. "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>MV.AV.4 :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_mov_avg_4' value='" . $TMPred2D->{'mov_avg_4'}. "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>MIN TMS LENGTH :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_TMLengthMin_for_valpred2' value='" . $TMPred2D->{'TMLengthMin_for_valpred2'}. "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>MAX TMS LENGTH :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_TMLengthMax_for_valpred2' value='" . $TMPred2D->{'TMLengthMax_for_valpred2'} . "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>MAX LENGTH NON-TMS :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_max_nonTMS_length' value='" . $TMPred2D->{'max_nonTMS_length'}. "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>MIN LENGTH NON-TMS :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_min_nonTMS_length' value='" . $TMPred2D->{'min_nonTMS_length'}. "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>MIN NON-TMS SCORE DIFFERENCE :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_min_nonTMS_score_diff' value='" . $TMPred2D->{'min_nonTMS_score_diff'} . "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>MIN NON-TMS AREA :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_min_nonTMS_area' value='" . $TMPred2D->{'min_nonTMS_area'}. "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>TWILIGHT TMS ZONE AREA PER RESIDUE - LOWER LIMIT :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_twilight_area_per_residue_lower_limit' value='" . $TMPred2D->{'twilight_area_per_residue_lower_limit'}. "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>TWILIGHT TMS ZONE AREA PER RESIDUE - UPPER LIMIT :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_twilight_area_per_residue_upper_limit' value='" . $TMPred2D->{'twilight_area_per_residue_upper_limit'}. "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>IGNORE TWILIGHT TMS ZONE :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<select name='option_ignore_twilight_area'><option value='1'>Yes</option><option value='0' selected='selected'>No</option></select>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>NO. RESIDUES SEARCH AREA FOR TMS NEIGHBOURS :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_search_area_for_neighbour_tms' value='" . $TMPred2D->{'search_area_for_neighbour_tms'}. "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>VALPRED2 MIN HELIX LENGTH :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_helix_length_min' value='" . $TMPred2D->{'helix_length_min'}. "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>MIN SCORE DIFF FOR TMS ENDS (MOV.AVG.1) :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_min_score_diff_for_TMS_ends' value='" . $TMPred2D->{'min_score_diff_for_TMS_ends'}. "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>MIN SCORE DIFF. MOV.AVG.1 :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_minScoreDiff_movavg1_for_valpred2' value='" . $TMPred2D->{'minScoreDiff_movavg1_for_valpred2'}. "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>MIN SCORE DIFF. MOV.AVG.2 :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_minScoreDiff_movavg2_for_valpred2' value='" . $TMPred2D->{'minScoreDiff_movavg2_for_valpred2'}. "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>MIN AV AREA MOV.AVG.1 :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_min_avg_area_movavg1' value='" . $TMPred2D->{'min_avg_area_movavg1'}. "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr>\n",
"						<td valign='top' align='left'>MIN AV AREA MOV.AVG.2 :&nbsp;</td>\n",
"						<td valign='top' align='left' colspan='2'>\n",
							"<input type='text' name='option_min_avg_area_movavg2' value='" . $TMPred2D->{'min_avg_area_movavg2'}. "'>\n",
"						</td>\n",
"					</tr>\n";
print "					<tr><td height='20' colspan='3'></td></tr>\n";



print "					<tr>\n",
"						<td colspan='3'>DISPLAY OPTIONS :</td>\n",
"					</tr>\n";
print "					<tr><td height='20' colspan='3'></td></tr>\n";

					my @option_lines = ('PROFILES3D-LINE','REPIMPS-LINE','BOTH-LINE','PROFILES3D-MOV-AVG-1','REPIMPS-MOV-AVG-1','TMS-MOV-AVG-1',
						'PROFILES3D-MOV-AVG-2','REPIMPS-MOV-AVG-2','TMS-MOV-AVG-2','PROFILES3D-MOV-AVG-3','REPIMPS-MOV-AVG-3','HYDROPHILIC-MOV-AVG-3',
						'PROFILES3D-MOV-AVG-4','REPIMPS-MOV-AVG-4','HYDROPHILIC-MOV-AVG-4',
						'TMS');
						# 'TMS','TRUE-TMS','TRUE-INSIDE','TRUE-OUTSIDE');
					my @option_attributes = ('showGraph','colour_name','legend_name','smooth');
					my @option_attributes_display_name = ('Show line','Line colour','Line label','Smooth line');

					for ( my $o = 0; $o < @option_lines; $o++ ) {

						my $option_name = $option_lines[$o];
						my $ln = $program_options->{'lnOpt_index'}->{$option_name};
						my $html_display_option_name = $program_options->{'lnOpt'}->[$ln]->{'legend_name'} . ' :&nbsp;';

						for ( my $n = 0; $n < @option_attributes; $n++ ) {

							my $option_attribute = $option_attributes[$n];
							my $param_name = 'option_' . $ln . '_' . $option_attribute;

							my $html_display_attribute_name = $option_attributes_display_name[$n];

							my $html_input = '';

							if ($option_attribute eq 'showGraph') {
								my $html_display_select_box = "<select name='$param_name'>";
								my @html_select_options_display_name = ('Yes','No');
								my @html_select_options_value = (1,0);
								for (my $v = 0; $v < @html_select_options_value; $v++) {
									my $html_selected_text = "";
									if ($program_options->{'lnOpt'}->[$ln]->{'showGraph'} == $html_select_options_value[$v]) {
										$html_selected_text = " selected='selected'";
									}
									$html_display_select_box .= "<option value='" . $html_select_options_value[$v] . "'" . $html_selected_text . ">" . $html_select_options_display_name[$v] . "</option>";
								}
								$html_display_select_box .= '</select>';
								$html_input = $html_display_select_box;
							} elsif ($option_attribute eq 'mov_avg') { # 'Moving average'
								my $html_input_value = $program_options->{'lnOpt'}->[$ln]->{'mov_avg'};
								$html_input = "<input type='text' name='$param_name' value='$html_input_value'>";
							} elsif ($option_attribute eq 'colour_name') {
								my $html_display_select_box = "<select name='$param_name'>";
								my @html_select_options_display_name = ('blue','red','black','green','purple','cyan','orange','yellow','pink','lred','lgreen','lblue','lyellow','lpurple','lorange');
								my @html_select_options_value = ('blue','red','black','green','purple','cyan','orange','yellow','pink','lred','lgreen','lblue','lyellow','lpurple','lorange');
								for (my $v = 0; $v < @html_select_options_value; $v++) {
									my $html_selected_text = "";
									if ($program_options->{'lnOpt'}->[$ln]->{'colour_name'} eq $html_select_options_value[$v]) {
										$html_selected_text = " selected='selected'";
									}
									$html_display_select_box .= "<option value='" . $html_select_options_value[$v] . "'" . $html_selected_text . ">" . $html_select_options_display_name[$v] . "</option>";
								}
								$html_display_select_box .= '</select>';
								$html_input = $html_display_select_box;
							} elsif ($option_attribute eq 'legend_name') {
								my $html_input_value = $program_options->{'lnOpt'}->[$ln]->{'legend_name'};
								$html_input = "<input type='text' name='$param_name' value='$html_input_value'>";
							} elsif ($option_attribute eq 'smooth') {
								my $html_display_select_box = "<select name='$param_name'>";
								my @html_select_options_display_name = ('Yes','No');
								my @html_select_options_value = (1,0);
								for (my $v = 0; $v < @html_select_options_value; $v++) {
									my $html_selected_text = "";
									if ($program_options->{'lnOpt'}->[$ln]->{'smooth'} == $html_select_options_value[$v]) {
										$html_selected_text = " selected='selected'";
									}
									$html_display_select_box .= "<option value='" . $html_select_options_value[$v] . "'" . $html_selected_text . ">" . $html_select_options_display_name[$v] . "</option>";
								}
								$html_display_select_box .= '</select>';
								$html_input = $html_display_select_box;
							}

		print "					<tr>\n",
		"						<td valign='top' align='left'><nobr>$html_display_option_name</nobr></td>\n",
		"						<td valign='top' align='left'><nobr>$html_display_attribute_name:</nobr></td>\n",
		"						<td valign='top' align='left'>\n",
									"$html_input\n",
		"						</td>\n",
		"					</tr>\n";

							$html_display_option_name = '';
						}
					}

print "				</table>\n",
				"<br><br><b>Notes about these parameters :</b><br>\n",
				"<br>\n",
				"VALPRED uses the area between the REPIMPS and PROFILES3D lines to predict transmembrane segments (TMS). ",
				"VALPRED performs well at predicting transmembrane helices whilst avoiding false positives. ",
				"VALPRED2 performs well at predicting predicting small helices, such as in potassium channels, but tends to predict some false positives. <br>\n",
				"<br>\n",
				"By default, VALPRED2 will question whether area falling within the TWILIGHT TMS ZONE AREA PER RESIDUE - LOWER LIMIT and TWILIGHT TMS ZONE AREA PER RESIDUE - UPPER LIMIT might be a hydrophobic helix in a soluble protein instead of being a membrane helix. If there are no other membrane helices close by, within NO. RESIDUES SEARCH AREA FOR TMS NEIGHBOURS, that have areas above the twilight zone, then it will not be marked as a membrane helix. To display all such hydrophobic helices in soluble proteins, set IGNORE TWILIGHT TMS ZONE to Yes.<br>\n",
				"<br>\n",

"			</td>\n",
"		</tr>\n",
"	</table>\n",
"	</form>\n";

	print "</td></tr></table>\n";
	print $q->end_html();
}



sub program_get_input {

	set_program_options_defaults();
	# read_program_options_defaults_file(); # defaults file from windows version
	adjust_options();

	$program_mode = 'command-line';
	if (defined($ENV{'REQUEST_METHOD'})) {
		if ($ENV{'REQUEST_METHOD'} eq 'POST') {
			my $params = $q->Vars;
			if (defined($params->{'chart_coords'})) {
				$program_mode = 'cgi-output-graph';
			} else {
				$program_mode = 'cgi-output';
			}
		} else {
			$program_mode = 'cgi-display-form';
		}
	} else {
		GetOptions( "infile=s" => \$input_infile, "fastachainfile=s" => \$input_fasta_chain_file, "fastafile=s" => \$input_fasta_file, "inscoresfile=s" => \$input_inscoresfile, 
				"predictions=s" => \$input_predictions_file, "benchmarkpredictions=s" => \$input_benchmarkpredictions_file, "graphpointsfile=s" => \$input_graph_points_file, 
				"version=s" => \$input_flag_version, "amphipathic" => \$input_flag_amphipathic, 
				"options=s" => \$input_options_file, "tms=s" => \$input_tms_file, "multiplegraphs" => \$input_flag_multiplegraphs, "benchmark" => \$input_flag_benchmark, 
				"help" => \$input_flag_help, "h" => \$input_flag_h, "debug" => \$input_debug_flag );
		if ((defined ($input_infile)) || (defined ($input_fasta_chain_file)) || (defined ($input_fasta_file)) || (defined ($input_graph_points_file))) {
			$program_mode = 'command-line';
		}
	}
}



sub program_output_command_line {

	my $got_command_line_error = 0;

	my $command_line_help_text = '';
	$command_line_help_text .= "\n";
	$command_line_help_text .= "This program is called VALPRED, and it predicts transmembrane segments (TMS) from amino acid sequences.\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "There are various ways to call this program in command line mode :\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "          perl perl_valpred_2d.pl -infile YOUR_FILE_NAME\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "or\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "          perl perl_valpred_2d.pl -fastachainfile YOUR_FILE_NAME\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "or\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "          perl perl_valpred_2d.pl -fastafile YOUR_FILE_NAME\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "or\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "          perl perl_valpred_2d.pl -inscoresfile YOUR_FILE_NAME\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "or\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "          perl perl_valpred_2d.pl -graphpointsfile YOUR_FILE_NAME\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "If the -infile option is used, then your file needs to contain one or more sequences, one line per sequence, in the following format :\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "          SEQ-ID:::::>FASTA-ID:::::SEQUENCE OF ONE-CODE AMINO ACIDS\n";
	$command_line_help_text .= "          SEQ-ID:::::>FASTA-ID:::::SEQUENCE OF ONE-CODE AMINO ACIDS\n";
	$command_line_help_text .= "          SEQ-ID:::::>FASTA-ID:::::SEQUENCE OF ONE-CODE AMINO ACIDS\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "If the -fastachainfile option is used, then your file needs to contain one sequence of one or more chains, over multiple lines, in the following fasta format :\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "          >FASTA-ID WITH OPTIONAL CHAIN-ID\n";
	$command_line_help_text .= "          SEQUENCE OF ONE-CODE AMINO ACIDS\n";
	$command_line_help_text .= "          CONTINUATION OF SEQUENCE OF ONE-CODE AMINO ACIDS\n";
	$command_line_help_text .= "          CONTINUATION OF SEQUENCE OF ONE-CODE AMINO ACIDS ETC.\n";
	$command_line_help_text .= "          >FASTA-ID WITH OPTIONAL CHAIN-ID\n";
	$command_line_help_text .= "          SEQUENCE OF ONE-CODE AMINO ACIDS\n";
	$command_line_help_text .= "          CONTINUATION OF SEQUENCE OF ONE-CODE AMINO ACIDS\n";
	$command_line_help_text .= "          CONTINUATION OF SEQUENCE OF ONE-CODE AMINO ACIDS ETC.\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "If the -fastafile option is used, then your file needs to contain one or more sequences, over multiple lines, in the following fasta format :\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "          >FASTA-ID\n";
	$command_line_help_text .= "          SEQUENCE OF ONE-CODE AMINO ACIDS\n";
	$command_line_help_text .= "          CONTINUATION OF SEQUENCE OF ONE-CODE AMINO ACIDS\n";
	$command_line_help_text .= "          CONTINUATION OF SEQUENCE OF ONE-CODE AMINO ACIDS ETC.\n";
	$command_line_help_text .= "          >FASTA-ID\n";
	$command_line_help_text .= "          SEQUENCE OF ONE-CODE AMINO ACIDS\n";
	$command_line_help_text .= "          CONTINUATION OF SEQUENCE OF ONE-CODE AMINO ACIDS\n";
	$command_line_help_text .= "          CONTINUATION OF SEQUENCE OF ONE-CODE AMINO ACIDS ETC.\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "If the -graphpointsfile option is used, then your file needs to contain one or more sequences, in the valpred graph-points output format.\n";
	$command_line_help_text .= "You need to use your valpred_graphpoints.txt output file from a previous run, as input here.\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "You can also add the -version option and/or the -options option to the command-line, as shown below.\n";
	$command_line_help_text .= "          perl perl_valpred_2d.pl -infile YOUR_FILE_NAME -version valpred1\n";
	$command_line_help_text .= "          perl perl_valpred_2d.pl -infile YOUR_FILE_NAME -version valpred2\n";
	$command_line_help_text .= "          perl perl_valpred_2d.pl -infile YOUR_FILE_NAME -options YOUR_OPTIONS_FILE_NAME\n";
	$command_line_help_text .= "The -version option will set the TMS-prediction algorithm to that specified, and use default options for that specified algorithm.\n";
	$command_line_help_text .= "The -options option file allows you to set options specified in the file, instead of using defaults. (Please see program code for format.)\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "You can also add the -multiplegraphs option to the command-line, as shown below, to produce 1 graph file for every input sequence.\n";
	$command_line_help_text .= "          perl perl_valpred_2d.pl -infile YOUR_FILE_NAME -multiplegraphs\n";
	$command_line_help_text .= "          perl perl_valpred_2d.pl -fastachainfile YOUR_FILE_NAME -multiplegraphs\n";
	$command_line_help_text .= "          perl perl_valpred_2d.pl -fastafile YOUR_FILE_NAME -multiplegraphs\n";
	$command_line_help_text .= "          perl perl_valpred_2d.pl -graphpointsfile YOUR_FILE_NAME -multiplegraphs\n";
	$command_line_help_text .= "Otherwise, a graph file will be produced for only the first sequence in the input file.\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "You can also add the -tms option to the command-line, as shown below.\n";
	$command_line_help_text .= "          perl perl_valpred_2d.pl -infile YOUR_FILE_NAME -tms YOUR_TMS_FILE_NAME\n";
	$command_line_help_text .= "The TMS file shows where TMS are thought to be, in the following format :\n";
	$command_line_help_text .= "          SEQ_16:::::>2MPR:A:::::HELIX;;;;;113-115,237-239,:::::TMS;;;;;113-115,237-239,:::::\n";
	$command_line_help_text .= "          SEQ_17:::::>1MAL:A:::::HELIX;;;;;113-115,160-162,398-400,:::::TMS;;;;;113-115,160-162,398-400,:::::\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "You can also add the -benchmark option to the command-line, as shown below, to produce the YOUR_FILE_NAME.valpred_benchmark.txt file that can be used as input to the Rost/Kernytsky TMH benchmark server.\n";
	$command_line_help_text .= "          perl perl_valpred_2d.pl -fastafile YOUR_FILE_NAME -benchmark\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "If the -predictions or -benchmarkpredictions option is used with the input file being the output file of a previous valpred run, ";
	$command_line_help_text .= "and with the -tms option and input file containing observed/known tms, then the program will calculate statistics for the accuracy of the predictions.\n";
	$command_line_help_text .= "          perl perl_valpred_2d.pl -predictions YOUR_FILE.valpred_segments.txt -tms YOUR_TMS_FILE\n";
	$command_line_help_text .= "          perl perl_valpred_2d.pl -benchmarkpredictions YOUR_FILE.valpred_benchmark.txt -tms YOUR_TMS_FILE\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "If the -debug option is present on the command-line, as shown below, then the valpred.pdb file will be produced, containing the helix built from sequence.\n";
	$command_line_help_text .= "          perl perl_valpred_2d.pl -fastafile YOUR_FILE_NAME -debug\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "The various output files from this program, depending on input options, are :\n";	$command_line_help_text .= "They will be plotted on the output graph, and compared to this program's predictions.\n";
	$command_line_help_text .= "          - *.valpred_summary.txt - contains the valpred TMS summary for the first input sequence\n";
	$command_line_help_text .= "          - *.valpred_graph.gif - contains the valpred TMS graph for the first sequence\n";
	$command_line_help_text .= "          - *.valpred_calculations.txt - contains the valpred calculation results of accessibility and hydropathy environment for each residue\n";
	$command_line_help_text .= "          - *.valpred_segments.txt - is the valpred TMS output file with one line per input sequence\n";
	$command_line_help_text .= "          - *.valpred_scores.txt - is the valpred scores output file with one line per input sequence, with the PROFILES3D and REPIMPS scores for each residue\n";
	$command_line_help_text .= "          - *.valpred_graphpoints.txt - is the valpred TMS segments and graph points output file with one line per input sequence\n";
	$command_line_help_text .= "          - *.valpred_graphs.html - is an html page containing multiple valpred graphs, 1 graph per input sequence\n";
	$command_line_help_text .= "          - *.valpred_graph_*.gif - are the graph output files, one for each input sequence\n";
	$command_line_help_text .= "          - *.valpred_evaluate.txt - compares the valpred predicted to actual known TMS found in the input TMS file\n";
	$command_line_help_text .= "          - *.valpred_evaluate_summary.txt - contains a summary of the comparison of the valpred predicted to actual known TMS found in the input TMS file\n";
	$command_line_help_text .= "          - *.valpred_benchmark.txt - is the output file of TMS predictions for the Rost TMH benchmark server at http://cubic.bioc.columbia.edu/services/tmh_benchmark/\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "VALPRED and VALPRED2 use the area between the REPIMPS and PROFILES3D lines to predict transmembrane segments (TMS). ";
	$command_line_help_text .= "VALPRED performs well at predicting transmembrane helices whilst avoiding false positives. \n";
	$command_line_help_text .= "VALPRED2 performs well at predicting predicting small helices, such as in potassium channels, but tends to predict some false positives. \n";
	#$command_line_help_text .= "VALPRED1 is the default and performs best on the Kernytsky/Rost TMS benchmark server (at http://cubic.bioc.columbia.edu/services/tmh_benchmark/). ";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "The parameters used by VALPRED are : MINIMUM AREA DISTANCE, MOVING AVERAGE, MIN TMS LENGTH (VALPRED), MAX TMS LENGTH (VALPRED), MIN AV PROFILES3D RANGE, MAX AV PROFILES3D RANGE, MIN AV REPIMPS RANGE, MAX AV REPIMPS RANGE, MIN LENGTH DIVIDED BY AREA, MAX LENGTH DIVIDED BY AREA, MIN SCORE DIFFERENCE, MIN AV AREA, Profiles 3D (mv.av.1), REPIMPS (mv.av.1), Transmembrane Segments.\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "The parameters used by VALPRED2 are : MV.AV.1, MV.AV.2, MV.AV.3, MV.AV.4, MIN TMS LENGTH (VALPRED2), MAX TMS LENGTH (VALPRED2), MAX LENGTH NON-TMS, MIN LENGTH NON-TMS, MIN NON-TMS SCORE DIFFERENCE, MIN NON-TMS AREA, MIN AREA PER RESIDUE, MAX AREA PER RESIDUE, NO. RESIDUES SEARCH AREA FOR TMS NEIGHBOURS, VALPRED2 MIN HELIX LENGTH, Profiles 3D, REPIMPS, Profiles 3D (mv.av.1), REPIMPS (mv.av.1), TMS (mv.av.1), Profiles 3D (mv.av.2), REPIMPS (mv.av.2), TMS (mv.av.2), Profiles 3D (mv.av.3), REPIMPS (mv.av.3), Hydrophilic area (mv.av.3), Profiles 3D (mv.av.4), REPIMPS (mv.av.4), Hydrophilic area (mv.av.4), Transmembrane Segments.\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "You can also add the -amphipathic option to the command-line, as shown below, to predict amphipathic helices from sequence.\n";
	$command_line_help_text .= "          perl perl_valpred_2d.pl -fastafile YOUR_FILE_NAME -amphipathic\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "Please note that it is also possible to run this program as a cgi-script on a web server.\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "There is a Windows version of this program. ";
	$command_line_help_text .= "you can download the Windows version of this program from http://www.canoz.com/valpred/valpred_execution_files_jun2009.zip ";
	$command_line_help_text .= "or from http://www.canoz.com/valpred/valpred_execution_files_jun2009.rar and run it on your Windows PC.\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "The Windows version has only VALPRED1 implemented, whereas the Perl version has both VALPRED1 and VALPRED2 implemented, and has the amphipathic helices prediction implemented in it too.\n";
	$command_line_help_text .= "The Perl command line version allows you to specify and display the true TMS in addition to the predicted TMS.\n";
	$command_line_help_text .= "The Perl command line version allows you to run multiple sequences in batch mode.\n";
	$command_line_help_text .= "The Windows version, however, has the implementation of the 3D verification of a protein structure with its REPIMPS and PROFILES3D environments, as described by Dastmalchi et al. 2001. (This implementation may need some final attention before general use.) This Perl verions does not have that implementation.\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "References :\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "          Lawrence K. Lee, Noeris K. Salam, Hong Wing Lee, Emma M. Rath and W. Bret Church (in preparation)\n";
	$command_line_help_text .= "          'Transmembrane helix analysis using threading approaches'\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "          James U. Bowie, Roland Lüthy, David Eisenberg (1991)\n";
	$command_line_help_text .= "          'A Method to Identify Protein Sequences that Fold into a Known Three-Dimensional Structure'\n";
	$command_line_help_text .= "          Science, 253, 164-170\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "          Roland Lüthy, James U. Bowie, David Eisenberg (1992)\n";
	$command_line_help_text .= "          'Assessment of protein models with three-dimensional profiles'\n";
	$command_line_help_text .= "          Nature, 356, 83-85\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "          Siavoush Dastmalchi, Michael B. Morris and W. Bret Church (2001)\n";
	$command_line_help_text .= "          'Modelling of the structural features of integral-membrane proteins using REPIMPS\n";
	$command_line_help_text .= "          (Reverse-Environment Prediction of Integral Membrane Protein Structure)'\n";
	$command_line_help_text .= "          Protein Science, 10, 1529-1538\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "Copyright :\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "          Valpred Windows Program Copyright © 2010 Lawrence K. Lee. All rights reserved.\n";
	$command_line_help_text .= "          Valpred Perl Program Copyright © 2010 Emma M. Rath. All rights reserved.\n";
	$command_line_help_text .= "          Soon this program will be released under an open source software license such as GNU General Public License or \n";
	$command_line_help_text .= "          Creative Commons license for Free Software Foundation's GNU General Public License at creativecommons.org\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "\n";
	$command_line_help_text .= "\n";

	my $num_input_files = 0;
	if (defined ($input_infile)) {
		$num_input_files++;
	}
	if (defined ($input_fasta_chain_file)) {
		$num_input_files++;
	}
	if (defined ($input_fasta_file)) {
		$num_input_files++;
	}
	if (defined ($input_predictions_file)) {
		$num_input_files++;
	}
	if (defined ($input_benchmarkpredictions_file)) {
		$num_input_files++;
	}
	if (defined ($input_graph_points_file)) {
		$num_input_files++;
	}
	if (defined ($input_inscoresfile)) {
		$num_input_files++;
	}

	if ( $num_input_files > 1 ) {
		print $command_line_help_text;
		print "Program was not run due to input errors. More than one input file was specied in command line mode. Don't know which input file to use.\n\n";
		$got_command_line_error = 1;
	}

	if ( $num_input_files == 0 ) {
		print $command_line_help_text;
		print "Program was not run due to input errors. No input file was specied in command line mode.\n\n";
		$got_command_line_error = 1;
	}

	if ($got_command_line_error == 0) {

		if (defined ($input_infile)) {
			$input_file_name_for_output_files = $input_infile;
		} elsif (defined ($input_fasta_file)) {
			$input_file_name_for_output_files = $input_fasta_file;
		} elsif (defined ($input_fasta_chain_file)) {
			$input_file_name_for_output_files = $input_fasta_chain_file;
		} elsif (defined ($input_predictions_file)) {
			$input_file_name_for_output_files = $input_predictions_file;
		} elsif (defined ($input_benchmarkpredictions_file)) {
			$input_file_name_for_output_files = $input_benchmarkpredictions_file;
		} elsif (defined ($input_graph_points_file)) {
			$input_file_name_for_output_files = $input_graph_points_file;
		} elsif (defined ($input_inscoresfile)) {
			$input_file_name_for_output_files = $input_inscoresfile;
		}

		$input_file_name_for_output_files = untaint($input_file_name_for_output_files);

		if (defined ($input_flag_benchmark)) {
			# open output file
			$benchmark_output_file = "$input_file_name_for_output_files.valpred_benchmark.txt";
			open( BENCHOUTFILE, ">$benchmark_output_file") or
				die "Cannot open $benchmark_output_file for writing : $!\n";
		}

		if (defined ($input_flag_multiplegraphs)) {
			$html_graphs_output_file = "$input_file_name_for_output_files.valpred_graphs.html";
		}

		if (defined ($input_tms_file)) {
			# open output file
			$evaluate_output_file = "$input_file_name_for_output_files.valpred_evaluate.txt";
			open( EVALOUTFILE, ">$evaluate_output_file") or
				die "Cannot open $evaluate_output_file for writing : $!\n";
		}

		if (defined ($input_infile)) {

			program_output_command_line_for_infile();

		} elsif (defined ($input_fasta_file)) {

			program_output_command_line_for_fastafile();

		} elsif (defined ($input_fasta_chain_file)) {

			program_output_command_line_for_fastachainfile();

		} elsif (defined ($input_predictions_file)) {

			program_output_command_line_for_predictions();

		} elsif (defined ($input_benchmarkpredictions_file)) {

			program_output_command_line_for_benchmarkpredictions();

		} elsif (defined ($input_graph_points_file)) {

			program_output_command_line_for_graphpointsfile();

		} elsif (defined ($input_inscoresfile)) {

			program_output_command_line_for_inscoresfile();

		} else {

			print "There was no input file specified in command-line mode.\n";

		}
	}

	if (defined ($input_flag_benchmark)) {
		close BENCHOUTFILE;
	}
	if (defined ($input_tms_file)) {
		close EVALOUTFILE;
		output_evaluate_summary_file();
	}
}



sub program_output_command_line_for_benchmarkpredictions {

	get_command_line_program_options();

	# open input file
	open PREDFILE, $input_benchmarkpredictions_file or die $!;
	my @input_lines = <PREDFILE>;

	my $seq_num = 0;

	my @new_input_lines;
	my $new_input_line = '';
	my $this_new_line_segments = '';
	my $finish_prev_seq = 0;

	foreach my $input_line (@input_lines) {

		chomp( $input_line );
		if ($input_line ne '') {
			$input_line = trim($input_line);
		}
		if ($input_line ne '') {

			if (substr($input_line,0,1) eq '>') {
				if ($finish_prev_seq == 1) {
					if ($this_new_line_segments eq '') {
						$this_new_line_segments = '-';
					}
					$new_input_line .= $this_new_line_segments . ':::::';
					push( @new_input_lines, $new_input_line );
				}
				$input_fasta_id = trim($input_line);
				$new_input_line = $input_fasta_id . ':::::' . $input_fasta_id . ':::::-:::::';
				$this_new_line_segments = '';
				$finish_prev_seq = 1;
			} else {
				my @bits2 = split(/\,/, $input_line);
				my $start = $bits2[0];
				my $end = $bits2[1];
				$this_new_line_segments .= $start . '-' . $end . ',';
			}
		}
	}
	if ($finish_prev_seq == 1) {
		if ($this_new_line_segments eq '') {
			$this_new_line_segments = '-';
		}
		$new_input_line .= $this_new_line_segments . ':::::';
		push( @new_input_lines, $new_input_line );
	}

	@input_lines = @new_input_lines;

	foreach my $input_line (@input_lines) {

		chomp( $input_line );
		if ($input_line ne '') {
			$input_line = trim($input_line);
		}
		if ($input_line ne '') {

			my @bits = split(/:::::/, $input_line);

			my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
			$year += 1900;
			$mon += 1;
			$mon = sprintf("%02d", $mon);
			$mday = sprintf("%02d", $mday);
			$hour = sprintf("%02d", $hour);
			$min = sprintf("%02d", $min);
			$sec = sprintf("%02d", $sec);
			my $time_str = "$year-$mon-$mday.$hour:$min:$sec";

			$seq_num++;
			$input_fasta_id = trim($bits[1]);
			$input_sequence_id = trim($bits[0]);
			#if ((defined($input_flag_multiplegraphs)) || (defined($input_flag_benchmark))) {
			print "$time_str : Processing sequence number $seq_num : $input_fasta_id\n";
			#}

			my @aa_seq_chains;
			$aa_seq_chains[0] = '-';
			$input_aa_sequence = '-';
			my $segments = trim($bits[3]);

			get_2nd_struc_info();

			if (defined($input_tms_file)) {
				my @predicted_tms_array_start;
				my @predicted_tms_array_end;
				if ($segments ne '-') {
					my @bits2 = split(/\,/, $segments);
					foreach my $segment (@bits2) {
						my @bits3 = split(/\-/, $segment);
						my $start = $bits3[0];
						my $end = $bits3[1];
						push( @predicted_tms_array_start, $start );
						push( @predicted_tms_array_end, $end );
					}
				}
				my $predicted_tms_array_start_ref = \@predicted_tms_array_start;
				my $predicted_tms_array_end_ref = \@predicted_tms_array_end;
				output_evaluate_file( $predicted_tms_array_start_ref, $predicted_tms_array_end_ref );
			}
		}
	}
}



sub program_output_command_line_for_fastachainfile {

	# open output file
	$segment_output_file = "$input_file_name_for_output_files.valpred_segments.txt";
	$scores_output_file = "$input_file_name_for_output_files.valpred_scores.txt";
	$amphipathic_output_file = "$input_file_name_for_output_files.valpred_amphipathic.txt";

	open( SEGOUTFILE, ">$segment_output_file") or
		die "Cannot open $segment_output_file for writing : $!\n";
	open( SCOREOUTFILE, ">$scores_output_file") or
		die "Cannot open $scores_output_file for writing : $!\n";
	if (defined($input_flag_amphipathic)) {
		open( AMPHIFILE, ">$amphipathic_output_file") or
			die "Cannot open $amphipathic_output_file for writing : $!\n";
	}

	if (defined($input_flag_multiplegraphs)) {

		my $segment_graphpoints_output_file = "$input_file_name_for_output_files.valpred_graphpoints.txt";
		open( SEGPTSOUTFILE, ">$segment_graphpoints_output_file") or
			die "Cannot open $segment_graphpoints_output_file for writing : $!\n";

		my $html_graphs_output_file = "$input_file_name_for_output_files.valpred_graphs.html";
		open( HTMLGRAPHSOUTFILE, ">$html_graphs_output_file") or
			die "Cannot open $html_graphs_output_file for writing : $!\n";
		create_output_graphs_html_file( 1 );
	}

	read_amino_acid_pdbs();
	get_command_line_program_options();

	# open input file
	open FASTACHAINFILE, $input_fasta_chain_file or die $!;
	my @fasta_lines = <FASTACHAINFILE>;

	my $converted_input_line = '';
	my $need_to_push_line = 0;

	my @converted_input_array;

	my $i = 0;
	foreach my $fasta_line (@fasta_lines) {

		chomp($fasta_line);
		if ($fasta_line ne '') {
			$fasta_line = trim($fasta_line);
		}

		if (substr($fasta_line,0,1) eq '>') {

			$i++;

			if ($need_to_push_line == 1) {
				push ( @converted_input_array, $converted_input_line );
				$need_to_push_line = 0;
			}

			my @bits = split(/\|/, $fasta_line);
			my $seqid = "SEQ_$i";
			if (defined $bits[3]) {
				my @bits2 = split(/\./, $bits[3]);
				$seqid = uc $bits2[0];
				$seqid =~ s/\W//g;
				$seqid =~ s/_//g;
				$seqid = "SEQ_$seqid";
			}
			$converted_input_line = $seqid . ':::::' . $fasta_line . ':::::';
			$need_to_push_line = 1;

		} elsif ($fasta_line eq '') {

			if ($need_to_push_line == 1) {
				push ( @converted_input_array, $converted_input_line );
				$need_to_push_line = 0;
			}

		} else {

			$fasta_line =~ s/\s//g;
			$converted_input_line .= $fasta_line;
		}
	}

	if ($need_to_push_line == 1) {
		push ( @converted_input_array, $converted_input_line );
		$need_to_push_line = 0;
	}

	my @aa_seq_chains;
	foreach my $converted_input_line (@converted_input_array) { # while(!(f.eof()) ){

		chomp( $converted_input_line );
		if ($converted_input_line ne '') {

			my @bits = split(/:::::/, $converted_input_line);
			if (($input_fasta_id eq '') or ($input_fasta_id eq '-')) {
				$input_fasta_id = trim($bits[1]);
			}
			if (($input_sequence_id eq '') or ($input_sequence_id eq '-')) {
				$input_sequence_id = trim($bits[0]);
			}
			my $chain_upto = @aa_seq_chains;
			$aa_seq_chains[$chain_upto] = trim($bits[2]);
			$input_aa_sequence .= trim($bits[2]);
		}
	}

	get_2nd_struc_info();

	my $ref_aa_seq_chains = \@aa_seq_chains;
	my $seq_num = 1;
	CValpredDoc_OnToolsTmpred2d( $ref_aa_seq_chains, $seq_num );

	close SEGOUTFILE;
	close SCOREOUTFILE;
	if (defined($input_flag_multiplegraphs)) {
		close SEGPTSOUTFILE;
		create_output_graphs_html_file( 3 );
		close HTMLGRAPHSOUTFILE;
	}
}



sub program_output_command_line_for_fastafile {

	# open output file
	$segment_output_file = "$input_file_name_for_output_files.valpred_segments.txt";
	$scores_output_file = "$input_file_name_for_output_files.valpred_scores.txt";
	$amphipathic_output_file = "$input_file_name_for_output_files.valpred_amphipathic.txt";

	open( SEGOUTFILE, ">$segment_output_file") or
		die "Cannot open $segment_output_file for writing : $!\n";
	open( SCOREOUTFILE, ">$scores_output_file") or
		die "Cannot open $scores_output_file for writing : $!\n";
	if (defined($input_flag_amphipathic)) {
		open( AMPHIFILE, ">$amphipathic_output_file") or
			die "Cannot open $amphipathic_output_file for writing : $!\n";
	}

	if (defined($input_flag_multiplegraphs)) {

		my $segment_graphpoints_output_file = "$input_file_name_for_output_files.valpred_graphpoints.txt";
		open( SEGPTSOUTFILE, ">$segment_graphpoints_output_file") or
			die "Cannot open $segment_graphpoints_output_file for writing : $!\n";

		my $html_graphs_output_file = "$input_file_name_for_output_files.valpred_graphs.html";
		open( HTMLGRAPHSOUTFILE, ">$html_graphs_output_file") or
			die "Cannot open $html_graphs_output_file for writing : $!\n";
		create_output_graphs_html_file( 1 );
	}

	read_amino_acid_pdbs();
	get_command_line_program_options();

	# open input file
	open FASTAFILE, $input_fasta_file or die $!;
	my @fasta_lines = <FASTAFILE>;

	my $converted_input_line = '';
	my $need_to_push_line = 0;

	my @converted_input_array;

	my $i = 0;
	foreach my $fasta_line (@fasta_lines) {

		chomp($fasta_line);
		if ($fasta_line ne '') {
			$fasta_line = trim($fasta_line);
		}

		if (substr($fasta_line,0,1) eq '>') {

			$i++;

			if ($need_to_push_line == 1) {
				push ( @converted_input_array, $converted_input_line );
				$need_to_push_line = 0;
			}

			my @bits = split(/\|/, $fasta_line);
			my $seqid = "SEQ_$i";
			if (defined $bits[3]) {
				my @bits2 = split(/\./, $bits[3]);
				$seqid = uc $bits2[0];
				$seqid =~ s/\W//g;
				$seqid =~ s/_//g;
				$seqid = "SEQ_$seqid";
			}
			$converted_input_line = $seqid . ':::::' . $fasta_line . ':::::';
			$need_to_push_line = 1;

		} elsif ($fasta_line eq '') {

			if ($need_to_push_line == 1) {
				push ( @converted_input_array, $converted_input_line );
				$need_to_push_line = 0;
			}

		} else {

			$fasta_line =~ s/\s//g;
			$converted_input_line .= $fasta_line;
		}
	}

	if ($need_to_push_line == 1) {
		push ( @converted_input_array, $converted_input_line );
		$need_to_push_line = 0;
	}

	my $seq_num = 0;

	foreach my $input_line (@converted_input_array) { # while(!(f.eof()) ){

		chomp( $input_line );
		if ($input_line ne '') {
			$input_line = trim($input_line);
		}
		if ($input_line ne '') {

			my @bits = split(/:::::/, $input_line);

			my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
			$year += 1900;
			$mon += 1;
			$mon = sprintf("%02d", $mon);
			$mday = sprintf("%02d", $mday);
			$hour = sprintf("%02d", $hour);
			$min = sprintf("%02d", $min);
			$sec = sprintf("%02d", $sec);
			my $time_str = "$year-$mon-$mday.$hour:$min:$sec";

			$seq_num++;
			$input_fasta_id = trim($bits[1]);
			$input_sequence_id = trim($bits[0]);
			#if ((defined($input_flag_multiplegraphs)) || (defined($input_flag_benchmark))) {
			print "$time_str : Processing sequence number $seq_num : $input_fasta_id\n";
			#}

			my @aa_seq_chains;
			$aa_seq_chains[0] = trim($bits[2]);
			$input_aa_sequence = trim($bits[2]);
			my $ref_aa_seq_chains = \@aa_seq_chains;

			$analyser = {}; # initialise this data structure for each new sequence to be read in and analysed

			get_2nd_struc_info();

			CValpredDoc_OnToolsTmpred2d( $ref_aa_seq_chains, $seq_num );
		}
	}

	close SEGOUTFILE;
	close SCOREOUTFILE;
	if (defined($input_flag_multiplegraphs)) {
		close SEGPTSOUTFILE;
		create_output_graphs_html_file( 3 );
		close HTMLGRAPHSOUTFILE;
	}
}



sub program_output_command_line_for_graphpointsfile {

	if (defined($input_flag_multiplegraphs)) {

		my $html_graphs_output_file = "$input_file_name_for_output_files.valpred_graphs.html";
		open( HTMLGRAPHSOUTFILE, ">$html_graphs_output_file") or
			die "Cannot open $html_graphs_output_file for writing : $!\n";
		create_output_graphs_html_file( 1 );
	}

	read_amino_acid_pdbs();
	get_command_line_program_options();

	# open input file
	open INGRAPHPTSFILE, $input_graph_points_file or die $!;
	my @input_lines = <INGRAPHPTSFILE>;

	my $seq_num = 0;

	foreach my $input_line (@input_lines) { # while(!(f.eof()) ){

		chomp( $input_line );
		if ($input_line ne '') {
			$input_line = trim($input_line);
		}
		if ($input_line ne '') {

			my @bits = split(/:::::/, $input_line);

			my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
			$year += 1900;
			$mon += 1;
			$mon = sprintf("%02d", $mon);
			$mday = sprintf("%02d", $mday);
			$hour = sprintf("%02d", $hour);
			$min = sprintf("%02d", $min);
			$sec = sprintf("%02d", $sec);
			my $time_str = "$year-$mon-$mday.$hour:$min:$sec";

			$seq_num++;
			$input_fasta_id = trim($bits[1]);
			$input_sequence_id = trim($bits[0]);
			#if ((defined($input_flag_multiplegraphs)) || (defined($input_flag_benchmark))) {
			print "$time_str : Processing sequence number $seq_num : $input_fasta_id\n";
			#}

			my @aa_seq_chains;
			$input_aa_sequence = trim($bits[2]);
			$input_aa_sequence =~ s/\s//g;
			$aa_seq_chains[0] = $input_aa_sequence;
			my $ref_aa_seq_chains = \@aa_seq_chains;

			# input file :
			# SEQ_1BL8A:::::>1BL8:A|PDBID|CHAIN|1BL8:A:::::ALHWRAAGAATVLLVIVLLAGSYLAVLAERGAPGAQLITYPRALWWSVETATTVGYGDLYPVTLWGRCVAVVVMVAGITSFGLVTAALATWFVGREQ:::::0-20,49-88,:::::PROFILES3D,,,,,-0.518,-0.358333333333333,-0.244285714285714,-0.135,-0.0711111111111112,-0.02,-0.084,-0.072,-0.112,-0.14,-0.059,-0.109,-0.278,-0.478,-0.568,-0.568,-0.485,-0.344,-0.353,-0.262,-0.248,-0.304,-0.225,-0.044,0.062,-0.033,-0.033,-0.00500000000000003,0.025,0.134,0.134,0.225,0.134,0.031,-0.049,-0.053,-0.141,-0.129,-0.06,-0.169,-0.322,-0.46,-0.361,-0.427,-0.45,-0.356,-0.287,-0.363,-0.427,-0.506,-0.334,-0.28,-0.179,-0.01,-0.013,-0.107,-0.176,-0.281,-0.281,-0.202,-0.374,-0.256,-0.37,-0.431,-0.355,-0.256,-0.356,-0.356,-0.461,-0.505,-0.521,-0.54,-0.426,-0.645,-0.695,-0.723,-0.779,-0.591,-0.603,-0.638,-0.533,-0.533,-0.552,-0.362,-0.298,-0.334,-0.262,-0.46,-0.448,-0.26,-0.291,-0.275,-0.29,-0.271111111111111,-0.36,-0.382857142857143,-0.265,:::::REPIMPS,,,,,0.206,0.298333333333333,0.364285714285714,0.26125,0.316666666666667,0.361,0.324,0.268,0.432,0.451,0.705,0.74,0.738,0.914,0.968,0.968,0.883,0.856,0.753,0.753,0.755,0.718,0.774,0.72,0.375,0.119,0.119,0.148,-0.0350000000000001,-0.211,-0.211,-0.423,-0.423,-0.388,-0.134,0.073,-0.037,-0.293,-0.061,0.115,0.15,0.399,0.316,0.279,0.025,0.037,0.269,0.488,0.451,0.395,0.238,0.154,0.061,-0.041,0.304,0.292,0.06,0.095,0.095,0.151,0.308,0.235,0.101,0.224,0.168,0.217,0.447,0.447,0.482,0.478,0.441,0.563,0.697,0.713,0.678,0.649,0.703,0.583,0.639,0.587,0.552,0.552,0.674,0.693,0.73,0.722,0.705,0.879,0.823,0.703,0.484,0.193,-0.021,-0.167777777777778,-0.28375,-0.38,-0.628333333333333,:::::TMS,,,,,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,:::::
			# $multiplegraphs_data : PROFILES3D,,,,,-0.518,-0.358333333333333,-0.244285714285714,-0.135,-0.0711111111111112,-0.02,-0.084,-0.072,-0.112,-0.14,-0.059,-0.109,-0.278,-0.478,-0.568,-0.568,-0.485,-0.344,-0.353,-0.262,-0.248,-0.304,-0.225,-0.044,0.062,-0.033,-0.033,-0.00500000000000003,0.025,0.134,0.134,0.225,0.134,0.031,-0.049,-0.053,-0.141,-0.129,-0.06,-0.169,-0.322,-0.46,-0.361,-0.427,-0.45,-0.356,-0.287,-0.363,-0.427,-0.506,-0.334,-0.28,-0.179,-0.01,-0.013,-0.107,-0.176,-0.281,-0.281,-0.202,-0.374,-0.256,-0.37,-0.431,-0.355,-0.256,-0.356,-0.356,-0.461,-0.505,-0.521,-0.54,-0.426,-0.645,-0.695,-0.723,-0.779,-0.591,-0.603,-0.638,-0.533,-0.533,-0.552,-0.362,-0.298,-0.334,-0.262,-0.46,-0.448,-0.26,-0.291,-0.275,-0.29,-0.271111111111111,-0.36,-0.382857142857143,-0.265,:::::REPIMPS,,,,,0.206,0.298333333333333,0.364285714285714,0.26125,0.316666666666667,0.361,0.324,0.268,0.432,0.451,0.705,0.74,0.738,0.914,0.968,0.968,0.883,0.856,0.753,0.753,0.755,0.718,0.774,0.72,0.375,0.119,0.119,0.148,-0.0350000000000001,-0.211,-0.211,-0.423,-0.423,-0.388,-0.134,0.073,-0.037,-0.293,-0.061,0.115,0.15,0.399,0.316,0.279,0.025,0.037,0.269,0.488,0.451,0.395,0.238,0.154,0.061,-0.041,0.304,0.292,0.06,0.095,0.095,0.151,0.308,0.235,0.101,0.224,0.168,0.217,0.447,0.447,0.482,0.478,0.441,0.563,0.697,0.713,0.678,0.649,0.703,0.583,0.639,0.587,0.552,0.552,0.674,0.693,0.73,0.722,0.705,0.879,0.823,0.703,0.484,0.193,-0.021,-0.167777777777778,-0.28375,-0.38,-0.628333333333333,:::::TMS,,,,,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,:::::TITLE,,,,,>1BL8:A|PDBID|CHAIN|1BL8:A:::::

			$multiplegraphs_data = $bits[4] . ':::::' . $bits[5] . ':::::' . $bits[6] . ':::::' . 'TITLE,,,,,' . $input_fasta_id . ':::::';

			get_2nd_struc_info();

			CValpredDoc_OnToolsTmpred2d( $ref_aa_seq_chains, $seq_num );
		}
	}

	if (defined($input_flag_multiplegraphs)) {
		create_output_graphs_html_file( 3 );
		close HTMLGRAPHSOUTFILE;
	}
}



sub program_output_command_line_for_infile {

	# open output file
	$segment_output_file = "$input_file_name_for_output_files.valpred_segments.txt";
	$scores_output_file = "$input_file_name_for_output_files.valpred_scores.txt";
	$amphipathic_output_file = "$input_file_name_for_output_files.valpred_amphipathic.txt";

	open( SEGOUTFILE, ">$segment_output_file") or
		die "Cannot open $segment_output_file for writing : $!\n";
	open( SCOREOUTFILE, ">$scores_output_file") or
		die "Cannot open $scores_output_file for writing : $!\n";
	if (defined($input_flag_amphipathic)) {
		open( AMPHIFILE, ">$amphipathic_output_file") or
			die "Cannot open $amphipathic_output_file for writing : $!\n";
	}

	if (defined($input_flag_multiplegraphs)) {

		my $segment_graphpoints_output_file = "$input_file_name_for_output_files.valpred_graphpoints.txt";
		open( SEGPTSOUTFILE, ">$segment_graphpoints_output_file") or
			die "Cannot open $segment_graphpoints_output_file for writing : $!\n";

		my $html_graphs_output_file = "$input_file_name_for_output_files.valpred_graphs.html";
		open( HTMLGRAPHSOUTFILE, ">$html_graphs_output_file") or
			die "Cannot open $html_graphs_output_file for writing : $!\n";
		create_output_graphs_html_file( 1 );
	}

	read_amino_acid_pdbs();
	get_command_line_program_options();

	# open input file
	open INFILE, $input_infile or die $!;
	my @input_lines = <INFILE>;

	my $seq_num = 0;

	foreach my $input_line (@input_lines) { # while(!(f.eof()) ){

		chomp( $input_line );
		if ($input_line ne '') {
			$input_line = trim($input_line);
		}
		if ($input_line ne '') {

			my @bits = split(/:::::/, $input_line);

			my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
			$year += 1900;
			$mon += 1;
			$mon = sprintf("%02d", $mon);
			$mday = sprintf("%02d", $mday);
			$hour = sprintf("%02d", $hour);
			$min = sprintf("%02d", $min);
			$sec = sprintf("%02d", $sec);
			my $time_str = "$year-$mon-$mday.$hour:$min:$sec";

			$seq_num++;
			$input_fasta_id = trim($bits[1]);
			$input_sequence_id = trim($bits[0]);
			#if ((defined($input_flag_multiplegraphs)) || (defined($input_flag_benchmark))) {
			print "$time_str : Processing sequence number $seq_num : $input_fasta_id\n";
			#}

			my @aa_seq_chains;
			$aa_seq_chains[0] = trim($bits[2]);
			$input_aa_sequence = trim($bits[2]);
			my $ref_aa_seq_chains = \@aa_seq_chains;

			$analyser = {}; # initialise this data structure for each new sequence to be read in and analysed

			get_2nd_struc_info();

			CValpredDoc_OnToolsTmpred2d( $ref_aa_seq_chains, $seq_num );
		}
	}

	close SEGOUTFILE;
	close SCOREOUTFILE;
	if (defined($input_flag_multiplegraphs)) {
		close SEGPTSOUTFILE;
		create_output_graphs_html_file( 3 );
		close HTMLGRAPHSOUTFILE;
	}
}



sub program_output_command_line_for_inscoresfile {

	# open output file
	$segment_output_file = "$input_file_name_for_output_files.valpred_segments.txt";
	$scores_output_file = "$input_file_name_for_output_files.valpred_scores.txt";
	$amphipathic_output_file = "$input_file_name_for_output_files.valpred_amphipathic.txt";

	open( SEGOUTFILE, ">$segment_output_file") or
		die "Cannot open $segment_output_file for writing : $!\n";
	open( SCOREOUTFILE, ">$scores_output_file") or
		die "Cannot open $scores_output_file for writing : $!\n";
	if (defined($input_flag_amphipathic)) {
		open( AMPHIFILE, ">$amphipathic_output_file") or
			die "Cannot open $amphipathic_output_file for writing : $!\n";
	}

	if (defined($input_flag_multiplegraphs)) {

		my $segment_graphpoints_output_file = "$input_file_name_for_output_files.valpred_graphpoints.txt";
		open( SEGPTSOUTFILE, ">$segment_graphpoints_output_file") or
			die "Cannot open $segment_graphpoints_output_file for writing : $!\n";

		my $html_graphs_output_file = "$input_file_name_for_output_files.valpred_graphs.html";
		open( HTMLGRAPHSOUTFILE, ">$html_graphs_output_file") or
			die "Cannot open $html_graphs_output_file for writing : $!\n";
		create_output_graphs_html_file( 1 );
	}

	read_amino_acid_pdbs();
	get_command_line_program_options();

	# open input file
	open INFILE, $input_inscoresfile or die $!;
	my @input_lines = <INFILE>;

	my $seq_num = 0;

	foreach my $input_line (@input_lines) { # while(!(f.eof()) ){

		chomp( $input_line );
		if ($input_line ne '') {
			$input_line = trim($input_line);
		}
		if ($input_line ne '') {

			my @bits = split(/:::::/, $input_line);

			my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
			$year += 1900;
			$mon += 1;
			$mon = sprintf("%02d", $mon);
			$mday = sprintf("%02d", $mday);
			$hour = sprintf("%02d", $hour);
			$min = sprintf("%02d", $min);
			$sec = sprintf("%02d", $sec);
			my $time_str = "$year-$mon-$mday.$hour:$min:$sec";

			$seq_num++;
			$input_fasta_id = trim($bits[1]);
			$input_sequence_id = trim($bits[0]);
			#if ((defined($input_flag_multiplegraphs)) || (defined($input_flag_benchmark))) {
			print "$time_str : Processing sequence number $seq_num : $input_fasta_id\n";
			#}

			my @aa_seq_chains;
			$aa_seq_chains[0] = trim($bits[2]);
			$input_aa_sequence = trim($bits[2]);
			my $ref_aa_seq_chains = \@aa_seq_chains;

			my $profiles3d_scores_and_hdr =  trim($bits[3]);
			my $repimps_scores_and_hdr =  trim($bits[4]);
			my @bits1 = split(/\;\;\;\;\;/, $profiles3d_scores_and_hdr);
			my $profiles3d_scores_string = $bits1[1];
			@bits1 = split(/\;\;\;\;\;/, $repimps_scores_and_hdr);
			my $repimps_scores_string = $bits1[1];

			$analyser = {}; # initialise this data structure for each new sequence to be read in and analysed

			get_2nd_struc_info();

			create_structure_from_input_scores( $profiles3d_scores_string, $repimps_scores_string, $input_aa_sequence );

			CValpredDoc_OnToolsTmpred2d( $ref_aa_seq_chains, $seq_num );
		}
	}

	close SEGOUTFILE;
	close SCOREOUTFILE;
	if (defined($input_flag_multiplegraphs)) {
		close SEGPTSOUTFILE;
		create_output_graphs_html_file( 3 );
		close HTMLGRAPHSOUTFILE;
	}
}



sub program_output_command_line_for_predictions {

	get_command_line_program_options();

	# open input file
	open PREDFILE, $input_predictions_file or die $!;
	my @input_lines = <PREDFILE>;

	my $seq_num = 0;

	foreach my $input_line (@input_lines) { # while(!(f.eof()) ){

		chomp( $input_line );
		if ($input_line ne '') {
			$input_line = trim($input_line);
		}
		if ($input_line ne '') {

			my @bits = split(/:::::/, $input_line);

			my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
			$year += 1900;
			$mon += 1;
			$mon = sprintf("%02d", $mon);
			$mday = sprintf("%02d", $mday);
			$hour = sprintf("%02d", $hour);
			$min = sprintf("%02d", $min);
			$sec = sprintf("%02d", $sec);
			my $time_str = "$year-$mon-$mday.$hour:$min:$sec";

			$seq_num++;
			$input_fasta_id = trim($bits[1]);
			$input_sequence_id = trim($bits[0]);
			#if ((defined($input_flag_multiplegraphs)) || (defined($input_flag_benchmark))) {
			print "$time_str : Processing sequence number $seq_num : $input_fasta_id\n";
			#}

			my @aa_seq_chains;
			$aa_seq_chains[0] = trim($bits[2]);
			$input_aa_sequence = trim($bits[2]);
			my $segments = trim($bits[3]);

			get_2nd_struc_info();

			if (defined($input_tms_file)) {
				my @predicted_tms_array_start;
				my @predicted_tms_array_end;
				if ($segments ne '-') {
					my @bits2 = split(/\,/, $segments);
					foreach my $segment (@bits2) {
						my @bits3 = split(/\-/, $segment);
						my $start = $bits3[0];
						my $end = $bits3[1];
						push( @predicted_tms_array_start, $start );
						push( @predicted_tms_array_end, $end );
					}
				}
				my $predicted_tms_array_start_ref = \@predicted_tms_array_start;
				my $predicted_tms_array_end_ref = \@predicted_tms_array_end;
				output_evaluate_file( $predicted_tms_array_start_ref, $predicted_tms_array_end_ref );
			}
		}
	}
}



sub program_output_html_graph {

	my $params = $q->Vars;
	my $chart_coords = '';
	if (defined($params->{'chart_coords'})) {
		$chart_coords = trim($params->{'chart_coords'});
	}

	set_program_options_defaults();
	# read_program_options_defaults_file(); # defaults file from windows version

	# OPTIONS,3,PROFILES3D-MOV-AVG-1,1,1,10,0,0,255,blue,Profiles 3D (mv.av.1),0,1,:::::
	# OPTIONS,4,REPIMPS-MOV-AVG-1,1,1,10,255,0,0,red,REPIMPS (mv.av.1),0,1,:::::
	# OPTIONS,15,TMS,1,1,0,0,0,0,black,Transmembrane Segments,1.7,1,:::::
	# NUM_RESIDUES,160,:::::NUM_CHAINS,1,:::::CHAIN,0,160,0,:::::
	# LINE,3,:::::
	# DATA,-0.185,-0.146,-0.197,-0.264,-0.167,-0.167,-0.229,-0.338,-0.229,-0.234,-0.261,-0.147,-0.085,0.021,-0.032,:::::
	# LINE,4,:::::
	# DATA,0.407,0.101,0.387,0.617,0.354,0.354,0.437,0.613,0.437,0.127,0.017,0.151,0.068,0.07,0.337,:::::
	# LINE,15,:::::
	# DATA,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,:::::
	# TITLE,,,,,>KCSA_STRLI Voltage-gated potassium channel:::::

	# get any program_options in the $chart_coords input string
	for ( my $o = 0; $o <= $program_options->{'option_max_index'}; $o++ ) {
		$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 0;
	}
	# my @option_attributes = ('line_id','showGraph','weight','mov_avg','colour0','colour1','colour2','colour_name','legend_name','y_value','smooth');
	my @option_attributes = @{$program_options->{'option_name'}};
	my @chart_coords_lines = split( /:::::/, $chart_coords );
	foreach my $chart_coords_line (@chart_coords_lines) {
		my @bits = split( /,/, $chart_coords_line );
		if ($bits[0] eq 'OPTIONS') {
			my $ln = '';
			if (defined($bits[1])) {
				$ln = trim($bits[1]);
				if ($ln =~ /^(\d+)$/) { # is it a +ve integer?
					for (my $o2 = 2; $o2 < @bits; $o2++ ) {
						my $o = $o2 - 2;
						my $value = trim($bits[$o2]);
						if ($value ne '') {
							my $option_attribute = $option_attributes[$o];
							$program_options->{'lnOpt'}->[$ln]->{$option_attribute} = $value;
						}
					}
				}
			}
		# } elsif ($bits[0] eq 'LINE') {
		#	my $ln = '';
		#	if (defined($bits[1])) {
		#		$ln = trim($bits[1]);
		#		if ($ln =~ /^(\d+)$/) { # is it a +ve integer?
		#			$program_options->{'lnOpt'}->[$ln]->{'showGraph'} = 1;
		#		}
		#	}
		}
	}

	draw_chart_from_coords( $chart_coords, 0 );
}



sub program_output_html_results {

	html_header();

	read_amino_acid_pdbs();
	get_cgi_script_program_options();

	my $ref_aa_seq_chains = get_html_aaseq();

	my @aa_seq_chains = @$ref_aa_seq_chains;
	if (defined($aa_seq_chains[0])) {
		if (length($aa_seq_chains[0]) > 0) {
			my $seq_num = 1;
			CValpredDoc_OnToolsTmpred2d( $ref_aa_seq_chains, $seq_num );
		} else {
			print "There was no sequence input.\n";
		}
	} else {
		print "There was no input.\n";
	}

	print "</td></tr></table>\n";
	print $q->end_html();
}



sub Range_new {

	my $r;
	$r->{ 'chainID' } = 0;
	$r->{ 'min' } = 0;
	$r->{ 'max' } = 0;

	# print 'in Range_new : $r = ' . Dumper( $r );

	return $r; #C++ Range
}


sub read_amino_acid_pdbs {

	my @amino_acids = ('ala','arg','asn','asp','cys','gln','glu','gly','his','ile','leu','lys','met','phe','pro','ser','thr','trp','tyr','val');

	foreach my $amino_acid (@amino_acids) {

		my $amino_acid_capitals = uc $amino_acid;

		my $filepath = $path_to_files_Amino_Acids . '/' . $amino_acid . '.pdb';	# filepath = currentDir + "/Amino Acids/" + r.name + ".pdb";

		my $file_exists = 1;
		open AAFILE, $filepath or $file_exists = 0;

		my @pdb_lines;
		if ($file_exists == 1) {
			@pdb_lines = <AAFILE>;
		} else {
			my $pdb_lines_ref = hard_coded_amino_acid_pdb( $amino_acid );
			@pdb_lines = @$pdb_lines_ref;
		}

		my $i = 0;
		foreach my $pdb_line (@pdb_lines) {

			chomp( $pdb_line );
			$pdb_line =~ s/\r//g; # windows file line end = 0D0A = \r\n   linux file line end = 0A = \n
			if (length($pdb_line) >= 66) { 				# if(buffer.length() < 4 || buffer.substr(0,4) != "ATOM"){
				if (substr($pdb_line,0,4) eq 'ATOM') {		# if(buffer.substr(0,4) == "ATOM"){

					$amino_acid_pdbs->{$amino_acid_capitals}->[$i]->{'recname'} = trim(substr($pdb_line,0,6));
					$amino_acid_pdbs->{$amino_acid_capitals}->[$i]->{'serial'} = trim(substr($pdb_line,6,6));
					$amino_acid_pdbs->{$amino_acid_capitals}->[$i]->{'name'} = trim(substr($pdb_line,12,4)); 	# string name = buffer.substr(12,4); eatWhite(name);
					$amino_acid_pdbs->{$amino_acid_capitals}->[$i]->{'x'} = trim(substr($pdb_line,30,8)); 		# float x = (float)(atof(buffer.substr(30, 8).c_str()));
					$amino_acid_pdbs->{$amino_acid_capitals}->[$i]->{'y'} = trim(substr($pdb_line,38,8));		# float y = (float)(atof(buffer.substr(38, 8).c_str()));
					$amino_acid_pdbs->{$amino_acid_capitals}->[$i]->{'z'} = trim(substr($pdb_line,46,8));		# float z = (float)(atof(buffer.substr(46, 8).c_str()));
				}
			}
			$i++;
		}
	}
}



sub read_program_options_defaults_file {

	# defaults file from windows version

	# CChartView::setLineParam, input in Valpred_Files/ChartOptions.opt
	# CChartView::loadChartOptions
			# outFile<<"PROFILES\n1\n1\n10\n0\n0\n255\n";
			# outFile<<"REPIMPS\n1\n1\n10\n255\n0\n0\n";
			# outFile<<"OVERALL\n0\n1\n10\n0\n0\n0\n";
			# outFile<<"XAXIS\n1\n20\n1\n0\n0\n \n \n0\n0\n";
			# outFile<<"YAXIS\n1\n1.0\n1\n0\n0\n \n \n0\n0\n";

	my $filepath = $path_to_files_Valpred_Files . '/ChartOptions.opt';

	my $file_exists = 1;
	open DEFOPTS, $filepath or $file_exists = 0;

	if ($file_exists == 1) {

		my @default_options_lines = <DEFOPTS>;

		for ( my $o = 0; $o < @default_options_lines; $o++ ) {
			my $default_options_line = $default_options_lines[$o];
			$default_options_line =~ s/\r//g; # windows file line end = 0D0A = \r\n
			$default_options_line =~ s/\n//g;
			$default_options_lines[$o] = $default_options_line;
		}

		my @option_name_in_ChartOptions_file = ('PROFILES','REPIMPS','OVERALL');
		my @option_name_in_this_program = ('PROFILES3D-MOV-AVG-1','REPIMPS-MOV-AVG-1','BOTH-LINE');

		for ( my $o = 0; $o < @option_name_in_ChartOptions_file; $o++ ) {

			my $ln_name = $option_name_in_this_program[$o];
			my $ln = $program_options->{'lnOpt_index'}->{$ln_name};

			my $f = $ln * 7;
			if (defined($default_options_lines[$f + 1])) {
				$program_options->{'lnOpt'}->[$ln]->{'showGraph'} = $default_options_lines[$f + 1];
			}
			if (defined($default_options_lines[$f + 2])) {
				$program_options->{'lnOpt'}->[$ln]->{'weight'} = $default_options_lines[$f + 2];
			}
			if (defined($default_options_lines[$f + 3])) {
				$program_options->{'lnOpt'}->[$ln]->{'mov_avg'} = $default_options_lines[$f + 3];
			}
			if (defined($default_options_lines[$f + 4])) {
				$program_options->{'lnOpt'}->[$ln]->{'colour0'} = $default_options_lines[$f + 4];
			}
			if (defined($default_options_lines[$f + 5])) {
				$program_options->{'lnOpt'}->[$ln]->{'colour1'} = $default_options_lines[$f + 5];
			}
			if (defined($default_options_lines[$f + 6])) {
				$program_options->{'lnOpt'}->[$ln]->{'colour2'} = $default_options_lines[$f + 6];
			}
		}
	}
}



sub read_program_options_input_file {

	# an example of the options input file for VALPRED1 :

	# SASARES:::::2:::::
	# GRAPHRESDIUEWIDTH:::::10:::::
	# GRAPHHEIGHT:::::600:::::
	# VALPRED1_OR_VALPRED2:::::1:::::
	# MINAREADIFF_FOR_VALPRED1:::::10:::::
	# MINAREADIFF_FOR_VALPRED2:::::10:::::
	# MVEAVE:::::10:::::
	# MOV_AVG_1:::::10:::::
	# MOV_AVG_2:::::5:::::
	# MOV_AVG_3:::::5:::::
	# MOV_AVG_4:::::6:::::
	# MINAVESASA:::::0:::::
	# TMLENGTHMIN_FOR_VALPRED1:::::12:::::
	# TMLENGTHMAX_FOR_VALPRED1:::::40:::::
	# TMLENGTHMIN_FOR_VALPRED2:::::15:::::
	# TMLENGTHMAX_FOR_VALPRED2:::::36:::::
	# AVEREPRANGEMIN:::::0.3:::::
	# AVEREPRANGEMAX:::::0.8:::::
	# AVEPROFRANGEMIN:::::-0.6:::::
	# AVEPROFRANGEMAX:::::0:::::
	# LENDIVAREAMIN:::::0.:::::
	# LENDIVAREAMAX:::::1.8:::::
	# MINSCOREDIFF_FOR_VALPRED1:::::0.2:::::
	# MINSCOREDIFF_MOVAVG1_FOR_VALPRED2:::::0.2:::::
	# MINSCOREDIFF_MOVAVG2_FOR_VALPRED2:::::0.2:::::
	# MIN_AVG_AREA_MOVAVG1:::::0.2:::::
	# MIN_AVG_AREA_MOVAVG2:::::0.2:::::
	# MAX_NONTMS_LENGTH:::::3:::::
	# MIN_NONTMS_LENGTH:::::1:::::
	# MIN_NONTMS_SCORE_DIFF:::::0:::::
	# MIN_NONTMS_AREA:::::1.2:::::
	# TWILIGHT_AREA_PER_RESIDUE_LOWER_LIMIT:::::0.8:::::
	# TWILIGHT_AREA_PER_RESIDUE_UPPER_LIMIT:::::1.3:::::
	# IGNORE_TWILIGHT_AREA:::::0:::::
	# SEARCH_AREA_FOR_NEIGHBOUR_TMS:::::60:::::
	# HELIX_LENGTH_MIN:::::6:::::
	# MIN_SCORE_DIFF_FOR_TMS_ENDS:::::0.6:::::
	# PROFILES3D-LINE:::::0:::::1:::::1:::::0:::::0:::::255:::::cyan:::::Profiles 3D:::::0:::::0:::::
	# REPIMPS-LINE:::::0:::::1:::::1:::::255:::::0:::::0:::::pink:::::REPIMPS:::::0:::::0:::::
	# OVERALL:::::0:::::1:::::10:::::0:::::0:::::0:::::black:::::Overall:::::0:::::1:::::
	# PROFILES3D-MOV-AVG-1:::::1:::::1:::::10:::::0:::::0:::::255:::::blue:::::Profiles 3D (mv.av.1):::::0:::::1:::::
	# REPIMPS-MOV-AVG-1:::::1:::::1:::::10:::::255:::::0:::::0:::::red:::::REPIMPS (mv.av.1):::::0:::::1:::::
	# TMS-MOV-AVG-1:::::0:::::1:::::0:::::0:::::0:::::0:::::red:::::TMS (mv.av.1):::::1.6:::::0:::::
	# PROFILES3D-MOV-AVG-2:::::0:::::1:::::5:::::0:::::0:::::255:::::lblue:::::Profiles 3D (mv.av.2):::::0:::::0:::::
	# REPIMPS-MOV-AVG-2:::::0:::::1:::::5:::::255:::::0:::::0:::::lred:::::REPIMPS (mv.av.2):::::0:::::0:::::
	# TMS-MOV-AVG-2:::::0:::::1:::::0:::::0:::::0:::::0:::::lred:::::TMS (mv.av.2):::::1.58:::::0:::::
	# PROFILES3D-MOV-AVG-3:::::0:::::1:::::5:::::0:::::0:::::255:::::lblue:::::Profiles 3D (mv.av.3):::::0:::::0:::::
	# REPIMPS-MOV-AVG-3:::::0:::::1:::::5:::::255:::::0:::::0:::::lred:::::REPIMPS (mv.av.3):::::0:::::0:::::
	# HYDROPHILIC-MOV-AVG-3:::::0:::::1:::::0:::::0:::::0:::::0:::::lblue:::::Hydrophilic area (mv.av.3):::::1.56:::::0:::::
	# PROFILES3D-MOV-AVG-4:::::0:::::1:::::10:::::0:::::0:::::255:::::cyan:::::Profiles 3D (mv.av.4):::::0:::::1:::::
	# REPIMPS-MOV-AVG-4:::::0:::::1:::::10:::::255:::::0:::::0:::::pink:::::REPIMPS (mv.av.4):::::0:::::1:::::
	# HYDROPHILIC-MOV-AVG-4:::::0:::::1:::::0:::::0:::::0:::::0:::::cyan:::::Hydrophilic area (mv.av.4):::::1.54:::::0:::::
	# TMS:::::1:::::1:::::0:::::0:::::0:::::0:::::black:::::Transmembrane Segments:::::1.7:::::0:::::
	# TRUE-TMS:::::1:::::1:::::10:::::0:::::0:::::0:::::green:::::true TMS:::::1.9:::::0:::::
	# TRUE-INSIDE:::::1:::::1:::::10:::::0:::::0:::::0:::::lgreen:::::true inside:::::1.86:::::0:::::
	# TRUE-OUTSIDE:::::1:::::1:::::10:::::0:::::0:::::0:::::dgreen:::::true outside:::::1.84:::::0:::::

	# an example of the options input file for VALPRED2 :

	# SASARES:::::2:::::
	# GRAPHRESDIUEWIDTH:::::10:::::
	# GRAPHHEIGHT:::::600:::::
	# VALPRED1_OR_VALPRED2:::::2:::::
	# MINAREADIFF_FOR_VALPRED1:::::10:::::
	# MINAREADIFF_FOR_VALPRED2:::::10:::::
	# MVEAVE:::::10:::::
	# MOV_AVG_1:::::10:::::
	# MOV_AVG_2:::::5:::::
	# MOV_AVG_3:::::5:::::
	# MOV_AVG_4:::::6:::::
	# MINAVESASA:::::0:::::
	# TMLENGTHMIN_FOR_VALPRED1:::::12:::::
	# TMLENGTHMAX_FOR_VALPRED1:::::40:::::
	# TMLENGTHMIN_FOR_VALPRED2:::::6:::::
	# TMLENGTHMAX_FOR_VALPRED2:::::38:::::
	# AVEREPRANGEMIN:::::0.3:::::
	# AVEREPRANGEMAX:::::0.8:::::
	# AVEPROFRANGEMIN:::::-0.6:::::
	# AVEPROFRANGEMAX:::::0:::::
	# LENDIVAREAMIN:::::0.:::::
	# LENDIVAREAMAX:::::1.8:::::
	# MINSCOREDIFF_FOR_VALPRED1:::::0.2:::::
	# MINSCOREDIFF_MOVAVG1_FOR_VALPRED2:::::0.2:::::
	# MINSCOREDIFF_MOVAVG2_FOR_VALPRED2:::::0.2:::::
	# MIN_AVG_AREA_MOVAVG1:::::0.2:::::
	# MIN_AVG_AREA_MOVAVG2:::::0.2:::::
	# MAX_NONTMS_LENGTH:::::3:::::
	# MIN_NONTMS_LENGTH:::::1:::::
	# MIN_NONTMS_SCORE_DIFF:::::0:::::
	# MIN_NONTMS_AREA:::::1.2:::::
	# TWILIGHT_AREA_PER_RESIDUE_LOWER_LIMIT:::::0.8:::::
	# TWILIGHT_AREA_PER_RESIDUE_UPPER_LIMIT:::::1.3:::::
	# IGNORE_TWILIGHT_AREA:::::0:::::
	# SEARCH_AREA_FOR_NEIGHBOUR_TMS:::::60:::::
	# HELIX_LENGTH_MIN:::::6:::::
	# MIN_SCORE_DIFF_FOR_TMS_ENDS:::::0.6:::::
	# PROFILES3D-LINE:::::1:::::1:::::1:::::0:::::0:::::255:::::cyan:::::Profiles 3D:::::0:::::0:::::
	# REPIMPS-LINE:::::1:::::1:::::1:::::255:::::0:::::0:::::pink:::::REPIMPS:::::0:::::0:::::
	# OVERALL:::::0:::::1:::::10:::::0:::::0:::::0:::::black:::::Overall:::::0:::::1:::::
	# PROFILES3D-MOV-AVG-1:::::1:::::1:::::10:::::0:::::0:::::255:::::blue:::::Profiles 3D (mv.av.1):::::0:::::1:::::
	# REPIMPS-MOV-AVG-1:::::1:::::1:::::10:::::255:::::0:::::0:::::red:::::REPIMPS (mv.av.1):::::0:::::1:::::
	# TMS-MOV-AVG-1:::::1:::::1:::::0:::::0:::::0:::::0:::::red:::::TMS (mv.av.1):::::1.6:::::0:::::
	# PROFILES3D-MOV-AVG-2:::::1:::::1:::::5:::::0:::::0:::::255:::::lblue:::::Profiles 3D (mv.av.2):::::0:::::0:::::
	# REPIMPS-MOV-AVG-2:::::1:::::1:::::5:::::255:::::0:::::0:::::lred:::::REPIMPS (mv.av.2):::::0:::::0:::::
	# TMS-MOV-AVG-2:::::1:::::1:::::0:::::0:::::0:::::0:::::lred:::::TMS (mv.av.2):::::1.58:::::0:::::
	# PROFILES3D-MOV-AVG-3:::::1:::::1:::::5:::::0:::::0:::::255:::::lblue:::::Profiles 3D (mv.av.3):::::0:::::0:::::
	# REPIMPS-MOV-AVG-3:::::1:::::1:::::5:::::255:::::0:::::0:::::lred:::::REPIMPS (mv.av.3):::::0:::::0:::::
	# HYDROPHILIC-MOV-AVG-3:::::1:::::1:::::0:::::0:::::0:::::0:::::lblue:::::Hydrophilic area (mv.av.3):::::1.56:::::0:::::
	# PROFILES3D-MOV-AVG-4:::::1:::::1:::::10:::::0:::::0:::::255:::::cyan:::::Profiles 3D (mv.av.4):::::0:::::1:::::
	# REPIMPS-MOV-AVG-4:::::1:::::1:::::10:::::255:::::0:::::0:::::pink:::::REPIMPS (mv.av.4):::::0:::::1:::::
	# HYDROPHILIC-MOV-AVG-4:::::1:::::1:::::0:::::0:::::0:::::0:::::cyan:::::Hydrophilic area (mv.av.4):::::1.54:::::0:::::
	# TMS:::::1:::::1:::::0:::::0:::::0:::::0:::::black:::::Transmembrane Segments:::::1.7:::::0:::::
	# TRUE-TMS:::::1:::::1:::::10:::::0:::::0:::::0:::::green:::::true TMS:::::1.9:::::0:::::
	# TRUE-INSIDE:::::1:::::1:::::10:::::0:::::0:::::0:::::lgreen:::::true inside:::::1.86:::::0:::::
	# TRUE-OUTSIDE:::::1:::::1:::::10:::::0:::::0:::::0:::::dgreen:::::true outside:::::1.84:::::0:::::

	# for the lines, the values represent :
	# line_id:::::showGraph:::::weight:::::mov_avg:::::colour0:::::colour1:::::colour2:::::colour_name:::::legend_name:::::y_value:::::smooth:::::

	# $program_options->{'option_name'}->[0] = 'line_id';
	# $program_options->{'option_name'}->[1] = 'showGraph';
	# $program_options->{'option_name'}->[2] = 'weight';
	# $program_options->{'option_name'}->[3] = 'mov_avg';
	# $program_options->{'option_name'}->[4] = 'colour0';
	# $program_options->{'option_name'}->[5] = 'colour1';
	# $program_options->{'option_name'}->[6] = 'colour2';
	# $program_options->{'option_name'}->[7] = 'colour_name';
	# $program_options->{'option_name'}->[8] = 'legend_name';
	# $program_options->{'option_name'}->[9] = 'y_value';
	# $program_options->{'option_name'}->[10] = 'smooth';

	my $o = $program_options->{'lnOpt_index'}->{'TRUE-TMS'};
	$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 0;
	$o = $program_options->{'lnOpt_index'}->{'TRUE-INSIDE'};
	$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 0;
	$o = $program_options->{'lnOpt_index'}->{'TRUE-OUTSIDE'};
	$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 0;

	if (defined($input_options_file)) {
		open OPTSFILE, $input_options_file or die $!;
		my @options_lines = <OPTSFILE>;
		foreach my $options_line (@options_lines) {
			my @bits = split(/:::::/, $options_line);
			my $o = -1;

			if ($bits[0] eq 'SASARES') {
				$program_options->{'SASARes'} = $bits[1];
			} elsif ($bits[0] eq 'GRAPHRESDIUEWIDTH') {
				$program_options->{'graph_residue_width'} = $bits[1];
			} elsif ($bits[0] eq 'GRAPHHEIGHT') {
				$program_options->{'graph_height'} = $bits[1];

			} elsif ($bits[0] eq 'VALPRED1_OR_VALPRED2') {
				$TMPred2D->{'valpred1_or_valpred2'} = $bits[1];
			} elsif ($bits[0] eq 'PREDICT_AMPHIPATHIC_HELICES') {
				$TMPred2D->{'predict_amphipathic_helices'} = $bits[1];
			} elsif ($bits[0] eq 'MINAREADIFF_FOR_VALPRED1') {
				$TMPred2D->{'minAreaDiff_for_valpred1'} = $bits[1];
			} elsif ($bits[0] eq 'MINAREADIFF_FOR_VALPRED2') {
				$TMPred2D->{'minAreaDiff_for_valpred2'} = $bits[1];
			} elsif ($bits[0] eq 'MVEAVE') {
				$TMPred2D->{'mveAve'} = $bits[1];
			} elsif ($bits[0] eq 'MOV_AVG_1') {
				$TMPred2D->{'mov_avg_1'} = $bits[1];
			} elsif ($bits[0] eq 'MOV_AVG_2') {
				$TMPred2D->{'mov_avg_2'} = $bits[1];
			} elsif ($bits[0] eq 'MOV_AVG_3') {
				$TMPred2D->{'mov_avg_3'} = $bits[1];
			} elsif ($bits[0] eq 'MOV_AVG_4') {
				$TMPred2D->{'mov_avg_4'} = $bits[1];
			} elsif ($bits[0] eq 'MINAVESASA') {
				$TMPred2D->{'minAveSasa'} = $bits[1];
			} elsif ($bits[0] eq 'TMLENGTHMIN_FOR_VALPRED1') {
				$TMPred2D->{'TMLengthMin_for_valpred1'} = $bits[1];
			} elsif ($bits[0] eq 'TMLENGTHMAX_FOR_VALPRED1') {
				$TMPred2D->{'TMLengthMax_for_valpred1'} = $bits[1];
			} elsif ($bits[0] eq 'TMLENGTHMIN_FOR_VALPRED2') {
				$TMPred2D->{'TMLengthMin_for_valpred2'} = $bits[1];
			} elsif ($bits[0] eq 'TMLENGTHMAX_FOR_VALPRED2') {
				$TMPred2D->{'TMLengthMax_for_valpred2'} = $bits[1];
			} elsif ($bits[0] eq 'AVEREPRANGEMIN') {
				$TMPred2D->{'aveRepRangeMin'} = $bits[1];
			} elsif ($bits[0] eq 'AVEREPRANGEMAX') {
				$TMPred2D->{'aveRepRangeMax'} = $bits[1];
			} elsif ($bits[0] eq 'AVEPROFRANGEMIN') {
				$TMPred2D->{'aveProfRangeMin'} = $bits[1];
			} elsif ($bits[0] eq 'AVEPROFRANGEMAX') {
				$TMPred2D->{'aveProfRangeMax'} = $bits[1];
			} elsif ($bits[0] eq 'LENDIVAREAMIN') {
				$TMPred2D->{'lenDivAreaMin'} = $bits[1];
			} elsif ($bits[0] eq 'LENDIVAREAMAX') {
				$TMPred2D->{'lenDivAreaMax'} = $bits[1];
			} elsif ($bits[0] eq 'MINSCOREDIFF_FOR_VALPRED1') {
				$TMPred2D->{'minScoreDiff_for_valpred1'} = $bits[1];
			} elsif ($bits[0] eq 'MINSCOREDIFF_MOVAVG1_FOR_VALPRED2') {
				$TMPred2D->{'minScoreDiff_movavg1_for_valpred2'} = $bits[1];
			} elsif ($bits[0] eq 'MINSCOREDIFF_MOVAVG2_FOR_VALPRED2') {
				$TMPred2D->{'minScoreDiff_movavg2_for_valpred2'} = $bits[1];
			} elsif ($bits[0] eq 'MIN_AVG_AREA_MOVAVG1') {
				$TMPred2D->{'min_avg_area_movavg1'} = $bits[1];
			} elsif ($bits[0] eq 'MIN_AVG_AREA_MOVAVG2') {
				$TMPred2D->{'min_avg_area_movavg2'} = $bits[1];
			} elsif ($bits[0] eq 'MAX_NONTMS_LENGTH') {
				$TMPred2D->{'max_nonTMS_length'} = $bits[1];
			} elsif ($bits[0] eq 'MIN_NONTMS_LENGTH') {
				$TMPred2D->{'min_nonTMS_length'} = $bits[1];
			} elsif ($bits[0] eq 'MIN_NONTMS_SCORE_DIFF') {
				$TMPred2D->{'min_nonTMS_score_diff'} = $bits[1];
			} elsif ($bits[0] eq 'MIN_NONTMS_AREA') {
				$TMPred2D->{'min_nonTMS_area'} = $bits[1];
			} elsif ($bits[0] eq 'TWILIGHT_AREA_PER_RESIDUE_LOWER_LIMIT') {
				$TMPred2D->{'twilight_area_per_residue_lower_limit'} = $bits[1];
			} elsif ($bits[0] eq 'TWILIGHT_AREA_PER_RESIDUE_UPPER_LIMIT') {
				$TMPred2D->{'twilight_area_per_residue_upper_limit'} = $bits[1];
			} elsif ($bits[0] eq 'IGNORE_TWILIGHT_AREA') {
				$TMPred2D->{'ignore_twilight_area'} = $bits[1];
			} elsif ($bits[0] eq 'SEARCH_AREA_FOR_NEIGHBOUR_TMS') {
				$TMPred2D->{'search_area_for_neighbour_tms'} = $bits[1];
			} elsif ($bits[0] eq 'HELIX_LENGTH_MIN') {
				$TMPred2D->{'helix_length_min'} = $bits[1];
			} elsif ($bits[0] eq 'MIN_SCORE_DIFF_FOR_TMS_ENDS') {
				$TMPred2D->{'min_score_diff_for_TMS_ends'} = $bits[1];

			} elsif ($bits[0] eq 'PROFILES3D-LINE') {
				$o = $program_options->{'lnOpt_index'}->{'PROFILES3D-LINE'};
			} elsif ($bits[0] eq 'REPIMPS-LINE') {
				$o = $program_options->{'lnOpt_index'}->{'REPIMPS-LINE'};
			} elsif ($bits[0] eq 'OVERALL') {
				$o = $program_options->{'lnOpt_index'}->{'BOTH-LINE'};
			} elsif ($bits[0] eq 'PROFILES3D-MOV-AVG-1') {
				$o = $program_options->{'lnOpt_index'}->{'PROFILES3D-MOV-AVG-1'};
			} elsif ($bits[0] eq 'REPIMPS-MOV-AVG-1') {
				$o = $program_options->{'lnOpt_index'}->{'REPIMPS-MOV-AVG-1'};
			} elsif ($bits[0] eq 'TMS-MOV-AVG-1') {
				$o = $program_options->{'lnOpt_index'}->{'TMS-MOV-AVG-1'};
			} elsif ($bits[0] eq 'PROFILES3D-MOV-AVG-2') {
				$o = $program_options->{'lnOpt_index'}->{'PROFILES3D-MOV-AVG-2'};
			} elsif ($bits[0] eq 'REPIMPS-MOV-AVG-2') {
				$o = $program_options->{'lnOpt_index'}->{'REPIMPS-MOV-AVG-2'};
			} elsif ($bits[0] eq 'TMS-MOV-AVG-2') {
				$o = $program_options->{'lnOpt_index'}->{'TMS-MOV-AVG-2'};
			} elsif ($bits[0] eq 'PROFILES3D-MOV-AVG-3') {
				$o = $program_options->{'lnOpt_index'}->{'PROFILES3D-MOV-AVG-3'};
			} elsif ($bits[0] eq 'REPIMPS-MOV-AVG-3') {
				$o = $program_options->{'lnOpt_index'}->{'REPIMPS-MOV-AVG-3'};
			} elsif ($bits[0] eq 'HYDROPHILIC-MOV-AVG-3') {
				$o = $program_options->{'lnOpt_index'}->{'HYDROPHILIC-MOV-AVG-3'};
			} elsif ($bits[0] eq 'PROFILES3D-MOV-AVG-4') {
				$o = $program_options->{'lnOpt_index'}->{'PROFILES3D-MOV-AVG-4'};
			} elsif ($bits[0] eq 'REPIMPS-MOV-AVG-4') {
				$o = $program_options->{'lnOpt_index'}->{'REPIMPS-MOV-AVG-4'};
			} elsif ($bits[0] eq 'HYDROPHILIC-MOV-AVG-4') {
				$o = $program_options->{'lnOpt_index'}->{'HYDROPHILIC-MOV-AVG-4'};
			} elsif ($bits[0] eq 'TMS') {
				$o = $program_options->{'lnOpt_index'}->{'TMS'};
			} elsif ($bits[0] eq 'TRUE-TMS') {
				$o = $program_options->{'lnOpt_index'}->{'TRUE-TMS'};
			} elsif ($bits[0] eq 'TRUE-INSIDE') {
				$o = $program_options->{'lnOpt_index'}->{'TRUE-INSIDE'};
			} elsif ($bits[0] eq 'TRUE-OUTSIDE') {
				$o = $program_options->{'lnOpt_index'}->{'TRUE-OUTSIDE'};
			}

			if ($o != -1) {
				if (defined($bits[1])) {
					if ($bits[1] ne '-') {
						$program_options->{'lnOpt'}->[$o]->{'showGraph'} = $bits[1];
					}
				}
				if (defined($bits[2])) {
					if ($bits[2] ne '-') {
						$program_options->{'lnOpt'}->[$o]->{'weight'} = $bits[2];
					}
				}
				if (defined($bits[3])) {
					if ($bits[3] ne '-') {
						$program_options->{'lnOpt'}->[$o]->{'mov_avg'} = $bits[3];
					}
				}
				if (defined($bits[4])) {
					if ($bits[4] ne '-') {
						$program_options->{'lnOpt'}->[$o]->{'colour0'} = $bits[4];
					}
				}
				if (defined($bits[5])) {
					if ($bits[5] ne '-') {
						$program_options->{'lnOpt'}->[$o]->{'colour1'} = $bits[5];
					}
				}
				if (defined($bits[6])) {
					if ($bits[6] ne '-') {
						$program_options->{'lnOpt'}->[$o]->{'colour2'} = $bits[6];
					}
				}
				if (defined($bits[7])) {
					if ($bits[7] ne '-') {
						$program_options->{'lnOpt'}->[$o]->{'colour_name'} = $bits[7];
					}
				}
				if (defined($bits[8])) {
					if ($bits[8] ne '-') {
						$program_options->{'lnOpt'}->[$o]->{'legend_name'} = $bits[8];
					}
				}
				if (defined($bits[9])) {
					if ($bits[9] ne '-') {
						$program_options->{'lnOpt'}->[$o]->{'y_value'} = $bits[9];
					}
				}
				if (defined($bits[10])) {
					if ($bits[10] ne '-') {
						$program_options->{'lnOpt'}->[$o]->{'smooth'} = $bits[10];
					}
				}
			}
		}
	}

	# read in data for plotting what the TMS should really be
	if (defined($input_tms_file)) {
		open TMSFILE, $input_tms_file or die $!;
		my @tms_lines = <TMSFILE>;
		foreach my $tms_line (@tms_lines) {
			$tms_line = trim($tms_line);
			if ($tms_line ne '') {
				my $tms_index = @input_tms_array_fasta_id;
				my @bits = split( /:::::/, $tms_line );
				my $tms_fasta_id = trim($bits[1]);
				$input_tms_array_fasta_id[$tms_index] = $tms_fasta_id;
				$input_tms_array_tms_line[$tms_index] = '-';
				$input_tms_array_inside_line[$tms_index] = '-';
				$input_tms_array_outside_line[$tms_index] = '-';
				$input_tms_array_HELIX_INSIDE_line[$tms_index] = '-';
				$input_tms_array_HELIX_OUTSIDE_line[$tms_index] = '-';
				$input_tms_array_HELIX_MEMBRANE_line[$tms_index] = '-';
				for ( my $t = 2; $t < @bits; $t++ ) {
					my $tms_2nd_struc_line = $bits[$t];
					my @bits2 = split( /;;;;;/, $tms_2nd_struc_line );
					if ($bits2[0] eq 'TMS') {
						$input_tms_array_tms_line[$tms_index] = $bits2[1];
					} elsif ($bits2[0] eq 'INSIDE') {
						$input_tms_array_inside_line[$tms_index] = $bits2[1];
					} elsif ($bits2[0] eq 'OUTSIDE') {
						$input_tms_array_outside_line[$tms_index] = $bits2[1];
					} elsif ($bits2[0] eq 'TMS') {
						$input_tms_array_TMS_line[$tms_index] = $bits2[1];
					} elsif ($bits2[0] eq 'TMS_HELIX') {
						$input_tms_array_TMS_HELIX_line[$tms_index] = $bits2[1];
					} elsif ($bits2[0] eq 'HELIX') {
						$input_tms_array_HELIX_line[$tms_index] = $bits2[1];
					} elsif ($bits2[0] eq 'HELIX_INSIDE') {
						$input_tms_array_HELIX_INSIDE_line[$tms_index] = $bits2[1];
					} elsif ($bits2[0] eq 'HELIX_OUTSIDE') {
						$input_tms_array_HELIX_OUTSIDE_line[$tms_index] = $bits2[1];
					} elsif ($bits2[0] eq 'HELIX_MEMBRANE') {
						$input_tms_array_HELIX_MEMBRANE_line[$tms_index] = $bits2[1];
					}
				}
			}
		}
	}
}



sub Residue_checkMainChain { #C++ bool Residue::checkMainChain(){
	my $r = shift;
	if (($r->{'C'} < 0) || ($r->{'CA'} < 0) || ($r->{'N'} < 0) || ($r->{'O'} < 0)) { # if (C < 0 || CA < 0 || N < 0 || O < 0){
		return 0;								# return false;
	}
	return 1;									# return true;
}



sub Residue_new {

	my $res;

	# ///structural information///
	$res->{ 'chainID' } = -1;		# char chainID;				# // chain Identifier
	$res->{ 'SS' } = '';			# SecStruct SS;				# // secondary structure
	$res->{ 'resSeq' } = -1;		# unsigned int resSeq;			# // residue sequence identifier
	$res->{ 'name' } = '';			# string name;				# // residue name
	$res->{ 'isTransMem' } = -1;		# bool isTransMem;
	$res->{ 'CA' } = -1;			# int CA;
	$res->{ 'C' } = -1;			# int C;
	$res->{ 'N' } = -1;			# int N;
	$res->{ 'O' } = -1;			# int O;				# // index to main chain atoms
	my @res_sideChains_array;
	$res->{ 'sideChains' } = \@res_sideChains_array; # vector<unsigned int> sideChains; # // vector of sidechain atoms
	$res->{ 'aaCode' } = '';		# ResCode aaCode;			# //enumeration residue type;

	# print 'in Residue_new : $res = ' . Dumper( $res );

	return $res; #C++ Residue
}



sub set_assignResType {
	$assignResType{'ALA'} = ALA;
	$assignResType{'ARG'} = ARG;
	$assignResType{'ASN'} = ASN;
	$assignResType{'CYS'} = CYS;
	$assignResType{'GLN'} = GLN;
	$assignResType{'GLU'} = GLU;
	$assignResType{'GLY'} = GLY;
	$assignResType{'HIS'} = HIS;
	$assignResType{'ILE'} = ILE;
	$assignResType{'LEU'} = LEU;
	$assignResType{'LYS'} = LYS;
	$assignResType{'MET'} = MET;
	$assignResType{'PHE'} = PHE;
	$assignResType{'PRO'} = PRO;
	$assignResType{'SER'} = SER;
	$assignResType{'THR'} = THR;
	$assignResType{'TRP'} = TRP;
	$assignResType{'TYR'} = TYR;
	$assignResType{'VAL'} = VAL;
	$assignResType{'UNK'} = UNK;
	return;
}



sub set_charToRes {

	$charToRes{'A'} = "ALA";
	$charToRes{'R'} = "ARG";
	$charToRes{'N'} = "ASN";
	$charToRes{'D'} = "ASP";
	$charToRes{'C'} = "CYS";
	$charToRes{'Q'} = "GLN";
	$charToRes{'E'} = "GLU";
	$charToRes{'G'} = "GLY";
	$charToRes{'H'} = "HIS";
	$charToRes{'I'} = "ILE";
	$charToRes{'L'} = "LEU";
	$charToRes{'K'} = "LYS";
	$charToRes{'M'} = "MET";
	$charToRes{'F'} = "PHE";
	$charToRes{'P'} = "PRO";
	$charToRes{'S'} = "SER";
	$charToRes{'T'} = "THR";
	$charToRes{'W'} = "TRP";
	$charToRes{'Y'} = "TYR";
	$charToRes{'V'} = "VAL";
	return;
}



sub set_program_options_defaults {

	$program_options->{'SASARes'} = 2;
	$program_options->{'graph_residue_width'} = 10; # 10 pixels to be shown per residue
	$program_options->{'graph_height'} = 600; # 600 pixels high

	$TMPred2D->{'valpred1_or_valpred2'} = 1;
	if (defined($input_flag_version)) {
		$input_flag_version = uc $input_flag_version;
		if ($input_flag_version eq 'VALPRED2') {
			$TMPred2D->{'valpred1_or_valpred2'} = 2;
		}
	}

	$TMPred2D->{'predict_amphipathic_helices'} = 0;
	if (defined($input_flag_amphipathic)) {
		$TMPred2D->{'predict_amphipathic_helices'} = 1;
	}

	# VALPRED1 parameters 

	$TMPred2D->{'minAreaDiff_for_valpred1'} = 10; # , minAreaDiff(10.0f)
	$TMPred2D->{'mveAve'} = 10;		# , mveAve(10)
	$TMPred2D->{'minAveSasa'} = 0;		# , minAveSasa(0) # ///not sure about this
	$TMPred2D->{'TMLengthMin_for_valpred1'} = 12; # , TMLengthMin(12)
	$TMPred2D->{'TMLengthMax_for_valpred1'} = 40; # , TMLengthMax(40)
	$TMPred2D->{'aveRepRangeMin'} = 0.3;	# , aveRepRangeMin(0.3f)
	$TMPred2D->{'aveRepRangeMax'} = 0.8;	# , aveRepRangeMax(0.8f)
	$TMPred2D->{'aveProfRangeMin'} = -0.6;	# , aveProfRangeMin(-0.6f)
	$TMPred2D->{'aveProfRangeMax'} = 0;	# , aveProfRangeMax(0.0f)
	$TMPred2D->{'lenDivAreaMin'} = 0.8;	# , lenDivAreaMin(0.8f)
	$TMPred2D->{'lenDivAreaMax'} = 1.8;	# , lenDivAreaMax(1.8f)
	$TMPred2D->{'minScoreDiff_for_valpred1'} = 0.2;	# , minScoreDiff(0.2f)

	# VALPRED2 parameters

	$TMPred2D->{'minAreaDiff_for_valpred2'} = 10; # 6
	$TMPred2D->{'mov_avg_1'} = 10;
	$TMPred2D->{'mov_avg_2'} = 5;
	$TMPred2D->{'mov_avg_3'} = 5;
	$TMPred2D->{'mov_avg_4'} = 6;
	$TMPred2D->{'TMLengthMin_for_valpred2'} = 6;
	$TMPred2D->{'TMLengthMax_for_valpred2'} = 38; # 45;
	$TMPred2D->{'minScoreDiff_movavg1_for_valpred2'} = 0; # 0;
	$TMPred2D->{'minScoreDiff_movavg2_for_valpred2'} = 0.2; # 0;
	$TMPred2D->{'min_avg_area_movavg1'} = 0.2; # 0.6;
	$TMPred2D->{'min_avg_area_movavg2'} = 0.2;
	$TMPred2D->{'max_nonTMS_length'} = 3;
	$TMPred2D->{'min_nonTMS_length'} = 1;
	$TMPred2D->{'min_nonTMS_score_diff'} = 0; # -0.4;
	$TMPred2D->{'min_nonTMS_area'} = 1.2;
	$TMPred2D->{'twilight_area_per_residue_lower_limit'} = 0.8;
	$TMPred2D->{'twilight_area_per_residue_upper_limit'} = 1.3;
	$TMPred2D->{'ignore_twilight_area'} = 0;
	$TMPred2D->{'search_area_for_neighbour_tms'} = 60;
	$TMPred2D->{'helix_length_min'} = 6;
	$TMPred2D->{'min_score_diff_for_TMS_ends'} = 0.1; # -10
	$TMPred2D->{'ASA_are_available'} = 1;

	$program_options->{'option_name'}->[0] = 'line_id';
	$program_options->{'option_name'}->[1] = 'showGraph';
	$program_options->{'option_name'}->[2] = 'weight';
	$program_options->{'option_name'}->[3] = 'mov_avg';
	$program_options->{'option_name'}->[4] = 'colour0';
	$program_options->{'option_name'}->[5] = 'colour1';
	$program_options->{'option_name'}->[6] = 'colour2';
	$program_options->{'option_name'}->[7] = 'colour_name';
	$program_options->{'option_name'}->[8] = 'legend_name';
	$program_options->{'option_name'}->[9] = 'y_value';
	$program_options->{'option_name'}->[10] = 'smooth';
	my $o;

	# PROFILES3D-LINE
	$o = 0;
	$program_options->{'lnOpt_index'}->{'PROFILES3D-LINE'} = $o;
	$program_options->{'lnOpt'}->[$o]->{'line_id'} = 'PROFILES3D-LINE';
	$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 0;
	if (defined($input_flag_version)) {
		if ($input_flag_version eq 'VALPRED2') {
			$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 1;
		}
	}
	$program_options->{'lnOpt'}->[$o]->{'weight'} = 1;
	$program_options->{'lnOpt'}->[$o]->{'mov_avg'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour0'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour1'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour2'} = 255;
	$program_options->{'lnOpt'}->[$o]->{'colour_name'} = 'cyan';
	$program_options->{'lnOpt'}->[$o]->{'legend_name'} = 'Profiles 3D';
	$program_options->{'lnOpt'}->[$o]->{'y_value'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'smooth'} = 0;

	# REPIMPS-LINE
	$o = 1;
	$program_options->{'lnOpt_index'}->{'REPIMPS-LINE'} = $o;
	$program_options->{'lnOpt'}->[$o]->{'line_id'} = 'REPIMPS-LINE';
	$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 0;
	if (defined($input_flag_version)) {
		if ($input_flag_version eq 'VALPRED2') {
			$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 1;
		}
	}
	$program_options->{'lnOpt'}->[$o]->{'weight'} = 1;
	$program_options->{'lnOpt'}->[$o]->{'mov_avg'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour0'} = 255;
	$program_options->{'lnOpt'}->[$o]->{'colour1'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour2'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour_name'} = 'pink';
	$program_options->{'lnOpt'}->[$o]->{'legend_name'} = 'REPIMPS';
	$program_options->{'lnOpt'}->[$o]->{'y_value'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'smooth'} = 0;

	# BOTH-LINE
	$o = 2;
	$program_options->{'lnOpt_index'}->{'BOTH-LINE'} = $o;
	$program_options->{'lnOpt'}->[$o]->{'line_id'} = 'BOTH-LINE';
	$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'weight'} = 1;
	$program_options->{'lnOpt'}->[$o]->{'mov_avg'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour0'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour1'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour2'} = 255;
	$program_options->{'lnOpt'}->[$o]->{'colour_name'} = 'black';
	$program_options->{'lnOpt'}->[$o]->{'legend_name'} = 'Both/Overall';
	$program_options->{'lnOpt'}->[$o]->{'y_value'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'smooth'} = 1;

	# PROFILES3D-MOV-AVG-1
	$o = 3;
	$program_options->{'lnOpt_index'}->{'PROFILES3D-MOV-AVG-1'} = $o;
	$program_options->{'lnOpt'}->[$o]->{'line_id'} = 'PROFILES3D-MOV-AVG-1';
	$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 1;
	$program_options->{'lnOpt'}->[$o]->{'weight'} = 1;
	$program_options->{'lnOpt'}->[$o]->{'mov_avg'} = 10;
	$program_options->{'lnOpt'}->[$o]->{'colour0'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour1'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour2'} = 255;
	$program_options->{'lnOpt'}->[$o]->{'colour_name'} = 'blue';
	$program_options->{'lnOpt'}->[$o]->{'legend_name'} = 'Profiles 3D (mv.av.1)';
	$program_options->{'lnOpt'}->[$o]->{'y_value'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'smooth'} = 1;

	# REPIMPS-MOV-AVG-1
	$o = 4;
	$program_options->{'lnOpt_index'}->{'REPIMPS-MOV-AVG-1'} = $o;
	$program_options->{'lnOpt'}->[$o]->{'line_id'} = 'REPIMPS-MOV-AVG-1';
	$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 1;
	$program_options->{'lnOpt'}->[$o]->{'weight'} = 1;
	$program_options->{'lnOpt'}->[$o]->{'mov_avg'} = 10;
	$program_options->{'lnOpt'}->[$o]->{'colour0'} = 255;
	$program_options->{'lnOpt'}->[$o]->{'colour1'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour2'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour_name'} = 'red';
	$program_options->{'lnOpt'}->[$o]->{'legend_name'} = 'REPIMPS (mv.av.1)';
	$program_options->{'lnOpt'}->[$o]->{'y_value'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'smooth'} = 1;

	# TMS-MOV-AVG-1
	$o = 5;
	$program_options->{'lnOpt_index'}->{'TMS-MOV-AVG-1'} = $o;
	$program_options->{'lnOpt'}->[$o]->{'line_id'} = 'TMS-MOV-AVG-1';
	$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 0;
	if (defined($input_flag_version)) {
		if ($input_flag_version eq 'VALPRED2') {
			$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 1;
		}
	}
	$program_options->{'lnOpt'}->[$o]->{'weight'} = 1;
	$program_options->{'lnOpt'}->[$o]->{'mov_avg'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour0'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour1'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour2'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour_name'} = 'red';
	$program_options->{'lnOpt'}->[$o]->{'legend_name'} = 'TMS (mv.av.1)';
	$program_options->{'lnOpt'}->[$o]->{'y_value'} = 1.6;
	$program_options->{'lnOpt'}->[$o]->{'smooth'} = 0;

	# PROFILES3D-MOV-AVG-2
	$o = 6;
	$program_options->{'lnOpt_index'}->{'PROFILES3D-MOV-AVG-2'} = $o;
	$program_options->{'lnOpt'}->[$o]->{'line_id'} = 'PROFILES3D-MOV-AVG-2';
	$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 0;
	if (defined($input_flag_version)) {
		if ($input_flag_version eq 'VALPRED2') {
			$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 1;
		}
	}
	$program_options->{'lnOpt'}->[$o]->{'weight'} = 1;
	$program_options->{'lnOpt'}->[$o]->{'mov_avg'} = 5;
	$program_options->{'lnOpt'}->[$o]->{'colour0'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour1'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour2'} = 255;
	$program_options->{'lnOpt'}->[$o]->{'colour_name'} = 'lblue';
	$program_options->{'lnOpt'}->[$o]->{'legend_name'} = 'Profiles 3D (mv.av.2)';
	$program_options->{'lnOpt'}->[$o]->{'y_value'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'smooth'} = 0;

	# REPIMPS-MOV-AVG-2
	$o = 7;
	$program_options->{'lnOpt_index'}->{'REPIMPS-MOV-AVG-2'} = $o;
	$program_options->{'lnOpt'}->[$o]->{'line_id'} = 'REPIMPS-MOV-AVG-2';
	$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 0;
	if (defined($input_flag_version)) {
		if ($input_flag_version eq 'VALPRED2') {
			$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 1;
		}
	}
	$program_options->{'lnOpt'}->[$o]->{'weight'} = 1;
	$program_options->{'lnOpt'}->[$o]->{'mov_avg'} = 5;
	$program_options->{'lnOpt'}->[$o]->{'colour0'} = 255;
	$program_options->{'lnOpt'}->[$o]->{'colour1'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour2'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour_name'} = 'lred';
	$program_options->{'lnOpt'}->[$o]->{'legend_name'} = 'REPIMPS (mv.av.2)';
	$program_options->{'lnOpt'}->[$o]->{'y_value'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'smooth'} = 0;

	# TMS-MOV-AVG-2
	$o = 8;
	$program_options->{'lnOpt_index'}->{'TMS-MOV-AVG-2'} = $o;
	$program_options->{'lnOpt'}->[$o]->{'line_id'} = 'TMS-MOV-AVG-2';
	$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 0;
	if (defined($input_flag_version)) {
		if ($input_flag_version eq 'VALPRED2') {
			$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 1;
		}
	}
	$program_options->{'lnOpt'}->[$o]->{'weight'} = 1;
	$program_options->{'lnOpt'}->[$o]->{'mov_avg'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour0'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour1'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour2'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour_name'} = 'lred';
	$program_options->{'lnOpt'}->[$o]->{'legend_name'} = 'TMS (mv.av.2)';
	$program_options->{'lnOpt'}->[$o]->{'y_value'} = 1.58;
	$program_options->{'lnOpt'}->[$o]->{'smooth'} = 0;

	# PROFILES3D-MOV-AVG-3
	$o = 9;
	$program_options->{'lnOpt_index'}->{'PROFILES3D-MOV-AVG-3'} = $o;
	$program_options->{'lnOpt'}->[$o]->{'line_id'} = 'PROFILES3D-MOV-AVG-3';
	$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 0;
	if (defined($input_flag_version)) {
		if ($input_flag_version eq 'VALPRED2') {
			$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 1;
		}
	}
	$program_options->{'lnOpt'}->[$o]->{'weight'} = 1;
	$program_options->{'lnOpt'}->[$o]->{'mov_avg'} = 5;
	$program_options->{'lnOpt'}->[$o]->{'colour0'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour1'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour2'} = 255;
	$program_options->{'lnOpt'}->[$o]->{'colour_name'} = 'purple';
	$program_options->{'lnOpt'}->[$o]->{'legend_name'} = 'Profiles 3D (mv.av.3)';
	$program_options->{'lnOpt'}->[$o]->{'y_value'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'smooth'} = 0;

	# REPIMPS-MOV-AVG-3
	$o = 10;
	$program_options->{'lnOpt_index'}->{'REPIMPS-MOV-AVG-3'} = $o;
	$program_options->{'lnOpt'}->[$o]->{'line_id'} = 'REPIMPS-MOV-AVG-3';
	$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 0;
	if (defined($input_flag_version)) {
		if ($input_flag_version eq 'VALPRED2') {
			$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 1;
		}
	}
	$program_options->{'lnOpt'}->[$o]->{'weight'} = 1;
	$program_options->{'lnOpt'}->[$o]->{'mov_avg'} = 5;
	$program_options->{'lnOpt'}->[$o]->{'colour0'} = 255;
	$program_options->{'lnOpt'}->[$o]->{'colour1'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour2'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour_name'} = 'orange';
	$program_options->{'lnOpt'}->[$o]->{'legend_name'} = 'REPIMPS (mv.av.3)';
	$program_options->{'lnOpt'}->[$o]->{'y_value'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'smooth'} = 0;

	# HYDROPHILIC-MOV-AVG-3
	$o = 11;
	$program_options->{'lnOpt_index'}->{'HYDROPHILIC-MOV-AVG-3'} = $o;
	$program_options->{'lnOpt'}->[$o]->{'line_id'} = 'HYDROPHILIC-MOV-AVG-3';
	$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 0;
	if (defined($input_flag_version)) {
		if ($input_flag_version eq 'VALPRED2') {
			$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 1;
		}
	}
	$program_options->{'lnOpt'}->[$o]->{'weight'} = 1;
	$program_options->{'lnOpt'}->[$o]->{'mov_avg'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour0'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour1'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour2'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour_name'} = 'purple';
	$program_options->{'lnOpt'}->[$o]->{'legend_name'} = 'Hydrophilic area (mv.av.3)';
	$program_options->{'lnOpt'}->[$o]->{'y_value'} = 1.56;
	$program_options->{'lnOpt'}->[$o]->{'smooth'} = 0;

	# PROFILES3D-MOV-AVG-4
	$o = 12;
	$program_options->{'lnOpt_index'}->{'PROFILES3D-MOV-AVG-4'} = $o;
	$program_options->{'lnOpt'}->[$o]->{'line_id'} = 'PROFILES3D-MOV-AVG-4';
	$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 0;
	if (defined($input_flag_version)) {
		if ($input_flag_version eq 'VALPRED2') {
			$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 1;
		}
	}
	$program_options->{'lnOpt'}->[$o]->{'weight'} = 1;
	$program_options->{'lnOpt'}->[$o]->{'mov_avg'} = 10;
	$program_options->{'lnOpt'}->[$o]->{'colour0'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour1'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour2'} = 255;
	$program_options->{'lnOpt'}->[$o]->{'colour_name'} = 'cyan';
	$program_options->{'lnOpt'}->[$o]->{'legend_name'} = 'Profiles 3D (mv.av.4)';
	$program_options->{'lnOpt'}->[$o]->{'y_value'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'smooth'} = 1;

	# REPIMPS-MOV-AVG-4
	$o = 13;
	$program_options->{'lnOpt_index'}->{'REPIMPS-MOV-AVG-4'} = $o;
	$program_options->{'lnOpt'}->[$o]->{'line_id'} = 'REPIMPS-MOV-AVG-4';
	$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 0;
	if (defined($input_flag_version)) {
		if ($input_flag_version eq 'VALPRED2') {
			$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 1;
		}
	}
	$program_options->{'lnOpt'}->[$o]->{'weight'} = 1;
	$program_options->{'lnOpt'}->[$o]->{'mov_avg'} = 10;
	$program_options->{'lnOpt'}->[$o]->{'colour0'} = 255;
	$program_options->{'lnOpt'}->[$o]->{'colour1'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour2'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour_name'} = 'pink';
	$program_options->{'lnOpt'}->[$o]->{'legend_name'} = 'REPIMPS (mv.av.4)';
	$program_options->{'lnOpt'}->[$o]->{'y_value'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'smooth'} = 1;

	# HYDROPHILIC-MOV-AVG-4
	$o = 14;
	$program_options->{'lnOpt_index'}->{'HYDROPHILIC-MOV-AVG-4'} = $o;
	$program_options->{'lnOpt'}->[$o]->{'line_id'} = 'HYDROPHILIC-MOV-AVG-4';
	$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 0;
	if (defined($input_flag_version)) {
		if ($input_flag_version eq 'VALPRED2') {
			$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 1;
		}
	}
	$program_options->{'lnOpt'}->[$o]->{'weight'} = 1;
	$program_options->{'lnOpt'}->[$o]->{'mov_avg'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour0'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour1'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour2'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour_name'} = 'cyan';
	$program_options->{'lnOpt'}->[$o]->{'legend_name'} = 'Hydrophilic area (mv.av.4)';
	$program_options->{'lnOpt'}->[$o]->{'y_value'} = 1.54;
	$program_options->{'lnOpt'}->[$o]->{'smooth'} = 0;

	# TMS (transmembrane segments)
	$o = 15;
	$program_options->{'lnOpt_index'}->{'TMS'} = $o;
	$program_options->{'lnOpt'}->[$o]->{'line_id'} = 'TMS';
	$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 1;
	$program_options->{'lnOpt'}->[$o]->{'weight'} = 1;
	$program_options->{'lnOpt'}->[$o]->{'mov_avg'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour0'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour1'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour2'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour_name'} = 'black';
	$program_options->{'lnOpt'}->[$o]->{'legend_name'} = 'Transmembrane Segments';
	$program_options->{'lnOpt'}->[$o]->{'y_value'} = 1.7;
	$program_options->{'lnOpt'}->[$o]->{'smooth'} = 0;

	# TRUE-TMS
	$o = 16;
	$program_options->{'lnOpt_index'}->{'TRUE-TMS'} = $o;
	$program_options->{'lnOpt'}->[$o]->{'line_id'} = 'TRUE-TMS';
	$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'weight'} = 1;
	$program_options->{'lnOpt'}->[$o]->{'mov_avg'} = 10;
	$program_options->{'lnOpt'}->[$o]->{'colour0'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour1'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour2'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour_name'} = 'green';
	$program_options->{'lnOpt'}->[$o]->{'legend_name'} = 'true TMS';
	$program_options->{'lnOpt'}->[$o]->{'y_value'} = 1.9;
	$program_options->{'lnOpt'}->[$o]->{'smooth'} = 0;

	# TRUE-INSIDE
	$o = 17;
	$program_options->{'lnOpt_index'}->{'TRUE-INSIDE'} = $o;
	$program_options->{'lnOpt'}->[$o]->{'line_id'} = 'TRUE-INSIDE';
	$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'weight'} = 1;
	$program_options->{'lnOpt'}->[$o]->{'mov_avg'} = 10;
	$program_options->{'lnOpt'}->[$o]->{'colour0'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour1'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour2'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour_name'} = 'lgreen';
	$program_options->{'lnOpt'}->[$o]->{'legend_name'} = 'true inside';
	$program_options->{'lnOpt'}->[$o]->{'y_value'} = 1.86;
	$program_options->{'lnOpt'}->[$o]->{'smooth'} = 0;

	# TRUE-OUTSIDE
	$o = 18;
	$program_options->{'lnOpt_index'}->{'TRUE-OUTSIDE'} = $o;
	$program_options->{'lnOpt'}->[$o]->{'line_id'} = 'TRUE-OUTSIDE';
	$program_options->{'lnOpt'}->[$o]->{'showGraph'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'weight'} = 1;
	$program_options->{'lnOpt'}->[$o]->{'mov_avg'} = 10;
	$program_options->{'lnOpt'}->[$o]->{'colour0'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour1'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour2'} = 0;
	$program_options->{'lnOpt'}->[$o]->{'colour_name'} = 'dgreen';
	$program_options->{'lnOpt'}->[$o]->{'legend_name'} = 'true outside';
	$program_options->{'lnOpt'}->[$o]->{'y_value'} = 1.84;
	$program_options->{'lnOpt'}->[$o]->{'smooth'} = 0;

	$program_options->{'option_max_index'} = $o;
}



sub set_Residue_constants {
	# constants set in Residue.h
	$Residue_SecStruct{'SHEET'} = 0;
	$Residue_SecStruct{'HELIX'} = 1;
	$Residue_SecStruct{'OTHER'} = 2;
	$Residue_ResCode{'ALA'} = 7;
	$Residue_ResCode{'ARG'} = 19;
	$Residue_ResCode{'ASN'} = 14;
	$Residue_ResCode{'ASP'} = 16;
	$Residue_ResCode{'CYS'} = 10;
	$Residue_ResCode{'GLN'} = 13;
	$Residue_ResCode{'GLU'} = 15;
	$Residue_ResCode{'GLY'} = 8;
	$Residue_ResCode{'HIS'} = 17;
	$Residue_ResCode{'ILE'} = 4;
	$Residue_ResCode{'LEU'} = 3;
	$Residue_ResCode{'LYS'} = 18;
	$Residue_ResCode{'MET'} = 6;
	$Residue_ResCode{'PHE'} = 1;
	$Residue_ResCode{'PRO'} = 9;
	$Residue_ResCode{'SER'} = 12;
	$Residue_ResCode{'THR'} = 11;
	$Residue_ResCode{'TRP'} = 0;
	$Residue_ResCode{'TYR'} = 2;
	$Residue_ResCode{'VAL'} = 5;
	$Residue_ResCode{'UNK'} = 21;
	$Residue_EnvCode{'B1a'} = 0;
	$Residue_EnvCode{'B1b'} = 1;
	$Residue_EnvCode{'B1'} = 2;
	$Residue_EnvCode{'B2a'} = 3;
	$Residue_EnvCode{'B2b'} = 4;
	$Residue_EnvCode{'B2'} = 5;
	$Residue_EnvCode{'B3a'} = 6;
	$Residue_EnvCode{'B3b'} = 7;
	$Residue_EnvCode{'B3'} = 8;
	$Residue_EnvCode{'P1a'} = 9;
	$Residue_EnvCode{'P1b'} = 10;
	$Residue_EnvCode{'P1'} = 11;
	$Residue_EnvCode{'P2a'} = 12;
	$Residue_EnvCode{'P2b'} = 13;
	$Residue_EnvCode{'P2'} = 14;
	$Residue_EnvCode{'Ea'} = 15;
	$Residue_EnvCode{'Eb'} = 16;
	$Residue_EnvCode{'E'} = 17;
	$Residue_EnvCode{'unknown'} = 18;
	return;
}



sub set_Value_constants {
	# constants set in Value.h
	# //initialisatiion of side chain GlyXGly values
							# const float SCGlyXGly[21] = {
	$SCGlyXGly->[TRP] = 245.263746; 		# 245.263746f,	//TRP
	$SCGlyXGly->[PHE] = 189.159849; 		# 189.159849f,	//PHE
	$SCGlyXGly->[TYR] = 208.929577; 		# 208.929577f,	//TYR
	$SCGlyXGly->[LEU] = 159.793969; 		# 159.793969f,	//LEU
	$SCGlyXGly->[ILE] = 156.162288; 		# 156.162288f,	//ILE
	$SCGlyXGly->[VAL] = 123.11399;	 		# 123.11399f,	//VAL
	$SCGlyXGly->[MET] = 172.169566; 		# 172.169566f,	//MET
	$SCGlyXGly->[ALA] = 60.892862;	 		# 60.892862f,	//ALA
	$SCGlyXGly->[GLY] = 45.032846;	 		# 45.032846f,	//GLY
	$SCGlyXGly->[PRO] = 136.55121; 			# 136.55121f,	//PRO
	$SCGlyXGly->[CYS] = 98.295329; 			# 98.295329f,	//CYS
	$SCGlyXGly->[THR] = 108.913991; 		# 108.913991f,	//THR
	$SCGlyXGly->[SER] = 95.026895; 			# 95.026895f,	//SER
	$SCGlyXGly->[GLN] = 144.73443; 			# 144.73443f,	//GLN
	$SCGlyXGly->[ASN] = 125.607472; 		# 125.607472f,	//ASN
	$SCGlyXGly->[GLU] = 152.160527; 		# 152.160527f,	//GLU
	$SCGlyXGly->[ASP] = 110.354411; 		# 110.354411f,	//ASP
	$SCGlyXGly->[HIS] = 152.417666; 		# 152.417666f,	//HIS
	$SCGlyXGly->[LYS] = 162.613548; 		# 162.613548f,	//LYS
	$SCGlyXGly->[ARG] = 199.132285; 		# 199.132285f,	//ARG
	$SCGlyXGly->[UNK] = 0; 				# 0.000000f,	//UNK
	return;
}



sub SpherePoints_generate { #C++ void SpherePoints::generate(vector<Vec3f>& spherePoints, string path)

	my $spherePoints = shift; #C++ vector<Vec3f>&
	my $path = shift; #C++ string

	my $file_exists = 1;
	open POINTSFILE, $path or $file_exists = 0;			# fstream pointsFile(path.c_str());

	my @lines;
	if ($file_exists == 1) {
		@lines = <POINTSFILE>;					# getline(pointsFile, line);
	} else {
		my $points_lines_ref = hard_coded_sphere_points();
		@lines = @$points_lines_ref;
	}

	$lines[0] =~ s/\r//g; # windows file line end = 0D0A = \r\n   linux file line end = 0A = \n
	my $numPoints = $lines[0];					# int numPoints = atoi(line.c_str());
									# spherePoints.resize(numPoints, Vec3f());
	my $l = 0;
	foreach my $line ( @lines ) {
		chomp( $line );
		$line =~ s/\r//g; # windows file line end = 0D0A = \r\n   linux file line end = 0A = \n
		$lines[$l] = $line;
		$l++;
	}

	my $line_upto = 0;
	for ( my $i = 0; $i < $numPoints; $i++ ) { 			# for (int i=0; i<numPoints; ++i)
		$line_upto++;						# getline(pointsFile, line);
		$spherePoints->[$i]->{'x'} = $lines[$line_upto];	# spherePoints[i].x = (float)atof(line.c_str());
		$line_upto++;						# getline(pointsFile, line);
		$spherePoints->[$i]->{'y'} = $lines[$line_upto];	# spherePoints[i].y = (float)atof(line.c_str());
		$line_upto++;						# getline(pointsFile, line);
		$spherePoints->[$i]->{'z'} = $lines[$line_upto];	# spherePoints[i].z = (float)atof(line.c_str());
	}

	return $spherePoints; #C++ vector<Vec3f>&
}



sub Structure_getResChain_indices {

	my $index = shift;

	my $res;
	$res->{'min'} = -1;
	$res->{'max'} = -1;
	my $chains_dot_size = $#{$structure->{'chains'}} + 1;
	if ($index < $chains_dot_size) {
		$res->{'min'} = $structure->{'chains'}->[$index]->{'min'};
		$res->{'max'} = $structure->{'chains'}->[$index]->{'max'};
		return $res;
	}
	print "ERROR: Internal error: Chain index not found!\n";
	return $res;
}



sub Structure_getResChain_int { #C++ vector<Residue*> Structure::getResChain(unsigned int index){

	my $index = shift; #C++ unsigned int

	my $res;							# vector<Residue*> res;
	my $chains_dot_size = $#{$structure->{'chains'}} + 1;
	if ($index < $chains_dot_size) {				# if (index < chains.size()){
		# for (unsigned int i = chains[index].min; i <= chains[index].max; i++){
		for ( my $i = $structure->{'chains'}->[$index]->{'min'}; $i <= $structure->{'chains'}->[$index]->{'max'}; $i++ ) {
			my $push_back_index = $#{$res} + 1;
			$res->[$push_back_index] = $structure->{'allResidues'}->[$i]; # res.push_back(&allResidues[i]);
		}
		return $res;						# return res;
	}
	print "ERROR: Internal error: Chain index not found!\n";	# AfxMessageBox("Internal error: Chain index not found!");
	return $res;							# return res;
}



sub Structure_getResChain_Range { #C++ vector<Residue*> Structure::getResChain(Range r){

	my $r = shift;
	my @res_array;
	my $res = \@res_array; #C++ vector<Residue*> res;

	for ( my $i = $r->{'min'}; $i <= $r->{'max'}; $i++ ) {		# for( unsigned int i = r.min; i <= r.max; i++){
		my $push_back_index = $#{$res} + 1;
		$res->[$push_back_index] = $structure->{'allResidues'}->[$i]; # res.push_back(&allResidues[i]);
	}
	return $res;							# return res;
}



sub Structure_getTMDomains { #C++ vector<Range> Structure::getTMDomains(){

	$tmDomains = (); # initialize this data structure for each new sequence

	my $allResidues_dot_size = $#{$structure->{'allResidues'}};
	my $i = 0;
	while ($i < $allResidues_dot_size) {						# for (unsigned int i = 0; i < allResidues.size(); i++){
		if ( $structure->{'allResidues'}->[$i]->{'isTransMem'} ) {		# if (allResidues[i].isTransMem){

			my %tm_hash;
			my $tm = \%tm_hash;						# Range tm;
			$tm->{'chainID'} = $structure->{'allResidues'}->[$i]->{'chainID'}; # tm.ChainID = allResidues[i].chainID;
			$tm->{'min'} = $structure->{'allResidues'}->[$i]->{'resSeq'}; 	# tm.min = allResidues[i].resSeq;

			my $at_end_of_a_tmh = 0;
			while ($at_end_of_a_tmh == 0) {

				$i++;
				my $i_plus_1 = $i + 1;

				if ($i_plus_1 >= $allResidues_dot_size) {
					$at_end_of_a_tmh = 1;
				} elsif ($structure->{'allResidues'}->[$i_plus_1]->{'isTransMem'} == 0) {
					$at_end_of_a_tmh = 1;
				} else {
					# for valpred1, if 2 predicted TMS are stuck together, separate them, don't predict 1 long TMS
					my $next_residue_is_a_tmDomain_start = 0;
					for ( my $a = 0; $a < @tmDomain_start; $a++ ) { # for valpred1, if 2 predicted TMS are stuck together, separate them, don't predict 1 long TMS
						if ($tmDomain_start[$a] == $i_plus_1) {
							$next_residue_is_a_tmDomain_start = 1;
						}
					}
					if ($next_residue_is_a_tmDomain_start == 1) {
						$at_end_of_a_tmh = 1;
					}
				}
			}

			$tm->{'max'} = $structure->{'allResidues'}->[$i]->{'resSeq'}; # tm.max = allResidues[i-1].resSeq;
			my $push_back_index = $#{$tmDomains} + 1;
			$tmDomains->[$push_back_index] = $tm;	
			# tmDomains.push_back(tm);

			$i++;
		} else {
			$i++;
		}
	}
	return;										# return tmDomains;
}



sub Structure_new {

	my @s_SheetInfo_array;
	$structure->{ 'SheetInfo' } = \@s_SheetInfo_array; 			#C++ vector<Range> SheetInfo;
	my @s_HelixInfo_array;
	$structure->{ 'HelixInfo' } = \@s_HelixInfo_array;	 		#C++ vector<Range> HelixInfo;
	my @s_allAtoms_array;
	$structure->{ 'allAtoms' } = \@s_allAtoms_array; 			#C++ vector<Atom> allAtoms;
	my @s_allResidues_array;
	$structure->{ 'allResidues' } = \@s_allResidues_array; 			#C++ vector<Residue> allResidues;
	my @s_chains_array;
	$structure->{ 'chains' } = \@s_chains_array; 				#C++ vector<Range> chains;
	my @s_hetAtoms_array;
	$structure->{ 'hetAtoms' } = \@s_hetAtoms_array; 			#C++ vector<HetAtom> hetAtoms;

	# // Construction/Destruction

	$structure->{'ASA'} = 0;						# ASA = 0.0f; # //Total SASA of structure
	$structure->{'aveResScore'} = 0;					# aveResScore = 0.0f;
	$structure->{'scoreTotal'} = 0;						# float scoreTotal;

	my @xyzExtremities_array;
	$structure->{'xyzExtremities'} = \@xyzExtremities_array;

	for (my $i = 0; $i < 3; $i++) {					# for (int i = 0; i < 3; i++){
									# struct vec2f{
									# float f1;
									# float f2;
		my $v;							# vec2f v;
		$v->{'f1'} = 0;						# v.f1 = 0.0f; # //minimum value
		$v->{'f2'} = 0;						# v.f2 = 0.0f; # //maximum value
									# vector<vec2f> xyzExtremities;	//the maximum and minimum xyz values
		my $push_back_index = $#{$structure->{'xyzExtremities'}} + 1;
		$structure->{'xyzExtremities'}->[$push_back_index] = $v;	# xyzExtremities.push_back(v);
	}

	return;
}



# Perl trim function to remove whitespace from the start and end of the string
sub trim {
        my $string = shift;
        $string =~ s/^\s+//;
        $string =~ s/\s+$//;
        return $string;
}



# //==========================================================================
# // Vec3.hpp : 3d vector class template. Works for any integer or real type.
# //==========================================================================
# // from GLUV
# // http://www.cs.unc.edu/~walk/software/glvu/html/vec4f_8hpp-source.html
# //Author: 
# //Kenny Hoff <hoff at cs.unc.edu> , Bill Baxter <baxter at billbaxter.com> ,
# //Wes Hunt <hunt at cs.unc.edu> , Mark Harris <harrism at cs.unc.edu> 



sub untaint {
	my $path = shift;
	# user_error("path cannot contain metacharacters") if $path=~/[\n|<>&!;\'\"]/;

	# extract only legal alphanumerics from the filename
	# $path =~ m!(/[a-zA-Z/0-9._~\-]+)!;
	# $path =~ /([\w"-"".""_"]+)/;
	$path =~ /([a-zA-Z0-9._~\-]+)/;
	if (defined($1)) {
		$path = $1;
	} else {
		$path = '';
	}

	# user_error("path cannot contain relative directories") 
	#	if $path=~m!\.\.!;

	return $path;
}



sub Vec3_add_Vec3 { #C++ inline Vec3& operator += (const Vec3& A) 

	my $vec3_1 = shift;
	my $vec3_2 = shift;
	my $vec3_new;
	# { x+=A.x; y+=A.y; z+=A.z;
	$vec3_new->{'x'} = $vec3_1->{'x'} + $vec3_2->{'x'};
	$vec3_new->{'y'} = $vec3_1->{'y'} + $vec3_2->{'y'};
	$vec3_new->{'z'} = $vec3_1->{'z'} + $vec3_2->{'z'};

	return $vec3_new; #C++ Vec3f or Vec3d	# return *this; }
}



sub Vec3_cross { #C++ inline Vec3 cross (const Vec3& A) const

	my $vec3_1 = shift;
	my $vec3_2 = shift;
	my $vec3_new;
	# Vec3 CrossProd(y*A.z-z*A.y, z*A.x-x*A.z, x*A.y-y*A.x);
	$vec3_new->{'x'} = ($vec3_1->{'y'} * $vec3_2->{'z'}) - ($vec3_1->{'z'} * $vec3_2->{'y'});
	$vec3_new->{'y'} = ($vec3_1->{'z'} * $vec3_2->{'x'}) - ($vec3_1->{'x'} * $vec3_2->{'z'});
	$vec3_new->{'z'} = ($vec3_1->{'x'} * $vec3_2->{'y'}) - ($vec3_1->{'y'} * $vec3_2->{'x'});

	return $vec3_new; #C++ Vec3f or Vec3d	# return(CrossProd); };
}



sub Vec3_dot_Vec3 { #C++ inline Type dot (const Vec3& A) const

	my $vec3_1 = shift;
	my $vec3_2 = shift;
	my $vec3_new;

	# { Type DotProd = x*A.x+y*A.y+z*A.z;
	$vec3_new = ($vec3_1->{'x'} * $vec3_2->{'x'}) + ($vec3_1->{'y'} * $vec3_2->{'y'}) + ($vec3_1->{'z'} * $vec3_2->{'z'});

	return $vec3_new; # return(DotProd); };
}



sub Vec3_length { #C++ inline Type length (void) const                  
    
	my $vec3 = shift;
	# { return ((Type)sqrt(x*x+y*y+z*z)); };
	my $vec3_new;
	$vec3_new = ($vec3->{'x'} * $vec3->{'x'}) + ($vec3->{'y'} * $vec3->{'y'}) + ($vec3->{'z'} * $vec3->{'z'});
	$vec3_new = sqrt($vec3_new);

	return $vec3_new;
}



sub Vec3_lengthSqr { #C++ inline Type lengthSqr (void) const                  
    
	my $vec3 = shift;
	# { return (x*x+y*y+z*z); };
	my $vec3_new;
	$vec3_new = ($vec3->{'x'} * $vec3->{'x'}) + ($vec3->{'y'} * $vec3->{'y'}) + ($vec3->{'z'} * $vec3->{'z'});

	return $vec3_new;
}



sub Vec3_multiply_number { #C++ inline Vec3& operator *= (const Type s) 

	my $vec3 = shift;
	my $s = shift;
	my $vec3_new;
	# { x*=s; y*=s; z*=s;
	$vec3_new->{'x'} = $vec3->{'x'} * $s;
	$vec3_new->{'y'} = $vec3->{'y'} * $s;
	$vec3_new->{'z'} = $vec3->{'z'} * $s;

	return $vec3_new; #C++ Vec3f or Vec3d	# return *this; }
}



# in C++ vec3h :
# static Vec3 ZERO;
# typedef Vec3<float> Vec3f;
# typedef Vec3<double> Vec3d;
# template<class Type> Vec3<Type> Vec3<Type>::ZERO = Vec3<Type>(0,0,0);
# in Perl, apparently all numbers are treated as double, and not as the less precise float,
# and the Perl compiler/interpreter is written in C/C++,
# so all numbers should end up being represented as C++ double.
sub Vec3_new {

	my $vec3_new;
	$vec3_new->{ 'x' } = 0;
	$vec3_new->{ 'y' } = 0;
	$vec3_new->{ 'z' } = 0;

	return $vec3_new; #C++ Vec3f or Vec3d
}



sub Vec3_new_set1 {

	my $xyx = shift;
	my $vec3_new;
	$vec3_new->{ 'x' } = $xyx->{ 'x' };
	$vec3_new->{ 'y' } = $xyx->{ 'y' };
	$vec3_new->{ 'z' } = $xyx->{ 'z' };

	return $vec3_new; #C++ Vec3f or Vec3d
}



sub Vec3_new_set3 {

	my $x = shift;
	my $y = shift;
	my $z = shift;
	my $vec3_new;
	$vec3_new->{ 'x' } = $x;
	$vec3_new->{ 'y' } = $y;
	$vec3_new->{ 'z' } = $z;

	return $vec3_new; #C++ Vec3f or Vec3d
}



sub Vec3_normalize { #C++ inline Vec3& normalize (void)                  
    
	my $vec3 = shift;
	my $vec3_new;
	$vec3_new = $vec3;
	my $L = Vec3_length( $vec3_new );	# { Type L = length();			# //  CALCULATE LENGTH
	if ($L > 0) {				# if (L>0) { x/=L; y/=L; z/=L; }	# //  DIV COMPONENTS BY LENGTH
		$vec3_new->{'x'} = $vec3_new->{'x'} / $L;
		$vec3_new->{'y'} = $vec3_new->{'y'} / $L;
		$vec3_new->{'z'} = $vec3_new->{'z'} / $L;
	}

	return $vec3_new; #C++ Vec3f or Vec3d	# return *this;
}



sub Vec3_set3 {

	my $vec3 = shift;
	my $x = shift;
	my $y = shift;
	my $z = shift;
	my $vec3_new;

	$vec3_new->{ 'x' } = $x;
	$vec3_new->{ 'y' } = $y;
	$vec3_new->{ 'z' } = $z;

	return $vec3_new; #C++ Vec3f or Vec3d
}



sub Vec3_subtract_Vec3 {

	my $vec3_1 = shift;
	my $vec3_2 = shift;
	my $vec3_new;
	$vec3_new->{ 'x' } = $vec3_1->{ 'x' } - $vec3_2->{ 'x' };
	$vec3_new->{ 'y' } = $vec3_1->{ 'y' } - $vec3_2->{ 'y' };
	$vec3_new->{ 'z' } = $vec3_1->{ 'z' } - $vec3_2->{ 'z' };

	return $vec3_new; #C++ Vec3f or Vec3d
}



sub write_output_amphipathic_helices {

	$structure = $analyser->{'st'};

	my $output_line = '';

	my $output_sequence_id = trim($input_sequence_id);
	my $output_fasta_id = trim($input_fasta_id);
	if ($output_fasta_id eq '') {
		$output_fasta_id = '-';
	}
	if ($output_sequence_id eq '') {
		$output_sequence_id = $output_fasta_id;
	}

	$output_line .= $output_sequence_id . ':::::';

	$output_line .= $output_fasta_id . ':::::';

	my %code_3to1;
	$code_3to1{'ALA'} = 'A';
	$code_3to1{'ARG'} = 'R';
	$code_3to1{'ASN'} = 'N';
	$code_3to1{'ASP'} = 'D';
	$code_3to1{'CYS'} = 'C';
	$code_3to1{'GLN'} = 'Q';
	$code_3to1{'GLU'} = 'E';
	$code_3to1{'GLY'} = 'G';
	$code_3to1{'HIS'} = 'H';
	$code_3to1{'ILE'} = 'I';
	$code_3to1{'LEU'} = 'L';
	$code_3to1{'LYS'} = 'K';
	$code_3to1{'MET'} = 'M';
	$code_3to1{'PHE'} = 'F';
	$code_3to1{'PRO'} = 'P';
	$code_3to1{'SER'} = 'S';
	$code_3to1{'THR'} = 'T';
	$code_3to1{'TRP'} = 'W';
	$code_3to1{'TYR'} = 'Y';
	$code_3to1{'VAL'} = 'V';
	$code_3to1{'UNK'} = 'X';

	my $output_residues = '';
	my $s_arrow_allResidues_dot_size = $#{$structure->{'allResidues'}} + 1;
	for ( my $i = 0; $i < $s_arrow_allResidues_dot_size; $i++ ) {
		my $code3 = $structure->{'allResidues'}->[$i]->{'name'};
		my $code1 = $code_3to1{$code3};
		$output_residues .= $code1;
	}
	if ($output_residues eq '') {
		$output_residues = '-';
	}

	$output_line .= $output_residues . ':::::';

	my $output_segments = '';
	get_amphipathic_helices();

	my $amphipathic_helices_dot_size = $#{$amphipathic_helices} + 1;
	for ( my $i = 0; $i < $amphipathic_helices_dot_size; $i++ ) {
		my $amphipathic_helix_start = $amphipathic_helices->[$i]->{'start_residue'};
		my $amphipathic_helix_end = $amphipathic_helices->[$i]->{'end_residue'};
		my $output_segment = $amphipathic_helix_start . '-' . $amphipathic_helix_end . ',';
		$output_segments .= $output_segment;
	}
	if ($output_segments eq '') {
		$output_segments = '-';
	}

	$output_line .= $output_segments . ':::::';

	print AMPHIFILE "$output_line\n";

	close AMPHIFILE;
	open( AMPHIFILE, ">>$amphipathic_output_file") or
		die "Cannot open $amphipathic_output_file for writing : $!\n";
}



sub write_output_segments_and_graphpoints_files {

	$structure = $analyser->{'st'};

	my $output_line = '';

	my $output_sequence_id = trim($input_sequence_id);
	my $output_fasta_id = trim($input_fasta_id);
	if ($output_fasta_id eq '') {
		$output_fasta_id = '-';
	}
	if ($output_sequence_id eq '') {
		$output_sequence_id = $output_fasta_id;
	}

	$output_line .= $output_sequence_id . ':::::';

	$output_line .= $output_fasta_id . ':::::';

	my %code_3to1;
	$code_3to1{'ALA'} = 'A';
	$code_3to1{'ARG'} = 'R';
	$code_3to1{'ASN'} = 'N';
	$code_3to1{'ASP'} = 'D';
	$code_3to1{'CYS'} = 'C';
	$code_3to1{'GLN'} = 'Q';
	$code_3to1{'GLU'} = 'E';
	$code_3to1{'GLY'} = 'G';
	$code_3to1{'HIS'} = 'H';
	$code_3to1{'ILE'} = 'I';
	$code_3to1{'LEU'} = 'L';
	$code_3to1{'LYS'} = 'K';
	$code_3to1{'MET'} = 'M';
	$code_3to1{'PHE'} = 'F';
	$code_3to1{'PRO'} = 'P';
	$code_3to1{'SER'} = 'S';
	$code_3to1{'THR'} = 'T';
	$code_3to1{'TRP'} = 'W';
	$code_3to1{'TYR'} = 'Y';
	$code_3to1{'VAL'} = 'V';
	$code_3to1{'UNK'} = 'X';

	my $output_residues = '';
	my $output_scores_profiles3d = '';
	my $output_scores_repimps = '';
	my $s_arrow_allResidues_dot_size = $#{$structure->{'allResidues'}} + 1;
	for ( my $i = 0; $i < $s_arrow_allResidues_dot_size; $i++ ) {
		my $code3 = $structure->{'allResidues'}->[$i]->{'name'};
		my $code1 = $code_3to1{$code3};
		$output_residues .= $code1;
		$output_scores_profiles3d .= $structure->{'allResidues'}->[$i]->{'scores'}->[0] . ',';
		$output_scores_repimps .= $structure->{'allResidues'}->[$i]->{'scores'}->[1] . ',';
	}
	if ($output_residues eq '') {
		$output_residues = '-';
	}

	$output_line .= $output_residues . ':::::';

	my $output_line_scores = $output_line . 'PROFILES3D;;;;;' . $output_scores_profiles3d . ':::::REPIMPS;;;;;' . $output_scores_repimps . ':::::';
	print SCOREOUTFILE "$output_line_scores\n";
	close SCOREOUTFILE;
	open( SCOREOUTFILE, ">>$scores_output_file") or
		die "Cannot open $scores_output_file for writing : $!\n";

	my $output_segments = '';
	Structure_getTMDomains();					# vector<Range> tmDomains = s->getTMDomains();
	my $tmDomains_dot_size = $#{$tmDomains} + 1;
	for ( my $i = 0; $i < $tmDomains_dot_size; $i++ ) {
		my $tms_start = $tmDomains->[$i]->{'min'} + 1;
		my $tms_end = $tmDomains->[$i]->{'max'} + 1;
		my $output_segment = $tms_start . '-' . $tms_end . ',';
		$output_segments .= $output_segment;
	}
	if ($output_segments eq '') {
		$output_segments = '-';
	}

	$output_line .= $output_segments . ':::::';

	print SEGOUTFILE "$output_line\n";

	close SEGOUTFILE;
	open( SEGOUTFILE, ">>$segment_output_file") or
		die "Cannot open $segment_output_file for writing : $!\n";

	if (defined($input_flag_multiplegraphs)) {
		$output_line .= $multiplegraphs_data;
		print SEGPTSOUTFILE "$output_line\n";
	}
}



