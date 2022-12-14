# CType_Motif
Looking for motifs in the cell position data of the cell types of the Mouse retina. 

This report  illustrates an approach to searching for potential  motifs in a set of retina cell position sequences, generated by random walk,  by identifying significantly over-represented ungapped words of fixed length. The discussion is based on the case study presented in Chapter 10 of "Cristianini, and Hahn, Introduction to Computational Genomics: A Case Studies Approach",2007,  and the example developed based on this case study in  [Matlab documentation](https://de.mathworks.com/help/bioinfo/ug/identifying-over-represented-regulatory-motifs.html).


We aim to evaluate possible motifs in the real cell position data in contrast to the shuffled one. The sequences are generated with a lazy random walk with laziness parameter(the probability that the random walk stays at the same node) set to $0.15$. We run the random walk process for $100000$ iterations to make sure that the mixing time of the walk has passed and the probability of hitting each node is close to the steady state solution. The graph itself is generated by Delaunay triangle algorithm. Long edges between far nodes are discarded, and the remaing edges  are assumed to be uniform.
