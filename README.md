# RNN-DBSCAN

Graduate Databases Fall 2018 Final Project.

Attempt to implement an RNN-DBSCAN based clusterer
and compare the results against other well known clustering algorithms.
In particular we attempted to confirm the results in the original paper describing the algorithm

# What is RNN-DBSCAN?

It is a variant of DBSCAN that uses the k Reverse Nearest Neighbors to estimate local density.
It is structured to be sensitive to changes in density, thereby handling nearby clusters of differing
densities without grouping them together incorrectly or labeling one as noise.

It also has the benefit of having only one parameter, k.
This is apposed to DBSCAN's two, minPts and epsilon (neighborhood size), making parameter space exploration much easier.

See the original paper for details.

## 3rd Party Libraries Used

* Hipparchus Clustering 1.3  (originally Apache Commons, ported)
* Java-graphs 0.41  (reverse nearest neighborhood calculations)
* jfreechart 1.5.0  (display 2d clusters, test correctness and graphics for paper)

## Programs

Excluding shell and gnuplot scripts

### BasicTest

  > java BasicTest <input_file> [eps minPts k]
>
* input-file must be in our format (see datasets below)
* eps and minPts for dbscan
* k for other clusterer (default is rnn-dbscan)

Displays two windows showing the results of clustering the given dataset
with DBSCAN and RNN-DBSCAN with the current parameter configuration.

To use other clusterers the code must be edited.
You may use any Hipparchus clusterer or our RNNDBSCAN, RECORD, or ISBDBSCAN clusterer classes.

### NMITest

Used for debugging our broken implementation of Normalized Mutual Information statistical score.
Takes no parameters, but you may edit the code to supply your own data for testing.

### ParameterScan

  > Usage: java ParameterScan [-v] (RNN | REC | IS | ISB | DBS | OPT) <input_file> (ARI | NMI)
  
  > Output: <input_file> \t ARI/NMI \t ARI/NMI_value \t parameter(s)

Program for determining the best value of k / minPts and eps for a given clusterer (indicated by abbreviation) and input file
(as determined by ARI or NMI as requested)

* Input file must be "original" (coordinates and label) format
* Unlabeled data *MUST NOT* be used.
* Label must be a positive integer (Label of 0 is verboten!).
* -v outputs a line for each iteration of the sweep

It performs a parameter sweep over the parameters accepted by the given clusterer.
The exact ranges are hard-coded.  Reflection is not used, the iteration to perform
is written for each type of clusterer (1 or 2 parameter)

Primarily invoked via a helper bash script that should automatically generate tables like Tbl 3-5


## Datasets

### Origin and Selection
We tried to find all the datasets originally used in the paper.  They appear to be relatively common sets used for testing
clustering algorithms.  We lifted most files from other projects and from https://cs.joensuu.fi/sipu/datasets/.

For the grid example we generated our own (much like the original paper) with the script 'generate_grid.sh'
found in the artificial datasets subdirectory

See Notes.txt in the resources directory for more info

convert.sh is used to convert from the original format of the files from joensuu.fi to the format of the rest of the datasets.
See BasicTest.java comments for more info.

### Our Format
Our format lists number of data points, n, and dimensionality, d, on the first line. 
Then the next n lines list d floating point numbers, the components of each point.


## Original Paper
A. Bryant and K. Cios, "RNN-DBSCAN: A Density-Based Clustering Algorithm Using Reverse Nearest Neighbor Density Estimates," in IEEE Transactions on Knowledge and Data Engineering, vol. 30, no. 6, pp. 1109-1121, 1 June 2018, doi: 10.1109/TKDE.2017.2787640.

