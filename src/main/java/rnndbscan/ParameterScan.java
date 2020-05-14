package rnndbscan;

import info.debatty.java.graphs.Graph;
import info.debatty.java.graphs.Neighbor;
import info.debatty.java.graphs.build.ThreadedBrute;
import org.hipparchus.clustering.Cluster;
import org.hipparchus.clustering.Clusterer;
import org.hipparchus.clustering.DBSCANClusterer;
import org.hipparchus.clustering.DoublePoint;
import org.hipparchus.clustering.distance.DistanceMeasure;
import org.hipparchus.clustering.distance.EuclideanDistance;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.*;

import static java.lang.Math.log;
import static org.hipparchus.util.CombinatoricsUtils.binomialCoefficient;

/***
 * Program for determining the best value of k / minPts and eps for a given clusterer and input file
 * as determined by ARI or NMI as requested
 * Input file must be "original" (coordinates and label) format
 * Unlabeled data MUST NOT be used
 * And label must be a positive integer (Label of 0 is verboten!).
 *
 * Use: java ParameterScan [-v] (RNN | REC | IS | ISB | DBS | OPT) <input_file> (ARI | NMI)
 * Output: Filename \t ARI/NMI \t ARI/NMI value \t Parameter(s)
 *
 * The v option outputs each iterations
 *
 * Primarily invoked via a helper bash script that should automatically generate tables like Tbl 3-5
 ***/

public class ParameterScan {
    private static String[] validAbbrev = {"RNN", "REC", "IS", "ISB", "DBS", "OPT"};

    static int verbose = 0;

    // valid k ranges
    static int kmin = 1;
    static int kmax = 101;

    // valid DBSCAN minPts values (DBSCAN eps values are more complicated)
    static int[] minPtsValues = {1, 5, 10, 20};

    // For storing the raw data
    private static ArrayList<DoublePoint> data = null;

    // Canonical clustering
    private static ArrayList<Cluster<DoublePoint>> canonical = null;

    // global data (here to prevent dumb/spaghetti-code splitting of print statements)
    static int bestK = -1;
    static int bestMinPts = 0;
    static double bestEps = 0.0;
    static double maxScore = 0.0;
    static int numClusters = 0;

    public static void main(final String[] args)  {

        if (args.length < 3) {
            System.err.println("Usage: java ParameterScan [-v] (RNN | REC | IS | ISB | DBS | OPT) <input_file> (ARI | NMI)");
            System.exit(10);
        }

        if (args[0].equals("-v")) {
            verbose = 1;
        }

        String clustAbbrev = args[verbose];
        if (!Arrays.asList(validAbbrev).contains(clustAbbrev)) {
            System.err.println("Invalid clusterer abbreviation given");
            System.exit(20);
        }

        String comparisonMethod = args[2+verbose];
        if (!(comparisonMethod.equalsIgnoreCase("ARI") || comparisonMethod.equalsIgnoreCase("NMI"))) {
            System.err.println("Invalid comparison class given");
            System.exit(30);
        }

        data = new ArrayList<>();
        canonical = new ArrayList<>();

        try {
            Scanner in = new Scanner(new BufferedReader(new FileReader(args[1+verbose])));

            // Number of fields is num tabs+1, number of coordinates is num of tabs
            String firstLine = in.nextLine();
            int dimension = Math.toIntExact(firstLine.codePoints().filter(ch -> ch == '\t').count());

            // move back to beginning and begin reading in data
            in.reset();

            while (in.hasNextFloat()) {
                double[] point = new double[dimension];
                for (int j=0; j<dimension; j++) {
                    point[j]=in.nextFloat();
                }

                // add to raw data
                DoublePoint dp = new DoublePoint(point);
                data.add(dp);

                // add to cluster of given number (label) [ -1 due to 0 indexing]
                // Usually an int, but if it is "noise" then don't add it to any cluster
                String finalField = in.next();
                if (!finalField.equalsIgnoreCase("noise")) {
                    int label = Integer.parseInt(finalField) - 1;
                    while (canonical.size() <= label) {
                        canonical.add(null);
                    }
                    Cluster<DoublePoint> cluster = canonical.get(label);
                    if (cluster == null) {
                        cluster = new Cluster<>();
                        canonical.set(label, cluster);
                    }
                    cluster.addPoint(dp);
                }
            }
            in.close();
        }
        catch (FileNotFoundException ex) {
            System.err.println("File does not exist");
            System.exit(40);
        }

        // Sort the data for ease of comparison, here for efficiency
        canonical.sort(Comparator.comparingInt(left -> left.getPoints().size()));

        if (clustAbbrev.equalsIgnoreCase("DBS") || clustAbbrev.equalsIgnoreCase("OPT")) {
            minepsLoop(clustAbbrev, comparisonMethod);
            System.out.println(args[1+verbose]+"\t"+comparisonMethod+"\t"+maxScore+"\tclu\t"+numClusters+"\teps="+bestEps+"\tminPts="+bestMinPts);
        }
        else {
            kLoop(clustAbbrev, comparisonMethod);
            System.out.println(args[1+verbose]+"\t"+comparisonMethod+"\t"+maxScore+"\tclu\t"+numClusters+"\tk="+bestK);
        }
        System.exit(0);
    }

    // Main loop for min_pts/eps clusterers (DBS and OPT)
    static void minepsLoop(String type, String comparison) {
        double curScore = 0.0;
        boolean ari = comparison.equalsIgnoreCase("ARI");

        for (int minPts : minPtsValues) {
            DistanceMeasure distance = new EuclideanDistance();
            ThreadedBrute<DoublePoint> builder = new ThreadedBrute<>();
            builder.setK(minPts);
            builder.setSimilarity((left,right) ->
                1.0 / (distance.compute(left.getPoint(), right.getPoint())));
            Graph<DoublePoint> graph = builder.computeGraph(data);

            TreeSet<Double> epsSet = new TreeSet<>();
            for (DoublePoint node : graph.getNodes()) {
                for (Object neigh : graph.getNeighbors(node)) {
                    epsSet.add(Math.abs(distance.compute(((Neighbor<DoublePoint>)neigh).getNode().getPoint(), node.getPoint())));
                }
            }

            for (double eps : epsSet) {
                DBSCANClusterer<DoublePoint> dbs = new DBSCANClusterer<>(eps, minPts);
                List<Cluster<DoublePoint>> clusters = dbs.cluster(data);
                clusters.sort(Comparator.comparingInt(left -> left.getPoints().size()));
                matchClusters(clusters, canonical);
                curScore = (ari ? ARI(clusters, canonical) : NMI(canonical, clusters));

                if (verbose!=0) {
                    System.out.println(Double.toString(eps)+"\t"+Integer.toString(minPts) + "\t" + clusters.size()+"\t"+curScore);
                }

                if (curScore > maxScore) {
                    maxScore = curScore;
                    bestMinPts = minPts;
                    bestEps = eps;
                    numClusters = clusters.size();
                }
            }
        }
    }

    // Main loop for k parameters
    static void kLoop(String type, String comparison) {
        double curScore = 0.0;
        boolean ari = comparison.equalsIgnoreCase("ARI");

        for (int k=kmin; k<kmax; k++) {
            Clusterer<DoublePoint> clusterer = factory(type, k);
            List<Cluster<DoublePoint>> clusters = (List<Cluster<DoublePoint>>) clusterer.cluster(data);

            clusters.sort(Comparator.comparingInt(left -> left.getPoints().size()));
            matchClusters(clusters, canonical);
            curScore = (ari ? ARI(clusters, canonical) : NMI(canonical, clusters));

            if (verbose!=0) {
                System.out.println(Integer.toString(k) + "\t" + clusters.size()+"\t"+curScore);
            }

            if (curScore > maxScore) {
                maxScore = curScore;
                bestK = k;
                numClusters = clusters.size();
            }
        }
    }

    // match up, position-wise, each canonical cluster with the most similar cluster (helps with correct ARI calc)
    // Basically insertion sort (possible to replace with sort call?)
    private static void matchClusters(List<Cluster<DoublePoint>> clusters, List<Cluster<DoublePoint>> canon) {
        for (int i=0; i<canon.size(); i++) {
            int bestIndex=-1;
            int bestMatches=0;
            for (int j=i; j<clusters.size(); j++) {
                ArrayList<DoublePoint> common = new ArrayList<>(canon.get(i).getPoints());
                common.retainAll(clusters.get(j).getPoints());
                int size = common.size();
                if (size > bestMatches) {
                    bestIndex=j;
                    bestMatches=size;
                }
            }

            // basic swap routine
            if (bestIndex>-1) {
                Cluster<DoublePoint> original = clusters.get(i);
                clusters.set(i, clusters.get(bestIndex));
                clusters.set(bestIndex, original);
            }
        }
    }

    static Clusterer<DoublePoint> factory(String type, int k) {
        switch (type.toUpperCase()) {
            case "RNN": return new RNNDBSCANClusterer<>(k);
            case "REC": return new RECORDClusterer<>(k);
            //case "IS" : return new ISDBSCANClusterer<>(k);
            case "ISB" : return new ISBDBSCANClusterer<>(k);
            default : throw new IllegalArgumentException();
        }
    }

    /** Adjusted Rand Index (Formula thanks to Dave Tang: https://davetang.org/muse/2017/09/21/adjusted-rand-index)
     *                      (Explanation of n parameter: http://faculty.washington.edu/kayee/pca/supp.pdf)
     *  WARNING: Highly inaccurate, depends on inputs being sorted by size in order to align the matrix
     *          correctly.  This fails if the sizes are off, or multiple clusters have same size.
     *
     *  Possible Solution: Match position-wise each canonical cluster with its most similar cluster
     *  Before sorting range of results was mostly -0.3 to -4, now it caps at 1.0 but can still dip below -1.0
     *  Correct range is -1.0 to 1.0
     *
     *  Yields in ballpark results for RNN, but too high for REC
     */
    private static double ARI(List<? extends Cluster<DoublePoint>> canonical,
                              List<? extends Cluster<DoublePoint>> clusters) {
        int cols = canonical.size();
        int rows = clusters.size();
        int matrix[][] = new int[rows+1][cols+1];

        for (int i=0; i<rows; i++) {
            for (int j=0; j<cols; j++) {
                for (DoublePoint rp : clusters.get(i).getPoints()) {
                    for (DoublePoint cp : canonical.get(j).getPoints()) {
                        if (rp.equals(cp)) {
                            matrix[i][j]++;
                            matrix[i][cols]++;
                            matrix[rows][j]++;
                        }
                    }
                }
            }
        }

        long diagonalPart=0;
        for (int i=0; i<Integer.min(rows,cols); i++) {
            if (matrix[i][i] >= 2)
                diagonalPart += binomialCoefficient(matrix[i][i], 2);
        }

        int rowsN = 0;
        long rowsPart=0;
        for (int i=0; i<rows; i++) {
            if (matrix[i][cols] >= 2)
                rowsPart += binomialCoefficient(matrix[i][cols], 2);
            rowsN += matrix[i][cols];
        }

        int colsN = 0;
        long colsPart=0;
        for (int i=0; i<cols; i++) {
            if (matrix[rows][i] >= 2)
                colsPart += binomialCoefficient(matrix[rows][i], 2);
            colsN += matrix[rows][i];
        }

        if (rowsN!=colsN)
            System.err.println("IMPOSSIBLE!!!");

        // n choose 2 part
        double nc2 = (rowsN>2) ? binomialCoefficient(rowsN, 2) : binomialCoefficient(2, rowsN);
        double common = (rowsPart*colsPart)/nc2;
        double top = diagonalPart-common;
        double bot = (0.5*(rowsPart+colsPart))-common;

        return top / bot;
    }

    /** Normalized Mutual Information (Broken Implementation)
     *  Sources of formula: https://course.ccs.neu.edu/cs6140sp15/7_locality_cluster/Assignment-6/NMI.pdf
     *                      https://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-clustering-1.html
     *
     *  This implementation is tested in NMITest, it passes both simple tests (lifted from ccs.neu.edu)
     *  It fails spectacularly on actual data, however it stays in the correct range.
     *
     *  A previous version produced closer results but behavior was either wrong (continually growing)
     *  or was outside of the proper range (negative values, and greater than 1.0) no matter the tweaks made.
     */
    static double NMI(List<Cluster<DoublePoint>> canon, List<Cluster<DoublePoint>> clust) {

        int    canonTotal = canon.stream().map(Cluster::getPoints).mapToInt(List::size).sum();
        int    clustTotal = clust.stream().map(Cluster::getPoints).mapToInt(List::size).sum();
        int    total = clustTotal + canonTotal;

        double HY = 0.0;
        for (int i=0; i<Integer.max(clust.size(), canon.size()); i++) {
            double classTotal = 0;
            if (i < canon.size()) {
               classTotal += canon.get(i).getPoints().size();
            }
            if (i < clust.size()) {
               classTotal += clust.get(i).getPoints().size();
            }
            double Pcj = classTotal / total;
            if (Pcj > 0) {
                HY += Pcj * (log(Pcj) / log(2));
            }
        }
        HY = -HY;

        double HC1frac = (double)canonTotal/total;
        double HC2frac = (double)clustTotal/total;
        double HC = -(HC1frac)*(log(HC1frac)/log(2)) - (HC2frac)*(log(HC2frac)/log(2));

        double HYC1 = 0.0;
        for (Cluster<DoublePoint> wk : clust) {
            double Pwk = (double)wk.getPoints().size() / clustTotal;
            HYC1 += Pwk * (log(Pwk)/log(2));
        }
        HYC1 = -((double)clustTotal/total)*HYC1;

        double HYC2 = 0.0;
        for (Cluster<DoublePoint> cj : canon) {
            double Pcj = (double)cj.getPoints().size() / canonTotal;
            HYC2 += Pcj * (log(Pcj)/log(2));
        }
        HYC2 = -((double)canonTotal/total)*HYC2;

        double I = HY - (HYC1 + HYC2);

        return 2*I/(HY+HC);
    }

/* Outdated version of NMI, here in case of an attempt to fix first method
        int    classesTotal = classes.stream().map(Cluster::getPoints).mapToInt(List::size).sum();
        int    clusterTotal = clusters.stream().map(Cluster::getPoints).mapToInt(List::size).sum();
        int    total = clusterTotal + classesTotal;
        //int total = data.size();

        double I = 0.0;
        for (Cluster<DoublePoint> wk : clusters) {
            for (Cluster<DoublePoint> cj : classes) {
                ArrayList<DoublePoint> intersection = new ArrayList<>(wk.getPoints());
                intersection.retainAll(cj.getPoints());
                double frac = (double)intersection.size() / (2*total);
                if (frac>0) {
                    double Pwk = (double) wk.getPoints().size() / total;
                    double Pcj = (double) cj.getPoints().size() / total;
                    I += frac * log(frac / (Pwk * Pcj));
                }
            }
        }

        double HO = 0.0;
        for (Cluster<DoublePoint> wk : clusters) {
            double Pwk = (double)wk.getPoints().size() / total;
            HO += Pwk * log(Pwk);
        }
        HO = -HO;

        double HC = 0.0;
        for (Cluster<DoublePoint> cj : classes) {
            double Pcj = (double)cj.getPoints().size() / total;
            HC += Pcj * log(Pcj);
        }
        HC = -HC;

        double HC = 0.0;
        for (int i=0; i<Integermin(clusters.size(), classes.size()); i++) {
            Cluster<DoublePoint> cj = classes.get(i);
            Cluster<DoublePoint> wi = clusters.get(i);
            int classtotal = cj.getPoints().size() + wi.getPoints().size();
            double Pcj = (double)cj.getPoints().size() / classtotal;
            HC += Pcj * log(Pcj);
        }
        HC = -HC;
*/
}
