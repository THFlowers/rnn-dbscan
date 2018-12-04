package rnndbscan;

import info.debatty.java.graphs.*;
import info.debatty.java.graphs.build.NNDescent;
import info.debatty.java.graphs.build.ThreadedBrute;
import info.debatty.java.graphs.build.ThreadedNNDescent;
import org.hipparchus.clustering.*;
import org.hipparchus.clustering.distance.DistanceMeasure;
import org.hipparchus.clustering.distance.EuclideanDistance;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;

import java.util.*;

public class RECORDClusterer<T extends Clusterable> extends Clusterer<T> implements SimilarityInterface<T> {
    private int k;
    private static final int UNCLASSIFIED=-1, NOISE=-2;

    public RECORDClusterer(int k) {
        super(new EuclideanDistance());
        this.k = k;
    }

    public RECORDClusterer(int k, DistanceMeasure measure) {
        super(measure);
        this.k = k;
    }

    public int getK() {
        return k;
    }

    // Shared Data between RNN-DBSCAN subroutines, here for convenience
    private HashMap<T,Integer> assign;
    private ArrayList<T> nodes;
    private Graph<T>   NN, rNN;

    @Override
    public List<Cluster<T>> cluster(Collection<T> points) throws MathIllegalArgumentException, MathIllegalStateException {
        ArrayList<Cluster<T>> clusters = new ArrayList<>();

        // Assign is used in SCC calculation, a bit like how it's used in RNNDBSCAN
        assign = new HashMap<>(points.size());
        for (Object point : points) {
            assign.put((T)point, UNCLASSIFIED);
        }
        nodes = new ArrayList<>(points.size());
        nodes.addAll(points);

        ThreadedBrute<T> builder = new ThreadedBrute<>();
        builder.setK(k);

        /*
        NNDescent<T> builder = new NNDescent<>();
        builder.setRho(0.1);
        builder.setDelta(0.001);
        builder.setK(k);
        builder.setMaxIterations(15);
        */

        builder.setSimilarity(this);

        // Maps to first few stages of record
        NN = builder.computeGraph(nodes);
        rNN = transposeGraph(NN);


        // Now apply any algorithm to get Strongly Connected Components
        // Here we use Kosaraju's algorithm (b/c its the first I found)
        // NOTE: Didn't know about graph.connectedComponents or stronglyConnectedComponents before writing this
        Cluster<T> cluster;
        for (T node : nodes) {
            cluster = visit(node, clusters.size());
            if (cluster != null) { //&& cluster.getPoints().size()>=k) {
                if (cluster.getPoints().size()>=k)
                    clusters.add(cluster);
                else {
                    for (T point : cluster.getPoints()) {
                        assign.put(point, NOISE);
                    }
                }
            }
        }

        // Local outlier incorporation
        int size = clusters.size();
        int numLocal[] = new int[size];
        int nodeAssignment;
        for (T node : nodes) {
            for (int i=0; i<size; i++) {
                numLocal[i]=0;
            }
            NeighborList neighborhood = rNN.getNeighbors(node);
            if (neighborhood!=null) {
                for (Neighbor<T> neigh : neighborhood) {
                    nodeAssignment = assign.get(neigh.getNode());
                    if (nodeAssignment != NOISE) {
                        numLocal[nodeAssignment]++;
                    }
                }

                // dimensionality, here b/c other obvious spots require too much abstraction
                int d = node.getPoint().length;
                for (int i = 0; i < size; i++) {
                    if (numLocal[i] > k / d) {
                        assign.put(node, i);
                        clusters.get(i).addPoint(node);
                        break;
                    }
                }
            }
        }

        // clear shared data, let garbage collector clean up memory
        assign = null;
        nodes  = null;
        NN = null;
        rNN = null;

        return clusters;
    }

    // Here UNCLASSIFIED=unvisited, NOISE or a clusterNum=visited
    Cluster<T> visit(T u, int clustNum, Cluster<T> accum) {
        if (assign.get(u)==UNCLASSIFIED) {
            NeighborList revNeigh=rNN.getNeighbors(u);
            if (revNeigh==null || revNeigh.size() < k) {
                assign.put(u, NOISE);
            }
            else {
                assign.put(u, clustNum);
                accum.addPoint(u);
                for (Neighbor<T> v : rNN.getNeighbors(u)) {
                    visit(v.getNode(), clustNum, accum);
                }
            }
        }
        return accum;
    }

    Cluster<T> visit(T u, int clustNum) {
        if (assign.get(u)==UNCLASSIFIED)
            return visit(u, clustNum, new Cluster<>());
        else
            return null;
    }

    public double similarity(T val1, T val2) {
        return 1.0 / (getDistanceMeasure().compute(val1.getPoint(), val2.getPoint()));
    }

    // transposes the graph, returns a graph (performs some reverse engineering to filter to k)
    private Graph<T> transposeGraph(Graph<T> g) {
        Graph<T> newGraph = new Graph<>(g.getNodes().size());
        //newGraph.setK(k);
        newGraph.setSimilarity(g.getSimilarity());

        HashMap<T, NeighborList> orig = g.getHashMap();

        for (T key : orig.keySet()) {
            for (Neighbor neighbor : orig.get(key)) {
                T node = (T)neighbor.getNode();
                NeighborList nl;
                if (!newGraph.containsKey(node)) {
                    // need to supply a size here, but k seems too restrictive. Assuming forall x in X revNN(x) <= k^2
                    nl = new NeighborList(k*k);
                    newGraph.put(node, nl);
                }
                else {
                    nl = newGraph.getNeighbors(node);
                }

                // note the order, should be transpose of original calculation
                nl.add(new Neighbor(key, similarity(node, key)));
            }
        }

        return newGraph;
    }
}
