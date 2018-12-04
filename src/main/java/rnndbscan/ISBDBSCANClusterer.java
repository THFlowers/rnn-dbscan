package rnndbscan;

import info.debatty.java.graphs.*;
import info.debatty.java.graphs.build.NNDescent;
import org.hipparchus.clustering.*;
import org.hipparchus.clustering.distance.DistanceMeasure;
import org.hipparchus.clustering.distance.EuclideanDistance;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;

import java.util.*;

public class ISBDBSCANClusterer<T extends Clusterable> extends Clusterer<T> implements SimilarityInterface<T> {
    private int k;
    private static final int UNCLASSIFIED = -1, NOISE = -2;

    public ISBDBSCANClusterer(int k) {
        super(new EuclideanDistance());
        this.k = k;
    }

    public ISBDBSCANClusterer(int k, DistanceMeasure measure) {
        super(measure);
        this.k = k;
    }

    public int getK() {
        return k;
    }

    public double similarity(T val1, T val2) {
        return 1.0 / (getDistanceMeasure().compute(val1.getPoint(), val2.getPoint()));
    }

    // Shared Data between IS-DBSCAN subroutines, here for convenience
    private HashMap<T, Integer> assign;
    private ArrayList<T> nodes;
    private Graph<T> NN, rNN;

    @Override
    public List<Cluster<T>> cluster(Collection<T> points) throws MathIllegalArgumentException, MathIllegalStateException {
        ArrayList<Cluster<T>> clusters = new ArrayList<>();

        assign = new HashMap<>(points.size());
        for (Object point : points) {
            assign.put((T) point, UNCLASSIFIED);
        }
        nodes = new ArrayList<T>(points.size());
        nodes.addAll(points);

        NNDescent<T> builder = new NNDescent<>();
        builder.setRho(0.1);
        builder.setDelta(0.001);
        builder.setK(k);
        builder.setMaxIterations(15);

        builder.setSimilarity(this::similarity);

        NN = builder.computeGraph(nodes);
        rNN = transposeGraph(NN);

        Cluster<T> cluster;
        for (T node : nodes) {
            if (assign.get(node)==UNCLASSIFIED) {
                cluster = ExpandCluster(node, clusters.size(), assign);
                if (cluster!=null && cluster.getPoints().size() > 0)
                    clusters.add(cluster);
            }
        }

        // Borderpoint handling, line 12 and onward in ISB-DBSCAN paper
        for (T x : nodes) {
            int clustnum = NOISE;
            double mindist = Double.POSITIVE_INFINITY;
            if (assign.get(x)==UNCLASSIFIED) {
                for (T p : ISk(x)) {
                    if (assign.get(p)!=NOISE) {
                        double distance = distance(x,p);
                        if (distance < mindist) {
                            mindist = distance;
                            clustnum = assign.get(p);
                        }
                    }
                }
            }
            assign.put(x, clustnum);
            if (clustnum >= 0)
                clusters.get(clustnum).addPoint(x);
        }

        // clear shared data, let garbage collector clean up memory
        assign = null;
        nodes = null;
        NN = null;
        rNN = null;

        return clusters;
    }

    // transposes the graph, returns a graph (performs some reverse engineering to filter to k)
    private Graph<T> transposeGraph(Graph<T> g) {
        Graph<T> newGraph = new Graph<>(g.getNodes().size());
        //newGraph.setK(k);
        newGraph.setSimilarity(g.getSimilarity());

        HashMap<T, NeighborList> orig = g.getHashMap();

        for (T key : orig.keySet()) {
            for (Neighbor<T> neighbor : orig.get(key)) {
                T node = neighbor.getNode();
                NeighborList nl;
                if (!newGraph.containsKey(node)) {
                    // need to supply a size here, but k seems too restrictive. Assuming forall x in X revNN(x) <= k^2
                    nl = new NeighborList(k * k);
                    newGraph.put(node, nl);
                } else {
                    nl = newGraph.getNeighbors(node);
                }

                // note the order, should be transpose of original calculation
                nl.add(new Neighbor(key, similarity(node, key)));
            }
        }

        return newGraph;
    }

    private Cluster<T> ExpandCluster(T x, int clusterNum, HashMap<T, Integer> assign) {
        Cluster<T> cluster = new Cluster<>();
        ArrayDeque<T> seeds = new ArrayDeque<>(ISk(x));
        if (seeds.size() > (2.0/3.0 * k)) {
            assign.put(x, clusterNum);
            cluster.addPoint(x);
        }
        else
            return null;

        while (!seeds.isEmpty()) {
            T y=seeds.removeFirst();
            ArrayList<T> n = ISk(y);
            if (n.size() > (2.0/3.0 * k)) {
                assign.put(y, clusterNum);
                cluster.addPoint(y);
                for (T z : ISk(y)) {
                    int status = assign.get(z);
                    if (status==UNCLASSIFIED || status==NOISE) {
                        if (!seeds.contains(z)) {
                            seeds.addLast(z);
                        }
                    }
                }
            }
        }
        return cluster;
    }

    ArrayList<T> ISk(T p) {
        ArrayList<T> ret = new ArrayList<>();
        NeighborList forward = NN.getNeighbors(p);
        NeighborList reverse = rNN.getNeighbors(p);

        if (forward!=null && reverse!=null) {
            for (Neighbor n : forward) {
                T node = (T) n.getNode();
                if (node != null && reverse.contains(n)) {
                    ret.add(node);
                }
            }
        }

        return ret;
    }
}