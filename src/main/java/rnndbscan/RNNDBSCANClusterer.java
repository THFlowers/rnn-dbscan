package rnndbscan;


import info.debatty.java.graphs.*;
import info.debatty.java.graphs.build.NNDescent;
import info.debatty.java.graphs.build.ThreadedBrute;
import org.hipparchus.clustering.*;
import org.hipparchus.clustering.distance.DistanceMeasure;
import org.hipparchus.clustering.distance.EuclideanDistance;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import java.util.*;

public class RNNDBSCANClusterer<T extends Clusterable> extends Clusterer<T> implements SimilarityInterface<T> {
    private int k;
    private static final int UNCLASSIFIED=-1, NOISE=-2;

    public RNNDBSCANClusterer(int k) {
        super(new EuclideanDistance());
        this.k = k;
    }

    public RNNDBSCANClusterer(int k, DistanceMeasure measure) {
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
    private Map<Cluster, Double> denCache;

    @Override
    public List<Cluster<T>> cluster(Collection<T> points) throws MathIllegalArgumentException, MathIllegalStateException {
        ArrayList<Cluster<T>> clusters = new ArrayList<>();
        denCache = new HashMap<>();

        assign = new HashMap<>(points.size());
        for (Object point : points) {
            assign.put((T)point, UNCLASSIFIED);
        }
        nodes = new ArrayList<T>(points.size());
        nodes.addAll(points);

        // Horribly inefficient, esp. in comparison to DBSCAN but results are pretty good
        // Use this for testing accuracy of implementation
        ThreadedBrute<T> builder = new ThreadedBrute<>();
        builder.setK(k);

        // Wickedly fast as long as maxIterations*k is low,
        // Use this for testing performance after brute is correct

        // Paper recommends rho=0.1, delta=0.001, and k=100, implies iterations of <= 12
        // too low k and results become crap (<25), too high maxIterations (>20) and time explodes even with good low k

        /*
        NNDescent<T> builder = new NNDescent<>();
        builder.setRho(0.1);
        builder.setDelta(0.001);
        builder.setK(k);
        builder.setMaxIterations(15); // 10, a little slower than DBSCAN any lower and wrong, any higher time explodes
        */

        builder.setSimilarity(this);

        NN = builder.computeGraph(nodes);
        rNN = transposeGraph(NN);

        Cluster<T> cluster;
        for (T node : nodes) {
            if (assign.get(node)==UNCLASSIFIED) {
                cluster = ExpandCluster(node, clusters.size(), assign);
                if (cluster!=null && cluster.getPoints().size() > 0) {
                    clusters.add(cluster);
                }
            }
        }

        ExpandClusters(nodes, assign, clusters);

        // clear shared data, let garbage collector clean up memory
        assign = null;
        nodes  = null;
        NN = null;
        rNN = null;
        denCache = null;

        return clusters;
    }

    private void ExpandClusters(ArrayList<T> nodes, HashMap<T, Integer> assign, List<Cluster<T>> clusters) {
        for (T x : nodes) {
            if (assign.get(x)==NOISE) {
                NeighborList neighbors = NN.getNeighbors(x);
                int mincluster=NOISE;
                double mindist = Double.POSITIVE_INFINITY;
                for (Neighbor neigh : neighbors) {
                    T n = (T)neigh.getNode();
                    int clusterNum = assign.get(n);
                    double d=distance(x,n);

                    // Hit a NOISE or
                    if (clusterNum < 0) {
                        continue;
                    }
                    if (clusterNum >= clusters.size()) {
                        System.err.println("Impossible cluster encountered");
                        throw new IllegalStateException();
                    }

                    Cluster<T> cluster = clusters.get(clusterNum);
                    if ((rNN.getNeighbors(n).size() >= k) && (d <= den(cluster)) && (d < mindist)) {
                        mincluster = clusterNum;
                        mindist = d;
                    }

                    if (mincluster>=0) {
                        assign.put(x, mincluster);
                        clusters.get(mincluster).addPoint(x);
                    }
                }
            }
        }
    }

    private double den(Cluster<T> cluster) {
        return denCache.computeIfAbsent(cluster, c -> {
           double max = 0;
           for (Object x : c.getPoints()) {
               for (Object y : c.getPoints()) {
                   int xsz = rNN.getNeighbors((T)x).size();
                   int ysz = rNN.getNeighbors((T)y).size();

                   if (xsz >= k && ysz >= k && directlyDensityReachable((T)x,(T)y)) {
                       double dist = distance((T)x, (T)y);
                       if (dist > max)
                           max = dist;
                   }
               }
           }
           return max;
        });
    }

    // Returns if x is reachable from y (might not be symmetric)
    private boolean directlyDensityReachable(T x, T y) {
         return NN.getNeighbors(y).containsNode(x) && rNN.getNeighbors(y).size() >= k;
    }

    public double similarity(T val1, T val2) {
            return 1.0 / (getDistanceMeasure().compute(val1.getPoint(), val2.getPoint()));
    }

    private Cluster<T> ExpandCluster(T x, int clusterNum, HashMap<T, Integer> assign) {
        NeighborList revNeigh = rNN.getNeighbors(x);
        if (revNeigh == null) {
            return null;
        }
        if (revNeigh.size() < k) {
            assign.put(x, NOISE);
            return null;
        }

        Cluster<T>    cluster = new Cluster<>();
        ArrayList<T>  nbhood  = neighborhood(x);
        ArrayDeque<T> seeds   = new ArrayDeque<>();

        assign.put(x,clusterNum);
        cluster.addPoint(x);
        for (T n : nbhood) {
            seeds.add(n);
            assign.put(n, clusterNum);
            cluster.addPoint(n);
        }

        while (seeds.size() != 0) {
            T y = seeds.remove();
            NeighborList revNeigh2 = rNN.getNeighbors(y);
            if (revNeigh2.size() >= k) {
                ArrayList<T> neighborhood2 = neighborhood(y);
                for (T z : neighborhood2) {
                    if (assign.get(z)==UNCLASSIFIED) {
                        seeds.add(z);
                        cluster.addPoint(z);
                        assign.put(z, clusterNum);
                    }
                    else if (assign.get(z)==NOISE) {
                        cluster.addPoint(z);
                        assign.put(z, clusterNum);
                    }
                }
            }
        }

        return cluster;
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

    // Generate Neighborhood(x,k)
    private ArrayList<T> neighborhood(T x) {
        ArrayList<T> neighbors = new ArrayList<>();

        for (Neighbor n : NN.getNeighbors(x)) {
            neighbors.add((T)n.getNode());
        }

        for (Neighbor y : rNN.getNeighbors(x)) {
            NeighborList revNb = rNN.getNeighbors((T)y.getNode());
            if (revNb != null && revNb.size() > k) {
                neighbors.add((T)y.getNode());
            }
        }

        return neighbors;
    }
}
