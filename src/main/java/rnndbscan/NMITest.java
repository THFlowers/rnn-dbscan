package rnndbscan;

import org.hipparchus.clustering.Cluster;
import org.hipparchus.clustering.DoublePoint;

import java.util.ArrayList;

public class NMITest {
    public static void main(String[] args) {
        int[][] canonArr1 = {{1,2,3}, {4,5,6}, {7,8,9,10}};
        int[][] clustArr1 = {{1,2}, {3,4,5,6,7,8,9}, {10}};

        double value = ParameterScan.NMI(convert(canonArr1), convert(clustArr1));
        System.out.println(value);

        int[][] canonArr2 = {{1,2,3}, {4,5,6,7,8,9,10}};
        int[][] clustArr2 = {{1,2}, {3,4,5}, {6,7,8,9,10}};

        value = ParameterScan.NMI(convert(canonArr2), convert(clustArr2));
        System.out.println(value);
    }

    public static ArrayList<Cluster<DoublePoint>> convert(int[][] input) {
        ArrayList<Cluster<DoublePoint>> ret = new ArrayList<>();
        for (int[] arr : input) {
            Cluster<DoublePoint> clu = new Cluster<>();
            for (int i : arr) {
                clu.addPoint(new DoublePoint(new int[]{i}));
            }
            ret.add(clu);
        }
        return ret;
    }
}
