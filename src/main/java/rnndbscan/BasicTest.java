package rnndbscan;

import org.hipparchus.clustering.*;

import javax.swing.*;
import java.io.*;
import java.util.ArrayList;
import java.util.InputMismatchException;
import java.util.List;
import java.util.Scanner;

/***
 * Simple executable for testing implementation of clusterers
 * Change rdb type / right hand side to change Class being tested
 *
 * INPUT DATA MUST BE TWO DIMENSIONAL (2D) FOR GUI TO WORK
 *
 * use: java BasicTest <input_file> [eps minPts k]
 * eps and minPts for dbscan, k for other clusterer
 *
 * Output:
 * Number of clusters in DBSCAN, number in other, and gui displayed images of both
 *
 * File format should be same as t4.8k dataset
 * That is, first line is number of  entries and dimensionality of data
 *
 * Use script convert.sh in resources to convert the other main format
 * (coordinates and label) into this format
 ***/

public class BasicTest {
    public static void main(final String[] args)  {

        double eps = 18.0;
        int minPts = 75;
        int k = 40;

        if (args.length < 1) {
            System.err.println("Arg1 should be text file containing set to cluster");
            System.exit(20);
        }
        if (args.length > 1) {
            eps = Double.parseDouble(args[1]);
            minPts = Integer.parseInt(args[2]);
            k = Integer.parseInt(args[3]);
        }

        DBSCANClusterer<DoublePoint> db = new DBSCANClusterer<>(eps, minPts);
        RNNDBSCANClusterer<DoublePoint> rdb = new RNNDBSCANClusterer<>(k);

        ArrayList<DoublePoint> data = null;
        try {
            Scanner in = new Scanner(new BufferedReader(new FileReader(args[0])));
            int rows=in.nextInt();
            int cols=in.nextInt();

            data = new ArrayList<>(rows);
            for (int i=0; i<rows; i++) {

                // in case of empty last line(s)
                if (!in.hasNextFloat())
                    break;

                double[] point = new double[cols];
                for (int j=0; j<cols; j++) {
                    point[j]=in.nextFloat();
                }
                data.add(new DoublePoint(point));
            }
            in.close();
        }
        catch (FileNotFoundException ex) {
            System.err.print("File not found: ");
            System.err.println(args[0]);
            System.exit(30);
        }
        catch (InputMismatchException ex) {
            System.err.println("Found non-floating point data");
            System.exit(30);
        }

        final List<Cluster<DoublePoint>> clusters = db.cluster(data);
        System.out.println(clusters.size());

        ArrayList<DoublePoint> finalData = (ArrayList<DoublePoint>)data.clone();
        SwingUtilities.invokeLater(() -> {
            xyScatterPlot dbplot = new xyScatterPlot("DBSCAN", clusters, finalData);
            dbplot.setSize(800, 400);
            dbplot.setLocationRelativeTo(null);
            dbplot.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
            dbplot.setVisible(true);
        });

        final List<Cluster<DoublePoint>> clusters2 = rdb.cluster(data);
        System.out.println(clusters2.size());

        ArrayList<DoublePoint> finalData2 = (ArrayList<DoublePoint>)data.clone();
        SwingUtilities.invokeLater(() -> {
            xyScatterPlot rnnplot = new xyScatterPlot("RNN-DBSCAN", clusters2, finalData2);
            rnnplot.setSize(800, 400);
            rnnplot.setLocationRelativeTo(null);
            rnnplot.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
            rnnplot.setVisible(true);
        });
    }
}
