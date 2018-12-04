package rnndbscan;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;
import javax.swing.JFrame;
import javax.swing.SwingUtilities;
import javax.swing.WindowConstants;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import org.hipparchus.clustering.Cluster;
import org.hipparchus.clustering.DoublePoint;

import static java.lang.String.format;

public class xyScatterPlot extends JFrame {
    public xyScatterPlot(List<Cluster<DoublePoint>> clusters, ArrayList<DoublePoint> data) {

        super("xy Results");

        XYSeriesCollection dataset = new XYSeriesCollection();
        XYSeries[] series = new XYSeries[clusters.size()];

        for (int i=0; i<clusters.size(); i++) {
            series[i] = new XYSeries(format("Cluster %d", i+1));
            Cluster<DoublePoint> cluster = clusters.get(i);
            for (DoublePoint point : cluster.getPoints()) {
                double[] raw = point.getPoint();
                series[i].add(raw[0], raw[1]);
            }
            data.removeAll(cluster.getPoints());
        }

        // data now only contains un-clustered points
        XYSeries unclustered = new XYSeries("Unclustered");
        for (DoublePoint point : data) {
            double[] raw = point.getPoint();
            unclustered.add(raw[0], raw[1]);
        }

        dataset.addSeries(unclustered);
        for (int i=0; i<clusters.size(); i++) {
            dataset.addSeries(series[i]);
        }

        JFreeChart chart = ChartFactory.createScatterPlot("", "", "", dataset);
        //XYPlot plot = (XYPlot)chart.getPlot();
        ChartPanel panel = new ChartPanel(chart);
        setContentPane(panel);
    }
}
