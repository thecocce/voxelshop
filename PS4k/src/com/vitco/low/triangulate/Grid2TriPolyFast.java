package com.vitco.low.triangulate;

import com.vitco.util.graphic.G2DUtil;
import com.vitco.util.misc.IntegerTools;
import gnu.trove.set.hash.TIntHashSet;
import org.poly2tri.Poly2Tri;
import org.poly2tri.geometry.polygon.Polygon;
import org.poly2tri.geometry.polygon.PolygonPoint;
import org.poly2tri.triangulation.TriangulationAlgorithm;
import org.poly2tri.triangulation.TriangulationContext;
import org.poly2tri.triangulation.TriangulationPoint;
import org.poly2tri.triangulation.delaunay.DelaunayTriangle;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Takes a 2D bit array and converts it into polygons with holes.
 *
 * Also implements access to the Poly2Tri algorithm to convert the created polygons into triangles.
 *
 * Note: The main difference between this algorithm and the JAI algorithm of "voxel -> polygon" is
 * that this one already combines holes that touch inside a polygon and furthermore already integrates
 * holes that touch the contour of the polygon. This is helpful as it is required for processing
 * with Poly2Tri, but could probably be changes by changing
 *    (edgeR[4] == 1 ? 1 : -1)
 * to
 *    (edgeR[4] == 1 ? -1 : 1)
 * two times (!). Verification still required.
 */
public class Grid2TriPolyFast {

    // interpolation value
    private static final double INTERP = 0.000001;

    private static final int MIN_ANGLE = 30;

    // ==============

    // helper - we need only one context for all conversion (faster)
    private final static TriangulationContext tcx = Poly2Tri.createContext(TriangulationAlgorithm.DTSweep);

    /**
     * Refresh a polygon. This is necessary after triangulation to be able to perform another triangulation.
     */
    private static void refresh_polygon(ArrayList<ArrayList<PolygonPoint>> polygon, List<TriangulationPoint> steinerPoints) {
        for (ArrayList<PolygonPoint> contour : polygon) {
            for (int i = 0; i < contour.size(); i++) {
                contour.set(i, new PolygonPoint(contour.get(i).getX(), contour.get(i).getY()));
            }
        }
        for (int i = 0; i < steinerPoints.size(); i++) {
            steinerPoints.set(i, new PolygonPoint(steinerPoints.get(i).getX(), steinerPoints.get(i).getY()));
        }
    }

    /**
     * Compute circle with radius r that has it center on the "left side" of the
     * line and has both points on it's border.
     */
    public static double[] find_line_circle_center(double x1, double y1, double x2, double y2, double r) {
        double[] middle_point = new double[] {(x1 + x2) / 2, (y1 + y2) / 2};
        double dist = Math.sqrt(Math.pow(x2 - x1, 2) + Math.pow(y2 - y1, 2));
        double[] direction = new double[] {(y1 - y2) / dist, (x2 - x1) / dist};
        return new double[] {
                middle_point[0] + Math.sqrt(Math.pow(r, 2) - Math.pow(dist / 2, 2)) * direction[0],
                middle_point[1] + Math.sqrt(Math.pow(r, 2) - Math.pow(dist / 2, 2)) * direction[1]
        };
    }

    /**
     * Do one step in the while loop of the referenced algorithm.
     * Reference: http://www.cs.berkeley.edu/~jrs/meshpapers/ErtenUngor.pdf
     */
    private static boolean degenerate_removal_step(DelaunayTriangle bad_tri, List<DelaunayTriangle> triangles, ArrayList<ArrayList<PolygonPoint>> polygon, List<TriangulationPoint> steinerPoints) {
        // find shortest side
        double[][] points = G2DUtil.get_triangle_points(bad_tri);
        double[] sides_length = G2DUtil.get_triangle_sides_length(points);
        double min_side_length;
        double[][] min_side;
        if (sides_length[0] < sides_length[1] && sides_length[0] < sides_length[2]) {
            min_side = new double[][] {points[1], points[2]};
            min_side_length = sides_length[0];
        } else if (sides_length[1] < sides_length[2]) {
            min_side = new double[][] {points[2], points[0]};
            min_side_length = sides_length[1];
        } else {
            min_side = new double[][] {points[0], points[1]};
            min_side_length = sides_length[2];
        }
        // compute circle center
        double r = 2 * Math.sin((MIN_ANGLE/180d) * Math.PI) * min_side_length;
        double[] center = find_line_circle_center(min_side[0][0], min_side[0][1], min_side[1][0], min_side[1][1], r);
        // collect all points
        ArrayList<TriangulationPoint> all_points = new ArrayList<TriangulationPoint>();
        for (ArrayList<PolygonPoint> contour : polygon) {
            all_points.addAll(contour);
        }
        all_points.addAll(steinerPoints);
        // this will hold all possible points
        ArrayList<double[]> possible_addition = new ArrayList<double[]>();
        for (DelaunayTriangle tri : triangles) {
            // add all triangle centers
            possible_addition.add(G2DUtil.computeCircumcenter(tri));
            // add all side centers
            possible_addition.add(G2DUtil.getCenter(tri.points[0].getX(), tri.points[0].getY(), tri.points[1].getX(), tri.points[1].getY()));
            possible_addition.add(G2DUtil.getCenter(tri.points[1].getX(), tri.points[1].getY(), tri.points[2].getX(), tri.points[2].getY()));
            possible_addition.add(G2DUtil.getCenter(tri.points[2].getX(), tri.points[2].getY(), tri.points[0].getX(), tri.points[0].getY()));
            // add tri line orthogonal intersections with circle
            double[][] line = G2DUtil.getOrthogonalLine(tri.points[0].getX(), tri.points[0].getY(), tri.points[1].getX(), tri.points[1].getY());
            Collections.addAll(possible_addition, G2DUtil.getIntersections(line[0][0], line[0][1], line[1][0], line[1][1], center[0], center[1], r));
            line = G2DUtil.getOrthogonalLine(tri.points[1].getX(), tri.points[1].getY(), tri.points[2].getX(), tri.points[2].getY());
            Collections.addAll(possible_addition, G2DUtil.getIntersections(line[0][0], line[0][1], line[1][0], line[1][1], center[0], center[1], r));
            line = G2DUtil.getOrthogonalLine(tri.points[2].getX(), tri.points[2].getY(), tri.points[0].getX(), tri.points[0].getY());
            Collections.addAll(possible_addition, G2DUtil.getIntersections(line[0][0], line[0][1], line[1][0], line[1][1], center[0], center[1], r));
            // add tri line intersections with circle
            Collections.addAll(possible_addition, G2DUtil.getIntersections(tri.points[0].getX(), tri.points[0].getY(), tri.points[1].getX(), tri.points[1].getY(), center[0], center[1], r));
            Collections.addAll(possible_addition, G2DUtil.getIntersections(tri.points[1].getX(), tri.points[1].getY(), tri.points[2].getX(), tri.points[2].getY(), center[0], center[1], r));
            Collections.addAll(possible_addition, G2DUtil.getIntersections(tri.points[2].getX(), tri.points[2].getY(), tri.points[0].getX(), tri.points[0].getY(), center[0], center[1], r));
        }
        // filter points
        for (int i = 0; i < possible_addition.size(); i++) {
            double[] p = possible_addition.get(i);
            // check if in circle
            if (!G2DUtil.contains(center[0], center[1], r, p[0], p[1])) {
                possible_addition.remove(i);
                i--;
                continue;
            }
            // check if a valid point
            boolean contained = false;
            for (DelaunayTriangle tri : triangles) {
                if (G2DUtil.onLine(tri.points[0].getX(), tri.points[0].getY(), tri.points[1].getX(), tri.points[1].getY(), p[0], p[1])) {
                    contained = true;
                    break;
                }
                if (G2DUtil.onLine(tri.points[1].getX(), tri.points[1].getY(), tri.points[2].getX(), tri.points[2].getY(), p[0], p[1])) {
                    contained = true;
                    break;
                }
                if (G2DUtil.onLine(tri.points[2].getX(), tri.points[2].getY(), tri.points[0].getX(), tri.points[0].getY(), p[0], p[1])) {
                    contained = true;
                    break;
                }
                if (G2DUtil.inTriangle(p[0], p[1], tri.points[0].getX(), tri.points[0].getY(), tri.points[1].getX(), tri.points[1].getY(), tri.points[2].getX(), tri.points[2].getY())) {
                    contained = true;
                    break;
                }
            }
            if (!contained) {
                possible_addition.remove(i);
                i--;
            }
        }
        // check which point is best
        double largest_min_dist = 0;
        double[] addition = null;
        for (double[] p : possible_addition) {
            // compute the distance
            double min_dist = Double.MAX_VALUE;
            for (TriangulationPoint vertex : all_points) {
                double dist = Math.sqrt(Math.pow(p[0] - vertex.getX(), 2) + Math.pow(p[1] - vertex.getY(), 2));
                if (dist < min_dist) {
                    min_dist = dist;
                }
            }
            if (largest_min_dist < min_dist) {
                largest_min_dist = min_dist;
                addition = p;
            }
        }
        if (largest_min_dist < 2) {
            addition = null;
        }
        if (addition != null) {
            // this might be a new point on a contour
            boolean added = false;
            edgeLoop: for (ArrayList<PolygonPoint> contour : polygon) {
                for (int i = 0; i < contour.size(); i++) {
                    int index2 = i + 1 < contour.size() ? i + 1 : 0;
                    PolygonPoint p1 = contour.get(i);
                    PolygonPoint p2 = contour.get(index2);
                    if (G2DUtil.onLine(p2.getX(), p2.getY(), p1.getX(), p1.getY(), addition[0], addition[1])) {
                        contour.add(index2, new PolygonPoint(addition[0], addition[1]));
                        added = true;
                        break edgeLoop;
                    }
                }
            }
            if (!added) {
                steinerPoints.add(new PolygonPoint(addition[0], addition[1]));
            }
            return true;
        }
        return false;
    }

    /**
     * Triangulate a list of polygons while removing degenerate triangles.
     */
    public static ArrayList<DelaunayTriangle> triangulate_degenerate_removal(short[][][] polys) {
        ArrayList<DelaunayTriangle> result = new ArrayList<DelaunayTriangle>();

        // loop over all polygon (a polygon consists of exterior and interior ring)
        for (short[][] poly : polys) {
            List<DelaunayTriangle> triangles; // stores an interpolated polygon ("0" entry is contour, others are holes)
            List<TriangulationPoint> steinerPoints = new ArrayList<TriangulationPoint>();
            boolean has_bad_triangles;
            ArrayList<ArrayList<PolygonPoint>> polygon = short_array_to_poly(poly);

            // ensure that there are no bad triangles
            do {
                has_bad_triangles = false;
                refresh_polygon(polygon, steinerPoints);
                triangles = triangulate_poly(polygon, steinerPoints);
                Collections.shuffle(triangles);

                for (DelaunayTriangle tri : triangles) {
                    if (G2DUtil.get_min_angle(tri) < MIN_ANGLE) {
                        if (degenerate_removal_step(tri, triangles, polygon, steinerPoints)) {
                            has_bad_triangles = true;
                            break;
                        }
                    }
                }

            } while (has_bad_triangles);

            result.addAll(triangles);

        }

        // return all triangles
        return result;
    }

    // =========================

    /**
     * Convert an contour (this could be a hole or a polygon) into a list of PolygonPoint. The PolygonPoints
     * are interpolated to avoid conflict with overlapping point. This is important for the poly2tri library.
     */
    private static ArrayList<PolygonPoint> parse_contour(short[] contour) {
        // stores and manages all seen points
        TIntHashSet indexer = new TIntHashSet();
        ArrayList<PolygonPoint> resultOutline = new ArrayList<PolygonPoint>();
        for (int i = 0, len = contour.length - 2; i < len; i+=2) {
            // check if we need to interpolate
            if (!indexer.add(IntegerTools.makeInt(contour[i], contour[i + 1]))) {
                resultOutline.add(new PolygonPoint(contour[i] - Math.signum(contour[i] - contour[i+2]) * INTERP, contour[i+1] - Math.signum(contour[i+1] - contour[i+3]) * INTERP));
            } else {
                resultOutline.add(new PolygonPoint(contour[i], contour[i+1]));
            }
        }
        return resultOutline;
    }

    /**
     * Convert a nested short polygon into a nest ArrayList of PolygonPoints
     */
    private static ArrayList<ArrayList<PolygonPoint>> short_array_to_poly(short[][] poly) {
        // stores an interpolated polygon ("0" entry is contour, others are holes)
        ArrayList<ArrayList<PolygonPoint>> polygon = new ArrayList<ArrayList<PolygonPoint>>();

        // loop over polygon contour (j=0) and holes (j>0)
        for (short[] contour : poly) {
            polygon.add(parse_contour(contour));
        }
        return polygon;
    }

    /**
     * Triangulate a polygon, optionally using Steiner points
     */
    private static List<DelaunayTriangle> triangulate_poly(ArrayList<ArrayList<PolygonPoint>> polygon, List<TriangulationPoint> steinerPoints) {
        // convert to polygon from raw data (zero is always the id that contains the exterior of the polygon)
        org.poly2tri.geometry.polygon.Polygon polyR = new Polygon(polygon.get(0));
        for (int j = 1, len = polygon.size(); j < len; j++) {
            polyR.addHole(new Polygon(polygon.get(j)));
        }
        if (steinerPoints != null) {
            polyR.addSteinerPoints(steinerPoints);
        }

        // do the triangulation and add the triangles for this polygon
        // Note: This needs to be synchronized to prevent multiple instances
        // from accessing the tcx context at once
        synchronized (tcx) {
            tcx.prepareTriangulation(polyR);
            Poly2Tri.triangulate(tcx);
            tcx.clear();
        }

        return polyR.getTriangles();
    }

    /**
     * Triangulate a polygon
     */
    private static List<DelaunayTriangle> triangulate_poly(ArrayList<ArrayList<PolygonPoint>> polygon) {
        return triangulate_poly(polygon, null);
    }

    // triangulate a polygon, the input data is interpolated to allow Poly2Tri to process it.
    // Hence the output data is slightly "off". This can be fixed by rounding the output data, don't use (int)
    // casting though as this might round down instead of up.
    public static ArrayList<DelaunayTriangle> triangulate(short[][][] polys) {
        ArrayList<DelaunayTriangle> result = new ArrayList<DelaunayTriangle>();

        // loop over all polygon (a polygon consists of exterior and interior ring)
        for (short[][] poly : polys) {
            // stores an interpolated polygon ("0" entry is contour, others are holes)
            ArrayList<ArrayList<PolygonPoint>> polygon = short_array_to_poly(poly);
            result.addAll(triangulate_poly(polygon));
        }

        // return all triangles
        return result;
    }

}
