package com.vitco.low.triangulate;

import com.vitco.util.graphic.Util3D;
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

    // ==============

    // helper - we need only one context for all conversion (faster)
    private final static TriangulationContext tcx = Poly2Tri.createContext(TriangulationAlgorithm.DTSweep);

    /**
     * Do one step in the while loop of the referenced algorithm.
     * Reference: http://www.cs.berkeley.edu/~jrs/meshpapers/ErtenUngor.pdf
     */
    private static void degenerate_removal_step(DelaunayTriangle tri, ArrayList<ArrayList<PolygonPoint>> polygon, List<TriangulationPoint> steinerPoints) {
        double[][] points = Util3D.get_triangle_points(tri);
        double[] sides_length = Util3D.get_triangle_sides_length(points);
        TriangulationPoint[] longest_side;
        if (sides_length[0] > sides_length[1] && sides_length[0] > sides_length[2]) {
            longest_side = new TriangulationPoint[] {tri.points[1], tri.points[2]};
        } else if (sides_length[1] > sides_length[2]) {
            longest_side = new TriangulationPoint[] {tri.points[0], tri.points[2]};
        } else {
            longest_side = new TriangulationPoint[] {tri.points[0], tri.points[1]};
        }
        PolygonPoint middle_point = new PolygonPoint(
                (longest_side[0].getX() + longest_side[1].getX()) / 2f,
                (longest_side[0].getY() + longest_side[1].getY()) / 2f,
                (longest_side[0].getZ() + longest_side[1].getZ()) / 2f
        );
        PolygonPoint center_point = new PolygonPoint(
                (points[0][0] + points[1][0] + points[2][0]) / 3f,
                (points[0][1] + points[1][1] + points[2][1]) / 3f,
                (points[0][2] + points[1][2] + points[2][2]) / 3f
        );

        boolean added = false;
        // check if side is contained (then we add a point in the middle)
        mainLoop: for (ArrayList<PolygonPoint> contour : polygon) {
            for (int i = 0, len = contour.size(); i < len; i++) {
                PolygonPoint p = contour.get(i);
                if (p == longest_side[0] || p == longest_side[1]) {
                    PolygonPoint p2 = i+1 == len ? contour.get(0) : contour.get(i+1);
                    if (p2 == longest_side[0] || p2 == longest_side[1]) {
                        contour.add(i + 1, middle_point);
                        added = true;
                        break mainLoop;
                    }
                }
            }
        }
        if (!added) {
            steinerPoints.add(center_point);
        }
    }

    /**
     * Refresh a polygon. This is necessary after triangulation to be able to perform another triangulation.
     */
    private static void refresh_polygon(ArrayList<ArrayList<PolygonPoint>> polygon) {
        for (ArrayList<PolygonPoint> contour : polygon) {
            for (int i = 0; i < contour.size(); i++) {
                contour.set(i, new PolygonPoint(contour.get(i).getX(), contour.get(i).getY()));
            }
        }
    }

    /**
     * Triangulate a list of polygons while removing degenerate triangles.
     */
    public static ArrayList<DelaunayTriangle> triangulate_degenerate_removal(short[][][] polys) {
        ArrayList<DelaunayTriangle> result = new ArrayList<DelaunayTriangle>();

        // loop over all polygon (a polygon consists of exterior and interior ring)
        for (short[][] poly : polys) {
            List<DelaunayTriangle> triangles;
            List<TriangulationPoint> steinerPoints = new ArrayList<TriangulationPoint>();
            boolean has_bad_triangles;
            ArrayList<ArrayList<PolygonPoint>> polygon = short_array_to_poly(poly);

            do {
                has_bad_triangles = false;
                // stores an interpolated polygon ("0" entry is contour, others are holes)
                refresh_polygon(polygon);
                triangles = triangulate_poly(polygon, steinerPoints);

                double angle = Float.MAX_VALUE;
                DelaunayTriangle worst_triangle = null;
                for (DelaunayTriangle tri : triangles) {
                    double new_angle = Util3D.get_min_angle(tri);
                    if (new_angle < angle) {
                        angle = new_angle;
                        worst_triangle = tri;
                    }
                }
                if (angle < 20) {
                    degenerate_removal_step(worst_triangle, polygon, steinerPoints);
                    has_bad_triangles = true;
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
