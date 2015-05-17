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
import java.util.HashMap;
import java.util.List;

/**
 * Takes a 2D bit array and converts it into polygons with holes.
 *
 * Also implements access to the Poly2Tri algorithm to convert the created polygons into triangles.
 *
 * Note: The main difference between this algorithm and the JAI algorithm of "voxel -> polygon" is
 * that this one already combines holes that touch inside a polygon and furthermore already integrates
 * holes that touch the outline of the polygon. This is helpful as it is required for processing
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
     * Convert an outline (this could be a hole or a polygon) into a list of PolygonPoint. The PolygonPoints
     * are interpolated to avoid conflict with overlapping point. This is important for the poly2tri library.
     */
    private static ArrayList<PolygonPoint> parse_contour(short[] outline) {
        // stores and manages all seen points
        TIntHashSet indexer = new TIntHashSet();
        ArrayList<PolygonPoint> resultOutline = new ArrayList<PolygonPoint>();
        for (int i = 0, len = outline.length - 2; i < len; i+=2) {
            // check if we need to interpolate
            if (!indexer.add(IntegerTools.makeInt(outline[i], outline[i + 1]))) {
                resultOutline.add(new PolygonPoint(outline[i] - Math.signum(outline[i] - outline[i+2]) * INTERP, outline[i+1] - Math.signum(outline[i+1] - outline[i+3]) * INTERP));
            } else {
                resultOutline.add(new PolygonPoint(outline[i], outline[i+1]));
            }
        }
        return resultOutline;
    }

    /**
     * Do one step in the while loop of the referenced algorithm.
     * Reference: http://www.cs.berkeley.edu/~jrs/meshpapers/ErtenUngor.pdf
     */
    private static void prevent_degenerate(DelaunayTriangle tri, HashMap<Integer, ArrayList<PolygonPoint>> polygon, List<TriangulationPoint> steinerPoints) {
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

        System.out.println("Triangle:" + tri.points[0] + " " + tri.points[1] + " " + tri.points[2]);

        boolean added = false;
        // check if side is contained (then we add a point in the middle)
        mainLoop: for (ArrayList<PolygonPoint> contour : polygon.values()) {
            for (int i = 0, len = contour.size(); i < len; i++) {
                PolygonPoint p = contour.get(i);
                if (p == longest_side[0] || p == longest_side[1]) {
                    PolygonPoint p2 = i+1 == len ? contour.get(0) : contour.get(i+1);
                    if (p2 == longest_side[0] || p2 == longest_side[1]) {
                        System.out.println("Edge " + middle_point);
                        contour.add(i + 1, middle_point);
                        added = true;
                        break mainLoop;
                    }
                }
            }
        }
        if (!added) {
            System.out.println("Steiner " + center_point);
            steinerPoints.add(center_point);
        }

        boolean first = true;
        for (ArrayList<PolygonPoint> contour : polygon.values()) {
            System.out.println(first ? "First" : "Hole");
            for (PolygonPoint p : contour) {
                System.out.println(p);
            }
            first = false;
        }
    }

    private static void remove_degenerate_triangles() {

    }

    // triangulate a polygon, the input data is interpolated to allow Poly2Tri to process it.
    // Hence the output data is slightly "off". This can be fixed by rounding the output data, don't use (int)
    // casting though as this might round down instead of up.
    public static ArrayList<DelaunayTriangle> triangulate(short[][][] polys) {
        ArrayList<DelaunayTriangle> result = new ArrayList<DelaunayTriangle>();

        // loop over all polygon (a polygon consists of exterior and interior ring)
        for (short[][] poly : polys) {

            // stores an interpolated polygon ("0" entry is outline, others are holes)
            HashMap<Integer, ArrayList<PolygonPoint>> polygon = new HashMap<Integer, ArrayList<PolygonPoint>>();

            // loop over polygon outline (j=0) and holes (j>0)
            for (int j = 0; j < poly.length; j++) {
                polygon.put(j, parse_contour(poly[j]));
            }

            List<DelaunayTriangle> triangles;
            boolean has_bad_triangles;
            List<TriangulationPoint> steinerPoints = new ArrayList<TriangulationPoint>();

            boolean first = true;
            for (ArrayList<PolygonPoint> contour : polygon.values()) {
                System.out.println(first ? "First" : "Hole");
                for (PolygonPoint p : contour) {
                    System.out.println(p);
                }
                first = false;
            }

            do {
                HashMap<Integer, ArrayList<PolygonPoint>> polygon2 = new HashMap<Integer, ArrayList<PolygonPoint>>();
                int i = 0;
                for (ArrayList<PolygonPoint> contour : polygon.values()) {
                    ArrayList<PolygonPoint> list = new ArrayList<PolygonPoint>();
                    for (PolygonPoint p : contour) {
                        list.add(new PolygonPoint(p.getX(), p.getY()));
                    }
                    polygon2.put(i, list);
                    i++;
                }
                polygon = polygon2;
                List<TriangulationPoint> steinerPoints2 = new ArrayList<TriangulationPoint>();
                for (TriangulationPoint p : steinerPoints) {
                    steinerPoints2.add(new PolygonPoint(p.getX(), p.getY()));
                }
                steinerPoints = steinerPoints2;

                has_bad_triangles = false;
                // convert to polygon from raw data (zero is always the id that contains the exterior of the polygon)
                org.poly2tri.geometry.polygon.Polygon polyR = null;
                for (ArrayList<PolygonPoint> contour : polygon.values()) {
                    Polygon outline = new Polygon(contour);
                    if (polyR == null) {
                        polyR = outline;
                    } else {
                        polyR.addHole(outline);
                    }
                }
                polyR.addSteinerPoints(steinerPoints);

                // do the triangulation and add the triangles for this polygon
                // Note: This needs to be synchronized to prevent multiple instances
                // from accessing the tcx context at once
                synchronized (tcx) {
                    tcx.prepareTriangulation(polyR);
                    Poly2Tri.triangulate(tcx);
                    tcx.clear();
                }
                triangles = polyR.getTriangles();

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
                    prevent_degenerate(worst_triangle, polygon, steinerPoints);
                    has_bad_triangles = true;
                }

            } while (has_bad_triangles);

            System.out.println("==========");

            result.addAll(triangles);

        }

        // return all triangles
        return result;
    }

}
