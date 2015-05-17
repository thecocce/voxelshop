package com.vitco.util.graphic;

import com.threed.jpct.SimpleVector;
import org.poly2tri.triangulation.delaunay.DelaunayTriangle;

/**
 * Helper class to perform 3D calculation tasks
 */
public final class Util3D {

    /**
     * Get the points of a DelaunayTriangle as a double list.
     */
    public static double[][] get_triangle_points(DelaunayTriangle tri) {
        return new double[][] {
                new double[] {tri.points[0].getX(), tri.points[0].getY(), tri.points[0].getZ()},
                new double[] {tri.points[1].getX(), tri.points[1].getY(), tri.points[1].getZ()},
                new double[] {tri.points[2].getX(), tri.points[2].getY(), tri.points[2].getZ()}
        };
    }

    /**
     * Get the side length of a DelaunayTriangle as a double array. The first side is the side opposite of the
     * first point and so on.
     */
    public static double[] get_triangle_sides_length(double[][] points) {
        double abs_a = Math.sqrt(Math.pow(points[1][0] - points[2][0], 2) + Math.pow(points[1][1] - points[2][1], 2) + Math.pow(points[1][2] - points[2][2], 2));
        double abs_b = Math.sqrt(Math.pow(points[0][0] - points[2][0], 2) + Math.pow(points[0][1] - points[2][1], 2) + Math.pow(points[0][2] - points[2][2], 2));
        double abs_c = Math.sqrt(Math.pow(points[0][0] - points[1][0], 2) + Math.pow(points[0][1] - points[1][1], 2) + Math.pow(points[0][2] - points[1][2], 2));
        return new double[] {abs_a, abs_b, abs_c};
    }

    /**
     * Obtain the minimal angle of a DelaunayTriangle.
     */
    public static double get_min_angle(DelaunayTriangle tri) {
        double[][] points = get_triangle_points(tri);
        double[] side_length = get_triangle_sides_length(points);
        double max_side_length = Math.max(Math.max(side_length[0], side_length[1]), side_length[2]);
        double med_side_length = Math.max(
                Math.min(side_length[0], side_length[1]),
                Math.min(Math.max(side_length[0], side_length[1]), side_length[2])
        );
        return Math.acos(med_side_length / max_side_length) * 180/Math.PI;
    }


    public static SimpleVector nearestPoint(SimpleVector[] line, SimpleVector point) {
        SimpleVector a = new SimpleVector(line[0]);
        a.sub(point);

        SimpleVector b = new SimpleVector(line[1]);
        b.sub(line[0]);

        float top = a.calcDot(b);

        double t = - top / Math.pow(line[1].distance(line[0]),2.0);

        return new SimpleVector(
                line[0].x + (line[1].x - line[0].x) * t,
                line[0].y + (line[1].y - line[0].y) * t,
                line[0].z + (line[1].z - line[0].z) * t
        );
    }

    public static SimpleVector[] nearestLinePoints(SimpleVector[] line1, SimpleVector[] line2) {
        SimpleVector d21 = new SimpleVector(line1[0]);
        d21.sub(line1[1]); // 3 sub x 3
        SimpleVector d34 = new SimpleVector(line2[1]);
        d34.sub(line2[0]);
        SimpleVector d13 = new SimpleVector(line2[0]);
        d13.sub(line1[0]);

        // m * u = x
        float a = d21.calcDot(d21); // (3 mul + 3 add ) x 5
        float b = d21.calcDot(d34);
        float c = d34.calcDot(d34);
        float d = -d13.calcDot(d21);
        float e = -d13.calcDot(d34);

        // Solve for u1 & u2
        float[] u = new float[2];
        u[0] = (d*c-e*b)/(c*a-b*b); // 4 mul, 2 sub, 1 div
        u[1] = (e - b * u[0]) / c; // 1 mul, 1 sub, 1 div

        return new SimpleVector[] {
                new SimpleVector(
                        line1[0].x + (line1[1].x - line1[0].x) * u[0],
                        line1[0].y + (line1[1].y - line1[0].y) * u[0],
                        line1[0].z + (line1[1].z - line1[0].z) * u[0]
                ),
                new SimpleVector(
                        line2[0].x + (line2[1].x - line2[0].x) * u[1],
                        line2[0].y + (line2[1].y - line2[0].y) * u[1],
                        line2[0].z + (line2[1].z - line2[0].z) * u[1]
                )
        };

    }

    // test if a ray intersects a triangle (single sided parameter indicates
    // whether only testing from the front should be performed)
    // Note: This assumes that the direction is normalized (!)
    // returns null if no hit is detected, otherwise the hit position is returned
    public static SimpleVector rayTriangleIntersects(
            SimpleVector v0, SimpleVector v1, SimpleVector v2,
            SimpleVector origin, SimpleVector dir,
            boolean isSingledSided
    ) {
        assert Math.abs(1 - dir.length()) < 0.0001;
        SimpleVector v0v1 = v1.calcSub(v0);
        SimpleVector v0v2 = v2.calcSub(v0);

        SimpleVector N = v0v1.calcCross(v0v2);

        float nDotRay = N.calcDot(dir);

        // ray parallel to triangle or hit from the back
        if (nDotRay == 0 || (isSingledSided && nDotRay > 0)) return null;

        float d = N.calcDot(v0);
        float t = -(N.calcDot(origin) + d) / nDotRay;

        if (t < 0) return null; // ray behind triangle
        // inside-out test
        SimpleVector dist = new SimpleVector(dir);
        dist.scalarMul(t);
        SimpleVector pos = origin.calcAdd(dist);

        // inside-out test edge0
        SimpleVector v0p = pos.calcSub(v0);
        float v = N.calcDot(v0v1.calcCross(v0p));
        if (v < 0) return null; // P outside triangle

        // inside-out test edge1
        SimpleVector v1p = pos.calcSub(v1);
        SimpleVector v1v2 = v2.calcSub(v1);
        float w = N.calcDot(v1v2.calcCross(v1p));
        if (w < 0) return null; // P outside triangle

        // inside-out test edge2
        SimpleVector v2p = pos.calcSub(v2);
        SimpleVector v2v0 = v0.calcSub(v2);
        float u = N.calcDot(v2v0.calcCross(v2p));
        if (u < 0) return null; // P outside triangle

        // Note: t contains the distance
        return pos;
    }
}
