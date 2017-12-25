import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import java.util.*;

public class PackRegularTetrahedra {

    public static final Vector3D origin = new Vector3D(new double[]{0, 0, 0});

    /**
     * Function to calculate the unit normal of the least-squares plane of a
     * given set of points in 3D space. Following the procedure outlined in
     * detail in:
     * Y. Kallus, V. Elser, and S. Gravel, "Method for dense packing
     * discovery," Physical Review E, vol. 82, no. 5, Nov. 2010.
     *
     * Short summary:
     * The algorithm proceeds to compute the centroid (r*) of the given set of
     * vectors as follows:
     * r* = sum_{r belongs to s} r / |s|, where r is the vector and |s| is the
     * size of the set.
     *
     * Then computes the symmetric matrix M as follows:
     * M = sum_{r belongs to s} (r - r*)^T (r - r*) , where ^T denotes the
     * transpose operation.
     *
     * Computes the eigenvalues and eigenvectors of M. The unit normal is the
     * eigenvector corresponding to the smallest eigenvalue.
     *
     * @param s A list containing the subset of vectors denoting the positions
     *         of vertices of the two intersecting regular tetrahedra.
     * @return An array containing the centroid of the given subset as the
     * first element and the unit normal of the least-squares plane as the
     * second element.
     */
    public static Vector3D[] findLeastSqNormal(List<Vector3D> s) {

        int sz = s.size();
        if (sz < 3) throw new RuntimeException("The size should be at least " +
                "3");

        Vector3D[] result = new Vector3D[2];

        // Compute centroid.
        Vector3D centroid = origin;
        for (Vector3D r : s) centroid = centroid.add(r);
        centroid = centroid.scalarMultiply(1.0 / sz);
        result[0] = centroid;

        // Compute matrix M.
        BlockRealMatrix M = new BlockRealMatrix(3, 3);
        BlockRealMatrix t = new BlockRealMatrix(3, 3);
        for (Vector3D r : s) {
            t.setColumn(0, r.subtract(centroid).toArray());
            M = M.add(t.multiply(t.transpose()));
        }

        // Compute eigenvalues and eigenvectors.
        EigenDecomposition eig = new EigenDecomposition(M);
        if (eig.hasComplexEigenvalues()) throw new RuntimeException("Oops, " +
                "system has complex eigenvalues!");

        // Compute minimum eigenvalue and corresponding eigenvector.
        double min = Double.POSITIVE_INFINITY;
        Vector3D normal = null;
        for (int i=0; i<3; i++) {
            double eValue = eig.getRealEigenvalue(i);
            if (eValue < min) {
                min = eValue;
                normal = new Vector3D(eig.getEigenvector(i).toArray());
            }
        }
        if (normal == null) throw new RuntimeException("Eigenvector is null!");
        if (normal.isNaN() || normal.isInfinite()) throw new RuntimeException
                ("Normal is NaN or infinity!");
        result[1] = normal;
        return result;
    }

    /**
     * Function to generate all subsets of a given list of vectors. We return
     * a list of lists here instead of a set of sets because two different
     * vertices of the two tetrahedra may coincide. If we used a set, we
     * would miss this case. Doesn't return empty set as a subset.
     *
     * Note: This algorithm only works if list size is less than 32 as we use
     * integer's binary representation to compute the subsets. For two
     * tetrahedra with maximum of 8 vectors in this list, this algorithm
     * should suffice.

     * @param v List of vectors.
     * @return The list of list containing all subsets of the given list.
     */
    public static List<List<Vector3D>> generateSubsets(List<Vector3D> v) {
        // Total size of list provided.
        int n = v.size();
        if (n <= 0) throw new IllegalArgumentException("Input list should " +
                "have at least 1 element!");

        // Generate subsets.
        List<List<Vector3D>> subsets = new ArrayList<>();
        for (int i=0; i<(1<<n); i++) {
            List<Vector3D> set = new ArrayList<>();
            for (int j=0; j<n; j++)
                // Check if bit is present or not.
                if ((i & (1 << j)) > 0) set.add(v.get(j));
            if (set.size() > 0) subsets.add(set);
        }
        return subsets;
    }

    /**
     * Function to compute whether two tetrahedra intersect or not. Based on
     * the algorithm explained in detail here:
     * Y. Kallus, V. Elser, and S. Gravel, "Method for dense packing
     * discovery," Physical Review E, vol. 82, no. 5, Nov. 2010.
     *
     * Short summary:
     * The algorithm proceeds to first compute all the subsets of the set
     * formed by the union operation on the sets of vertices of the two
     * tetrahedra (X1 and X2, say). From this, all subsets of length 3 are
     * chosen and that contain at least one element each from X1 and X2. A
     * least-squares plane is found for each of these subsets (S, say). The
     * least-squares plane is defined as the plane to which the sum of
     * squared distances of all vertices in a set S is the lowest. Based on
     * this, a measure of interpenetration \Delta^2 is computed. If \Delta^2
     * is greater than zero, the two tetrahedra intersect. If \Delta^2 = 0,
     * they don't.
     *
     * @param rt1 Regular tetrahedron instance 1.
     * @param rt2 Regular tetrahedron instance 2.
     * @return Whether or not they intersect.
     */
    public static boolean isOverlap(RegularTetrahedron rt1,
                                    RegularTetrahedron rt2) {

        // List of vertices of X1.
        List<Vector3D> X1 = Arrays.asList(rt1.getAllVertices());
        // All possible subsets.
        List<List<Vector3D>> subsetsX1 = generateSubsets(X1);

        // List of vertices of X2.
        List<Vector3D> X2 = Arrays.asList(rt2.getAllVertices());
        // All possible subsets.
        List<List<Vector3D>> subsetsX2 = generateSubsets(X2);

        // List of vectors formed by the union of X1 and X2.
        List<Vector3D> X1UX2 = new ArrayList<>();
        X1UX2.addAll(X1);
        X1UX2.addAll(X2);
        // All possible subsets.
        List<List<Vector3D>> subsetsX1UX2 = generateSubsets(X1UX2);

        int d = 3;
        // Of all possible subsets of X1 union X2, select those with size of 3.
        List<List<Vector3D>> S = new ArrayList<>();
        for (List<Vector3D> v : subsetsX1UX2) if (v.size() == d) S.add(v);

        List<List<Vector3D>> x1 = new ArrayList<>();
        for (List<Vector3D> v : subsetsX1) if (v.size() == d) x1.add(v);
        List<List<Vector3D>> x2 = new ArrayList<>();
        for (List<Vector3D> v : subsetsX2) if (v.size() == d) x2.add(v);

        // From the list S, remove those that don't contain at least one
        // element each of X1 and X2.
        for (List<Vector3D> l1 : x1) S.remove(l1);
        for (List<Vector3D> l2 : x2) S.remove(l2);

        double minDS = Double.POSITIVE_INFINITY;
        for (List<Vector3D> s : S) {
            Vector3D[] result = findLeastSqNormal(s);
            Vector3D centroid = result[0];
            Vector3D normal = result[1];
            double h = normal.dotProduct(centroid);
            double dsp = 0;
            double dsm = 0;
            for (Vector3D x : X1) {
                double dp = normal.dotProduct(x);
                double value = (dp - h) * (dp - h);
                if (dp > h) dsp += value;
                else dsm += value;
            }
            for (Vector3D x : X2) {
                double dp = normal.dotProduct(x);
                double value = (dp - h) * (dp - h);
                if (dp > h) dsm += value;
                else dsp += value;
            }

            double ds = Math.min(dsp, dsm);
            minDS = Math.min(minDS, ds);
        }
        // If minDS is zero, then a separating plane exists, which implies
        // that the input tetrahedra don't overlap.
        return minDS > 0;
    }

    /**
     * Function to move the two intersecting tetrahedra by the smallest
     * displacement of their vertices that resolves the overlap. Following
     * the procedure outline in detail in:
     * Y. Kallus, V. Elser, and S. Gravel, "Method for dense packing
     * discovery," Physical Review E, vol. 82, no. 5, Nov. 2010.
     *
     * The set S is defined to be the union of the sets containing the
     * vertices of the two tetrahedra (X1, X2). The algorithm proceeds to find
     * the unit normal of the least-squares plane of all subsets of S that
     * have at least 4 elements in it and are made up of at least one element
     * of each set. This is used to compute the magnitude of the displacement
     * needed to resolve the collision. We iterate through all the possible
     * subsets and find the minimum such displacement and apply that to the
     * particular subset of S.
     *
     * @param rt1 Regular tetrahedron instance 1.
     * @param rt2 Regular tetrahedron instance 2.
     */
    public static void resolveOverlap(RegularTetrahedron rt1,
                                      RegularTetrahedron rt2) {

        // To compute all the subsets of S that have at least d elements (d =
        // 4, here) in it and have at least one element from X1 and X2, we
        // simply get all the subsets of S that have at least d elements in
        // it and subtract from it the subsets of X1 and X2 that have exactly
        // d elements in it. Since X1 and X2 have exactly 4 elements, we just
        // have to subtract them from the list of all subsets.

        assert (!isOverlap(rt1, rt2));
        // List of vertices of X1.
        List<Vector3D> X1 = Arrays.asList(rt1.getAllVertices());
        // All possible subsets.
        List<List<Vector3D>> subsetsX1 = generateSubsets(X1);

        // List of vertices of X2.
        List<Vector3D> X2 = Arrays.asList(rt2.getAllVertices());
        // All possible subsets.
        List<List<Vector3D>> subsetsX2 = generateSubsets(X2);

        // List of vectors formed by the union of X1 and X2.
        List<Vector3D> X1UX2 = new ArrayList<>();
        X1UX2.addAll(X1);
        X1UX2.addAll(X2);
        // All possible subsets.
        List<List<Vector3D>> subsetsX1UX2 = generateSubsets(X1UX2);

        int d = 3;
        // Of all possible subsets of X1 union X2, select those with size
        // greater than 3.
        List<List<Vector3D>> S = new ArrayList<>();
        for (List<Vector3D> v : subsetsX1UX2) if (v.size() > d) S.add(v);

        // From the list S, remove those that don't contain at least one
        // element each of X1 and X2. Since for tetrahedra, X1 and X2 have
        // only 4 vertices each, we only need to remove two lists from S -
        // one each containing all elements of X1 and X2.
        S.remove(X1);
        S.remove(X2);

        double minDS = Double.POSITIVE_INFINITY;
        List<Vector3D> T = null;
        double hT = Double.POSITIVE_INFINITY;
        Vector3D nT = null;
        for (List<Vector3D> s : S) {
            Vector3D[] result = findLeastSqNormal(s);
            Vector3D centroid = result[0];
            Vector3D normal = result[1];
            double h = normal.dotProduct(centroid);

            // Get the list of vectors in X1 and X2 that are not in S.
            List<Vector3D> X1S = new ArrayList<>(X1);
            List<Vector3D> X2S = new ArrayList<>(X2);
            for (Vector3D v : s) {
                X1S.remove(v);
                X2S.remove(v);
            }

            // Check if the least-squares plane separates the two sets of
            // vertices X1S and X2S.
            double ds2 = Double.POSITIVE_INFINITY;
            if (doesPlaneSeparate(normal, X1S, X2S)) ds2 =
                    sumSquaredDist(normal, s);
            if (ds2 < minDS) {
                minDS = ds2;
                hT = h;
                nT = normal;
                T = s;
            }
        }

        if (nT == null) throw new RuntimeException("Minimization went wrong!");

        for (Vector3D x : T) {
            Vector3D displacement = new Vector3D((hT - x.dotProduct(nT)), nT);
            int x1Index = X1.indexOf(x);
            int x2Index = X2.indexOf(x);
            if (x1Index != -1) rt1.translateVertex(x1Index, displacement);
            else if (x2Index != -1) rt2.translateVertex(x2Index, displacement);
            else throw new RuntimeException("Vertex not found in either X1 or" +
                        " X2!!!");
        }
    }

    /**
     * Function to check if a plane separates two given sets of vectors.
     * Compares the sign of the dot product of the given sets of vectors to
     * return true only when the two sets of vectors have opposite signs.
     * @param normal Unit normal vector of the given plane.
     * @param set1 First set of vectors.
     * @param set2 Second set of vectors.
     * @return Whether the plane separates the two sets of vectors or not.
     */
    public static boolean doesPlaneSeparate(Vector3D normal, List<Vector3D>
            set1, List<Vector3D> set2) {
        boolean positive = false;
        boolean negative = false;
        boolean flag = false;
        boolean separates = true;
        for (Vector3D v : set1) {
            double dot = v.dotProduct(normal);
            if (!flag) {
                // Flag has not been set yet. Meaning we have not decided if
                // this set has positive or negative sign. This will occur
                // until we encounter the first non-zero valued dot product.
                flag = true;
                if (dot > 0) positive = true;
                else if (dot < 0) negative = true;
                else flag = false;
            }
            else {
                // If one of the vectors changes sign, no need to check the
                // rest. We can conclude that the plane doesn't separate.
                if ((positive && dot < 0) || (negative && dot > 0)) {
                    separates = false;
                    break;
                }
            }
        }

        if (separates) {
            for (Vector3D v : set2) {
                double dot = v.dotProduct(normal);
                if (!flag) {
                    flag = true;
                    if (dot > 0) positive = true;
                    else if (dot < 0) negative = true;
                    else flag = false;
                } else {
                    if ((positive && dot > 0) || (negative && dot < 0)) {
                        separates = false;
                        break;
                    }
                }
            }
        }
        return separates;
    }

    /**
     * Function to compute the sum of squared distances of a list of vectors
     * to a plane.
     * @param normal Unit normal vector of the given plane.
     * @param s List of desired vectors.
     * @return Sum of squared distances.
     */
    public static double sumSquaredDist(Vector3D normal, List<Vector3D> s) {
        double squaredDist = 0.0;
        for (Vector3D v : s) {
            double dot = normal.dotProduct(v);
            squaredDist += dot * dot;
        }
        return squaredDist;
    }

    public static void main(String[] args) {
//        System.out.println("Integer.parseInt(\"1110\") = " + Integer.parseInt
//                ("11100000", 2));
//        System.out.println("Integer.parseInt(\"01101\") = " + Integer.parseInt
//                ("01101000", 2));

//        System.out.println("Math.min(1e-31, 0) = " + Math.min(1e-31, 0));
//        System.exit(1);
//        List<Vector3D> l1 = new ArrayList<>();
//        Vector3D v1 = new Vector3D(new double[]{0, 0, 2});
//        Vector3D v2 = new Vector3D(new double[]{2, 0, 2});
//        Vector3D v3 = new Vector3D(new double[]{2, 3, 2});
//        l1.add(v1);
//        l1.add(v2);
//        l1.add(v3);
//        List<Vector3D> l2 = new ArrayList<>();
//        Vector3D w1 = new Vector3D(new double[]{0, 0, 4});
//        Vector3D w2 = new Vector3D(new double[]{-2, 0, 4});
//        Vector3D w3 = new Vector3D(new double[]{-4, -2, 4});
//        l2.add(w1);
//        l2.add(w2);
//        l2.add(w3);
//        Vector3D n = origin;
//        List<Vector3D> l3 = new ArrayList<>();
//        Vector3D t1 = new Vector3D(new double[]{0, 0, 0});
//        Vector3D t2 = new Vector3D(new double[]{-2, 0, 2});
//        Vector3D t3 = new Vector3D(new double[]{-2, -5, -9});
//        l3.add(t1);
//        l3.add(t2);
//        l3.add(t3);
//        n = n.add(Vector3D.PLUS_I);
//        n = n.add(Vector3D.PLUS_J);
//        n = n.add(Vector3D.PLUS_K);
//        n.normalize();
//        double d = doesPlaneSeparate(n, l1, l2, l3);
//        System.out.println(d);
//        l.sort(new Comparator<Vector3D>() {
//            @Override
//            public int compare(Vector3D v1, Vector3D v2) {
//                return Double.compare(v1.getNorm(), v2.getNorm());
//            }
//        });
//        for (Vector3D v : l) System.out.println(v.toString());
//        v.remove(v1);
//        for (List<Vector3D> s : generateSubsets(1, setSize
//                .AT_LEAST, v))
//            System.out.println(s.toString());
//        System.exit(1);

        RegularTetrahedron r1 = new RegularTetrahedron();
        r1.vertices[0] = new Vector3D(new double[]{2328.824,1256.559,
                9965.406});//0.0, 0.0, 0.0});
//        List<Vector3D> l1 = new ArrayList<>();
//        l1.add(r1.vertices[0]);

//        System.out.println("r1.vertices[0].equals(new Vector3D(new " +
//                "double[]{0, 0, 0})) = " + r1.vertices[0].equals(new Vector3D
//                (new double[]{0, 0, 0})));
//        System.exit(1);
        r1.vertices[1] = new Vector3D(new double[]{2328.824, 1256.559,
                9889.206});//50.0, 50.0, 0.0});
        r1.vertices[2] = new Vector3D(new double[]{2208.937, 1187.342,
                9965.406});//100.0, 0.0, 0.0});
        r1.vertices[3] = new Vector3D(new double[]{2132.737,1319.324,
                9965.406});//50.0, 25.0, 100.0});
//
        RegularTetrahedron r2 = new RegularTetrahedron();
        r2.vertices[0] = new Vector3D(new double[]{2137.673 ,1234.186 ,
                10003.130});//0.0, 0.0, 0.0});
//        List<Vector3D> l2 = new ArrayList<>();
//        l2.add(r2.vertices[0]);
//        System.out.println("l2.equals(l1) = " + l2.equals(l1));
//        System.exit(1);
        r2.vertices[1] = new Vector3D(new double[]{2060.557 ,1281.337 ,
                9973.415});//-50.0, 50.0, 0.0});
        r2.vertices[2] = new Vector3D(new double[]{2096.147, 1303.468,
                9925.334});//-200.0, 0.0, 20.0});
        r2.vertices[3] = new Vector3D(new double[]{2101.378,1213.231,
                9973.415});//150.0, 25.0, 100.0});
        System.out.println("isOverlap(r1, r2) = " + isOverlap(r1, r2));
        resolveOverlap(r1, r2);
        System.out.println("isOverlap(r1, r2) = " + isOverlap(r1, r2));
        for (Vector3D v : r1.vertices) System.out.print(v.toString()+" ");
        System.out.println();
        for (Vector3D v : r2.vertices) System.out.print(v.toString()+" ");
        System.out.println();
//        Plane p0_012 = new Plane(r1.vertices[0], r1.vertices[1], r1.vertices[2],
//                1e-10);
//        Plane p0_123 = new Plane(r1.vertices[1], r1.vertices[2], r1.vertices[3],
//                1e-10);
//        Plane p0_230 = new Plane(r1.vertices[2], r1.vertices[3], r1.vertices[0],
//                1e-10);
//        Plane p0_301 = new Plane(r1.vertices[3], r1.vertices[0], r1.vertices[1],
//                1e-10);
//        Plane p1_012 = new Plane(r2.vertices[0], r2.vertices[1], r2.vertices[2],
//                1e-10);
//        Plane p1_123 = new Plane(r2.vertices[1], r2.vertices[2], r2.vertices[3],
//                1e-10);
//        Plane p1_230 = new Plane(r2.vertices[2], r2.vertices[3], r2.vertices[0],
//                1e-10);
//        Plane p1_301 = new Plane(r2.vertices[3], r2.vertices[0], r2.vertices[1],
//                1e-10);
//        System.out.println("p0_012 = " + p0_012.intersection(p1_012));
//        System.out.println("p0_012 = " + p0_012.intersection(p1_123));
//        System.out.println("p0_012 = " + p0_012.intersection(p1_230));
//        System.out.println("p0_012 = " + p0_012.intersection(p1_301));
//
//        System.out.println("p0_123 = " + p0_123.intersection(p1_012));
//        System.out.println("p0_123 = " + p0_123.intersection(p1_123));
//        System.out.println("p0_123 = " + p0_123.intersection(p1_230));
//        System.out.println("p0_123 = " + p0_123.intersection(p1_301));
//
//        System.out.println("p0_230 = " + p0_230.intersection(p1_012));
//        System.out.println("p0_230 = " + p0_230.intersection(p1_123));
//        System.out.println("p0_230 = " + p0_230.intersection(p1_230));
//        System.out.println("p0_230 = " + p0_230.intersection(p1_301));
//
//        System.out.println("p0_301 = " + p0_301.intersection(p1_012));
//        System.out.println("p0_301 = " + p0_301.intersection(p1_123));
//        System.out.println("p0_301 = " + p0_301.intersection(p1_230));
//        System.out.println("p0_301 = " + p0_301.intersection(p1_301));
//        System.out.println("isOverlapLazy(r1, r2) = " + isOverlapLazy(r1,
//                r2));
//        System.out.println("isIntersecting2(r1, r2) = " + isIntersecting2(r1,
//                r2));

//        resolveCollision(r1, r2);
//        r2.vertices[0] = new Vector3D(new double[]{100.0001, 0.0, 0.0});
//        r2.vertices[1] = new Vector3D(new double[]{150.0, 50.0, 0.0});
//        r2.vertices[2] = new Vector3D(new double[]{200.0, 0.0, 0.0});
//        r2.vertices[3] = new Vector3D(new double[]{150.0, 25.0, 10.0});

//        System.out.println("isOverlapLazy(r1, r2) = " + isOverlapLazy(r1,
//                r2));
//        r1.vertices[0] = new Vector3D(new double[]{2328.824, 1256.559,
//                9965.406});
//        r1.vertices[1] = new Vector3D(new double[]{2328.824, 1256.559,
//                9889.206});
//        r1.vertices[2] = new Vector3D(new double[]{2208.937, 1187.342,
//                9965.406});
//        r1.vertices[3] = new Vector3D(new double[]{2132.737, 1319.324,
//                9965.406});
//        r2.vertices[0] = new Vector3D(new double[]{2137.673, 1234.186,
//                10003.130});
//        r2.vertices[1] = new Vector3D(new double[]{2060.557, 1281.337,
//                9973.415});
//        r2.vertices[2] = new Vector3D(new double[]{2096.147, 1303.468,
//                9925.334});
//        r2.vertices[3] = new Vector3D(new double[]{2101.378, 1213.231,
//                9973.415});
//        System.out.println("isOverlapLazy(rt1, rt2) = " + isOverlapLazy(r1,
//                r2));
//        resolveOverlap(r1, r2);
//        System.out.println("resolveOverlap(r1, r2) = " + resolveOverlap(r1,
//                r2));
//        System.out.println("isIntersecting2(r1, r2) = " + isIntersecting2(r1,
//                r2));
//        BlockRealMatrix m = new BlockRealMatrix(3, 3);
//        m.setRow(1, new double[]{1, 2, 3});
//        double[][] x = m.getData();
//        for (double[] x1 : x) {
//            for (double x2 : x1) System.out.print(x2+" ");
//            System.out.println();
//        }
//        (def p [[2328.824 1256.559 9965.406] [2328.824 1256.559 9889.206] [2208.937 1187.342 9965.406] [2132.737 1319.324 9965.406]])
//        (def q [[2137.673 1234.186 10003.130] [2060.557 1281.337 9973.415] [2096.147 1303.468 9925.334] [2101.378 1213.231 9973.415]])
//        (def p [[0.0 0.0 0.0] [50.0 50.0 0.0] [100.0 0.0 0.0] [50.0 25.0 100.0]])
//        ;(def q [[0.0 0.0 0.0] [-50.0 50.0 0.0] [-200.0 0.0 20.0] [150.0 25.0 100.0]])
    }
}

