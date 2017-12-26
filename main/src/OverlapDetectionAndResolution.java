package src;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import java.util.*;

/**
 * Class to detect and resolve overlaps between two tetrahedra. The methods
 * implemented in this class are adapted from:
 * Y. Kallus, V. Elser, and S. Gravel, "Method for dense packing discovery,"
 * Physical Review E, vol. 82, no. 5, Nov. 2010.
 *
 * More specifically, this class is an implementation of the divide
 * projection as part of the divide and concur algorithm.
 *
 * @author Ram
 */
public class OverlapDetectionAndResolution {

    // List of subsets of all vertices of both the tetrahedra.
    private List<List<Vector3D>> subsetsX1UX2 = null;

    // List of subsets of all vertices of tetrahedron 1.
    private List<List<Vector3D>> subsetsX1 = null;

    // List of subsets of all vertices of tetrahedron 2.
    private List<List<Vector3D>> subsetsX2 = null;

    // List of vectors in the set S which correspond to the best
    // least-squares plane compared to all subsets considered.
    private List<Vector3D> T = null, X1 = null, X2 = null;

    // Unit normal of the best least-squares plane.
    private Vector3D nT = null;

    // Projection of the centroid of the set of vectors T on the nT.
    private double hT = Double.POSITIVE_INFINITY;

    // Regular tetrahedra instances.
    private RegularTetrahedron r1, r2;

    // Threshold used to decide if the two tetrahedra overlap or not. If
    // \Delta^2, the measure of interpenetration is above this threshold, we
    // conclude that the tetrahedra are overlapping. Default value: 1e-30.
    private double threshold = 1e-30;

    /**
     * Constructs the OverlapDetectionAndResolution object which can be
     * queried using various methods and fields as described in the contents
     * below.
     * @param rt1 Regular Tetrahedron instance 1.
     * @param rt2 Regular Tetrahedron instance 2.
     */
    public OverlapDetectionAndResolution(RegularTetrahedron rt1,
                                         RegularTetrahedron rt2) {
        r1 = rt1;
        r2 = rt2;
        updateSubsets();
    }

    /**
     * Empty constructor for convenience.
     */
    public OverlapDetectionAndResolution(){
    }

    /**
     * Function to update the list of vectors in the subsets of various sets.
     * Need to update subsets every time we update rt1 and/or rt2.
     */
    public void updateSubsets() {
        // List of vertices of X1.
        X1 = Arrays.asList(r1.getAllVertices());
        // All possible subsets.
        subsetsX1 = generateSubsets(X1);

        // List of vertices of X2.
        X2 = Arrays.asList(r2.getAllVertices());
        // All possible subsets.
        subsetsX2 = generateSubsets(X2);

        // List of vectors formed by the union of X1 and X2.
        List<Vector3D> X1UX2 = new ArrayList<>();
        X1UX2.addAll(X1);
        X1UX2.addAll(X2);
        // All possible subsets.
        subsetsX1UX2 = generateSubsets(X1UX2);
    }

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
    private Vector3D[] findLeastSqNormal(List<Vector3D> s) {

        int sz = s.size();
        if (sz < 3) throw new RuntimeException("The size should be at least " +
                "3");

        Vector3D[] result = new Vector3D[2];

        // Compute centroid.
        Vector3D centroid = Vector3D.ZERO;
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
    private List<List<Vector3D>> generateSubsets(List<Vector3D> v) {
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
     * @return Whether or not they intersect.
     */
    public boolean isOverlap() {
        if (r1 == null || r2 == null) throw new RuntimeException("Regular " +
                "tetrahedra r1 and r2 haven't been initialized. Please " +
                "initialize them by calling the appropriate constructor or " +
                "setter methods.");
        // Check if the two tetrahedra are exactly the same.
        if (r1.equals(r2)) return true;

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
        if ((minDS == 0) || (minDS <= threshold)) return false;
        return true;
    }

    /**
     * Function to compute the best least-squares plane that results in the
     * smallest displacement of the vertices that resolves the overlap.
     * Following the procedure outlined in detail in:
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
     * @return Vector of displacements, one for each tetrahedron.
     */
    public Vector3D[] findDisplacementsToResolve(boolean apply) {
        if (r1 == null || r2 == null) throw new RuntimeException("Regular " +
                "tetrahedra r1 and r2 haven't been initialized. Please " +
                "initialize them by calling the appropriate constructor or " +
                "setter methods.");
        // To compute all the subsets of S that have at least d elements (d =
        // 4, here) in it and have at least one element from X1 and X2, we
        // simply get all the subsets of S that have at least d elements in
        // it and subtract from it the subsets of X1 and X2 that have exactly
        // d elements in it. Since X1 and X2 have exactly 4 elements, we just
        // have to subtract them from the list of all subsets.

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
        if (r1.isRegular() && r2.isRegular()) return computeDisplacements
                (apply);
        else return compDisplacements(apply);
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
    private boolean doesPlaneSeparate(Vector3D normal, List<Vector3D>
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
    private double sumSquaredDist(Vector3D normal, List<Vector3D> s) {
        double squaredDist = 0.0;
        for (Vector3D v : s) {
            double dot = normal.dotProduct(v);
            squaredDist += dot * dot;
        }
        return squaredDist;
    }

    /**
     * Function to compute the displacements of the vertices that belong to
     * T, the set which has the best least-squares plane.
     * @param apply Flag to apply the displacements directly to the
     *              tetrahedra. Mainly useful for debugging purposes.
     * @return Computed displacements for the two tetrahedra as an array
     * of vectors.
     */
    public Vector3D[] compDisplacements(boolean apply) {
        if (T == null || nT == null || !Double.isFinite(hT)) throw new
                RuntimeException("Overlap hasn't been resolved yet! Resolve " +
                "overlap and then make displacements");

        Vector3D[][] displacements = new Vector3D[2][4];
        for (int i=0; i<2; i++) {
            for(int j=0; j<4; j++) displacements[i][j] = Vector3D.ZERO;
        }
        for (Vector3D x : T) {
            int x1Index = X1.indexOf(x);
            int x2Index = X2.indexOf(x);
            if (x1Index != -1) {
                displacements[0][x1Index] = new Vector3D(hT - x
                        .dotProduct(nT), nT);
                if (apply) r1.translateVertex(x1Index,
                        displacements[0][x1Index]);
            }
            else {
                displacements[1][x2Index] = new Vector3D(hT - x
                        .dotProduct(nT), nT);
                if (apply) r2.translateVertex(x2Index,
                        displacements[1][x2Index]);
            }

        }
        return new Vector3D[0];
    }

    /**
     * Function to compute the displacements of the two tetrahedra in order
     * to resolve the overlap. The original algorithm quite often results in
     * irregular tetrahedra shapes as part of this step, which gets fixed in
     * the 'concur' part of the algorithm. However, since this implementation
     * doesn't involve the 'concur' part, the sum of displacements for all
     * the vertices of the tetrahedra that correspond to the best
     * least-squares plane is applied to all of its vertices instead.
     * @param apply Flag to apply the displacements directly to the
     *              tetrahedra. Mainly useful for debugging purposes.
     * @return Computed displacements for the two tetrahedra as an array
     * of vectors.
     */
    public Vector3D[] computeDisplacements(boolean apply) {
        if (T == null || nT == null || !Double.isFinite(hT)) throw new
                RuntimeException("Overlap hasn't been resolved yet! Resolve " +
                "overlap and then make displacements");

        Vector3D[] displacements = new Vector3D[2];
        for (int i=0; i<2; i++) displacements[i] = Vector3D.ZERO;

        for (Vector3D x : T) {
            int x1Index = X1.indexOf(x);
            int x2Index = X2.indexOf(x);
            if (x1Index != -1) displacements[0] = displacements[0].add(new Vector3D(hT - x
                        .dotProduct(nT), nT));
            else displacements[1] = displacements[1].add(new Vector3D(hT - x
                        .dotProduct(nT), nT));
        }
        if (apply) {
            r1.translateAllVertices(displacements[0]);
            r2.translateAllVertices(displacements[1]);
        }
        return displacements;
    }


    /**
     * Function to get hT.
     * @return hT.
     */
    public double getHT() {
        return hT;
    }

    /**
     * Function to get nT.
     * @return nT.
     */
    public Vector3D getNT() {
        return nT;
    }

    /**
     * Function to get T.
     * @return T.
     */
    public List<Vector3D> getT() {
        return T;
    }

    /**
     * Function to get R1.
     * @return R1.
     */
    public RegularTetrahedron getR1() {
        return r1;
    }

    /**
     * Function to get R2.
     * @return R2.
     */
    public RegularTetrahedron getR2() {
        return r2;
    }

    /**
     * Function to set R1R2.
     */
    public void setR1R2(RegularTetrahedron rt1, RegularTetrahedron rt2) {
        r1 = rt1;
        r2 = rt2;
        updateSubsets();
    }

    /**
     * Function to set R1.
     */
    public void setR1(RegularTetrahedron rt1) {
        r1 = rt1;
        updateSubsets();
    }

    /**
     * Function to set R2.
     */
    public void setR2(RegularTetrahedron rt2) {
        r2 = rt2;
        updateSubsets();
    }

    /**
     * Function to set the threshold used by this class.
     * @param x Desired threshold.
     */
    public void setThreshold(double x) {
        threshold = x;
    }
}