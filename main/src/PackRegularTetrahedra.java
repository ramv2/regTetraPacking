package src;

import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.util.random.IRandom;
import etomica.util.random.RandomMersenneTwister;
import etomica.util.random.RandomNumberGeneratorUnix;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import java.util.*;

/**
 * Class to pack regular tetrahedra as densely as possible. Uses hard walls
 * instead of periodic boundary conditions. The packing algorithm is
 * summarized as follows:
 *
 * 1) Obtain input 'N' denoting the number of regular tetrahedra to pack,
 * from the user.
 * 2) Initialize the possible length of edges to an array containing few values:
 * [0.01, 0.2, 0.1, 1, 2, 4, 8, 16, 100].
 * 3) Loop through the array of edge lengths and let the current edge length
 * be l. Do the following:
 *      a) Randomize the positions of N regular tetrahedra of length l in 3D
 *      space.
 *      b) Fix the maximum number of steps needed for convergence.
 *      c) Repeated pair-wise overlap detection and resolution based on
 *      computing the divide projections for both tetrahedra. See
 *      {@linkplain OverlapDetectionAndResolution}
 *      d) After convergence, loop through all position vectors to find the
 *      range (i.e., max - min) of x-, y-, z- coordinates among all
 *      tetrahedra vertices. This gives us the basis vectors and the volume
 *      of the cell.
 *      e) Evaluate the density and store the position vectors of all
 *      tetrahedra as well as the basis vectors.
 * 4) Choose the l that results in the highest packing density as the result.
 *
 * @author Ram
 */
public class PackRegularTetrahedra {

    // Number of regular tetrahedra to be packed.
    public int n;

    // Maximum number of steps for convergence.
    public int maxSteps = 10000;

    // Array of tetrahedra.
    public RegularTetrahedron[] rts, bestPackingRTs;

    // Overlap detection and resolution tool.
    public OverlapDetectionAndResolution od;

    // Maximum density possible for a given N.
    public double maxDensity = -1.0;

    // Unit cell vectors corresponding to maximum density.
    public Vector3D[] unitCell;

    // Vertex positions corresponding to maximum density.
    public Vector3D[][] vertexPositions;

    // Edge length corresponding to maximum density.
    public double maxEdgeLength = Double.NaN;

    /**
     * Constructor that takes the number of regular tetrahedra to pack as an
     * argument.
     * @param N Desired number.
     */
    public PackRegularTetrahedra(int N) {
        if (N < 2) throw new IllegalArgumentException("N must be at least 2");
        n = N;
        rts = new RegularTetrahedron[n];
        od = new OverlapDetectionAndResolution();
    }

    /**
     * Function to randomize the positions of N regular tetrahedra. The
     * algorithm proceeds by generating a uniform random number on the
     * surface of a sphere (equivalent to generating a random unit vector
     * uniformly). Uses this as axis for the random rotation as well as the
     * basis for the random translation.
     * @param length Edge length of the regular tetrahedra.
     */
    private void randomizePositions(double length) {
        Space space = Space3D.getInstance();
        int[] seeds = RandomNumberGeneratorUnix.getRandSeedArray();
//        System.out.print("[");
//        for (int s : seeds) System.out.print(s+", ");
//        System.out.print("]");
//        System.out.println();
        IRandom r = new RandomMersenneTwister(seeds);
        for (int i=0; i<n; i++) {
            // Create regular tetrahedra.
            rts[i] = new RegularTetrahedron(length);

            // Rotation part.
            Vector ax = space.makeVector();
            ax.setRandomSphere(r);
            Vector3D axis = new Vector3D(ax.toArray());
            double angle = 2 * Math.PI * r.nextDouble();

            // Perform rotation.
            rts[i].rotate(angle, axis);

            // Translation part.
            Vector displacement = space.makeVector();
            displacement.setRandomSphere(r);
            Vector3D disp = new Vector3D(displacement.toArray());
            double stepSize = length * r.nextDouble();
            disp = disp.scalarMultiply(stepSize);

            // Perform translation.
            rts[i].translateAllVertices(disp);
        }
    }

    /**
     * Function to pack N regular tetrahedra as densely as possible. Choose
     * edge lengths of different scales to try. For a fixed length, the
     * algorithm proceeds until it is converged or has reached maximum number
     * of steps. If the number of overlaps has either increased or stayed the
     * same compared to the previous iteration, the algorithm performs a
     * random move (50% probability of rotation and translation). This helps
     * speed up the convergence rate. If it has converged, we compute the
     * unit cell parameters, density and store the vertex positions.
     */
    public void pack() {
        double[] lengths = new double[] {0.01, 0.2, 0.1, 1, 2, 4};
        for (double l : lengths) {
            randomizePositions(l);
            boolean isConverged = false;
            int steps = 0;
            double density = Double.POSITIVE_INFINITY;
            int prevOverlaps = 0;
            while (steps++ <= maxSteps) {
                // Compute pairs of overlapping tetrahedra.
                Map<Integer, List<Integer>> overlaps = getCurrentOverlaps();
                int currentOverlaps = overlaps.size();
                if (currentOverlaps == 0) {
                    isConverged = true;
                    break;
                }
                for (Map.Entry<Integer, List<Integer>> entry : overlaps
                        .entrySet()) {

                    int index1 = entry.getKey();
                    for (int index2 : entry.getValue()) {
                        if (prevOverlaps >= currentOverlaps) {
                            performRandomMove(index2);
                        }
                        od.setR1R2(rts[index1], rts[index2]);
                        if (od.isOverlap()) od.findDisplacementsToResolve(true);
                    }
                }
                prevOverlaps = currentOverlaps;
            }
            if (isConverged) {
                unitCell = computeUnitCellParams();
                density = computePackingDensity(rts[0].getVolume());
                if (density > maxDensity) {
                    maxEdgeLength = l;
                    maxDensity = density;
                    vertexPositions = getPositions();
                }
            }
        }
        if (maxDensity != 0.0) {
            System.out.println("Successfully converged!");
            printSolution();
        }
        else {
            System.out.println("Failed to converge for all l!");
        }
    }

    /**
     * Performs a random translation or rotation of the regular tetrahedron
     * with the specified index. This is used to improve the convergence rate.
     * @param index2 Index of tetrahedron.
     */
    private void performRandomMove(int index2) {
        IRandom r = new RandomMersenneTwister(RandomNumberGeneratorUnix
                .getRandSeedArray());
        Space space = Space3D.getInstance();

        // Perform rotation half the time.
        if (r.nextInt(1) == 0) {
            // Rotation part.
            Vector ax = space.makeVector();
            ax.setRandomSphere(r);
            Vector3D axis = new Vector3D(ax.toArray());
            double angle = 2 * Math.PI * r.nextDouble();

            // Perform rotation.
            rts[index2].rotate(angle, axis);
        }
        else {
            // Perform translation the other half.
            // Translation part.
            Vector displacement = space.makeVector();
            displacement.setRandomSphere(r);
            Vector3D disp = new Vector3D(displacement.toArray());
            double stepSize = rts[index2].getEdgeLength() * r.nextDouble() /
                    20;
            disp = disp.scalarMultiply(stepSize);
            rts[index2].translateAllVertices(disp);
        }
    }

    /**
     * Function to get the unit cell basis vectors.
     * @return Unit cell basis vectors.
     */
    public Vector3D[] getUnitCell() {
        return unitCell;
    }

    /**
     * Function to get the tetrahedra vertex positions.
     * @return Vertex positions of all tetrahedra.
     */
    public Vector3D[][] getVertexPositions() {
        return vertexPositions;
    }

    /**
     * Function to get the maximum density.
     * @return Maximum density.
     */
    public double getMaxDensity() {
        return maxDensity;
    }

    /**
     * Function to print the solution in the desired format.
     */
    private void printSolution() {
        System.out.println("Solution: Packed "+n+" regular tetrahedra with " +
                "edge length: "+maxEdgeLength+" (arbitrary length units) and " +
                "packing density of: " +maxDensity);
        System.out.println("Basis:");
        for (Vector3D v : unitCell) System.out.println(v.toString());
        System.out.println("Tetrahedra positions:");
        for (int i=0; i<n; i++) {
            for (Vector3D v: vertexPositions[i]) System.out.print(v.toString
                    ()+" ");
            System.out.println();
        }
    }

    /**
     * Function to store the vertex positions of all tetrahedra.
     * @return Vertex positions of all tetrahedra.
     */
    private Vector3D[][] getPositions() {
        vertexPositions = new Vector3D[n][4];
        for (int i=0; i<n; i++) vertexPositions[i] = rts[i].getAllVertices();
        return vertexPositions;
    }

    /**
     * Function to compute the packing density, given the volume of a single
     * tetrahedron.
     * @param tetVolume Volume of a single tetrahedron.
     * @return Packing density.
     */
    private double computePackingDensity(double tetVolume) {
        double cellVolume = unitCell[0].getNorm() * unitCell[1].getNorm() *
                unitCell[2].getNorm();
        return n * tetVolume / cellVolume;
    }

    /**
     * Function to compute unit cell basis vectors.
     * @return Unit cell basis vectors.
     */
    private Vector3D[] computeUnitCellParams() {
        double minX = Double.POSITIVE_INFINITY;
        double minY = Double.POSITIVE_INFINITY;
        double minZ = Double.POSITIVE_INFINITY;
        double maxX = Double.NEGATIVE_INFINITY;
        double maxY = Double.NEGATIVE_INFINITY;
        double maxZ = Double.NEGATIVE_INFINITY;
        for (int i=0; i<n; i++) {
            Vector3D[] vertices = rts[i].getAllVertices();
            for (int j=0; j<4; j++) {
                double x = vertices[j].getX();
                double y = vertices[j].getY();
                double z = vertices[j].getZ();
                if (x < minX) minX = x;
                if (x > maxX) maxX = x;
                if (y < minY) minY = y;
                if (y > maxY) maxY = y;
                if (z < minZ) minZ = z;
                if (z > maxZ) maxZ = z;
            }
        }
        Vector3D[] lattice = new Vector3D[3];
        lattice[0] = new Vector3D(maxX - minX, Vector3D.PLUS_I);
        lattice[1] = new Vector3D(maxY - minY, Vector3D.PLUS_J);
        lattice[2] = new Vector3D(maxZ - minZ, Vector3D.PLUS_K);
        return lattice;
    }

    /**
     * Function to compute the map containing the index of one tetrahedron as
     * the key and a list of indices of the tetrahedra with which it
     * overlaps, as the value.
     * @return Map of overlaps.
     */
    private Map<Integer, List<Integer>> getCurrentOverlaps() {
        Map<Integer, List<Integer>> overlaps = new TreeMap<>();
        for (int i=0; i<n; i++) {
            for (int j=0; j<i; j++) {
                od.setR1R2(rts[i], rts[j]);
                if (od.isOverlap()) {
                    if (!overlaps.containsKey(j)) overlaps.put(j, new
                            ArrayList<>());
                    overlaps.get(j).add(i);
                }
            }
        }
        return overlaps;
    }

    /**
     * Function to set the number of maximum steps to be used in the algorithm.
     * @param mSteps Desired number. Default is 10,000.
     */
    public void setMaxSteps(int mSteps) {
        maxSteps = mSteps;
    }

    public static void main(String[] args) {
        int N = 20;
        if (args.length >= 1) {
            try {
                N = Integer.parseInt(args[0]);
            } catch (Exception e) {
                System.out.println("First argument must be of type integer.");
                System.out.println("Please enter a valid integer. Or to run with " +

                        "a default value of 20, simply run it without any " +
                        "arguments.");
                System.exit(1);
            }
        }
        PackRegularTetrahedra p = new PackRegularTetrahedra(N);
        p.pack();

    }
}