package src;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

/**
 * Class that defines a regular tetrahedron using 4 vectors to denote the
 * positions of its vertices in 3D space. Always places the vertices in the
 * same exact location chosen arbitrarily as specified here:
 * http://mathworlvertices[3].wolfram.com/RegularTetrahedron.html
 * Users can modify the orientation of it a later time by performing
 * rotations or translations as needed.
 *
 * @author Ram
 */
public class RegularTetrahedron {

    // Denotes the edge length of the regular tetrahedron.
    public double edgeLength;

    // Array of position vectors for the four vertices.
    public Vector3D[] vertices = new Vector3D[4];

    // Threshold used to check if this is a regular tetrahedra or not.
    public double threshold = 1e-20;

    /**
     * Empty constructor to initialize a regular tetrahedron with edge length
     * 1.0.
     */
    public RegularTetrahedron() {
        this(1.0);
    }

    /**
     * Constructs a regular tetrahedron with given edge length.
     * @param length Desired edge length.
     */
    public RegularTetrahedron(double length) {
        if (length <= 0) throw new IllegalArgumentException("Edge length must" +
                " be positive");

        edgeLength = length;

        vertices[0] = new Vector3D(new double[]{Math.sqrt(3) / 3.0, 0, 0});
        vertices[0] = vertices[0].scalarMultiply(edgeLength);
        vertices[1] = new Vector3D(new double[]{-Math.sqrt(3) / 6.0, -0.5, 0});
        vertices[1] = vertices[1].scalarMultiply(edgeLength);
        vertices[2] = new Vector3D(new double[]{-Math.sqrt(3) / 6.0, 0.5, 0});
        vertices[2] = vertices[2].scalarMultiply(edgeLength);
        vertices[3] = new Vector3D(new double[]{0, 0, Math.sqrt(6) / 3.0});
        vertices[3] = vertices[3].scalarMultiply(edgeLength);
    }

    /**
     * Function to translate all the vertices by the same amount.
     * @param disp Desired displacement vector.
     */
    public void translateAllVertices(Vector3D disp) {
        for (int i=0; i<4; i++) translateVertex(i, disp);
    }

    /**
     * Function to translate a specific vertex of the regular tetrahedron along
     * the given displacement vector.
     * @param index Index of the vertex to be translated.
     * @param displacement Vector defining the amounts by which to translate.
     */
    public void translateVertex(int index, Vector3D displacement) {
        if (displacement.isInfinite() || displacement.isNaN()) throw new
                IllegalArgumentException("Vector for translation must be " +
                "finite.");
        if (index < 0 || index > 3) throw new IllegalArgumentException("Index" +
                " of the vertex to be translated should be between [0, 3]!");
        vertices[index] = vertices[index].add(displacement);
    }

    /**
     * Function to rotate all the vertices about a given axis by a given
     * angle. Code adapted from molecular simulation package Etomica:
     * https://github.com/etomica/etomica/
     * @param angle Angle (in radians) by which to rotate.
     * @param axis Vector denoting the axis.
     */
    public void rotate(double angle, Vector3D axis) {
        // consider a circle on the surface of the unit sphere.  The given axis
        // passes through the center of the circle.  The circle passes through
        // the current direction vector and the vector v4 defined below.  We
        // rotate the direction by the given angle (dt) around the circle.

        // v1 is the projection of direction onto axis
        // v2 is the component of direction perpendicular to axis
        // v3 has the same magnitude as v2 and is perpendicular to both
        //    direction and axis
        // v4 is a unit vector whose components are v1 and v3

        // v1 = v1overAxis * axis
        double cdt = Math.cos(angle);
        double sdt = Math.sin(angle);
        for (int i=0; i<4; i++) {
            double v1overAxis = axis.dotProduct(vertices[i]);
            Vector3D temp = new Vector3D(-v1overAxis, axis);
            temp = temp.add(vertices[i]);
            // now temp = v2
            Vector3D temp2 = axis;
            temp2 = temp2.crossProduct(vertices[i]);
            // now temp2 = v3
            vertices[i] = new Vector3D(cdt, temp);
            vertices[i] = vertices[i].add(new Vector3D(sdt, temp2));
            vertices[i] = vertices[i].add(new Vector3D(v1overAxis, axis));
        }
    }

    /**
     * Function to get a copy of the position vectors of all vertices in this
     * tetrahedron instance.
     * @return A copy of the variable vertices.
     */
    public Vector3D[] getAllVertices() {
        Vector3D[] v = new Vector3D[4];
        for (int i=0; i<4; i++) v[i] = new Vector3D(vertices[i].toArray());
        return v;
    }

    /**
     * Function to compute the volume of the regular tetrahedron.
     * @return Volume.
     */
    public double getVolume() {
        return edgeLength * edgeLength * edgeLength * Math.sqrt(2) / 12.0;
    }

    /**
     * Function to set the positions of all the vertices.
     * @param pos Desired positions as a 4x3 double array.
     */
    public void setAllVertices(double[][] pos) {
        if (pos.length <= 0) throw new IllegalArgumentException("Input array " +
                "can't be empty!!");
        if (pos.length != 4 || pos[0].length != 3) throw new
                IllegalArgumentException("Expected 4x3 array. Received "+pos
                .length+"x"+pos[0].length+" instead!");
        for (int i=0; i<4; i++) vertices[i] = new Vector3D(pos[i]);
        computeEdgeLength();
    }

    private void computeEdgeLength() {
        edgeLength = vertices[1].subtract(vertices[0]).getNorm();
    }

    /**
     * Function to check if the current regular tetrahedron instance is
     * regular (i.e., all the edges have the same length).
     * @return Whether the current instance is a regular tetrahedron or not.
     */
    public boolean isRegular() {
        // TODO: More rigorous tests to check if it is regular.
        // Get the three edges from vertex A.
        Vector3D AB = vertices[1].subtract(vertices[0]);
        Vector3D AC = vertices[2].subtract(vertices[0]);
        Vector3D AD = vertices[3].subtract(vertices[0]);

        // Check if they are coplanar.
        double dp = AB.dotProduct(AC.crossProduct(AD));
        if ((dp * dp) <= threshold) return false;

        // Check if all edge lengths are the same.
        double l01 = vertices[0].subtract(vertices[1]).getNorm();
        double l02 = vertices[0].subtract(vertices[2]).getNorm();
        double l03 = vertices[0].subtract(vertices[3]).getNorm();
        double l12 = vertices[1].subtract(vertices[2]).getNorm();
        double l13 = vertices[1].subtract(vertices[3]).getNorm();
        double l23 = vertices[2].subtract(vertices[3]).getNorm();

        double diffSqSum = (l01 - l02) * (l01 - l02) + (l02 - l03) * (l02 -
                l03) + (l03 - l12) * (l03 - l12) + (l12 - l13) * (l12 - l13)
                + (l13 - l23) * (l13 - l23);
        return (diffSqSum <= threshold) && ((edgeLength - l01) * (edgeLength -
                l01) <= threshold);
    }

    /**
     * Function to check if the current tetrahedron instance is equal to
     * another instance.
     * @param o Another tetrahedron to compare.
     * @return True if all the vertices of the given tetrahedron match with
     * the current instance, else, false.
     */
    @Override
    public boolean equals(Object o) {
        if (o instanceof RegularTetrahedron) {
            for (int i=0; i<4; i++) if (!vertices[i].equals((
                    (RegularTetrahedron) o).vertices[i])) return false;
            return true;
        }
        return false;
    }

    /**
     * Function to get the edge length of this regular tetrahedron instance.
     * @return Edge length.
     */
    public double getEdgeLength() {
        return edgeLength;
    }
}