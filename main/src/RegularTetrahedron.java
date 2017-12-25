import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

import java.util.ArrayList;
import java.util.List;

/**
 * Class that defines a regular tetrahedron using 4 vectors to denote the
 * positions of its vertices in 3D space. Always places the vertices in the
 * same exact location chosen arbitrarily as specified here:
 * http://mathworlvertices[3].wolfram.com/RegularTetrahedron.html
 * Users can modify the orientation of it a later time by performing
 * rotations or translations as needevertices[3].
 *
 */
public class RegularTetrahedron {

    /**
     * Denotes the edge length of the regular tetrahedron.
     */
    public double edgeLength;

    /**
     * Array of position vectors for the four vertices.
     */
    public Vector3D[] vertices = new Vector3D[4];

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
        vertices[0].scalarMultiply(edgeLength);
        vertices[1] = new Vector3D(new double[]{-Math.sqrt(3) / 6.0, -0.5, 0});
        vertices[1].scalarMultiply(edgeLength);
        vertices[2] = new Vector3D(new double[]{-Math.sqrt(3) / 6.0, 0.5, 0});
        vertices[2].scalarMultiply(edgeLength);
        vertices[3] = new Vector3D(new double[]{0, 0, Math.sqrt(6) / 3.0});
        vertices[3].scalarMultiply(edgeLength);
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

    public void rotate(double angle, Vector3D axis) {

    }

    /**
     * Function to get the position vectors of a specific vertex in this
     * tetrahedron instance.
     * @return The position vector of the vertex at the desired index.
     */
    public Vector3D getVertex(int index) {
        if (index < 0 || index > 3) throw new IllegalArgumentException("Index" +
                " of the vertex to be translated should be between [0, 3]!");
        return vertices[index];
    }

    /**
     * Function to get a copy of the position vectors of all vertices in this
     * tetrahedron instance.
     * @return A copy of the variable vertices.
     */
    public Vector3D[] getAllVertices() {
        return vertices.clone();
    }

    public double getVolume() {
        return edgeLength * edgeLength * edgeLength * Math.sqrt(2) / 12.0;
    }
}