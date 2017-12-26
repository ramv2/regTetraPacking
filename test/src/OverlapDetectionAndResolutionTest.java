package src;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.util.random.IRandom;
import etomica.util.random.RandomMersenneTwister;
import etomica.util.random.RandomNumberGeneratorUnix;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * @author Ram
 */
public class OverlapDetectionAndResolutionTest {

    private RegularTetrahedron r1, r2;
    private OverlapDetectionAndResolution odTool;

    @Before
    public void setUp() {
        r1 = new RegularTetrahedron();
        r2 = new RegularTetrahedron();
        odTool = new OverlapDetectionAndResolution(r1, r2);
    }

    @After
    public void tearDown() {
        r1 = null;
        r2 = null;
        odTool = null;
    }

    @Test
    public void testIsOverlap() {
        // Simple cases:
        // Test overlap with itself.
        // Create regular tetrahedron of unit edge length.

        assertTrue(odTool.isOverlap());

        // Translate one of them by a displacement vector whose norm is
        // greater than 1 so that they don't overlap.
        Vector3D displacement = new Vector3D(3, Vector3D.PLUS_I);
        for (int i=0; i<4; i++) r2.translateVertex(i, displacement);
        odTool.updateSubsets();
        assertFalse(odTool.isOverlap());

        // More complicated cases. Set vertices of tetrahedron to make it
        // temporarily irregular. Cases 1 and 2 are taken from the example
        // provided here:
        // https://gist.github.com/postspectacular/9021724
        // Case 3 is from one of the comments provided in that link.
        // Note: The clojure algorithm provided in the link above results in
        // an incorrect prediction. These two tetrahedra do not overlap. All
        // cases were confirmed by plotting the convex hull of these vertices
        // using python's matplotlib library.

        // Case 1:
        double[][] pos1 = new double[][] {{0.0, 0.0, 0.0}, {50.0, 50.0, 0.0},
                {100.0, 0.0, 0.0}, {50.0, 25.0, 100.0}};
        r1.setAllVertices(pos1);
        double[][] pos2 = new double[][] {{0.0, 0.0, 0.0}, {-50.0, 50.0, 0.0},
                {-200.0, 0.0, 20.0}, {150.0, 25.0, 100.0}};
        r2.setAllVertices(pos2);
        odTool.updateSubsets();
        assertTrue(odTool.isOverlap());

        // Case 2:
        pos1 = new double[][] {{0.0, 0.0, 0.0}, {50.0, 50.0, 0.0},
                {100.0, 0.0, 0.0}, {50.0, 25.0, 100.0}};
        r1.setAllVertices(pos1);
        pos2 = new double[][] {{100.001, 0.0, 0.0}, {150.0, 50.0,
                0.0}, {200.0, 0.0, 0.0}, {150.0, 25.0, 10.0}};
        r2.setAllVertices(pos2);
        odTool.updateSubsets();
        assertFalse(odTool.isOverlap());

        // Case 3:
        pos1 = new double[][] {{2328.824, 1256.559, 9889.406}, {2328.824,
                1256.559, 9889.206}, {2208.937, 1187.342, 9965.406},
                {2132.737,1319.324, 9965.406}};
        r1.setAllVertices(pos1);
        pos2 = new double[][] {{2137.673, 1234.186, 10003.130}, {2060.557,
                1281.337, 9973.415}, {2096.147, 1303.468, 9925.334},
                {2101.378, 1213.231, 9973.415}};
        r2.setAllVertices(pos2);
        odTool.updateSubsets();
        assertFalse(odTool.isOverlap());
    }

    @Test
    public void testOverlapResolution() {
        // Test for irregular tetrahedra.
        double[][] pos1 = new double[][] {{0.0, 0.0, 0.0}, {50.0, 50.0, 0.0},
                {100.0, 0.0, 0.0}, {50.0, 25.0, 100.0}};
        r1.setAllVertices(pos1);
        double[][] pos2 = new double[][] {{0.0, 0.0, 0.0}, {-50.0, 50.0, 0.0},
                {-200.0, 0.0, 20.0}, {150.0, 25.0, 100.0}};
        r2.setAllVertices(pos2);
        odTool.updateSubsets();
        assertTrue(odTool.isOverlap());

        odTool.findDisplacementsToResolve(true);
        odTool.updateSubsets();
        assertFalse(odTool.isOverlap());

        // Test for regular tetrahedra.
        r1 = new RegularTetrahedron();
        r2 = new RegularTetrahedron();
        Space space = Space3D.getInstance();
        IRandom r = new RandomMersenneTwister(RandomNumberGeneratorUnix
                .getRandSeedArray());
        Vector ax = space.makeVector();
        ax.setRandomSphere(r);
        Vector3D axis = new Vector3D(ax.toArray());
        double angle = 2 * Math.PI * r.nextDouble();

        // Perform rotation.
        r2.rotate(angle, axis);
        odTool.setR1R2(r1, r2);
        assertTrue(odTool.isOverlap());
        odTool.findDisplacementsToResolve(true);
        odTool.updateSubsets();
        assertFalse(odTool.isOverlap());
    }
}
