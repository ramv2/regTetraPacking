package src;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

/**
 * @author Ram
 */
public class PackingRegularTetrahedraTest {
    @Test
    public void testPack() {
        // Simple test.
        int N = 20;
        PackRegularTetrahedra p = new PackRegularTetrahedra(N);
        p.pack();


        Vector3D[][] vPos = p.getVertexPositions();
        RegularTetrahedron[] rt = new RegularTetrahedron[N];
        for (int i=0; i<N; i++) {
            rt[i] = new RegularTetrahedron();
            double[][] riPos = new double[4][3];
            for (int k=0; k<4; k++) riPos[k] = vPos[i][k].toArray();
            rt[i].setAllVertices(riPos);
            assertEquals(p.maxEdgeLength, rt[i].getEdgeLength(), 1e-10);
        }

        OverlapDetectionAndResolution odTool = new
                OverlapDetectionAndResolution();
        for (int i=0; i<N; i++) {
            for (int j=0; j<N; j++) {
                if (i != j) {
                    odTool.setR1R2(rt[i], rt[j]);
                    odTool.updateSubsets();
                    assertFalse(odTool.isOverlap());
                }
            }
        }
    }
}