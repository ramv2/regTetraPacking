<h2><u>Regular tetrahedron packing</u></h2>
An optimization scheme for packing regular tetrahedra. Implements the divide projection of the divide and concur algorithm discussed in detail here:
Y. Kallus, V. Elser, and S. Gravel, "Method for dense packing discovery", Physical Review E, vol. 82, no. 5, Nov. 2010.

Dependencies: 
1. org.apache.commons:commons-math3:3.6.1 (can be downloaded from Maven repository)
2. etomica.jar (see: https://github.com/etomica/etomica)
3. JUnit4 (can be downloaded from Maven repository)

Usage:
java -cp <path_to_regTetraPacking.jar>:<path_to_commons-math3-3.6.1.jar>:<path_to_etomica.jar> src.PackRegularTerahedra <int_value>

Here <...> denotes user input. For instance, the user must replace <path_to_commons-math3-3.6.1.jar> with the path to commons-math3-3.6.1.jar in her/his local machine. <int_value> denotes any input in the integer format.

<h2><u>Divide and Concur</u></h2>
My implementation of the packing algorithm is based on the divide and concur algorithm reported by Kallus et al. [1]. The main idea behind the divide and concur algorithm is to solve overlaps pair-wise and combine them efficiently to obtain a packing in which no two shapes overlap. The divide projection is a two-step process involving overlap detection and resolution.

Overlap detection:
In a d-dimensional space, let the set of vertices of the two tetrahedra be X1 and X2 respectively. Let S denote all the subsets of X1 U X2 (here U denotes the union operation) of size d that contain at least one vertex from each tetrahedra. For every set in S, we find the least-squares plane (i.e., a plane that minimizes the sum of squared distances from a set of points to it) and check to see if separates the sets X1\S and X2\S. Here, X1\S and X2\S are the sets that are formed by removing the vertices present in S from X1 and X2 respectively. If there exists such a separating plane, then the two tetrahedra don’t overlap.

Overlap resolution:
If no such separating plane exists, then we can find one. Let S denote all the subsets of X1 U X2 (here U denotes the union operation) of size greater than d that contain at least one vertex from each tetrahedra. For every set in S, we find the least-squares plane (i.e., a plane that minimizes the sum of squared distances from a set of points to it) and check to see if separates the sets X1\S and X2\S. Here, X1\S and X2\S are the sets that are formed by removing the vertices present in S from X1 and X2 respectively. If there exists such a separating plane, then the minimum displacement is set to the sum of squared distances of the set S to the plane. Otherwise, we simply ignore this set S. The overall minimum displacement required to resolve the overlap is obtained by computing the minimum of all such ‘sum of squared distances’.

In the original algorithm, the overlap resolution step usually results in non-congruent tetrahedra. To maintain congruency, we simply translate the entire tetrahedra by a displacement which is equal to the sum of displacements for each of its vertices in the set S that results in the best least-squares plane.

<h2><u>Summary of approach</u></h2>
<ol>
<li> Obtain input 'N' denoting the number of regular tetrahedra to pack, from the user.</li>
<li> Initialize the possible length of edges to an array containing few values:
[0.01, 0.2, 0.1, 1, 2, 4, 8, 16, 100].</li>
<li> Loop through the array of edge lengths and let the current edge length
be l. Do the following:
<ul>
<li> Randomize the positions of N regular tetrahedra of length l in 3D
space.</li>
<li> Repeated pair-wise overlap detection and resolution based on
computing the divide projections for both tetrahedra. See
{@linkplain OverlapDetectionAndResolution}.</li>
<li> If the number of overlaps from the previous step is greater than
or equal to the number of overlaps in the current step, perform a random
rotation or translation of one of the tetrahedra. This is useful to get
out of problematic configurations and improves convergence rate.
<li> After convergence, loop through all position vectors to find the
range (i.e., max - min) of x-, y-, z- coordinates among all
tetrahedra vertices. This gives us the basis vectors and the volume
of the cell. </li>
<li> Evaluate the density and store the position vectors of all
tetrahedra as well as the basis vectors.</li>
</ul>
</li>
<li> Choose the l that results in the highest packing density as the result.</li>
</ol>

<h3><u>References</u></h3>
[1]	Y. Kallus, V. Elser, and S. Gravel, "Method for dense packing discovery", Physical Review E, vol. 82, no. 5, Nov. 2010.
