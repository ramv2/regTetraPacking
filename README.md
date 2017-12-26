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
