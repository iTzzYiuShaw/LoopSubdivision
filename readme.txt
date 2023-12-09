To compile on feng-linux / feng-gps:

module purge
module add legacy-eng
module load qt/5.15.2
qmake -project QT+=opengl
qmake
make

To run the project, please input command AS FOLLOW:
./LoopSubdivisionRelease ../path_to/model.diredgenormal loopTimes

The EXAMPLE of input command is (loop subdivision for 3 times):
./LoopSubdivisionRelease ../LoopSubdivisionRelease/diredgenormals/cube.diredgenormal 3


After the project runs, a path that indicates the output file is shown in the command line as follows:

[sc22y2x@feng-linux-07 LoopSubdivisionRelease]$ ./LoopSubdivisionRelease ../LoopSubdivisionRelease/diredgenormals/cube.diredgenormal 3
filePath: ../LoopSubdivisionRelease/diredgenormals/cube_editted.diredgenormal


