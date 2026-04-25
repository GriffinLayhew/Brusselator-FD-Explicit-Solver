To run the code you can do the following:

1) make run
    -This will use the default number of processors that are set in the make file

2) mpirun -np <np> bin/Program


NOTES:
-----
a) Note that many of my make files include OpenMP. This is because I am using a makefile from my research.
b) I would like to learn more on how to order output. This was not done here, but could be added in the future.
c) I make a note on to how I think the ring could be better using modulo arithmetic. This is a possible future improvment.
