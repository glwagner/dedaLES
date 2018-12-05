from mpi4py import MPI
import dedaLES

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

print(rank)
