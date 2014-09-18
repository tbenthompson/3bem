#include "mpi.h"
#include "UnitTest++.h"

int main(int, char const *[])
{
    MPI_Init(NULL, NULL);
    int retval = UnitTest::RunAllTests();
    MPI_Finalize();
    return retval;
}
