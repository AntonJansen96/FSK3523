#include "md.h"

void MD::writeFrame(size_t step) const
{
    auto file = fopen("traj.pdb", "a");

    // Print the header.
    fprintf(file, "MODEL        %zu\n", step);

    if (d_usePBC)
        fprintf(file, "CRYST1  %7.3f  %7.3f  %7.3f  90.00  90.00  90.00 P 1           1\n", 10 * d_boxsize[0], 10 * d_boxsize[1], 10 * d_boxsize[2]);

    for (Atom const &atom : d_AtomList)
        fprintf(file, "%-6s%5zu %-4s %-4s %4zu      %-8.3f%-8.3f%-8.3f\n", "ATOM", atom.idx, atom.name.c_str(), atom.name.c_str(), atom.idx, 10 * atom.x[0], 10 * atom.x[1], 10 * atom.x[2]);

    fprintf(file, "TER\nENDMDL\n");
    fclose(file);
}
