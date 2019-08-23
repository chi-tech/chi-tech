#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <iostream>
#include <fstream>
#include <cmath>

using namespace dealii;

void first_grid()
{
    Triangulation<2> triangulation;

    GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(4);

    std::ofstream out("grid-1.eps");
    GridOut       grid_out;
    grid_out.write_eps(triangulation, out);
    std::cout << "Grid written to grid-1.eps" << std::endl;
}


void second_grid()
{
    Triangulation<2> triangulation;

    const Point<2> center(1, 0);
    const double   inner_radius = 0.5, outer_radius = 1.0;
    GridGenerator::hyper_shell(
            triangulation, center, inner_radius, outer_radius, 10);

    for (unsigned int step = 0; step < 5; ++step)
    {

        for (auto cell : triangulation.active_cell_iterators())
        {

            for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_cell; ++v)
            {

                const double distance_from_center =
                        center.distance(cell->vertex(v));

                if (std::fabs(distance_from_center - inner_radius) < 1e-10)
                {
                    cell->set_refine_flag();
                    break;
                }
            }
        }

        triangulation.execute_coarsening_and_refinement();
    }


    std::ofstream out("grid-2.eps");
    GridOut       grid_out;
    grid_out.write_eps(triangulation, out);

    std::cout << "Grid written to grid-2.eps" << std::endl;

    triangulation.reset_manifold(0);

}



// @sect3{The main function}

// Finally, the main function. There isn't much to do here, only to call the
// two subfunctions, which produce the two grids.
int dealtest()
{
    first_grid();
    second_grid();
}
