#include "grissom_custom_quadrature.h"
#include "ChiLua/chi_lua.h"

#include <fstream>
#include <iomanip>

#include "chi_log.h"
extern ChiLog& chi_log;

#include "ChiMath/chi_math.h"
extern ChiMath& chi_math_handler;

void chi_math::GrissomCustomQuadrature::BuildDiscreteToMomentOperator(int scatt_order, bool oneD)
{
    chi_log.Log() << "Hello";

    int num_angles = abscissae.size();
    int num_moms = 0;
    float a;

    d2m_op.clear();

    std::ifstream file;
    file.open(filename);
    if (not file.is_open())
    {
        chi_log.Log(LOG_ALLERROR) << "Failed to open Quadrature File in call to " << __FUNCTION__ << ".";
        exit(EXIT_FAILURE);
    }

    int mc=-1; //moment count
    for (int ell=0; ell<=scatt_order; ell++)
    {
        for (int m = -ell; m <= ell; m++)
        {
            std::vector<double> cur_mom; mc++;
            num_moms++;

            for (int n = 0; n < num_angles; n++)
            {
                file >> a;
                double value = a;
                cur_mom.push_back(value);
//                chi_log.Log() << value;
            }

            d2m_op.push_back(cur_mom);
        }//for m
    }//for ell

    std::stringstream outs;
    outs
            << "\nQuadrature d2m operator:\n";
    for (int n=0; n<num_angles; n++)
    {
        outs << std::setw(5) << n;
        for (int m=0; m<num_moms; m++)
        {
            outs
                    << std::setw(15) << std::left << std::fixed
                    << std::setprecision(10) << d2m_op[m][n] << " ";
        }
        outs << "\n";
    }
    chi_log.Log() << outs.str();

    file.close();

    chi_log.Log() << "Goodbye";
}
void chi_math::GrissomCustomQuadrature::BuildMomentToDiscreteOperator(int scatt_order, bool oneD)
{
    chi_log.Log() << "Hello";

    int num_angles = abscissae.size();
    int num_moms = 0;
    float a;

    m2d_op.clear();

    std::ifstream file;
    file.open(filename);
    if (not file.is_open())
    {
        chi_log.Log(LOG_ALLERROR) << "Failed to open Quadrature File in call to " << __FUNCTION__ << ".";
        exit(EXIT_FAILURE);
    }

    int mc=-1; //moment count
    for (int ell=0; ell<=scatt_order; ell++)
    {
        for (int m = -ell; m <= ell; m++)
        {
            std::vector<double> cur_mom; mc++;
            num_moms++;

            for (int n = 0; n < num_angles; n++)
            {
                file >> a;
                double value = a;
                cur_mom.push_back(value);
//                chi_log.Log() << value;
            }

            m2d_op.push_back(cur_mom);
        }//for m
    }//for ell

    std::stringstream outs;
    outs
            << "\nQuadrature d2m operator:\n";
    for (int n=0; n<num_angles; n++)
    {
        outs << std::setw(5) << n;
        for (int m=0; m<num_moms; m++)
        {
            outs
                    << std::setw(15) << std::left << std::fixed
                    << std::setprecision(10) << m2d_op[m][n] << " ";
        }
        outs << "\n";
    }
    chi_log.Log() << outs.str();

    file.close();

    chi_log.Log() << "Goodbye";
}

// #########################################################################################
// #########################################################################################
// #########################################################################################

int chiCreateGrissomCustomQuadrature(lua_State* L)
{

//    printf("Hello World");

    int num_args = lua_gettop(L);
    if (num_args != 1)
        LuaPostArgAmountError(__FUNCTION__ , 1, num_args);

    const char* file_name_raw = lua_tostring(L, 1);
    std::string file_name(file_name_raw);
    std::ifstream file;
    file.open(file_name);
    if (not file.is_open())
    {
        chi_log.Log(LOG_ALLERROR) << "Failed to open Quadrature File in call to " << __FUNCTION__ << ".";
        exit(EXIT_FAILURE);
    }

    auto new_quad = std::make_shared<chi_math::GrissomCustomQuadrature>(file_name);
    chi_math::QuadraturePointPhiTheta q_points;

    float a;
    while(file >> a)
    {
        q_points.phi = a;
        file >> a;
        q_points.theta = a;
//        chi_log.Log() << "Here";
//        chi_log.Log() << q_points.theta;
        double x = sin(q_points.theta)*cos(q_points.phi);
        double y = sin(q_points.theta)*sin(q_points.phi);
        double z = cos(q_points.theta);

        new_quad->abscissae.push_back(q_points);
        new_quad->weights.push_back(1.0);
        new_quad->omegas.emplace_back(x, y, z);
    }
//    q_points.phi = 0.5;
//    q_points.theta = 0.5;
//    double x = sin(q_points.theta)*cos(q_points.phi);
//    double y = sin(q_points.theta)*sin(q_points.phi);
//    double z = cos(q_points.theta);
//
//    new_quad->abscissae.push_back(q_points);
//    new_quad->weights.push_back(1.0);
//    new_quad->omegas.emplace_back(x, y, z);

//    std::string line;
//    bool not_eof = bool(std::getline(file, line));
//
//    while (not_eof)
//    {
//        chi_log.Log() << line << "\n";
//        not_eof = bool(std::getline(file, line));
//    }

    file.close();

    chi_math_handler.angular_quadratures.push_back(new_quad);
    lua_pushnumber(L, chi_math_handler.angular_quadratures.size()-1);

    return 1;

}