// Meson electroproduction kinematics.
// Original Fortran version written by Garth Huber in 2000.
// Email: Garth.Huber@uregina.ca
// Converted to C++ by Wenliang (Bill) Li in 2017.
// Email: wenliang.li@msstate.edu

#ifndef KINEMATICS_H_INCLUDED
#define KINEMATICS_H_INCLUDED

#include <cstddef>

// Units throughout this program:
//   masses and energies: MeV, unless the name ends in _g (GeV)
//   momentum components: MeV/c, with c = 1 in the calculations
//   q_2_g, t_g, u_g: GeV^2; t and u: MeV^2
//   trigonometric angles: radians; scan angles and printed angles: degrees
// CM denotes the photon-proton center-of-mass frame; lab denotes the
// frame in which the target proton is initially at rest.

// Scan degree-valued proton-q angles, then calculate and print each row.
void test_out(const double array[], std::size_t limit);
void Cal_print();

// Rest masses in MeV. The _g proton mass is used in GeV-based formulas.
const double m_e = 0.511;
const double m_pi = 139.57;
const double m_pi_g = m_pi / 1000.0;
const double m_p = 938.27;
const double m_p_g = m_p / 1000.0;
const double PI = 3.14159265358979323846;

// Incoming electron kinetic energies in MeV. Add entries to scan more beams;
// main() derives the element count automatically rather than assuming a size.
const int beam_energy[] = {10730};

// Angle conversions occur at the boundary between degree-valued settings
// or output and the radian-valued C++ trigonometric functions.
inline double Rad_to_Deg(double radian) {
    return radian * 180.0 / PI;
}

inline double Deg_to_Rad(double degree) {
    return degree * PI / 180.0;
}

// Historical angular grids retained as alternative configurations.
// Every entry is a proton lab angle relative to q, measured in degrees.
const double theqp_Q1p6_low_eps[] = {-9.534, -10.5, -12.54};
const double theqp_Q1p6_high_eps[] = {-13.29, -10.55, -16.29};
const double theqp_Q2p45_low_eps[] = {-9.19, -10.54, -12.2};
const double theqp_Q2p45_high_eps[] = {-13.61, -10.61, -16.61};
const double theqp_generic[] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20};

#endif // KINEMATICS_H_INCLUDED
