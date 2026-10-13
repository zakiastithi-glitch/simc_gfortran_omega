// Meson electroproduction kinematics.
// Original Fortran version written by Garth Huber in 2000.
// Email: Garth.Huber@uregina.ca
// Converted to C++ by Wenliang (Bill) Li in 2017.
// Email: wenliang.li@msstate.edu

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <iomanip>
#include <iostream>

#include "kinematics.h"

using namespace std;

// Per-run state shared by main(), test_out(), and Cal_print(). Keep these
// definitions private to this translation unit rather than in the header.
// Each beam overwrites electron state; each scan row overwrites final-state state.
namespace {
// User inputs: mass and W in GeV, positive Q^2 in GeV^2.
double m_meson_g, q_2_g, W_g;
double m_meson, W; // MeV equivalents of mass and W.

// Photon energy and momentum, and the boost along its direction.
double nu_g, nu, p_gamma_g, p_gamma;
double p_gammax, p_gammaz;
double beta_cm, gamma_cm; // Dimensionless boost speed and Lorentz factor.
double p_cm, e_meson_cm, e_proton_cm;

// Incoming and scattered electron energies, momenta, and scattering angle.
double t_in, e_in, p_in, p_inx, p_inz;
double e_scat, t_scat, p_scat, p_scatx, p_scatz;
double angle; // Cosine of the electron scattering angle, not an angle itself.
double theta_angle, theta_gamma; // Radians.
double epsilon; // Dimensionless virtual-photon polarization.

// Scan and final-state angles; thpq_lab and thmq_lab are relative to q,
// while thp_lab is a beam-relative, in-plane angle estimate.
double phix, thpq_lab, thp_lab; // Degrees.
double thp_cm, thmq_lab; // Radians.
double a, b, c, cos_thp; // Quadratic coefficients and selected root.
double p_proton_lab, e_proton, p_meson_lab, e_meson;

// Mandelstam invariants and the dimensionless Bjorken scaling variable.
double t, u; // MeV^2.
double t_g, u_g; // GeV^2.
double x;
bool change_beam_e; // Request one header for the next beam's table.
} // namespace

// Calculate p(e,e'p)meson kinematics with a stationary target proton.
// Natural units (c = 1) are used: masses and energies share the same units.
int main() {
    // The proton azimuth about the virtual-photon direction is in degrees.
    // Zero selects the scattering plane.
    phix = 0;

    printf("*--------------------------------------------------*/\n");
    printf("* omega meson mass: 0.78265 GeV \n");
    printf("* phi meson mass: 1.019461 GeV \n");
    printf("* J/psi meson mass: 3.096916 GeV \n");
    printf("*--------------------------------------------------*/\n");

    printf(" Enter meson mass, Q^2, W in GeV\n");

    // Input mass and W are in GeV; Q^2 is in GeV^2, not GeV.
    // This legacy input path still assumes valid, physically allowed values.
    cin >> m_meson_g >> q_2_g >> W_g;
    cout << m_meson_g << "  " << q_2_g << "  " << W_g << endl;

    // Work in MeV for particle masses and the CM two-body calculation.
    // Quantities ending in _g retain GeV-based units.
    m_meson = m_meson_g * 1000.;
    W = W_g * 1000.;

    // W^2 = M_p^2 + 2 M_p nu - Q^2 for a proton at rest.
    // Solving for nu gives the energy transferred from the electron.
    nu_g = (pow(W_g, 2) + q_2_g - pow(m_p_g, 2) ) / (2 * m_p_g);
    nu = nu_g * 1000.;

    // The virtual photon is spacelike: q^2 = nu^2 - |q|^2 = -Q^2.
    // Thus its three-momentum magnitude is sqrt(nu^2 + Q^2).
    p_gamma_g = sqrt(pow(nu_g, 2) + q_2_g);
    p_gamma   = p_gamma_g * 1000.;

    // The photon-proton system has lab energy nu + M_p and momentum |q|.
    // These define the boost from its center-of-mass (CM) frame to the lab.
    beta_cm = p_gamma / (nu + m_p);
    gamma_cm = (nu + m_p) / W;

    // Two-body final states have equal and opposite CM momenta.
    // This is sqrt(lambda(W^2, m_meson^2, M_p^2)) / (2 W), where
    // lambda(x,y,z) = x^2 + y^2 + z^2 - 2xy - 2xz - 2yz.
    // Physical production requires W >= M_p + m_meson; the current
    // angular solver is singular at the exact zero-momentum threshold.
    p_cm = sqrt(pow((pow(W, 2) + pow(m_meson, 2) - pow(m_p, 2)), 2) - 4. * pow((W * m_meson), 2) ) / (2.0 * W);

    // Each final-state particle obeys E^2 = p^2 + m^2 in the CM frame.
    e_meson_cm = sqrt(pow(p_cm, 2) + pow(m_meson, 2));

    e_proton_cm = sqrt(pow(p_cm, 2) + pow(m_p, 2));

    // Determine the number of configured beam energies from the array.
    const std::size_t size = sizeof(beam_energy) / sizeof(beam_energy[0]);

    for (std::size_t i = 0; i < size; ++i) {
        // beam_energy contains kinetic energy, so add the electron rest
        // mass to obtain total energy. The incoming beam defines the +z axis.
        t_in = beam_energy[i];

        e_in = t_in + m_e;

        p_in = sqrt(pow(e_in, 2) - pow(m_e, 2));
        p_inz = p_in;
        p_inx = 0;

        // Energy conservation fixes the scattered electron total energy.
        e_scat = e_in - nu;

        cout << "Electron scattering energy: " << e_scat << endl;

        // Skip beam energies that cannot supply the requested transfer.
        // This skips one beam configuration rather than exiting main().
        if (e_scat < m_e) {
            cerr << "Scattering energy is less than electron mass!" << endl
                << "Program ends!" << endl;
            continue;
        }
        t_scat = e_scat - m_e;

        p_scat = sqrt(pow(e_scat, 2) - pow(m_e, 2));

        // Apply the cosine rule to q = p_in - p_scat. Here angle is
        // cos(theta_e), while theta_angle below is measured in radians.
        // The current implementation does not validate the acos domain.
        angle = ( pow(p_in, 2) + pow(p_scat, 2) - pow(p_gamma, 2) ) / (2. * p_in * p_scat);

        theta_angle = acos(angle);

        // Choose the outgoing electron in the +x half of the x-z plane.
        p_scatz = p_scat * cos(theta_angle);
        p_scatx = p_scat * sin(theta_angle);

        // Virtual-photon polarization in the usual electron-scattering
        // approximation: epsilon = 1 / [1 + 2 |q|^2 / Q^2 tan^2(theta_e / 2)].
        // This expression requires Q^2 > 0.
        epsilon = 1.0 / (1.0 + 2.0 * pow(p_gamma_g, 2) / q_2_g * pow((tan(theta_angle / 2.)), 2));

        // Momentum conservation determines the photon direction; its x
        // component is negative with the electron-plane convention above.
        p_gammaz = p_inz - p_scatz;
        p_gammax = p_inx - p_scatx;

        // Photon angle relative to the beam, in radians. This legacy atan
        // expression assumes a positive longitudinal photon momentum.
        theta_gamma = atan(p_gammax / p_gammaz);

        // Print one table header per beam, then scan proton lab angles
        // relative to the virtual-photon axis (theqp_generic is in degrees).
        change_beam_e = true;

        test_out(theqp_generic, sizeof(theqp_generic) / sizeof(theqp_generic[0]));

        cout << "-----------------------------" << endl;
    }

    return 0;
}

// Solve the proton CM angle for each requested lab angle and print a row.
// Only the +sqrt branch of the quadratic is selected; the alternate branch
// may describe another allowed two-body solution for the same lab angle.
void test_out(const double array[], std::size_t limit) {
    for (std::size_t ii = 0; ii < limit; ++ii) {
        thpq_lab = array[ii];

        // Convert the photon angle before combining it with a degree angle.
        thp_lab = -(abs(Rad_to_Deg(theta_gamma)) + cos(Deg_to_Rad(phix)) * thpq_lab);

        // Under a boost along q, transverse momentum is unchanged and
        // p_parallel_lab = gamma (p_cm cos(theta_cm) + beta E_proton_cm).
        // tan(theta_lab) = p_transverse / p_parallel_lab gives a quadratic
        // a C^2 + b C + c = 0 in C = cos(theta_proton_cm).
        a = pow((tan(Deg_to_Rad(thpq_lab)) * gamma_cm * p_cm), 2) + pow(p_cm, 2);
        b = 2.0 * p_cm * e_proton_cm * beta_cm * pow(tan(Deg_to_Rad(thpq_lab)) * gamma_cm, 2);
        c = pow((tan(Deg_to_Rad(thpq_lab)) * gamma_cm * beta_cm * e_proton_cm), 2) - pow(p_cm, 2);

        // A negative discriminant means no real solution for this angle.
        // Stop the increasing-angle scan at its first inaccessible angle.
        if (pow(b, 2) - 4.0 * a * c < 0) {
            cout << "b^2-4*a*c value is less than 0! " << endl
                << "For loop breaks! " << endl;
            break;
        }

        // Recover the CM polar angle in radians from the selected root.
        cos_thp = (-b + sqrt(pow(b, 2) - 4.0 * a * c)) / (2.0 * a);

        if (abs(cos_thp) <= 1) {
            thp_cm = acos(cos_thp);
        } else {
            cerr << "Cos of the proton angle is larger than 1 ! " << endl
                << "Program end! " << endl;
            break;
        }

        Cal_print();
    }
}

// Boost the final-state particles, evaluate invariants, and format one row.
void Cal_print() {
    // Historical fixed momentum settings, disabled in favor of the boost.
    // if (q_2_g == 1.6) {
    //     p_proton_lab = 2931;
    // } else if (q_2_g == 2.45) {
    //     p_proton_lab = 3336;
    // }

    // The boost gives the longitudinal lab momentum; divide by
    // cos(theta_lab) to obtain its magnitude for the selected branch.
    // Apply this kinematic solution for every Q^2 setting.
    p_proton_lab = gamma_cm * (p_cm * cos(thp_cm) + beta_cm * e_proton_cm) / cos(Deg_to_Rad(thpq_lab));

    e_proton = sqrt(pow(p_proton_lab, 2) + pow(m_p, 2));

    // Proton and meson are back-to-back in the CM frame. The meson has
    // longitudinal momentum -p_cm*cos(thp_cm) and transverse magnitude
    // p_cm*sin(thp_cm). A boost along q changes only the longitudinal
    // component. Both components below are in MeV/c.
    const double meson_transverse = p_cm * sin(thp_cm);
    const double meson_longitudinal =
        gamma_cm * (beta_cm * e_meson_cm - p_cm * cos(thp_cm));

    // atan2 retains the quadrant: negative longitudinal momentum gives
    // a backward angle between pi/2 and pi, including pi on the -q axis.
    // It also handles zero longitudinal momentum without division by zero.
    // The polar angle is in radians here and converted to degrees for output.
    thmq_lab = atan2(meson_transverse, meson_longitudinal);

    // A momentum magnitude is nonnegative regardless of travel direction.
    // hypot avoids dividing by cos(theta), which vanishes at 90 degrees.
    p_meson_lab = hypot(meson_transverse, meson_longitudinal);

    e_meson = sqrt(pow(p_meson_lab, 2) + pow(m_meson, 2));

    // With metric (+,-,-,-), t = (p_target - p_final_proton)^2
    // = 2 M_p (M_p - E_proton_lab). Convert MeV^2 to GeV^2.
    t = 2. * m_p * (m_p-e_proton);

    t_g = t / 1.e6;

    // Similarly, u = (p_target - p_meson)^2. Both t and u are
    // Lorentz invariants; the table reports their negatives, -t and -u.
    u = pow((m_p - e_meson), 2) - pow(p_meson_lab, 2);
    u_g = u / 1.e6;

    // Bjorken x is dimensionless. Use GeV units throughout this ratio.
    x = q_2_g / (2.0 * m_p_g * nu_g);

    // Emit the configuration and units only before the first output row.
    if (change_beam_e) {
        cout << "            p(e,e*p)meson kinematics for Q^2 = " << q_2_g << " GeV^2" << "    W = " << W_g
            << " GeV    " << endl
            << "               nu = omega = " << nu << "MeV    Pgamma = " << p_gamma << " MeV/c     "
            << endl << "                   Meson Mass = " << m_meson << " MeV    " << "Phi = " << phix
            << endl;

        cout << "  Tinc       Te*      THe*    epsilon       -t     "
            << "THq     Pp       THpq     Pm        Thmq    -u         x " << endl
            << "  MeV        MeV      deg                  GeV^2   "
            << "deg    MeV/c     deg     MeV/c       deg   GeV^2  " << endl;

        change_beam_e = false;
    }

    // All table angles are degrees; C++ trigonometric functions use radians.
    // THpq already comes from a degree-valued scan array.
    const double electron_angle_deg = Rad_to_Deg(theta_angle);
    const double photon_angle_deg = Rad_to_Deg(theta_gamma);
    const double proton_angle_deg = thpq_lab;
    const double meson_angle_deg = Rad_to_Deg(thmq_lab);

    // Column order: beam kinetic energy, scattered kinetic energy, electron
    // angle, epsilon, -t, photon angle, proton momentum, proton-q angle,
    // meson momentum, meson-q angle, -u, and Bjorken x.
    // Energy: MeV; momentum: MeV/c; invariant squares: GeV^2; angles: deg.
    cout << setprecision(7) << setw(8) << t_in << " " << setw(9) << t_scat << "   " << setw(6)
        << setprecision(5) << electron_angle_deg << setprecision(4)
        << setw(5) << " " << epsilon << "  " << fixed << setprecision(3) << setw(9) << -t_g << " " << setw(8)
        << photon_angle_deg << "  " << setprecision(3) << p_proton_lab;

    cout << " " << setprecision(1) << setw(4) << proton_angle_deg
        << " " << setprecision(2) << setw(10) << p_meson_lab << " " << setprecision(4)
        << setw(8) << fixed << meson_angle_deg;

    // Restore the default float format before formatting the invariant columns.
    cout.unsetf(ios_base::floatfield);

    cout << " " << fixed << setprecision(4) << setw(7) << -u_g << "    " << setprecision(4) << setw(7) << x << endl;
    // Leave the stream ready for significant-digit formatting on the next row.
    cout.unsetf(ios_base::floatfield);
}
