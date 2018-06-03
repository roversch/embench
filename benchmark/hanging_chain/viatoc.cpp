
#include "mpcexport.h"

int main()
{
    constexpr int num_masses = 3;
    constexpr int num_free_masses = num_masses - 2;
    constexpr int nx = (2*num_free_masses + 1)*3;
    constexpr int nu = 3;
    constexpr int N = 40;
    constexpr double Ts = 0.2;

    constexpr double L = 0.033;
    constexpr double D = 1.0;
    constexpr double m = 0.03;
    constexpr double g = 9.81;

    DifferentialState px;
    DifferentialState py;
    DifferentialState pz;
    DifferentialState vx;
    DifferentialState vy;
    DifferentialState vz;
    DifferentialState p_endx;
    DifferentialState p_endy;
    DifferentialState p_endz;
    Control ux;
    Control uy;
    Control uz;
    
    IntermediateState scale1, scale2, F1x, F1y, F1z, F2x, F2y, F2z;
    
    scale1 = -D / m * (1 - L / sqrt(1e-4+px*px+py*py+pz*pz));
    scale2 = D / m * (1 - L / sqrt(1e-4+(p_endx-px)*(p_endx-px)+(p_endy-py)*(p_endy-py)+(p_endz-pz)*(p_endz-pz)));
    
    F1x = scale1 * px;
    F1y = scale1 * py;
    F1z = scale1 * pz;

    F2x = scale2 * (p_endx - px);
    F2y = scale2 * (p_endy - py);
    F2z = scale2 * (p_endz - pz);

    DifferentialEquation f;
    f << dot(px) == vx;
    f << dot(py) == vy;
    f << dot(pz) == vz;
    f << dot(vx) == F1x + F2x;
    f << dot(vy) == F1y + F2y;
    f << dot(vz) == F1z + F2z - g;
    f << dot(p_endx) == ux;
    f << dot(p_endy) == uy;
    f << dot(p_endz) == uz;

    Matrix Q = eye(nx);
    Q(0, 0) = 25;
    Q(1, 1) = 25;
    Q(2, 2) = 25;
    Q(3, 3) = 1;
    Q(4, 4) = 1;
    Q(5, 5) = 1;
    Q(6, 6) = 25;
    Q(7, 7) = 25;
    Q(8, 8) = 25;

    Matrix R = eye(nu)*0.01;

    MPCexport mpc(0.0, Ts*N, N);
    mpc.minimizeLSQ(Q, R);
    mpc.minimizeLSQEndTerm(Q);

    mpc.subjectTo(f);
    mpc.subjectTo(-1.0 <= ux <= 1.0);
    mpc.subjectTo(-1.0 <= uy <= 1.0);
    mpc.subjectTo(-1.0 <= uz <= 1.0);
    
    mpc.set(INTEGRATOR_TYPE, INT_RK4);
    mpc.set(EXPORT_SFUNCTION, YES);
    mpc.set(EXPORT_CPP, YES);

    mpc.exportCode("_viatoc");
    
    return 0;
}
