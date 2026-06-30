#ifndef OP_DIFFBASIN_H
#define OP_DIFFBASIN_H

/**
 * Diffusion basin-orientation constraint (single 180 in the middle).
 *
 *     mean(g, pre-180 window)  >=  +eps
 *     mean(g, post-180 window) <=  -eps        eps = eps_factor * gmax
 *
 * Essentially forces a gradient flip around the 180 (or that they stay the same sign if same_sign=true). 
 * The mean is taken over a window of length window_time on each side of the 180.
 */

#include "Eigen/Dense"
#include <string>

#include "op_main.hpp"

namespace Gropt {

class Op_DiffBasin : public Operator {
  protected:
    double window_time; // window on each side of the 180 [s]
    double eps_factor;  // min |mean window gradient| as a fraction of gmax
    double gmax;        // gradient limit [T/m], scales eps
    double eps = 0.0;   // = eps_factor * gmax (set in init)
    bool same_sign;     // false: post<=-eps (opposite/flip); true: post>=+eps (same sign / no-flip)

    Eigen::VectorXd a_pre;  // mean-gradient functional over the pre-180 window (length Ntot)
    Eigen::VectorXd a_post; // mean-gradient functional over the post-180 window (length Ntot)

  public:
    Op_DiffBasin(const ProblemData &_pdata, double _window_time, double _eps_factor, double _gmax,
                 double _weight_mod, bool _same_sign = false);
    virtual void init();

    virtual void forward(Eigen::VectorXd &X, Eigen::VectorXd &out);
    virtual void transpose(Eigen::VectorXd &X, Eigen::VectorXd &out);
    virtual void prox(Eigen::VectorXd &X);
    virtual void check(Eigen::VectorXd &X);
};

} // namespace Gropt

#endif
