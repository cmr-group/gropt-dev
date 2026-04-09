#ifndef OP_ACOUSTIC_H
#define OP_ACOUSTIC_H

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include "Eigen/Dense"

#include "op_main.hpp"
#include "fft_helper.hpp"

namespace Gropt {

class Op_Acoustic : public Operator
{
    protected:
        std::vector<double> freqs;
        std::vector<double> bws;
        double bw_scale = 1.0;

        Eigen::VectorXcd H;
        int N_pad = 32768;
        std::unique_ptr<FFT_Helper> ffth;

    public:
        Op_Acoustic(const ProblemData &_pdata, std::vector<double> _freqs, std::vector<double> _bws, double _weight_mod,
                    double _bw_scale);
        Op_Acoustic(const Op_Acoustic& other);
        virtual void init();

        virtual void forward(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void transpose(Eigen::VectorXd &X, Eigen::VectorXd &out);
        virtual void prox(Eigen::VectorXd &X);
        virtual void check(Eigen::VectorXd &X);
};

}

#endif
