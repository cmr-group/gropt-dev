#include <iostream> 
#include <string>
#include <math.h>  
#include "Eigen/Dense"

#include "fft_helper.hpp" 
// #include "mkl.h"

namespace Gropt {

FFT_Helper::FFT_Helper(int N) 
{
    N_FT = N;

    #if defined(FFTMODE_POCKETFFT)
        // FFT params
        shape.push_back(N_FT);
        stride_cd.push_back(sizeof(std::complex<double>));
        axes.push_back(0);

        ft_c_vec0.setZero(N_FT);
        ft_c_vec1.setZero(N_FT);
    #elif defined(FFTMODE_FFTW)
        // Threading doesn't really help for our small-ish 1D FFTs used here

        // omp_set_num_threads(24);

        // fftw_init_threads();
        // fftw_plan_with_nthreads(2); 

        ft_c_vec0.setZero(N_FT);
        ft_c_vec1.setZero(N_FT);

        plan_forward = fftw_plan_dft_1d(N_FT, reinterpret_cast<fftw_complex*>(ft_c_vec0.data()), reinterpret_cast<fftw_complex*>(ft_c_vec1.data()), FFTW_FORWARD, FFTW_ESTIMATE);
        plan_backward = fftw_plan_dft_1d(N_FT, reinterpret_cast<fftw_complex*>(ft_c_vec0.data()), reinterpret_cast<fftw_complex*>(ft_c_vec1.data()), FFTW_BACKWARD, FFTW_ESTIMATE);
    #endif

}

FFT_Helper::~FFT_Helper() {
    #if defined(FFTMODE_FFTW)
        fftw_destroy_plan(plan_forward);
        fftw_destroy_plan(plan_backward);
    #endif
}


void FFT_Helper::fft(Eigen::VectorXd &X, Eigen::VectorXcd &out) 
{
    #if defined(FFTMODE_POCKETFFT)
        ft_c_vec0.setZero();
        ft_c_vec1.setZero();

        ft_c_vec0.real() = X;

        pocketfft::c2c(shape, stride_cd, stride_cd, axes, pocketfft::FORWARD,
                        ft_c_vec0.data(), ft_c_vec1.data(), 1.);

        out = ft_c_vec1;
    #elif defined(FFTMODE_FFTW)
        ft_c_vec0.setZero();
        ft_c_vec0.real() = X/sqrtf(N_FT);
        fftw_execute_dft(plan_forward, reinterpret_cast<fftw_complex*>(ft_c_vec0.data()), reinterpret_cast<fftw_complex*>(out.data()));
    #endif

}




void FFT_Helper::ifft(Eigen::VectorXcd &X, Eigen::VectorXd &out) 
{
    #if defined(FFTMODE_POCKETFFT)
        ft_c_vec0.setZero();
        ft_c_vec1.setZero();

        ft_c_vec1 = X;

        pocketfft::c2c(shape, stride_cd, stride_cd, axes, pocketfft::BACKWARD,
                        ft_c_vec1.data(), ft_c_vec0.data(), 1.0/N_FT);

        out = ft_c_vec0.real();
    #elif defined(FFTMODE_FFTW)
        fftw_execute_dft(plan_backward, reinterpret_cast<fftw_complex*>(X.data()), reinterpret_cast<fftw_complex*>(ft_c_vec1.data()));
        out = ft_c_vec1.real()/sqrtf(N_FT);
    #endif

}


void FFT_Helper::fft_convolve(Eigen::VectorXd &X, Eigen::VectorXd &out, Eigen::VectorXcd &H, bool do_shift, bool transpose) 
{
    #if defined(FFTMODE_POCKETFFT)
        ft_c_vec0.setZero();
        ft_c_vec1.setZero();

        ft_c_vec0.real() = X;

        pocketfft::c2c(shape, stride_cd, stride_cd, axes, pocketfft::FORWARD,
                        ft_c_vec0.data(), ft_c_vec1.data(), 1.);

        if (!transpose) {
            if (!do_shift) {
                ft_c_vec1.array() *= H.array();
            } else {
                int center = (int)ft_c_vec1.size()/2;
                ft_c_vec1.segment(center,center).array() *= H.segment(0, center).array();
                ft_c_vec1.segment(0,center).array() *= H.segment(center, center).array();
            }
        } else {
            if (!do_shift) {
                ft_c_vec1.array() *= H.conjugate().array();
            } else {
                int center =  (int)ft_c_vec1.size()/2;
                ft_c_vec1.segment(center,center).array() *= H.conjugate().segment(0, center).array();
                ft_c_vec1.segment(0,center).array() *= H.conjugate().segment(center, center).array();
            }
        }

        pocketfft::c2c(shape, stride_cd, stride_cd, axes, pocketfft::BACKWARD,
                        ft_c_vec1.data(), ft_c_vec0.data(), 1.0/N_FT);

        out = ft_c_vec0.real();
    #elif defined(FFTMODE_FFTW)

        ft_c_vec0.setZero();
        ft_c_vec1.setZero();
        ft_c_vec0.real() = X;

        fftw_execute_dft(plan_forward, reinterpret_cast<fftw_complex*>(ft_c_vec0.data()), reinterpret_cast<fftw_complex*>(ft_c_vec1.data()));

        if (!transpose) {
            if (!do_shift) {
                ft_c_vec1.array() *= H.array();
            } else {
                int center = ft_c_vec1.size()/2;
                ft_c_vec1.segment(center,center).array() *= H.segment(0, center).array();
                ft_c_vec1.segment(0,center).array() *= H.segment(center, center).array();
            }
        } else {
            if (!do_shift) {
                ft_c_vec1.array() *= H.conjugate().array();
            } else {
                int center = ft_c_vec1.size()/2;
                ft_c_vec1.segment(center,center).array() *= H.conjugate().segment(0, center).array();
                ft_c_vec1.segment(0,center).array() *= H.conjugate().segment(center, center).array();
            }
        }

        fftw_execute_dft(plan_backward, reinterpret_cast<fftw_complex*>(ft_c_vec1.data()), reinterpret_cast<fftw_complex*>(ft_c_vec0.data()));

        
        out = ft_c_vec0.real()/N_FT;

    #endif

}

}




