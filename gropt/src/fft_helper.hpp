#ifndef FFT_HELPER_H
#define FFT_HELPER_H

/**
 * This class is a helper that is mostly to let you switch FFT libraries more easily.
 * PocketFFT is easier to bundle, but is slower, ideally you will have FFTW installed
 * and available to use for this.
 * This also works with MKL, with their drop in FFTW replacements
 */

#include <iostream> 
#include <string>
#include <vector>
#include "Eigen/Dense" 

// Here is where you pick which library you want to use, uncomment only one
// TODO: Move to a build option or similar, so you can pick at build time

// #define FFTMODE_POCKETFFT
#define FFTMODE_POCKETFFT 

#if defined(FFTMODE_POCKETFFT)
    #include "pocketfft_hdronly.hpp"
#elif defined(FFTMODE_FFTW)
    #include "fftw3.h" 
#endif 

namespace Gropt {

// Note the FFT_Helper class is for a single size FFT, in order to use FFTW wisdom and
// pre-allocate some vectors that get re-used.
class FFT_Helper
{
    public:
        unsigned long long N_FT;

        #if defined(FFTMODE_POCKETFFT)
            Eigen::VectorXcd ft_c_vec0;
            Eigen::VectorXcd ft_c_vec1;

            pocketfft::shape_t shape;
            pocketfft::shape_t shape_resize;
            pocketfft::stride_t stride_cd;
            pocketfft::shape_t axes;
        #elif defined(FFTMODE_FFTW)
            Eigen::VectorXcd ft_c_vec0;
            Eigen::VectorXcd ft_c_vec1;

            fftw_plan plan_forward;
            fftw_plan plan_backward;
        #endif

        FFT_Helper(int N);
        virtual ~FFT_Helper();
        virtual void fft(Eigen::VectorXd &X, Eigen::VectorXcd &out);
        virtual void ifft(Eigen::VectorXcd &X, Eigen::VectorXd &out);
        virtual void fft_convolve(Eigen::VectorXd &X, Eigen::VectorXd &out, Eigen::VectorXcd &H, bool do_shift, bool transpose);
};

}  // end namespace Gropt 

#endif