#include<arm_neon.h>
#include<stdio.h>

void fir_reset(FIR_T *fir, UINT32 n_taps, const COMPLEX_R64 coeff[])
{
    UINT32 i = 0;
    fir->fir_n_taps = n_taps;
    fir->fir_com_idx = 0;
    for (i = 0; i < n_taps; i +=2)
    {
        // load 4 complex numbers from coeff array
        float64x2_t coeff_real = vld1q_f64(&coeff[i].real);
        float64x2_t coeff_imag = vld1q_f64(&coeff[i].imag);
        // store 4 complex numbers in fir_h array
        vst1q_f64(&fir->fir_h[i].real, coeff_real);
        vst1q_f64(&fir->fir_h[i].imag, coeff_imag);
        // set 4 complex numbers in fir_state array to 0
        vst1q_f64(&fir->fir_state[i].real, vdupq_n_f64(0));
        vst1q_f64(&fir->fir_state[i].imag, vdupq_n_f64(0));
        vst1q_f64(&fir->fir_state[n_taps + i].real, vdupq_n_f64(0));
        vst1q_f64(&fir->fir_state[n_taps + i].imag, vdupq_n_f64(0));
    }
}

void fir_filter(FIR_T *fir, COMPLEX_R64 *x, UINT32 n, COMPLEX_R64 *y)
{
    UINT32 com_idx = fir->fir_com_idx;
    UINT32 n_taps = fir->fir_n_taps;
    UINT32 i, j,k;
    float64_t m11[2],m21,m13[2],m12;
    float64x2_t m31,m41;    
    float64x2_t x_real, x_imag, h_real, h_imag;
    float64x2x2_t h_val,x_val;
    double zero[] = {0.0};
    x_val = vld2q_dup_f64(zero);  h_val = vld2q_dup_f64(zero);
    COMPLEX_R64 *h = fir->fir_h,x_v[n_taps];
            
    for (j = 0; j < n; j++)
    {
        fir->fir_state[com_idx] = x[j];
        fir->fir_state[com_idx + n_taps] = x[j];
        float64x2_t vaccum_real = vdupq_n_f64(0);
        float64x2_t vaccum_imag = vdupq_n_f64(0);
        
        for (i = 0; i < n_taps; i++)
        {
            k = (n_taps + com_idx - i);        
            x_v[i] = fir->fir_state[k];                 
        }
        // vectorize the inner loop
        for (i = 0; i < n_taps; i+=2)
        {   
            //interleaved load of x coefficients
            x_val = vld2q_f64(x_v + i);    
            x_real = x_val.val[0];          
            x_imag = x_val.val[1];          
            // interleaved load of filter coefficients
            h_val = vld2q_f64(h + i);
            h_real = h_val.val[0];         
            h_imag = h_val.val[1];         
            //Multiply (a+ib)(x+iy)=(ax-by)+i(ay+bx)
            vaccum_real = vmlaq_f64(vaccum_real, h_real, x_real);     
            vaccum_real = vmlsq_f64(vaccum_real, h_imag, x_imag);  
            vaccum_imag = vmlaq_f64(vaccum_imag, h_real, x_imag);   
            vaccum_imag = vmlaq_f64(vaccum_imag, h_imag, x_real);
            vst1q_f64(m11,vaccum_real); m31=vld1q_f64(m11); m21=vaddvq_f64(m31);   
            vst1q_f64(m13,vaccum_imag); m41=vld1q_f64(m13); m12=vaddvq_f64(m41);    
        }
            com_idx = com_idx + 1;
            if (com_idx == n_taps)
            {
                com_idx = 0;
            }
            y[j].real=(double)m21;
            y[j].imag=(double)m12;
                       
   }
        
        fir->fir_com_idx = com_idx;
           
}
