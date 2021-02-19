/* uvg_common.hpp

   B. Bird - 07/02/2020
*/

#ifndef UVG_COMMON_HPP
#define UVG_COMMON_HPP

#include <vector>
#include <string>

struct MotionVector{
    int x;
    int y;
    int code;
};

//Convenience function to wrap around the nasty notation for 2d vectors
template<typename T>
std::vector<std::vector<T> > create_2d_vector(unsigned int outer, unsigned int inner){
    std::vector<std::vector<T> > V {outer, std::vector<T>(inner,T() )};
    return V;
}

//The floating point calculations we use while converting between 
//RGB and YCbCr can occasionally yield values slightly out of range
//for an unsigned char (e.g. -1 or 255.9).
//Furthermore, we want to ensure that any conversion uses rounding
//and not truncation (to improve accuracy).
inline unsigned char round_and_clamp_to_char(double v){
    //Round to int 
    int i = (int)(v+0.5);
    //Clamp to the range [0,255]
    if (i < 0)
        return 0;
    else if (i > 255)
        return 255;
    return i;
}

/* The exact RGB <-> YCbCr conversion formula used here is the "JPEG style"
   conversion (there is some debate over the best conversion formula) */
struct PixelYCbCr;
struct PixelRGB{
    unsigned char r, g, b;
    PixelYCbCr to_ycbcr(); //Implementation is below (since the PixelYCbCr type has to exist before we can fully define this function)
};

struct PixelYCbCr{
    unsigned char Y, Cb, Cr;
    inline PixelRGB to_rgb(){
        return {
            round_and_clamp_to_char(Y + 1.402*(Cr-128.0)),
            round_and_clamp_to_char(Y-0.344136*(Cb-128.0)-0.714136*(Cr-128.0)),
            round_and_clamp_to_char(Y+1.772*(Cb-128.0))
        };
    }
};


inline PixelYCbCr PixelRGB::to_ycbcr(){
    return {
        round_and_clamp_to_char(0.299*r + 0.587*g + 0.114*b),
        round_and_clamp_to_char(128 + -0.168736*r + -0.331264*g + 0.5*b),
        round_and_clamp_to_char(128 + 0.5*r + -0.418688*g + -0.081312*b)
    };
}


#endif