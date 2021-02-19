/* uvid_decompress.cpp
   CSC 485B/578B/SENG 480B - Data Compression - Summer 2020
   Starter code for Assignment 5
   
   This placeholder code reads the (basically uncompressed) data produced by
   the uvid_compress starter code and outputs it in the uncompressed 
   YCbCr (YUV) format used for the sample video input files. To play the 
   the decompressed data stream directly, you can pipe the output of this
   program to the ffplay program, with a command like 
     ffplay -f rawvideo -pixel_format yuv420p -framerate 30 -video_size 352x288 - 2>/dev/null
   (where the resolution is explicitly given as an argument to ffplay).
   B. Bird - 07/15/2020
*/

#include <iostream>
#include <fstream>
#include <array>
#include <string>
#include <cassert>
#include <math.h>
#include <cstdint>
#include <tuple>
#include "input_stream.hpp"
#include "yuv_stream.hpp"
#include "uvg_common.hpp"
int num_frames;
std::vector<std::vector<int>> luminance(8);
std::vector<std::vector<int>> chrominance(8);

// Initializes the quantization matrices
void init(){
    luminance = create_2d_vector<int>(8,8);
    chrominance = create_2d_vector<int>(8,8);

    luminance = {
        {16, 11, 10, 16, 24, 40, 51, 61},
        {12, 12, 14, 19, 26, 58, 60, 55}, 
        {14, 13, 16, 24, 40, 57, 69, 56},
        {14, 17, 22, 29, 51, 87, 80, 62},
        {18, 22, 37, 56, 68, 109, 103, 77},
        {24, 35, 55, 64, 81, 104, 113, 92},
        {49, 64, 78, 87, 103, 121, 120, 101},
        {72, 92, 95, 98, 112, 100, 103, 99}
    };
    
    for (int i = 0; i < 8; i++){
        for (int j = 0; j < 8; j++){
            chrominance.at(i).at(j) = 99;
        }
    }
    chrominance.at(0).at(0) = 17;
    chrominance.at(0).at(1) = 18;
    chrominance.at(0).at(2) = 24;
    chrominance.at(0).at(3) = 47;
    chrominance.at(1).at(0) = 18;
    chrominance.at(1).at(1) = 21;
    chrominance.at(1).at(2) = 26;
    chrominance.at(1).at(3) = 66;
    chrominance.at(2).at(0) = 24;
    chrominance.at(2).at(1) = 26;
    chrominance.at(2).at(2) = 56;
    chrominance.at(3).at(0) = 47;
    chrominance.at(3).at(1) = 66;
}

// Reads a value from stream written in unary
int read_delta(InputBitStream &stream){
    int start = stream.read_bit();
    if (start == 0){
        return 0;
    }else{
        int sign = stream.read_bit() == 0 ? 1 : -1;
        int value = 0;
        while (stream.read_bit() == 1){
            value++;
        }
        return value * sign;
    }
}

// Reads in a p-frame from stream
std::vector<std::vector<int>> get_quantized_pframe(InputBitStream &stream, int width, int height, std::vector<std::vector<int>> &motion_codes){
    auto quantized = create_2d_vector<int>(height, width);
    int y = 0;
    int x = 0;
    while (y <= height - 8){
        while (x <= width - 8){
            int sign = stream.read_bit() == 1 ? 1 : -1;
            int code = stream.read_bits(2);
            quantized.at(y).at(x) = stream.read_u16() * sign;
            int prev = quantized.at(y).at(x);
            motion_codes.at(y).at(x) = code;
            for (int j = 0; j < 8; j++){
                for (int i = 0; i < 8; i++){
                    if (j != 0 || i != 0){
                        int delta = read_delta(stream);
                        quantized.at(y + j).at(x + i) = prev + delta;
                        motion_codes.at(y + j).at(x + i) = code;
                        prev = quantized.at(y + j).at(x + i);
                    }
                }
            }
            x += 8;
        }
        y += 8;
        x = 0;
    }
    return quantized;
}

// Reads in an i-frame from stream
std::vector<std::vector<int>> get_quantized_iframe(InputBitStream &stream, int width, int height){
    auto quantized = create_2d_vector<int>(height, width);
    int y = 0;
    int x = 0;
    while (y <= height - 8){
        while (x <= width - 8){
            quantized.at(y).at(x) = stream.read_u16();
            int prev = quantized.at(y).at(x);
            for (int j = 0; j < 8; j++){
                for (int i = 0; i < 8; i++){
                    if (j != 0 || i != 0){
                        quantized.at(y + j).at(x + i) = prev + read_delta(stream);
                        prev = quantized.at(y + j).at(x + i);
                    }
                }
            }
            x += 8;
        }
        y += 8;
        x = 0;
    }
    return quantized;
}

// Inverts the DCT for the matrix @q
std::vector<std::vector<int>> inv_dct(const std::vector<std::vector<int>> &q){
    std::vector<std::vector<int>> dct_vals = create_2d_vector<int>(8,8);
    float ci, cj, sum;
    for (int y = 0; y < 8; y++){
        for (int x = 0; x < 8; x++){
            sum = 0;
            for (int j = 0; j < 8; j++){
                for (int i = 0; i < 8; i++){
                    ci = i == 0 ? 1. / sqrt(2) : 1.;
                    cj = j == 0 ? 1. / sqrt(2) : 1.;
                    sum += ci * cj * q.at(j).at(i) * cos(M_PI * i * (2.0 * x + 1)/(2*8)) * cos(M_PI * j * (2.0 * y + 1)/(2*8));
                }
            }
            dct_vals.at(y).at(x) = round(1.0/4 * sum);
        }    
    }

    return dct_vals;
}

// Inverts a quantized matrix
void inv_quantized(std::vector<std::vector<int>> &vals, const std::vector<std::vector<int>> &factors, unsigned int width, unsigned int height, float quality){
    unsigned int x = 0;
    unsigned int y = 0;
    while (y <= height - 8){
        while (x <= width - 8){
            for (int j = 0; j < 8; j++){
                for (int i = 0; i < 8; i++){
                    vals.at(y + j).at(x + i) *= (factors.at(j).at(i) / quality);
                }
            }
            x += 8;
        }
        y += 8;
        x = 0;
    }
}

// Parses through one plane of a frame in 8x8 blocks
// Inverts the DCT of each block
void untransform(std::vector<std::vector<int>> &quant, unsigned int padded_height, unsigned int padded_width){
    unsigned int x = 0;
    unsigned int y = 0;
    while (y <= padded_height - 8){
        while (x <= padded_width - 8){
            std::vector<std::vector<int>> vals = create_2d_vector<int>(8,8);
            for (int j = 0; j < 8; j++){
                for (int i = 0; i < 8; i++){
                    vals.at(j).at(i) = quant.at(y + j).at(x + i);
                }
            }
            std::vector<std::vector<int>> inverted = inv_dct(vals);
            for (int j = 0; j < 8; j++){
                for (int i = 0; i < 8; i++){
                    quant.at(y + j).at(x + i) = inverted.at(j).at(i);
                }
            }
            x += 8;
        }
        y += 8;
        x = 0;
    }
}

// Same as reading an i-frame except we have to add the difference afterwards
void read_pframe(InputBitStream &input_stream, YUVFrame420 &prev_frame, YUVStreamWriter &writer, u32 height, u32 width, float quality){
    YUVFrame420& frame = writer.frame();

    unsigned int padded_height = height + (8 - (height % 8));
    unsigned int padded_width = width + (8 - (width % 8));

    unsigned int scaled_height = (height+1)/2;
    unsigned int scaled_width = (width+1)/2;
    unsigned int padded_scaled_height = scaled_height + 8 - (scaled_height % 8);
    unsigned int padded_scaled_width = scaled_width + 8 - (scaled_width % 8);

    auto predicted_blocks = create_2d_vector<bool>(padded_height, padded_width);
    auto motion_codes = create_2d_vector<int>(padded_height, padded_width);

    int x_mv[4] = {0, 8, 0, -8};
    int y_mv[4] = {8, 0, -8, 0};

    // Read Y values first
    auto y_quant = get_quantized_pframe(input_stream, padded_width, padded_height, motion_codes);
    inv_quantized(y_quant, luminance, padded_width, padded_height, quality);

    // Invert the DCT for the Y vals
    untransform(y_quant, padded_height, padded_width);

    // Now add the differences from the previous frame
    for (u32 j = 0; j < height; j++){
        for (u32 i = 0; i < width; i++){
            int x_motion = x_mv[motion_codes.at(j).at(i)];
            int y_motion = y_mv[motion_codes.at(j).at(i)];
            frame.Y(i,j) = round_and_clamp_to_char(y_quant.at(j).at(i) + prev_frame.Y(i + x_motion,j + y_motion));
        }
    }

    // Then get the scaled colour planes
    auto cb_quant = get_quantized_pframe(input_stream, padded_scaled_width, padded_scaled_height, motion_codes);
    inv_quantized(cb_quant, chrominance, padded_scaled_width, padded_scaled_height, quality);

    // Invert the DCT for the cb plane
    untransform(cb_quant, padded_scaled_height, padded_scaled_width);

    // Add the differences from the previous frame
    for (u32 j = 0; j < (height+1)/2; j++){
        for (u32 i = 0; i < (width+1)/2; i++){
            int x_motion = x_mv[motion_codes.at(j).at(i)];
            int y_motion = y_mv[motion_codes.at(j).at(i)];
            frame.Cb(i,j) = round_and_clamp_to_char(cb_quant.at(j).at(i) + prev_frame.Cb(i + x_motion, j + y_motion));
        }
    }

    // Get the values for the Cr plane
    auto cr_quant = get_quantized_pframe(input_stream, padded_scaled_width, padded_scaled_height, motion_codes);
    inv_quantized(cr_quant, chrominance, padded_scaled_width, padded_scaled_height, quality);

    // Invert the DCT for the colour planes
    untransform(cr_quant, padded_scaled_height, padded_scaled_width);

    // Add the differences from the previous frame
    for (u32 j = 0; j < (height+1)/2; j++){
        for (u32 i = 0; i < (width+1)/2; i++){
            int x_motion = x_mv[motion_codes.at(j).at(i)];
            int y_motion = y_mv[motion_codes.at(j).at(i)];
            frame.Cr(i,j) = round_and_clamp_to_char(cr_quant.at(j).at(i) + prev_frame.Cr(i + x_motion,j + y_motion));
        }
    }
    prev_frame = frame;
    writer.write_frame();
}

// Reads an i-frame from stream
void read_iframe(InputBitStream &input_stream, YUVFrame420 &prev_frame, YUVStreamWriter &writer, u32 height, u32 width, float quality){
    YUVFrame420& frame = writer.frame();

    unsigned int padded_height = height + (8 - (height % 8));
    unsigned int padded_width = width + (8 - (width % 8));

    unsigned int scaled_height = (height+1)/2;
    unsigned int scaled_width = (width+1)/2;
    unsigned int padded_scaled_height = scaled_height + 8 - (scaled_height % 8);
    unsigned int padded_scaled_width = scaled_width + 8 - (scaled_width % 8);

    // Y values first
    auto y_quant = get_quantized_iframe(input_stream, padded_width, padded_height);
    inv_quantized(y_quant, luminance, padded_width, padded_height, quality);

    // Then scale colour planes
    auto cb_quant = get_quantized_iframe(input_stream, padded_scaled_width, padded_scaled_height);
    inv_quantized(cb_quant, chrominance, padded_scaled_width, padded_scaled_height, quality);
    auto cr_quant = get_quantized_iframe(input_stream, padded_scaled_width, padded_scaled_height);
    inv_quantized(cr_quant, chrominance, padded_scaled_width, padded_scaled_height, quality);

    // Invert the DCT for the Y vals
    untransform(y_quant, padded_height, padded_width);

    // Invert the DCT for the colour planes
    untransform(cb_quant, padded_scaled_height, padded_scaled_width);
    untransform(cr_quant, padded_scaled_height, padded_scaled_width);

    for (u32 j = 0; j < height; j++)
        for (u32 i = 0; i < width; i++)
            frame.Y(i,j) = round_and_clamp_to_char(y_quant.at(j).at(i));
    for (u32 j = 0; j < (height+1)/2; j++)
        for (u32 i = 0; i < (width+1)/2; i++)
            frame.Cb(i,j) = round_and_clamp_to_char(cb_quant.at(j).at(i));
    for (u32 j = 0; j < (height+1)/2; j++)
        for (u32 i = 0; i < (width+1)/2; i++)
            frame.Cr(i,j) = round_and_clamp_to_char(cr_quant.at(j).at(i));
    prev_frame = frame;
    writer.write_frame();
}

int main(int argc, char** argv){
    init();
    //Note: This program must not take any command line arguments. (Anything
    //      it needs to know about the data must be encoded into the bitstream)
    
    InputBitStream input_stream {std::cin};

    u32 height {input_stream.read_u32()};
    u32 width {input_stream.read_u32()};
    unsigned int quality_bits = input_stream.read_bits(2);
    float quality = quality_bits == 0 ? 0.2 : (quality_bits == 1 ? 1 : 4);
    unsigned int quality_freq = quality_bits == 0 ? 12 : (quality_bits == 1 ? 10 : 5);

    YUVStreamWriter writer {std::cout, width, height};
    num_frames = 0;
    YUVFrame420 prev_frame(height, width);
    if (input_stream.read_byte()){
        read_iframe(input_stream, prev_frame, writer, height, width, quality);
        num_frames++;
    }

    while (input_stream.read_byte()){
        num_frames++;
        if (num_frames % quality_freq != 0){
            read_pframe(input_stream, prev_frame, writer, height, width, quality);
        }else{
            read_iframe(input_stream, prev_frame, writer, height, width, quality);
        }
    }

    return 0;
}