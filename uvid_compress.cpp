/* uvid_compress.cpp
   CSC 485B/578B/SENG 480B - Data Compression - Summer 2020
   Starter code for Assignment 5
   Reads video data from stdin in uncompresed YCbCr (YUV) format 
   (With 4:2:0 subsampling). To produce this format from 
   arbitrary video data in a popular format, use the ffmpeg
   tool and a command like 
     ffmpeg -i videofile.mp4 -f rawvideo -pixel_format yuv420p - 2>/dev/null | ./this_program <width> height>
   Note that since the width/height of each frame is not encoded into the raw
   video stream, those values must be provided to the program as arguments.
   B. Bird - 07/15/2020
*/

#include <iostream>
#include <fstream>
#include <array>
#include <math.h>
#include <string>
#include <cassert>
#include <cstdint>
#include <tuple>
#include <cmath>
#include "output_stream.hpp"
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

// Pads the 2D array so that the width and height are divisible by 8
void pad(std::vector<std::vector<int>> &frame, u32 &height, u32 &width){
    unsigned int right_padding = 8 - (width % 8);
    unsigned int lower_padding = 8 - (height % 8);
    for (unsigned int y = 0; y < height; y++){
        for (unsigned int x = 0; x < right_padding; x++){
            frame.at(y).at(x + width) = frame.at(y).at(width - 1);
        }
    }

    for (unsigned int y = 0; y < lower_padding; y++){
        for (unsigned int x = 0; x < width; x++){
            frame.at(y + height).at(x) = frame.at(height - 1).at(x);
        }
    }
}

// Applies the 2-dimensional discrete cosine transform to the "matrix" @g
std::vector<std::vector<float>> dct(const std::vector<std::vector<int>> &g){
    std::vector<std::vector<float>> dct_vals = create_2d_vector<float>(8,8);
    float cx, cy, sum;
    for (int y = 0; y < 8; y++){
        for (unsigned int x = 0; x < 8; x++){
            cx = x == 0 ? 1. / sqrt(2) : 1.;
            cy = y == 0 ? 1. / sqrt(2) : 1.;
            sum = 0;
            for (int j = 0; j < 8; j++){
                for (int i = 0; i < 8; i++){
                    sum += cx * cy * g.at(j).at(i) * cos(M_PI * x * (2.0 * i + 1)/(2*8)) * cos(M_PI * y * (2.0 * j + 1)/(2*8));
                }
            }
            dct_vals.at(y).at(x) = 1.0/4 * sum;
        }    
    }

    return dct_vals;
}

// Quantizes the matrix @dct_vals by dividing elementwise by each value in @quant_vals
std::vector<std::vector<int>> quantize(const std::vector<std::vector<float>> &dct_vals, const std::vector<std::vector<int>> &quant_vals, float quality){
    std::vector<std::vector<int>> q = create_2d_vector<int>(8,8);
    for (int j = 0; j < 8; j++){
        for (int i = 0; i < 8; i++){
            q.at(j).at(i) = std::round(dct_vals.at(j).at(i) / (quant_vals.at(j).at(i) / quality));
        }
    }
    return q;
}

// 10 -> positive
// 11 -> negative
// Followed by @difference 1's and a 0 to indicate a stop
void delta(OutputBitStream &stream, int difference){
    if (difference != 0){
        stream.push_bits(difference > 0 ? 1 : 3, 2);
        for (int i = 0; i < abs(difference); i++){
            stream.push_bit(1);
        }
    }
    stream.push_bit(0);
}

// Uses delta compression to output values to stream
void write_codes(OutputBitStream &output_stream, const std::vector<std::vector<int>> &q){
    int prev_val = q.at(0).at(0);
    output_stream.push_u16(abs(prev_val));
    for (int j = 0; j < 8; j++){
        for (int i = 0; i < 8; i++){
            if (j != 0 || i != 0){
                int difference = q.at(j).at(i) - prev_val;
                delta(output_stream, difference);
                prev_val = q.at(j).at(i);
            }
        }
    }
}

// Applies the DCT and outputs an i-frame
void transform_and_output(OutputBitStream &output_stream, std::vector<std::vector<int>> &values, unsigned int padded_height, unsigned int padded_width, float quality_val, bool use_luminance){
    unsigned int x = 0;
    unsigned int y = 0;
    auto eight_grid = create_2d_vector<int>(8,8);
    while (y <= padded_height - 8){
        while (x <= padded_width - 8){
            for (unsigned int j = y; j < y + 8; j++){
                for (unsigned int i = x; i < x + 8; i++){
                    eight_grid.at(j % 8).at(i % 8) = values.at(j).at(i);
                }
            }
            std::vector<std::vector<float>> dct_vals =  dct(eight_grid);
            std::vector<std::vector<int>> quant = quantize(dct_vals, (use_luminance ? luminance : chrominance), quality_val);

            write_codes(output_stream, quant);
            x += 8;
        }
        y += 8;
        x = 0;
    }
}

// Performs local search to get motion vector which minimizes the average absolute difference
// cur_y is the position of the topmost pixel in the current block
// cur_x is the position of the leftmost pixel in the current block
MotionVector get_motion_vector(std::vector<std::vector<int>> &cur_frame, std::vector<std::vector<int>> &prev_frame, int cur_y, int cur_x, int height, int width){
    int upper = INT32_MAX;
    int lower = INT32_MAX;
    int right = INT32_MAX;
    int left = INT32_MAX;
    MotionVector m;
    if (cur_y + 16 < height){
        // Get AAD of the block 8 pixels below
        lower = 0;
        for (int j = 0; j < 8; j++){
            for (int i = 0; i < 8; i++){
                lower += abs(cur_frame.at(cur_y + j).at(cur_x + i) - prev_frame.at(cur_y + j + 8).at(cur_x + i));
            }
        }
    }
    if (cur_y - 8 >= 0){
        // Get AAD of the block 8 pixels above
        upper = 0;
        for (int j = 0; j < 8; j++){
            for (int i = 0; i < 8; i++){
                upper += abs(cur_frame.at(cur_y + j).at(cur_x + i) - prev_frame.at(cur_y + j - 8).at(cur_x + i));
            }
        }
    }
    if (cur_x + 16 < width){
        // Get AAD of the block 8 pixels to the right
        right = 0;
        for (int j = 0; j < 8; j++){
            for (int i = 0; i < 8; i++){
                right += abs(cur_frame.at(cur_y + j).at(cur_x + i) - prev_frame.at(cur_y + j).at(cur_x + i + 8));
            }
        }
    }
    if (cur_x - 8 >= 0){
        // Get AAD of the block 8 pixels to the left
        left = 0;
        for (int j = 0; j < 8; j++){
            for (int i = 0; i < 8; i++){
                left += abs(cur_frame.at(cur_y + j).at(cur_x + i) - prev_frame.at(cur_y + j).at(cur_x + i - 8));
            }
        }
    }
    
    int ultimate_min = std::min(std::min(upper, lower), std::min(right, left));
    if (ultimate_min == lower){
        m.code = 0;
        m.x = 0;
        m.y = 8;
    }else if (ultimate_min == upper){
        m.code = 2;
        m.x = 0;
        m.y = -8;
    }else if (ultimate_min == right){
        m.code = 1;
        m.x = 8;
        m.y = 0;
    }else{
        m.code = 3;
        m.x = -8;
        m.y = 0;
    }
    return m;
}

// For each block, sum up how many bits are needed to output it as a predicted block
// If the number of bits for a predicted block > number of bits for the original data, use a predicted block
void transform_and_output_pframe(OutputBitStream &output_stream, std::vector<std::vector<int>> &cur_frame, std::vector<std::vector<int>> &prev_frame, unsigned int padded_height, unsigned int padded_width, float quality_val, bool use_luminance){
    unsigned int x = 0;
    unsigned int y = 0;
    auto eight_grid = create_2d_vector<int>(8,8);

    while (y <= padded_height - 8){
        while (x <= padded_width - 8){
            // Sum up the values in each matrix
            MotionVector m = get_motion_vector(cur_frame, prev_frame, y, x, padded_height, padded_width);
            for (unsigned int j = y; j < y + 8; j++){
                for (unsigned int i = x; i < x + 8; i++){
                    eight_grid.at(j % 8).at(i % 8) = cur_frame.at(j).at(i) - prev_frame.at(j + m.y).at(i + m.x);
                }
            }
            std::vector<std::vector<float>> dct_vals =  dct(eight_grid);
            std::vector<std::vector<int>> quant = quantize(dct_vals, (use_luminance ? luminance : chrominance), quality_val);
            // Push a 1 if the block is a predicted block, 0 otherwise
            output_stream.push_bit(quant.at(0).at(0) > 0);
            output_stream.push_bits(m.code, 2);
            write_codes(output_stream, quant);
            x += 8;
        }
        y += 8;
        x = 0;
    }
}

// Parses the frame in 8x8 chunks
// Applies the DCT, quantizes, then writes to stream
void output_iframe(OutputBitStream &output_stream, YUVFrame420& frame, u32 &height, u32 &width, float quality_val){
    // Get the frame data into 3 2D vectors each with padded width and height
    // Start with Y values
    unsigned int padded_height = height + (8 - (height % 8));
    unsigned int padded_width = width + (8 - (width % 8));
    auto y_values = create_2d_vector<int>(padded_height, padded_width);
    for (u32 y = 0; y < height; y++)
        for (u32 x = 0; x < width; x++)
            y_values.at(y).at(x) = frame.Y(x,y);
    pad(y_values, height, width);

    transform_and_output(output_stream, y_values, padded_height, padded_width, quality_val, true);

    // Pad the scaled values
    unsigned int scaled_height = (height + 1)/2;
    unsigned int scaled_width = (width + 1)/2;
    unsigned int right_padding = 8 - (scaled_width % 8);
    unsigned int lower_padding = 8 - (scaled_height % 8);
    unsigned int padded_scaled_height = scaled_height + lower_padding;
    unsigned int padded_scaled_width = scaled_width + right_padding;
    auto cb_values = create_2d_vector<int>(padded_scaled_height, padded_scaled_width);
    auto cr_values = create_2d_vector<int>(padded_scaled_height, padded_scaled_width);
    for (u32 y = 0; y < scaled_height; y++){
        for (u32 x = 0; x < scaled_width; x++){
            cb_values.at(y).at(x) = frame.Cb(x,y);
            cr_values.at(y).at(x) = frame.Cr(x,y);
        }
    }
    pad(cb_values, scaled_height, scaled_width);
    pad(cr_values, scaled_height, scaled_width);

    transform_and_output(output_stream, cb_values, padded_scaled_height, padded_scaled_width, quality_val, false);
    transform_and_output(output_stream, cr_values, padded_scaled_height, padded_scaled_width, quality_val, false);
}

// Writes a p-frame
// Differences are computed before the DCT
void make_pframe(OutputBitStream &output_stream, YUVFrame420& frame, YUVFrame420 &prev_frame, u32 &height, u32 &width, float quality_val){
    // Y values first
    unsigned int padded_height = height + (8 - (height % 8));
    unsigned int padded_width = width + (8 - (width % 8));
    auto cur_y = create_2d_vector<int>(padded_height, padded_width);
    auto prev_y = create_2d_vector<int>(padded_height, padded_width);
    for (int y = 0; y < height; y++){
        for (int x = 0; x < width; x++){
            cur_y.at(y).at(x) = (int)frame.Y(x,y); 
            prev_y.at(y).at(x) = (int)prev_frame.Y(x,y);
        }
    }
    pad(cur_y, height, width);
    pad(prev_y, height, width);

    transform_and_output_pframe(output_stream, cur_y, prev_y, padded_height, padded_width, quality_val, true);

    // Pad the scaled values
    unsigned int scaled_height = (height + 1)/2;
    unsigned int scaled_width = (width + 1)/2;
    unsigned int right_padding = 8 - (scaled_width % 8);
    unsigned int lower_padding = 8 - (scaled_height % 8);
    unsigned int padded_scaled_height = scaled_height + lower_padding;
    unsigned int padded_scaled_width = scaled_width + right_padding;
    auto cur_cb = create_2d_vector<int>(padded_scaled_height, padded_scaled_width);
    auto cur_cr = create_2d_vector<int>(padded_scaled_height, padded_scaled_width);
    auto prev_cb = create_2d_vector<int>(padded_scaled_height, padded_scaled_width);
    auto prev_cr = create_2d_vector<int>(padded_scaled_height, padded_scaled_width);
    for (u32 y = 0; y < scaled_height; y++){
        for (u32 x = 0; x < scaled_width; x++){
            cur_cb.at(y).at(x) = (int)frame.Cb(x,y);
            prev_cb.at(y).at(x) = (int)prev_frame.Cb(x,y);
            cur_cr.at(y).at(x) = (int)frame.Cr(x,y);
            prev_cr.at(y).at(x) = (int)prev_frame.Cr(x,y);
        }
    }
    pad(cur_cb, scaled_height, scaled_width);
    pad(cur_cr, scaled_height, scaled_width);
    pad(prev_cb, scaled_height, scaled_width);
    pad(prev_cr, scaled_height, scaled_width);

    transform_and_output_pframe(output_stream, cur_cb, prev_cb, padded_scaled_height, padded_scaled_width, quality_val, false);
    transform_and_output_pframe(output_stream, cur_cr, prev_cr, padded_scaled_height, padded_scaled_width, quality_val, false);
}

int main(int argc, char** argv){
    init();
    if (argc < 4){
        std::cerr << "Usage: " << argv[0] << " <width> <height> <low/medium/high>" << std::endl;
        return 1;
    }
    u32 width = std::stoi(argv[1]);
    u32 height = std::stoi(argv[2]);
    std::string quality{argv[3]};
    float quality_val = quality == "low" ? 0.2 : (quality == "medium" ? 1 : 4);
    unsigned int quality_freq = quality == "low" ? 12 : (quality == "medium" ? 10 : 5);

    YUVStreamReader reader {std::cin, width, height};
    OutputBitStream output_stream {std::cout};

    output_stream.push_u32(height);
    output_stream.push_u32(width);
    output_stream.push_bits(quality == "low" ? 0 : (quality == "medium" ? 1 : 2), 2);
    num_frames = 0;
    YUVFrame420 prev_frame(width, height);
    // Read first frame and output as an i-frame
    // Then loop to have option of p-frame
    if (reader.read_next_frame()){
        output_stream.push_byte(1);
        YUVFrame420& frame = reader.frame();
        num_frames++;
        output_iframe(output_stream, frame, height, width, quality_val);
        prev_frame = frame;
    }

    while (reader.read_next_frame()){
        output_stream.push_byte(1); //Use a one byte flag to indicate whether there is a frame here
        YUVFrame420& frame = reader.frame();
        num_frames++;
        if (num_frames % quality_freq != 0){
            make_pframe(output_stream, frame, prev_frame, height, width, quality_val);
        }
        else{
            output_iframe(output_stream, frame, height, width, quality_val);
        }
        prev_frame = frame;
    }

    output_stream.push_byte(0); //Flag to indicate end of data
    output_stream.flush_to_byte();

    return 0;
}