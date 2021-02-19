# Video Compression Program
Author: Ashley Van Spankeren

Description
-----------
This program uses the discrete cosine transform, quantization, predicted frames, and motion vectors to achieve compression.

In the compressor:
* Height, width, and quality arguments are read from command line and immediately written to stream.
* The quality specified determines the frequency of I-frames and P-frames (lower quality allows for more p-frames which leads to better compression).
* Frames are first padded to the height and width at the next multiple of 8.
* Frames are then parsed in 8x8 blocks, applying the following transformations to each block sequentially as follows: 
    * For P-frames, the next step is to determine the motion vector which minimizes the absolute average distance (see Architecture section for how motion vectors are determined). The compressor then writes the 2-bit "motion code" to stream to specify which motion vector was chosen. I-frames skip this step.
    * The discrete cosine transform is then applied to the block.
    * The block is then quantized according to some scalar multiple of the JPEG recommendation. The scalar multiple is determined from the quality argument.
    * Delta compression is applied to write the block to stream.

In the decompressor:
* Height, width, and quality are read from stream.
* Based on the quality that was read from the stream, the decompressor knows the frequency of I-frames compared to P-frames, and hence alternates decompressing each type according to this frequency.
* For P-frames, the motion code is read in for each block. This is stored in the 2D vector motion_codes.
* The delta values are read in and decoded to retrieve the frames quantized values.
* Next, the planes are parsed in 8x8 blocks. For each block, quantization is inverted and the inverse DCT is applied.
* For P-frames, the planes are again parsed in 8x8 blocks to apply motion compensation based on the motion vectors retrieved from the motion codes.
* The frame is reconstructed using the 3 planes and then output to stream.


Architecture
--------------
Algorithms:
* Computing and inverting the DCT:
    * This algorithm is taken directly from the slides on the Video Compression lecture
    * In the compressor, the function dct() takes an 8x8 matrix @g and returns an 8x8 matrix containing the DCT of g.
    * In the decompressor, the function inv_dct() takes an 8x8 matrix @q and returns an 8x8 matrix containing the inverse DCT of q.
* Motion vector computation:
    * Local search is used for finding motion vectors. For each block, the motion vector corresponding to the blocks directly above, below, to the right, and to the left are analyzed. For each potential choice, the absolute average distance is taken (formula taken directly from the Video Compression slides). The block producing the minimum AAD is selected. Given that there are 4 potential motion vectors, the choice can be conveyed to the decompressor in 2 bits:
        * Code 0 represents <0, 8> 
        * Code 1 represents <8, 0>
        * Code 2 represents <0, -8>
        * Code 3 represents <-8, 0>

Data structures:
* MotionVector:
    * Used to represent the motion vector producing the minimum AAD
    * Contains the vertical offset (y), horizontal offset (x), and the motion code described above


Bitstream
-----------
The first 64 bits convey the height and width of the frame (each using 32 bits).
Next, 2 bits are used to convey the quality of the stream:
* 00 --> low
* 01 --> medium
* 10 --> high

Each frame is then encoded as follows:
* A byte of value 1 is pushed at the start to signal that the stream has not finished. 
* The Y plane is encoded first, followed by the Cb then Cr planes.
* Each plane is encoded in 8x8 blocks.
* For P-frames, the sign of the first pixel (ie. (0,0)) is encoded in 1 bit (0 for negative, 1 otherwise). Next, the motion code is encoded using 2 bits (as described above).
* For both I-frames and P-frames, the first pixel of each block (ie (0,0)) is encoded using 16 bits representing the absolute value.
* Every following pixel is encoded as its difference from the previous pixel:
  * Differences of 0 are simply represented by a single 0 bit
  * For nonzero values, the sign is pushed first using two bits (10 for positive, 11 for negative) followed by the unary representation of the absolute value (terminated with a 0 to indicate the stop)
Finally, a 0 byte is pushed to indicate the end of the stream.
    

Sources
-----------
* My submission for Assignment 4
* Video Compression slides
